#!/usr/bin/env python3
"""
make_shelter_mesh.py -- SIH26051 parametric shelter mesh (gmsh -> CGNS/Fluent).

Builds the TEACHING/DEMO topology of docs/SIH26051_Fluent_Setup_Guide.md, fully conformal
(shared nodes -> conjugate heat transfer with zero contact resistance), and writes:

    shelter.cgns     -- read into Fluent:  File > Read > Mesh  (or PyFluent script)
    shelter.msh      -- gmsh format backup / ParaView inspection

Coordinate system (matches the guide & PyFluent script):
    x = East, y = Up (gravity -y), z = North.  South facade at z = 0.

LAYOUT CONVENTION
    LX, WZ, H are the INTERIOR clear dimensions.
    External footprint: x in [0, FX], z in [0, FZ] with FX = LX + 2*T_WALL, FZ = WZ + 2*T_WALL.
    Interior air: x in [T_WALL, T_WALL+LX], y in [T_FLOOR, T_FLOOR+H], z in [T_WALL, T_WALL+WZ].
    Roof top at y = FY = T_FLOOR + H + T_ROOF; grade plane at y = 0.

Zones created (physical groups -> Fluent cell zones / BCs):
    fluids : air_interior, win1_gain, win2_gain, trombe_gap, air_exterior
    solids : wall_steel, wall_gap, wall_puf, wall_mass, wall_lining,
             roof_steel_in, roof_puf, roof_steel_out, floor_ply, floor_puf, slab_conc,
             steel_stud, glass_win1, glass_win2, glass_trombe, soil
    BCs    : wind_inlet, outlet_all, ground_ext, facade_ext, int_envelope,
             int_glass_win1, int_glass_win2, glass_trombe_out, trombe_absorber, soil_far

DEMO simplifications (production guidance in the guide):
    * steel skins modelled as 5 mm solids (production: 2-D shell conduction for <= 2 mm)
    * glazing fills the whole wall cut (production: 40 mm IGU in a framed reveal)
    * the south steel skin doubles as the Trombe absorber plate (paint matte black)
    * no boundary-layer inflation (production: y+ ~ 1 inflation, see guide Sec. 5)

Requires:  pip install gmsh numpy       Run:  python make_shelter_mesh.py
"""
import os
import sys

import gmsh

# ============================================================ PARAMETERS (edit here) =======
LX, WZ, H = 5.6, 3.8, 2.6            # interior clear dims (E-W, N-S, height)
T_LINING, T_MASS, T_PUF, T_GAP, T_STEEL = 0.012, 0.100, 0.080, 0.025, 0.005
T_WALL = T_LINING + T_MASS + T_PUF + T_GAP + T_STEEL                  # 0.222
ROOF = [("roof_steel_in", 0.005), ("roof_puf", 0.100), ("roof_steel_out", 0.005)]
FLOOR = [("floor_ply", 0.018), ("floor_puf", 0.050), ("slab_conc", 0.100)]
T_FLOOR = sum(t for _, t in FLOOR)                                    # 0.168
T_ROOF = sum(t for _, t in ROOF)                                      # 0.110
T_SOIL = 2.0                                                          # subsoil block depth

IX0, IX1 = T_WALL, T_WALL + LX                                        # interior air x-range
IY0, IY1 = T_FLOOR, T_FLOOR + H                                       # interior air y-range
IZ0, IZ1 = T_WALL, T_WALL + WZ                                        # interior air z-range
FX, FY, FZ = LX + 2 * T_WALL, T_FLOOR + H + T_ROOF, WZ + 2 * T_WALL   # external footprint

# openings ------------------------------------------------------------------
WIN1 = dict(x=(3.3, 4.5), y=(0.9, 1.9), axis="z")    # south facade (facing -z)
WIN2 = dict(z=(1.3, 2.5), y=(0.9, 1.9), axis="x")    # east facade  (facing +x)
GAIN_T = 0.050                                        # interior 'solar gain sheet' depth
TROMBE = dict(x=(0.8, 2.8), y=(0.4, 2.4))             # on south facade
T_TGAP, T_GLAZE = 0.100, 0.010                         # Trombe air gap, glazing thickness
VENTS = [dict(x=(1.6, 2.0), y=(0.6, 0.9)),            # bottom vent (through south wall)
         dict(x=(1.6, 2.0), y=(2.0, 2.3))]            # top vent
STUD_W, STUD_PITCH = 0.10, 0.60                       # C-stud width / spacing (thermal bridges)

# exterior air and soil blocks ----------------------------------------------
EXT = dict(x=(-8.0, FX + 8.0), y=(0.0, 8.0), z=(-6.0, FZ + 6.0))
SOIL = dict(x=(-3.0, FX + 3.0), y=(-T_SOIL, 0.0), z=(-3.0, FZ + 3.0))

# mesh sizes per material ---------------------------------------------------
SIZE = dict(air_interior=0.22, win1_gain=0.08, win2_gain=0.08, trombe_gap=0.07,
            air_exterior=0.60, soil=0.45, steel_stud=0.05,
            wall_steel=0.15, wall_gap=0.12, wall_puf=0.12, wall_mass=0.15, wall_lining=0.10,
            roof_steel_in=0.15, roof_puf=0.15, roof_steel_out=0.15,
            floor_ply=0.12, floor_puf=0.12, slab_conc=0.20,
            glass_win1=0.10, glass_win2=0.10, glass_trombe=0.10)
SIZE_DEFAULT = 0.5
TOL = 1e-6

# ============================================================ geometry helpers =============

def box(name, x, y, z):
    """Add an OCC box; returns (dim3 tag, name)."""
    xmin, xmax = x; ymin, ymax = y; zmin, zmax = z
    tag = gmsh.model.occ.addBox(xmin, ymin, zmin, xmax - xmin, ymax - ymin, zmax - zmin)
    return tag, name

WALL_SPEC = [("wall_lining", T_LINING), ("wall_mass", T_MASS), ("wall_puf", T_PUF),
             ("wall_gap", T_GAP), ("wall_steel", T_STEEL)]

def layer_boxes_side(axis, side):
    """Wall stack boxes for one side ('S','N','W','E'), inner face -> outer face."""
    out, d = [], 0.0
    span_x = (0.0, FX)
    span_z = (0.0, FZ)
    if axis == "z":                                   # S: outward = -z ; N: outward = +z
        d0 = IZ0 if side == "S" else IZ1
        s = -1.0 if side == "S" else 1.0
        for nm, t in WALL_SPEC:
            za, zb = d0 + s * d, d0 + s * (d + t)
            out.append(box(nm, span_x, (0.0, FY), (min(za, zb), max(za, zb))))
            d += t
    else:                                             # W: outward = -x ; E: outward = +x
        d0 = IX0 if side == "W" else IX1
        s = -1.0 if side == "W" else 1.0
        for nm, t in WALL_SPEC:
            xa, xb = d0 + s * d, d0 + s * (d + t)
            out.append(box(nm, (min(xa, xb), max(xa, xb)), (0.0, FY), span_z))
            d += t
    return out

def layer_boxes_floorroof(kind):
    out, d = [], 0.0
    spec = FLOOR if kind == "floor" else ROOF
    if kind == "floor":                               # downward from IY0
        for nm, t in spec:
            out.append(box(nm, (0.0, FX), (IY0 - d - t, IY0 - d), (0.0, FZ)))
            d += t
    else:                                             # upward from IY1
        for nm, t in spec:
            out.append(box(nm, (0.0, FX), (IY1 + d, IY1 + d + t), (0.0, FZ)))
            d += t
    return out

def in_rect(p, r, axis):
    """Point-in-openings test; r has x/y and one of them keyed by the wall axis."""
    x, y, z = p
    if axis == "z":                                   # opening in a z-normal (S/N) wall
        return r["x"][0] - TOL < x < r["x"][1] + TOL and r["y"][0] - TOL < y < r["y"][1] + TOL
    return r["z"][0] - TOL < z < r["z"][1] + TOL and r["y"][0] - TOL < y < r["y"][1] + TOL

STUD_POSITIONS = {
    "S": [0.3 + i * STUD_PITCH for i in range(int((FX - 0.6) / STUD_PITCH) + 1)],
    "N": [0.3 + i * STUD_PITCH for i in range(int((FX - 0.6) / STUD_PITCH) + 1)],
    "W": [0.3 + i * STUD_PITCH for i in range(int((FZ - 0.6) / STUD_PITCH) + 1)],
    "E": [0.3 + i * STUD_PITCH for i in range(int((FZ - 0.6) / STUD_PITCH) + 1)],
}

# ============================================================ zone classification ==========
def classify_volume(p, bb):
    """Material name for a fragment volume from its centroid p and bbox bb
    (xmin,ymin,zmin,xmax,ymax,zmax). Fenestration rules additionally require the
    piece bounding box to FIT inside the cut region -> no centroid false-positives."""
    x, y, z = p
    x0, y0, z0, x1, y1, z1 = bb
    fits = lambda a, lo, hi: a[0] >= lo - TOL and a[1] <= hi + TOL   # (min,max) inside band
    # 1 fenestration cut regions (checked before anything else).
    # NOTE: wall layers intersected by a cut survive as solids WITH HOLES whose centroid
    # still lies inside the cut -> every axis of the bbox must fit inside the cut region.
    if in_rect(p, WIN1, "z") and fits((z0, z1), -TOL, IZ0 + TOL) \
       and fits((x0, x1), *WIN1["x"]) and fits((y0, y1), *WIN1["y"]):
        return "glass_win1"
    if in_rect(p, WIN2, "x") and fits((x0, x1), IX1 - TOL, IX1 + T_WALL + TOL) \
       and fits((z0, z1), *WIN2["z"]) and fits((y0, y1), *WIN2["y"]):
        return "glass_win2"
    # 2 interior solar-gain sheets (strictly inside the interior air box)
    if in_rect(p, WIN1, "z") and fits((z0, z1), IZ0 - TOL, IZ0 + GAIN_T + TOL) \
       and fits((x0, x1), IX0, IX1) and fits((y0, y1), IY0, IY1):
        return "win1_gain"
    if in_rect(p, WIN2, "x") and fits((x0, x1), IX1 - GAIN_T, IX1 + TOL) \
       and fits((y0, y1), IY0, IY1) and fits((z0, z1), IZ0, IZ1):
        return "win2_gain"
    # 3 Trombe vent ducts + gap + glazing
    for v in VENTS:
        if v["x"][0] - TOL < x < v["x"][1] + TOL and v["y"][0] - TOL < y < v["y"][1] + TOL \
           and fits((z0, z1), -T_TGAP, IZ0 + 0.03):
            return "trombe_gap"
    if TROMBE["x"][0] - TOL < x < TROMBE["x"][1] + TOL and TROMBE["y"][0] - TOL < y < TROMBE["y"][1] + TOL:
        if -T_TGAP - T_GLAZE - TOL < z < -T_TGAP + TOL:
            return "glass_trombe"
        if -T_TGAP - TOL < z < TOL:
            return "trombe_gap"
    # 4 interior / exterior fluids, soil
    if IX0 < x < IX1 and IY0 < y < IY1 and IZ0 < z < IZ1:
        return "air_interior"
    if y < -TOL:
        return "soil"
    above_building = y > FY - TOL
    if y > -TOL and (above_building or not (0 - TOL < x < FX + TOL and 0 - TOL < z < FZ + TOL)):
        return "air_exterior"
    # 5 steel studs (thermal bridges) inside the insulation of the 4 walls
    for side, axis in (("S", "z"), ("N", "z"), ("W", "x"), ("E", "x")):
        for coord in STUD_POSITIONS[side]:
            lo, hi = coord - STUD_W / 2, coord + STUD_W / 2
            if axis == "z" and lo < x < hi and (
               (side == "S" and IZ0 - T_MASS - T_PUF - T_TGAP < z < IZ0 - T_LINING) or
               (side == "N" and IZ1 + T_LINING < z < IZ1 + T_MASS + T_PUF + T_GAP)):
                return "steel_stud"
            if axis == "x" and lo < z < hi and (
               (side == "W" and IX0 - T_MASS - T_PUF - T_TGAP < x < IX0 - T_LINING) or
               (side == "E" and IX1 + T_LINING < x < IX1 + T_MASS + T_PUF + T_GAP)):
                return "steel_stud"
    # 6 layered shell: roof / floor / walls (measure from the interior face, inner -> outer)
    if y > IY1 - TOL:
        d, cum = y - IY1, 0.0
        for nm, t in ROOF:
            if d <= cum + t + TOL:
                return nm
            cum += t
        return ROOF[-1][0]
    if y < IY0 + TOL:
        d, cum = IY0 - y, 0.0
        for nm, t in FLOOR:
            if d <= cum + t + TOL:
                return nm
            cum += t
        return FLOOR[-1][0]
    for d0, sgn, ax in ((IZ0, -1.0, "z"), (IZ1, 1.0, "z"), (IX0, -1.0, "x"), (IX1, 1.0, "x")):
        d = sgn * ((z - d0) if ax == "z" else (x - d0))
        if 0 <= d <= T_WALL + TOL:
            cum = 0.0
            for nm, t in WALL_SPEC:
                if d <= cum + t + TOL:
                    return nm
                cum += t
            return "wall_steel"
    return "air_exterior"   # unreachable fallback

def classify_point(p):
    """Point-only variant (used for mesh-size assignment at geometry vertices)."""
    return classify_volume(p, (p[0], p[1], p[2], p[0], p[1], p[2]))

# ============================================================ boundary classification ======
def classify_face(p):
    """BC name for a boundary-face centroid (external faces only; internal faces -> None)."""
    x, y, z = p
    on = lambda v, t: abs(v - t) < 1e-4
    # far-field BCs on the exterior-air and soil boxes
    if on(x, EXT["x"][0]):
        return "wind_inlet"
    if on(x, EXT["x"][1]) or on(y, EXT["y"][1]) or on(z, EXT["z"][0]) or on(z, EXT["z"][1]):
        return "outlet_all"
    if on(y, 0.0) and not (0 - TOL < x < FX + TOL and 0 - TOL < z < FZ + TOL):
        return "ground_ext"                          # ambient ground plane (set symmetry)
    # soil far-field
    if on(y, SOIL["y"][0]) or on(x, SOIL["x"][0]) or on(x, SOIL["x"][1]) \
       or on(z, SOIL["z"][0]) or on(z, SOIL["z"][1]):
        return "soil_far"
    # Trombe glazing outer face
    if on(z, -T_TGAP - T_GLAZE):
        return "glass_trombe_out"
    # facade: outer surface of the opaque envelope (4 sides + roof top), except Trombe zone
    in_trombe = TROMBE["x"][0] - TOL < x < TROMBE["x"][1] + TOL and \
                TROMBE["y"][0] - TOL < y < TROMBE["y"][1] + TOL
    if in_trombe and on(z, 0.0):
        return "trombe_absorber"                     # steel skin = absorber plate
    if (on(z, 0.0) or on(z, FZ) or on(x, 0.0) or on(x, FX)) and -0.3 < y < FY + 0.3:
        return "facade_ext"
    if on(y, FY) and -TOL < x < FX + TOL and -TOL < z < FZ + TOL:
        return "facade_ext"
    # interior envelope faces
    if on(x, IX0) or on(x, IX1) or on(z, IZ0) or on(z, IZ1) or on(y, IY0) or on(y, IY1):
        if in_rect(p, WIN1, "z") and on(z, IZ0):
            return "int_glass_win1"
        if in_rect(p, WIN2, "x") and on(x, IX1):
            return "int_glass_win2"
        return "int_envelope"
    return None                                       # internal CHT interface -> auto-interior

# ============================================================ build ========================
def main():
    out_dir = os.path.dirname(os.path.abspath(__file__))
    cgns_path = os.path.join(out_dir, "shelter.cgns")
    msh_path = os.path.join(out_dir, "shelter.msh")

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add("sih26051_shelter")

    occ = gmsh.model.occ
    parts = []
    # fluids & ground ---------------------------------------------------------
    parts.append(box("air_interior", (IX0, IX1), (IY0, IY1), (IZ0, IZ1)))
    parts.append(box("air_exterior", **EXT))
    parts.append(box("soil", **SOIL))
    # walls, roof, floor ------------------------------------------------------
    for side in ("S", "N"):
        parts += layer_boxes_side("z", side)
    for side in ("W", "E"):
        parts += layer_boxes_side("x", side)
    parts += layer_boxes_floorroof("floor")
    parts += layer_boxes_floorroof("roof")
    # studs (thermal bridges through insulation) ------------------------------
    for side in ("S", "N"):
        for c in STUD_POSITIONS[side]:
            lo, hi = c - STUD_W / 2, c + STUD_W / 2
            if side == "S":
                parts.append(box("steel_stud", (lo, hi), (0.0, FY),
                                 (IZ0 - T_MASS - T_PUF - T_TGAP, IZ0 - T_LINING)))
            else:
                parts.append(box("steel_stud", (lo, hi), (0.0, FY),
                                 (IZ1 + T_LINING, IZ1 + T_MASS + T_PUF + T_GAP)))
    for side in ("W", "E"):
        for c in STUD_POSITIONS[side]:
            lo, hi = c - STUD_W / 2, c + STUD_W / 2
            if side == "W":
                parts.append(box("steel_stud", (IX0 - T_MASS - T_PUF - T_TGAP, IX0 - T_LINING),
                                 (0.0, FY), (lo, hi)))
            else:
                parts.append(box("steel_stud", (IX1 + T_LINING, IX1 + T_MASS + T_PUF + T_GAP),
                                 (0.0, FY), (lo, hi)))
    # windows: cut through the wall stack + interior gain sheets ---------------
    parts.append(box("glass_win1", (WIN1["x"][0], WIN1["x"][1]), (WIN1["y"][0], WIN1["y"][1]),
                     (0.0, IZ0)))
    parts.append(box("win1_gain", (WIN1["x"][0], WIN1["x"][1]), (WIN1["y"][0], WIN1["y"][1]),
                     (IZ0, IZ0 + GAIN_T)))
    parts.append(box("glass_win2", (IX1, IX1 + T_WALL), (WIN2["y"][0], WIN2["y"][1]),
                     (WIN2["z"][0], WIN2["z"][1])))
    parts.append(box("win2_gain", (IX1 - GAIN_T, IX1), (WIN2["y"][0], WIN2["y"][1]),
                     (WIN2["z"][0], WIN2["z"][1])))
    # Trombe gap + glazing (on the south facade) -------------------------------
    parts.append(box("trombe_gap", (TROMBE["x"][0], TROMBE["x"][1]),
                     (TROMBE["y"][0], TROMBE["y"][1]), (-T_TGAP, 0.0)))
    parts.append(box("glass_trombe", (TROMBE["x"][0], TROMBE["x"][1]),
                     (TROMBE["y"][0], TROMBE["y"][1]), (-T_TGAP - T_GLAZE, -T_TGAP)))
    # vent ducts (gap <-> interior, through the south wall) --------------------
    for v in VENTS:
        parts.append(box("trombe_gap", (v["x"][0], v["x"][1]), (v["y"][0], v["y"][1]),
                         (-T_TGAP, IZ0 + 0.03)))

    print(f"fragmenting {len(parts)} boxes ...")
    frags, _ = occ.fragment([(3, t) for t, _ in parts], [])
    occ.removeAllDuplicates()
    occ.synchronize()

    # ---- classify fragment volumes by centroid ------------------------------
    zones = {}
    for dim, tag in frags:
        if dim != 3:
            continue
        mass = gmsh.model.occ.getMass(3, tag)
        if mass < 1e-9:                               # sliver: attach to exterior air
            zones.setdefault("air_exterior", []).append(tag)
            continue
        cx, cy, cz = gmsh.model.occ.getCenterOfMass(3, tag)
        bb = gmsh.model.getBoundingBox(3, tag)
        mat = classify_volume((cx, cy, cz), bb)
        zones.setdefault(mat, []).append(tag)

    print("\nvolume zones:")
    total = 0.0
    for name in sorted(zones):
        vtot = sum(gmsh.model.occ.getMass(3, t) for t in zones[name])
        total += vtot
        print(f"  {name:<15s} n={len(zones[name]):4d}  V={vtot:10.3f} m3")
    print(f"  {'TOTAL':<15s}        V={total:10.3f} m3")

    # ---- physical groups (volumes) ------------------------------------------
    for name, tags in zones.items():
        p = gmsh.model.addPhysicalGroup(3, tags)
        gmsh.model.setPhysicalName(3, p, name)

    # ---- physical groups (external boundary faces) --------------------------
    bcs = {}
    for dim, tag in gmsh.model.getEntities(2):
        com = gmsh.model.occ.getCenterOfMass(2, tag)
        bc = classify_face(com)
        if bc:
            bcs.setdefault(bc, []).append(tag)
    print("\nboundary groups:")
    for name in sorted(bcs):
        print(f"  {name:<18s} faces={len(bcs[name])}")
    for name, tags in bcs.items():
        p = gmsh.model.addPhysicalGroup(2, tags)
        gmsh.model.setPhysicalName(2, p, name)

    # ---- mesh size field from material classification at geometry points ----
    pts = gmsh.model.getEntities(0)
    by_size = {}
    for dim, tag in pts:
        xyz = gmsh.model.getValue(0, tag, [])
        mat = classify_point(tuple(xyz))
        by_size.setdefault(SIZE.get(mat, SIZE_DEFAULT), []).append((dim, tag))
    for s, tags in by_size.items():
        gmsh.model.mesh.setSize(tags, s)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.Algorithm", 6)        # 2-D frontal-Delaunay
    gmsh.option.setNumber("Mesh.Algorithm3D", 1)      # 3-D Delaunay (HXT=10 can crash on
    gmsh.option.setNumber("Mesh.MaxNumThreads1D", 1)  # heavily fragmented multi-material
    gmsh.option.setNumber("Mesh.MaxNumThreads2D", 1)  # OCC assemblies; keep single-threaded
    gmsh.option.setNumber("Mesh.MaxNumThreads3D", 1)
    gmsh.option.setNumber("Mesh.Optimize", 1)
    gmsh.option.setNumber("Mesh.OptimizeThreshold", 0.3)

    print("\ngenerating mesh ...")
    for algo, name in ((1, "Delaunay"), (4, "Frontal"), (10, "HXT")):
        try:
            gmsh.option.setNumber("Mesh.Algorithm3D", algo)
            gmsh.model.mesh.generate(3)
            break
        except Exception as exc:
            print(f"3-D algorithm {name} failed ({exc}); trying next ...")
            gmsh.model.mesh.clear()
    else:
        sys.exit("all 3-D algorithms failed - coarsen sizes or simplify the assembly")
    gmsh.model.mesh.removeDuplicateNodes()

    n_el = len(gmsh.model.mesh.getElementsByType(4)[0])          # tets
    print(f"tetrahedra: {n_el}")

    # ---- export --------------------------------------------------------------
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.write(msh_path)
    try:
        gmsh.write(cgns_path)                          # needs CGNS support in the gmsh build
        print(f"wrote {cgns_path}")
        cgns_ok = True
    except Exception as exc:                           # pragma: no cover
        print(f"CGNS export failed ({exc}); use shelter.msh with a gmsh->Fluent converter "
              "or mesh via SpaceClaim/Fluent Meshing.")
        cgns_ok = False
    print(f"wrote {msh_path}")
    gmsh.finalize()

    print("\nNEXT STEPS")
    print(" 1. Fluent: File > Read > Mesh -> shelter.cgns" if cgns_ok else
          " 1. Convert shelter.msh (gmsh) to CGNS via meshio, or rebuild in SpaceClaim")
    print(" 2. automation/pyfluent_run.py  (full setup + transient run)")
    print(" 3. After loading, read zone IDs (TUI /mesh/zones) -> paste into climate_data.h")

if __name__ == "__main__":
    main()
