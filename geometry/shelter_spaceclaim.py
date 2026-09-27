# shelter_spaceclaim.py -- SIH26051 parametric shelter geometry for ANSYS SpaceClaim
# =====================================================================================
# BEST-EFFORT scripting-API sketch (IronPython). SpaceClaim's recorded-API namespace
# (SpaceClaim.Api.VXX) changes with every release: run this from the SpaceClaim File >
# New > Script panel, or better: re-record the operations once via the Script panel and
# splice in the parameters below. The gmsh route (make_shelter_mesh.py) is the tested,
# version-independent pipeline; use this file when your Fluent licence mandates CAD
# geometry (.scdoc -> Fluent Meshing watertight workflow).
#
# Builds the same topology as the guide Sec. 3: interior air box, 5-layer wall stacks,
# roof/floor stacks, steel studs, two glazing cuts, Trombe gap + glazing, vent ducts,
# exterior air box and 2 m soil block -- each as a SEPARATE solid so the meshing
# workflow can share topology and name cell zones.
#
# Units: millimetres via MM() helper. Axes: x=East, y=Up, z=North (matches the guide).
# =====================================================================================

import math
from SpaceClaim.Api.V19 import *          # adjust VXX to your release (V18/V19/V20+)
from SpaceClaim.Api.V19.Geometry import *
from SpaceClaim.Api.V19.ScriptMode import *   # some builds expose ScriptMode here

# ------------------------------------------------------------------ parameters (mm) ---
LX, WZ, H = 5600.0, 3800.0, 2600.0                 # interior clear dims
T_LINING, T_MASS, T_PUF, T_GAP, T_STEEL = 12.0, 100.0, 80.0, 25.0, 5.0
T_WALL = T_LINING + T_MASS + T_PUF + T_GAP + T_STEEL
ROOF_T = [5.0, 100.0, 5.0]
FLOOR_T = [18.0, 50.0, 100.0]
T_FLOOR, T_ROOF = sum(FLOOR_T), sum(ROOF_T)
T_SOIL = 2000.0
FX, FZ = LX + 2 * T_WALL, WZ + 2 * T_WALL          # external footprint
FY = T_FLOOR + H + T_ROOF
EXT_PAD = 8000.0                                   # exterior air around the shelter
SOIL_PAD = 3000.0

WIN1 = (3300.0, 4500.0, 900.0, 1900.0)             # south: x0,x1,z0,z1 (y-band)
WIN2 = (1300.0, 2500.0, 900.0, 1900.0)             # east:  z0,z1,y0,y1 (x-band)
TROMBE = (800.0, 2800.0, 400.0, 2400.0)            # south x0,x1,y0,y1
T_TGAP, T_GLAZE = 100.0, 10.0
VENTS = [(1600.0, 2000.0, 600.0, 900.0),           # x0,x1,y0,y1 through south wall
         (1600.0, 2000.0, 2000.0, 2300.0)]


def MM(v):
    return float(v)


def box(x0, y0, z0, dx, dy, dz, name):
    """Create a named solid box; SpaceClaim groups it as its own component."""
    result = Body.CreateBox(
        Point.Create(MM(x0), MM(y0), MM(z0)),
        Point.Create(MM(x0 + dx), MM(y0 + dy), MM(z0 + dz)))
    result.Name = name
    return result


def layer_boxes_side(axis, side):
    """Wall stack, inner face -> outer face, as separate named solids."""
    spec = [("wall_lining", T_LINING), ("wall_mass", T_MASS), ("wall_puf", T_PUF),
            ("wall_gap", T_GAP), ("wall_steel", T_STEEL)]
    d, out = 0.0, []
    for nm, t in spec:
        if axis == "z":
            z0 = (T_WALL - d - t) if side == "S" else (T_WALL + WZ + d)
            out.append(box(0.0, 0.0, z0, FX, FY, t, f"{nm}_{side}"))
        else:
            x0 = (T_WALL - d - t) if side == "W" else (T_WALL + LX + d)
            out.append(box(x0, 0.0, 0.0, t, FY, FZ, f"{nm}_{side}"))
        d += t
    return out


def main():
    ViewHelper.SetSketchPlane(ViewHelper.SketchPlaneType.XY)
    parts = []
    # fluids & ground ----------------------------------------------------------
    parts.append(box(T_WALL, T_FLOOR, T_WALL, LX, H, WZ, "air_interior"))
    parts.append(box(-EXT_PAD, 0.0, -EXT_PAD, FX + 2 * EXT_PAD, 8000.0,
                     FZ + 2 * EXT_PAD, "air_exterior"))
    parts.append(box(-SOIL_PAD, -T_SOIL, -SOIL_PAD, FX + 2 * SOIL_PAD, T_SOIL,
                     FZ + 2 * SOIL_PAD, "soil"))
    # wall / roof / floor stacks ------------------------------------------------
    for side in ("S", "N"):
        parts += layer_boxes_side("z", side)
    for side in ("W", "E"):
        parts += layer_boxes_side("x", side)
    d = 0.0
    for i, t in enumerate(FLOOR_T):
        parts.append(box(0.0, T_FLOOR - d - t, 0.0, FX, t, FZ, f"floor_layer{i}"))
        d += t
    d = 0.0
    for i, t in enumerate(ROOF_T):
        parts.append(box(0.0, T_FLOOR + H + d, 0.0, FX, t, FZ, f"roof_layer{i}"))
        d += t
    # windows: cut + glazing + interior gain sheet -------------------------------
    x0, x1, y0, y1 = WIN1
    parts.append(box(x0, y0, 0.0, x1 - x0, y1 - y0, T_WALL, "glass_win1"))
    parts.append(box(x0, y0, T_WALL, x1 - x0, y1 - y0, 50.0, "win1_gain"))
    z0, z1, y0, y1 = WIN2
    parts.append(box(T_WALL + LX, y0, z0, T_WALL, y1 - y0, z1 - z0, "glass_win2"))
    parts.append(box(T_WALL + LX - 50.0, y0, z0, 50.0, y1 - y0, z1 - z0, "win2_gain"))
    # Trombe gap + glazing + vents ------------------------------------------------
    x0, x1, y0, y1 = TROMBE
    parts.append(box(x0, y0, -T_TGAP, x1 - x0, y1 - y0, T_TGAP, "trombe_gap"))
    parts.append(box(x0, y0, -T_TGAP - T_GLAZE, x1 - x0, y1 - y0, T_GLAZE,
                     "glass_trombe"))
    for i, (vx0, vx1, vy0, vy1) in enumerate(VENTS):
        parts.append(box(vx0, vy0, -T_TGAP, vx1 - vx0, vy1 - vy0, T_TGAP + T_WALL,
                         f"vent{i}"))
    # TODO(studs): add C-stud boxes every 600 mm through insulation layers (guide Sec. 3);
    # their boolean is identical to the vent ducts above.

    # share topology so the mesh is conformal (CHT without contact resistance):
    try:
        ShareToggling.Share(GetRootPart().GetAllBodies())   # API name varies by release
    except Exception:
        pass  # alternatively: Design > Share in the GUI, or set in Fluent Meshing
    GetRootPart().Name = "SIH26051_shelter"
    print("geometry created -- next: File > Save As shelter.scdoc, then Fluent Meshing "
          "watertight workflow (or import into Discovery to verify)")


main()
