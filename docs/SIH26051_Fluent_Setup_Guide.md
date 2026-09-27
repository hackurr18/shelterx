# SIH26051 — Transient CFD Model of an Area-Specific Passive Shelter (ANSYS Fluent)
## Complete Setup Guide · DRDO High-Altitude Cold Region (Ladakh Baseline)

**Scope:** Software-based model development for design of an area-specific shelter for thermal comfort
maintenance (SIH26051, DRDO / Dept. of Defence Production–iDEX, Software Edition).
This document is the engineering ground truth for the simulation pipeline: geometry, physics,
solver, climate-data ingestion, meshing, thermal-loss accounting, comfort metrics, and the
PyFluent automation that turns the model into the required *user-friendly predictive tool*.

**Companion artifacts in this package**

| File | Purpose |
|---|---|
| `udf/climate_shelter_udf.c` | UDF suite: EPW ingestion, wind profile, sky temp, dynamic-SHGC solar gain, run logger |
| `udf/climate_data.h` | Generated hourly climate arrays (example: Leh design day, 21 Dec) |
| `tools/make_inputs.py` | EPW → `climate_data.h` + Fluent transient profile + transmitted-gain tables |
| `geometry/make_shelter_mesh.py` | Parametric gmsh script → conformal multi-zone mesh (CGNS for Fluent) |
| `geometry/shelter_spaceclaim.py` | Parametric SpaceClaim (IronPython) geometry alternative |
| `automation/pyfluent_run.py` | PyFluent: launch → setup → UDF hook → DOE batch → KPI ranking (CSV) |

---

## 0. Requirement → Model Mapping (problem-statement traceability)

| PS requirement | Where implemented |
|---|---|
| Prediction of shelter inside temperature from user-defined inputs | Transient CFT run; `T_in` monitors; DOE loop in `pyfluent_run.py` |
| Prediction of thermal energy generated from solar radiation | Solar Ray Tracing + `DEFINE_SOLAR_INTENSITY` (EPW DNI/DHI); transmitted-gain tables (dynamic SHGC) |
| Heat-flow details vs (T_amb − T_shelter) over time | Flux report definitions per zone; UDF hourly CSV audit; §9 |
| Multi-layer walls, thermal mass, composite multi-material | 5-layer wall stack + PCM option, §3–4 |
| Effect of openings, size, shape, orientation | Parametric geometry (gmsh/SpaceClaim scripts); orientation-aware gain tables |
| Comparative material analysis → most efficient combination | `automation/pyfluent_run.py` design-of-experiments + ranking CSV |
| Minimize energy use / passive self-sufficiency | Comfort-hours-below-setpoint KPI + net heating demand proxy, §12 |

**3-day fast-track (submission 30 Sep 2026):**
1. Day 1 — run `make_shelter_mesh.py` (demo topology, sealed interior), `pyfluent_run.py RUN_MODE="demo"` → one 24-h cycle.
2. Day 2 — validate vs hand-calc U-values (§13), produce 3-day spin-up + recording day, charts.
3. Day 3 — DOE of 4–6 wall/glazing variants, ranking table, slide deck. Full external-wind domain is a stretch goal; the sealed-envelope model with convective/radiative external BCs (§6.4, Path B) is engineering-acceptable for the PS.

---

## 1. Workflow Overview

```
EPW file ──► make_inputs.py ──┬─► udf/climate_data.h          (hourly T, DNI, DHI, wind, gain tables)
                              ├─► profiles/amb_profile.pro    (Fluent transient profile, UDF-free option)
                              └─► KPI constants (window areas, SHGC schedule)

gmsh / SpaceClaim ──► conformal mesh (interior air | exterior air | Trombe gap |
                                       5-layer walls | glazing | studs | slab | 2 m soil)
        │
        ▼
PyFluent: models (energy, Realizable k-ε + EWT, DO + Solar Ray Tracing)
          BCs ← UDF/profiles  ·  initialization  ·  spin-up ≥5 diurnal cycles
        │
        ▼
Recording day(s): hourly CSV (T_in stats, heat-loss split) + fields (T, q", MRT)
        │
        ▼
Post: time lag φ, decrement factor f, MRT_solar, PMV/PPD, loss breakdown → ranking
```

Coordinate convention (used everywhere below): **x = East, y = Up (gravity −y), z = North.**
Shelter south face at z = 0. Fluent solar-calculator mesh orientation: North (0,0,1), East (1,0,0).

---

## 2. Governing Physics & Key Equations

### 2.1 Flow (transient RANS, Boussinesq interior)
Continuity and momentum (i-component):

```
∂ρ/∂t + ∇·(ρV) = 0
∂(ρu_i)/∂t + ∂(ρ u_i u_j)/∂x_j = −∂p/∂x_i + ∂τ_ij/∂x_j + (ρ − ρ_ref) g_i
```

Boussinesq closure for the interior air (|ΔT| ≲ 15 K):

```
ρ = ρ_ref [1 − β (T − T_ref)],   β = 1/T_ref   (ideal gas at T_ref)
```

Exterior air: incompressible-ideal-gas with **operating density = ρ(T_amb)** so hydrostatic
pressure cancels — this avoids spurious buoyancy from the reference-pressure offset.

### 2.2 Energy (conjugate heat transfer)
Fluid: `∂(ρE)/∂t + ∇·(V(ρE+p)) = ∇·(k_eff ∇T) + S_h`
Solid (each wall/soil layer): `∂(ρ c T)/∂t = ∇·(k ∇T)`
Effective conductivity: `k_eff = k + c_p μ_t/Pr_t`, Pr_t = 0.85.

### 2.3 Realizable k-ε + Enhanced Wall Treatment (EWT)

```
∂(ρk)/∂t + ∇·(ρkV) = ∇·[(μ + μ_t/σ_k)∇k] + G_k + G_b − ρε − Y_M
∂(ρε)/∂t + ∇·(ρεV) = ∇·[(μ + μ_t/σ_ε)∇ε] + ρC₁Sε − ρC₂ ε²/(k + √(νε)) + C₁ε ε/k C₃ε G_b
μ_t = ρ C_μ k²/ε,  C_μ = 1/(A₀ + A_s U* k/ε),  C₁ = max(0.43, η/(η+5)),  C₂ = 1.9
```
EWT: two-layer formulation with blending function Re_y; valid y+ ≈ 1 (also tolerates
y+ 30–300 through automatic blending). Target **y+ ≈ 1** on all interior and exterior
surfaces because natural-convection Nusselt numbers dominate the answer.

### 2.4 Discrete Ordinates radiation (RTE) + Solar Ray Tracing

```
dI(r,s)/ds + (κ + σ_s) I(r,s) = κ I_b(T) + σ_s/4π ∫ I(r,s′) Φ(s·s′) dΩ′
```
- Shelter air is treated as **non-participating** (κ = σ_s = 0); radiation is surface-to-surface + solar.
- Glazing uses **two-band** semi-transparent DO properties: solar band 0.25–2.7 µm (τ_vis, ρ, α) and
  IR band > 2.7 µm (τ_IR ≈ 0, ε ≈ 0.9) → reproduces the visible/IR split that drives SHGC.
- Solar Ray Tracing overlays the direct beam: shadow rays project sun-load through openings and onto
  facades; diffuse sky + ground-reflected applied per ISO-type view factors.
- Solar geometry (used by Fluent's solar calculator, and by `make_inputs.py`):

```
declination   δ = 23.45° sin[360°(284+n)/365]
hour angle    ω = 15°(t_solar − 12)
altitude      α_s:  sin α_s = sin φ sin δ + cos φ cos δ cos ω     (φ = latitude)
incidence on vertical wall of azimuth ψ:  cosθ = cos α_s · cos(ψ_sun − ψ_wall)
```

### 2.5 Envelope conduction & the numbers to beat
Series thermal resistance (ISO 6946; R_si = 0.13, R_se = 0.04 m²K/W):

```
U = 1 / ( R_si + Σ t_i/k_i + R_se )            [W/m²K]
Sol-air temperature:  T_sa = T_amb + α_abs G/h_o − ε ΔR/h_o
Exterior film (McAdams):  h_o = 5.7 + 3.8 V_wind   [W/m²K]
Sky temperature (Swinbank, clear):  T_sky = 0.0552 · T_amb[K]^1.5
```

Baseline wall (GI 1.5 mm / air gap 25 mm / PUF 80 mm / concrete 100 mm / plywood 12 mm):
R = 0.13 + 0.03 + 0.18 + 3.64 + 0.071 + 0.092 + 0.04 = **4.19 → U = 0.239 W/m²K**.
This hand value is the §13 validation anchor for the CFD-conjugate result.

### 2.6 Thermal bridges (ISO 10211/14683)
```
linear:   ψ = ( Q_2D − Σ_i U_i l_i ΔT ) / ΔT        [W/mK]
point:    χ = ( Q_3D − Σ U_i A_i ΔT − Σ ψ_j L_j ΔT ) / ΔT   [W/K]
```
In Fluent: integrate total surface heat flux over the framed wall zone, subtract the clear-wall
U·A·ΔT product (same ΔT window) → per-frame ψ contribution. Steel C-studs at 600 mm c/c typically
add 20–40 % to an 80 mm PUF wall's transmission loss — a headline result for the DRDO report.

### 2.7 Ground coupling
Soil column R = d/k_soil = 2.0/1.5 = 1.33 m²K/W. Slab + soil in series → effective floor loss
coefficient ≈ 0.35 W/m²K referenced to deep-ground temperature. The 2 m depth is adequate because
the annual wave decays as exp(−z/dₐ), dₐ = √(2α/ω) ≈ 2.8 m for α = 8×10⁻⁷ m²/s: at 2 m the
amplitude is damped ~50 % and lagged ~2 months; the residual swing is captured since the soil is
meshed as a transient solid, with only the *deep* face driven by data (constant or annual sinusoid).

### 2.8 Occupant comfort metrics
Mean radiant temperature from view factors (exact, post-processed):
```
T_MRT⁴ = Σ_i F_i T_i⁴ ,   Σ F_i = 1
```
Field approximation inside the Fluent domain, from DO incident radiation G (isotropic enclosure):
```
T_MRT ≈ ( G / (4σ) )^0.25 ,   σ = 5.670e-8
```
Solar-adjusted MRT (ASHRAE-55 / Arens et al., as implemented in the CBE tool):
```
MRT_sol = MRT + ΔT_sw
ΔT_sw = F_sun · (α_sw/ε_lw) · f_sw · I_person / h_r
h_r = ε σ (T_MRT + T_a)(T_MRT² + T_a²) ≈ 4.6 W/m²K
α_sw ≈ 0.7 (skin/clothing shortwave absorptivity), ε_lw ≈ 0.95,
F_sun = fraction of body in sun (0–1), f_sw ≈ 0.22 (projected-area factor, seated),
I_person = shortwave irradiance at occupant position from DO/solar model
```
Time lag & decrement factor (walls, roof — the PS's "diurnal time lag"):
```
φ = t_peak(T_interior surface) − t_peak(T_sol-air)          [hours]
f  = (T_i,max − T̄_i) / (T_o,max − T̄_o)                      (want φ large, f small)
```
PMV/PPD per ISO 7730 — Fluent's Thermal Comfort post-processing computes these directly once DO
is on (needs metabolic rate M ≈ 70–120 W, clothing I_cl ≈ 1.0–1.4 clo, humidity input):
```
PMV = (0.303 e^(−0.036M) + 0.028){M−W − 3.05e−3[5733 − 6.99(M−W) − p_a]
      − 0.42[(M−W) − 58.15] − 1.7e−5 M(5867 − p_a) − 0.0014 M(34 − t_a)
      − 3.96e−8 f_cl[(t_cl+273)⁴ − (T_MRT+273)⁴] − f_cl h_c (t_cl − t_a)}
```
Comfort band for the DRDO brief: **PMV ∈ [−0.7, +0.7]** ⇔ ~18–24 °C operative for 1.2 clo.

---

## 3. Geometry & Computational Domain (parametric defaults)

### 3.1 Dimensions (edit in one place — scripts read the same dict)
```
Shelter internal: 5.6 (E-W, x) × 3.8 (N-S, z) × 2.6 m high   → 21.3 m² floor, 55 m³
Wall stack (outer→inner): GI steel 1.5 mm | air gap 25 mm | PUF 80 mm | concrete 100 mm | plywood 12 mm
Roof:  steel 1.5 | air 50 | PUF 100 | steel 1.5                      U ≈ 0.21
Floor: ply 18 | PUF 50 | concrete slab 100 (on steel joists)         U_slab ≈ 0.32
Window south (win1): 1.2 × 1.0 m; window east (win2): 1.2 × 1.0 m;
       triple-glazed low-e, U_g = 1.0, SHGC 0.25–0.60 dynamic (electrochromic)
Door (north): 0.9 × 2.0 m insulated, U = 1.5, airlock excluded (leakage BC instead)
Trombe wall (south, x∈[0.6,3.6]): glazing | 100 mm gap | concrete mass 150 mm, black;
       top & bottom vents 0.4 × 0.2 m with dampers
Shading louver: horizontal, 0.6 m projection over south windows (geometry toggle)
Steel C-studs: 100 mm, 600 mm c/c, embedded through insulation (thermal-bridge solids)
Exterior air domain: x∈[−12,18], z∈[−7,11], y∈[0,10]  (~2–3 H around, 3 H above)
Ground block: footprint + 3 m margin, 2.0 m deep (y∈[−2,0])
```

### 3.2 ASCII section (south face, x–z cut not to scale)
```
                sky  ·  DO/solar rays  ☀
   ┌────────────────────────────────────────────────────────┐ exterior air (wind inlet W)
   │   ┌────────────────────────────────────────────┐      │
   │   │ roof stack (steel|gap|PUF|steel)           │      │
   │   ├────────────────────────────────────────────┤      │
   │ W │                                            │ E    │
   │ a │            INTERIOR AIR                    │ i    │
   │ l │        (Boussinesq, 55 m³)                 │ n    │
   │ l │   T(x,y,z,t)  ·  MRT plane y = 1.1 m       │ d    │
   │   ├──────┬──────────┬─────────────────────────┤      │
   │   │window│ Trombe:  │ south wall stack        │      │
   │   │glass │ gap+mass │ (5 layers + studs)      │      │
   └───┴──────┴──────────┴─────────────────────────┴──────┘
  ═══════════ floor slab (ply|PUF|concrete) ═══════════
  ▓▓▓▓▓▓▓▓▓▓▓▓▓ soil block, 2.0 m deep ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
       deep face: T_ground,deep (UDF, from data, e.g. +8 °C)
```

### 3.3 Zone-naming contract (scripts depend on these exact names)
| Zone (cell zone unless noted) | Type | Material |
|---|---|---|
| `air_interior` | fluid | air-boussinesq |
| `air_exterior` | fluid | air-incompressible-ideal-gas |
| `trombe_gap` | fluid | air-boussinesq |
| `wall_steel_out`, `wall_gap`, `wall_puf`, `wall_mass`, `wall_lining` | solids | per §4 |
| `roof_*`, `floor_ply`, `floor_puf`, `slab_conc` | solids | per §4 |
| `soil` | solid | soil-wet |
| `stud_steel` | solid | steel (thermal bridge) |
| `glass_win1`, `glass_win2`, `glass_trombe` | solids | glazing |
| `trombe_mass` | solid | concrete-black |
| `win1_gain`, `win2_gain` (thin fluid sheets inside `air_interior`) | fluid | air (solar-gain source zones) |

Boundary walls (surfaces): `out_wall_s/e/w/n`, `out_roof`, `out_soil_sides/bottom`,
`in_wall_*`, `in_roof`, `in_floor`, `in_glass_*`, `trombe_absorber`, `door_n`, `vent_top`, `vent_bottom`.

---

## 4. Material Database (SI units, paste into Fluent)

| Name | ρ [kg/m³] | c_p [J/kgK] | k [W/mK] | Notes |
|---|---|---|---|---|
| steel-gi | 7850 | 460 | 50 | sheets ≤2 mm → shell conduction in production |
| puf-80 | 35 | 1400 | 0.022 | spray PUF |
| vip-30 | 200 | 1000 | 0.007 | vacuum panel (DOE variant) |
| eps-100 | 25 | 1400 | 0.034 | lightweight variant |
| concrete-mass | 2300 | 880 | 1.4 | thermal mass / Trombe |
| pcm-rt24 | 800 | 2000 | 0.2 | **Enthalpy–porosity**: T_liq 24 °C, T_sol 22 °C, L = 190 kJ/kg |
| plywood | 550 | 2100 | 0.13 | lining |
| glazing-lowE | 2500 | 840 | 1.0 | DO two-band semi-transparent |
| air-boussinesq | 1.34 (T_ref) | 1006 | 0.024 | β = 1/T_ref |
| soil-wet | 1900 | 1000 | 1.5 | frozen crust variant k = 2.0 |

DOE wall variants (comparative analysis demanded by the PS): W1 = EPS-100 (light), W2 = baseline
PUF+concrete, W3 = PUF + PCM-RT24, W4 = VIP-30 + concrete-60, W5 = adobe/stone 300 mm
(k=1.05, ρ=1900, c=850 — Ladakh vernacular benchmark).

---

## 5. Meshing Guidelines

### 5.1 Topology & connectivity
- **Conformal mesh across every solid–solid and solid–fluid interface** (shared topology in SpaceClaim;
  OCC `fragment` in the gmsh script). Non-conformal interfaces add resistance & noise — avoid.
- Fluid regions: poly-hexcore or trimmed hex. Solids (layers, slab, soil): pure hex, swept/graded.
- Thin sheets ≤ 2 mm (GI skins): **do not mesh in 3D** — use Fluent **shell conduction** (2-D wall
  layers) in production. The demo gmsh mesh models them as 5 mm solids only to keep a fully
  watertight teaching mesh.

### 5.2 Sizes (production targets)
| Region | Min size | Max size | Cells through thickness |
|---|---|---|---|
| Wall layers | 20 mm | 150 mm | PUF ≥ 4, concrete ≥ 4, lining ≥ 2, air gap ≥ 2 |
| Steel studs | 15 mm | 40 mm | ≥ 4 across web |
| Glazing | 10 mm | 100 mm | ≥ 2 |
| Interior air | 40 mm | 250 mm | inflation 1st layer 1 mm, 12 layers, GR 1.2 (y+≈1) |
| Exterior air | 100 mm | 600 mm | inflation 1st layer 0.1 mm, 20 layers, GR 1.2 (y+≈1 at V=3 m/s) |
| Trombe gap | 20 mm | 60 mm | ≥ 5 across 100 mm gap |
| Soil | 25 mm at slab → 200 mm | — | ≥ 12 graded layers over 2 m |

y+ sizing check (exterior, U = 3 m/s, L = 5.6 m, air at −10 °C): Re_L ≈ 1.1×10⁶ →
c_f ≈ 0.0035, u_τ = U√(c_f/2) ≈ 0.126 m/s → **Δ₁ = y⁺ μ/(ρ u_τ) ≈ 0.10 mm**.
Interior natural convection (u ≈ 0.3 m/s): Δ₁ ≈ 0.6–1.0 mm.
**Always re-check with `/mesh/check-mesh` + y⁺ histograms after the first 100 iterations.**

### 5.3 Refinement must-haves
Local sizing (body-of-influence or curvature) at: window reveals (10 mm), Trombe vents (10 mm),
louver gaps, door seal line, stud/wall junctions, and the win*_gain source sheets (20 mm).
Snow ground reflectivity ≠ geometry, but remember exterior near-ground cells see high albedo (§6.6).

### 5.4 Quality gates (Fluent Meshing / gmsh `gmsh.model.mesh.quality`)
- Skewness ≤ 0.95 (target ≤ 0.90 poly), orthogonal quality ≥ 0.10 (target 0.15)
- Max aspect ratio: fluid ≤ 100 (inflation 1000 acceptable if orthogonal), solids ≤ 50
- Zero negative volumes; 100 % interface coverage; `mesh/check` warnings all cleared
- Volume budget: production 5–8 M cells; demo mesh from `make_shelter_mesh.py` ≤ 2 M cells.

---

## 6. Step-by-Step Solver Setup (2023R2+; TUI given where GUI menus drift)

**S1. Read mesh.** `File → Read → Mesh` (`shelter.cgns` from the gmsh script, or your .msh).
`/mesh/check` → fix scale (CGNS arrives in metres), `/mesh/check-verbosity` for interface stats.

**S2. General.** Solver: **Pressure-Based, Transient**. Gravity ON, `Y = −9.81`.
`/define/operating-conditions operating-density` = ρ(T_amb) ≈ 1.34 kg/m³ (exterior reference).

**S3. Models.**
- Energy: ON (`/define/models/energy yes`).
- Viscous: **k-epsilon → Realizable**, Near-Wall Treatment: **Enhanced Wall Treatment**
  (`/define/models/viscous/model/ke-realizable yes` → near-wall `?` enhanced). Enable
  **viscous heating** OFF, thermal effects: include buoyancy effects in k & ε production (G_b term ON).
- Radiation: **DO** (`/define/models/radiation/model/do`): θ-divisions 4, φ-divisions 8
  (raise to 8×16 for MRT-quality fields), iterations-per-radiation 1, convergence 1e-6
  energy-residual-linked. **Solar Load: Solar Ray Tracing** (`/define/models/radiation/solar-load`).
  In the Solar Calculator: Latitude 34.16, Longitude 77.58, TZ +5.5, North (0,0,1), East (1,0,0),
  irradiation method **Fair Weather** (values overridden by UDF, §7), spectral fraction 0.5.
- Set **Time-steps-per-solar-load-update = 1** so the sun walks across the sky each step.

**S4. Materials.** Create §4 table via GUI or the `pyfluent_run.py` `MATERIALS` dict. Enable
enthalpy–porosity (Solidification/Melting) **only** if PCM variant is active.

**S5. Cell zones.** Assign materials (§3.3). `air_interior` & `trombe_gap`: Boussinesq,
T_ref = initial interior T (e.g. 263 K). `air_exterior`: incompressible-ideal-gas,
operating density as S2. `win1_gain`/`win2_gain`: same Boussinesq air; these thin sheets only host
the SHGC energy source (S8).

**S6. Climate ingestion — UDF.** Compile & load (`/define/user-defined/compiled-functions compile
"libudf" → load`). Hook, per §7 table: `ambient_temp_profile`, `wind_speed_profile`, `sky_temp_profile`,
`ground_deep_profile`, `direct_from_epw` + `diffuse_from_epw` (Radiation Model dialog, Direct/Diffuse
Solar Irradiation dropdowns → "udf"), `win1_solar_gain`, `win2_solar_gain` (cell-zone sources).
*Profile-table alternative (no compiler on judges' laptop):* `profiles/amb_profile.pro` loaded via
`File → Read → Profile`, transient piecewise-linear, attached to the same BC fields.

**S7. Boundary conditions.**
| Boundary | Thermal BC | Values |
|---|---|---|
| `out_wall_*`, `out_roof` | **Mixed** (h_ext + T_ext + T_sky) | h_ext = `wind_conv_profile` UDF (5.7+3.8V), T_ext = `ambient_temp_profile`, T_sky = `sky_temp_profile`, ε = 0.9 |
| `out_soil_sides` | Heat flux ≈ 0 (adiabatic) | margin ≥ 3 m already decouples |
| `out_soil_bottom` | Temperature | `ground_deep_profile` (e.g. 281 K) |
| Wind inlet (west face) | Velocity-inlet, magnitude + direction | magnitude = `wind_speed_profile`, T = `ambient_temp_profile` |
| Opposite/top faces | Pressure-outlet / symmetry | T backflow = `ambient_temp_profile` |
| Glazing walls (DO two-band) | Semi-transparent | τ_vis 0.6/0.12 (clear/tinted), τ_IR 0, ε 0.84; direct/diffuse absorb per band |
| `trombe_absorber` | Opaque, participates | α_solar 0.94, ε 0.9 |
| Vents | Interior openings | porous-jump (K=10, Δp via Darcy–Forchheimer) or clear fluid |
| `door_n` | Wall U=1.5 + infiltration | leakage via §6.5 ACH source, or resolve 5 mm perimeter gap (Δ₁ mesh 2 mm) |

**Path A vs Path B (choose per fenestration — never both):**
- **Path A (dynamic SHGC, default):** glazing is *opaque to solar ray tracing* (set solar transmissivity
  0) and transmitted energy enters through `win*_gain` volumetric source = UDF `win*_solar_gain`,
  which applies P(t) = SHGC(t)·A_w·G_incident,orientation(t) pre-computed from EPW by `make_inputs.py`.
  This is what makes electrochromic/scheduled SHGC exact and trims DO cost.
- **Path B (fixed SHGC):** semi-transparent two-band DO glazing; solar ray tracing carries the beam
  into the room (correct sun-patch geometry, validates Path A).

**S8. SHGC & Trombe controls.** Hook `win1_solar_gain`/`win2_solar_gain` on the gain zones.
Trombe dampers: schedule via `vent_porosity_profile` UDF (porous-jump permeability × open-fraction)
— e.g. vents open 09:00–17:00 winter, closed at night.

**S9. Solver numerics.** Pressure–velocity coupling: **Coupled** (or PISO with skewness correction);
spatial: gradient least-squares cell-based, pressure **PRESTO!** (buoyancy), momentum/energy/k/ε
**Second-Order Upwind**; DO second order. Temporal: **Second-Order Implicit**.
URFs (coupled): 0.3–0.5 initial → relax to 0.7. Flows with vents may need PISO fallback
(pressure 0.3, momentum 0.7) if coupled stalls.

**S10. Initialization & spin-up.** Standard-init whole domain at T_amb(t₀) = 02:00 hr value;
patch `soil` to T_ground,deep, walls to linear gradient. Then run **5–7 full diurnal cycles**
(Δt = 120 s, ≤ 25 iter/step) with periodic climate data before the recording day. Periodicity is
declared when peak-interior-T drifts < 0.5 % between consecutive days.

**S11. Monitors & autosave.** Report definitions (§12) plotted each step; autosave case+data every
2 h simulation time; convergence: energy residual < 1e-7, continuity/k/ε < 1e-4 *and* flatlined
monitors (residual-only is insufficient with DO + buoyancy).

---

## 7. Climate-Data Ingestion (UDF contract)

`make_inputs.py --epw <file.epw> --day 355` regenerates `udf/climate_data.h` with hourly arrays:

```
TA[]      dry-bulb ambient        [°C]   → ambient_temp_profile  (face temp, K)
WS[]      wind speed @10 m        [m/s]  → wind_speed_profile (power-law α=0.16) + wind_conv_profile
DNI[]     direct normal irradiance[W/m²] → DEFINE_SOLAR_INTENSITY direct_from_epw
DHI[]     diffuse horizontal      [W/m²] → DEFINE_SOLAR_INTENSITY diffuse_from_epw
PWIN1[].. transmitted gain per window incl. SHGC(t) & orientation [W] → win*_solar_gain (source)
TG        deep-ground temp        [°C]   → ground_deep_profile
T_START_HR, N_DATA_HOURS, zone constants (areas, volumes)
```

Key macros in `climate_shelter_udf.c` (full, compilable file in `udf/`):

```c
DEFINE_PROFILE(ambient_temp_profile, thread, nv)     /* hourly → K, linear interp, periodic */
{
    face_t f; real h = sim_hour();
    begin_f_loop(f, thread)
    {
        real x[ND_ND]; F_CENTROID(x, f, thread);
        F_PROFILE(f, thread, nv) = C2K(interpolate(TA, h));
    }
    end_f_loop(f, thread)
}

DEFINE_SOLAR_INTENSITY(direct_from_epw, sun_x, sun_y, sun_z, hour, minute)
{   return MAX(0.0, interpolate(DNI, hour + minute/60.0)); }   /* W/m², solar time from solver */

DEFINE_SOURCE(win1_solar_gain, cell, thread, dS, eqn)          /* dynamic-SHGC path (Path A) */
{
    real P = interpolate(PWIN1, sim_hour());                   /* W, SHGC(t)-weighted */
    dS[eqn] = 0.0;
    return P / (real)N_CELLS_WIN1 / C_VOLUME(cell, thread);    /* W/m³, uniform */
}
```
`DEFINE_SOLAR_INTENSITY` is the documented Fluent hook for solar-load irradiation (3-D only);
sun *position* still comes from the Solar Calculator (date/time/latitude), so shadows stay physical.

Wind profile: `V(z) = V₁₀ (z/10)^0.16` (open terrain). Sky: Swinbank on TA.
Logger: `DEFINE_EXECUTE_AT_END(hourly_audit)` appends `results/audit.csv`
(`t[h], T_in_avg, T_in_min, T_in_max`) once per simulated hour — this CSV feeds time-lag/decrement
post-processing in the automation script.

---

## 8. Openings, Fenestration & Trombe Modeling Notes
- **Trombe vents**: resolving blades costs cells; porous-jump with inertial resistance
  `K = Δp/(ρ v²/2)/t` calibrated from damper open-area fraction is the standard economy. Buoyant
  circulation is driven by Δρ·g·H_stack between gap and room — keep ≥5 cells across the gap width.
- **Shading louvers**: geometry on (resolved) for south-summer study; equivalent-τ porous wall off.
- **Night shutters/curtains** (big lever in Ladakh): emulate as transient glazing τ schedule —
  i.e., switch Path-A PWIN tables to night values (τ_night ≈ 0.05).
- **Infiltration**: Q_inf = ρ V ACH c_p ΔT/3600; implement as negative energy source on
  `air_interior` via `DEFINE_SOURCE` with ACH from user data (or resolve door gap at 2 mm mesh).

---

## 9. Thermal-Loss Accounting (what to report, exactly)

Every quantity below is a Fluent **surface integral of Total Surface Heat Flux** (W), logged hourly
(report definitions; the automation script writes them to `results/losses.csv`):

1. **Opaque conduction**: `in_wall_s/e/w/n + in_roof` interiors (conduction out > 0 = loss).
2. **Fenestration**: `in_glass_win1/2` total flux split into its **convective** and **radiation**
   parts (Fluent reports both on wall flux panel) → window U-loss vs solar-gain balance.
3. **Floor/ground**: `in_floor` (slab↔room) and `out_soil_bottom` (deep sink); their difference is
   soil thermal-storage — the "free battery" the PS wants quantified.
4. **Thermal bridges**: `stud_steel` wetted flux − (U_clear·A_stud_share·ΔT) → ψ-matrix per §2.6;
   report % of total transmission loss (expect 20–40 %).
5. **Ventilation/infiltration**: mass-flow-weighted enthalpy flux through vents/door gap:
   Σ ṁ c_p (T_out − T_in) over vent openings.
6. **Radiative exterior penalty**: Σ (radiation flux) on `out_wall_*` to T_sky — motivates
   low-ε skin/roof选项 in variants.
7. **Solar input ledger**: sunload report (TUI `/report/solar-data`) + `in_glass_*` absorbed +
   `win*_gain` sources = "thermal energy generated from solar radiation" (PS Task 2).

**Closure check every 24 h:** |Σ losses − Σ gains − ΔU_storage| / Σ gains < 2 %.

---

## 10. Comfort Evaluation & Post-Processing
1. **Air temperature field**: volume stats on `air_interior` (mean/min/max) + planes y = 0.1/1.1/1.7 m;
   report stratification ΔT(1.7 − 0.1) < 3 K target.
2. **Operative temperature** T_op = (T_a + T_MRT)/2 at occupant points (bed/cot & desk positions).
3. **MRT plane**: custom field function `((incident-radiation)/(4*5.670e-8))^0.25` on y = 1.1 m;
   export incident radiation + surface temps → exact Σ F_i T_i⁴ in the post script.
4. **MRT_solar**: add ΔT_sw (§2.8) using DO shortwave at occupant position; report both day & night.
5. **PMV/PPD**: Fluent Thermal Comfort post-processing (needs DO) or offline ISO-7730 calc from
   T_a, T_MRT, V; report comfort-hours on the recording day (target ≥ 90 % light-activity hours
   ≥ 18 °C with 1.2 clo).
6. **Time lag / decrement factor** per envelope element: extract inner-surface T and sol-air series
   from `audit.csv` + flux reports; φ and f per §2.7 formulas; done in `pyfluent_run.py` post stage.
7. **Sun patches**: isosurfaces of incident radiation > 200 W/m² to visualize beam entry (Path B).

---

## 11. Run Plan & Compute
- Δt = 120 s (60 s during 08:00–16:00 transients if stability demands), 2nd-order implicit.
- Demo mesh (~1.5 M cells, 16 cores): ≈ 4–8 h per 24-h cycle. Production (6 M, 64 cores): ~1–2 days
  for 7-day spin-up + 2 recording days. Use checkpoint-restart (autosave + `Read → Data`).
- DOE mode: short 48-h runs (1 spin + 1 record) per variant on the demo mesh ranks variants
  correctly; confirm the top-2 on the production mesh.

## 12. Outputs → PS "Tasks 1–3"
| PS task | Artifact |
|---|---|
| 1. Inside-temperature prediction | `audit.csv` T_in series; contour movie; KPI table (T_min, T_mean, comfort-hours) |
| 2. Solar energy generated | Daily kWh ledger per §9.7 (direct+diffuse→interior), per variant |
| 3. Heat-flow details vs ΔT | Hourly loss-split table (§9), Q vs (T_in − T_amb) scatter with h_eff fit |

`pyfluent_run.py` DOE ranking CSV columns: `variant, U_wall, T_in_min, T_in_mean, comfort_hours,
Q_solar_in_kWh_day, Q_loss_kWh_day, time_lag_h, decrement_factor` → the "most efficient combination
of materials, shape, size" the PS expects, with CFD evidence.

## 13. Validation & QA Checklist
1. U-value: steady two-zone checkerboard run (no solar, ΔT = 30 K fixed) vs hand §2.5 (±5 %).
2. y+ histograms ≤ 1 (walls), ≤ 5 max. 3. 24-h energy closure < 2 %.
4. Grid: GCI-style re-run with 1.5× coarser fluid size — peak T shift < 0.3 K.
5. Time step: 60 s vs 120 s — T_in_min shift < 0.2 K.
6. DO discretization 4×8 vs 8×16 — interior flux < 3 % shift.
7. Solar sanity: noon vertical-south incidence on 21 Dec ≈ cos(34.16°+23.45°-tilt geometry) per §2.4;
   compare ray-traced facade flux vs hand DNI·cosθ + DHI/2.
8. Known pitfalls: double-counting solar (Path A + semi-transparent glazing); operating density not
   set (spurious −1 to +2 K drift); Boussinesq beyond ΔT ≈ 15 K; shell-conduction walls not marked
   "participates in radiation" (zero MRT contribution); forgetting "Time-steps-per-solar-update = 1"
   (sun frozen at t₀).

## 14. References
Fluent User's Guide §Solar Load Model, §Discrete Ordinates, §Enhanced Wall Treatment; Fluent UDF
Manual §DEFINE_SOLAR_INTENSITY, §DEFINE_PROFILE/SOURCE; ISO 6946, ISO 10211, ISO 14683, ISO 13370,
ISO 7730, ASHRAE Fundamentals Ch. 18 (sol-air, Swinbank), ASHRAE 55 Annex solar-MRT (Arens et al.
2015); EnergyPlus EPW format documentation; IMD/NIWE Ladakh climate normals & MNRE solar handbooks
(GHI 1900–2100 kWh/m²·yr, sunshine 7.9 h/day — PS figures).
