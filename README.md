# SIH26051 — Area-Specific Passive Shelter Thermal Model (ANSYS Fluent)

Software package for **SIH26051** (DRDO, Smart India Hackathon 2026): a transient
conjugate-heat-transfer model that predicts inside temperature, solar gains and heat
losses of an area-specific passive shelter for high-altitude cold regions (Ladakh
baseline), driven by real EPW climate data.

## Quick start (3 commands)

```bash
pip install gmsh numpy                      # geometry/mesh toolchain
python tools/make_inputs.py --demo          # EPW -> UDF header + Fluent profile
python geometry/make_shelter_mesh.py        # -> geometry/shelter.cgns (~0.7 M cells)
```

Then in Fluent (or `automation/pyfluent_run.py`):

```
File > Read > Mesh...        geometry/shelter.cgns
User-Defined > Functions > Compiled...   udf/libudf   (build from climate_shelter_udf.c)
Run the transient case per docs/SIH26051_Fluent_Setup_Guide.md  Sec. 6 (S1–S11)
```

`python automation/pyfluent_run.py --run-mode demo|full|doe` automates steps 2–12
(models, materials, UDF hook, solar calculator, BCs, run, KPIs, DOE ranking).

## Package map

| Path | What it is |
|---|---|
| `docs/SIH26051_Fluent_Setup_Guide.md` | **Start here** — full engineering setup guide: equations, geometry, materials, meshing numbers, solver settings, BCs, loss accounting, comfort/MRT, validation checklist |
| `udf/climate_shelter_udf.c` | UDF suite: hourly EPW ingestion (`DEFINE_PROFILE`), `DEFINE_SOLAR_INTENSITY` (DNI/DHI), dynamic-SHGC window sources, infiltration, sky temperature, hourly audit logger. Syntax-checked for serial + parallel |
| `udf/climate_data.h` | Generated hourly climate arrays (demo: Leh, 21 Dec). Regenerate from your EPW |
| `tools/make_inputs.py` | EPW → `climate_data.h` + Fluent transient profile + per-window transmitted-gain tables (orientation-aware, dynamic SHGC, snow albedo) |
| `geometry/make_shelter_mesh.py` | Parametric gmsh model: 5-layer walls, studs, windows + gain sheets, Trombe + vents, exterior air, 2 m soil → conformal **CGNS for Fluent**. Tested: 712k cells, all zones named |
| `geometry/shelter_spaceclaim.py` | SpaceClaim/IronPython alternative (recorded-API sketch — verify against your release) |
| `automation/pyfluent_run.py` | PyFluent end-to-end: launch → setup → UDF → run → KPIs → DOE wall-variant ranking (`results/doe_ranking.csv`) |
| `profiles/amb_profile.pro` | Fluent transient profile (UDF-free fallback) |
| `profiles/gain_table.csv` | Hourly per-window transmitted solar power (W) with SHGC state |
| `geometry/shelter.cgns`, `geometry/shelter.msh` | Generated demo meshes (0.7 M cells) |

## How it answers the problem statement

1. **Inside-temperature prediction** — transient CHT run; `results/audit.csv` hourly
   interior mean/min/max from the UDF logger.
2. **Solar energy generated** — `DEFINE_SOLAR_INTENSITY` feeds measured DNI/DHI into
   Solar Ray Tracing (physical shadows); transmitted gains with **dynamic SHGC** per
   window from `make_inputs.py` (ledger in `gain_table.csv`).
3. **Heat-flow details vs ΔT** — per-surface heat-flux report definitions (envelope,
   glazing, slab↔soil, studs→thermal bridges) + hourly energy-balance closure check.
4. **Comparative analysis** — `--run-mode doe` ranks wall stacks (EPS / PUF+concrete /
   PCM / VIP / adobe) by comfort hours, minimum inside temperature, time lag and
   decrement factor → the "most efficient combination of materials, shape, size".

## Field data you should substitute (placeholders marked in code)

- EPW for Leh (ISHRAE/MetDotClim) — `make_inputs.py --epw ... --day 355`
- Deep-ground temperature at 2 m (calibrate from IMD/site data; demo: +8 °C)
- Infiltration ACH from the actual shelter door/seal spec (demo: 0.5)
- Damper loss coefficient for Trombe vents (demo placeholder 200 1/m)

## Version notes

- Fluent 2023R2+ assumed; every version-fragile PyFluent call degrades to a printed
  MANUAL STEPS checklist instead of crashing.
- gmsh ≥ 4.11 with CGNS support (pip wheel). If CGNS export fails on your build,
  Fluent also reads `.msh` via File > Import > Mesh (or convert with meshio).
- Deadline triage: docs guide §0 has a 3-day fast-track plan (submission 30 Sep 2026).
