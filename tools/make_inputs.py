#!/usr/bin/env python3
"""
make_inputs.py -- SIH26051 shelter model: climate-data pipeline (EPW -> Fluent/UDF inputs).

Reads an EnergyPlus Weather (EPW) file, extracts a design day (or multi-day average),
models orientation-resolved solar incidence + dynamic-SHGC window gains, and emits:

  1. climate_data.h            -- C header consumed by udf/climate_shelter_udf.c
  2. amb_profile.pro           -- Fluent transient profile (UDF-free fallback for judges' laptops)
  3. gain_table.csv            -- per-window transmitted power, human-checkable

Usage
-----
  # Real data (recommended): ISHRAE/MetDotClim EPW for Leh (e.g. IN_Leh.420xxx_ISHRAE.epw)
  python make_inputs.py --epw path/to/leh.epw --day 355 --out ../udf

  # No EPW at hand? Built-in Leh 21-Dec design day:
  python make_inputs.py --demo

Window definition (edit WINDOWS below or pass --windows):
  WINDOWS = [dict(name="win1", azimuth_deg=180, tilt_deg=90, area_m2=1.2, gain_zone_volume_m3=0.6), ...]
  azimuth: 180 = south, 90 = east, 270 = west, 0 = north.  tilt: 90 = vertical.

Solar model: NOAA declination/hour-angle, isotropic-sky diffuse + ground-reflected:
  I_vert = DNI*cos(theta) + DHI*0.5 + GHI*rho_g*0.5
Dynamic SHGC (electrochromic / smart glass): SHGC = SHGC_TINTED if I_vert > THRESHOLD else SHGC_CLEAR.
"""
import argparse
import csv
import math
import os
import sys

# ---------------------------------------------------------------- user parameters ---------
WINDOWS = [
    dict(name="win1", azimuth_deg=180.0, tilt_deg=90.0, area_m2=1.2),   # south facade
    dict(name="win2", azimuth_deg=90.0, tilt_deg=90.0, area_m2=1.2),    # east facade
]
SHGC_CLEAR, SHGC_TINTED, SHGC_THRESHOLD_W = 0.60, 0.25, 350.0
GROUND_ALBEDO = 0.60          # snow-covered high-altitude ground (0.2 bare, up to 0.8 fresh snow)
TG_DEEP_C = 8.0               # deep-ground temperature at 2 m (calibrate with IMD/site data)
T_START_HR = 2.0              # simulation starts 02:00 (coldest hour) for clean spin-up
WIND_ALPHA = 0.16             # open-terrain power-law exponent
ACH = 0.5                     # infiltration air changes per hour (door/seal spec)
V_ROOM_M3 = 55.0              # interior air volume
V_SOURCE_ZONE_M3 = 0.6        # volume of each thin 'win*_gain' fluid sheet
VENT_K_OPEN = 200.0           # Trombe vent inertial resistance, 1/m (calibrate from damper)
ZONE_ID_AIR_INTERIOR, ZONE_ID_WIN1_GAIN, ZONE_ID_WIN2_GAIN = 5, 6, 7  # <-- FIX AFTER MESHING

# ---------------------------------------------------------------- solar geometry ----------
def solar_position(doy, hour_solar, lat_deg):
    """Returns (altitude_rad, azimuth_compass_rad) for solar time hour_solar (0-24)."""
    phi = math.radians(lat_deg)
    decl = math.radians(23.45) * math.sin(math.radians(360.0 * (284 + doy) / 365.0))
    omega = math.radians(15.0 * (hour_solar - 12.0))
    sin_alt = math.sin(phi) * math.sin(decl) + math.cos(phi) * math.cos(decl) * math.cos(omega)
    alt = math.asin(max(-1.0, min(1.0, sin_alt)))
    # azimuth measured from SOUTH (NOAA form): psi = atan2(sin w, cos w sin phi - tan d cos phi)
    psi = math.atan2(math.sin(omega), math.cos(omega) * math.sin(phi) - math.tan(decl) * math.cos(phi))
    az_compass = math.pi + psi          # 0=N, 90=E, 180=S, 270=W  (psi<0 => morning/east)
    return alt, az_compass

def cos_incidence(alt, az_compass, surf_az_deg, surf_tilt_deg):
    dpsi = az_compass - math.radians(surf_az_deg)
    cos_beam = math.cos(alt) * math.cos(dpsi)                       # vertical-surface beam term
    if abs(surf_tilt_deg - 90.0) > 1e-3:                            # non-vertical: generic form
        ct = math.sin(alt) * math.cos(math.radians(surf_tilt_deg)) + \
             math.cos(alt) * math.cos(dpsi) * math.sin(math.radians(surf_tilt_deg))
        return max(0.0, ct)
    return max(0.0, cos_beam)

# ---------------------------------------------------------------- demo climate ------------
def demo_leh_design_day():
    """Leh (34.16N, 77.58E), 21 December, ~300 cloud-free days: cold, clear, calm morning."""
    ta  = [-13.5, -14.0, -14.2, -14.3, -14.0, -13.2, -12.0, -10.5, -9.0, -7.2, -5.5, -4.0,
           -2.8, -2.0, -2.2, -3.0, -4.5, -6.5, -8.5, -10.0, -11.2, -12.0, -12.6, -13.1]
    ws  = [1.2, 1.0, 1.0, 1.0, 1.1, 1.3, 1.6, 2.0, 2.4, 2.8, 3.2, 3.6,
           3.8, 3.6, 3.2, 3.0, 3.1, 3.4, 3.0, 2.5, 2.0, 1.6, 1.4, 1.3]
    dni = [0, 0, 0, 0, 0, 0, 40, 300, 520, 650, 720, 760, 780, 770, 720, 610, 380, 90, 0, 0, 0, 0, 0, 0]
    dhi = [0, 0, 0, 0, 0, 0, 25, 55, 75, 90, 100, 105, 110, 108, 100, 88, 65, 30, 0, 0, 0, 0, 0, 0]
    return ta, ws, dni, dhi

# ---------------------------------------------------------------- EPW reader --------------
def read_epw_day(epw_path, doy, avg_days=1):
    """Returns (ta, ws, dni, dhi) hourly arrays for the day centered on DOY (avg of avg_days)."""
    rows = []
    with open(epw_path, newline="") as fh:
        for _ in range(8):
            next(fh)
        rd = csv.reader(fh)
        for r in rd:
            rows.append(r)
    def doy_of(r):
        return int(r[1]) * 100 + int(r[2])
    from datetime import date, timedelta
    base = date(2001, 1, 1) + timedelta(days=doy - 1)
    idxs = []
    for off in range(-(avg_days // 2), avg_days // 2 + 1):
        d = base + timedelta(days=off)
        key = d.month * 100 + d.day
        idxs.extend(i for i, r in enumerate(rows) if doy_of(r) == key)
    if not idxs:
        sys.exit(f"EPW: no rows for day-of-year {doy}")
    def f(r, col):
        try:
            return float(r[col])
        except (ValueError, IndexError):
            return 0.0
    n = len(idxs) // 24
    ta  = [sum(f(rows[i], 6)  for i in idxs[h*24:(h+1)*24]) / n for h in range(24)]
    ws  = [sum(f(rows[i], 21) for i in idxs[h*24:(h+1)*24]) / n for h in range(24)]
    dni = [sum(f(rows[i], 14) for i in idxs[h*24:(h+1)*24]) / n for h in range(24)]
    dhi = [sum(f(rows[i], 15) for i in idxs[h*24:(h+1)*24]) / n for h in range(24)]
    return ta, ws, dni, dhi

# ---------------------------------------------------------------- physics -----------------
def window_gains(doy, lat, ta, dni, dhi):
    """Per-window transmitted power [W] with dynamic SHGC; returns dict name->(list24, states24)."""
    out = {}
    for w in WINDOWS:
        p, states = [], []
        for h in range(24):
            alt, az = solar_position(doy, h + 0.5, lat)
            ghi = dni[h] * max(0.0, math.sin(alt)) + dhi[h]
            iv = dni[h] * cos_incidence(alt, az, w["azimuth_deg"], w["tilt_deg"]) \
                 + dhi[h] * 0.5 * (1.0 + math.cos(math.radians(w["tilt_deg"]))) \
                 + ghi * GROUND_ALBEDO * 0.5 * (1.0 - math.cos(math.radians(w["tilt_deg"])))
            tinted = iv > SHGC_THRESHOLD_W
            states.append("tinted" if tinted else ("clear" if iv > 1.0 else "night"))
            p.append((SHGC_TINTED if tinted else SHGC_CLEAR) * iv * w["area_m2"])
        out[w["name"]] = (p, states)
    return out

# ---------------------------------------------------------------- writers -----------------
HDR_TMPL = """/* climate_data.h -- GENERATED by tools/make_inputs.py ({src}). DO NOT EDIT BY HAND.
 * Hourly climate inputs for udf/climate_shelter_udf.c  ({desc})
 * Regenerate with your EPW:  python tools/make_inputs.py --epw <file.epw> --day <DOY>
 */
#ifndef CLIMATE_DATA_H
#define CLIMATE_DATA_H

#define N_DATA_HOURS {n}
#define TG   {tg}      /* deep-ground temperature at 2 m  [C]  */
#define T_START_HR {tstart}  /* simulation start hour (2.0 = coldest hour) */
#define WIND_ALPHA {walpha}    /* wind power-law exponent */
#define ACH  {ach}      /* infiltration air changes per hour */
#define V_ROOM {vroom}     /* interior air volume [m3] */
#define V_SOURCE_ZONE {vsrc}   /* win*_gain fluid-sheet volume [m3] */
#define VENT_K_OPEN {ventk}  /* Trombe vent inertial resistance when open [1/m] */
#define AIR_RHO {airrho}   /* air density at T_ref [kg/m3] */
#define AIR_CP  {aircp}    /* air specific heat [J/kgK] */

/* cell-zone IDs assigned by Fluent (see TUI: /mesh/zones or GUI Outlines->Cell Zones) */
#define ZONE_ID_AIR_INTERIOR {z_air}
#define ZONE_ID_WIN1_GAIN   {z_w1}
#define ZONE_ID_WIN2_GAIN   {z_w2}
#define N_CELLS_WIN1 {nc1}   /* update from hourly_audit first-call printout */
#define N_CELLS_WIN2 {nc2}

/* hour-locked climate tables [N_DATA_HOURS] */
static const real TA[{n}]   = {{ {ta} }};
static const real WS[{n}]   = {{ {ws} }};
static const real DNI[{n}]  = {{ {dni} }};
static const real DHI[{n}]  = {{ {dhi} }};
/* transmitted solar power per window incl. dynamic SHGC & orientation [W] */
static const real PWIN1[{n}] = {{ {p1} }};
static const real PWIN2[{n}] = {{ {p2} }};

#endif /* CLIMATE_DATA_H */
"""

def write_climate_header(path, ta, ws, dni, dhi, p1, p2, src, desc):
    def fmt(vals, nd=2):
        return ", ".join(f"{v:.{nd}f}" for v in vals)
    with open(path, "w") as fh:
        fh.write(HDR_TMPL.format(
            src=src, desc=desc, n=len(ta), tg=TG_DEEP_C, tstart=T_START_HR, walpha=WIND_ALPHA,
            ach=ACH, vroom=V_ROOM_M3, vsrc=V_SOURCE_ZONE_M3, ventk=VENT_K_OPEN,
            airrho=1.34, aircp=1006.0,
            z_air=ZONE_ID_AIR_INTERIOR, z_w1=ZONE_ID_WIN1_GAIN, z_w2=ZONE_ID_WIN2_GAIN,
            nc1=1000, nc2=1000, ta=fmt(ta, 2), ws=fmt(ws, 2), dni=fmt(dni, 1), dhi=fmt(dhi, 1),
            p1=fmt(p1, 1), p2=fmt(p2, 1)))
    print(f"wrote {path}")

def write_fluent_profile(path, ta_k, ws):
    """Fluent transient profile: ((name transient n periodic) time n-vals name n-vals ...)."""
    n = len(ta_k)
    with open(path, "w") as fh:
        fh.write(f"((amb_climate transient {n} periodic)\n")
        fh.write("time\n")
        fh.write("\n".join(f"{h*3600.0:.1f}" for h in range(n)) + "\n")
        fh.write("ambient_temp_k\n")
        fh.write("\n".join(f"{v:.2f}" for v in ta_k) + "\n")
        fh.write("wind_speed_ms\n")
        fh.write("\n".join(f"{v:.2f}" for v in ws) + "\n")
        fh.write(")\n")
    print(f"wrote {path}")

def write_gain_csv(path, gains):
    p1, s1 = gains[WINDOWS[0]["name"]]
    p2, s2 = (gains[WINDOWS[1]["name"]] if len(WINDOWS) > 1 else ([0.0] * 24, s1))
    with open(path, "w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(["hour", "P_win1_W", "SHGC_win1", "P_win2_W", "SHGC_win2"])
        for h in range(len(p1)):
            wr.writerow([h, f"{p1[h]:.1f}", s1[h], f"{p2[h]:.1f}", s2[h]])
    print(f"wrote {path}")

# ---------------------------------------------------------------- main --------------------
def main():
    ap = argparse.ArgumentParser(description="EPW -> Fluent/UDF climate inputs (SIH26051)")
    ap.add_argument("--epw", help="path to .epw file")
    ap.add_argument("--day", type=int, default=355, help="design day-of-year (355 = 21 Dec)")
    ap.add_argument("--avg-days", type=int, default=1, help="average ±N/2 days around --day")
    ap.add_argument("--lat", type=float, default=34.16, help="latitude (Leh = 34.16)")
    ap.add_argument("--out", default=os.path.join(os.path.dirname(__file__), "..", "udf"))
    ap.add_argument("--profiles-out", default=os.path.join(os.path.dirname(__file__), "..", "profiles"))
    ap.add_argument("--demo", action="store_true", help="use built-in Leh 21-Dec design day")
    a = ap.parse_args()

    if a.demo or not a.epw:
        ta, ws, dni, dhi = demo_leh_design_day()
        src, desc = "built-in demo", "Leh design day, 21 Dec (demo arrays)"
    else:
        ta, ws, dni, dhi = read_epw_day(a.epw, a.day, a.avg_days)
        src, desc = os.path.basename(a.epw), f"EPW day {a.day} (avg {a.avg_days} d), lat {a.lat}"

    gains = window_gains(a.day if not a.demo else 355, a.lat, ta, dni, dhi)
    p1 = gains[WINDOWS[0]["name"]][0]
    p2 = gains[WINDOWS[1]["name"]][0] if len(WINDOWS) > 1 else [0.0] * 24

    out = os.path.abspath(a.out)
    prof = os.path.abspath(a.profiles_out)
    os.makedirs(out, exist_ok=True)
    os.makedirs(prof, exist_ok=True)

    write_climate_header(os.path.join(out, "climate_data.h"), ta, ws, dni, dhi, p1, p2, src, desc)
    write_fluent_profile(os.path.join(prof, "amb_profile.pro"), [t + 273.15 for t in ta], ws)
    write_gain_csv(os.path.join(prof, "gain_table.csv"), gains)

    peak = max(p1)
    print(f"\nKPI: peak transmitted gain win1 = {peak:.0f} W | daily energy = "
          f"{sum(p1)/1000.0:.2f} kWh | T range {min(ta):.1f}..{max(ta):.1f} C | "
          f"peak DNI {max(dni):.0f} W/m2")
    print("REMINDER: after meshing, paste Fluent cell-zone IDs and N_CELLS_* counts into "
          "climate_data.h (hourly_audit prints them on first call).")

if __name__ == "__main__":
    main()
