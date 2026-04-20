#!/usr/bin/env python3
"""
createGeo_height.py

Parametri supportati in &geoSettings:

Heights (mm):
- primary_height, secondary_height
  oppure in alternativa:
- shell_height            -> usato per entrambi (primary e secondary)

Starts (mm):
- primary_start, secondary_start   (default 0.0)

Thickness law:
- thickness_intercept, thickness_slope
  (accetta anche thickQ / ThickM per compatibilità)

Altri:
- nShells
- focal_length                    (m)
- maxDiam                         (mm)
- angular_shell_separation_deg
- minThick                        (mm)
- WallDensity                     (kg/m^3)
- profileType                     (es. ph)

Uso:
    python createGeo_height.py geoSettings.txt
    python createGeo_height.py geoSettings.txt -d diameters.txt -o shellStruct.dat
"""

#TODO: verify varibable names (i.e. if include units in the name)
#TODO: verify what is the use case for -d: one might want to provide diameters at IP, or at entrance pupil, or maybe it is never useful.
#TODO: verify handling of IP gap and partial illumination considering hyperbola and parabola separately. Add test cases.
#TODO: document how the geometry file is used in EA.py (how the relevant column is selected).

from __future__ import annotations

import argparse
import math
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

PI = math.pi


def _strip_fortran_comments(line: str) -> str:
    return line.split("!")[0]


def replace_fortran_exponent(s: str) -> str:
    return re.sub(r"([0-9])([dD])([+\-]?[0-9]+)", r"\1e\3", s)


def _parse_scalar(val: str) -> Any:
    v = val.strip()
    if not v:
        return v

    vl = v.lower()
    if vl in (".true.", "true"):
        return True
    if vl in (".false.", "false"):
        return False

    if (v.startswith("'") and v.endswith("'")) or (v.startswith('"') and v.endswith('"')):
        return v[1:-1]

    if "d" in vl:
        v = replace_fortran_exponent(v)

    try:
        if any(c in v for c in ".eE"):
            return float(v)
        return int(v)
    except ValueError:
        return v


def _split_top_level_commas(s: str) -> List[str]:
    """Split by commas, ignoring commas inside quotes."""
    out, buf = [], []
    in_sq = in_dq = False
    for ch in s:
        if ch == "'" and not in_dq:
            in_sq = not in_sq
        elif ch == '"' and not in_sq:
            in_dq = not in_dq
        if ch == "," and not in_sq and not in_dq:
            token = "".join(buf).strip()
            if token:
                out.append(token)
            buf = []
        else:
            buf.append(ch)
    token = "".join(buf).strip()
    if token:
        out.append(token)
    return out


def read_geo_settings(path: Path, group: str = "geoSettings") -> Dict[str, Any]:
    """
    Legge un file stile namelist Fortran:
        &geoSettings
        key = value
        ...
        /
    Supporta assegnamenti separati per riga e/o separati da virgola.
    """
    txt = path.read_text(encoding="utf-8", errors="ignore").splitlines()
    lines = [_strip_fortran_comments(l).strip() for l in txt]

    start = None
    gpat = re.compile(rf"^\s*&\s*{re.escape(group)}\b", re.IGNORECASE)
    for i, l in enumerate(lines):
        if gpat.search(l):
            start = i
            break
    if start is None:
        raise ValueError(f"Namelist group '&{group}' not found in {path}")

    out: Dict[str, Any] = {}
    for l in lines[start + 1 :]:
        if re.match(r"^\s*/\s*$", l) or re.match(r"^\s*&\s*end\b", l, re.IGNORECASE):
            break
        if not l:
            continue
        for token in _split_top_level_commas(l):
            if "=" not in token:
                continue
            k, v = token.split("=", 1)
            out[k.strip()] = _parse_scalar(v)
    return out


def _get_ci(d: Dict[str, Any], key: str, default: Any = None) -> Any:
    kl = key.lower()
    for k, v in d.items():
        if k.lower() == kl:
            return v
    return default


def poly_coeff(profile_type: str, radius_ip: float, f_length_mm: float) -> Tuple[Tuple[float, float, float, float, float],
                                                                                Tuple[float, float, float, float, float]]:
    if not profile_type or len(profile_type) < 2:
        raise ValueError("profileType must be a 2-character string (e.g., 'cc', 'cp', 'ph').")

    p1 = profile_type[0].lower()
    p2 = profile_type[1].lower()

    angle_one = math.atan2(radius_ip, f_length_mm) / 4.0
    a1 = -2.0 * math.tan(angle_one)
    a2 = a3 = a4 = a5 = 0.0

    if p1 == "c":
        a2 = math.tan(angle_one) ** 2
    elif p1 == "p":
        a2 = 0.0
    elif p1 == "u":
        a2 = 0.0002005
        a3 = 0.0003318
        a4 = 0.0002711
        a5 = 0.00009652
    elif p1 == "h":
        raise NotImplementedError("Hyperbolic primary surface (profileType[0]='h') is not implemented.")
    else:
        raise ValueError(f"Primary profile type '{profile_type[0]}' not recognized (valid: c, h, p, u).")

    angle_two = 3.0 * angle_one
    b1 = -2.0 * math.tan(angle_two)
    b2 = b3 = b4 = b5 = 0.0

    if p2 == "c":
        b2 = math.tan(angle_two) ** 2
    elif p2 == "p":
        b2 = 0.0
    elif p2 == "u":
        b2 = 0.001
        b3 = 0.0004591
        b4 = -0.0004248
        b5 = 0.0001651
    elif p2 == "h":
        b2 = 2.0 * (radius_ip * math.tan(angle_two)) / (f_length_mm + radius_ip / math.tan(2.0 * angle_one))
    else:
        raise ValueError(f"Secondary profile type '{profile_type[1]}' not recognized (valid: c, h, p, u).")

    return (a1, a2, a3, a4, a5), (b1, b2, b3, b4, b5)


def poly_invert(z: float, r1: float, c1: float, c2: float) -> float:
    delta = (c1 * z) ** 2 - 4.0 * c2 * (z ** 2) + 4.0 * (r1 ** 2)
    if delta < 0:
        delta = 0.0
    return (-c1 * z + math.sqrt(delta)) / 2.0


def poly_profile(z: float, r0: float, c: Tuple[float, float, float, float, float]) -> float:
    c1, c2, c3, c4, c5 = c
    t = z / r0
    inside = 1.0 + c1 * t + c2 * (t ** 2) + c3 * (t ** 3) + c4 * (t ** 4) + c5 * (t ** 5)
    if inside < 0:
        inside = 0.0
    return r0 * math.sqrt(inside)


@dataclass
class GeoSettings:
    n_shells: int
    max_diam_mm: float
    focal_length_m: float
    angular_shell_separation_deg: float

    primary_height_mm: float
    secondary_height_mm: float
    primary_start_mm: float = 0.0
    secondary_start_mm: float = 0.0

    thickness_intercept_mm: float = 0.0
    thickness_slope_per_mm: float = 0.0

    min_thick_mm: float = 0.0
    wall_density_kg_m3: float = 0.0
    profile_type: str = "cc"

    @staticmethod
    def from_dict(d: Dict[str, Any]) -> "GeoSettings":
        shell_h = _get_ci(d, "shell_height", None)
        if shell_h is not None and str(shell_h).strip() != "":
            L1 = float(shell_h)
            L2 = float(shell_h)
        else:
            L1 = float(_get_ci(d, "primary_height", 0.0))
            L2 = float(_get_ci(d, "secondary_height", L1))

        s1 = float(_get_ci(d, "primary_start", 0.0))
        s2 = float(_get_ci(d, "secondary_start", 0.0))

        q = _get_ci(d, "thickness_intercept", None)
        if q is None:
            q = _get_ci(d, "thickQ", 0.0)
        m = _get_ci(d, "thickness_slope", None)
        if m is None:
            m = _get_ci(d, "ThickM", 0.0)

        return GeoSettings(
            n_shells=int(_get_ci(d, "nShells", 0)),
            max_diam_mm=float(_get_ci(d, "maxDiam", 0.0)),
            focal_length_m=float(_get_ci(d, "focal_length", 0.0)),
            angular_shell_separation_deg=float(_get_ci(d, "angular_shell_separation_deg", 0.0)),
            primary_height_mm=float(L1),
            secondary_height_mm=float(L2),
            primary_start_mm=float(s1),
            secondary_start_mm=float(s2),
            thickness_intercept_mm=float(q),
            thickness_slope_per_mm=float(m),
            min_thick_mm=float(_get_ci(d, "minThick", 0.0)),
            wall_density_kg_m3=float(_get_ci(d, "WallDensity", 0.0)),
            profile_type=str(_get_ci(d, "profileType", "cc")).strip(),
        )


def compute_geometry(settings: GeoSettings, diam_file: Optional[Path] = None) -> Dict[str, Any]:
    nw = settings.n_shells
    if nw <= 0:
        raise ValueError("nShells must be > 0.")

    f_length_cm = settings.focal_length_m * 100.0
    f_length_mm = f_length_cm * 10.0

    shell_length1 = settings.primary_height_mm
    shell_length2 = settings.secondary_height_mm

    alkw = shell_length1 * math.tan(settings.angular_shell_separation_deg * PI / 180.0)

    size = nw + 2
    rmed = [0.0] * size
    rmax = [0.0] * size
    rmin = [0.0] * size
    thick = [0.0] * size
    theta = [0.0] * size
    mass = [0.0] * size
    area = [0.0] * size

    if diam_file is not None:
        vals: List[float] = []
        with diam_file.open("r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                line = _strip_fortran_comments(line).strip()
                if not line:
                    continue
                vals.append(float(line.split()[0]))
                if len(vals) >= nw:
                    break
        if len(vals) < nw:
            raise ValueError(f"Expected {nw} diameters in {diam_file}, got {len(vals)}")
        vals = [v / 2.0 for v in vals]
        if vals[0] < vals[-1]:
            vals = list(reversed(vals))
        for i in range(1, nw + 1):
            rmed[i] = vals[i - 1]
        rstart = max(rmed[1:nw + 1])
    else:
        if settings.max_diam_mm <= 0:
            raise ValueError("maxDiam not provided and no diameters file supplied.")
        rstart = settings.max_diam_mm / 2.0

    h1_start = settings.primary_start_mm
    h2_start = settings.secondary_start_mm

    rtoll = 1e-10

    r = rstart
    for i in range(1, nw + 2):
        if (diam_file is None) or (i == nw + 1):
            if i == nw + 1:
                spes = 0.0
            else:
                m = settings.thickness_slope_per_mm
                q = settings.thickness_intercept_mm
                spes = (q + m * 2.0 * r) / (1.0 - m * 2.0)
                spes = max(spes, settings.min_thick_mm)

            r = r - spes

            a_c, b_c = poly_coeff(settings.profile_type, r, f_length_mm)
            r0 = poly_invert(-shell_length1, r, a_c[0], a_c[1])

            while True:
                a_c, b_c = poly_coeff(settings.profile_type, r0, f_length_mm)
                r1 = poly_invert(-shell_length1, r, a_c[0], a_c[1])
                if abs(r1 - r0) <= rtoll:
                    break
                r0 = r1

            rmed[i] = r0

        r = rmed[i] - alkw

        theta[i] = math.atan2(rmed[i], f_length_mm) / 4.0
        a_c, b_c = poly_coeff(settings.profile_type, rmed[i], f_length_mm)

        rmax[i] = poly_profile(-shell_length1, rmed[i], a_c)
        rmin[i] = poly_profile(shell_length2, rmed[i], b_c)

        m = settings.thickness_slope_per_mm
        q = settings.thickness_intercept_mm
        spes = (q + m * 2.0 * rmax[i]) / (1.0 - m * 2.0)
        spes = max(spes, settings.min_thick_mm)
        thick[i] = spes

        if i != 1:
            if (rmin[i] + spes) > rmin[i - 1]:
                raise ValueError(
                    "Shell overlapping (at the exit pupil) for shell nr. {}:\n"
                    "  shell {} inner: {}\n"
                    "  shell {} outer: {}".format(i, i - 1, rmin[i - 1], i, rmin[i] + spes)
                )

        if i <= nw:
            apot_p = shell_length1 / math.cos(theta[i])
            apot_h = shell_length2 / math.cos(3.0 * theta[i])
            volume = spes * PI * ((rmed[i] + rmax[i]) * apot_p + (rmed[i] + rmin[i]) * apot_h)
            mass[i] = volume * settings.wall_density_kg_m3 / 1_000_000.0

        _ = poly_profile(h1_start, rmed[i], a_c)
        _ = poly_profile(-h2_start, rmed[i], a_c)

    total_area = 0.0
    for i in range(1, nw + 1):
        inner = max(rmed[i], rmax[i + 1] + thick[i])
        a = PI * (rmax[i] ** 2 - inner ** 2) / 100.0
        area[i] = a
        total_area += a

    return dict(
        rmax=rmax,
        rmed=rmed,
        rmin=rmin,
        thick=thick,
        theta=theta,
        area=area,
        mass=mass,
        f_length_cm=f_length_cm,
        f_length_mm=f_length_mm,
        shell_length1=shell_length1,
        shell_length2=shell_length2,
        total_area=total_area,
        total_mass=sum(mass[1:nw + 1]),
    )


def write_shell_struct(out_path: Path, geom: Dict[str, Any], n_shells: int) -> None:
    header = (
        "Nshell   Dmax(mm)   Dmid(mm)   Dmin(mm)   thickness(mm)   "
        "Angle(rad)   Area(cm^2)  Mass(kg)   primary_height(mm)   secondary_height(mm)"
    )
    with out_path.open("w", encoding="utf-8") as f:
        f.write(header + "\n")
        for k in range(1, n_shells + 1):
            f.write(
                f"{k:4d}"
                f"{geom['rmax'][k]*2:15.6f}"
                f"{geom['rmed'][k]*2:15.6f}"
                f"{geom['rmin'][k]*2:15.6f}"
                f"{geom['thick'][k]:15.6f}"
                f"{geom['theta'][k]:15.6f}"
                f"{geom['area'][k]:15.6f}"
                f"{geom['mass'][k]:15.6f}"
                f"{geom['shell_length1']:15.6f}"
                f"{geom['shell_length2']:15.6f}\n"
            )


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(description="Create shell geometry file (Python port of createGeo.f90).")
    ap.add_argument("settings", nargs="?", default="geoSettings.txt",
                    help="File namelist con &geoSettings (default: geoSettings.txt).")
    ap.add_argument("-d", "--diameters", default=None,
                    help="Optional file with diameters at the intersection plane (one per line).")
    ap.add_argument("-o", "--output", default="shellStruct.dat",
                    help="Output geometry file (default: shellStruct.dat).")
    args = ap.parse_args(argv)

    settings_path = Path(args.settings)
    if not settings_path.exists():
        print(f"Error: settings file not found: {settings_path}", file=sys.stderr)
        return 2

    d = read_geo_settings(settings_path, group="geoSettings")
    settings = GeoSettings.from_dict(d)

    diam_path = Path(args.diameters) if args.diameters else None
    if diam_path is not None and not diam_path.exists():
        print(f"Error: diameters file not found: {diam_path}", file=sys.stderr)
        return 2

    print("Starting Python createGeo (height-based settings)")
    print("Settings:", settings_path.resolve())
    if diam_path:
        print("Diameters:", diam_path.resolve())
    print(f"F_length (cm) = {settings.focal_length_m*100:.6f}")
    print(f"primary_height (mm) = {settings.primary_height_mm:.6f}")
    print(f"secondary_height (mm) = {settings.secondary_height_mm:.6f}")
    print(f"Shell number: {settings.n_shells}")

    geom = compute_geometry(settings, diam_file=diam_path)
    out_path = Path(args.output)
    write_shell_struct(out_path, geom, settings.n_shells)

    print(f"Total mass (kg): {geom['total_mass']:.6f}")
    print(f"Total collecting area (cm^2): {geom['total_area']:.6f}")
    print(f"Geometry file successfully created: {out_path.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
