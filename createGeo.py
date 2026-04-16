#!/usr/bin/env python3
"""
createGeo.py — Python translation of the Fortran program `createGeo.f90`.

What it does
------------
- Reads a Fortran namelist file (default: geoSettings.txt) containing the `geoSettings` group.
- Optionally reads a text file with diameters at the intersection plane (one per line).
- Computes shell geometry (rmax/rmed/rmin, thickness, angle, collecting area, mass).
- Writes `shellStruct.dat` in the current working directory (or a user-specified path).

Notes about fidelity
--------------------
- The Fortran code computes "collecting area" using rmax(i+1) (next shell). In the original loop this
  is referenced before it is known; this Python version computes area in a second pass after all shells
  are generated, which matches the apparent intent of the 2019/01/17 fix.
- The Fortran variables `h1_start` / `h2_start` are never assigned. Here they are mapped to
  the namelist entries `shell_H1_start` and `shell_H2_start` (default 0.0).
"""

from __future__ import annotations

import argparse
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any

PI = math.pi


# -------------------- Fortran namelist parsing --------------------

def _strip_fortran_comments(line: str) -> str:
    # Fortran comments typically start with '!' and go to end of line.
    return line.split("!")[0]


def _split_top_level_commas(s: str) -> List[str]:
    """Split by commas, ignoring commas inside single/double quotes."""
    out, buf = [], []
    in_sq = False
    in_dq = False
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


def _parse_fortran_scalar(val: str) -> Any:
    """Parse a scalar from a Fortran namelist assignment RHS."""
    v = val.strip()
    if not v:
        return v

    # logicals
    vl = v.lower()
    if vl in (".true.", "true"):
        return True
    if vl in (".false.", "false"):
        return False

    # strings
    if (v.startswith("'") and v.endswith("'")) or (v.startswith('"') and v.endswith('"')):
        return v[1:-1]

    # Fortran double exponent uses D or d
    if any(c in v for c in "dD"):
        try:
            return float(replace_fortran_exponent(v))
        except ValueError:
            return v

    # numbers
    try:
        if any(c in v for c in ".eE"):
            return float(v)
        return int(v)
    except ValueError:
        return v


def replace_fortran_exponent(s: str) -> str:
    # 1.0D+03 -> 1.0e+03
    return re.sub(r"([0-9])([dD])([+\-]?[0-9]+)", r"\1e\3", s)


def read_namelist(filepath: Path, group: str = "geoSettings") -> Dict[str, Any]:
    """
    Minimal Fortran namelist reader for scalar values.

    Looks for:
        &geoSettings
            key = value,
            ...
        /

    Returns a dict with case-preserving keys as found; lookups should be case-insensitive.
    """
    raw_lines = filepath.read_text(encoding="utf-8", errors="ignore").splitlines()
    cleaned = [_strip_fortran_comments(l).rstrip() for l in raw_lines]

    # Find start of group
    start_idx = None
    group_pat = re.compile(rf"^\s*&\s*{re.escape(group)}\b", re.IGNORECASE)
    for i, l in enumerate(cleaned):
        if group_pat.search(l):
            start_idx = i
            break
    if start_idx is None:
        raise ValueError(f"Namelist group '{group}' not found in {filepath}")

    # Collect until '/' or '&end'
    block_lines: List[str] = []
    for l in cleaned[start_idx:]:
        if re.search(r"^\s*/\s*$", l) or re.search(r"^\s*&\s*end\b", l, re.IGNORECASE):
            break
        block_lines.append(l)

    # Drop the leading "&group" line
    if block_lines:
        block_lines[0] = group_pat.sub("", block_lines[0])

    block = " ".join(block_lines)
    # Split into assignments by comma
    parts = _split_top_level_commas(block)

    out: Dict[str, Any] = {}
    for p in parts:
        if "=" not in p:
            continue
        k, v = p.split("=", 1)
        key = k.strip()
        val = _parse_fortran_scalar(v)
        out[key] = val
    return out


def _get_ci(d: Dict[str, Any], key: str, default: Any = None) -> Any:
    """Case-insensitive dict lookup."""
    kl = key.lower()
    for k, v in d.items():
        if k.lower() == kl:
            return v
    return default


# -------------------- Geometry model --------------------

def poly_coeff(profile_type: str, radius_ip: float, f_length_mm: float) -> Tuple[Tuple[float, float, float, float, float],
                                                                                Tuple[float, float, float, float, float]]:
    """
    Return (a_coeffs, b_coeffs) for the primary (surface 1) and secondary (surface 2).

    profile_type: 2-char string:
      c -> cone, p -> parabola, h -> hyperbola, u -> user-defined (fixed polynomial coefficients)
    """
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
        # as in the Fortran code: fixed "user" polynomial coefficients
        a2 = 0.0002005
        a3 = 0.0003318
        a4 = 0.0002711
        a5 = 0.00009652
    elif p1 == "h":
        raise NotImplementedError("Hyperbolic primary surface (profileType[0]='h') is not implemented (as in Fortran).")
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
        # Matches the Fortran expression (secondary hyperbola coefficient)
        b2 = 2.0 * (radius_ip * math.tan(angle_two)) / (f_length_mm + radius_ip / math.tan(2.0 * angle_one))
    else:
        raise ValueError(f"Secondary profile type '{profile_type[1]}' not recognized (valid: c, h, p, u).")

    return (a1, a2, a3, a4, a5), (b1, b2, b3, b4, b5)


def poly_invert(z: float, r1: float, c1: float, c2: float) -> float:
    """
    Return the radius at intersection plane, given a (mostly) quadratic polynomial profile
    and a value r1 at coordinate z.

    Fortran:
      deltaroot = sqrt((c1*z)^2 - 4*c2*z^2 + 4*r1^2)
      r0 = (-c1*z + deltaroot)/2
    """
    delta = (c1 * z) ** 2 - 4.0 * c2 * (z ** 2) + 4.0 * (r1 ** 2)
    if delta < 0:
        # Numerical guard
        delta = 0.0
    return (-c1 * z + math.sqrt(delta)) / 2.0


def poly_profile(z: float, r0: float, c: Tuple[float, float, float, float, float]) -> float:
    """Return the radius at coordinate z from the polynomial profile coefficients."""
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
    shell_length_mm: float
    thickness_slope: float
    thickness_intercept: float
    min_thick_mm: float
    wall_density_kg_m3: float
    profile_type: str
    shell_h1_start_mm: float = 0.0
    shell_h2_start_mm: float = 0.0

    @staticmethod
    def from_namelist(d: Dict[str, Any]) -> "GeoSettings":
        return GeoSettings(
            n_shells=int(_get_ci(d, "nShells", 0)),
            max_diam_mm=float(_get_ci(d, "maxDiam", 0.0)),
            focal_length_m=float(_get_ci(d, "focal_length", 0.0)),
            angular_shell_separation_deg=float(_get_ci(d, "angular_shell_separation_deg", 0.0)),
            shell_length_mm=float(_get_ci(d, "shell_length", 0.0)),
            thickness_slope=float(_get_ci(d, "thickness_slope", 0.0)),
            thickness_intercept=float(_get_ci(d, "thickness_intercept", 0.0)),
            min_thick_mm=float(_get_ci(d, "minThick", 0.0)),
            wall_density_kg_m3=float(_get_ci(d, "WallDensity", 0.0)),
            profile_type=str(_get_ci(d, "profileType", "cc")).strip(),
            shell_h1_start_mm=float(_get_ci(d, "shell_H1_start", 0.0)),
            shell_h2_start_mm=float(_get_ci(d, "shell_H2_start", 0.0)),
        )


def compute_geometry(settings: GeoSettings, diam_file: Optional[Path] = None) -> Dict[str, List[float]]:
    """
    Compute geometry arrays (1-based indexing; index 0 unused to resemble Fortran).
    Returns dict with keys: rmax, rmed, rmin, thick, theta, area, mass.
    """
    nw = settings.n_shells
    if nw <= 0:
        raise ValueError("nShells must be > 0 (read from geoSettings namelist).")

    # Focal length conversions as in Fortran:
    #   F_LENGTH (cm) = focal_length * 100
    #   F_LENGTHmm = F_LENGTH * 10
    f_length_cm = settings.focal_length_m * 100.0
    f_length_mm = f_length_cm * 10.0

    shell_length1 = settings.shell_length_mm
    shell_length2 = settings.shell_length_mm

    alkw = shell_length1 * math.tan(settings.angular_shell_separation_deg * PI / 180.0)

    # allocate (nw + 1) shells + 1 extra for central block => total indices 1..nw+1
    size = nw + 2
    rmed = [0.0] * size
    rmax = [0.0] * size
    rmin = [0.0] * size
    thick = [0.0] * size
    theta = [0.0] * size
    mass = [0.0] * size
    area = [0.0] * size

    # determine starting radius at entrance pupil
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

        # convert diameter->radius and ensure decreasing (outer-to-inner)
        vals = [v / 2.0 for v in vals]
        if vals[0] < vals[-1]:
            vals = list(reversed(vals))
        for i in range(1, nw + 1):
            rmed[i] = vals[i - 1]
        rstart = max(rmed[1:nw + 1])
    else:
        if settings.max_diam_mm <= 0:
            raise ValueError("maxDiam not provided in geoSettings and no diameter file was supplied.")
        rstart = settings.max_diam_mm / 2.0

    # Map Fortran's unassigned h1_start/h2_start to the namelist parameters.
    h1_start = settings.shell_h1_start_mm
    h2_start = settings.shell_h2_start_mm

    # iterative tolerance
    rtoll = 1e-10

    # build shells
    r = rstart
    for i in range(1, nw + 2):  # 1..nw+1
        # If diameters were NOT provided, compute rmed(i) recursively from outer radius
        if (diam_file is None) or (i == nw + 1):
            if i == nw + 1:
                spes = 0.0
            else:
                m = settings.thickness_slope
                q = settings.thickness_intercept
                spes = (q + m * 2.0 * r) / (1.0 - m * 2.0)
                spes = max(spes, settings.min_thick_mm)

            # radius at mirror surface (internal)
            r = r - spes

            # recursive inversion (as in Fortran)
            a_c, b_c = poly_coeff(settings.profile_type, r, f_length_mm)
            r0 = poly_invert(-shell_length1, r, a_c[0], a_c[1])

            while True:
                a_c, b_c = poly_coeff(settings.profile_type, r0, f_length_mm)
                r1 = poly_invert(-shell_length1, r, a_c[0], a_c[1])
                if abs(r1 - r0) <= rtoll:
                    break
                r0 = r1

            rmed[i] = r0

        # external radius at exit pupil for next cycle
        r = rmed[i] - alkw

        theta[i] = math.atan2(rmed[i], f_length_mm) / 4.0
        a_c, b_c = poly_coeff(settings.profile_type, rmed[i], f_length_mm)

        rmax[i] = poly_profile(-shell_length1, rmed[i], a_c)
        rmin[i] = poly_profile(shell_length2, rmed[i], b_c)

        # thickness as a function of the entrance pupil radius (rmax)
        m = settings.thickness_slope
        q = settings.thickness_intercept
        spes = (q + m * 2.0 * rmax[i]) / (1.0 - m * 2.0)
        spes = max(spes, settings.min_thick_mm)
        thick[i] = spes

        # overlap check at exit pupil
        if i != 1:
            if (rmin[i] + spes) > rmin[i - 1]:
                raise ValueError(
                    "Shell overlapping (at the exit pupil) for shell nr. {}:\n"
                    "  shell {} inner: {}\n"
                    "  shell {} outer: {}".format(i, i - 1, rmin[i - 1], i, rmin[i] + spes)
                )

        # mass (only real shells, not the extra block)
        if i <= nw:
            apot_p = shell_length1 / math.cos(theta[i])
            apot_h = shell_length2 / math.cos(3.0 * theta[i])
            volume = spes * PI * ((rmed[i] + rmax[i]) * apot_p + (rmed[i] + rmin[i]) * apot_h)
            mass[i] = volume * settings.wall_density_kg_m3 / 1_000_000.0  # kg

        # (rip1/rip2/obstruction terms are computed in Fortran but not used later;
        #  keep them here for completeness / future extension.)
        _ = poly_profile(h1_start, rmed[i], a_c)
        _ = poly_profile(-h2_start, rmed[i], a_c)

    # Collecting area (second pass, after rmax(i+1) is known)
    total_area = 0.0
    for i in range(1, nw + 1):
        inner = max(rmed[i], rmax[i + 1] + thick[i])  # matches the Fortran expression (uses current shell thickness)
        a = PI * (rmax[i] ** 2 - inner ** 2) / 100.0  # cm^2
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
        alkw=alkw,
        total_area=total_area,
        total_mass=sum(mass[1:nw + 1]),
    )


def write_shell_struct(out_path: Path, geom: Dict[str, List[float]], n_shells: int) -> None:
    header = (
        "Nshell   Dmax(mm)   Dmid(mm)   Dmin(mm)   thickness(mm)   "
        "Angle(rad)   Area(cm^2)  Mass(kg)\tShell_length1\tShell_length2"
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
    ap.add_argument("diameters", nargs="?", default=None,
                    help="Optional file with diameters at the intersection plane (one per line).")
    ap.add_argument("-s", "--settings", default="geoSettings.txt",
                    help="Fortran namelist file containing &geoSettings (default: geoSettings.txt).")
    ap.add_argument("-o", "--output", default="shellStruct.dat",
                    help="Output geometry file (default: shellStruct.dat).")
    args = ap.parse_args(argv)

    settings_path = Path(args.settings)
    if not settings_path.exists():
        print(f"Error: settings file not found: {settings_path}", file=sys.stderr)
        return 2

    d = read_namelist(settings_path, group="geoSettings")
    settings = GeoSettings.from_namelist(d)

    diam_path = Path(args.diameters) if args.diameters else None
    if diam_path is not None and not diam_path.exists():
        print(f"Error: diameters file not found: {diam_path}", file=sys.stderr)
        return 2

    print("Starting Python createGeo (port of createGeo.f90)")
    print("Settings:", settings_path.resolve())
    if diam_path:
        print("Diameters:", diam_path.resolve())
    print(f"F_length (cm) = {settings.focal_length_m*100:.6f}")
    print(f"shell_length (mm) = {settings.shell_length_mm:.6f}")
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
