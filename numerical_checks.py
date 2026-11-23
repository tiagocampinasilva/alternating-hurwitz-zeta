"""Reproduce Table 1 and the sampled bounds in equations (6), (16) and (17)."""
import argparse
import csv
import os
from decimal import Decimal
from fractions import Fraction
from pathlib import Path
from tempfile import NamedTemporaryFile

from mpmath import mp

# (q, x at the maximum, printed q|R|), from Table 1.
TABLE = (
    ("100", "-1", "0.083681979528"),
    ("100.001", "-1", "0.083681976829"),
    ("100.5", "-1", "0.083680636233"),
    ("101", "-1", "0.083679304576"),
    ("150", "-1", "0.083587445357"),
    ("400", "-1", "0.083445093725"),
    ("1600", "-1", "0.083350517492"),
)
GRID_X = tuple(Decimal(j) / 20 for j in range(-20, 21))
MAIN = Fraction(43, 10)
SPLIT = Fraction(371609064358084350147257, 86597871701497153200000)
GLOBAL = Fraction(14258445238530472531037, 3463914868059886128000)
APPENDIX = Fraction(356452924783991052355217, 86597871701497153200000)
FIELDS = ("q", "x", "working_dps", "Z_re", "Z_im", "F_re", "F_im", "G_re", "G_im",
          "approx_re", "approx_im", "R_re", "R_im", "abs_one_term", "abs_R", "q_abs_R",
          "main_bound", "split_bound", "global_bound", "appendix_bound")


def rational_value(fraction):
    return mp.mpf(fraction.numerator) / fraction.denominator


def alternating_zeta(s, q):
    """Equation (2), with the paired series at its removable pole."""
    if q <= 0 or mp.re(s) <= 0:
        raise ValueError("Equation (2) requires q > 0 and Re(s) > 0.")
    if s == 1:
        return mp.nsum(lambda k: 1 / ((2*k + q) * (2*k + q + 1)), [0, mp.inf])
    return (mp.zeta(s, q / mp.mpf(2)) - mp.zeta(s, (q + 1) / mp.mpf(2))) / mp.power(2, s)


def normalized_z(q, x):
    t = mp.pi * q + mp.sqrt(q) * x
    s = mp.mpf("0.5") + mp.j*t
    return mp.exp(mp.j*t*mp.log(q)) * alternating_zeta(s, q)


def profile(x):
    angle = mp.exp(-mp.j*mp.pi/4)
    phase = mp.exp(-mp.j*x*x/(2*mp.pi))
    return (1 + mp.j)/2 * phase * mp.erfc(-angle*x/mp.sqrt(2*mp.pi))


def correction(x, f):
    return mp.mpf("0.5") + mp.j/(6*mp.pi) - x*x*(1 - mp.j*x*f)/(6*mp.pi**2)


def evaluate(q_text, x_text, digits):
    with mp.workdps(digits):
        q, x = mp.mpf(q_text), mp.mpf(str(x_text))
        z = normalized_z(q, x)
        f = profile(x)
        g = correction(x, f)
        approximation = f + g / mp.sqrt(q)
        residual = z - approximation
        absolute = abs(residual)
        scaled = q*absolute
        bounds = tuple(rational_value(c)/q for c in (MAIN, SPLIT, GLOBAL, APPENDIX))
        if not all(mp.isfinite(v) and absolute < v for v in bounds):
            raise ArithmeticError(f"Sampled bound failed at q={q_text}, x={x_text}.")
        values = (z.real, z.imag, f.real, f.imag, g.real, g.imag,
                  approximation.real, approximation.imag, residual.real, residual.imag,
                  abs(z-f), absolute, scaled, *bounds)
        return (q_text, str(x_text), digits, *(mp.nstr(v, digits) for v in values)), scaled


def manuscript_grid():
    rows = []
    print("q         x     q|R| at the maximizing point (75 working digits)")
    for q, expected_x, printed in TABLE:
        samples = [evaluate(q, x, 55) for x in GRID_X]
        rows.extend(row for row, _ in samples)
        best, peak55 = max(samples, key=lambda sample: sample[1])
        x_peak = best[1]
        repeat, peak75 = evaluate(q, x_peak, 75)
        rows.append(repeat)
        with mp.workdps(75):
            if Decimal(x_peak) != Decimal(expected_x) or abs(peak75 - mp.mpf(printed)) > mp.mpf("5e-13"):
                raise ArithmeticError(f"Table 1 mismatch for q={q}.")
            if abs(peak55 - peak75) >= mp.mpf("1e-44"):
                raise ArithmeticError(f"55/75-digit comparison failed for q={q}.")
            print(f"{q:<9} {x_peak:<5} {mp.nstr(peak75, 61)}")
    return rows


def abel_plana_check():
    """Normalized equation (7), integrated on the positive real axis at q=100."""
    for x in (-1, 0, 1):
        with mp.workdps(60):
            q = mp.mpf(100)
            root = mp.sqrt(q)
            s = mp.mpf("0.5") + mp.j*(mp.pi*q + x*root)

            def integrand(y):
                if y == 0:
                    return -mp.j*s/(mp.pi*q)
                plus = -s*mp.log1p(mp.j*y/q)
                minus = -s*mp.log1p(-mp.j*y/q)
                return mp.exp(plus-mp.pi*y)*mp.expm1(minus-plus)/mp.expm1(-2*mp.pi*y)

            integral = (mp.mpf("0.5") + mp.j*mp.quad(integrand, [0, 1, root, 4*root, mp.inf])) / root
        with mp.workdps(70):
            difference = abs(integral - normalized_z(q, x))
            print(f"Abel-Plana, x={x}: absolute difference {mp.nstr(difference, 8)}")
            if not difference < mp.mpf("1e-35"):
                raise ArithmeticError(f"Equation (7) comparison failed at x={x}.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path,
                        default=Path(__file__).resolve().with_name("Hurwitz_Numerical_Results.csv"))
    parser.add_argument("--abel-plana", action="store_true", help="Also integrate equation (7).")
    args = parser.parse_args()
    rows = manuscript_grid()
    if args.abel_plana:
        abel_plana_check()
    with NamedTemporaryFile(mode="w", encoding="utf-8", newline="", dir=args.output.parent,
                            delete=False) as stream:
        temporary = Path(stream.name)
        try:
            writer = csv.writer(stream)
            writer.writerow(FIELDS)
            writer.writerows(rows)
            stream.close()
            os.replace(temporary, args.output)
        finally:
            temporary.unlink(missing_ok=True)
    print(f"{len(rows)} evaluations saved to {args.output}")


if __name__ == "__main__":
    main()
