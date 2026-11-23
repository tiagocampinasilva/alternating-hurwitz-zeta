# Explicit Uniform Remainder Estimate for the Alternating Hurwitz Zeta Function in a Joint Large-Parameter Regime

Computational supplement to the manuscript by **Tiago Manuel Serafim Campina da Silva**. The manuscript gives 24 November 2026 as its version date. The reproduction recorded here was performed on 22 November 2026.

## Run the three programs

From this directory, create a virtual environment using a Python installation with the standard venv module. These commands were tested with CPython 3.14.4 on macOS ARM64. On the project Mac, the interpreter used was `/Library/Frameworks/Python.framework/Versions/3.14/bin/python3`. If necessary, substitute that path for `python3` in the first command.

```sh
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python exact_arithmetic.py
.venv/bin/python identity_checks.py
.venv/bin/python numerical_checks.py --abel-plana --output reproduced.csv
```

The only external Python dependencies are mpmath 1.3.0 (numerical evaluation) and SymPy 1.14.0 (symbolic checks). The exact program uses only the Python standard library, including `fractions.Fraction`, and can run before installing the two packages. Each program returns a nonzero exit status on a failed check. The example writes a separate CSV. Omitting `--output` replaces `Hurwitz_Numerical_Results.csv` beside the numerical script after the checks succeed. Omitting `--abel-plana` skips the three quadratures.

## Regenerate the Python files with Mathematica

Open `Hurwitz_Python_Code_Generation.nb` in Wolfram Mathematica and select **Evaluation → Evaluate Notebook**. The notebook contains the Wolfram definitions, the numbers printed in Table 1, the Python source templates and the export code. It constructs the scripts in its own directory. Regeneration requires write access there. The Python programs can run separately from Mathematica and from the manuscript files.

The notebook evaluates independent Wolfram versions of the algebraic residuals, rational constants, numerical grid and quadratures. Its final cell reads the CSV with the distributed filename, matches rows by exact `(q, x, working_dps)` keys, and compares all 17 numerical fields. For a fresh calculation and comparison in Mathematica, rerun the numerical program with its default output name, then evaluate the final cell after the Wolfram grid cells. The recorded Wolfram kernel version was **15.0.1 for Mac OS X ARM (64-bit), July 2, 2026**.

## Files, equations and precision

| File | Calculation |
| --- | --- |
| numerical_checks.py | Equations (1)–(6), the exact constants in (16), (17) and Appendix 1, and optional numerical quadrature of (7). |
| identity_checks.py | Prints each SymPy residual. It checks the Abel–Plana weight and origin limit, the two kernel forms, Taylor coefficients, angular derivatives, differential identities for F, the reduction of G, the polynomial in (14) and four kernel moments. |
| exact_arithmetic.py | Derives the polynomial coefficients, rational moment bounds, residual constants and strict margins. It also checks both transition estimates, the 4,096-term floor sum and the rational comparisons used in the Appendix 1 counterexample. |
| Hurwitz_Numerical_Results.csv | 287 main-grid evaluations at 55 working digits and seven reevaluations at 75 digits. Each group of 41 rows is followed by its maximizing point. |
| Hurwitz_Python_Code_Generation.nb | Self-contained program generation and independent Mathematica comparison. |
| requirements.txt | Versions of the two external Python libraries used in this reproduction. |

The grid uses the seven real shifts in Table 1 and the 41 rationals x = −1 + j/20 for j = 0,…,40. The implementation uses the principal logarithm, t = πq + x√q and s = 1/2 + it. At s = 1 the numerical zeta routine evaluates the paired series from equation (2) to avoid the individual Hurwitz poles. Every precision run reconstructs its inputs from decimal strings representing exact rational coordinates.

The CSV suffixes `_re` and `_im` denote real and imaginary parts. `approx` is F + G/√q, `R` is Z minus that approximation, `abs_one_term` is |Z − F|, `abs_R` is |R| and `q_abs_R` is q|R|. The four bound columns are absolute-error bounds or comparison values in the same units as |R|. They give 4.3/q, the exact coefficient from (16) divided by q, the exact coefficient from (17) divided by q, and the Appendix 1 coefficient divided by q, respectively. CSV numbers use a decimal point. Read their text directly into an arbitrary-precision type to keep the stored digits. Working precision is the number of decimal digits used in a calculation. Cancellation in the remainder reduces its effective accuracy.

## Recorded comparisons

- Recompiling the complete LaTeX source reproduced the text and layout of all nine PDF pages. Independent page renders at 144 dpi were pixel-identical, including equations, both figures, Table 1 and the appendix.
- All seven Table 1 maxima occurred at x = −1. Each agreed with the printed value to 12 decimal places. The largest absolute difference from print was 4.021 × 10⁻¹³ or less. The highest sampled q|R| was 0.083681979527711946562461945089… at q = 100. The seven 75-digit Python maxima and their Mathematica counterparts agree when rounded to 60 decimal places.
- Mathematica matched all 294 Python rows. Across 17 quantities per row, the largest absolute differences were approximately 3.0081 × 10⁻⁴⁶ at 55 working digits and 3.8116 × 10⁻⁶⁷ at 75 digits. The absolute tolerances in the notebook are 10⁻⁴⁴ and 10⁻⁶⁴. A second Python run in a fresh virtual environment reproduced the distributed CSV byte for byte.
- The 41 plotted samples in Figure 1 and both sets of 123 markers in Figure 2 matched the coordinates in the LaTeX source to their printed precision. At q = 100 the ratio of the one-term to two-term absolute error ranged from 58.94136 to 471.01661 (rounded outwards).
- All 21 individual symbolic residuals simplified to zero in SymPy and Mathematica. The manuscript describes nine identities. The program separates individual coefficients and differentiated consequences and includes the polynomial and moment evaluations in its count. The rational checks ran successfully, including the strict final margins.
- For equation (7) both implementations integrated at q = 100 and x = −1, 0, 1 at 60 working digits and compared with the Hurwitz evaluation at 70 digits. Mathematica used the manuscript's accuracy and precision goals of 35 and reproduced its three printed discrepancies. The Python and Wolfram quadrature routines give different absolute discrepancies from their respective zeta evaluations.

| x | Python absolute discrepancy | Mathematica absolute discrepancy |
| ---: | ---: | ---: |
| −1 | 6.216 × 10⁻⁶¹ | 7.531 × 10⁻⁵⁸ |
| 0 | 1.047 × 10⁻⁶⁰ | 1.672 × 10⁻⁵⁴ |
| 1 | 1.567 × 10⁻⁶⁰ | 1.187 × 10⁻⁵³ |

## Exact constants and scope

| Quantity | Exact value |
| --- | ---: |
| Split residual bound B_s | 689147/1575000 |
| Whole-axis residual bound B_g | 22061/84000 |
| Multiplier-3 residual bound B₃ | 275663/1050000 |
| Integral upper bound J̄ | 533948559110614545113/69278297361197722560 |
| J̄/2 + B_s, equation (16) | 371609064358084350147257/86597871701497153200000 |
| Margin below 43/10 | 761783958353408612743/86597871701497153200000 |
| J̄/2 + B_g, equation (17) | 14258445238530472531037/3463914868059886128000 |
| J̄/2 + B₃, Appendix 1 | 356452924783991052355217/86597871701497153200000 |

The even Gaussian recurrence propagates the upper bound 67/80 for g₀. Equation (17) uses a whole-axis kernel estimate. Appendix 1 uses a multiplier of 3 on the short interval and presents this refinement as provisional. Both exact alternative coefficients lie below 4117/1000. The appendix mentions 24.8/q and 11.3/q as historical working bounds.

The numerical comparisons cover the sampled points. The uniform assertion for every real q ≥ 100 and |x| ≤ 1 additionally uses the manuscript's analytical arguments for contour deformation, differentiation and domination.

## Manuscript files used

SHA-256 digests identify the source pair used for this reproduction.

```text
327ca3b7a501f7ce5bff239c377dd3e12b3956f9369e8a05471bd97652ea74ee  HurwitzZetaManuscript.tex
6999455d483bdcc79310452d3d28b9bada62eaa37d4f7c6eb99eaef2c579e1ad  HurwitzZetaManuscript.pdf
```
