# Explicit Uniform Remainder Estimate for the Alternating Hurwitz Zeta Function in a Joint Large-Parameter Regime

**Tiago Manuel Serafim Campina da Silva**

The [manuscript](HurwitzZetaManuscript.pdf) links to the individual data files, figures and verification scripts in release v1.0.0 of the [GitHub repository](https://github.com/tiagocampinasilva/alternating-hurwitz-zeta). The computational supplement has DOI [10.5281/zenodo.24483843](https://doi.org/10.5281/zenodo.24483843). Open [Supplement.pdf](Supplement.pdf) to read the article with its data, figures and verification code. On pages 6 and 7, click a figure caption, the Table 1 caption or a linked verification method to jump to the corresponding appendix. Use “Return to the manuscript” to go back.

| Manuscript item | Data and code |
| --- | --- |
| Table 1 | [Seven maxima, printed values and absolute discrepancies](Verification/table_1.csv) · [Full Wolfram grid and reevaluations](Verification/verification_results.csv) |
| Figure 1, scaled error at q = 100 | [Vector figure](Figures/figure_1.pdf) · [41 plotted points](Figures/figure_1_data.csv) |
| Figure 2, errors for one and two terms at q = 100, 400 and 1600 | [Three-panel vector figure](Figures/figure_2.pdf) · [123 pairs of errors](Figures/figure_2_data.csv) |
| Numerical comparisons and quadrature | [Wolfram Language source](Verification/verification.wl) · [Python/mpmath source](Python/numerical_checks.py) |
| Symbolic identities | [SymPy checks](Python/identity_checks.py) |
| Uniform bound and rational refinements | [Exact arithmetic in Python](Python/exact_arithmetic.py) · [Exact coefficients and strict margins](Verification/exact_coefficients.csv) |
| Verification record | [Check results from the package build](Verification/Check_Results.csv) · [Archived notebook](Verification/Hurwitz_Final_Verification.nb) · [Archived report](Verification/Hurwitz_Final_Verification.pdf) |

## Reproduce the calculations

1. Make a working copy of this directory. Open [Reproduce.nb](Reproduce.nb) in Mathematica and evaluate its input cells in order. The first calculation runs the verification against the archived LaTeX table, reference grid and original notebook. It writes the Wolfram results CSV inside the working copy. The second calculation exports the two figures and their accompanying CSVs.
2. For a separate Python calculation, follow the [environment and run instructions](Python/README.md). The package includes the [dependency versions](Python/requirements.txt) and the [Mathematica notebook that exports the Python sources](Python/Hurwitz_Python_Code_Generation.nb).
3. The [editable manuscript](HurwitzZetaManuscript.tex) links to the individual files in the versioned GitHub repository. Compile it from the package root. The [archived source pair](Verification/Manuscrito/) identifies the version used for the numerical checks.

## What is checked

The main estimate is **|R| < 4.3/q** for real q ≥ 100 and |x| ≤ 1, with the normalization and functions defined in the manuscript. The proof supplies estimates on the whole domain. The exact arithmetic checks the rational coefficient in (16), the refinement in (17) and the separate coefficient in Appendix 1. The manuscript describes the appendix refinement as provisional.

Table 1 reports the largest sampled q|R| at each of seven q values, with x = −1 + j/20 for j = 0,…,40. The calculation uses 55 working digits and repeats the seven maximizing points at 75 working digits. Both figure CSVs use the 55-digit evaluations and retain the stored decimal strings. Figure 2 uses a common logarithmic scale for its three panels. The five [additional plots](Verification/Plots/) belong to the archived verification run.

The Wolfram verification directory and the Python supplement retain their respective input grids, which have different columns. The [archived log](Verification/Verification_Log.txt) contains an unevaluated formatting expression in the line about the Abel–Plana discrepancy. The [three numerical discrepancies](Verification/Abel_Plana.csv) are exported separately from the fresh calculation.

Both PDFs were compiled from LaTeX with Tectonic through Wolfram Mathematica. The notebook, standalone figures and verification data accompany the source. The manuscript version is dated 24 November 2026. SHA-256 digests are listed in [SHA256SUMS.txt](SHA256SUMS.txt).
