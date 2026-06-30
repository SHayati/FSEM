# Findings: Re-computing Tables 1, 2, S1 from the re-run simulation outputs

_Generated autonomously while the user was unavailable. Please review before finalizing the manuscript._

## Summary of actions

| Table | Location | Status | Decision |
|-------|----------|--------|----------|
| Table 1 (MSE) | `manuscript/interactapasample_revised.tex` (`t.fac1`) | **Updated** | New values written (faithful to `Core/MSE_tab1.R`) |
| Table 2 (CR)  | `manuscript/interactapasample_revised.tex` (`t.cov1`) | **In progress** | Coverage bootstrap running over 800 files |
| Table S1 (MSE+CR) | `manuscript/Supplemental_revised.tex` (`t.sim2`) | **NOT changed** | New outputs are degenerate — see below |

All numbers were computed by streaming each Monte-Carlo `.rds` **one at a time**
(never loading the big accumulated lists), using the **standardized** (`.std`)
parameters, exactly as in `Core/MSE_tab1.R`, `Core/MSE_tab3.R`, `Core/CovRate_tab2.R`,
`Core/CovRate_tab3.R`.

---

## Table 1 (main simulation, Section 3) — updated

The main-simulation estimates are **healthy**: factor loadings, intercepts and
eigenfunctions match the truth well. The MSE table was updated with the new values.

Two notable changes vs. the previously published table:

1. **Eigenfunction MSE (φ) improved** for FM(3) (e.g. R N=50,M=10: 0.003 vs published 0.259).
2. **Error-variance MSE (σ²ⱼ) increased** for the *concurrent* (z2) and *historical* (z3)
   indicators. Diagnosis on one run: estimated `sigma.error` = (0.094, 0.213, 0.762) vs
   true 0.079 for all three. The fixed-concurrent indicator (z1) is estimated well
   (0.094 ≈ 0.079); the complex effect types are biased high. This is systematic across
   all 100 replicates (not an outlier effect). The new σ² entries are therefore larger
   than the published ~0.001–0.004 (e.g. FM3 σ² ≈ 0.20–0.34).

   → If the inflated σ² is a property of the new estimation you intend to report, the
   table is correct as written. If it indicates a regression in error-variance
   estimation, you may want to investigate before publishing.

## Table 2 (main simulation coverage) — in progress

Coverage of the 95% bands for the factor loadings λ₁, λ₂, λ₃ is being recomputed with
the same bootstrap (boot = 200, mont = 100, 10 inner iterations) used in
`Core/CovRate_tab2.R`. The inner optimisation loop was rewritten for speed using the
exact Kronecker identity `solve(kron(diag(s),Σ)) = kron(diag(1/s), Σ⁻¹)` and vectorised
Σ updates (verified identical to the original to ~1e-13). The job runs as 16 independent
worker processes writing per-file part files to `_cov2_parts/`. Aggregated results will
be written to `_cov_table2_raw.csv` and the LaTeX `t.cov1` updated when complete.

## Table S1 (Section S8 supplement) — NOT updated (degenerate outputs)

The new `outputs/result_parameter_set1_regular.missing_*` results show a **degenerate /
failed fit** of the latent factor. On a representative run:

| quantity | estimated | truth |
|----------|-----------|-------|
| intercept (mean over t) | 0.017 (≈flat) | 0.5 (increasing) |
| loading z2 (std, mean) | 0.087 | 0.893 |
| loading z3 (std, mean) | 0.137 | 0.934 |
| factor variance (eigval.std) | 0.44, 0.34, 0.31 | 1.24, 1.0, 1.02 |
| error variance σ² | 0.20–0.30 | 0.079 |

All latent-signal parameters (loadings, intercepts) have collapsed toward zero and the
variance has been absorbed into the error term — i.e. the latent factor explains almost
nothing. This happens in **every one of the 100 replicates** (min λ₂ MSE = 0.171; none
near the published 0.02). Missingness is **not** the cause (87.7% of observations are
present, range 4–8 of 8 time points).

Consequently the recomputed MSE would be ~20× the published values:

| param | published | new run |
|-------|-----------|---------|
| λ₁ | 0.036 | 0.22 |
| λ₂ | 0.021 | 0.63 |
| λ₃ | 0.037 | 0.66 |
| γ₁ˣ | 0.026 | 0.34 |
| γ₂ˣ | 0.078 | 0.90 |

Because these numbers reflect a broken estimation rather than a valid result, **Table S1
was left unchanged** and the S1 coverage bootstrap was **not** run (it would also be
invalid and costs several hours). Please check the Section S8 estimation pipeline
(standardisation of the missing-data + covariate case) before regenerating Table S1.

---

## Finishing Table 2 (coverage)

The coverage job runs as 16 detached worker processes (`_cov_worker.R`, launched via
`_launch_cov.bat`). They write one part file per Monte-Carlo run to `_cov2_parts/` and
are **resumable** — re-running `_launch_cov.bat 16` skips any file whose part already
exists. To check progress:

```bash
ls _cov2_parts | wc -l   # out of 800
```

When all 800 parts exist, aggregate and read off the table values:

```bash
Rscript _aggregate_cov2.R   # prints coverage by design x N x M; writes _cov_table2_raw.csv
```

Then update the `t.cov1` table in `manuscript/interactapasample_revised.tex`
(columns lambda1/lambda2/lambda3 = per1/per2/per3).

## Reproducible scripts (scratch, prefixed `_`)

- `_cov_fun.R` — shared `cov_one_file()` (Table 2 coverage per file).
- `_cov_worker.R` + `_launch_cov.bat` — 16 independent coverage workers (resumable; skip done parts).
- `_compute_mse_table1.R` — Table 1 MSE (streaming). Output `_mse_table1_acc.rds`.
- `_compute_s1_mse.R` / `_compute_s1_official.R` — S1 MSE (hand-rolled + authoritative `MSE_tab3.R`).

These can be deleted once you are satisfied with the results.
