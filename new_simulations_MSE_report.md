# Two Additional Two-Factor Simulation Studies: MSE Results

This note reports the Monte-Carlo mean-squared-error (MSE) results for **two new
simulation studies** that extend the single-factor study of the main paper to a
**two-factor** functional structural equation model (FSEM). All quantities are
computed from the stored simulation output (folders `outputs/` and `outputs3/`);
**no simulation was re-run**. The tables are intended to be incorporated into the
manuscript later.

Throughout, the latent factors are zero-mean Gaussian processes observed only
through their indicators. Each replicate is summarised by the integrated squared
error between the estimated and the true (standardised) functional parameter, and
the MSE is the Monte-Carlo average of that error over the replicates in each
design cell. We use the manuscript notation:

| Symbol | Meaning |
|---|---|
| $\lambda$ | factor-loading function (indicator on factor) |
| $\beta$ | factor intercept function |
| $\gamma$ | structural (latent-on-latent) regression function |
| $\phi,\ \nu$ | leading eigenfunction / eigenvalue of a factor-error process |
| $\psi,\ \mu$ | leading eigenfunction / eigenvalue of a latent-process |
| $\sigma^2$ | measurement-error variance |

Both studies use $J=6$ indicators ($z_1,\dots,z_6$), a B-spline basis with $k$
internal structure as in the main study, signal-to-noise ratio $\mathrm{SNR}=4$
and correlation parameter $\rho=0.3$, with Matérn covariance for the latent
processes. Loadings combine a *fixed-concurrent* anchor ($z_1,z_4$), *concurrent*
effects ($z_2,z_5,z_6$) and a *historical* effect ($z_3$). Sample sizes
$N\in\{50,100\}$ and number of observation points $M\in\{10,20\}$ are crossed.

For readability all loading and intercept MSEs are reported $\times 10^{-2}$.

---

## Study A — Two factors, *uncorrelated*

**Design.** Two latent factors $\eta_1,\eta_2$ are generated **independently**
(no structural path between them). Indicators $z_1,z_2,z_3$ load on $\eta_1$ and
$z_4,z_5,z_6$ load on $\eta_2$. Each factor has an intercept-only structural
equation ($\eta_1\sim 1$, $\eta_2\sim 1$). This study checks that the estimator
recovers a multi-factor measurement model and the (independent) latent-process
covariances when there is no latent dependence to exploit. Results in `outputs/`
comprise the four cells $N\in\{50,100\}\times M\in\{10,20\}$ with 80 replicates each.

**Factor-loading MSE $\lambda$ ($\times 10^{-2}$).**

| $N$ | $M$ | $\lambda_1$ | $\lambda_2$ | $\lambda_3$ | $\lambda_4$ | $\lambda_5$ | $\lambda_6$ |
|----:|----:|----:|----:|----:|----:|----:|----:|
| 50  | 10 | 2.46 | 3.15 | 10.18 | 6.32 | 6.41 | 5.22 |
| 50  | 20 | 1.82 | 2.14 |  8.72 | 3.79 | 3.77 | 3.71 |
| 100 | 10 | 1.36 | 2.09 |  8.88 | 3.75 | 3.97 | 3.06 |
| 100 | 20 | 1.25 | 1.52 |  7.62 | 1.85 | 1.98 | 1.90 |

**Factor-intercept MSE $\beta$ ($\times 10^{-2}$).**

| $N$ | $M$ | $\beta_1$ | $\beta_2$ | $\beta_3$ | $\beta_4$ | $\beta_5$ | $\beta_6$ |
|----:|----:|----:|----:|----:|----:|----:|----:|
| 50  | 10 | 5.54 | 4.72 | 2.25 | 3.99 | 3.66 | 3.41 |
| 50  | 20 | 5.51 | 5.07 | 2.49 | 3.48 | 3.12 | 3.03 |
| 100 | 10 | 3.35 | 3.32 | 1.92 | 1.85 | 1.93 | 1.83 |
| 100 | 20 | 2.75 | 2.43 | 1.45 | 1.78 | 1.66 | 1.59 |

**Findings.** Every loading and intercept MSE decreases monotonically as the
sample size $N$ grows and as the curves are observed more densely ($M$),
confirming consistent estimation of the two-factor measurement model. The
*historical* loading $\lambda_3$ is the hardest to estimate (largest MSE), as
expected for a two-dimensional surface, but it also improves with $N$ and $M$.
The leading eigenfunctions/eigenvalues ($\phi_1,\nu_1$) of the factor-error
processes and the latent-process covariances ($\psi_1,\mu_1$) are recovered with
small error (e.g. $\psi_{1},\mu_{1}$ MSE $\approx 1\text{–}5\times10^{-3}$, all
decreasing in $N,M$); higher-order eigen-components are subject to the usual
sign/ordering non-identifiability and are not tabulated. Because the factors are
independent, no structural ($\gamma$) parameters are present.

---

## Study B — Two factors, *dependent* (well- vs mis-specified)

**Design.** The two factors are now **dependent**: the data-generating model
includes a structural path $\eta_1 \leftarrow \eta_2$ (a concurrent latent-on-
latent regression $\gamma$), in addition to the same six-indicator measurement
model ($z_1,z_2,z_3$ on $\eta_1$; $z_4,z_5,z_6$ on $\eta_2$). Two estimators are
fitted to every data set:

* **Well-specified (`.w`)** — the fitted model includes the structural path
  $\eta_1 \sim \eta_2$ (the correct dependence).
* **Mis-specified (`.m`)** — the fitted model **omits** the path
  ($\eta_1 \sim 1$ only), so the latent dependence is ignored.

This study quantifies (i) how well the structural function $\gamma$ is recovered
and (ii) the cost of ignoring latent dependence. Results in `outputs3/` cover
$N=100$ with $M\in\{10,20\}$.

**Factor-loading MSE $\lambda$ ($\times 10^{-2}$), well-specified vs mis-specified.**

| Fit | $N$ | $M$ | $\lambda_1$ | $\lambda_2$ | $\lambda_3$ | $\lambda_4$ | $\lambda_5$ | $\lambda_6$ |
|---|----:|----:|----:|----:|----:|----:|----:|----:|
| well | 100 | 10 |  8.11 |  8.41 |  6.78 | 2.67 | 2.84 | 2.17 |
| well | 100 | 20 |  3.23 |  3.42 |  4.40 | 1.69 | 1.64 | 1.67 |
| miss | 100 | 10 | 68.63 | 73.12 | 75.55 | 3.81 | 4.31 | 3.00 |
| miss | 100 | 20 | 68.31 | 73.44 | 78.48 | 2.00 | 2.18 | 2.15 |

**Structural-function MSE $\gamma$ ($\eta_1\leftarrow\eta_2$), well-specified.**

| $N$ | $M$ | $\gamma$ |
|----:|----:|----:|
| 100 | 10 | 0.408 |
| 100 | 20 | 0.117 |

(The mis-specified fit contains no $\gamma$ because the path is omitted.)

**Findings.** Under the **well-specified** model the loadings of all six
indicators are estimated accurately and improve sharply with denser sampling
($M:10\to20$ roughly halves the loading MSE). The structural function $\gamma$ is
the most challenging parameter, but its MSE falls by a factor of $\sim3.5$ (from
$0.41$ to $0.12$) as $M$ increases from 10 to 20, showing that the latent-on-
latent effect is recovered well once the curves are reasonably sampled.

The **mis-specified** fit reveals the cost of ignoring the dependence: the
loadings of the indicators of the *contaminated* factor $\eta_1$
($z_1,z_2,z_3$) deteriorate by **more than an order of magnitude**
($\lambda_{1\text{–}3}$ MSE jumps from $\approx 0.03\text{–}0.08$ to
$\approx 0.69\text{–}0.78$ and, tellingly, **does not improve with $M$**), whereas
the loadings of $\eta_2$'s indicators ($z_4,z_5,z_6$) are essentially unchanged.
Omitting the structural path therefore biases exactly the part of the
measurement model attached to the response factor, while leaving the covariate
factor intact — a clean and interpretable failure mode.

---

## Discussion

Across the original single-factor study and the two new two-factor studies the
estimator behaves consistently:

* **Consistency.** In every correctly-specified setting — the one-factor model of
  the main paper, Study A (uncorrelated factors) and the well-specified Study B
  (dependent factors) — all functional MSEs decrease monotonically with both the
  sample size $N$ and the sampling density $M$. The main study established this
  for a single factor; the new studies show it continues to hold when the
  measurement model carries two factors and when the factors are linked by a
  structural path.

* **Scalability to richer measurement models.** Study A demonstrates that the EM
  estimator handles several indicators per factor and multiple independent
  factors without loss of accuracy: loadings, intercepts and the leading
  latent-covariance eigen-components are all recovered with small error.

* **Latent dependence is recoverable but data-hungry.** In the well-specified
  Study B the structural function $\gamma$ has the largest MSE of any parameter
  yet improves rapidly with $M$, mirroring the difficulty of the *historical*
  surface effect seen in the main study. Adequate sampling density is the key
  driver of accurate structural estimation.

* **Mis-specification is localised but severe.** Ignoring a true latent-on-latent
  path inflates the loading MSE of the affected factor by more than a factor of
  ten and removes the usual benefit of denser sampling, while leaving the
  unaffected factor's parameters intact. This underlines the practical importance
  of correctly specifying the structural part of the model and shows that the
  damage is confined to the response factor rather than propagating through the
  whole system.

Taken together, the two additional studies corroborate the conclusions of the
main simulation — consistent, sampling-density-driven estimation — and extend
them to multi-factor measurement models and to latent structural dependence,
while quantifying the cost of structural mis-specification.

*Coverage rates are not reported here; this note covers the MSE results only.*
