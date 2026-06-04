# GPU EM estimation: `em.estimation.optimised.gpu`

`Core/FSEM_gpu.R` adds a GPU (CUDA / `torch`) re-implementation of
`em.estimation.optimised.revised`, named **`em.estimation.optimised.gpu`**.
It targets the RTX 5090 (Blackwell, compute capability `sm_120`) but runs on
any CUDA GPU and falls back to CPU `torch` tensors if CUDA is unavailable.

## Same interface as the CPU `revised` estimator

```r
em.estimation.optimised.gpu(model, data, x.data = NULL,
                            initial.parameter, n.b,
                            n.em = 100, n.monte = 100,
                            s.p = c("min", "sd"),
                            range.min = NULL, range.max = NULL,
                            design = c("regular","irregular",
                                       "regular.truncated","regular.missing"),
                            plot.progress = FALSE, save.every = NULL,
                            gamma.init.sd = 0.3, seed = NULL,
                            gpu.device = c("cuda","cpu"), gpu.verbose = TRUE)
```

* All arguments of `em.estimation.optimised.revised` are accepted verbatim,
  including the `gamma.init.sd` / `seed` jitter that fixes the coupled-latent
  EM degeneracy (see `Core/FSEM_optimised_revised.R`).
* Two extra optional arguments: `gpu.device` (`"cuda"` default, or `"cpu"` to
  force CPU tensors) and `gpu.verbose`.

## What runs where

* **GPU (resident across all inner iterations).** The compute-dominant M-step:
  penalised weighted least squares (IRLS), the 4-point `delta` cross
  validation, the `sigma` updates, and the CV likelihood/bound arithmetic —
  for every factor model and every structural regression. The stacked
  `(monte, sample, basis, p)` design tensor is uploaded **once** per block;
  the block forms `t(F) kron(I, Sigma^-1) F` and `t(F) kron(I, Sigma^-1) z`
  are single batched `torch_einsum` reductions; `linalg_solve` / `linalg_inv`
  and the `sigma` einsum stay on the device. No host↔device transfer happens
  inside the IRLS loop. Computation is in `float64` for fidelity.
* **CPU (reused from the protected helpers).** The E-step posterior moments
  and Monte-Carlo latent draws (`eta.distribution`, `params.z1`, `params.z`,
  `generator.z`, `MASS::mvrnorm` rejection sampling), the model-coupled design
  assembly, and the final standardization / convergence pipeline
  (`ex.estimations`, `parameter.estimated.evaluted`). These are model-structure
  bookkeeping, not numeric hot spots, and live in `Core/FSEM.R`, which is left
  untouched.

A literal "100% on GPU" port is not possible without rewriting the protected
`Core/FSEM.R` helpers; everything that is arithmetic-heavy is on the GPU.

## Fidelity

Results are statistically equivalent to `em.estimation.optimised.revised` but
**not bit-identical** (GPU reduction order and CPU RNG differ). Use the same
`seed` to reproduce a GPU run.

Validation (`verify_gpu.R`, N=50, M=10, n.b=6, n.em=8, n.monte=20, seed=1):

| coefficient | corr(CPU, GPU) | RMSE ratio |
|---|---|---|
| loadings z1–z6 ← eta | +0.994 … +0.9998 | 0.97 … 1.02 |
| structural eta1 ← eta2 | +0.998 | 1.07 |

GPU time 119 s vs CPU 154 s on this tiny problem; the GPU advantage grows with
`n.b`, `n.monte`, and `N`.

## One-time backend install (RTX 50-series / Blackwell)

The R `torch` package auto-detects the system CUDA, but the system CUDA may be
newer than `torch` supports (e.g. 13.x), and `sm_120` requires CUDA ≥ 12.8.
Pin CUDA 12.8 and raise the download timeout (the backend is ~3.4 GB, and R's
default 60 s `download.file` timeout aborts it):

```r
options(timeout = 7200)
Sys.setenv(CUDA = "12.8")
torch::install_torch(type = "cuda")
```

Confirm a real CUDA kernel executes (not just `cuda_is_available()`):

```r
Sys.setenv(CUDA = "12.8"); library(torch)
x <- torch_randn(2000, 2000, device = "cuda", dtype = torch_float64())
as.numeric(torch_mm(x, x)[1, 1]$cpu())        # finite number => sm_120 kernels run
```

## Run

```bash
"C:/Program Files/R/R-4.5.1/bin/Rscript.exe" run_gpu_example.R   # end-to-end demo
"C:/Program Files/R/R-4.5.1/bin/Rscript.exe" verify_gpu.R        # CPU vs GPU check
```

Source order (required): `Core/FSEM.R` → `Core/simualtion_utilities.R` →
`Core/FSEM_optimised.R` → (`Core/fsem_kernels_fix.R`) →
`Core/FSEM_optimised_revised.R` → `Core/FSEM_gpu.R`.
