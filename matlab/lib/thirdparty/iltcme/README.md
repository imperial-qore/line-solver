# iltcme - Inverse Laplace Transform via CME

Pre-computed Talbot-contour parameters for the Concentrated Matrix Exponential
(CME) method of numerical Laplace transform inversion, in the Abate-Whitt
framework.

`iltcme.json` contains 332 optimal parameter sets for `n` from 1 to 1000. Each
entry stores the coefficients `a`, `b`, `c`, `mu1` and `omega` that define the
contour shape, plus `cv2` (the squared coefficient of variation of the
approximation error) used to select a parameter set.

## Usage from LINE

`CME.m` reads this table through `CME.table()`, which caches the decoded struct
after the first call:

```matlab
d = CME(1.0, 7);   % order-7 concentrated matrix exponential, mean 1
```

`lineStart` puts this directory on the MATLAB path via `addpath(genpath(...))`,
so `fileread('iltcme.json')` in `CME.table` resolves without an explicit path.

## Provenance

The same table is vendored for the other codebases at
`jar/src/main/resources/iltcme.json` (read as a classpath resource) and
`python/line_solver/api/lti/iltcme.json`. All copies are byte-identical; keep
them so.

Upstream: Horvath, Horvath, Almousa and Telek, "Numerical inverse Laplace
transformation using concentrated matrix exponential distributions",
Performance Evaluation 137 (2020), https://inverselaplace.org/
