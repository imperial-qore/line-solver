# iltcme - Inverse Laplace Transform via CME

This is a derivative work ported from MATLAB to Python.

The original MATLAB implementation provides numerical inverse Laplace transform
using the Abate-Whitt framework with three method variants:

- **CME (Concentrated Matrix Exponential / Talbot contour)**: Uses
  pre-computed optimal parameters from `iltcme.json` to place quadrature
  nodes along a Talbot-style contour in the complex plane. The parameters
  were obtained by optimizing the contour shape to minimise the coefficient
  of variation of the approximation error for a given number of function
  evaluations.

- **Euler**: Binomial-coefficient weighted Euler summation acceleration of
  the Bromwich inversion integral.

- **Gaver-Stehfest**: Factorial-based weights with real logarithmic
  abscissae.

## Pre-computed Parameters

The file `iltcme.json` contains 332 pre-computed optimal Talbot contour
parameter sets for `n` ranging from 1 to 1000. Each entry stores the
coefficients `a`, `b`, `c`, `mu1`, `omega` that define the contour shape,
along with `cv2` (squared coefficient of variation of the approximation
error) used for parameter set selection.

## Usage

```python
from line_solver.lib.thirdparty.iltcme import matlab_ilt
import numpy as np

# Inverse Laplace transform of 1/s^2 (expected: f(t) = t)
T = np.array([1.0, 2.0, 3.0])
result = matlab_ilt(lambda s: 1.0 / s**2, T, maxFnEvals=32)
```
