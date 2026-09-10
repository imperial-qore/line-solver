"""Inverse Laplace Transform via CME (Concentrated Matrix Exponential).

This package provides numerical inverse Laplace transform using three methods:
  - CME/Talbot contour (default) with pre-computed optimal parameters
  - Euler summation
  - Gaver-Stehfest

Ported from the MATLAB iltcme library.
"""

from .matlab_ilt import matlab_ilt, matlab_ilt_matrix, cme_parameters

__all__ = ["matlab_ilt", "matlab_ilt_matrix", "cme_parameters"]
