"""
Perm: Matrix Permanent Computation.

Native Python implementations of matrix permanent algorithms, the twin of the
JAR package jline.lib.perm and of MATLAB perm.m / perm_heur.m.

The permanent of a matrix is similar to the determinant but uses only additions
(no subtractions). Computing the permanent is #P-complete, so exact computation
is expensive for large matrices and the package offers three families:

Exact (exact submodule):
    Permanent: inclusion-exclusion exploiting repeated rows or columns (perm.m)
    RyzerPermanent: Ryser's formula, Gray-code and naive variants
    NaivePermanent: enumeration of all permutations, O(n!)

Deterministic approximations (approx submodule):
    BethePermanent: sum-product (belief propagation) Bethe approximation
    HeuristicPermanent: Sinkhorn scaling with mean-field and Gurvits bounds
    SaddlePointPermanent: Laplace expansion of the coefficient integral,
        the homogeneous variant of cache_spm (perm_spm.m)

Randomized approximations (sampling submodule):
    AdaPartSampler: adaptive partitioning rejection sampling on Soules bounds
    HuberLawSampler: Huber-Law acceptance-rejection importance sampling

Permanent-based network marginals (network submodule):
    preprocessing_ds, NetworkNoThink, NetworkThink

The samplers take a seed and draw from a numpy Generator, so their runs are
reproducible; the JAR twins use an unseeded java.util.Random and therefore
agree in distribution only.
"""

from .approx import (
    BethePermanent,
    HeuristicPermanent,
    SaddlePointPermanent,
    perm_bethe,
    perm_heur,
    perm_spm,
)
from .base import PermResult, PermSolver
from .exact import (
    NaivePermanent,
    Permanent,
    RyzerPermanent,
    compute_permanent,
    perm,
    permanent,
    _binomial_coefficient,
    _find_unique_columns_with_multiplicities,
    _permanent_naive,
    _permanent_ryser,
    _permanent_with_multiplicities,
    _pprod_next,
)
from .network import MarginalResult, NetworkNoThink, NetworkThink, preprocessing_ds
from .sampling import AdaPartSampler, HuberLawSampler

__all__ = [
    'compute_permanent',
    'permanent',
    'perm',
    'perm_heur',
    'perm_bethe',
    'perm_spm',
    'Permanent',
    'PermResult',
    'PermSolver',
    'NaivePermanent',
    'RyzerPermanent',
    'BethePermanent',
    'HeuristicPermanent',
    'SaddlePointPermanent',
    'AdaPartSampler',
    'HuberLawSampler',
    'preprocessing_ds',
    'NetworkNoThink',
    'NetworkThink',
    'MarginalResult',
]
