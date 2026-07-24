"""
Moment Conversion Algorithms for Discrete Distributions.

Native Python implementations of the one-to-one conversions between the power
(raw), factorial, upward-factorial, binomial, negative-binomial and central
moments of a discrete random variable, together with the combinatorial
triangles that underpin them.

Key algorithms:
    moment_stirling1: Signed Stirling numbers of the first kind
    moment_stirling2: Stirling numbers of the second kind
    moment_stirlingcycle: Stirling cycle numbers
    moment_lah: Lah numbers
    moment_binotrans: Binomial transform of a sequence
    moment_factorial_from_raw: Raw to factorial moments
    moment_upfactorial_from_raw: Raw to upward-factorial moments
    moment_central_from_raw: Raw to central moments
    moment_cumulant_from_raw: Raw moments to cumulants
    moment_factcumulant_from_factorial: Factorial moments to factorial cumulants
    moment_joint_*: the same house for the joint moments of a random vector
    moment_joint_marking: per-class counts of a multinomially marked count
    moment_joint_aggregate: factorial moments of a sum from the joint ones
"""

from .moment import (
    moment_stirling1,
    moment_stirling2,
    moment_stirlingcycle,
    moment_lah,
    moment_binotrans,
    moment_binotransinv,
    moment_factorial_from_raw,
    moment_raw_from_factorial,
    moment_upfactorial_from_raw,
    moment_raw_from_upfactorial,
    moment_binomial_from_factorial,
    moment_factorial_from_binomial,
    moment_negbinomial_from_upfactorial,
    moment_upfactorial_from_negbinomial,
    moment_binomial_from_negbinomial,
    moment_negbinomial_from_binomial,
    moment_factorial_from_upfactorial,
    moment_upfactorial_from_factorial,
    moment_central_from_raw,
    moment_raw_from_central,
    moment_cumulant_from_raw,
    moment_raw_from_cumulant,
    moment_factcumulant_from_factorial,
    moment_factorial_from_factcumulant,
    moment_tensortrans,
    moment_housematrix,
    moment_jointtrans,
    moment_joint_factorial_from_raw,
    moment_joint_raw_from_factorial,
    moment_joint_upfactorial_from_raw,
    moment_joint_raw_from_upfactorial,
    moment_joint_binomial_from_factorial,
    moment_joint_factorial_from_binomial,
    moment_joint_negbinomial_from_upfactorial,
    moment_joint_upfactorial_from_negbinomial,
    moment_joint_factorial_from_upfactorial,
    moment_joint_upfactorial_from_factorial,
    moment_joint_negbinomial_from_binomial,
    moment_joint_binomial_from_negbinomial,
    moment_joint_central_from_raw,
    moment_joint_central_from_raw_mean,
    moment_joint_raw_from_central,
    moment_joint_cumulant_from_raw,
    moment_joint_raw_from_cumulant,
    moment_joint_factcumulant_from_factorial,
    moment_joint_factorial_from_factcumulant,
    moment_binomial_from_tail,
    moment_tail_from_binomial,
    moment_joint_binomial_from_tail,
    moment_joint_tail_from_binomial,
    moment_joint_central_from_tail,
    moment_joint_marking,
    moment_joint_aggregate,
)

__all__ = [
    'moment_stirling1',
    'moment_stirling2',
    'moment_stirlingcycle',
    'moment_lah',
    'moment_binotrans',
    'moment_binotransinv',
    'moment_factorial_from_raw',
    'moment_raw_from_factorial',
    'moment_upfactorial_from_raw',
    'moment_raw_from_upfactorial',
    'moment_binomial_from_factorial',
    'moment_factorial_from_binomial',
    'moment_negbinomial_from_upfactorial',
    'moment_upfactorial_from_negbinomial',
    'moment_binomial_from_negbinomial',
    'moment_negbinomial_from_binomial',
    'moment_factorial_from_upfactorial',
    'moment_upfactorial_from_factorial',
    'moment_central_from_raw',
    'moment_raw_from_central',
    'moment_cumulant_from_raw',
    'moment_raw_from_cumulant',
    'moment_factcumulant_from_factorial',
    'moment_factorial_from_factcumulant',
    'moment_tensortrans',
    'moment_housematrix',
    'moment_jointtrans',
    'moment_joint_factorial_from_raw',
    'moment_joint_raw_from_factorial',
    'moment_joint_upfactorial_from_raw',
    'moment_joint_raw_from_upfactorial',
    'moment_joint_binomial_from_factorial',
    'moment_joint_factorial_from_binomial',
    'moment_joint_negbinomial_from_upfactorial',
    'moment_joint_upfactorial_from_negbinomial',
    'moment_joint_factorial_from_upfactorial',
    'moment_joint_upfactorial_from_factorial',
    'moment_joint_negbinomial_from_binomial',
    'moment_joint_binomial_from_negbinomial',
    'moment_joint_central_from_raw',
    'moment_joint_central_from_raw_mean',
    'moment_joint_raw_from_central',
    'moment_joint_cumulant_from_raw',
    'moment_joint_raw_from_cumulant',
    'moment_joint_factcumulant_from_factorial',
    'moment_joint_factorial_from_factcumulant',
    'moment_binomial_from_tail',
    'moment_tail_from_binomial',
    'moment_joint_binomial_from_tail',
    'moment_joint_tail_from_binomial',
    'moment_joint_central_from_tail',
    'moment_joint_marking',
    'moment_joint_aggregate',
]
