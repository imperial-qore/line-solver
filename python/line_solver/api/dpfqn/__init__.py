"""
Discrete-time (slotted) product-form queueing network algorithms.

Native Python counterparts of :mod:`line_solver.api.pfqn` for models that live
on a lattice of unit slots rather than on the continuous time axis. In each
slot a busy server completes with probability p and an arrival occurs with
probability b, both recorded at the end of the slot with the departure
resolved before the arrival (Daduna's LA rule and D/A rule).

Key algorithms:
    dpfqn_nc: normalizing constants of a closed cycle of Bernoulli servers
        with state independent service probabilities (Buzen-style recursions)
    dpfqn_ncld: the same constants with state dependent service probabilities,
        by truncated convolution, plus the complement and arrival weights

Reference: H. Daduna, Queueing Networks with Discrete Time Scale, LNCS 2046,
Springer, 2001.
"""

from .nc import dpfqn_nc, dpfqn_ncld

__all__ = ['dpfqn_nc', 'dpfqn_ncld']
