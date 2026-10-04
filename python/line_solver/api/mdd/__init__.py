"""
Decision-diagram state-space storage and approximate aggregation.

Multi-valued decision diagram (MDD) storage of a CTMC reachable set, and the
Miner-Ciardo-Donatelli level aggregation that solves the chain without ever
forming the |S|-state generator.

Key algorithms:
    MDD: quasi-reduced ordered multi-valued decision diagram
    mdd_reachset: reachability set generation into an MDD
    mdd_closedqn: exact solve of a closed exponential QN over the diagram
    mdd_mcd: approximate stationary measures by level aggregation
    mdd_rec: exact normalising constant of a product form over the diagram
    mdd_descriptor: Kronecker descriptor, count-plus-phase (non-preemptive)
    mdd_ps: Kronecker descriptor, per-phase counts (shared servers)

References:
    A.S. Miner, G. Ciardo, "Efficient Reachability Set Generation and Storage
    Using Decision Diagrams", ICATPN 1999, LNCS 1639, pp.6-25.
    A.S. Miner, G. Ciardo, S. Donatelli, "Using the exact state space of a
    Markov model to compute approximate stationary measures", SIGMETRICS 2000.
    S. Balsamo, A. Marin, I. Stojic, "Computation of the normalising constant
    for product-form models of distributed systems with synchronisation",
    Future Generation Computer Systems 111 (2020) 475-490.
"""

from .mdd import MDD, MDDStruct, TERM_FALSE, TERM_TRUE
from .reachset import mdd_reachset
from .closedqn import mdd_closedqn
from .mcd import mdd_mcd
from .rec import mdd_rec, mdd_rec_masked, mdd_rec_marginal
from .descriptor import mdd_descriptor
from .ps import mdd_ps

__all__ = [
    'MDD',
    'MDDStruct',
    'TERM_TRUE',
    'TERM_FALSE',
    'mdd_reachset',
    'mdd_closedqn',
    'mdd_mcd',
    'mdd_rec',
    'mdd_rec_masked',
    'mdd_rec_marginal',
    'mdd_descriptor',
    'mdd_ps',
]
