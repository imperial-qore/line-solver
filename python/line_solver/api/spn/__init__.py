"""
Stochastic Petri net analysis.

Algorithms specific to the SPN formalism, as opposed to the formalism-agnostic
decision-diagram machinery in api/mdd.

Key algorithms:
    spn_mdd: decision-diagram reachable set and Kronecker rate descriptor of a
        stochastic Petri net, consumed by mdd_mcd
    spn_pf: decide the product form of a net and derive the per-level factors
        g_l that mdd_rec and spn_metrics take as input
    spn_rec_enabled: enabling-degree distribution of one mode, by masked MDD-rec
    spn_metrics: stationary token, utilization and throughput measures
    spn_sinvariants: minimal-support S-invariants and the load vector S m0
    spn_lpbnd: moment-relaxation linear program bounding the mean marking and
        the throughputs, without building the reachable set
    spn_conv: normalising constant by convolution over the load vector

References:
    A.S. Miner, G. Ciardo, "Efficient Reachability Set Generation and Storage
    Using Decision Diagrams", ICATPN 1999, LNCS 1639, pp.6-25.
    S. Balsamo, A. Marin, I. Stojic, "Computation of the normalising constant
    for product-form models of distributed systems with synchronisation",
    Future Generation Computer Systems 111 (2020) 475-490.
    J. Coleman, W. Henderson, P. Taylor, "Product form equilibrium
    distributions and a convolution algorithm for stochastic Petri nets",
    Performance Evaluation 26(3), 1996, 159-180.
    Z. Liu, "Performance Analysis of Stochastic Timed Petri Nets Using Linear
    Programming Approach", IEEE Trans. Software Engineering 24(11), 1998,
    1014-1030.
"""

from .mdd import spn_mdd
from .pf import spn_pf
from .rec_enabled import spn_rec_enabled
from .metrics import spn_metrics
from .sinvariants import spn_sinvariants
from .lpbnd import spn_lpbnd
from .conv import spn_conv

__all__ = [
    'spn_mdd',
    'spn_pf',
    'spn_rec_enabled',
    'spn_metrics',
    'spn_sinvariants',
    'spn_lpbnd',
    'spn_conv',
]
