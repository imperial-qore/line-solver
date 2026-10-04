"""
MAP distance measures for continuous and discrete-time MAPs.

Continuous-time distances based on:
    G. Horvath, "Measuring the distance between MAPs and some
    applications," in Proc. ASMTA 2015, LNCS 9081, pp. 95-109.
    https://link.springer.com/chapter/10.1007/978-3-319-18579-8_8

Discrete-time extensions by QORE Lab (https://qore.doc.ic.ac.uk/)
"""

from .continuous import (
    map_exp_mul_int,
    map_dist,
    map_dist_lag1,
    map_geo_mul_sum,
    map_dist_acf,
    map_optim_dist,
    map_optim_dist_acf,
)

from .discrete import (
    dmap_geo_mul_sum,
    dmap_dist,
    dmap_dist_lag1,
    dmap_geo_mul_sum_acf,
    dmap_dist_acf,
    dmap_optim_dist,
    dmap_optim_dist_acf,
)
