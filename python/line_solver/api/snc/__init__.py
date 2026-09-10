"""Native Python implementations for stochastic network calculus (SNC).

The only api family in LINE whose deliverable is a TAIL QUANTILE with a
certified violation probability rather than a mean. MGF envelopes describe the
arrivals and the service, the min-plus operations compose them across a
feed-forward network, and the bounds are Chernoff infima over the free
parameter theta.

Key algorithms:
    arrival envelopes: snc_env_poisson, snc_env_cpoisson, snc_env_map,
        snc_env_tokenbucket
    service envelopes: snc_srv_rate (constant rate), snc_srv_exp (Exp server,
        JOB units -- the one that composes across hops)
    min-plus calculus: snc_leftover, snc_conv, snc_output
    tail bounds: snc_bound_backlog, snc_bound_delay
    quantiles: snc_perc_backlog, snc_perc_delay
    mean bounds: snc_mean_delay, snc_mean_backlog
    Chernoff search: snc_thetaopt

Wired into SolverBA as the method 'snc.upper'.
"""

from .envelopes import (
    snc_env_poisson, snc_env_cpoisson, snc_env_tokenbucket, snc_env_map,
    snc_srv_rate, snc_srv_exp,
    snc_leftover, snc_conv, snc_output,
)
from .bounds import (
    snc_thetaopt,
    snc_bound_backlog, snc_bound_delay,
    snc_perc_backlog, snc_perc_delay,
    snc_mean_delay, snc_mean_backlog,
)

__all__ = [
    'snc_env_poisson', 'snc_env_cpoisson', 'snc_env_tokenbucket', 'snc_env_map',
    'snc_srv_rate', 'snc_srv_exp',
    'snc_leftover', 'snc_conv', 'snc_output',
    'snc_thetaopt',
    'snc_bound_backlog', 'snc_bound_delay',
    'snc_perc_backlog', 'snc_perc_delay',
    'snc_mean_delay', 'snc_mean_backlog',
]
