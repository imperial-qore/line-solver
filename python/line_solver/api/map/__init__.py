"""
MAP Queue Analysis Algorithms.

Native Python implementations for analyzing MAP-driven queues
including MAP/M/1-PS response time distribution.

Key algorithms:
    map_m1ps_cdf_respt: Complementary CDF of sojourn time in MAP/M/1-PS queue
    map_m1ps_sojourn: Alias for map_m1ps_cdf_respt
    map_compute_R: Compute rate matrix R for MAP/M/1 queue
    map_m1ps_h_recursive: Recursive h_{n,k} coefficients

References:
    Masuyama, H., & Takine, T. (2003). Sojourn time distribution in a
    MAP/M/1 processor-sharing queue. Operations Research Letters, 31(6), 406-412.
"""

from ..mam.mapm1ps import (
    map_m1ps_cdf_respt,
    map_compute_R,
    map_m1ps_h_recursive,
    map_m1ps_sojourn,
)

__all__ = [
    'map_m1ps_cdf_respt',
    'map_compute_R',
    'map_m1ps_h_recursive',
    'map_m1ps_sojourn',
]
