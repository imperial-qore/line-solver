"""
Utilities for SolverMAM.

Provides helper functions for network analysis, traffic computation, and
result building.
"""

from .network_adapter import (
    extract_mam_params,
    identify_station_types,
    build_routing_matrix,
    get_service_distribution,
    extract_visit_counts,
    check_singleclass_network,
    check_closed_network,
    check_product_form,
    is_fork_join_network,
)

__all__ = [
    'extract_mam_params',
    'identify_station_types',
    'build_routing_matrix',
    'get_service_distribution',
    'extract_visit_counts',
    'check_singleclass_network',
    'check_closed_network',
    'check_product_form',
    'is_fork_join_network',
]
