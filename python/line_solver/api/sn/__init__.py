"""
Service Network (SN) utilities.

Native Python implementations for stochastic network structure
analysis, validation, and parameter extraction.

Key classes:
    NetworkStruct: Data structure summarizing network characteristics
    SnGetDemandsResult: Result of sn_get_demands_chain calculation

Key functions:
    sn_get_demands_chain: Aggregate class-level parameters into chain-level
    sn_has_*: Network property predicates
    sn_is_*: Model type checks
    sn_get_*: Parameter extraction
    sn_validate: Network validation
"""

from .network_struct import (
    MatrixArray,
    NetworkStruct,
    NodeType,
    SchedStrategy,
    RoutingStrategy,
    DropStrategy,
)

from .demands import (
    sn_get_demands_chain,
    SnGetDemandsResult,
)

from .deaggregate import (
    sn_deaggregate_chain_results,
    SnDeaggregateResult,
)

from .transforms import (
    ProductFormParams,
    sn_get_product_form_params,
    sn_get_residt_from_respt,
    sn_get_state_aggr,
    sn_set_arrival,
    sn_set_service,
    sn_set_servers,
    sn_set_population,
    sn_set_priority,
    sn_set_routing,
    sn_refresh_visits,
    sn_set_fork_fanout,
    sn_set_service_batch,
    sn_nonmarkov_toph,
)

from .getters import (
    ChainParams,
    sn_get_arvr_from_tput,
    sn_get_node_arvr_from_tput,
    sn_get_node_tput_from_tput,
    sn_get_product_form_chain_params,
    sn_set_routing_prob,
)

from .utils import (
    sn_print,
    sn_print_routing_matrix,
    sn_refresh_process_fields,
    sn_rtnodes_to_rtorig,
)

from .predicates import (
    # Model type predicates
    sn_is_closed_model,
    sn_is_open_model,
    sn_is_mixed_model,
    sn_is_population_model,
    # Class predicates
    sn_has_closed_classes,
    sn_has_open_classes,
    sn_has_mixed_classes,
    sn_has_single_class,
    sn_has_multi_class,
    sn_has_multiple_closed_classes,
    # Chain predicates
    sn_has_single_chain,
    sn_has_multi_chain,
    # Scheduling predicates
    sn_has_fcfs,
    sn_has_ps,
    sn_has_inf,
    sn_has_lcfs,
    sn_has_lcfspr,
    sn_has_lcfs_pi,
    sn_has_siro,
    sn_has_dps,
    sn_has_dps_prio,
    sn_has_gps,
    sn_has_gps_prio,
    sn_has_ps_prio,
    sn_has_hol,
    sn_has_lps,
    sn_has_setf,
    sn_has_sept,
    sn_has_lept,
    sn_has_sjf,
    sn_has_ljf,
    sn_has_polling,
    sn_has_homogeneous_scheduling,
    sn_has_multi_class_fcfs,
    sn_has_multi_class_heter_fcfs,
    sn_has_multi_class_heter_exp_fcfs,
    # Server predicates
    sn_has_multi_server,
    # Load dependence predicates
    sn_has_load_dependence,
    sn_has_joint_dependence,
    # Structure predicates
    sn_has_fork_join,
    sn_has_priorities,
    sn_has_class_switching,
    sn_has_fractional_populations,
    # Product form predicates
    sn_has_sd_routing,
    sn_has_product_form,
    sn_has_product_form_not_het_fcfs,
    sn_has_product_form_except_multi_class_heter_exp_fcfs,
    # State predicates
    sn_is_state_valid,
)

__all__ = [
    # Core classes
    'MatrixArray',
    'NetworkStruct',
    'NodeType',
    'SchedStrategy',
    'RoutingStrategy',
    'DropStrategy',
    # Chain functions
    'sn_get_demands_chain',
    'SnGetDemandsResult',
    'sn_deaggregate_chain_results',
    'SnDeaggregateResult',
    # Transform functions
    'ProductFormParams',
    'sn_get_product_form_params',
    'sn_get_residt_from_respt',
    'sn_get_state_aggr',
    'sn_set_arrival',
    'sn_set_service',
    'sn_set_servers',
    'sn_set_population',
    'sn_set_priority',
    'sn_set_routing',
    'sn_refresh_visits',
    'sn_set_fork_fanout',
    'sn_set_service_batch',
    'sn_nonmarkov_toph',
    # Getter functions
    'ChainParams',
    'sn_get_arvr_from_tput',
    'sn_get_node_arvr_from_tput',
    'sn_get_node_tput_from_tput',
    'sn_get_product_form_chain_params',
    'sn_set_routing_prob',
    # Model type predicates
    'sn_is_closed_model',
    'sn_is_open_model',
    'sn_is_mixed_model',
    'sn_is_population_model',
    # Class predicates
    'sn_has_closed_classes',
    'sn_has_open_classes',
    'sn_has_mixed_classes',
    'sn_has_single_class',
    'sn_has_multi_class',
    'sn_has_multiple_closed_classes',
    # Chain predicates
    'sn_has_single_chain',
    'sn_has_multi_chain',
    # Scheduling predicates
    'sn_has_fcfs',
    'sn_has_ps',
    'sn_has_inf',
    'sn_has_lcfs',
    'sn_has_lcfspr',
    'sn_has_lcfs_pi',
    'sn_has_siro',
    'sn_has_dps',
    'sn_has_dps_prio',
    'sn_has_gps',
    'sn_has_gps_prio',
    'sn_has_ps_prio',
    'sn_has_hol',
    'sn_has_lps',
    'sn_has_setf',
    'sn_has_sept',
    'sn_has_lept',
    'sn_has_sjf',
    'sn_has_ljf',
    'sn_has_polling',
    'sn_has_homogeneous_scheduling',
    'sn_has_multi_class_fcfs',
    'sn_has_multi_class_heter_fcfs',
    'sn_has_multi_class_heter_exp_fcfs',
    # Server predicates
    'sn_has_multi_server',
    # Load dependence predicates
    'sn_has_load_dependence',
    'sn_has_joint_dependence',
    # Structure predicates
    'sn_has_fork_join',
    'sn_has_priorities',
    'sn_has_class_switching',
    'sn_has_fractional_populations',
    # Product form predicates
    'sn_has_sd_routing',
    'sn_has_product_form',
    'sn_has_product_form_not_het_fcfs',
    'sn_has_product_form_except_multi_class_heter_exp_fcfs',
    # State predicates
    'sn_is_state_valid',
    # Utility functions
    'sn_print',
    'sn_print_routing_matrix',
    'sn_refresh_process_fields',
    'sn_rtnodes_to_rtorig',
]
