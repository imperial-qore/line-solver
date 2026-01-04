
"""
Stochastic network (SN) utility functions.

This module provides utility functions for analyzing stochastic networks,
including parameter extraction, model introspection, result processing,
and various transformations.

Key function categories:
- Model introspection: sn_has_* functions to check model properties
- Parameter extraction: sn_get_* functions for demands, visits, etc.
- Result processing: sn_deaggregate_chain_results  
- Model classification: sn_is_* functions (closed, open, mixed models)
- Utility functions: state validation, routing matrix operations

These functions support the internal workings of LINE solvers by providing
common operations on stochastic network models and results.
"""

import jpype
import numpy as np
from line_solver import jlineMatrixToArray, jlineMatrixFromArray


def sn_deaggregate_chain_results(sn, Lchain, ST, STchain, Vchain, alpha, Qchain, Uchain, Rchain, Tchain, Cchain, Xchain):
    """
    Deaggregate chain-level results back to individual node and class results.

    Takes aggregated performance metrics at the chain level and disaggregates them
    to provide detailed node-level and class-level performance metrics.

    Args:
        sn: The stochastic network model
        Lchain: Chain service demands
        ST: Station service times  
        STchain: Chain station service times
        Vchain: Chain visit ratios
        alpha: Chain weights/proportions
        Qchain: Chain queue lengths
        Uchain: Chain utilizations
        Rchain: Chain response times
        Tchain: Chain throughputs
        Cchain: Chain capacities
        Xchain: Chain rates

    Returns:
        dict: Dictionary containing disaggregated results with keys:
            - 'QN': Node queue lengths
            - 'UN': Node utilizations  
            - 'RN': Node response times
            - 'TN': Node throughputs
            - 'CN': Node capacities
            - 'XN': Node rates
            - 'lG': Log normalizing constant
    """
    ST_matrix = jlineMatrixFromArray(ST) if ST is not None else None
    Qchain_matrix = jlineMatrixFromArray(Qchain) if Qchain is not None else None
    Uchain_matrix = jlineMatrixFromArray(Uchain) if Uchain is not None else None
    Cchain_matrix = jlineMatrixFromArray(Cchain) if Cchain is not None else None

    result = jpype.JPackage('jline').api.sn.SnDeaggregateChainResultsKt.sn_deaggregate_chain_results(
        sn,
        jlineMatrixFromArray(Lchain),
        ST_matrix,
        jlineMatrixFromArray(STchain),
        jlineMatrixFromArray(Vchain),
        jlineMatrixFromArray(alpha),
        Qchain_matrix,
        Uchain_matrix,
        jlineMatrixFromArray(Rchain),
        jlineMatrixFromArray(Tchain),
        Cchain_matrix,
        jlineMatrixFromArray(Xchain)
    )

    return {
        'QN': jlineMatrixToArray(result.QN),
        'UN': jlineMatrixToArray(result.UN),
        'RN': jlineMatrixToArray(result.RN),
        'TN': jlineMatrixToArray(result.TN),
        'CN': jlineMatrixToArray(result.CN),
        'XN': jlineMatrixToArray(result.XN),
        'lG': result.lG
    }


def sn_get_arvr_from_tput(sn, TN=None, TH=None):
    """
    Calculate arrival rates from throughput values.

    Computes arrival rates based on given throughput measurements,
    considering the network topology and routing probabilities.

    Args:
        sn: The stochastic network model
        TN: Node throughputs (optional)
        TH: System throughput (optional)

    Returns:
        numpy.ndarray: Arrival rates for each node and class
    """
    TN_matrix = jlineMatrixFromArray(TN) if TN is not None else None

    result = jpype.JPackage('jline').api.sn.SnGetArvRFromTputKt.sn_get_arv_r_from_tput(
        sn, TN_matrix, TH
    )

    return jlineMatrixToArray(result)


def sn_get_demands_chain(sn):
    """
    Extract chain-level service demands and parameters.

    Aggregates individual class service demands into chain-level parameters
    for analysis in the chain-decomposed representation of the model.

    Args:
        sn: The stochastic network model

    Returns:
        dict: Dictionary containing chain parameters:
            - 'Lchain': Chain service demands
            - 'Nchain': Chain populations
            - 'Zchain': Chain think times
            - 'refstat': Reference station indices
            - 'alpha': Chain mixing probabilities
            - 'sn': Updated network model
    """
    result = jpype.JPackage('jline').api.sn.SnGetDemandsChainKt.sn_get_demands_chain(sn)

    return {
        'Lchain': jlineMatrixToArray(result.Lchain),
        'Nchain': jlineMatrixToArray(result.Nchain),
        'Zchain': jlineMatrixToArray(result.Zchain),
        'refstat': result.refstat,
        'alpha': jlineMatrixToArray(result.alpha),
        'sn': result.sn
    }


def sn_get_node_arvr_from_tput(sn, TN, TH=None, AN=None):
    """
    Calculate node arrival rates from node throughput values.

    Computes arrival rates at specific nodes based on measured throughput,
    accounting for network routing and feedback loops.

    Args:
        sn: The stochastic network model
        TN: Node throughputs
        TH: System throughput (optional)
        AN: Node arrivals (optional)

    Returns:
        numpy.ndarray: Node arrival rates
    """
    AN_matrix = jlineMatrixFromArray(AN) if AN is not None else None

    result = jpype.JPackage('jline').api.sn.SnGetNodeArvRFromTputKt.snGetNodeArvRFromTput(
        sn, jlineMatrixFromArray(TN), TH, AN_matrix
    )

    return jlineMatrixToArray(result)


def sn_get_node_tput_from_tput(sn, TN, TH=None, ANn=None):
    """
    Calculate node throughput from system or other node throughputs.

    Derives individual node throughput values from overall system throughput
    or other node throughput measurements using visit ratios.

    Args:
        sn: The stochastic network model
        TN: Node throughputs
        TH: System throughput (optional)
        ANn: Node-specific arrivals (optional)

    Returns:
        numpy.ndarray: Calculated node throughputs
    """
    ANn_matrix = jlineMatrixFromArray(ANn) if ANn is not None else None

    result = jpype.JPackage('jline').api.sn.SnGetNodeTputFromTputKt.snGetNodeTputFromTput(
        sn, jlineMatrixFromArray(TN), TH, ANn_matrix
    )

    return jlineMatrixToArray(result)


def sn_get_product_form_chain_params(sn):
    """
    Extract product-form chain parameters for MVA and related algorithms.

    Retrieves the parameters needed for product-form analysis at the chain level,
    including service demands, populations, and server characteristics.

    Args:
        sn: The stochastic network model

    Returns:
        dict: Chain-level product-form parameters:
            - 'L': Chain service demands
            - 'N': Chain populations
            - 'Z': Chain think times
            - 'mu': Chain service rates
            - 'phi': Chain service time distributions
            - 'nservers': Number of servers per station
            - 'schedid': Scheduling discipline identifiers
            - 'refstat': Reference station index
            - 'sn': Updated network model
    """
    result = jpype.JPackage('jline').api.sn.SnGetProductFormChainParamsKt.snGetProductFormChainParams(sn)

    return {
        'L': jlineMatrixToArray(result.L),
        'N': jlineMatrixToArray(result.N),
        'Z': jlineMatrixToArray(result.Z),
        'mu': jlineMatrixToArray(result.mu),
        'phi': jlineMatrixToArray(result.phi),
        'nservers': jlineMatrixToArray(result.nservers),
        'schedid': jlineMatrixToArray(result.schedid),
        'refstat': result.refstat,
        'sn': result.sn
    }


def sn_get_product_form_params(sn):
    """
    Extract product-form parameters for class-level analysis.

    Retrieves parameters needed for product-form queueing network analysis
    algorithms like MVA, AMVA, and Convolution algorithms.

    Args:
        sn: The stochastic network model

    Returns:
        dict: Product-form parameters:
            - 'L': Service demands matrix
            - 'N': Population vector
            - 'Z': Think times
            - 'mu': Service rates
            - 'phi': Service time coefficients
            - 'nservers': Number of servers per station
            - 'schedid': Scheduling discipline identifiers
            - 'refstat': Reference station index
            - 'sn': Updated network model
    """
    result = jpype.JPackage('jline').api.sn.SnGetProductFormParamsKt.snGetProductFormParams(sn)

    return {
        'L': jlineMatrixToArray(result.L),
        'N': jlineMatrixToArray(result.N),
        'Z': jlineMatrixToArray(result.Z),
        'mu': jlineMatrixToArray(result.mu),
        'phi': jlineMatrixToArray(result.phi),
        'nservers': jlineMatrixToArray(result.nservers),
        'schedid': jlineMatrixToArray(result.schedid),
        'refstat': result.refstat,
        'sn': result.sn
    }


def sn_get_residt_from_respt(sn, RNclass, WH=None):
    """
    Calculate residence times from response times.

    Computes node residence times based on measured response times,
    useful for performance analysis and bottleneck identification.

    Args:
        sn: The stochastic network model
        RNclass: Response times by class
        WH: Waiting times (optional)

    Returns:
        numpy.ndarray: Residence times matrix
    """
    result = jpype.JPackage('jline').api.sn.SnGetResidTFromRespTKt.snGetResidTFromRespT(
        sn, jlineMatrixFromArray(RNclass), WH
    )

    return jlineMatrixToArray(result)


def sn_get_state_aggr(sn):
    """
    Get state space aggregation information for the network.

    Returns aggregated state information for each node, which is useful
    for state-dependent analysis and Markov chain representations.

    Args:
        sn: The stochastic network model

    Returns:
        dict: State aggregation mapping for each node
    """
    result = jpype.JPackage('jline').api.sn.SnGetStateAggrKt.snGetStateAggr(sn)

    state_dict = {}
    for entry in result.entrySet():
        node_key = entry.getKey()
        state_matrix = jlineMatrixToArray(entry.getValue())
        state_dict[node_key] = state_matrix

    return state_dict


def sn_is_state_valid(sn):
    """
    Check if the current network state is valid.

    Validates that the current state configuration satisfies all
    network constraints and population requirements.

    Args:
        sn: The stochastic network model

    Returns:
        bool: True if the state is valid, False otherwise
    """
    return jpype.JPackage('jline').api.sn.SnIsStateValidKt.snIsStateValid(sn)


def sn_refresh_visits(sn, chains=None, rt=None, rtnodes=None):
    """
    Refresh visit ratios for the network model.

    Updates visit ratio calculations based on routing probabilities,
    ensuring consistency with the current network configuration.

    Args:
        sn: The stochastic network model
        chains: Chain specifications (optional)
        rt: Routing table (optional)
        rtnodes: Routing nodes (optional)

    Returns:
        The updated stochastic network model
    """
    chains_matrix = jlineMatrixFromArray(chains) if chains is not None else None
    rt_matrix = jlineMatrixFromArray(rt) if rt is not None else None
    rtnodes_matrix = jlineMatrixFromArray(rtnodes) if rtnodes is not None else None

    return jpype.JPackage('jline').api.sn.SnRefreshVisitsKt.snRefreshVisits(
        sn, chains_matrix, rt_matrix, rtnodes_matrix
    )


def sn_has_class_switching(sn):
    """
    Check if the network has class switching.

    Determines whether jobs can change class while moving through the network.

    Args:
        sn: The stochastic network model

    Returns:
        bool: True if the network has class switching, False otherwise
    """
    return jpype.JPackage('jline').api.sn.SnHasClassSwitchingKt.snHasClassSwitching(sn)


def sn_has_fork_join(sn):
    """
    Check if the network contains fork-join nodes.

    Identifies the presence of fork-join synchronization points in the network.

    Args:
        sn: The stochastic network model

    Returns:
        bool: True if fork-join nodes are present, False otherwise
    """
    return jpype.JPackage('jline').api.sn.SnHasForkJoinKt.snHasForkJoin(sn)


def sn_has_load_dependence(sn):
    """
    Check if the network has load-dependent service rates.

    Determines whether any stations have service rates that depend on
    the number of customers at the station.

    Args:
        sn: The stochastic network model

    Returns:
        bool: True if load-dependent service is present, False otherwise
    """
    return jpype.JPackage('jline').api.sn.SnHasLoadDependenceKt.snHasLoadDependence(sn)


def sn_has_multi_server(sn):
    """
    Check if the network contains multi-server stations.

    Identifies stations with more than one server.

    Args:
        sn: The stochastic network model

    Returns:
        bool: True if multi-server stations exist, False otherwise
    """
    return jpype.JPackage('jline').api.sn.SnHasMultiServerKt.snHasMultiServer(sn)


def sn_has_priorities(sn):
    """
    Check if the network uses priority scheduling.

    Determines whether any stations implement priority-based job scheduling.

    Args:
        sn: The stochastic network model

    Returns:
        bool: True if priority scheduling is used, False otherwise
    """
    return jpype.JPackage('jline').api.sn.SnHasPrioritiesKt.snHasPriorities(sn)


def sn_has_product_form(sn):
    """
    Check if the network has product-form solution.

    Determines whether the network satisfies the conditions for a
    product-form equilibrium distribution, enabling exact analysis.

    Args:
        sn: The stochastic network model

    Returns:
        bool: True if the network has product-form, False otherwise
    """
    return jpype.JPackage('jline').api.sn.SnHasProductFormKt.sn_has_product_form(sn)


def sn_has_closed_classes(sn):
    """
    Check if the network contains closed job classes.

    Identifies the presence of closed classes where jobs circulate
    indefinitely without external arrivals or departures.

    Args:
        sn: The stochastic network model

    Returns:
        bool: True if closed classes exist, False otherwise
    """
    return jpype.JPackage('jline').api.sn.SnHasClosedClassesKt.snHasClosedClasses(sn)


def sn_has_open_classes(sn):
    """
    Check if the network contains open job classes.

    Identifies the presence of open classes with external arrivals and departures.

    Args:
        sn: The stochastic network model

    Returns:
        bool: True if open classes exist, False otherwise
    """
    return jpype.JPackage('jline').api.sn.SnHasOpenClassesKt.snHasOpenClasses(sn)


def sn_has_mixed_classes(sn):
    """
    Check if the network has both open and closed classes.

    Determines whether the network is a mixed model with both
    open and closed job classes.

    Args:
        sn: The stochastic network model

    Returns:
        bool: True if both open and closed classes exist, False otherwise
    """
    return jpype.JPackage('jline').api.sn.SnHasMixedClassesKt.snHasMixedClasses(sn)


def sn_has_multi_chain(sn):
    """
    Check if the network has multiple chains.

    Determines whether the network has more than one routing chain.

    Args:
        sn: The stochastic network model

    Returns:
        bool: True if multiple chains exist, False otherwise
    """
    return jpype.JPackage('jline').api.sn.SnHasMultiChainKt.snHasMultiChain(sn)


def sn_is_closed_model(sn):
    """
    Check if the network is a closed queueing model.

    Determines whether all job classes are closed (no external arrivals/departures).

    Args:
        sn: The stochastic network model

    Returns:
        bool: True if the model is closed, False otherwise
    """
    return jpype.JPackage('jline').api.sn.SnIsClosedModelKt.snIsClosedModel(sn)


def sn_is_open_model(sn):
    """
    Check if the network is an open queueing model.

    Determines whether all job classes are open (with external arrivals/departures).

    Args:
        sn: The stochastic network model

    Returns:
        bool: True if the model is open, False otherwise
    """
    return jpype.JPackage('jline').api.sn.SnIsOpenModelKt.snIsOpenModel(sn)


def sn_is_mixed_model(sn):
    """
    Check if the network is a mixed queueing model.

    Determines whether the model contains both open and closed job classes.

    Args:
        sn: The stochastic network model

    Returns:
        bool: True if the model is mixed, False otherwise
    """
    return jpype.JPackage('jline').api.sn.SnIsMixedModelKt.snIsMixedModel(sn)


def sn_has_product_form_not_het_fcfs(sn):
    """
    Check if network has product-form but not heterogeneous FCFS.

    Determines whether the network has product-form properties but does not
    contain heterogeneous First-Come-First-Served stations.

    Args:
        sn: The stochastic network model

    Returns:
        bool: True if product-form without heterogeneous FCFS, False otherwise
    """
    return jpype.JPackage('jline').api.sn.SnHasProductFormNotHetFCFSKt.sn_has_product_form_not_het_fcfs(sn)


def sn_print_routing_matrix(sn, onlyClassName=None):
    """
    Print the routing matrix for the network.

    Displays the routing probability matrix showing how jobs move between stations.

    Args:
        sn: The stochastic network model
        onlyClassName: Restrict output to specific class name (str, optional)
    """
    jpype.JPackage('jline').api.sn.SnPrintRoutingMatrixKt.snPrintRoutingMatrix(sn, onlyClassName)


def sn_has_fcfs(sn):
    """
    Check if the network has First-Come-First-Served scheduling.

    Determines whether any stations use FCFS scheduling discipline.

    Args:
        sn: The stochastic network model

    Returns:
        bool: True if FCFS scheduling is present, False otherwise
    """
    return bool(jpype.JPackage('jline').api.sn.SnHasFCFSKt.snHasFCFS(sn))


def sn_has_lcfs(sn):
    """
    Check if the network has Last-Come-First-Served scheduling.

    Determines whether any stations use LCFS scheduling discipline.

    Args:
        sn: The stochastic network model

    Returns:
        bool: True if LCFS scheduling is present, False otherwise
    """
    return bool(jpype.JPackage('jline').api.sn.SnHasLCFSKt.snHasLCFS(sn))


def sn_has_lcfspr(sn):
    """
    Check if the network has LCFS with preemptive resume scheduling.

    Determines whether any stations use LCFS-PR scheduling discipline.

    Args:
        sn: The stochastic network model

    Returns:
        bool: True if LCFS-PR scheduling is present, False otherwise
    """
    return bool(jpype.JPackage('jline').api.sn.SnHasLCFSPRKt.snHasLCFSPR(sn))


def sn_has_ps(sn):
    """
    Check if the network has Processor Sharing scheduling.

    Determines whether any stations use PS scheduling discipline.

    Args:
        sn: The stochastic network model

    Returns:
        bool: True if PS scheduling is present, False otherwise
    """
    return bool(jpype.JPackage('jline').api.sn.SnHasPSKt.snHasPS(sn))


def sn_has_dps(sn):
    """
    Check if the network has Discriminatory Processor Sharing scheduling.

    Determines whether any stations use DPS scheduling discipline.

    Args:
        sn: The stochastic network model

    Returns:
        bool: True if DPS scheduling is present, False otherwise
    """
    return bool(jpype.JPackage('jline').api.sn.SnHasDPSKt.snHasDPS(sn))


def sn_has_gps(sn):
    """
    Check if the network has Generalized Processor Sharing scheduling.

    Determines whether any stations use GPS scheduling discipline.

    Args:
        sn: The stochastic network model

    Returns:
        bool: True if GPS scheduling is present, False otherwise
    """
    return bool(jpype.JPackage('jline').api.sn.SnHasGPSKt.snHasGPS(sn))


def sn_has_inf(sn):
    """
    Check if the network has infinite server stations.

    Determines whether any stations have unlimited server capacity.

    Args:
        sn: The stochastic network model

    Returns:
        bool: True if infinite server stations exist, False otherwise
    """
    return bool(jpype.JPackage('jline').api.sn.SnHasINFKt.snHasINF(sn))


def sn_has_hol(sn):
    """
    Check if the network has Head-of-Line priority scheduling.

    Determines whether any stations use HOL priority discipline.

    Args:
        sn: The stochastic network model

    Returns:
        bool: True if HOL priority is present, False otherwise
    """
    return bool(jpype.JPackage('jline').api.sn.SnHasHOLKt.snHasHOL(sn))


def sn_has_sjf(sn):
    """
    Check if the network has Shortest Job First scheduling.

    Determines whether any stations use SJF scheduling discipline.

    Args:
        sn: The stochastic network model

    Returns:
        bool: True if SJF scheduling is present, False otherwise
    """
    return bool(jpype.JPackage('jline').api.sn.SnHasSJFKt.snHasSJF(sn))


def sn_has_ljf(sn):
    """
    Check if the network has Longest Job First scheduling.

    Determines whether any stations use LJF scheduling discipline.

    Args:
        sn: The stochastic network model

    Returns:
        bool: True if LJF scheduling is present, False otherwise
    """
    return bool(jpype.JPackage('jline').api.sn.SnHasLJFKt.snHasLJF(sn))


def sn_has_sept(sn):
    """
    Check if the network has Shortest Expected Processing Time scheduling.

    Determines whether any stations use SEPT scheduling discipline.

    Args:
        sn: The stochastic network model

    Returns:
        bool: True if SEPT scheduling is present, False otherwise
    """
    return bool(jpype.JPackage('jline').api.sn.SnHasSEPTKt.snHasSEPT(sn))


def sn_has_lept(sn):
    """
    Check if the network has Longest Expected Processing Time scheduling.

    Determines whether any stations use LEPT scheduling discipline.

    Args:
        sn: The stochastic network model

    Returns:
        bool: True if LEPT scheduling is present, False otherwise
    """
    return bool(jpype.JPackage('jline').api.sn.SnHasLEPTKt.snHasLEPT(sn))


def sn_has_siro(sn):
    """
    Check if the network has Service In Random Order scheduling.

    Determines whether any stations use SIRO scheduling discipline.

    Args:
        sn: The stochastic network model

    Returns:
        bool: True if SIRO scheduling is present, False otherwise
    """
    return bool(jpype.JPackage('jline').api.sn.SnHasSIROKt.snHasSIRO(sn))


def sn_has_dps_prio(sn):
    """
    Check if the network has DPS with priority scheduling.

    Determines whether any stations use DPS with priority discipline.

    Args:
        sn: The stochastic network model

    Returns:
        bool: True if DPS-PRIO scheduling is present, False otherwise
    """
    return bool(jpype.JPackage('jline').api.sn.SnHasDPSPRIOKt.snHasDPSPRIO(sn))


def sn_has_gps_prio(sn):
    """
    Check if the network has GPS with priority scheduling.

    Determines whether any stations use GPS with priority discipline.

    Args:
        sn: The stochastic network model

    Returns:
        bool: True if GPS-PRIO scheduling is present, False otherwise
    """
    return bool(jpype.JPackage('jline').api.sn.SnHasGPSPRIOKt.snHasGPSPRIO(sn))


def sn_has_ps_prio(sn):
    """
    Check if the network has PS with priority scheduling.

    Determines whether any stations use PS with priority discipline.

    Args:
        sn: The stochastic network model

    Returns:
        bool: True if PS-PRIO scheduling is present, False otherwise
    """
    return bool(jpype.JPackage('jline').api.sn.SnHasPSPRIOKt.snHasPSPRIO(sn))


def sn_has_single_class(sn):
    """
    Check if the network has only a single job class.

    Determines whether the network contains exactly one job class.

    Args:
        sn: The stochastic network model

    Returns:
        bool: True if single class, False otherwise
    """
    return bool(jpype.JPackage('jline').api.sn.SnHasSingleClassKt.snHasSingleClass(sn))


def sn_has_single_chain(sn):
    """
    Check if the network has only a single routing chain.

    Determines whether the network contains exactly one routing chain.

    Args:
        sn: The stochastic network model

    Returns:
        bool: True if single chain, False otherwise
    """
    return bool(jpype.JPackage('jline').api.sn.SnHasSingleChainKt.snHasSingleChain(sn))


def sn_has_fractional_populations(sn):
    """
    Check if the network has non-integer job populations.

    Determines whether any job classes have fractional population values.

    Args:
        sn: The stochastic network model

    Returns:
        bool: True if fractional populations exist, False otherwise
    """
    return bool(jpype.JPackage('jline').api.sn.SnHasFractionalPopulationsKt.snHasFractionalPopulations(sn))


def sn_has_multiple_closed_classes(sn):
    """
    Check if the network has multiple closed job classes.

    Determines whether the network contains more than one closed job class.

    Args:
        sn: The stochastic network model

    Returns:
        bool: True if multiple closed classes exist, False otherwise
    """
    return bool(jpype.JPackage('jline').api.sn.SnHasMultipleClosedClassesKt.snHasMultipleClosedClasses(sn))


def sn_has_multiclass_fcfs(sn):
    """
    Check if the network has multi-class FCFS stations.

    Determines whether any FCFS stations serve multiple job classes.

    Args:
        sn: The stochastic network model

    Returns:
        bool: True if multi-class FCFS stations exist, False otherwise
    """
    return bool(jpype.JPackage('jline').api.sn.SnHasMultiClassFCFSKt.snHasMultiClassFCFS(sn))


def sn_has_multiclass_heter_fcfs(sn):
    """
    Check if the network has heterogeneous multi-class FCFS stations.

    Determines whether any FCFS stations serve multiple classes with
    different service time distributions.

    Args:
        sn: The stochastic network model

    Returns:
        bool: True if heterogeneous multi-class FCFS exists, False otherwise
    """
    return bool(jpype.JPackage('jline').api.sn.SnHasMultiClassHeterFCFSKt.snHasMultiClassHeterFCFS(sn))


def sn_has_multiclass_heter_exp_fcfs(sn):
    """
    Check if the network has heterogeneous exponential multi-class FCFS stations.

    Determines whether any FCFS stations serve multiple classes with
    different exponential service time distributions.

    Args:
        sn: The stochastic network model

    Returns:
        bool: True if heterogeneous exponential multi-class FCFS exists, False otherwise
    """
    return bool(jpype.JPackage('jline').api.sn.SnHasMultiClassHeterExpFCFSKt.snHasMultiClassHeterExpFCFS(sn))


def sn_has_homogeneous_scheduling(sn, strategy):
    """
    Check if the network has homogeneous scheduling for a given strategy.

    Determines whether all stations use the same scheduling discipline.

    Args:
        sn: The stochastic network model
        strategy: The scheduling strategy to check for

    Returns:
        bool: True if all stations use the specified strategy, False otherwise
    """
    return bool(jpype.JPackage('jline').api.sn.SnHasHomogeneousSchedulingKt.snHasHomogeneousScheduling(sn, strategy))


def sn_has_multi_class(sn):
    """
    Check if the network has multiple job classes.

    Determines whether the network contains more than one job class.

    Args:
        sn: The stochastic network model

    Returns:
        bool: True if multiple job classes exist, False otherwise
    """
    return bool(jpype.JPackage('jline').api.sn.SnHasMultiClassKt.snHasMultiClass(sn))


def sn_chain_analysis(sn, options=None):
    """
    Perform chain-level analysis of the stochastic network.

    Analyzes the network at the routing chain level, providing insights
    into chain behavior and performance characteristics.

    Args:
        sn: The stochastic network model
        options: Analysis options (optional)

    Returns:
        dict: Chain analysis results containing:
            - 'chain_info': Information about each chain
            - 'analysis_result': Performance analysis results
    """
    if options is None:
        options = {}

    java_result = jpype.JPackage('jline').api.sn.SnChainAnalysisKt.snChainAnalysis(sn, options)

    return {
        'chain_info': dict(java_result.getChainInfo()),
        'analysis_result': dict(java_result.getAnalysisResult())
    }


def sn_get_demands(sn, options=None):
    """
    Extract service demands from the stochastic network.

    Computes service demand matrix, service times, and visit ratios
    for all stations and job classes.

    Args:
        sn: The stochastic network model
        options: Extraction options (optional)

    Returns:
        tuple: (D, ST, V) where:
            - D: Service demands matrix
            - ST: Service times matrix
            - V: Visit ratios matrix
    """
    if options is None:
        options = {}

    from .. import jlineMatrixToArray

    java_result = jpype.JPackage('jline').api.sn.SnGetDemandsKt.snGetDemands(sn, options)

    D = jlineMatrixToArray(java_result.getD())
    ST = jlineMatrixToArray(java_result.getST())
    V = jlineMatrixToArray(java_result.getV())

    return D, ST, V


def sn_get_visits_chain(sn, options=None):
    """
    Get visit ratios at the chain level.

    Computes visit ratios aggregated at the routing chain level
    rather than individual class level.

    Args:
        sn: The stochastic network model
        options: Computation options (optional)

    Returns:
        numpy.ndarray: Chain-level visit ratios
    """
    if options is None:
        options = {}

    from .. import jlineMatrixToArray

    java_result = jpype.JPackage('jline').api.sn.SnGetVisitsChainKt.snGetVisitsChain(sn, options)

    return jlineMatrixToArray(java_result)


def sn_check_balance(sn, options=None):
    """
    Check traffic balance in the stochastic network.

    Verifies that flow conservation equations are satisfied at all nodes,
    ensuring the routing probabilities are consistent.

    Args:
        sn: The stochastic network model
        options: Check options (optional)

    Returns:
        dict: Balance check results:
            - 'is_balanced': Whether traffic is balanced
            - 'violations': List of balance violations
            - 'details': Detailed balance information
    """
    if options is None:
        options = {}

    java_result = jpype.JPackage('jline').api.sn.SnCheckBalanceKt.snCheckBalance(sn, options)

    return {
        'is_balanced': bool(java_result.isBalanced()),
        'violations': list(java_result.getViolations()),
        'details': dict(java_result.getDetails())
    }


def sn_check_consistency(sn, options=None):
    """
    Check model consistency in the stochastic network.

    Validates that the network model is internally consistent and
    all parameters are properly defined.

    Args:
        sn: The stochastic network model
        options: Check options (optional)

    Returns:
        dict: Consistency check results:
            - 'is_consistent': Whether the model is consistent
            - 'errors': List of consistency errors
            - 'warnings': List of warnings
            - 'details': Detailed consistency information
    """
    if options is None:
        options = {}

    java_result = jpype.JPackage('jline').api.sn.SnCheckConsistencyKt.snCheckConsistency(sn, options)

    return {
        'is_consistent': bool(java_result.isConsistent()),
        'errors': list(java_result.getErrors()),
        'warnings': list(java_result.getWarnings()),
        'details': dict(java_result.getDetails())
    }


def sn_check_feasibility(sn, options=None):
    """
    Check model feasibility in the stochastic network.

    Determines whether the network model represents a feasible queueing system
    that can be analyzed or simulated.

    Args:
        sn: The stochastic network model
        options: Check options (optional)

    Returns:
        dict: Feasibility check results:
            - 'is_feasible': Whether the model is feasible
            - 'issues': List of feasibility issues
            - 'recommendations': Suggested fixes
    """
    if options is None:
        options = {}

    java_result = jpype.JPackage('jline').api.sn.SnCheckFeasibilityKt.snCheckFeasibility(sn, options)

    return {
        'is_feasible': bool(java_result.isFeasible()),
        'issues': list(java_result.getIssues()),
        'recommendations': list(java_result.getRecommendations())
    }


def sn_has_blocking(sn):
    """
    Check if the network has blocking (finite capacity) stations.

    Determines whether any stations have finite buffer capacity that
    can block arriving customers.

    Args:
        sn: The stochastic network model

    Returns:
        bool: True if blocking stations exist, False otherwise
    """
    return bool(jpype.JPackage('jline').api.sn.SnHasBlockingKt.snHasBlocking(sn))


def sn_has_caches(sn):
    """
    Check if the network contains cache stations.

    Determines whether any stations implement caching behavior.

    Args:
        sn: The stochastic network model

    Returns:
        bool: True if cache stations exist, False otherwise
    """
    return bool(jpype.JPackage('jline').api.sn.SnHasCachesKt.snHasCaches(sn))


def sn_has_delays(sn):
    """
    Check if the network contains delay (infinite server) stations.

    Determines whether any stations are delay stations with no queueing.

    Args:
        sn: The stochastic network model

    Returns:
        bool: True if delay stations exist, False otherwise
    """
    return bool(jpype.JPackage('jline').api.sn.SnHasDelaysKt.snHasDelays(sn))


def sn_has_finite_capacity(sn):
    """
    Check if the network has stations with finite capacity.

    Determines whether any stations have capacity constraints.

    Args:
        sn: The stochastic network model

    Returns:
        bool: True if finite capacity stations exist, False otherwise
    """
    return bool(jpype.JPackage('jline').api.sn.SnHasFiniteCapacityKt.snHasFiniteCapacity(sn))


def sn_has_loadindep(sn):
    """
    Check if the network has load-independent service rates.

    Determines whether any stations have service rates that do not
    depend on the queue length.

    Args:
        sn: The stochastic network model

    Returns:
        bool: True if load-independent stations exist, False otherwise
    """
    return bool(jpype.JPackage('jline').api.sn.SnHasLoadindepKt.snHasLoadindep(sn))


def sn_has_state_dependent(sn):
    """
    Check if the network has state-dependent service rates.

    Determines whether any stations have service rates that depend
    on the system state.

    Args:
        sn: The stochastic network model

    Returns:
        bool: True if state-dependent stations exist, False otherwise
    """
    return bool(jpype.JPackage('jline').api.sn.SnHasStateDependentKt.snHasStateDependent(sn))


def sn_validate_model(sn, options=None):
    """
    Perform comprehensive validation of the stochastic network model.

    Runs all validation checks including consistency, feasibility,
    and balance checks to ensure the model is ready for analysis.

    Args:
        sn: The stochastic network model
        options: Validation options (optional)

    Returns:
        dict: Complete validation results:
            - 'is_valid': Overall validity status
            - 'consistency': Consistency check results
            - 'feasibility': Feasibility check results
            - 'balance': Traffic balance results
            - 'summary': Validation summary
    """
    if options is None:
        options = {}

    java_result = jpype.JPackage('jline').api.sn.SnValidateModelKt.snValidateModel(sn, options)

    return {
        'is_valid': bool(java_result.isValid()),
        'consistency': dict(java_result.getConsistency()),
        'feasibility': dict(java_result.getFeasibility()),
        'balance': dict(java_result.getBalance()),
        'summary': dict(java_result.getSummary())
    }

# Additional SN functions for complete API coverage

def sn_has_polling(sn):
    """
    Check if stochastic network has polling stations.
    
    Args:
        sn: NetworkStruct object
        
    Returns:
        bool: True if network has polling stations
    """
    result = jpype.JPackage('jline').api.sn.Sn_has_pollingKt.sn_has_polling(sn.obj)
    return bool(result)


def sn_has_rr(sn):
    """
    Check if stochastic network has round-robin scheduling.
    
    Args:
        sn: NetworkStruct object
        
    Returns:
        bool: True if network has RR stations
    """
    result = jpype.JPackage('jline').api.sn.Sn_has_rrKt.sn_has_rr(sn.obj)
    return bool(result)


def sn_has_slc(sn):
    """
    Check if stochastic network has self-looping classes.
    
    Args:
        sn: NetworkStruct object
        
    Returns:
        bool: True if network has self-looping classes
    """
    result = jpype.JPackage('jline').api.sn.Sn_has_slcKt.sn_has_slc(sn.obj)
    return bool(result)


def sn_has_snc(sn):
    """
    Check if stochastic network has state-dependent node capacities.
    
    Args:
        sn: NetworkStruct object
        
    Returns:
        bool: True if network has state-dependent node capacities
    """
    result = jpype.JPackage('jline').api.sn.Sn_has_sncKt.sn_has_snc(sn.obj)
    return bool(result)


def sn_has_srpt(sn):
    """
    Check if stochastic network has shortest remaining processing time scheduling.
    
    Args:
        sn: NetworkStruct object
        
    Returns:
        bool: True if network has SRPT stations
    """
    result = jpype.JPackage('jline').api.sn.Sn_has_srptKt.sn_has_srpt(sn.obj)
    return bool(result)


def sn_has_state_dependence(sn):
    """
    Check if stochastic network has state-dependent behavior.
    
    Args:
        sn: NetworkStruct object
        
    Returns:
        bool: True if network has state-dependent behavior
    """
    result = jpype.JPackage('jline').api.sn.Sn_has_state_dependenceKt.sn_has_state_dependence(sn.obj)
    return bool(result)


def sn_print(sn):
    """
    Print stochastic network structure information.
    
    Args:
        sn: NetworkStruct object
        
    Returns:
        None
    """
    jpype.JPackage('jline').api.sn.Sn_printKt.sn_print(sn.obj)


def sn_summary(sn):
    """
    Get summary information about stochastic network.
    
    Args:
        sn: NetworkStruct object
        
    Returns:
        dict: Summary information
    """
    result = jpype.JPackage('jline').api.sn.Sn_summaryKt.sn_summary(sn.obj)
    return {
        'stations': result.getStations(),
        'classes': result.getClasses(),
        'chains': result.getChains(),
        'features': str(result.getFeatures())
    }


def sn_validate(sn):
    """
    Validate stochastic network structure.
    
    Args:
        sn: NetworkStruct object
        
    Returns:
        tuple: (is_valid, error_messages)
    """
    result = jpype.JPackage('jline').api.sn.Sn_validateKt.sn_validate(sn.obj)
    return result.isValid(), list(result.getErrors())


def sn_get_arv_r_from_tput(sn, tput):
    """
    Computes arrival rates from throughput values.

    Note: This is an alias for sn_get_arvr_from_tput with corrected name.

    Args:
        sn: Stochastic network object
        tput: Throughput matrix

    Returns:
        numpy.ndarray: Arrival rates
    """
    return sn_get_arvr_from_tput(sn, tput)


def sn_get_node_arv_r_from_tput(sn, tput):
    """
    Computes node arrival rates from throughput values.

    Note: This is an alias for sn_get_node_arvr_from_tput with corrected name.

    Args:
        sn: Stochastic network object
        tput: Throughput matrix

    Returns:
        numpy.ndarray: Node arrival rates
    """
    return sn_get_node_arvr_from_tput(sn, tput)


def sn_has_lcfs_pr(sn):
    """
    Checks if network has LCFS-PR (Last Come First Served - Preemptive Resume) nodes.

    Note: This is an alias for sn_has_lcfspr.

    Args:
        sn: Stochastic network object

    Returns:
        bool: True if network has LCFS-PR nodes
    """
    return sn_has_lcfspr(sn)


def sn_has_multi_class_fcfs(sn):
    """
    Checks if network has multi-class FCFS nodes.

    Args:
        sn: Stochastic network object

    Returns:
        bool: True if network has multi-class FCFS nodes
    """
    result = jpype.JPackage('jline').api.sn.Sn_has_multi_class_fcfsKt.sn_has_multi_class_fcfs(sn.obj)
    return bool(result)


def sn_has_multi_class_heter_exp_fcfs(sn):
    """
    Checks if network has multi-class heterogeneous exponential FCFS nodes.

    Args:
        sn: Stochastic network object

    Returns:
        bool: True if network has multi-class heterogeneous exponential FCFS nodes
    """
    result = jpype.JPackage('jline').api.sn.Sn_has_multi_class_heter_exp_fcfsKt.sn_has_multi_class_heter_exp_fcfs(sn.obj)
    return bool(result)


def sn_has_multi_class_heter_fcfs(sn):
    """
    Checks if network has multi-class heterogeneous FCFS nodes.

    Args:
        sn: Stochastic network object

    Returns:
        bool: True if network has multi-class heterogeneous FCFS nodes
    """
    result = jpype.JPackage('jline').api.sn.Sn_has_multi_class_heter_fcfsKt.sn_has_multi_class_heter_fcfs(sn.obj)
    return bool(result)


def sn_is_population_model(sn):
    """
    Checks if the network is a population model.

    Args:
        sn: Stochastic network object

    Returns:
        bool: True if the network is a population model
    """
    result = jpype.JPackage('jline').api.sn.Sn_is_population_modelKt.sn_is_population_model(sn.obj)
    return bool(result)


def sn_rtnodes_to_rtorig(sn, rt_nodes):
    """
    Converts response times at nodes to original response times.

    Maps response times from individual nodes back to the original
    station/class structure of the network.

    Args:
        sn: Stochastic network object
        rt_nodes: Response times at nodes

    Returns:
        numpy.ndarray: Original response times
    """
    result = jpype.JPackage('jline').api.sn.Sn_rtnodes_to_rtorigKt.sn_rtnodes_to_rtorig(
        sn.obj, jlineMatrixFromArray(rt_nodes)
    )
    return jlineMatrixToArray(result)


# =============================================================================
# NetworkStruct Direct Modification Methods
# =============================================================================

def sn_set_service(sn, station_idx, class_idx, rate, scv=1.0, auto_refresh=False):
    """
    Set service rate at a specific station and class.

    Directly modifies the service rate in NetworkStruct without rebuilding
    the full Network object. Useful for fast parameter updates in optimization.

    Args:
        sn: The stochastic network structure
        station_idx: Station index (0-based)
        class_idx: Class index (0-based)
        rate: New service rate (must be positive)
        scv: Squared coefficient of variation (default 1.0 for exponential)
        auto_refresh: If True, refresh process fields after modification

    Returns:
        The modified NetworkStruct
    """
    ModifyMode = jpype.JPackage('jline').api.sn.ModifyMode
    ValidationLevel = jpype.JPackage('jline').api.sn.ValidationLevel

    return jpype.JPackage('jline').api.sn.SnSetServiceKt.snSetService(
        sn, station_idx, class_idx, rate, scv,
        ModifyMode.IN_PLACE, ValidationLevel.MINIMAL, auto_refresh
    )


def sn_set_service_batch(sn, rates, scvs=None, auto_refresh=False):
    """
    Set service rates for multiple station-class pairs.

    Batch update of service rates. NaN values are skipped.
    More efficient than calling sn_set_service multiple times.

    Args:
        sn: The stochastic network structure
        rates: Matrix of new rates (nstations x nclasses), NaN = skip
        scvs: Matrix of new SCVs (optional), NaN = skip
        auto_refresh: If True, refresh process fields after modification

    Returns:
        The modified NetworkStruct
    """
    ModifyMode = jpype.JPackage('jline').api.sn.ModifyMode
    ValidationLevel = jpype.JPackage('jline').api.sn.ValidationLevel

    rates_matrix = jlineMatrixFromArray(np.array(rates))
    scvs_matrix = jlineMatrixFromArray(np.array(scvs)) if scvs is not None else None

    return jpype.JPackage('jline').api.sn.SnSetServiceKt.snSetServiceBatch(
        sn, rates_matrix, scvs_matrix,
        ModifyMode.IN_PLACE, ValidationLevel.MINIMAL, auto_refresh
    )


def sn_set_arrival(sn, class_idx, rate, scv=1.0, auto_refresh=False):
    """
    Set arrival rate for a class at the Source station.

    Args:
        sn: The stochastic network structure
        class_idx: Class index (0-based)
        rate: New arrival rate (lambda)
        scv: Squared coefficient of variation (default 1.0)
        auto_refresh: If True, refresh process fields after modification

    Returns:
        The modified NetworkStruct
    """
    ModifyMode = jpype.JPackage('jline').api.sn.ModifyMode
    ValidationLevel = jpype.JPackage('jline').api.sn.ValidationLevel

    return jpype.JPackage('jline').api.sn.SnSetArrivalKt.snSetArrival(
        sn, class_idx, rate, scv,
        ModifyMode.IN_PLACE, ValidationLevel.MINIMAL, auto_refresh
    )


def sn_set_servers(sn, station_idx, n_servers):
    """
    Set the number of servers at a station.

    Args:
        sn: The stochastic network structure
        station_idx: Station index (0-based)
        n_servers: Number of servers (positive, or float('inf'))

    Returns:
        The modified NetworkStruct
    """
    ModifyMode = jpype.JPackage('jline').api.sn.ModifyMode
    ValidationLevel = jpype.JPackage('jline').api.sn.ValidationLevel

    return jpype.JPackage('jline').api.sn.SnSetServersKt.snSetServers(
        sn, station_idx, float(n_servers),
        ModifyMode.IN_PLACE, ValidationLevel.MINIMAL
    )


def sn_set_population(sn, class_idx, n_jobs, auto_refresh=False):
    """
    Set the number of jobs for a closed class.

    Args:
        sn: The stochastic network structure
        class_idx: Class index (0-based)
        n_jobs: Number of jobs (non-negative)
        auto_refresh: If True, refresh visit ratios after modification

    Returns:
        The modified NetworkStruct
    """
    ModifyMode = jpype.JPackage('jline').api.sn.ModifyMode
    ValidationLevel = jpype.JPackage('jline').api.sn.ValidationLevel

    return jpype.JPackage('jline').api.sn.SnSetPopulationKt.snSetPopulation(
        sn, class_idx, float(n_jobs),
        ModifyMode.IN_PLACE, ValidationLevel.MINIMAL, auto_refresh
    )


def sn_set_priority(sn, class_idx, priority):
    """
    Set the priority for a class.

    Args:
        sn: The stochastic network structure
        class_idx: Class index (0-based)
        priority: Priority value

    Returns:
        The modified NetworkStruct
    """
    ModifyMode = jpype.JPackage('jline').api.sn.ModifyMode
    ValidationLevel = jpype.JPackage('jline').api.sn.ValidationLevel

    return jpype.JPackage('jline').api.sn.SnSetPriorityKt.snSetPriority(
        sn, class_idx, float(priority),
        ModifyMode.IN_PLACE, ValidationLevel.MINIMAL
    )


def sn_set_routing(sn, rt, auto_refresh=False):
    """
    Set the routing matrix for stateful nodes.

    Args:
        sn: The stochastic network structure
        rt: New routing matrix (nstateful*nclasses x nstateful*nclasses)
        auto_refresh: If True, refresh visit ratios after modification

    Returns:
        The modified NetworkStruct
    """
    ModifyMode = jpype.JPackage('jline').api.sn.ModifyMode
    ValidationLevel = jpype.JPackage('jline').api.sn.ValidationLevel

    rt_matrix = jlineMatrixFromArray(np.array(rt))

    return jpype.JPackage('jline').api.sn.SnSetRoutingKt.snSetRoutingMatrix(
        sn, rt_matrix,
        ModifyMode.IN_PLACE, ValidationLevel.MINIMAL, auto_refresh
    )


def sn_set_routing_prob(sn, from_stateful, from_class, to_stateful, to_class, prob, auto_refresh=False):
    """
    Set a routing probability between two stateful node-class pairs.

    Args:
        sn: The stochastic network structure
        from_stateful: Source stateful node index (0-based)
        from_class: Source class index (0-based)
        to_stateful: Destination stateful node index (0-based)
        to_class: Destination class index (0-based)
        prob: Routing probability [0, 1]
        auto_refresh: If True, refresh visit ratios after modification

    Returns:
        The modified NetworkStruct
    """
    ModifyMode = jpype.JPackage('jline').api.sn.ModifyMode
    ValidationLevel = jpype.JPackage('jline').api.sn.ValidationLevel

    return jpype.JPackage('jline').api.sn.SnSetRoutingKt.snSetRoutingProb(
        sn, from_stateful, from_class, to_stateful, to_class, float(prob),
        ModifyMode.IN_PLACE, ValidationLevel.MINIMAL, auto_refresh
    )


def sn_set_fork_fanout(sn, fork_node_idx, fan_out):
    """
    Set fork fanout (tasksPerLink) for a Fork node.

    Args:
        sn: The stochastic network structure
        fork_node_idx: Node index of the Fork node (0-based)
        fan_out: Number of tasks per output link (>= 1)

    Returns:
        The modified NetworkStruct
    """
    ModifyMode = jpype.JPackage('jline').api.sn.ModifyMode
    ValidationLevel = jpype.JPackage('jline').api.sn.ValidationLevel

    return jpype.JPackage('jline').api.sn.SnSetForkFanoutKt.snSetForkFanout(
        sn, fork_node_idx, int(fan_out),
        ModifyMode.IN_PLACE, ValidationLevel.MINIMAL
    )


def sn_refresh_process_fields(sn, station_idx, class_idx):
    """
    Refresh process fields for a station-class pair based on rate and SCV.

    Updates mu, phi, proc, pie, phases based on current rate and SCV values.

    Args:
        sn: The stochastic network structure
        station_idx: Station index (0-based)
        class_idx: Class index (0-based)

    Returns:
        The modified NetworkStruct
    """
    return jpype.JPackage('jline').api.sn.SnRefreshProcessFieldsKt.snRefreshProcessFields(
        sn, station_idx, class_idx
    )


def sn_refresh_all_process_fields(sn):
    """
    Refresh process fields for all station-class pairs.

    Args:
        sn: The stochastic network structure

    Returns:
        The modified NetworkStruct
    """
    return jpype.JPackage('jline').api.sn.SnRefreshProcessFieldsKt.snRefreshAllProcessFields(sn)
