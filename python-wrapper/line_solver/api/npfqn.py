
"""
Non-product-form queueing network (NPFQN) approximations.

This module provides approximation algorithms for queueing networks
that do not satisfy product-form conditions, including traffic
merging and splitting operations for non-exponential service times.

These approximations enable analysis of more general queueing
networks beyond the product-form assumption.
"""

import jpype
from line_solver import jlineMatrixToArray, jlineMatrixFromArray


def npfqn_nonexp_approx(method, sn, ST, V, SCV, Tin, Uin, gamma, nservers):
    """
    Non-exponential approximation for non-product-form queueing networks.
    
    Applies approximation methods to analyze queueing networks with
    non-exponential service time distributions.
    
    Args:
        method: Approximation method identifier.
        sn: Service network structure.
        ST: Service time matrix.
        V: Visit count matrix.
        SCV: Squared coefficient of variation matrix.
        Tin: Input think times.
        Uin: Input utilizations.
        gamma: Load balancing parameters.
        nservers: Number of servers per station.
        
    Returns:
        list: [ST_new, gamma_new, nservers_new, rho, scva, scvs, eta].
    """
    method_str = str(method)

    ret = jpype.JPackage('jline').api.npfqn.Npfqn_nonexp_approxKt.npfqn_nonexp_approx(
        method_str,
        sn,
        jlineMatrixFromArray(ST),
        jlineMatrixFromArray(V),
        jlineMatrixFromArray(SCV),
        jlineMatrixFromArray(Tin),
        jlineMatrixFromArray(Uin),
        jlineMatrixFromArray(gamma),
        jlineMatrixFromArray(nservers)
    )

    ST_new = jlineMatrixToArray(ret.get(0))
    gamma_new = jlineMatrixToArray(ret.get(1))
    nservers_new = jlineMatrixToArray(ret.get(2))
    rho = jlineMatrixToArray(ret.get(3))
    scva = jlineMatrixToArray(ret.get(4))
    scvs = jlineMatrixToArray(ret.get(5))
    eta = jlineMatrixToArray(ret.get(6))

    return [ST_new, gamma_new, nservers_new, rho, scva, scvs, eta]


def npfqn_traffic_merge(lambda_rates, scv_rates):
    """
    Merge traffic streams with different arrival rates and variabilities.
    
    Combines multiple traffic streams into a single merged stream,
    computing the resulting arrival rate and variability.
    
    Args:
        lambda_rates: Array of arrival rates for each stream.
        scv_rates: Array of squared coefficients of variation.
        
    Returns:
        tuple: (lambda_merged, scv_merged) - Merged traffic parameters.
    """
    result = jpype.JPackage('jline').api.npfqn.Npfqn_traffic_mergeKt.npfqn_traffic_merge(
        jlineMatrixFromArray(lambda_rates),
        jlineMatrixFromArray(scv_rates)
    )

    return result.lambda_merged, result.scv_merged


def npfqn_traffic_merge_cs(lambda_rates, scv_rates, P):
    """
    Merge traffic streams with class switching probabilities.
    
    Combines multiple traffic streams considering class switching
    behavior and routing probabilities.
    
    Args:
        lambda_rates: Array of arrival rates for each stream.
        scv_rates: Array of squared coefficients of variation.
        P: Class switching probability matrix.
        
    Returns:
        tuple: (lambda_merged, scv_merged) - Merged traffic with class switching.
    """
    result = jpype.JPackage('jline').api.npfqn.Npfqn_traffic_merge_csKt.npfqn_traffic_merge_cs(
        jlineMatrixFromArray(lambda_rates),
        jlineMatrixFromArray(scv_rates),
        jlineMatrixFromArray(P)
    )

    return jlineMatrixToArray(result.lambda_merged), jlineMatrixToArray(result.scv_merged)


def npfqn_traffic_split_cs(lambda_rate, scv_rate, P):
    """
    Split traffic stream with class switching probabilities.
    
    Divides a single traffic stream into multiple streams based
    on class switching probabilities.
    
    Args:
        lambda_rate: Input arrival rate.
        scv_rate: Input squared coefficient of variation.
        P: Class switching probability matrix.
        
    Returns:
        tuple: (lambda_split, scv_split) - Split traffic streams.
    """
    result = jpype.JPackage('jline').api.npfqn.Npfqn_traffic_split_csKt.npfqn_traffic_split_cs(
        jpype.JDouble(lambda_rate),
        jpype.JDouble(scv_rate),
        jlineMatrixFromArray(P)
    )

    return jlineMatrixToArray(result.lambda_split), jlineMatrixToArray(result.scv_split)