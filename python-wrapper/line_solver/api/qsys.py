
"""
Single queueing system analysis functions.

This module provides analytical formulas and approximations for single
queueing systems, including classical models like M/M/1, M/M/k, M/G/1,
and various G/G/1 approximations.

Supported queueing systems:
- qsys_mm1: M/M/1 queue (Poisson arrivals, exponential service)
- qsys_mmk: M/M/k queue (Poisson arrivals, k exponential servers)
- qsys_gm1: G/M/1 queue (general arrivals, exponential service)  
- qsys_mg1: M/G/1 queue (Poisson arrivals, general service)
- Various G/G/1 approximations (Whitt, Allen-Cunneen, Kingman, etc.)

These functions compute exact results where available, or high-quality
approximations for more general cases.
"""

import jpype
import numpy as np
from line_solver import jlineMatrixToArray, jlineMatrixFromArray


def qsys_mm1(lambda_val, mu):
    """
    Analyze M/M/1 queue (Poisson arrivals, exponential service).

    Args:
        lambda_val (float): Arrival rate.
        mu (float): Service rate.

    Returns:
        dict: Performance measures including:
            - L: Mean number in system
            - Lq: Mean number in queue
            - W: Mean response time
            - Wq: Mean waiting time
            - rho: Utilization (λ/μ)
    """
    result = jpype.JPackage('jline').api.qsys.Qsys_mm1Kt.qsys_mm1(
        jpype.JDouble(lambda_val), jpype.JDouble(mu)
    )

    # Java returns W and rho; derive other metrics using queueing formulas
    W = float(result.W)
    rho = float(result.rho)
    Wq = W - 1.0 / mu  # Mean waiting time = response time - service time
    L = lambda_val * W  # Little's law
    Lq = lambda_val * Wq  # Little's law for queue

    return {
        'L': L,
        'Lq': Lq,
        'W': W,
        'Wq': Wq,
        'rho': rho
    }


def qsys_mmk(lambda_val, mu, k):
    """
    Analyze M/M/k queue (Poisson arrivals, k exponential servers).

    Args:
        lambda_val (float): Arrival rate.
        mu (float): Service rate per server.
        k (int): Number of parallel servers.

    Returns:
        dict: Performance measures including L, Lq, W, Wq, rho.
    """
    result = jpype.JPackage('jline').api.qsys.Qsys_mmkKt.qsys_mmk(
        jpype.JDouble(lambda_val), jpype.JDouble(mu), jpype.JInt(k)
    )

    # Java returns W and rho; derive other metrics
    W = float(result.W)
    rho = float(result.rho)
    Wq = W - 1.0 / mu  # Mean waiting time = response time - service time
    L = lambda_val * W  # Little's law
    Lq = lambda_val * Wq  # Little's law for queue

    return {
        'L': L,
        'Lq': Lq,
        'W': W,
        'Wq': Wq,
        'rho': rho
    }


def qsys_gm1(lambda_val, mu, sigma_s_squared):
    """
    Analyze G/M/1 queue (general arrivals, exponential service).

    Args:
        lambda_val (float): Arrival rate.
        mu (float): Service rate.
        sigma_s_squared (float): Variance of inter-arrival times.

    Returns:
        dict: Performance measures including L, Lq, W, Wq, rho.
    """
    result = jpype.JPackage('jline').api.qsys.Qsys_gm1Kt.qsys_gm1(
        jpype.JDouble(lambda_val), jpype.JDouble(mu), jpype.JDouble(sigma_s_squared)
    )

    # Java returns W and rho; derive other metrics
    W = float(result.W)
    rho = float(result.rho)
    Wq = W - 1.0 / mu
    L = lambda_val * W
    Lq = lambda_val * Wq

    return {
        'L': L,
        'Lq': Lq,
        'W': W,
        'Wq': Wq,
        'rho': rho
    }


def qsys_mg1(lambda_val, mu, sigma_s_squared):
    """
    Analyze M/G/1 queue (Poisson arrivals, general service).

    Args:
        lambda_val (float): Arrival rate.
        mu (float): Mean service rate.
        sigma_s_squared (float): Variance of service times.

    Returns:
        dict: Performance measures including L, Lq, W, Wq, rho.
    """
    result = jpype.JPackage('jline').api.qsys.Qsys_mg1Kt.qsys_mg1(
        jpype.JDouble(lambda_val), jpype.JDouble(mu), jpype.JDouble(sigma_s_squared)
    )

    # Java returns W and rho; derive other metrics
    W = float(result.W)
    rho = float(result.rho)
    Wq = W - 1.0 / mu
    L = lambda_val * W
    Lq = lambda_val * Wq

    return {
        'L': L,
        'Lq': Lq,
        'W': W,
        'Wq': Wq,
        'rho': rho
    }


def qsys_gig1_approx_lin(lambda_val, mu, ca_squared, cs_squared):
    """
    G/G/1 queue approximation using linear interpolation.

    Args:
        lambda_val (float): Arrival rate.
        mu (float): Service rate.
        ca_squared (float): Squared coefficient of variation of arrivals.
        cs_squared (float): Squared coefficient of variation of service.

    Returns:
        dict: Approximate performance measures.
    """
    result = jpype.JPackage('jline').api.qsys.Qsys_gig1_approx_linKt.qsys_gig1_approx_lin(
        jpype.JDouble(lambda_val), jpype.JDouble(mu),
        jpype.JDouble(ca_squared), jpype.JDouble(cs_squared)
    )

    W = float(result.W)
    rho = float(result.rho)
    Wq = W - 1.0 / mu
    L = lambda_val * W
    Lq = lambda_val * Wq

    return {
        'L': L,
        'Lq': Lq,
        'W': W,
        'Wq': Wq,
        'rho': rho
    }


def qsys_gig1_approx_kk(lambda_val, mu, ca_squared, cs_squared):
    """
    G/G/1 queue approximation using Kraemer-Langenbach-Belz method.

    Args:
        lambda_val (float): Arrival rate.
        mu (float): Service rate.
        ca_squared (float): Squared coefficient of variation of arrivals.
        cs_squared (float): Squared coefficient of variation of service.

    Returns:
        dict: Approximate performance measures.
    """
    result = jpype.JPackage('jline').api.qsys.Qsys_gig1_approx_kkKt.qsys_gig1_approx_kk(
        jpype.JDouble(lambda_val), jpype.JDouble(mu),
        jpype.JDouble(ca_squared), jpype.JDouble(cs_squared)
    )

    W = float(result.W)
    rho = float(result.rho)
    Wq = W - 1.0 / mu
    L = lambda_val * W
    Lq = lambda_val * Wq

    return {
        'L': L,
        'Lq': Lq,
        'W': W,
        'Wq': Wq,
        'rho': rho
    }


def qsys_gig1_approx_whitt(lambda_val, mu, ca_squared, cs_squared):
    """
    G/G/1 queue approximation using Whitt's method.

    Args:
        lambda_val (float): Arrival rate.
        mu (float): Service rate.
        ca_squared (float): Squared coefficient of variation of arrivals.
        cs_squared (float): Squared coefficient of variation of service.

    Returns:
        dict: Approximate performance measures.
    """
    result = jpype.JPackage('jline').api.qsys.Qsys_gig1_approx_whittKt.qsys_gig1_approx_whitt(
        jpype.JDouble(lambda_val), jpype.JDouble(mu),
        jpype.JDouble(ca_squared), jpype.JDouble(cs_squared)
    )

    W = float(result.W)
    rho = float(result.rho)
    Wq = W - 1.0 / mu
    L = lambda_val * W
    Lq = lambda_val * Wq

    return {
        'L': L,
        'Lq': Lq,
        'W': W,
        'Wq': Wq,
        'rho': rho
    }

def qsys_gig1_approx_allencunneen(lambda_val, mu, ca_squared, cs_squared):
    """
    G/G/1 queue approximation using Allen-Cunneen method.

    Args:
        lambda_val: Arrival rate
        mu: Service rate
        ca_squared: Squared coefficient of variation of arrivals
        cs_squared: Squared coefficient of variation of service

    Returns:
        dict: Performance measures (W, rho)
    """
    result = jpype.JPackage('jline').api.qsys.Qsys_gig1_approx_allencunneenKt.qsys_gig1_approx_allencunneen(
        jpype.JDouble(lambda_val), jpype.JDouble(mu),
        jpype.JDouble(ca_squared), jpype.JDouble(cs_squared)
    )

    return {
        'W': result.W,
        'rho': result.rho
    }


def qsys_gig1_approx_heyman(lambda_val, mu, ca_squared, cs_squared):
    """
    G/G/1 queue approximation using Heyman method.

    Args:
        lambda_val: Arrival rate
        mu: Service rate
        ca_squared: Squared coefficient of variation of arrivals
        cs_squared: Squared coefficient of variation of service

    Returns:
        dict: Performance measures (W, rho)
    """
    result = jpype.JPackage('jline').api.qsys.Qsys_gig1_approx_heymanKt.qsys_gig1_approx_heyman(
        jpype.JDouble(lambda_val), jpype.JDouble(mu),
        jpype.JDouble(ca_squared), jpype.JDouble(cs_squared)
    )

    return {
        'W': result.W,
        'rho': result.rho
    }


def qsys_gig1_approx_kobayashi(lambda_val, mu, ca_squared, cs_squared):
    """
    G/G/1 queue approximation using Kobayashi method.

    Args:
        lambda_val: Arrival rate
        mu: Service rate
        ca_squared: Squared coefficient of variation of arrivals
        cs_squared: Squared coefficient of variation of service

    Returns:
        dict: Performance measures (W, rho)
    """
    result = jpype.JPackage('jline').api.qsys.Qsys_gig1_approx_kobayashiKt.qsys_gig1_approx_kobayashi(
        jpype.JDouble(lambda_val), jpype.JDouble(mu),
        jpype.JDouble(ca_squared), jpype.JDouble(cs_squared)
    )

    return {
        'W': result.W,
        'rho': result.rho
    }


def qsys_gig1_approx_marchal(lambda_val, mu, ca_squared, cs_squared):
    """
    G/G/1 queue approximation using Marchal method.

    Args:
        lambda_val: Arrival rate
        mu: Service rate
        ca_squared: Squared coefficient of variation of arrivals
        cs_squared: Squared coefficient of variation of service

    Returns:
        dict: Performance measures (W, rho)
    """
    result = jpype.JPackage('jline').api.qsys.Qsys_gig1_approx_marchalKt.qsys_gig1_approx_marchal(
        jpype.JDouble(lambda_val), jpype.JDouble(mu),
        jpype.JDouble(ca_squared), jpype.JDouble(cs_squared)
    )

    return {
        'W': result.W,
        'rho': result.rho
    }


def qsys_gig1_ubnd_kingman(lambda_val, mu, ca_squared, cs_squared):
    """
    G/G/1 queue upper bound using Kingman method.

    Args:
        lambda_val: Arrival rate
        mu: Service rate
        ca_squared: Squared coefficient of variation of arrivals
        cs_squared: Squared coefficient of variation of service

    Returns:
        dict: Upper bound performance measures (W, rho)
    """
    result = jpype.JPackage('jline').api.qsys.Qsys_gig1_ubnd_kingmanKt.qsys_gig1_ubnd_kingman(
        jpype.JDouble(lambda_val), jpype.JDouble(mu),
        jpype.JDouble(ca_squared), jpype.JDouble(cs_squared)
    )

    return {
        'W': result.W,
        'rho': result.rho
    }


def qsys_gigk_approx(lambda_val, mu, ca_squared, cs_squared, k):
    """
    G/G/k queue approximation for k servers.

    Args:
        lambda_val: Arrival rate
        mu: Service rate per server
        ca_squared: Squared coefficient of variation of arrivals
        cs_squared: Squared coefficient of variation of service
        k: Number of servers

    Returns:
        dict: Performance measures (W, rho)
    """
    result = jpype.JPackage('jline').api.qsys.Qsys_gigk_approxKt.qsys_gigk_approx(
        jpype.JDouble(lambda_val), jpype.JDouble(mu),
        jpype.JDouble(ca_squared), jpype.JDouble(cs_squared), jpype.JInt(k)
    )

    return {
        'W': result.W,
        'rho': result.rho
    }


def qsys_gigk_approx_kingman(lambda_val, mu, ca_squared, cs_squared, k):
    """
    G/G/k queue approximation using Kingman method.

    Args:
        lambda_val: Arrival rate
        mu: Service rate per server
        ca_squared: Squared coefficient of variation of arrivals
        cs_squared: Squared coefficient of variation of service
        k: Number of servers

    Returns:
        dict: Performance measures (W, rho)
    """
    result = jpype.JPackage('jline').api.qsys.Qsys_gigk_approx_kingmanKt.qsys_gigk_approx_kingman(
        jpype.JDouble(lambda_val), jpype.JDouble(mu),
        jpype.JDouble(ca_squared), jpype.JDouble(cs_squared), jpype.JInt(k)
    )

    return {
        'W': result.W,
        'rho': result.rho
    }


def qsys_gig1_approx_klb(lambda_val, mu, ca_squared, cs_squared):
    """
    G/G/1 queue approximation using Kraemer-Langenbach-Belz method.

    Args:
        lambda_val: Arrival rate
        mu: Service rate
        ca_squared: Squared coefficient of variation of arrivals
        cs_squared: Squared coefficient of variation of service

    Returns:
        dict: Performance measures (L, Lq, W, Wq, rho)
    """
    result = jpype.JPackage('jline').api.qsys.Qsys_gig1_approx_klbKt.qsys_gig1_approx_klb(
        jpype.JDouble(lambda_val), jpype.JDouble(mu),
        jpype.JDouble(ca_squared), jpype.JDouble(cs_squared)
    )

    return {
        'L': result.L,
        'Lq': result.Lq,
        'W': result.W,
        'Wq': result.Wq,
        'rho': result.rho
    }



# Additional QSYS functions for complete API coverage

def qsys_gg1(arrival_mean, arrival_scv, service_mean, service_scv):
    """
    Analyze G/G/1 queueing system.
    
    Args:
        arrival_mean: Mean inter-arrival time
        arrival_scv: Squared coefficient of variation of inter-arrival times
        service_mean: Mean service time
        service_scv: Squared coefficient of variation of service times
        
    Returns:
        dict: Performance metrics including utilization, queue length, and waiting time
    """
    result = jpype.JPackage('jline').api.qsys.Qsys_gg1Kt.qsys_gg1(
        float(arrival_mean), float(arrival_scv), 
        float(service_mean), float(service_scv)
    )
    return {
        'utilization': float(result.getUtilization()),
        'queue_length': float(result.getQueueLength()),
        'waiting_time': float(result.getWaitingTime()),
        'response_time': float(result.getResponseTime())
    }


def qsys_gig1_approx_gelenbe(arrival_mean, arrival_scv, service_mean, service_scv):
    """
    G/G/1 approximation using Gelenbe's method.
    
    Args:
        arrival_mean: Mean inter-arrival time
        arrival_scv: Squared coefficient of variation of inter-arrival times
        service_mean: Mean service time
        service_scv: Squared coefficient of variation of service times
        
    Returns:
        dict: Approximate performance metrics
    """
    result = jpype.JPackage('jline').api.qsys.Qsys_gig1_approx_gelenbeKt.qsys_gig1_approx_gelenbe(
        float(arrival_mean), float(arrival_scv),
        float(service_mean), float(service_scv)
    )
    return {
        'queue_length': float(result.getQueueLength()),
        'waiting_time': float(result.getWaitingTime())
    }


def qsys_gig1_approx_kimura(arrival_mean, arrival_scv, service_mean, service_scv):
    """
    G/G/1 approximation using Kimura's method.
    
    Args:
        arrival_mean: Mean inter-arrival time
        arrival_scv: Squared coefficient of variation of inter-arrival times
        service_mean: Mean service time
        service_scv: Squared coefficient of variation of service times
        
    Returns:
        dict: Approximate performance metrics
    """
    result = jpype.JPackage('jline').api.qsys.Qsys_gig1_approx_kimuraKt.qsys_gig1_approx_kimura(
        float(arrival_mean), float(arrival_scv),
        float(service_mean), float(service_scv)
    )
    return {
        'queue_length': float(result.getQueueLength()),
        'waiting_time': float(result.getWaitingTime())
    }


def qsys_gig1_approx_klb(arrival_mean, arrival_scv, service_mean, service_scv):
    """
    G/G/1 approximation using Kraemer-Langenbach-Belz formula.
    
    Args:
        arrival_mean: Mean inter-arrival time
        arrival_scv: Squared coefficient of variation of inter-arrival times
        service_mean: Mean service time
        service_scv: Squared coefficient of variation of service times
        
    Returns:
        dict: Approximate performance metrics
    """
    result = jpype.JPackage('jline').api.qsys.Qsys_gig1_approx_klbKt.qsys_gig1_approx_klb(
        float(arrival_mean), float(arrival_scv),
        float(service_mean), float(service_scv)
    )
    return {
        'queue_length': float(result.getQueueLength()),
        'waiting_time': float(result.getWaitingTime())
    }


def qsys_gig1_approx_myskja(arrival_mean, arrival_scv, service_mean, service_scv):
    """
    G/G/1 approximation using Myskja's method.
    
    Args:
        arrival_mean: Mean inter-arrival time
        arrival_scv: Squared coefficient of variation of inter-arrival times
        service_mean: Mean service time
        service_scv: Squared coefficient of variation of service times
        
    Returns:
        dict: Approximate performance metrics
    """
    result = jpype.JPackage('jline').api.qsys.Qsys_gig1_approx_myskjaKt.qsys_gig1_approx_myskja(
        float(arrival_mean), float(arrival_scv),
        float(service_mean), float(service_scv)
    )
    return {
        'queue_length': float(result.getQueueLength()),
        'waiting_time': float(result.getWaitingTime())
    }


def qsys_gig1_approx_myskja2(arrival_mean, arrival_scv, service_mean, service_scv):
    """
    G/G/1 approximation using Myskja's refined method.
    
    Args:
        arrival_mean: Mean inter-arrival time
        arrival_scv: Squared coefficient of variation of inter-arrival times
        service_mean: Mean service time
        service_scv: Squared coefficient of variation of service times
        
    Returns:
        dict: Approximate performance metrics
    """
    result = jpype.JPackage('jline').api.qsys.Qsys_gig1_approx_myskja2Kt.qsys_gig1_approx_myskja2(
        float(arrival_mean), float(arrival_scv),
        float(service_mean), float(service_scv)
    )
    return {
        'queue_length': float(result.getQueueLength()),
        'waiting_time': float(result.getWaitingTime())
    }


def qsys_gig1_lbnd(arrival_mean, arrival_scv, service_mean, service_scv):
    """
    Lower bound for G/G/1 queue performance.
    
    Args:
        arrival_mean: Mean inter-arrival time
        arrival_scv: Squared coefficient of variation of inter-arrival times
        service_mean: Mean service time
        service_scv: Squared coefficient of variation of service times
        
    Returns:
        dict: Lower bounds on performance metrics
    """
    result = jpype.JPackage('jline').api.qsys.Qsys_gig1_lbndKt.qsys_gig1_lbnd(
        float(arrival_mean), float(arrival_scv),
        float(service_mean), float(service_scv)
    )
    return {
        'queue_length_lower': float(result.getQueueLengthLower()),
        'waiting_time_lower': float(result.getWaitingTimeLower())
    }


def qsys_gigk_approx_cosmetatos(arrival_mean, arrival_scv, service_mean, service_scv, num_servers):
    """
    G/G/k approximation using Cosmetatos method.
    
    Args:
        arrival_mean: Mean inter-arrival time
        arrival_scv: Squared coefficient of variation of inter-arrival times
        service_mean: Mean service time
        service_scv: Squared coefficient of variation of service times
        num_servers: Number of servers
        
    Returns:
        dict: Approximate performance metrics
    """
    result = jpype.JPackage('jline').api.qsys.Qsys_gigk_approx_cosmetatosKt.qsys_gigk_approx_cosmetatos(
        float(arrival_mean), float(arrival_scv),
        float(service_mean), float(service_scv),
        int(num_servers)
    )
    return {
        'queue_length': float(result.getQueueLength()),
        'waiting_time': float(result.getWaitingTime()),
        'utilization': float(result.getUtilization())
    }


def qsys_gigk_approx_whitt(arrival_mean, arrival_scv, service_mean, service_scv, num_servers):
    """
    G/G/k approximation using Whitt's method.
    
    Args:
        arrival_mean: Mean inter-arrival time
        arrival_scv: Squared coefficient of variation of inter-arrival times
        service_mean: Mean service time
        service_scv: Squared coefficient of variation of service times
        num_servers: Number of servers
        
    Returns:
        dict: Approximate performance metrics
    """
    result = jpype.JPackage('jline').api.qsys.Qsys_gigk_approx_whittKt.qsys_gigk_approx_whitt(
        float(arrival_mean), float(arrival_scv),
        float(service_mean), float(service_scv),
        int(num_servers)
    )
    return {
        'queue_length': float(result.getQueueLength()),
        'waiting_time': float(result.getWaitingTime()),
        'utilization': float(result.getUtilization())
    }


def qsys_mg1k_loss(arrival_rate, service_mean, service_scv, buffer_size):
    """
    M/G/1/k queue with finite buffer - compute loss probability.
    
    Args:
        arrival_rate: Arrival rate
        service_mean: Mean service time
        service_scv: Squared coefficient of variation of service times
        buffer_size: Buffer capacity
        
    Returns:
        dict: Performance metrics including loss probability
    """
    result = jpype.JPackage('jline').api.qsys.Qsys_mg1k_lossKt.qsys_mg1k_loss(
        float(arrival_rate), float(service_mean),
        float(service_scv), int(buffer_size)
    )
    return {
        'loss_prob': float(result.getLossProbability()),
        'queue_length': float(result.getQueueLength()),
        'utilization': float(result.getUtilization())
    }


def qsys_mmkk(arrival_rate, service_rate, num_servers, buffer_size):
    """
    M/M/k/k Erlang loss system.
    
    Args:
        arrival_rate: Arrival rate
        service_rate: Service rate per server
        num_servers: Number of servers
        buffer_size: System capacity (including servers)
        
    Returns:
        dict: Performance metrics including blocking probability
    """
    result = jpype.JPackage('jline').api.qsys.Qsys_mmkkKt.qsys_mmkk(
        float(arrival_rate), float(service_rate),
        int(num_servers), int(buffer_size)
    )
    return {
        'blocking_prob': float(result.getBlockingProbability()),
        'utilization': float(result.getUtilization()),
        'mean_customers': float(result.getMeanCustomers())
    }


def qsys_mmm(arrival_rate, service_rate, num_servers):
    """
    M/M/m infinite buffer multi-server queue.
    
    Args:
        arrival_rate: Arrival rate
        service_rate: Service rate per server
        num_servers: Number of servers
        
    Returns:
        dict: Performance metrics
    """
    result = jpype.JPackage('jline').api.qsys.Qsys_mmmKt.qsys_mmm(
        float(arrival_rate), float(service_rate), int(num_servers)
    )
    return {
        'utilization': float(result.getUtilization()),
        'queue_length': float(result.getQueueLength()),
        'waiting_time': float(result.getWaitingTime()),
        'response_time': float(result.getResponseTime())
    }


def qsys_mminf(arrival_rate, service_rate):
    """
    M/M/∞ infinite server queue.
    
    Args:
        arrival_rate: Arrival rate
        service_rate: Service rate per server
        
    Returns:
        dict: Performance metrics
    """
    result = jpype.JPackage('jline').api.qsys.Qsys_mminfKt.qsys_mminf(
        float(arrival_rate), float(service_rate)
    )
    return {
        'mean_customers': float(result.getMeanCustomers()),
        'variance_customers': float(result.getVarianceCustomers())
    }


def qsys_mginf(arrival_rate, service_mean, service_scv):
    """
    Analyzes M/G/∞ queue (infinite servers).

    Computes performance metrics for a queue with Poisson arrivals,
    general service times, and infinite servers.

    Args:
        arrival_rate: Arrival rate.
        service_mean: Mean service time.
        service_scv: Squared coefficient of variation of service time.

    Returns:
        dict: Performance metrics.
    """
    result = jpype.JPackage('jline').api.qsys.Qsys_mginfKt.qsys_mginf(
        float(arrival_rate), float(service_mean), float(service_scv)
    )
    return {
        'mean_customers': float(result.getMeanCustomers()),
        'variance_customers': float(result.getVarianceCustomers()),
        'throughput': float(result.getThroughput())
    }


def qsys_mm1k_loss(arrival_rate, service_rate, buffer_size):
    """
    Analyzes M/M/1/K queue with finite buffer.

    Computes performance metrics and loss probability for a single-server
    queue with Poisson arrivals, exponential service, and finite capacity.

    Args:
        arrival_rate: Arrival rate.
        service_rate: Service rate.
        buffer_size: Maximum system capacity (including server).

    Returns:
        dict: Performance metrics including loss probability.
    """
    result = jpype.JPackage('jline').api.qsys.Qsys_mm1k_lossKt.qsys_mm1k_loss(
        float(arrival_rate), float(service_rate), int(buffer_size)
    )
    return {
        'loss_probability': float(result.getLossProbability()),
        'mean_queue_length': float(result.getMeanQueueLength()),
        'mean_response_time': float(result.getMeanResponseTime()),
        'effective_arrival_rate': float(result.getEffectiveArrivalRate()),
        'utilization': float(result.getUtilization())
    }


def qsys_mg1k_loss_mgs(arrival_rate, service_mean, service_scv, buffer_size):
    """
    Analyzes M/G/1/K queue using modified generating function method.

    Computes loss probability and performance metrics for finite capacity
    queues with general service time distributions using the MGS
    (Modified Generating Function with Spectral) method.

    Args:
        arrival_rate: Arrival rate.
        service_mean: Mean service time.
        service_scv: Squared coefficient of variation of service time.
        buffer_size: Maximum system capacity.

    Returns:
        dict: Performance metrics including loss probability.
    """
    result = jpype.JPackage('jline').api.qsys.Qsys_mg1k_loss_mgsKt.qsys_mg1k_loss_mgs(
        float(arrival_rate), float(service_mean), float(service_scv), int(buffer_size)
    )
    return {
        'loss_probability': float(result.getLossProbability()),
        'mean_queue_length': float(result.getMeanQueueLength()),
        'mean_waiting_time': float(result.getMeanWaitingTime()),
        'effective_arrival_rate': float(result.getEffectiveArrivalRate())
    }


# =============================================================================
# M/G/1 Queue with Non-Preemptive Priorities
# =============================================================================

def qsys_mg1_prio(lambda_vals, mu_vals, cs_vals):
    """
    Analyze M/G/1 queue with non-preemptive (Head-of-Line) priorities.

    Computes per-class mean response times and overall utilization for a system
    with multiple priority classes using the Pollaczek-Khinchine formula extended
    for non-preemptive priority scheduling.

    For K priority classes (class 1 = highest priority), the mean response time
    for class k is:
        W_k = B_0 / ((1 - sum_{i=1}^{k-1} rho_i) * (1 - sum_{i=1}^{k} rho_i)) + 1/mu_k

    Where:
        - rho_i = lambda_i / mu_i (utilization of class i)
        - B_0 = sum_{i=1}^{K} lambda_i * E[S_i^2] / 2
        - E[S_i^2] = (1 + cs_i^2) / mu_i^2

    Args:
        lambda_vals: List or array of arrival rates per priority class.
                     Class 0 = highest priority.
        mu_vals: List or array of service rates per priority class.
        cs_vals: List or array of coefficients of variation per priority class.

    Returns:
        dict: Performance measures including:
            - W: Array of mean response times per priority class
            - rho: Overall system utilization (rhohat)

    Raises:
        ValueError: If inputs have different lengths or non-positive values.
        RuntimeError: If system is unstable (rho >= 1).

    Example:
        >>> # Two priority classes: high priority (class 0) and low priority (class 1)
        >>> result = qsys_mg1_prio([0.3, 0.2], [1.0, 0.8], [1.0, 1.5])
        >>> print(f"Response times: {result['W']}")
        >>> print(f"Utilization: {result['rho']:.4f}")
    """
    # Convert inputs to Java Matrix objects
    lambda_matrix = jlineMatrixFromArray(np.array(lambda_vals).flatten())
    mu_matrix = jlineMatrixFromArray(np.array(mu_vals).flatten())
    cs_matrix = jlineMatrixFromArray(np.array(cs_vals).flatten())

    result = jpype.JPackage('jline').api.qsys.Qsys_mg1_prioKt.qsys_mg1_prio(
        lambda_matrix, mu_matrix, cs_matrix
    )

    return {
        'W': jlineMatrixToArray(result.W).flatten(),
        'rho': float(result.rho)
    }


# =============================================================================
# MAP/PH Queue Analyzers (BUTools integration)
# =============================================================================

def qsys_mapg1(D0, D1, service_moments, num_ql_moms=3, num_ql_probs=100, num_st_moms=3):
    """
    Analyze MAP/G/1 queue (Markovian arrival process, general service).

    Uses BUTools MMAPPH1FCFS by fitting the general service time distribution
    to a Phase-Type distribution via moment matching.

    Args:
        D0: MAP hidden transition matrix (n x n). Can be list of lists or numpy array.
        D1: MAP arrival transition matrix (n x n). Can be list of lists or numpy array.
        service_moments: First k raw moments of service time [E[S], E[S^2], ...].
                        At least 2 moments required for good fitting.
        num_ql_moms: Number of queue length moments to compute (default 3).
        num_ql_probs: Number of queue length probabilities to compute (default 100).
        num_st_moms: Number of sojourn time moments to compute (default 3).

    Returns:
        dict: Performance metrics including:
            - mean_queue_length: Mean number of customers in system
            - mean_waiting_time: Mean waiting time in queue
            - mean_sojourn_time: Mean total time in system
            - utilization: Server utilization
            - queue_length_dist: Queue length distribution (if computed)
            - queue_length_moments: Raw moments of queue length
            - sojourn_time_moments: Raw moments of sojourn time
            - analyzer: Name of analyzer used

    Example:
        >>> # MMPP/G/1 queue
        >>> D0 = [[-2, 1], [1, -3]]  # Hidden transitions
        >>> D1 = [[1, 0], [0, 2]]    # Arrival transitions
        >>> moments = [1.0, 2.5, 10.0]  # First 3 moments of service time
        >>> result = qsys_mapg1(D0, D1, moments)
    """
    D0_matrix = jlineMatrixFromArray(np.array(D0, dtype=float))
    D1_matrix = jlineMatrixFromArray(np.array(D1, dtype=float))

    result = jpype.JPackage('jline').api.qsys.Qsys_mapg1Kt.qsys_mapg1(
        D0_matrix, D1_matrix,
        jpype.JArray(jpype.JDouble)(service_moments),
        jpype.JInt(num_ql_moms), jpype.JInt(num_ql_probs), jpype.JInt(num_st_moms)
    )

    return _parse_map_ph_result(result)


def qsys_mapg1_cv(D0, D1, mean_service, cv_service):
    """
    Analyze MAP/G/1 queue with service specified by mean and CV.

    Convenience function that converts mean and coefficient of variation
    to moments for the MAP/G/1 analysis.

    Args:
        D0: MAP hidden transition matrix.
        D1: MAP arrival transition matrix.
        mean_service: Mean service time.
        cv_service: Coefficient of variation of service time.

    Returns:
        dict: Performance metrics (see qsys_mapg1).
    """
    D0_matrix = jlineMatrixFromArray(np.array(D0, dtype=float))
    D1_matrix = jlineMatrixFromArray(np.array(D1, dtype=float))

    result = jpype.JPackage('jline').api.qsys.Qsys_mapg1Kt.qsys_mapg1(
        D0_matrix, D1_matrix,
        jpype.JDouble(mean_service), jpype.JDouble(cv_service)
    )

    return _parse_map_ph_result(result)


def qsys_mapmap1(C0, C1, D0, D1, num_ql_moms=3, num_ql_probs=100, num_st_moms=3):
    """
    Analyze MAP/MAP/1 queue (MAP arrivals, MAP service).

    Uses BUTools MMAPPH1FCFS by converting the service MAP to an equivalent
    Phase-Type representation.

    Args:
        C0: Arrival MAP hidden transition matrix (n x n).
        C1: Arrival MAP arrival transition matrix (n x n).
        D0: Service MAP hidden transition matrix (m x m).
        D1: Service MAP observable transition matrix (m x m).
        num_ql_moms: Number of queue length moments to compute (default 3).
        num_ql_probs: Number of queue length probabilities to compute (default 100).
        num_st_moms: Number of sojourn time moments to compute (default 3).

    Returns:
        dict: Performance metrics (see qsys_mapg1).

    Example:
        >>> # MAP/MAP/1 with bursty arrivals and correlated service
        >>> C0 = [[-2, 1], [1, -3]]
        >>> C1 = [[1, 0], [0, 2]]
        >>> D0 = [[-1.5, 0.5], [0.3, -1.2]]
        >>> D1 = [[1, 0], [0, 0.9]]
        >>> result = qsys_mapmap1(C0, C1, D0, D1)
    """
    C0_matrix = jlineMatrixFromArray(np.array(C0, dtype=float))
    C1_matrix = jlineMatrixFromArray(np.array(C1, dtype=float))
    D0_matrix = jlineMatrixFromArray(np.array(D0, dtype=float))
    D1_matrix = jlineMatrixFromArray(np.array(D1, dtype=float))

    result = jpype.JPackage('jline').api.qsys.Qsys_mapmap1Kt.qsys_mapmap1(
        C0_matrix, C1_matrix, D0_matrix, D1_matrix,
        jpype.JInt(num_ql_moms), jpype.JInt(num_ql_probs), jpype.JInt(num_st_moms)
    )

    return _parse_map_ph_result(result)


def qsys_mapmc(D0, D1, mu, c, mode="SylvesCR", max_num_comp=1000, verbose=0):
    """
    Analyze MAP/M/c queue (MAP arrivals, exponential service, c servers).

    Uses Q-MAM qCtMapMC which implements the Gaver-Jacobs-Latouche
    level-dependent QBD approach for MAP/M/c/FCFS queues.

    Args:
        D0: MAP hidden transition matrix (n x n).
        D1: MAP arrival transition matrix (n x n).
        mu: Exponential service rate.
        c: Number of servers.
        mode: Solver mode - one of "SylvesCR", "DirectCR", "SylvesFI", "DirectFI".
        max_num_comp: Maximum number of queue length components (default 1000).
        verbose: Verbosity level (default 0).

    Returns:
        dict: Performance metrics (see qsys_mapg1).

    Example:
        >>> # MMPP arrivals to a 3-server queue
        >>> D0 = [[-2, 1], [1, -3]]
        >>> D1 = [[1, 0], [0, 2]]
        >>> result = qsys_mapmc(D0, D1, mu=2.0, c=3)
    """
    D0_matrix = jlineMatrixFromArray(np.array(D0, dtype=float))
    D1_matrix = jlineMatrixFromArray(np.array(D1, dtype=float))

    result = jpype.JPackage('jline').api.qsys.Qsys_mapmcKt.qsys_mapmc(
        D0_matrix, D1_matrix,
        jpype.JDouble(mu), jpype.JInt(c),
        str(mode), jpype.JInt(max_num_comp), jpype.JInt(verbose)
    )

    return _parse_map_ph_result(result)


def qsys_mapm1(D0, D1, mu):
    """
    Analyze MAP/M/1 queue (MAP arrivals, exponential service).

    Convenience function for single-server case with exponential service.
    Delegates to BUTools via qsys_mapph1 for better performance metrics.

    Args:
        D0: MAP hidden transition matrix.
        D1: MAP arrival transition matrix.
        mu: Exponential service rate.

    Returns:
        dict: Performance metrics (see qsys_mapg1).
    """
    D0_matrix = jlineMatrixFromArray(np.array(D0, dtype=float))
    D1_matrix = jlineMatrixFromArray(np.array(D1, dtype=float))

    result = jpype.JPackage('jline').api.qsys.Qsys_mapmcKt.qsys_mapm1(
        D0_matrix, D1_matrix, jpype.JDouble(mu)
    )

    return _parse_map_ph_result(result)


def qsys_mapdc(D0, D1, s, c, max_num_comp=1000, num_steps=1, verbose=0):
    """
    Analyze MAP/D/c queue (MAP arrivals, deterministic service, c servers).

    Uses Q-MAM qCtMapDC which implements Non-Skip-Free (NSF) Markov chain
    analysis for MAP/D/c/FCFS queues.

    Args:
        D0: MAP hidden transition matrix (n x n).
        D1: MAP arrival transition matrix (n x n).
        s: Deterministic service time (positive scalar).
        c: Number of servers.
        max_num_comp: Maximum number of queue length components (default 1000).
        num_steps: Number of waiting time distribution points per service interval (default 1).
        verbose: Verbosity level (default 0).

    Returns:
        dict: Performance metrics including:
            - mean_queue_length: Mean number of customers in system
            - mean_waiting_time: Mean waiting time in queue
            - mean_sojourn_time: Mean sojourn time (waiting + service)
            - utilization: Server utilization (per server)
            - queue_length_dist: Queue length distribution P(Q=n)
            - waiting_time_dist: Waiting time CDF at discrete points
            - analyzer: Analyzer identifier

    Example:
        >>> # MMPP arrivals to a 2-server deterministic queue
        >>> D0 = [[-2, 1], [1, -3]]
        >>> D1 = [[1, 0], [0, 2]]
        >>> result = qsys_mapdc(D0, D1, s=0.5, c=2)
    """
    D0_matrix = jlineMatrixFromArray(np.array(D0, dtype=float))
    D1_matrix = jlineMatrixFromArray(np.array(D1, dtype=float))

    result = jpype.JPackage('jline').api.qsys.Qsys_mapdcKt.qsys_mapdc(
        D0_matrix, D1_matrix,
        jpype.JDouble(s), jpype.JInt(c),
        jpype.JInt(max_num_comp), jpype.JInt(num_steps), jpype.JInt(verbose)
    )

    return _parse_map_dc_result(result)


def qsys_mapd1(D0, D1, s, max_num_comp=1000, num_steps=1):
    """
    Analyze MAP/D/1 queue (MAP arrivals, deterministic service).

    Convenience function for single-server case with deterministic service.

    Args:
        D0: MAP hidden transition matrix.
        D1: MAP arrival transition matrix.
        s: Deterministic service time.
        max_num_comp: Maximum number of queue length components (default 1000).
        num_steps: Number of waiting time distribution points per service interval (default 1).

    Returns:
        dict: Performance metrics (see qsys_mapdc).
    """
    return qsys_mapdc(D0, D1, s, 1, max_num_comp, num_steps)


def qsys_phdc(alpha, T, s, c):
    """
    Analyze PH/D/c queue (Phase-Type arrivals, deterministic service, c servers).

    Converts Phase-Type arrival process to equivalent MAP representation
    and delegates to qsys_mapdc.

    The conversion is:
        D0 = T (generator matrix)
        D1 = (-T * e) * alpha (exit rate outer product with initial vector)

    Args:
        alpha: PH arrival initial probability vector (1 x n).
        T: PH arrival generator matrix (n x n).
        s: Deterministic service time (positive scalar).
        c: Number of servers.

    Returns:
        dict: Performance metrics (see qsys_mapdc).

    Example:
        >>> # Erlang-2 arrivals to a 2-server deterministic queue
        >>> alpha = [1.0, 0.0]  # Start in phase 1
        >>> T = [[-2, 2], [0, -2]]  # Erlang-2 with rate 2
        >>> result = qsys_phdc(alpha, T, s=0.5, c=2)
    """
    alpha_matrix = jlineMatrixFromArray(np.array(alpha, dtype=float).reshape(1, -1))
    T_matrix = jlineMatrixFromArray(np.array(T, dtype=float))

    result = jpype.JPackage('jline').api.qsys.Qsys_mapdcKt.qsys_phdc(
        alpha_matrix, T_matrix,
        jpype.JDouble(s), jpype.JInt(c)
    )

    return _parse_map_dc_result(result)


def qsys_mapph1(D0, D1, sigma, S, num_ql_moms=3, num_ql_probs=100, num_st_moms=3):
    """
    Analyze MAP/PH/1 queue (MAP arrivals, Phase-Type service).

    Uses BUTools MMAPPH1FCFS which computes performance measures for a
    single-server FCFS queue with Markovian Arrival Process arrivals
    and Phase-Type service times.

    Args:
        D0: MAP hidden transition matrix (n x n).
        D1: MAP arrival transition matrix (n x n).
        sigma: PH service initial probability vector (1 x m).
        S: PH service generator matrix (m x m).
        num_ql_moms: Number of queue length moments to compute (default 3).
        num_ql_probs: Number of queue length probabilities to compute (default 100).
        num_st_moms: Number of sojourn time moments to compute (default 3).

    Returns:
        dict: Performance metrics (see qsys_mapg1).

    Example:
        >>> # MAP arrivals with Erlang-2 service
        >>> D0 = [[-2, 1], [1, -3]]
        >>> D1 = [[1, 0], [0, 2]]
        >>> sigma = [1.0, 0.0]  # Start in phase 1
        >>> S = [[-2, 2], [0, -2]]  # Erlang-2 with rate 2
        >>> result = qsys_mapph1(D0, D1, sigma, S)
    """
    D0_matrix = jlineMatrixFromArray(np.array(D0, dtype=float))
    D1_matrix = jlineMatrixFromArray(np.array(D1, dtype=float))
    sigma_matrix = jlineMatrixFromArray(np.array(sigma, dtype=float).reshape(1, -1))
    S_matrix = jlineMatrixFromArray(np.array(S, dtype=float))

    result = jpype.JPackage('jline').api.qsys.Qsys_mapph1Kt.qsys_mapph1(
        D0_matrix, D1_matrix, sigma_matrix, S_matrix,
        jpype.JInt(num_ql_moms), jpype.JInt(num_ql_probs), jpype.JInt(num_st_moms)
    )

    return _parse_map_ph_result(result)


def qsys_phph1(alpha, T, beta, S, num_ql_moms=3, num_ql_probs=100, num_st_moms=3):
    """
    Analyze PH/PH/1 queue (Phase-Type arrivals and service).

    Uses BUTools MMAPPH1FCFS by converting the arrival PH to an equivalent
    Markovian Arrival Process.

    Args:
        alpha: Arrival PH initial probability vector (1 x n).
        T: Arrival PH generator matrix (n x n).
        beta: Service PH initial probability vector (1 x m).
        S: Service PH generator matrix (m x m).
        num_ql_moms: Number of queue length moments to compute (default 3).
        num_ql_probs: Number of queue length probabilities to compute (default 100).
        num_st_moms: Number of sojourn time moments to compute (default 3).

    Returns:
        dict: Performance metrics (see qsys_mapg1).

    Example:
        >>> # Erlang-2 arrivals, hyperexponential service
        >>> alpha = [1.0, 0.0]
        >>> T = [[-1, 1], [0, -1]]  # Erlang-2 with rate 1
        >>> beta = [0.6, 0.4]
        >>> S = [[-2, 0], [0, -0.5]]  # Hyperexp
        >>> result = qsys_phph1(alpha, T, beta, S)
    """
    alpha_matrix = jlineMatrixFromArray(np.array(alpha, dtype=float).reshape(1, -1))
    T_matrix = jlineMatrixFromArray(np.array(T, dtype=float))
    beta_matrix = jlineMatrixFromArray(np.array(beta, dtype=float).reshape(1, -1))
    S_matrix = jlineMatrixFromArray(np.array(S, dtype=float))

    result = jpype.JPackage('jline').api.qsys.Qsys_phph1Kt.qsys_phph1(
        alpha_matrix, T_matrix, beta_matrix, S_matrix,
        jpype.JInt(num_ql_moms), jpype.JInt(num_ql_probs), jpype.JInt(num_st_moms)
    )

    return _parse_map_ph_result(result)


def qsys_phm1(alpha, T, mu):
    """
    Analyze PH/M/1 queue (Phase-Type arrivals, exponential service).

    Convenience function for queues with PH renewal arrivals and
    exponential service.

    Args:
        alpha: Arrival PH initial probability vector.
        T: Arrival PH generator matrix.
        mu: Exponential service rate.

    Returns:
        dict: Performance metrics (see qsys_mapg1).
    """
    alpha_matrix = jlineMatrixFromArray(np.array(alpha, dtype=float).reshape(1, -1))
    T_matrix = jlineMatrixFromArray(np.array(T, dtype=float))

    result = jpype.JPackage('jline').api.qsys.Qsys_phph1Kt.qsys_phm1(
        alpha_matrix, T_matrix, jpype.JDouble(mu)
    )

    return _parse_map_ph_result(result)


def qsys_phmc(alpha, T, mu, c):
    """
    Analyze PH/M/c queue (Phase-Type arrivals, exponential service, c servers).

    Converts PH arrival process to equivalent MAP and uses Q-MAM.

    Args:
        alpha: Arrival PH initial probability vector.
        T: Arrival PH generator matrix.
        mu: Exponential service rate per server.
        c: Number of servers.

    Returns:
        dict: Performance metrics (see qsys_mapg1).
    """
    alpha_matrix = jlineMatrixFromArray(np.array(alpha, dtype=float).reshape(1, -1))
    T_matrix = jlineMatrixFromArray(np.array(T, dtype=float))

    result = jpype.JPackage('jline').api.qsys.Qsys_mapmcKt.qsys_phmc(
        alpha_matrix, T_matrix, jpype.JDouble(mu), jpype.JInt(c)
    )

    return _parse_map_ph_result(result)


def _parse_map_ph_result(result):
    """
    Parse QsysMapPhResult Java object into Python dict.

    Args:
        result: Java QsysMapPhResult object.

    Returns:
        dict: Parsed performance metrics.
    """
    response = {
        'mean_queue_length': float(result.getMeanQueueLength()),
        'mean_waiting_time': float(result.getMeanWaitingTime()),
        'mean_sojourn_time': float(result.getMeanSojournTime()),
        'utilization': float(result.getUtilization()),
        'analyzer': str(result.getAnalyzer())
    }

    # Optional fields - may be null
    if result.getQueueLengthDist() is not None:
        response['queue_length_dist'] = jlineMatrixToArray(result.getQueueLengthDist()).flatten()
    if result.getQueueLengthMoments() is not None:
        response['queue_length_moments'] = jlineMatrixToArray(result.getQueueLengthMoments()).flatten()
    if result.getSojournTimeMoments() is not None:
        response['sojourn_time_moments'] = jlineMatrixToArray(result.getSojournTimeMoments()).flatten()

    return response


def _parse_map_dc_result(result):
    """
    Parse QsysMapDcResult Java object into Python dict.

    Args:
        result: Java QsysMapDcResult object.

    Returns:
        dict: Parsed performance metrics for MAP/D/c queue.
    """
    response = {
        'mean_queue_length': float(result.getMeanQueueLength()),
        'mean_waiting_time': float(result.getMeanWaitingTime()),
        'mean_sojourn_time': float(result.getMeanSojournTime()),
        'utilization': float(result.getUtilization()),
        'analyzer': str(result.getAnalyzer())
    }

    # Queue length and waiting time distributions
    if result.getQueueLengthDist() is not None:
        response['queue_length_dist'] = jlineMatrixToArray(result.getQueueLengthDist()).flatten()
    if result.getWaitingTimeDist() is not None:
        response['waiting_time_dist'] = jlineMatrixToArray(result.getWaitingTimeDist()).flatten()

    return response


# =============================================================================
# MX/M/1 Queue (Batch Markovian Arrivals, Exponential Service)
# =============================================================================

def qsys_mxm1(lambda_batch, mu, mean_batch_size, second_moment_batch_size):
    """
    Analyze MX/M/1 queue with batch Markovian arrivals and exponential service.

    The MX/M/1 queue models a system with:
    - Batch Poisson arrivals at rate lambda_batch
    - Random batch sizes with mean E[X] and second moment E[X²]
    - Single exponential server with rate mu

    The formula accounts for both queueing delay and internal batch delay:
        Wq = rho/(mu*(1-rho)) + (E[X²] - E[X])/(2*mu*E[X]*(1-rho))

    Args:
        lambda_batch (float): Batch arrival rate (batches per unit time)
        mu (float): Service rate (jobs per unit time)
        mean_batch_size (float): Mean batch size E[X]
        second_moment_batch_size (float): Second moment of batch size E[X²]

    Returns:
        dict: Performance metrics with keys:
            - 'W': Mean time in system
            - 'Wq': Mean waiting time in queue
            - 'U': Server utilization
            - 'Q': Mean queue length (including service)
            - 'lambda': Effective job arrival rate (lambda_batch * E[X])

    Raises:
        RuntimeError: If system is unstable (rho >= 1)

    Examples:
        >>> # Deterministic batch size of 3
        >>> result = qsys_mxm1(lambda_batch=1.0, mu=5.0,
        ...                     mean_batch_size=3.0,
        ...                     second_moment_batch_size=9.0)
        >>> print(f"Mean time in system: {result['W']:.4f}")

        >>> # Variable batch sizes with mean=2, variance=1
        >>> E_X = 2.0
        >>> Var_X = 1.0
        >>> E_X2 = Var_X + E_X**2
        >>> result = qsys_mxm1(lambda_batch=1.5, mu=8.0,
        ...                     mean_batch_size=E_X,
        ...                     second_moment_batch_size=E_X2)
    """
    # Call JAR implementation
    java_result = jpype.JPackage('jline').api.qsys.Qsys_mxm1.qsys_mxm1(
        jpype.JDouble(lambda_batch),
        jpype.JDouble(mu),
        jpype.JDouble(mean_batch_size),
        jpype.JDouble(second_moment_batch_size)
    )

    # Extract results from Java Matrix
    # JAR returns Matrix with [W, Wq, U, Q]
    W = float(java_result.get(0))
    Wq = float(java_result.get(1))
    U = float(java_result.get(2))
    Q = float(java_result.get(3))

    # Compute effective job arrival rate
    lambda_eff = lambda_batch * mean_batch_size

    return {
        'W': W,
        'Wq': Wq,
        'U': U,
        'Q': Q,
        'lambda': lambda_eff
    }


# =============================================================================
# M/G/1 Queue with Size-Based Preemptive Scheduling
# =============================================================================

def qsys_mg1_srpt(lambda_vals, mu_vals, cs_vals):
    """
    Analyze M/G/1 queue with SRPT (Shortest Remaining Processing Time) scheduling.

    SRPT gives priority to jobs with the shortest remaining processing time.
    It is optimal for minimizing mean response time.

    For M/G/1/SRPT, the mean response time for a job of size x is:
        E[T(x)] = integral from 0 to x of dt/(1-rho(t))^2 + x/(1-rho(x))

    where rho(x) = lambda * integral from 0 to x of t*f(t)dt

    SRPT is "Sometimes Unfair": fair when rho <= 0.5, unfair for higher loads.

    Args:
        lambda_vals: List or array of arrival rates per class.
        mu_vals: List or array of service rates per class (1/mean service time).
        cs_vals: List or array of coefficients of variation per class.

    Returns:
        dict: Performance measures including:
            - W: Array of mean response times per class
            - rho: Overall system utilization

    Reference:
        A. Wierman and M. Harchol-Balter, "Classifying scheduling policies with
        respect to unfairness in an M/GI/1", SIGMETRICS 2003.
    """
    lambda_matrix = jlineMatrixFromArray(np.array(lambda_vals).flatten())
    mu_matrix = jlineMatrixFromArray(np.array(mu_vals).flatten())
    cs_matrix = jlineMatrixFromArray(np.array(cs_vals).flatten())

    result = jpype.JPackage('jline').api.qsys.Qsys_mg1_srptKt.qsys_mg1_srpt(
        lambda_matrix, mu_matrix, cs_matrix
    )

    return {
        'W': jlineMatrixToArray(result.W).flatten(),
        'rho': float(result.rho)
    }


def qsys_mg1_psjf(lambda_vals, mu_vals, cs_vals):
    """
    Analyze M/G/1 queue with PSJF (Preemptive Shortest Job First) scheduling.

    PSJF gives priority based on *original* job size (not remaining time).
    Jobs with smaller original sizes always preempt jobs with larger sizes.

    For PSJF, the mean response time for a job of size x is:
        E[T(x)] = (lambda * integral from 0 to x of t^2*f(t)dt) / (2*(1-rho(x))^2)
                  + x / (1 - rho(x))

    PSJF is "Always Unfair": some job sizes are treated unfairly under all loads.

    Args:
        lambda_vals: List or array of arrival rates per class.
        mu_vals: List or array of service rates per class (1/mean service time).
        cs_vals: List or array of coefficients of variation per class.

    Returns:
        dict: Performance measures including:
            - W: Array of mean response times per class
            - rho: Overall system utilization

    Reference:
        A. Wierman and M. Harchol-Balter, "Classifying scheduling policies with
        respect to unfairness in an M/GI/1", SIGMETRICS 2003.
    """
    lambda_matrix = jlineMatrixFromArray(np.array(lambda_vals).flatten())
    mu_matrix = jlineMatrixFromArray(np.array(mu_vals).flatten())
    cs_matrix = jlineMatrixFromArray(np.array(cs_vals).flatten())

    result = jpype.JPackage('jline').api.qsys.Qsys_mg1_psjfKt.qsys_mg1_psjf(
        lambda_matrix, mu_matrix, cs_matrix
    )

    return {
        'W': jlineMatrixToArray(result.W).flatten(),
        'rho': float(result.rho)
    }


def qsys_mg1_fb(lambda_vals, mu_vals, cs_vals):
    """
    Analyze M/G/1 queue with FB/LAS (Feedback/Least Attained Service) scheduling.

    FB/LAS gives priority to jobs with the least attained service (smallest age).
    This policy requires no knowledge of job sizes - only how much service
    a job has received.

    For FB, the mean response time for a job of size x is:
        E[T(x)] = (lambda * integral from 0 to x of t*F_bar(t)dt) / (1-rho_x)^2
                  + x / (1 - rho_x)

    FB is "Always Unfair" but approximates SRPT for heavy-tailed distributions.

    Args:
        lambda_vals: List or array of arrival rates per class.
        mu_vals: List or array of service rates per class (1/mean service time).
        cs_vals: List or array of coefficients of variation per class.

    Returns:
        dict: Performance measures including:
            - W: Array of mean response times per class
            - rho: Overall system utilization

    Reference:
        A. Wierman and M. Harchol-Balter, "Classifying scheduling policies with
        respect to unfairness in an M/GI/1", SIGMETRICS 2003.
    """
    lambda_matrix = jlineMatrixFromArray(np.array(lambda_vals).flatten())
    mu_matrix = jlineMatrixFromArray(np.array(mu_vals).flatten())
    cs_matrix = jlineMatrixFromArray(np.array(cs_vals).flatten())

    result = jpype.JPackage('jline').api.qsys.Qsys_mg1_fbKt.qsys_mg1_fb(
        lambda_matrix, mu_matrix, cs_matrix
    )

    return {
        'W': jlineMatrixToArray(result.W).flatten(),
        'rho': float(result.rho)
    }


def qsys_mg1_lrpt(lambda_vals, mu_vals, cs_vals):
    """
    Analyze M/G/1 queue with LRPT (Longest Remaining Processing Time) scheduling.

    LRPT gives priority to jobs with the longest remaining processing time.
    This prioritizes large jobs over small jobs.

    For LRPT, the expected response time for a job of size x is:
        E[T(x)] = x/(1-rho) + lambda*E[X^2] / (2*(1-rho)^2)

    The slowdown E[S(x)] converges to 1/(1-rho) as x -> infinity.

    LRPT is "Always Unfair": E[S(y)] > 1/(1-rho) for all finite job sizes y.

    Args:
        lambda_vals: List or array of arrival rates per class.
        mu_vals: List or array of service rates per class (1/mean service time).
        cs_vals: List or array of coefficients of variation per class.

    Returns:
        dict: Performance measures including:
            - W: Array of mean response times per class
            - rho: Overall system utilization

    Reference:
        A. Wierman and M. Harchol-Balter, "Classifying scheduling policies with
        respect to unfairness in an M/GI/1", SIGMETRICS 2003.
    """
    lambda_matrix = jlineMatrixFromArray(np.array(lambda_vals).flatten())
    mu_matrix = jlineMatrixFromArray(np.array(mu_vals).flatten())
    cs_matrix = jlineMatrixFromArray(np.array(cs_vals).flatten())

    result = jpype.JPackage('jline').api.qsys.Qsys_mg1_lrptKt.qsys_mg1_lrpt(
        lambda_matrix, mu_matrix, cs_matrix
    )

    return {
        'W': jlineMatrixToArray(result.W).flatten(),
        'rho': float(result.rho)
    }
