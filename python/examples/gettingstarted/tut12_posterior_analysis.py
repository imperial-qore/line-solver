"""
Tutorial 12: Posterior analysis with uncertain parameters

This tutorial demonstrates how to analyze queueing models with parameter
uncertainty using the Posterior solver. The Posterior solver works with
Prior distributions that represent uncertainty about model parameters,
computing posterior distributions of performance metrics.

Scenario: An M/M/1 queue where the service rate is uncertain. We model
this uncertainty using a Prior distribution with 30 alternatives,
creating a Gaussian-like distribution of possible service rates.

Copyright (c) 2012-2025, Imperial College London
All rights reserved.
"""

from line_solver import *
import numpy as np

# Block 1: Create model with uncertain service rate
model = Network('UncertainServiceModel')

# Create nodes
source = Source(model, 'Source')
queue = Queue(model, 'Queue', SchedStrategy.FCFS)
sink = Sink(model, 'Sink')

# Create job class
job_class = OpenClass(model, 'Jobs')

# Set arrival rate (lambda = 0.5)
arrival_rate = 0.5
source.setArrival(job_class, Exp(arrival_rate))

# Block 2: Define Prior distribution for uncertain service rate
# Use many alternatives to create a smooth, continuous-looking PDF
# Service rates range from 0.7 to 2.5 with a Gaussian-like prior
num_alternatives = 30
service_rates = np.linspace(0.7, 2.5, num_alternatives)

# Create Gaussian-like prior probabilities centered at mu=1.3
prior_mean = 1.3
prior_std = 0.4
prior_probs = np.exp(-0.5 * ((service_rates - prior_mean) / prior_std)**2)
prior_probs = prior_probs / np.sum(prior_probs)  # Normalize to sum to 1

alternatives = [Exp(rate) for rate in service_rates]
prior = Prior(alternatives, prior_probs.tolist())
queue.setService(job_class, prior)

# Block 3: Complete model topology
model.link(Network.serial_routing([source, queue, sink]))

print('Model: M/M/1 with uncertain service rate')
print(f'Arrival rate: lambda = {arrival_rate:.1f}')
print(f'Number of service rate alternatives: {num_alternatives}')
print(f'Service rate range: mu in [{service_rates.min():.2f}, {service_rates.max():.2f}]')
print(f'Prior: Gaussian-like centered at mu={prior_mean:.1f} with std={prior_std:.1f}\n')

# Block 4: Solve with Posterior wrapper using MVA
post = Posterior(model, MVA)
post.run_analyzer()

# Block 5: Get prior-weighted average results
avg_table = post.get_avg_table()
print('Prior-weighted average performance metrics:')
print(avg_table)
print()

# Block 6: Get posterior table with per-alternative results
post_table = post.get_posterior_table()
print('Posterior table (showing per-alternative results):')
print(post_table)
print()

# Block 7: Extract posterior distributions for different metrics
# Get posterior distribution of response time at the queue
resp_dist = post.get_posterior_dist('R', queue, job_class)

# Get posterior distribution of queue length at the queue
qlen_dist = post.get_posterior_dist('Q', queue, job_class)

# Get posterior distribution of utilization at the queue
util_dist = post.get_posterior_dist('U', queue, job_class)

# Block 8: Extract values and probabilities from the EmpiricalCDF objects
# For response time
resp_cdf_data = resp_dist.data  # [CDF, Value]
resp_values = resp_cdf_data[:, 1]
resp_cdf = resp_cdf_data[:, 0]
resp_probs = np.concatenate([[resp_cdf[0]], np.diff(resp_cdf)])  # Convert CDF to PMF

# For queue length
qlen_cdf_data = qlen_dist.data
qlen_values = qlen_cdf_data[:, 1]
qlen_cdf = qlen_cdf_data[:, 0]
qlen_probs = np.concatenate([[qlen_cdf[0]], np.diff(qlen_cdf)])

# For utilization
util_cdf_data = util_dist.data
util_values = util_cdf_data[:, 1]
util_cdf = util_cdf_data[:, 0]
util_probs = np.concatenate([[util_cdf[0]], np.diff(util_cdf)])

# Block 9: Print posterior distribution statistics
# Response time statistics
expected_R = np.sum(resp_values * resp_probs)
mode_idx_R = np.argmax(resp_probs)
print('Response Time (R) at Queue:')
print(f'  Expected Value E[R]: {expected_R:.4f}')
print(f'  Mode (most likely): {resp_values[mode_idx_R]:.4f}\n')

# Queue length statistics
expected_Q = np.sum(qlen_values * qlen_probs)
mode_idx_Q = np.argmax(qlen_probs)
print('Queue Length (Q) at Queue:')
print(f'  Expected Value E[Q]: {expected_Q:.4f}')
print(f'  Mode (most likely): {qlen_values[mode_idx_Q]:.4f}\n')

# Utilization statistics
expected_U = np.sum(util_values * util_probs)
mode_idx_U = np.argmax(util_probs)
print('Utilization (U) at Queue:')
print(f'  Expected Value E[U]: {expected_U:.4f}')
print(f'  Mode (most likely): {util_values[mode_idx_U]:.4f}\n')

# Optional: Find median from CDF
median_idx = np.where(resp_cdf >= 0.5)[0]
if len(median_idx) > 0:
    median_R = resp_values[median_idx[0]]
    print(f'Response Time Median: {median_R:.4f}')
