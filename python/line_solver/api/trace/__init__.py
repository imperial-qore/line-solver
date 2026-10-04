"""
Trace Analysis Functions.

Native Python implementations for statistical analysis of empirical trace data
including means, variances, correlations, and index of dispersion.

Key functions:
    trace_mean: Mean of trace data
    trace_var: Variance of trace data
    trace_scv: Squared coefficient of variation
    trace_acf: Autocorrelation function
    trace_summary: Comprehensive summary statistics
"""

from .trace_analysis import (
    # Single trace analysis
    trace_mean,
    trace_var,
    trace_scv,
    trace_acf,
    trace_gamma,
    trace_iat2counts,
    trace_idi,
    trace_idc,
    trace_pmf,
    trace_shuffle,
    trace_joint,
    trace_iat2bins,
    trace_summary,
    trace_bicov,
    # Multi-class trace analysis
    mtrace_mean,
    mtrace_var,
    mtrace_count,
    mtrace_sigma,
    mtrace_sigma2,
    mtrace_cross_moment,
    mtrace_forward_moment,
    mtrace_backward_moment,
    mtrace_cov,
    mtrace_pc,
    mtrace_summary,
    mtrace_split,
    mtrace_merge,
    mtrace_joint,
    mtrace_moment,
    mtrace_moment_simple,
    mtrace_bootstrap,
    mtrace_iat2counts,
)

__all__ = [
    # Single trace analysis
    'trace_mean',
    'trace_var',
    'trace_scv',
    'trace_acf',
    'trace_gamma',
    'trace_iat2counts',
    'trace_idi',
    'trace_idc',
    'trace_pmf',
    'trace_shuffle',
    'trace_joint',
    'trace_iat2bins',
    'trace_summary',
    'trace_bicov',
    # Multi-class trace analysis
    'mtrace_mean',
    'mtrace_var',
    'mtrace_count',
    'mtrace_sigma',
    'mtrace_sigma2',
    'mtrace_cross_moment',
    'mtrace_forward_moment',
    'mtrace_backward_moment',
    'mtrace_cov',
    'mtrace_pc',
    'mtrace_summary',
    'mtrace_split',
    'mtrace_merge',
    'mtrace_joint',
    'mtrace_moment',
    'mtrace_moment_simple',
    'mtrace_bootstrap',
    'mtrace_iat2counts',
]
