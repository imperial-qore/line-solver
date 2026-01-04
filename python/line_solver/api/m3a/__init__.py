"""
M3A: Markovian Arrival Process with 3-moment Approximation.

Native Python implementations of M3A compression algorithms for MMAPs.
"""

from .utils import (
    compute_autocorrelation,
    compute_idc,
    compute_coeff_var,
    compute_moments,
    compute_spectral_gap,
    validate_mmap,
    compute_kl_divergence,
    sample_interarrival_times,
)

from .compressor import (
    M3aCompressor,
    compress_mmap,
)

from .fit import (
    m3a_fit,
    m3a_fit_from_trace,
)

from .m3pp import (
    m3pp_rand,
    m3pp2m_interleave,
    m3pp2m_fitc_approx_ag_multiclass,
    m3pp2m_fitc_approx,
    m3pp2m_fitc_theoretical,
    m3pp_superpos,
)

__all__ = [
    # Utils
    'compute_autocorrelation',
    'compute_idc',
    'compute_coeff_var',
    'compute_moments',
    'compute_spectral_gap',
    'validate_mmap',
    'compute_kl_divergence',
    'sample_interarrival_times',
    # Compressor
    'M3aCompressor',
    'compress_mmap',
    # Fit
    'm3a_fit',
    'm3a_fit_from_trace',
    # M3PP functions
    'm3pp_rand',
    'm3pp2m_interleave',
    'm3pp2m_fitc_approx_ag_multiclass',
    'm3pp2m_fitc_approx',
    'm3pp2m_fitc_theoretical',
    'm3pp_superpos',
]
