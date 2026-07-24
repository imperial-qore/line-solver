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

from .amap2 import (
    amap2_assemble,
    amap2_adjust_gamma,
    amap2_fit_gamma,
    amap2_fit_gamma_map,
    amap2_fit_gamma_trace,
    amap2_fitall_gamma,
)

from .maph2m import (
    maph2m_fit,
    maph2m_fit_mmap,
    maph2m_fit_multiclass,
    maph2m_fit_trace,
)

from .mmap3k import (
    marking_block,
    marking_orders,
    mmap3k_fit,
)

from .mmap2k import (
    marking_inverse,
    mmap2k_fit,
)

from .mamap2m import (
    mamap2m_fit,
    mamap2m_fit_fb_multiclass,
    mamap2m_fit_gamma_fb,
    mamap2m_fit_mmap,
    mamap2m_fit_trace,
)

from .m3pp import (
    m3pp_rand,
    m3pp2m_interleave,
    m3pp2m_fitc_approx_ag_multiclass,
    m3pp2m_fitc_approx,
    m3pp2m_fitc,
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
    'm3pp2m_fitc',
    'm3pp2m_fitc_theoretical',
    'm3pp_superpos',

    'amap2_assemble',
    'amap2_adjust_gamma',
    'amap2_fit_gamma',
    'amap2_fit_gamma_map',
    'amap2_fit_gamma_trace',
    'amap2_fitall_gamma',
    'maph2m_fit',
    'maph2m_fit_mmap',
    'maph2m_fit_multiclass',
    'maph2m_fit_trace',
    'mmap3k_fit',
    'marking_block',
    'marking_orders',
    'mmap2k_fit',
    'marking_inverse',
    'mamap2m_fit',
    'mamap2m_fit_fb_multiclass',
    'mamap2m_fit_gamma_fb',
    'mamap2m_fit_mmap',
    'mamap2m_fit_trace',
]
