"""
Decomposition-aggregation (DA) toolkit.

Shared primitives for solution methods that decompose a model into isolated
submodels, exchange flows or rates, and iterate to a fixed point.
"""

from .fpi import da_fpi
from .superpos import da_traffic_superpos
from .cacheqn_retrieval import da_cacheqn_retrieval

__all__ = ['da_fpi', 'da_traffic_superpos', 'da_cacheqn_retrieval']
