"""
Statistical inference on measured data: fitting and testing model assumptions.
"""

from .nhpp_ks import infer_nhpp_ks

__all__ = [
    'infer_nhpp_ks',
]
