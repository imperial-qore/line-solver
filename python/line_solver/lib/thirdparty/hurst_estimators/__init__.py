"""
Hurst parameter estimators for long-range dependent time series.

Ported from MATLAB to Python. Original implementation by Chu Chen
(Version 1.0, 03/10/2008, chen-chu@163.com).

Reference: Murad S. Taqqu, Vadim Teverovsky and Walter Willinger,
"Estimators for long-range dependence: an empirical study".
"""

from .hurst_estimate import hurst_estimate
from .rs import rs
from .absval import absval
from .aggvar import aggvar
from .boxper import boxper
from .diffvar import diffvar
from .higuchi import higuchi
from .peng import peng
from .per import per

__all__ = [
    'hurst_estimate',
    'rs',
    'absval',
    'aggvar',
    'boxper',
    'diffvar',
    'higuchi',
    'peng',
    'per',
]
