"""
Main dispatcher for Hurst parameter estimation methods.

Ported from MATLAB implementation by Chu Chen (Version 1.0, 03/10/2008).

Reference: Murad S. Taqqu, Vadim Teverovsky and Walter Willinger,
"Estimators for long-range dependence: an empirical study".
"""

from .rs import rs
from .absval import absval
from .aggvar import aggvar
from .boxper import boxper
from .diffvar import diffvar
from .higuchi import higuchi
from .peng import peng
from .per import per

_METHODS = {
    'rs': rs,
    'RS': rs,
    'absval': absval,
    'aggvar': aggvar,
    'boxper': boxper,
    'diffvar': diffvar,
    'higuchi': higuchi,
    'peng': peng,
    'per': per,
}


def hurst_estimate(sequence, method, isplot=False, opt=None):
    """
    Estimate the Hurst parameter of a given sequence with an appointed method.

    The algorithms of the methods can be found in Murad S. Taqqu, Vadim
    Teverovsky and Walter Willinger's paper "Estimators for long-range
    dependence: an empirical study" or other related papers.

    Parameters
    ----------
    sequence : array_like
        The input sequence for estimation.
    method : str
        The name of the estimation method. One of:
        'aggvar', 'RS', 'rs', 'per', 'absval', 'boxper', 'diffvar',
        'higuchi', 'peng'.
    isplot : bool, optional
        Whether to display the plot (default False).
    opt : optional
        An optional parameter for some methods (e.g., ``moment`` for absval,
        ``boxnumber`` for boxper).

    Returns
    -------
    H : float
        The estimated Hurst parameter of the input sequence.

    Examples
    --------
    >>> import numpy as np
    >>> H = hurst_estimate(np.random.randn(10000), 'aggvar')
    >>> H = hurst_estimate(np.random.randn(10000), 'peng', isplot=True)
    >>> H = hurst_estimate(np.random.randn(10000), 'absval', opt=1)
    """
    if method not in _METHODS:
        raise ValueError(
            f"Unknown method '{method}'. Choose from: {list(_METHODS.keys())}")

    func = _METHODS[method]

    if opt is not None:
        return func(sequence, isplot, opt)
    else:
        return func(sequence, isplot)
