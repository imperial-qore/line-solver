"""
Laplace Transform Inversion (LTI) utilities.

Provides multiple methods for inverting Laplace transforms numerically:
- Talbot contour integration
- Gaverstehfest algorithm
- Euler method
- Abate-Whitt formula
- Laguerre expansion
- Custom Romberg integration

These methods are essential for moment fitting and distribution analysis.
"""

import numpy as np
from line_solver import native_to_array


def lib_lti_talbot(laplace_func, t, M=32):
    """
    Invert Laplace transform using Talbot's contour integration method.

    This method provides good numerical stability and accuracy.

    Args:
        laplace_func (callable): Function s -> F(s) (Laplace transform)
        t (float or array-like): Time points to evaluate
        M (int): Number of integration points (default 32)

    Returns:
        float or ndarray: Values of f(t)
    """
    t = np.atleast_1d(t)

    try:
        results = []
        for t_val in t:
            result = TalbotLTI.invert(laplace_func, float(t_val), M)
            results.append(float(result))
        return results[0] if len(results) == 1 else np.array(results)
    except Exception as e:
        raise RuntimeError(f"Talbot LTI failed: {str(e)}")


def lib_lti_gaverstehfest(laplace_func, t):
    """
    Invert Laplace transform using Gaverstehfest algorithm.

    Good for general-purpose inversion, especially for early times.

    Args:
        laplace_func (callable): Function s -> F(s)
        t (float or array-like): Time points

    Returns:
        float or ndarray: Values of f(t)
    """
    t = np.atleast_1d(t)

    try:
        results = []
        for t_val in t:
            result = GS.invert(laplace_func, float(t_val))
            results.append(float(result))
        return results[0] if len(results) == 1 else np.array(results)
    except Exception as e:
        raise RuntimeError(f"Gaverstehfest LTI failed: {str(e)}")


def lib_lti_euler(laplace_func, t, M=16):
    """
    Invert Laplace transform using Euler method.

    Stable method with good convergence for many functions.

    Args:
        laplace_func (callable): Function s -> F(s)
        t (float or array-like): Time points
        M (int): Number of terms (default 16)

    Returns:
        float or ndarray: Values of f(t)
    """
    t = np.atleast_1d(t)

    try:
        results = []
        for t_val in t:
            result = Euler.invert(laplace_func, float(t_val), M)
            results.append(float(result))
        return results[0] if len(results) == 1 else np.array(results)
    except Exception as e:
        raise RuntimeError(f"Euler LTI failed: {str(e)}")


def lib_lti_abatewhitt(laplace_func, t, M=16):
    """
    Invert Laplace transform using Abate-Whitt formula.

    Efficient method with controlled accuracy.

    Args:
        laplace_func (callable): Function s -> F(s)
        t (float or array-like): Time points
        M (int): Number of terms (default 16)

    Returns:
        float or ndarray: Values of f(t)
    """
    t = np.atleast_1d(t)

    try:
        results = []
        for t_val in t:
            result = AW.invert(laplace_func, float(t_val), M)
            results.append(float(result))
        return results[0] if len(results) == 1 else np.array(results)
    except Exception as e:
        raise RuntimeError(f"Abate-Whitt LTI failed: {str(e)}")


def lib_lti_laguerre(laplace_func, t, M=32):
    """
    Invert Laplace transform using Laguerre expansion.

    Good for probability distributions and densities.

    Args:
        laplace_func (callable): Function s -> F(s)
        t (float or array-like): Time points
        M (int): Number of Laguerre terms (default 32)

    Returns:
        float or ndarray: Values of f(t)
    """
    t = np.atleast_1d(t)

    try:
        results = []
        for t_val in t:
            result = Laguerre.invert(laplace_func, float(t_val), M)
            results.append(float(result))
        return results[0] if len(results) == 1 else np.array(results)
    except Exception as e:
        raise RuntimeError(f"Laguerre LTI failed: {str(e)}")


def lib_lti_custom_romberg(laplace_func, t, tol=1e-6):
    """
    Invert Laplace transform using custom Romberg integration.

    Adaptive method with user-specified tolerance.

    Args:
        laplace_func (callable): Function s -> F(s)
        t (float or array-like): Time points
        tol (float): Tolerance for convergence

    Returns:
        float or ndarray: Values of f(t)
    """
    t = np.atleast_1d(t)

    try:
        results = []
        for t_val in t:
            result = CR.invert(laplace_func, float(t_val), float(tol))
            results.append(float(result))
        return results[0] if len(results) == 1 else np.array(results)
    except Exception as e:
        raise RuntimeError(f"Custom Romberg LTI failed: {str(e)}")


def lib_lti_cme(laplace_func, t, M=16):
    """
    Invert Laplace transform using Crump-type method with exponential.

    Good for wide range of function types.

    Args:
        laplace_func (callable): Function s -> F(s)
        t (float or array-like): Time points
        M (int): Number of terms (default 16)

    Returns:
        float or ndarray: Values of f(t)
    """
    t = np.atleast_1d(t)

    try:
        results = []
        for t_val in t:
            result = CME.invert(laplace_func, float(t_val), M)
            results.append(float(result))
        return results[0] if len(results) == 1 else np.array(results)
    except Exception as e:
        raise RuntimeError(f"CME LTI failed: {str(e)}")


def lib_lti_compare_methods(laplace_func, t, methods=None):
    """
    Compare multiple LTI methods at a point.

    Args:
        laplace_func (callable): Laplace transform function
        t (float): Time point to evaluate
        methods (list): Which methods to use. Default: all available

    Returns:
        dict: Results from each method
    """
    if methods is None:
        methods = ['talbot', 'gaverstehfest', 'euler', 'abatewhitt', 'laguerre']

    results = {}

    for method in methods:
        try:
            if method == 'talbot':
                results['talbot'] = lib_lti_talbot(laplace_func, t)
            elif method == 'gaverstehfest':
                results['gaverstehfest'] = lib_lti_gaverstehfest(laplace_func, t)
            elif method == 'euler':
                results['euler'] = lib_lti_euler(laplace_func, t)
            elif method == 'abatewhitt':
                results['abatewhitt'] = lib_lti_abatewhitt(laplace_func, t)
            elif method == 'laguerre':
                results['laguerre'] = lib_lti_laguerre(laplace_func, t)
            elif method == 'romberg':
                results['romberg'] = lib_lti_custom_romberg(laplace_func, t)
            elif method == 'cme':
                results['cme'] = lib_lti_cme(laplace_func, t)
        except Exception as e:
            results[method] = f'Error: {str(e)}'

    return results
