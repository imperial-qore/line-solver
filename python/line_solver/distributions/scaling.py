"""Rate-scaled copy of a distribution, preserving its shape.

Port of matlab/src/lang/processes/dist_scale_rate.m.
"""

import numpy as np

from .base import Distribution


def dist_scale_rate(distrib, factor):
    """Return a new distribution whose rate is FACTOR times the rate of
    DISTRIB, i.e. the time-scaled variable X/FACTOR. The scaling is exact:
    every moment of order n is divided by FACTOR^n, so the mean is divided by
    FACTOR while the SCV, the skewness and the whole shape of the distribution
    are preserved.

    The scaled object is rebuilt from the parameters of the original, rather
    than by rescaling its Markovian representation, so that the parameter list
    stays coherent with the distribution family. Solvers that serialize the
    model (JMT, LDES) read the parameters, and would otherwise export the
    unscaled process.

    This is the perturbation primitive of the finite-difference branch of
    getSensitivityTable: scaling the rate at a station-class by (1+h) is
    exactly the perturbation the derivative d(.)/d(rate) is taken along.

    Args:
        distrib: A Distribution object.
        factor: Positive scaling factor for the rate.

    Returns:
        A new Distribution of the same family with rate*FACTOR.
    """
    from .continuous import (Exp, Det, Immediate, Erlang, HyperExp, Gamma,
                             Lognormal, Pareto, Uniform, Weibull)
    from .markovian import PH, APH, Coxian, Cox2, MAP, MMPP2

    if not isinstance(distrib, Distribution):
        raise ValueError("dist_scale_rate expects a Distribution object.")
    if not np.isscalar(factor) or not np.isfinite(factor) or factor <= 0:
        raise ValueError("The scaling factor must be a positive finite scalar.")
    factor = float(factor)

    # see _kb/03-api-layer.md (dist_scale_rate) for rationale
    cls = type(distrib)
    name = cls.__name__

    if name == 'Exp':
        return Exp(distrib.rate * factor)
    elif name == 'Erlang':
        # Erlang(phase_rate, nphases): the phase count is dimensionless.
        return Erlang(distrib._phase_rate * factor, distrib.phases)
    elif name == 'HyperExp':
        probs = np.asarray(distrib.probs, dtype=float)
        rates = np.asarray(distrib.rates, dtype=float)
        if len(rates) == 2:
            return HyperExp(probs[0], rates[0] * factor, rates[1] * factor)
        return HyperExp(probs, rates * factor)
    elif name in ('Coxian', 'Cox2'):
        # see _kb/03-api-layer.md (dist_scale_rate) for rationale
        rates = 1.0 / np.asarray(distrib.means, dtype=float)
        probs = np.asarray(distrib.probs, dtype=float)
        if name == 'Cox2':
            return Cox2(rates[0] * factor, rates[1] * factor, probs[0])
        return Coxian(rates * factor, probs)
    elif name == 'APH':
        return APH(distrib.alpha, distrib.T * factor)
    elif name == 'PH':
        return PH(distrib.alpha, distrib.T * factor)
    elif name == 'MAP':
        return MAP(distrib.D0 * factor, distrib.D1 * factor)
    elif name == 'MMPP2':
        # Every rate of the modulating chain and of the arrival process is
        # scaled, which time-scales the whole process.
        return MMPP2(distrib.lambda0 * factor, distrib.lambda1 * factor,
                     distrib.sigma0 * factor, distrib.sigma1 * factor)
    elif name == 'Det':
        return Det(distrib.value / factor)
    elif name == 'Uniform':
        return Uniform(distrib.min_val / factor, distrib.max_val / factor)
    elif name == 'Gamma':
        # Gamma(shape, scale): the shape is dimensionless.
        return Gamma(distrib.shape, distrib.scale / factor)
    elif name == 'Pareto':
        # Pareto(shape, scale): the scale is the minimum of the support.
        return Pareto(distrib.alpha, distrib.scale / factor)
    elif name == 'Weibull':
        # Weibull(shape, scale): the shape is dimensionless.
        return Weibull(distrib.shape, distrib.scale / factor)
    elif name == 'Lognormal':
        # X/factor is lognormal with mu - log(factor) and the same sigma.
        return Lognormal(distrib.mu - np.log(factor), distrib.sigma)
    elif name == 'NHPP':
        # see _kb/03-api-layer.md (dist_scale_rate) for rationale
        from .continuous import NHPP
        return NHPP(distrib.getBreakpoints() / factor,
                    distrib.getRates() * factor, distrib.isCyclic())
    elif name == 'Replayer':
        # see _kb/03-api-layer.md (dist_scale_rate) for rationale
        from .discrete import Replayer
        if (getattr(distrib, '_file_path', None) is not None
                and not getattr(distrib, '_file_path_is_temp', False)):
            raise ValueError(
                "Rate scaling is not defined for a file-backed Replayer: the "
                "trace is the parameter. Pass the samples as an array, "
                "Replayer(data), to differentiate it.")
        return Replayer(np.asarray(distrib.trace, dtype=float) / factor,
                        getattr(distrib, '_loop', True))
    elif name == 'Immediate':
        return Immediate()
    else:
        raise ValueError(
            "Rate scaling is not defined for a %s process. Supported: Exp, "
            "Erlang, HyperExp, Coxian, Cox2, APH, PH, MAP, MMPP2, Det, "
            "Uniform, Gamma, Pareto, Weibull, Lognormal, NHPP, Replayer, "
            "Immediate." % name)


# CamelCase alias
distScaleRate = dist_scale_rate
