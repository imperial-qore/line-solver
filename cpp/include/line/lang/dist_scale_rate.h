/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_LANG_DIST_SCALE_RATE_H
#define LINE_LANG_DIST_SCALE_RATE_H

/**
 * Rate-scaled copy of a distribution, preserving its shape.
 *
 * Port of `matlab/src/api/dist/dist_scale_rate.m`. The result is the law of
 * X / factor: every moment of order n is divided by factor^n, so the mean moves
 * and the SCV, the skewness and the whole shape do not.
 *
 * REBUILT FROM THE PARAMETERS, not by rescaling (D0, D1). Both routes give the
 * same law, but only the first leaves the parameter list coherent with the
 * family, and the parameters are what a serializing consumer writes out. The
 * reference gives the same reason.
 *
 * THIS IS THE PERTURBATION PRIMITIVE of the finite-difference branch of
 * getSensitivityTable: scaling the rate at a (station, class) by (1 + h) is
 * exactly the perturbation d(.)/d(rate) is taken along, which is why the
 * scaling has to be exact in every moment and not merely in the mean.
 *
 * REFUSED BY NAME, as the reference's `otherwise` arm refuses them: MMAP, BMAP,
 * MAPt, PHt and Prior. A marked or time-inhomogeneous process carries more than
 * one time scale (the per-mark blocks, the breakpoint schedule), and scaling
 * only the one the caller had in mind would silently change the others.
 */

#include <cmath>
#include <cstddef>
#include <vector>

#include "line/lang/lang_types.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace lang {

/** The law of X / factor, in the same family as `d`. */
template <class T>
Distrib<T> dist_scale_rate(const Distrib<T>& d, const T& factor) {
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    if (!(factor > zero))
        throw InputError("dist_scale_rate: the scaling factor must be positive");
    if (d.is_prior())
        throw UnsupportedError(
            "dist_scale_rate: a Prior is a set of alternative models, not one law with a time "
            "scale; perturb the alternatives instead");
    if (d.has_schedule())
        throw UnsupportedError(
            "dist_scale_rate: a MAPt / PHt / NHPP schedule carries a second time scale in its "
            "breakpoints, which scaling the rates alone would leave inconsistent");
    if (d.disabled) return d;

    switch (d.type) {
        case ProcessType::IMMEDIATE:
            return Distrib<T>::immediate();
        case ProcessType::EXP:
            return Distrib<T>::exp_rate(T(d.params[0] * factor));
        case ProcessType::ERLANG:
            return Distrib<T>::erlang(
                T(d.params[0] * factor),
                static_cast<std::size_t>(num_traits<T>::to_double(d.params[1]) + 0.5));
        case ProcessType::HYPEREXP:
            return Distrib<T>::hyperexp(d.params[0], T(d.params[1] * factor),
                                        T(d.params[2] * factor));
        case ProcessType::COXIAN:
        case ProcessType::COX2: {
            // params = [mu(1..n), phi(1..n)]. The completion probabilities are
            // dimensionless, so only the phase rates carry the time scale.
            const std::size_t n = d.params.size() / 2;
            std::vector<T> mu(n), phi(n);
            for (std::size_t i = 0; i < n; ++i) {
                mu[i] = T(d.params[i] * factor);
                phi[i] = d.params[n + i];
            }
            return Distrib<T>::coxian(mu, phi);
        }
        case ProcessType::APH:
        case ProcessType::PH: {
            // params carries alpha; the subgenerator is D0.
            std::vector<T> alpha(d.params.begin(), d.params.end());
            Matrix<T> A(d.D0.rows(), d.D0.cols(), zero);
            for (std::size_t i = 0; i < d.D0.rows(); ++i)
                for (std::size_t j = 0; j < d.D0.cols(); ++j) A(i, j) = T(d.D0(i, j) * factor);
            return Distrib<T>::phase_type(alpha, A, d.type == ProcessType::APH);
        }
        case ProcessType::MAP:
        case ProcessType::MMPP2: {
            // Every rate of the modulating chain and of the arrival process is
            // scaled, which time-scales the whole process.
            Matrix<T> D0(d.D0.rows(), d.D0.cols(), zero), D1(d.D1.rows(), d.D1.cols(), zero);
            for (std::size_t i = 0; i < d.D0.rows(); ++i)
                for (std::size_t j = 0; j < d.D0.cols(); ++j) D0(i, j) = T(d.D0(i, j) * factor);
            for (std::size_t i = 0; i < d.D1.rows(); ++i)
                for (std::size_t j = 0; j < d.D1.cols(); ++j) D1(i, j) = T(d.D1(i, j) * factor);
            Distrib<T> out = Distrib<T>::map_dist(D0, D1, d.type);
            for (const T& p : d.params) out.params.push_back(T(p * factor));
            out.mean = T(d.mean / factor);
            out.scv = d.scv;
            return out;
        }
        case ProcessType::DET:
            return Distrib<T>::det(T(d.params[0] / factor));
        case ProcessType::UNIFORM:
            return Distrib<T>::uniform(T(d.params[0] / factor), T(d.params[1] / factor));
        case ProcessType::GAMMA:
            // Gamma(shape, scale): the shape is dimensionless.
            return Distrib<T>::gamma_dist(d.params[0], T(d.params[1] / factor));
        case ProcessType::PARETO:
            // Pareto(shape, scale): the scale is the minimum of the support.
            return Distrib<T>::pareto(d.params[0], T(d.params[1] / factor));
        case ProcessType::WEIBULL:
            return Distrib<T>::weibull(T(d.params[0] / factor), d.params[1]);
        case ProcessType::LOGNORMAL: {
            if constexpr (!num_traits<T>::has_transcendental) {
                throw UnsupportedError(
                    "dist_scale_rate: scaling a Lognormal shifts its log-mean by log(factor), "
                    "which exact arithmetic has no representation for");
            } else {
                const double lf = std::log(num_traits<T>::to_double(factor));
                return Distrib<T>::lognormal(T(d.params[0] - num_traits<T>::from_double(lf)),
                                             d.params[1]);
            }
        }
        case ProcessType::REPLAYER: {
            // A trace is scaled sample by sample; the samples ARE the parameter.
            std::vector<T> s(d.trace.size(), zero);
            for (std::size_t i = 0; i < d.trace.size(); ++i) s[i] = T(d.trace[i] / factor);
            return Distrib<T>::replayer(s);
        }
        default:
            break;
    }
    (void)one;
    throw UnsupportedError(
        "dist_scale_rate: rate scaling is not defined for this process family; supported are "
        "Exp, Erlang, HyperExp, Coxian, Cox2, APH, PH, MAP, MMPP2, Det, Uniform, Gamma, Pareto, "
        "Weibull, Lognormal, Replayer and Immediate");
}

}  // namespace lang
}  // namespace line

#endif  // LINE_LANG_DIST_SCALE_RATE_H
