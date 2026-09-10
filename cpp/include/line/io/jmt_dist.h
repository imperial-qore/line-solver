/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_IO_JMT_DIST_H
#define LINE_IO_JMT_DIST_H

/**
 * The distribution subtree of a JMT `.jsimg` file, and the scalar formatting
 * every JMT writer shares.
 *
 * ONE EMITTER, SEVEN CALLERS. The reference spells this block out again in
 * `saveServiceStrategy`, `saveArrivalStrategy`, `saveDelayOffStrategy`,
 * `saveSwitchoverStrategy`, `saveTimingStrategies`, `saveRetrialDistributions`
 * and `saveImpatience` -- seven near-copies, which is why those files are the
 * seven largest in `@@JMTIO`. They differ only in the WRAPPER element that
 * carries the pair, never in the pair itself, so this port factors the pair out
 * and leaves each caller its own wrapper. A divergence between two copies of
 * the same emitter is a wrong service time in one place and not the other, and
 * it cannot be seen by reading either copy.
 *
 * WHAT JMT IS GIVEN IS THE FITTED PROCESS, NOT THE DECLARED ONE, for the
 * families whose parameters the reference re-derives from the first two
 * moments: Gamma, Pareto, Weibull, Lognormal and Uniform are all reconstructed
 * from `sn.rates` and `sn.scv` rather than from the constructor arguments. That
 * is the reference's behaviour and is kept, because the moments are what the
 * struct carries after a refresh; a model whose Weibull was fitted by moments
 * therefore round-trips, and one whose shape was set directly is exported as
 * the Justus-approximation Weibull with the same mean and SCV.
 */

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "line/lang/distribution.h"
#include "line/lang/lang_types.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"
#include "line/util/xml.h"

namespace line {
namespace io {

using lang::Distrib;
using lang::ProcessType;

/** `sprintf('%.12f', x)`, the numeric format of every JMT `<value>` body. */
inline std::string jmt_fmt(double x) {
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%.12f", x);
    return std::string(buf);
}

/** `int2str(x)`: round to nearest, print as a decimal integer. */
inline std::string jmt_int(double x) {
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%.0f", x);
    return std::string(buf);
}

/** `num2str(x)`: MATLAB's default five-significant-digit form. */
inline std::string jmt_num(double x) {
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%.5g", x);
    return std::string(buf);
}

/** `num2str(x, 2)`: the two-significant-digit form used for alpha/precision. */
inline std::string jmt_sig2(double x) {
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%.2g", x);
    return std::string(buf);
}

/** `<subParameter classPath="java.lang.Double" name="NAME"><value>V</value>` */
inline xml::Element& jmt_scalar(xml::Element& parent, const char* class_path, const char* name,
                                const std::string& value) {
    xml::Element& p = parent.add_child("subParameter");
    p.set_attr("classPath", class_path);
    p.set_attr("name", name);
    p.add_text_child("value", value);
    return p;
}

/** A `java.lang.Double` scalar subparameter. */
inline xml::Element& jmt_double(xml::Element& parent, const char* name, double v) {
    return jmt_scalar(parent, "java.lang.Double", name, jmt_fmt(v));
}

/**
 * The (D0, D1) pair a JMT MAP or phase-type parameter block is built from,
 * plus the moments the analytic families are reconstructed from.
 *
 * Assembled by `jmt_dist_view` so that a caller holding only a `Distrib` and a
 * caller holding a refreshed `(rate, scv)` row of `sn` reach the same emitter.
 */
template <class T>
struct JmtDistView {
    ProcessType type = ProcessType::DISABLED;
    std::size_t phases = 1;
    double rate = 0.0;  ///< `sn.rates(i,r)`
    double scv = 1.0;   ///< `sn.scv(i,r)`
    std::vector<double> pie;
    Matrix<double> D0, D1;
    std::string trace_file;  ///< Replayer / Trace file name, empty otherwise
};

/**
 * Lower a `Distrib` to the view above.
 *
 * The rate is `1/mean` and NOT `d.params[0]`: `sn.rates` is what the reference
 * reads, and for a fitted family the two differ. A DISABLED or IMMEDIATE
 * distribution is reported by type alone and carries no pair, exactly as
 * `sn.proc` leaves it empty there.
 */
template <class T>
JmtDistView<T> jmt_dist_view(const Distrib<T>& d) {
    JmtDistView<T> v;
    v.type = d.type;
    if (d.disabled || d.type == ProcessType::DISABLED) {
        v.type = ProcessType::DISABLED;
        return v;
    }
    if (d.type == ProcessType::IMMEDIATE) return v;
    const double mean = num_traits<T>::to_double(d.mean);
    v.rate = mean > 0.0 ? 1.0 / mean : 0.0;
    v.scv = num_traits<T>::to_double(d.scv);
    v.phases = d.phases();
    if (d.type == ProcessType::REPLAYER) {
        // The reference reads `sn.nodeparam{ind}{r}.fileName`; this port carries
        // the same path on the distribution (`Distrib::trace_file`). A Replayer
        // built from samples in memory has none, and is refused by the emitter
        // rather than exported as a Replayer reading nothing.
        v.trace_file = d.trace_file;
        return v;
    }
    // Only the Markovian families reach JMT as a matrix pair. For the rest the
    // pair is the refresh's Erlang approximation, which the analytic branches
    // below never look at, so building it would be wasted work and, under exact
    // arithmetic, would refuse on a Gamma that the analytic branch handles.
    switch (d.type) {
        case ProcessType::PH:
        case ProcessType::APH:
        case ProcessType::COXIAN:
        case ProcessType::COX2:
        case ProcessType::ERLANG:
        case ProcessType::HYPEREXP:
        case ProcessType::MAP:
        case ProcessType::MMPP2: {
            const mam::Map<T> m = lang::dist_to_map(d);
            v.D0 = Matrix<double>(m.D0.rows(), m.D0.cols(), 0.0);
            v.D1 = Matrix<double>(m.D1.rows(), m.D1.cols(), 0.0);
            for (std::size_t a = 0; a < m.D0.rows(); ++a)
                for (std::size_t b = 0; b < m.D0.cols(); ++b) {
                    v.D0(a, b) = num_traits<T>::to_double(m.D0(a, b));
                    v.D1(a, b) = num_traits<T>::to_double(m.D1(a, b));
                }
            const std::vector<T> pie = lang::dist_pie(d);
            v.pie.resize(pie.size());
            for (std::size_t k = 0; k < pie.size(); ++k) v.pie[k] = num_traits<T>::to_double(pie[k]);
            v.phases = m.D0.rows();
            break;
        }
        default:
            break;
    }
    return v;
}

/**
 * True when the reference routes the process through JMT's `PhaseTypeDistr`
 * rather than through a named analytic law.
 *
 * COX2 IS ADDED TO THE REFERENCE'S LIST. `saveServiceStrategy` sends PH, APH,
 * COXIAN and a HyperExp of more than two phases down this branch, and lets
 * COX2 fall into the analytic switch -- which has no COX2 case, so MATLAB
 * errors on an undefined `javaClass`. A two-phase Coxian has an exact
 * phase-type representation, and this port's `ProcessType` keeps COX2 distinct
 * from COXIAN, so it is exported through the same exact branch rather than
 * through an error.
 */
inline bool jmt_is_phase_type(ProcessType t, std::size_t phases) {
    return t == ProcessType::PH || t == ProcessType::APH || t == ProcessType::COXIAN ||
           t == ProcessType::COX2 || (phases > 2 && t == ProcessType::HYPEREXP);
}

/** `<subParameter array="true" classPath="java.lang.Object" name="NAME">` */
inline xml::Element& jmt_object_array(xml::Element& parent, const char* name) {
    xml::Element& e = parent.add_child("subParameter");
    e.set_attr("array", "true");
    e.set_attr("classPath", "java.lang.Object");
    e.set_attr("name", name);
    return e;
}

/** One `vector` of `entry` doubles, the JMT encoding of a matrix row. */
inline void jmt_append_row(xml::Element& parent, const Matrix<double>& M, std::size_t row,
                           std::size_t n, bool negate_diagonal) {
    xml::Element& vec = jmt_object_array(parent, "vector");
    for (std::size_t j = 0; j < n; ++j) {
        const double raw = M(row, j);
        const double v = negate_diagonal ? (row == j ? -std::fabs(raw) : std::fabs(raw)) : raw;
        jmt_scalar(vec, "java.lang.Double", "entry", jmt_fmt(v));
    }
}

/**
 * Append the `distribution` and `distrPar` pair for one process.
 *
 * @param wrapper the element the pair is appended to -- the caller's
 *                ServiceTimeStrategy, or a Timing / Impatience wrapper
 * @param v       the lowered process
 * @param who     the caller, for the refusal message of an unexportable family
 *
 * A family JMT has no law for is REFUSED BY NAME. Falling through to the
 * nearest available law would export a model that simulates cleanly and answers
 * a different question, and the featset gate is what is supposed to catch this:
 * a refusal here means the featset and this emitter have drifted apart.
 */
template <class T>
void jmt_append_distribution(xml::Element& wrapper, const JmtDistView<T>& v, const char* who) {
    if (jmt_is_phase_type(v.type, v.phases)) {
        xml::Element& dist = wrapper.add_child("subParameter");
        dist.set_attr("classPath", "jmt.engine.random.PhaseTypeDistr");
        dist.set_attr("name", "Phase-Type");
        xml::Element& par = wrapper.add_child("subParameter");
        par.set_attr("classPath", "jmt.engine.random.PhaseTypePar");
        par.set_attr("name", "distrPar");

        xml::Element& alpha = jmt_object_array(par, "alpha");
        xml::Element& alphavec = jmt_object_array(alpha, "vector");
        for (std::size_t k = 0; k < v.phases; ++k) {
            const double a = k < v.pie.size() ? std::fabs(v.pie[k]) : 0.0;
            jmt_scalar(alphavec, "java.lang.Double", "entry", jmt_fmt(a));
        }
        xml::Element& Tblk = jmt_object_array(par, "T");
        for (std::size_t k = 0; k < v.phases; ++k) jmt_append_row(Tblk, v.D0, k, v.phases, true);
        return;
    }

    if (v.type == ProcessType::MAP || v.type == ProcessType::MMPP2) {
        // MMPP2 reaches BOTH this branch and the analytic MMPP2Distr branch in
        // the reference: the phase-type test above excludes it, and the MAP
        // test here catches it first, so the `MMPP2Par` case of the analytic
        // switch is unreachable. Kept unreachable here too rather than
        // "corrected": the MAP form is exact and the two would otherwise
        // disagree on which one a Burst process is exported as.
        xml::Element& dist = wrapper.add_child("subParameter");
        dist.set_attr("classPath", "jmt.engine.random.MAPDistr");
        dist.set_attr("name", "Burst (MAP)");
        xml::Element& par = wrapper.add_child("subParameter");
        par.set_attr("classPath", "jmt.engine.random.MAPPar");
        par.set_attr("name", "distrPar");
        xml::Element& d0 = jmt_object_array(par, "D0");
        for (std::size_t k = 0; k < v.phases; ++k) jmt_append_row(d0, v.D0, k, v.phases, false);
        xml::Element& d1 = jmt_object_array(par, "D1");
        for (std::size_t k = 0; k < v.phases; ++k) jmt_append_row(d1, v.D1, k, v.phases, false);
        return;
    }

    const char* java_class = nullptr;
    const char* java_par_class = nullptr;
    const char* display = process_to_text(v.type);
    switch (v.type) {
        case ProcessType::DET:
            java_class = "jmt.engine.random.DeterministicDistr";
            java_par_class = "jmt.engine.random.DeterministicDistrPar";
            break;
        case ProcessType::ERLANG:
            java_class = "jmt.engine.random.Erlang";
            java_par_class = "jmt.engine.random.ErlangPar";
            break;
        case ProcessType::EXP:
            java_class = "jmt.engine.random.Exponential";
            java_par_class = "jmt.engine.random.ExponentialPar";
            display = "Exponential";
            break;
        case ProcessType::GAMMA:
            java_class = "jmt.engine.random.GammaDistr";
            java_par_class = "jmt.engine.random.GammaDistrPar";
            break;
        case ProcessType::HYPEREXP:
            java_class = "jmt.engine.random.HyperExp";
            java_par_class = "jmt.engine.random.HyperExpPar";
            display = "Hyperexponential";
            break;
        case ProcessType::PARETO:
            java_class = "jmt.engine.random.Pareto";
            java_par_class = "jmt.engine.random.ParetoPar";
            break;
        case ProcessType::WEIBULL:
            java_class = "jmt.engine.random.Weibull";
            java_par_class = "jmt.engine.random.WeibullPar";
            break;
        case ProcessType::LOGNORMAL:
            java_class = "jmt.engine.random.Lognormal";
            java_par_class = "jmt.engine.random.LognormalPar";
            break;
        case ProcessType::UNIFORM:
            java_class = "jmt.engine.random.Uniform";
            java_par_class = "jmt.engine.random.UniformPar";
            break;
        case ProcessType::REPLAYER:
            java_class = "jmt.engine.random.Replayer";
            java_par_class = "jmt.engine.random.ReplayerPar";
            display = "Replayer";
            break;
        default:
            throw UnsupportedError(std::string(who) + ": JMT has no distribution for '" +
                                   process_to_text(v.type) + "'");
    }

    xml::Element& dist = wrapper.add_child("subParameter");
    dist.set_attr("classPath", java_class);
    dist.set_attr("name", display);
    xml::Element& par = wrapper.add_child("subParameter");
    par.set_attr("classPath", java_par_class);
    par.set_attr("name", "distrPar");

    switch (v.type) {
        case ProcessType::DET:
            jmt_double(par, "t", v.rate > 0.0 ? 1.0 / v.rate : 0.0);
            break;
        case ProcessType::EXP:
            jmt_double(par, "lambda", v.rate);
            break;
        case ProcessType::HYPEREXP:
            jmt_double(par, "p", v.pie.empty() ? 0.0 : v.pie[0]);
            jmt_double(par, "lambda1", v.D0.rows() > 0 ? -v.D0(0, 0) : 0.0);
            jmt_double(par, "lambda2", v.D0.rows() > 1 ? -v.D0(1, 1) : 0.0);
            break;
        case ProcessType::ERLANG:
            jmt_double(par, "alpha", v.rate * static_cast<double>(v.phases));
            jmt_scalar(par, "java.lang.Long", "r", jmt_int(static_cast<double>(v.phases)));
            break;
        case ProcessType::GAMMA:
            jmt_double(par, "alpha", 1.0 / v.scv);
            jmt_double(par, "beta", v.scv / v.rate);
            break;
        case ProcessType::PARETO: {
            const double shape = std::sqrt(1.0 + 1.0 / v.scv) + 1.0;
            const double scale = (1.0 / v.rate) * (shape - 1.0) / shape;
            jmt_double(par, "alpha", shape);
            jmt_double(par, "k", scale);
            break;
        }
        case ProcessType::WEIBULL: {
            // Justus (1976) approximation of the shape from the SCV, as the
            // reference uses; `alpha` is JMT's scale and `r` its shape.
            const double c = std::sqrt(v.scv);
            const double rval = std::pow(c, -1.086);
            const double alpha = (1.0 / v.rate) / std::tgamma(1.0 + 1.0 / rval);
            jmt_double(par, "alpha", alpha);
            jmt_double(par, "r", rval);
            break;
        }
        case ProcessType::LOGNORMAL: {
            const double c = std::sqrt(v.scv);
            const double mu = std::log((1.0 / v.rate) / std::sqrt(c * c + 1.0));
            const double sigma = std::sqrt(std::log(c * c + 1.0));
            jmt_double(par, "mu", mu);
            jmt_double(par, "sigma", sigma);
            break;
        }
        case ProcessType::UNIFORM: {
            const double maxVal =
                (std::sqrt(12.0 * v.scv / (v.rate * v.rate)) + 2.0 / v.rate) / 2.0;
            const double minVal = 2.0 / v.rate - maxVal;
            jmt_double(par, "min", minVal);
            jmt_double(par, "max", maxVal);
            break;
        }
        case ProcessType::REPLAYER:
            if (v.trace_file.empty())
                throw UnsupportedError(std::string(who) +
                                       ": a Replayer must name the trace file JMT is to read; "
                                       "this port carries the samples on the distribution and the "
                                       "caller has not staged them beside the model");
            jmt_scalar(par, "java.lang.String", "fileName", v.trace_file);
            break;
        default:
            break;
    }
}

}  // namespace io
}  // namespace line

#endif  // LINE_IO_JMT_DIST_H
