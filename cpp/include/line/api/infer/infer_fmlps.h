/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_INFER_INFER_FMLPS_H
#define LINE_API_INFER_INFER_FMLPS_H

/**
 * Fluid response-time likelihood, and the FMLPS demand estimator built on it.
 *
 * Port of matlab/src/api/infer/infer_fluid_ps_rt_likelihood.m and
 * infer_fmlps.m. Both are MATLAB-ONLY: no JAR, no native-Python twin.
 *
 * Reference: Casale et al., "Fluid Analysis of Queueing in Processor Sharing
 * Systems".
 *
 * WHAT THE LIKELIHOOD IS. MLPS builds an exact CTMC per observation and reads a
 * phase-type density; that is exact and costs a state space. FMLPS replaces it
 * with the FLUID limit: mark one unit of fluid at the reference station in the
 * tagged class, let the deterministic drift carry it, and the passage-time
 * density at the observed response time is `-d/dt (marked mass) / mass0`. The
 * marked mass is nonincreasing, so that derivative is the density and the
 * likelihood needs no separate normalization.
 *
 * THE AUGMENTATION IS THE SAME ONE `fluid_passage_time` PERFORMS, and this file
 * deliberately reuses it rather than re-deriving the drift. The reference
 * expands `sn.mu`, `sn.phi`, `sn.proc` and the routing table from K to K+1
 * classes by hand and hands the loose arrays to `solver_fluid_odes`; the C++
 * fluid solver takes a `NetworkStruct`, and the marked-fluid construction it
 * already carries is exactly the tagged class K+1 -- a departure of the marked
 * block UNMARKS the fluid and delivers it wherever the class routes, which is
 * precisely the reference's "absorption at refIdx: the tagged class switches
 * back to the original classes". Re-transcribing the routing expansion would
 * duplicate the one part of the reference where the index arithmetic
 * (`l:Kc:end` against `l:K:end`) is easiest to get wrong.
 *
 * TWO DIFFERENCES FROM THE CTMC ESTIMATOR THAT A CALLER MUST KNOW:
 *
 *  - the fluid density is a LIMIT, so on a small population it is an
 *    approximation where MLPS is exact. It is what buys a model whose state
 *    space MLPS cannot enumerate;
 *  - the reference returns a likelihood of exactly ZERO when the integration
 *    terminates before the observed time (the marked fluid is gone). That is
 *    kept: a zero likelihood is information, and replacing it with a floor
 *    would make an impossible observation look merely unlikely.
 *
 * ARITHMETIC: double, following the ODE.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

#include "line/api/infer/infer_mlps.h"
#include "line/lang/qn/network_struct.h"
#include "line/num/number.h"
#include "line/solvers/fluid/fluid_odes.h"
#include "line/util/error.h"
#include "line/util/matrix.h"
#include "line/solvers/fluid/fluid_stiff.h"
#include "line/util/lsoda.h"
#include "line/util/ode.h"

namespace line {
namespace api {

/** What the fluid likelihood call reports. */
struct FluidRtLikelihood {
    double like = 0.0;          ///< the passage-time density at the observed time
    double marked0 = 0.0;       ///< the marked mass placed at t = 0
    double markedT = 0.0;       ///< what is left of it at the observed time
    std::size_t nstates = 0;    ///< size of the augmented system
};

namespace fmlpsdetail {

/**
 * The augmented drift: the untagged system plus one marked block at (i, c).
 *
 * This is `fluid_passage_time`'s construction. The station's service share is
 * computed on the FOLDED state -- marked mass added back into the block it came
 * from -- because a processor-sharing server does not know which of its jobs is
 * tagged, and computing the share without folding would give the marked job a
 * larger share than it has.
 */
struct TaggedDrift {
    const fluid::FluidOdeSystem* sys;
    std::size_t base, tag0, P, n, nt;
    struct TagEvent {
        std::size_t minus, plus, event_idx;
        double rate_base;
    };
    std::vector<TagEvent> tev;
};

inline TaggedDrift build_tagged(const fluid::FluidOdeSystem& sys, std::size_t i, std::size_t c) {
    TaggedDrift td;
    td.sys = &sys;
    td.P = sys.layout.kic[i][c];
    td.base = sys.layout.qidx[i][c];
    td.n = sys.layout.nstates;
    td.tag0 = td.n;
    td.nt = td.n + td.P;
    for (std::size_t e = 0; e < sys.events.size(); ++e) {
        const fluid::FluidEvent& ev = sys.events[e];
        if (ev.event_idx < td.base || ev.event_idx >= td.base + td.P) continue;
        const std::size_t k = ev.event_idx - td.base;
        TaggedDrift::TagEvent t2;
        t2.event_idx = td.tag0 + k;
        t2.rate_base = ev.rate_base;
        t2.minus = td.tag0 + k;
        if (e < sys.n_departures) {
            // A departure UNMARKS the fluid: it leaves the measured block and
            // arrives untagged wherever the class routes. That is the
            // reference's absorption.
            t2.plus = ev.plus;
        } else {
            t2.plus = td.tag0 + (ev.plus - td.base);  // a phase change stays marked
        }
        td.tev.push_back(t2);
    }
    return td;
}

}  // namespace fmlpsdetail

/**
 * The fluid passage-time density at one observed response time.
 *
 * @param sn       the network
 * @param ist      1-based reference station where the tagged job sits
 * @param cls      1-based tagged class
 * @param levels   (nstates) initial fluid state of the UNAUGMENTED system
 * @param rsampled the observed response time
 * @param marked   the marked mass, the reference's `newFluid`; one job
 */
template <class T>
FluidRtLikelihood infer_fluid_ps_rt_likelihood(const qn::NetworkStruct<T>& sn, std::size_t ist,
                                               std::size_t cls, const std::vector<double>& levels,
                                               double rsampled, double marked = 1.0) {
    const std::size_t M = sn.nstations, K = sn.nclasses;
    if (ist == 0 || ist > M) throw InputError("infer_fluid_ps_rt_likelihood: station out of range");
    if (cls == 0 || cls > K) throw InputError("infer_fluid_ps_rt_likelihood: class out of range");
    if (!(rsampled > 0.0))
        throw InputError("infer_fluid_ps_rt_likelihood: the response time must be positive");
    if (!(marked > 0.0))
        throw InputError("infer_fluid_ps_rt_likelihood: the marked mass must be positive");

    const std::size_t i = ist - 1, c = cls - 1;
    const fluid::FluidOdeSystem sys = fluid::fluid_ode_system(sn);
    const fluid::FluidLayout& L = sys.layout;
    if (levels.size() != L.nstates)
        throw InputError(
            "infer_fluid_ps_rt_likelihood: the initial fluid state has the wrong length for this "
            "model");

    const fmlpsdetail::TaggedDrift td = fmlpsdetail::build_tagged(sys, i, c);
    FluidRtLikelihood out;
    out.nstates = td.nt;
    out.marked0 = marked;
    if (td.P == 0) {
        // The class is not served at this station, so there is no passage to
        // measure and no density to report.
        out.like = 0.0;
        return out;
    }

    // The observed job is MOVED out of its class into the marked block, so the
    // total fluid is unchanged and the state the job found is what the drift
    // starts from.
    std::vector<double> y0(td.nt, 0.0);
    for (std::size_t s = 0; s < td.n; ++s) y0[s] = levels[s];
    double avail = 0.0;
    for (std::size_t k = 0; k < td.P; ++k) avail += y0[td.base + k];
    if (avail + 1e-12 < marked)
        throw InputError(
            "infer_fluid_ps_rt_likelihood: the initial state holds less fluid in the tagged block "
            "than the observation marks, so the observed job is not in the state it arrived to");
    // Remove the marked mass proportionally over the block's phases, which is
    // what leaves the rest of the block undisturbed.
    for (std::size_t k = 0; k < td.P; ++k) y0[td.base + k] -= marked * (y0[td.base + k] / avail);
    y0[td.tag0] = marked;  // all of it enters in phase one

    const std::size_t base = td.base, tag0 = td.tag0, P = td.P, n = td.n, nt = td.nt;
    const std::vector<fmlpsdetail::TaggedDrift::TagEvent> tev = td.tev;
    auto drift = [&sys, &tev, base, tag0, P, n, nt](double, const double* x, double* dx) {
        std::vector<double> xb(x, x + n);
        std::vector<double> g(xb);
        for (std::size_t k = 0; k < P; ++k) g[base + k] += x[tag0 + k];
        std::vector<double> gg(g);
        fluid::fluid_rates_closing(sys, g.data(), gg);
        double blk = 0.0, gblk = 0.0;
        for (std::size_t k = 0; k < P; ++k) {
            blk += g[base + k];
            gblk += gg[base + k];
        }
        const double share = (blk > 0.0) ? gblk / blk : 1.0;

        for (std::size_t s = 0; s < nt; ++s) dx[s] = 0.0;
        for (std::size_t e = 0; e < sys.events.size(); ++e) {
            const fluid::FluidEvent& ev = sys.events[e];
            double drive = gg[ev.event_idx];
            if (ev.event_idx >= base && ev.event_idx < base + P)
                drive = xb[base + (ev.event_idx - base)] * share;
            const double r = ev.rate_base * drive;
            if (r == 0.0) continue;
            dx[ev.minus] -= r;
            dx[ev.plus] += r;
        }
        for (std::size_t e = 0; e < tev.size(); ++e) {
            const double r = tev[e].rate_base * x[tev[e].event_idx] * share;
            if (r == 0.0) continue;
            dx[tev[e].minus] -= r;
            dx[tev[e].plus] += r;
        }
    };

    std::vector<double> grid;
    grid.push_back(0.0);
    grid.push_back(rsampled);
    LsodaOptions lopt;
    lopt.rtol = 1e-5;
    lopt.atol = 1e-8;
    const LsodaSolution sol = fluid::fluid_integrate_grid(drift, y0, grid, lopt);
    if (sol.y.empty()) throw InputError("infer_fluid_ps_rt_likelihood: the ODE returned no state");

    const std::vector<double>& yT = sol.y.back();
    double left = 0.0;
    for (std::size_t k = 0; k < P; ++k) left += std::max(0.0, yT[tag0 + k]);
    out.markedT = left;

    // The density is the rate at which the marked mass is leaving, normalized
    // by the mass placed. Reading the derivative rather than differencing the
    // trajectory is what the reference does and is what stays accurate when the
    // curve is nearly flat.
    std::vector<double> dy(nt, 0.0);
    drift(rsampled, yT.data(), dy.data());
    double ddt = 0.0;
    for (std::size_t k = 0; k < P; ++k) ddt += dy[tag0 + k];
    out.like = -ddt / marked;
    // The marked mass is nonincreasing, so a negative density is round-off at
    // an exhausted block, not a model result.
    if (out.like < 0.0) out.like = 0.0;
    return out;
}

/**
 * FMLPS: the fluid analogue of MLPS.
 *
 * Same objective as `infer_mlps` -- the negative log-likelihood of the observed
 * response times -- with the exact phase-type density replaced by the fluid
 * passage-time density above. The initial fluid state of each observation is
 * the per-class queue length it found, spread over the phases of each block.
 *
 * @param sn      a template network whose station `ist` carries the demands
 *                being estimated; only its service RATES are varied
 * @param ist     1-based PS station
 * @param samples the observations
 */
template <class T>
std::vector<double> infer_fmlps(const qn::NetworkStruct<T>& sn, std::size_t ist,
                                const std::vector<MlpsSample>& samples) {
    const std::size_t M = sn.nstations, K = sn.nclasses;
    if (ist == 0 || ist > M) throw InputError("infer_fmlps: station out of range");
    if (samples.empty()) throw InputError("infer_fmlps: no observations");
    for (std::size_t s = 0; s < samples.size(); ++s) {
        if (samples[s].cls < 1 || samples[s].cls > K)
            throw InputError("infer_fmlps: a sample names a class outside 1..K");
        if (samples[s].ql.size() != K)
            throw InputError("infer_fmlps: a sample's queue length has the wrong width");
        if (!(samples[s].rt > 0.0))
            throw InputError("infer_fmlps: a response time must be positive");
    }

    // The reference's starting point and box, as in MLPS.
    double meanQL = 0.0, rtmax = 0.0;
    for (std::size_t s = 0; s < samples.size(); ++s) {
        for (std::size_t r = 0; r < K; ++r) meanQL += samples[s].ql[r];
        rtmax = std::max(rtmax, samples[s].rt);
    }
    meanQL /= static_cast<double>(samples.size());
    const double nCores = sn.stations[ist - 1].nservers;
    const double Vtilde = std::min(meanQL, std::isinf(nCores) ? meanQL : nCores);
    std::vector<double> x0(K, 1e-3);
    for (std::size_t r = 0; r < K; ++r) {
        double sum = 0.0;
        std::size_t cnt = 0;
        for (std::size_t s = 0; s < samples.size(); ++s)
            if (samples[s].cls == r + 1) {
                sum += samples[s].rt;
                ++cnt;
            }
        if (cnt > 0 && meanQL > 0.0) x0[r] = Vtilde * (sum / static_cast<double>(cnt)) / meanQL;
    }

    const double TOL = 1e-6;
    auto objective = [&](const std::vector<double>& x) {
        for (std::size_t r = 0; r < K; ++r)
            if (!(x[r] > 0.0)) return std::numeric_limits<double>::infinity();

        qn::NetworkStruct<T> mod = sn;
        for (std::size_t r = 0; r < K; ++r)
            mod.set_service(ist, r + 1,
                            lang::Distrib<T>::exp_rate(num_traits<T>::from_double(1.0 / x[r])));
        mod.refresh_rates();

        const fluid::FluidOdeSystem sysx = fluid::fluid_ode_system(mod);
        double f = 0.0;
        for (std::size_t s = 0; s < samples.size(); ++s) {
            // The observed queue length, spread over each block's phases.
            std::vector<double> lev(sysx.layout.nstates, 0.0);
            for (std::size_t r = 0; r < K; ++r) {
                const std::size_t P = sysx.layout.kic[ist - 1][r];
                if (P == 0) continue;
                for (std::size_t k = 0; k < P; ++k)
                    lev[sysx.layout.qidx[ist - 1][r] + k] =
                        samples[s].ql[r] / static_cast<double>(P);
            }
            double like = 0.0;
            try {
                like = infer_fluid_ps_rt_likelihood(mod, ist, samples[s].cls, lev, samples[s].rt)
                           .like;
            } catch (const Error&) {
                // An observation the model cannot host contributes the floor
                // rather than aborting the whole fit.
                like = 0.0;
            }
            f -= std::log(TOL + std::max(0.0, like));
        }
        return f;
    };

    std::vector<Bound<double>> bounds(K);
    for (std::size_t r = 0; r < K; ++r) {
        bounds[r].has_lo = true;
        bounds[r].has_hi = true;
        bounds[r].lo = 0.0;
        bounds[r].hi = rtmax;
    }
    return nelder_mead_box(objective, x0, bounds).x;
}

}  // namespace api
}  // namespace line

#endif  // LINE_API_INFER_INFER_FMLPS_H
