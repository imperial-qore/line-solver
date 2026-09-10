/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_INFER_INFER_MLPS_H
#define LINE_API_INFER_INFER_MLPS_H

/**
 * Maximum-likelihood service-demand estimation at a processor-sharing queue.
 *
 * Port of matlab/src/api/infer/infer_mlps.m and infer_minps.m. THESE ARE
 * MATLAB-ONLY: there is no JAR or native-Python twin, so MATLAB is not merely
 * the ground truth here, it is the only prior art.
 *
 * WHAT MLPS IS. Each observation is a response time `rt`, the class of the
 * tagged job, and the per-class queue length seen ON ARRIVAL. For a given
 * vector of mean demands, the sojourn of a tagged job that arrives into a known
 * queue state is the absorption time of a small CTMC: build the same PS queue
 * with one EXTRA class carrying the tagged job, mark every transition that is
 * the tagged job departing, and the remaining sub-generator is the phase-type
 * representation of that sojourn. The likelihood of one sample is that
 * phase-type density at `rt`, and the estimate maximizes their product.
 *
 * WHY THE AUGMENTED MODELS ARE BUILT ONCE. The state space, the departure event
 * indices and the absorbing subset depend only on the (tagged class, arrival
 * queue length) PAIR, not on the demands being optimized. The reference caches
 * them per distinct pair and rebuilds only the rates inside the objective; so
 * does this port, because the enumeration is the expensive part and it would
 * otherwise be repeated once per likelihood evaluation.
 *
 * TWO SUBSTITUTIONS FOR MATLAB:
 *
 * 1. `fmincon` with BOX bounds only and no other constraint is
 *    `nelder_mead_box`. The reference passes empty A, b, Aeq, beq and no
 *    nonlinear constraint, so the interior-point machinery is doing nothing an
 *    ordinary box-constrained minimizer does not; the objective is a smooth
 *    log-likelihood in a handful of variables.
 * 2. `sn_set_service_coc` HAS NO C++ COUNTERPART AND NEEDS NONE. It exists only
 *    to write MATLAB's `sn.mu{i}{k}` cell-of-cells without reshaping the cell
 *    array, a storage quirk of that language; here the service law is
 *    `sn.service[i][r]`, an ordinary vector, and setting it is setting it.
 *
 * ARITHMETIC: double. The optimizer and the phase-type density are floating
 * point, and the reference's tolerances are absolute in double.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <string>
#include <vector>

#include "line/api/infer/infer_rps.h"
#include "line/api/mam/map_pdf.h"
#include "line/lang/qn/network_builder.h"
#include "line/lang/qn/network_struct.h"
#include "line/num/number.h"
#include "line/solvers/ctmc/solver_ctmc_getters.h"
#include "line/util/error.h"
#include "line/util/matrix.h"
#include "line/util/neldermead.h"

namespace line {
namespace api {

/** One observation: a response time, the tagged class, and the arrival state. */
struct MlpsSample {
    double rt = 0.0;
    std::size_t cls = 0;             ///< 1-based class of the tagged job
    std::vector<double> ql;          ///< per-class queue length seen on arrival
};

namespace mlpsdetail {

/**
 * Everything about one (tagged class, arrival queue length) pair that does NOT
 * move with the demands: the augmented model, the absorbing subset and the
 * transitions that are the tagged job departing.
 */
struct MlpsPrebuilt {
    std::size_t tagClass = 0;         ///< 1-based, in the ORIGINAL class space
    std::vector<double> N;            ///< augmented population, length R+1
    std::vector<std::size_t> depSync; ///< indices into the sync list
    std::vector<std::size_t> subset;  ///< states with the tagged job at the queue
    Matrix<double> SSqueue;           ///< the queue block of those states
    qn::Network<double> model;        ///< the augmented model, rates still placeholders
    std::size_t queueNode = 0, delayNode = 0;

    explicit MlpsPrebuilt(const qn::Network<double>& m) : model(m) {}
};

/** The row of `SSqueue` equal to `N`, or npos. */
inline std::size_t match_row(const Matrix<double>& S, const std::vector<double>& N) {
    for (std::size_t i = 0; i < S.rows(); ++i) {
        bool eq = true;
        for (std::size_t j = 0; j < S.cols() && j < N.size(); ++j)
            if (std::fabs(S(i, j) - N[j]) > 1e-9) eq = false;
        if (eq) return i;
    }
    return static_cast<std::size_t>(-1);
}

}  // namespace mlpsdetail

/**
 * MLPS demand estimation at a PS queue.
 *
 * @param muZ     (R) think-time rates of the delay station
 * @param nCores  servers at the PS queue
 * @param samples the observations
 * @return (R) estimated mean service demands
 */
inline std::vector<double> infer_mlps(const std::vector<double>& muZ, double nCores,
                                      const std::vector<MlpsSample>& samples) {
    using namespace mlpsdetail;
    const std::size_t R = muZ.size();
    if (R == 0) throw InputError("infer_mlps: at least one class is required");
    if (samples.empty()) throw InputError("infer_mlps: no observations");
    for (std::size_t i = 0; i < samples.size(); ++i) {
        if (samples[i].cls < 1 || samples[i].cls > R)
            throw InputError("infer_mlps: a sample names a class outside 1..R");
        if (samples[i].ql.size() != R)
            throw InputError("infer_mlps: a sample's queue length has the wrong width");
        if (!(samples[i].rt > 0.0))
            throw InputError("infer_mlps: a response time must be positive");
    }

    // ---- the initial point, as the reference forms it -------------------
    double meanQL = 0.0;
    for (std::size_t i = 0; i < samples.size(); ++i)
        for (std::size_t r = 0; r < R; ++r) meanQL += samples[i].ql[r];
    meanQL /= static_cast<double>(samples.size());
    const double Vtilde = std::min(meanQL, nCores);
    std::vector<double> x0(R, 1e-3), lo(R, 0.0), hi(R, 0.0);
    double rtmax = 0.0;
    for (std::size_t i = 0; i < samples.size(); ++i) rtmax = std::max(rtmax, samples[i].rt);
    for (std::size_t r = 0; r < R; ++r) {
        double sum = 0.0;
        std::size_t cnt = 0;
        for (std::size_t i = 0; i < samples.size(); ++i)
            if (samples[i].cls == r + 1) {
                sum += samples[i].rt;
                ++cnt;
            }
        if (cnt > 0 && meanQL > 0.0)
            x0[r] = Vtilde * (sum / static_cast<double>(cnt)) / meanQL;
        hi[r] = rtmax;
    }

    // ---- one augmented model per distinct (tagged class, arrival state) --
    const std::size_t newR = R + 1;
    std::vector<std::string> keys;
    std::vector<MlpsPrebuilt> pre;
    std::map<std::string, std::size_t> keyIndex;
    auto key_of = [](std::size_t tc, const std::vector<double>& q) {
        std::string k = std::to_string(tc);
        for (std::size_t i = 0; i < q.size(); ++i) k += "," + std::to_string(q[i]);
        return k;
    };

    for (std::size_t i = 0; i < samples.size(); ++i) {
        const std::string k = key_of(samples[i].cls, samples[i].ql);
        if (keyIndex.count(k)) continue;

        const std::size_t tc = samples[i].cls;
        // The tagged job is MOVED out of its own class into the extra one, so
        // the total population is unchanged and the queue state observed on
        // arrival is exactly what the augmented chain starts from.
        std::vector<double> N(newR, 0.0);
        for (std::size_t r = 0; r < R; ++r) N[r] = samples[i].ql[r];
        N[tc - 1] -= 1.0;
        N[newR - 1] = 1.0;
        if (N[tc - 1] < 0.0)
            throw InputError(
                "infer_mlps: a sample reports its own class empty on arrival, so the tagged job "
                "cannot be moved into the auxiliary class");

        qn::Network<double> m("mlps_aug");
        const std::size_t d = m.add_delay("Think");
        const std::size_t q = m.add_queue("Queue1", lang::SchedStrategy::PS);
        m.set_number_of_servers(q, nCores);
        std::vector<std::size_t> cls(newR, 0);
        for (std::size_t r = 0; r < newR; ++r) {
            cls[r] = m.add_closed_class("Class" + std::to_string(r + 1), N[r], d);
            // The auxiliary class thinks at the tagged class's own rate.
            const double mz = (r + 1 == newR) ? muZ[tc - 1] : muZ[r];
            m.set_service(d, cls[r], lang::Distrib<double>::exp_rate(mz));
            m.set_service(q, cls[r], lang::Distrib<double>::exp_rate(1.0));  // placeholder
        }
        qn::RoutingMatrix<double> P;
        for (std::size_t r = 0; r < newR; ++r) {
            P.set(cls[r], cls[r], d, q, 1.0);
            P.set(cls[r], cls[r], q, d, 1.0);
        }
        m.link(P);

        MlpsPrebuilt pb(m);
        pb.tagClass = tc;
        pb.N = N;
        pb.queueNode = q;
        pb.delayNode = d;

        const qn::NetworkStruct<double> asn = m.get_struct();
        ctmc::CtmcOptions co;
        const ctmc::CtmcGenerator<double> gen = ctmc::ctmc_get_generator(asn, co);
        const Matrix<double> aggr = ctmc::ctmc_get_state_space_aggr(asn, co);

        // The transitions that ARE the tagged job departing the queue. They are
        // invariant under a rate change, which is why they are cached.
        for (std::size_t e = 0; e < gen.sync.size(); ++e)
            if (gen.sync[e].active.node == q && gen.sync[e].active.cls == newR &&
                gen.sync[e].active.event == lang::EventType::DEP)
                pb.depSync.push_back(e);
        if (pb.depSync.empty())
            throw InputError(
                "infer_mlps: the augmented chain has no departure of the auxiliary class at the "
                "queue, so the sojourn has no absorbing event");

        // The states in which the tagged job is AT the queue: those are the
        // ones the sojourn runs over.
        const std::size_t qst = asn.stations[0].name == "Queue1" ? 1 : 2;
        const std::size_t taggedCol = (qst - 1) * newR + (newR - 1);
        for (std::size_t s = 0; s < aggr.rows(); ++s)
            if (std::fabs(aggr(s, taggedCol) - 1.0) < 1e-9) pb.subset.push_back(s);
        if (pb.subset.empty())
            throw InputError("infer_mlps: no state has the tagged job at the queue");

        pb.SSqueue = Matrix<double>(pb.subset.size(), newR, 0.0);
        for (std::size_t s = 0; s < pb.subset.size(); ++s)
            for (std::size_t r = 0; r < newR; ++r)
                pb.SSqueue(s, r) = aggr(pb.subset[s], (qst - 1) * newR + r);

        keyIndex[k] = pre.size();
        keys.push_back(k);
        pre.push_back(pb);
    }

    // ---- the negative log-likelihood ------------------------------------
    const double TOL = 1e-6;
    auto objective = [&](const std::vector<double>& x) {
        // The demands are means; the model carries rates.
        std::vector<double> rates(R, 0.0);
        for (std::size_t r = 0; r < R; ++r)
            rates[r] = (x[r] > 0.0) ? 1.0 / x[r] : std::numeric_limits<double>::infinity();
        for (std::size_t r = 0; r < R; ++r)
            if (!std::isfinite(rates[r])) return std::numeric_limits<double>::infinity();

        std::vector<Matrix<double>> A(pre.size());
        for (std::size_t p = 0; p < pre.size(); ++p) {
            qn::Network<double> m = pre[p].model;
            for (std::size_t r = 0; r < newR; ++r) {
                const double rate = (r + 1 == newR) ? rates[pre[p].tagClass - 1] : rates[r];
                m.set_service(pre[p].queueNode, r + 1, lang::Distrib<double>::exp_rate(rate));
            }
            const qn::NetworkStruct<double> asn = m.get_struct();
            ctmc::CtmcOptions co;
            const ctmc::CtmcGenerator<double> gen = ctmc::ctmc_get_generator(asn, co);

            // Q minus the tagged departures is the sub-generator of the sojourn.
            Matrix<double> Q = gen.Q;
            for (std::size_t di = 0; di < pre[p].depSync.size(); ++di) {
                const Matrix<double>& F = gen.filt[pre[p].depSync[di]];
                for (std::size_t i = 0; i < Q.rows(); ++i)
                    for (std::size_t j = 0; j < Q.cols(); ++j) Q(i, j) -= F(i, j);
            }
            const std::size_t ns = pre[p].subset.size();
            Matrix<double> S(ns, ns, 0.0);
            for (std::size_t i = 0; i < ns; ++i)
                for (std::size_t j = 0; j < ns; ++j) S(i, j) = Q(pre[p].subset[i], pre[p].subset[j]);
            A[p] = S;
        }

        double f = 0.0;
        for (std::size_t i = 0; i < samples.size(); ++i) {
            const std::size_t p = keyIndex.find(key_of(samples[i].cls, samples[i].ql))->second;
            const Matrix<double>& S = A[p];
            const std::size_t ns = S.rows();
            // The chain starts in the state the sample OBSERVED.
            const std::size_t idx = match_row(pre[p].SSqueue, pre[p].N);
            std::vector<double> pie(ns, 0.0);
            if (idx != static_cast<std::size_t>(-1)) pie[idx] = 1.0;
            // The absorbing MAP: D0 = S, D1 = (-S 1) pie, i.e. every absorption
            // restarts in the observed state, which makes the density of the
            // first passage the phase-type density this needs.
            mam::Map<double> mp;
            mp.D0 = S;
            mp.D1 = Matrix<double>(ns, ns, 0.0);
            for (std::size_t a = 0; a < ns; ++a) {
                double row = 0.0;
                for (std::size_t b = 0; b < ns; ++b) row += S(a, b);
                for (std::size_t b = 0; b < ns; ++b) mp.D1(a, b) = -row * pie[b];
            }
            std::vector<double> at(1, samples[i].rt);
            const double like = mam::map_pdf(mp, at)[0];
            f -= std::log(TOL + std::max(0.0, like));
        }
        return f;
    };

    std::vector<Bound<double>> bounds(R);
    for (std::size_t r = 0; r < R; ++r) {
        bounds[r].has_lo = true;
        bounds[r].has_hi = true;
        bounds[r].lo = lo[r];
        bounds[r].hi = hi[r];
    }
    const NelderMeadResult<double> res = nelder_mead_box(objective, x0, bounds);
    return res.x;
}

/**
 * MINPS: run MLPS and RPS and keep whichever gives the smaller mean demand.
 *
 * The reference's rule verbatim. It is a selection, not a blend: taking the
 * elementwise minimum instead would mix two estimators' class assignments and
 * report a demand vector neither of them produced.
 */
inline std::vector<double> infer_minps(const std::vector<double>& muZ, double nCores,
                                       const std::vector<MlpsSample>& samples) {
    const std::vector<double> mlps = infer_mlps(muZ, nCores, samples);

    std::vector<double> rt;
    std::vector<std::size_t> cls;
    Matrix<double> ql(samples.size(), muZ.size(), 0.0);
    for (std::size_t i = 0; i < samples.size(); ++i) {
        rt.push_back(samples[i].rt);
        // `infer_rps` indexes its classes from ZERO while an MlpsSample carries
        // the reference's 1-based label; passing the label through makes every
        // class below the maximum look empty and the estimator refuses.
        cls.push_back(samples[i].cls - 1);
        for (std::size_t r = 0; r < muZ.size(); ++r) ql(i, r) = samples[i].ql[r];
    }
    const std::vector<double> rps = infer::infer_rps<double>(rt, cls, ql, static_cast<long>(nCores));

    double ma = 0.0, mb = 0.0;
    for (std::size_t r = 0; r < mlps.size(); ++r) ma += mlps[r];
    for (std::size_t r = 0; r < rps.size(); ++r) mb += rps[r];
    if (!mlps.empty()) ma /= static_cast<double>(mlps.size());
    if (!rps.empty()) mb /= static_cast<double>(rps.size());
    return (ma < mb) ? mlps : rps;
}

}  // namespace api
}  // namespace line

#endif  // LINE_API_INFER_INFER_MLPS_H
