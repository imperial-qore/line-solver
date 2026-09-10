/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_FES_FES_AGGREGATE_H
#define LINE_API_FES_FES_AGGREGATE_H

/**
 * Flow-equivalent-server aggregation: replace a station subset by one station.
 *
 * Templated port of matlab/src/io/@ModelAdapter/aggregateFES.m (the JAR twin is
 * the deprecated `jline.api.fes.FESAggregator`, reached through
 * `ModelAdapter.aggregateFES`; MATLAB is the reference here).
 *
 * This is Chandy-Herzog-Woo's Norton theorem in its state-dependent form. The
 * subset is short-circuited and solved in isolation at every population on the
 * lattice; the resulting per-class throughputs X_r(n) become the service rates
 * of a single limited-class-dependent station, and the complement is rewired to
 * route through it. For a product-form closed network the replacement is EXACT
 * at every population, which is what makes the transform worth doing rather
 * than an approximation to be tuned.
 *
 * THREE THINGS THE ARITHMETIC DEPENDS ON, each easy to get wrong:
 *
 * 1. The routing seen by the subset is the STOCHASTIC COMPLEMENT of the full
 *    chain on the subset's rows, not the raw submatrix. A job that leaves the
 *    subset and comes back through the complement must re-enter with the right
 *    probability, and `dtmc_stochcomp` is what folds those excursions back in.
 *
 * 2. The isolated throughput is a rate INSIDE the subnetwork, whose visit
 *    ratios are normalized to the subset. Turning it into a rate the outside
 *    sees needs the per-class ESCAPE factor -- the visit-weighted probability
 *    of leaving the subset per subset visit. Without it the FES completes jobs
 *    at the subnetwork's internal circulation rate, which is too fast by
 *    exactly the number of internal hops per escape.
 *
 * 3. The class dependence carries `beta_r(n) = X_r(n) |n| / n_r`, not X_r(n).
 *    The |n|/n_r cancels the processor-sharing split the convolution applies
 *    (Sauer 1983, eq. 40), so the aggregate really completes class r at X_r(n).
 *    `fes_beta_handle` owns that factor; this file must not apply it twice.
 *
 * ONE DELIBERATE DEPARTURE FROM MATLAB, and it is a correctness fix rather than
 * a convention: the reference rebuilds each complement station's service law
 * from `sn.proc{i}{k}` as `APH(ones(1,n)/n, T)` whenever the process has more
 * than one phase, i.e. it DISCARDS the true initial phase vector and substitutes
 * a uniform one. That silently changes the distribution of every non-exponential
 * complement station (an Erlang(k) becomes a mixture starting in a random
 * phase, with a different mean and a much larger SCV). The port copies the
 * station's own `Distrib` verbatim instead, so a complement station keeps
 * exactly the law it had. On exponential service -- the case the reference's own
 * tests exercise -- the two agree, since a one-phase process is rebuilt as
 * Exp(rate) either way.
 *
 * ARITHMETIC: transcendental, inherited from the convolution behind
 * fes_compute_throughputs.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

#include "line/api/fes/fes_beta_handle.h"
#include "line/api/fes/fes_build_isolated.h"
#include "line/api/fes/fes_compute_throughputs.h"
#include "line/api/fes/fes_validate.h"
#include "line/api/mc/dtmc_stochcomp.h"
#include "line/api/pfqn/cd_peak_scaling.h"
#include "line/lang/qn/network_builder.h"
#include "line/lang/qn/network_struct.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace fes {

using lang::CdScaling;
using lang::SchedStrategy;

/** `options` of the reference; the solver field is implied by the convolution. */
struct FesOptions {
    std::vector<int> cutoffs;  ///< per-class population bound; empty means sn.njobs
    bool verbose = false;
};

/** Everything needed to map an FES result back onto the original model. */
template <class T>
struct FesDeaggInfo {
    std::vector<std::size_t> subsetIndices;      ///< 1-based, as given
    std::vector<std::size_t> complementIndices;  ///< 1-based
    std::vector<std::vector<T>> throughputTable; ///< per class, linearized on the lattice
    std::vector<int> cutoffs;
    Matrix<T> stochCompSubset, stochCompComplement;
    Matrix<T> isolatedDemands;  ///< (M_sub x K)
    Matrix<T> isolatedVisits;   ///< (M_sub x K)
    std::vector<int> isolatedServers;
    std::vector<bool> isolatedIsDelay;
    std::vector<T> escape;         ///< per-class Norton escape factor
    std::size_t fesNode = 0;       ///< 1-based node index of the FES in the new model
};

/** What fes_aggregate returns. */
template <class T>
struct FesAggregateResult {
    explicit FesAggregateResult(const qn::Network<T>& m) : model(m) {}
    qn::Network<T> model;
    std::size_t fesNode = 0;  ///< 1-based node index of the FES station
    FesDeaggInfo<T> deagg;
};

/**
 * @param sn            the original closed product-form network
 * @param subsetIndices 1-BASED station indices to aggregate, the reference's base
 * @param options       cutoffs and verbosity
 */
template <class T>
FesAggregateResult<T> fes_aggregate(const qn::NetworkStruct<T>& sn,
                                    const std::vector<std::size_t>& subsetIndices,
                                    const FesOptions& options = FesOptions()) {
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    const double fineTol = lang::GlobalConstants::FineTol;
    const std::size_t M = sn.nstations, K = sn.nclasses;

    const FesValidateResult ok = fes_validate(sn, subsetIndices);
    if (!ok.isValid) throw InputError("aggregateFES: " + ok.errorMsg);

    std::vector<std::size_t> complementIndices;
    for (std::size_t i = 1; i <= M; ++i) {
        bool inSubset = false;
        for (std::size_t a = 0; a < subsetIndices.size(); ++a)
            if (subsetIndices[a] == i) inSubset = true;
        if (!inSubset) complementIndices.push_back(i);
    }

    const std::vector<double> N = sn.njobs();
    std::vector<int> cutoffs = options.cutoffs;
    if (cutoffs.empty())
        for (std::size_t r = 0; r < K; ++r) cutoffs.push_back(static_cast<int>(N[r]));
    if (cutoffs.size() != K)
        throw InputError("aggregateFES: cutoffs must carry one entry per class");

    // ---- the two stochastic complements ---------------------------------
    // `rt` is indexed (stateful-1)*K + class, so a station contributes the K
    // consecutive rows of its own stateful index.
    std::vector<std::size_t> subsetRt, complementRt;
    for (std::size_t a = 0; a < subsetIndices.size(); ++a) {
        const std::size_t isf = sn.stateful_of_station(subsetIndices[a]);
        for (std::size_t r = 0; r < K; ++r) subsetRt.push_back((isf - 1) * K + r);
    }
    for (std::size_t a = 0; a < complementIndices.size(); ++a) {
        const std::size_t isf = sn.stateful_of_station(complementIndices[a]);
        for (std::size_t r = 0; r < K; ++r) complementRt.push_back((isf - 1) * K + r);
    }
    const Matrix<T> stochCompSubset = mc::dtmc_stochcomp(sn.rt, subsetRt);
    const Matrix<T> stochCompComplement = mc::dtmc_stochcomp(sn.rt, complementRt);

    // ---- the isolated subnetwork ----------------------------------------
    const std::size_t nSub = subsetIndices.size();
    Matrix<T> subRates(nSub, K, zero);
    std::vector<int> mi(nSub, 1);
    std::vector<bool> isDelay(nSub, false);
    for (std::size_t a = 0; a < nSub; ++a) {
        const qn::Station<T>& st = sn.stations[subsetIndices[a] - 1];
        for (std::size_t r = 0; r < K; ++r) subRates(a, r) = sn.rates(subsetIndices[a] - 1, r);
        isDelay[a] = (st.nodetype == qn::NodeType::Delay || st.sched == SchedStrategy::INF);
        mi[a] = std::isinf(st.nservers) ? 1 : static_cast<int>(st.nservers);
    }
    const FesIsolated<T> iso = fes_build_isolated(subRates, stochCompSubset);
    std::vector<std::vector<T>> scalingTable =
        fes_compute_throughputs(iso.L, mi, isDelay, cutoffs);

    // ---- the Norton escape factor ---------------------------------------
    // Per subset visit, the visit-weighted probability of leaving the subset.
    // The isolated throughput counts every internal hop, so without this the
    // FES would complete jobs at the internal circulation rate.
    std::vector<T> escape(K, zero);
    for (std::size_t r = 0; r < K; ++r) {
        for (std::size_t a = 0; a < nSub; ++a) {
            const std::size_t j = subsetIndices[a];
            const T Vjr = T(iso.L(a, r) * sn.rates(j - 1, r));
            const std::size_t isf_j = sn.stateful_of_station(j);
            T pexit = zero;
            for (std::size_t b = 0; b < complementIndices.size(); ++b) {
                const std::size_t isf_i = sn.stateful_of_station(complementIndices[b]);
                pexit += sn.rt((isf_j - 1) * K + r, (isf_i - 1) * K + r);
            }
            if (std::isfinite(num_traits<T>::to_double(Vjr))) escape[r] += T(Vjr * pexit);
        }
        if (num_traits<T>::to_double(escape[r]) > fineTol)
            for (std::size_t idx = 0; idx < scalingTable[r].size(); ++idx)
                scalingTable[r][idx] *= escape[r];
    }

    // ---- the FES model ---------------------------------------------------
    qn::Network<T> fesModel(sn.name.empty() ? std::string("FES") : sn.name + "_FES");

    // The complement stations, in their original order, then the FES.
    std::vector<std::size_t> complementNode(M + 1, 0);  // 1-based station -> new node
    for (std::size_t b = 0; b < complementIndices.size(); ++b) {
        const std::size_t i = complementIndices[b];
        const qn::Station<T>& st = sn.stations[i - 1];
        std::size_t nd = 0;
        if (st.nodetype == qn::NodeType::Delay) {
            nd = fesModel.add_delay(st.name);
        } else if (st.nodetype == qn::NodeType::Queue) {
            nd = fesModel.add_queue(st.name, st.sched);
            if (!std::isinf(st.nservers)) fesModel.set_number_of_servers(nd, st.nservers);
            if (st.cap > 0.0 && std::isfinite(st.cap))
                fesModel.set_capacity(nd, static_cast<int>(st.cap));
        } else {
            throw InputError("aggregateFES: unsupported station type in the complement");
        }
        complementNode[i] = nd;
    }
    const std::size_t fesNode = fesModel.add_queue("FES", SchedStrategy::PS);
    fesModel.set_number_of_servers(fesNode, 1.0);

    // The reference anchors the classes at the first complement station, and
    // at the FES only when the complement is empty (which fes_validate rules
    // out, since the subset must be proper).
    const std::size_t refNode =
        complementIndices.empty() ? fesNode : complementNode[complementIndices[0]];
    std::vector<std::size_t> newClass(K, 0);
    for (std::size_t r = 0; r < K; ++r)
        newClass[r] = fesModel.add_closed_class(sn.classes[r].name, N[r], refNode);

    // The complement keeps its own service laws, verbatim; see the header note
    // on why the reference's APH rebuild is not reproduced.
    for (std::size_t b = 0; b < complementIndices.size(); ++b) {
        const std::size_t i = complementIndices[b];
        for (std::size_t r = 0; r < K; ++r)
            fesModel.set_service(complementNode[i], newClass[r], sn.service[i - 1][r]);
    }

    // The FES serves at rate one and is scaled entirely by the class dependence.
    for (std::size_t r = 0; r < K; ++r)
        fesModel.set_service(fesNode, newClass[r], lang::Distrib<T>::exp_rate(one));

    // A zero entry would stall the recurrence, so the reference floors the
    // table rather than letting the handle return zero.
    const T floorVal = num_traits<T>::from_double(fineTol);
    for (std::size_t r = 0; r < K; ++r)
        for (std::size_t idx = 0; idx < scalingTable[r].size(); ++idx)
            if (num_traits<T>::to_double(scalingTable[r][idx]) < fineTol)
                scalingTable[r][idx] = floorVal;

    const FesBetaFun<T> beta = fes_beta_handle(scalingTable, cutoffs);
    const T peak = pfqn::cd_peak_scaling<T>(beta, cutoffs);
    // `sn.cdscaling` is evaluated on a population vector in the WORKING
    // arithmetic while the FES table is indexed by integer counts, so the
    // handle is wrapped rather than re-tabulated. The rounding is the
    // reference's own -- fes_beta_handle rounds before linearizing.
    const CdScaling<T> cd = [beta](const std::vector<T>& n) {
        std::vector<int> ni(n.size(), 0);
        for (std::size_t r = 0; r < n.size(); ++r)
            ni[r] = static_cast<int>(std::lround(num_traits<T>::to_double(n[r])));
        return beta(ni);
    };
    fesModel.set_class_dependence(fesNode, cd, std::vector<T>(1, peak));

    // ---- the routing -----------------------------------------------------
    qn::RoutingMatrix<T> P;
    for (std::size_t r = 0; r < K; ++r) {
        const std::size_t I = complementIndices.size() + 1;  // the complement plus the FES
        std::vector<std::vector<T>> Pk(I + 1, std::vector<T>(I + 1, zero));

        for (std::size_t b = 0; b < complementIndices.size(); ++b) {
            const std::size_t i = complementIndices[b];
            const std::size_t iNode = complementNode[i];
            const std::size_t isf_i = sn.stateful_of_station(i);

            // Complement to complement: the ORIGINAL probability, since those
            // paths do not pass through the aggregated subset.
            for (std::size_t c = 0; c < complementIndices.size(); ++c) {
                const std::size_t j = complementIndices[c];
                const std::size_t isf_j = sn.stateful_of_station(j);
                const T p = sn.rt((isf_i - 1) * K + r, (isf_j - 1) * K + r);
                if (num_traits<T>::to_double(p) > fineTol) Pk[iNode][complementNode[j]] = p;
            }
            // Complement to the subset: every such path now ends at the FES.
            for (std::size_t a = 0; a < nSub; ++a) {
                const std::size_t isf_j = sn.stateful_of_station(subsetIndices[a]);
                const T p = sn.rt((isf_i - 1) * K + r, (isf_j - 1) * K + r);
                if (num_traits<T>::to_double(p) > fineTol) Pk[iNode][fesNode] += p;
            }
        }

        // The FES leaves for the complement with the VISIT-WEIGHTED exit
        // probability of the subset: which subset station a job departs from is
        // no longer represented, so its visit ratio is what stands in.
        for (std::size_t c = 0; c < complementIndices.size(); ++c) {
            const std::size_t j = complementIndices[c];
            const std::size_t isf_j = sn.stateful_of_station(j);
            T probSum = zero;
            for (std::size_t a = 0; a < nSub; ++a) {
                const std::size_t isf_i = sn.stateful_of_station(subsetIndices[a]);
                probSum += T(iso.visits(a, r) * sn.rt((isf_i - 1) * K + r, (isf_j - 1) * K + r));
            }
            if (num_traits<T>::to_double(probSum) > fineTol)
                Pk[fesNode][complementNode[j]] = probSum;
        }

        // No FES self-loop: the internal circulation is already inside the
        // state-dependent rate, and adding it would count those hops twice.
        for (std::size_t nd = 1; nd <= I; ++nd) {
            T rowSum = zero;
            for (std::size_t md = 1; md <= I; ++md) rowSum += Pk[nd][md];
            if (num_traits<T>::to_double(rowSum) > fineTol)
                for (std::size_t md = 1; md <= I; ++md) Pk[nd][md] /= rowSum;
        }
        for (std::size_t nd = 1; nd <= I; ++nd)
            for (std::size_t md = 1; md <= I; ++md)
                if (num_traits<T>::to_double(Pk[nd][md]) != 0.0)
                    P.set(newClass[r], newClass[r], nd, md, Pk[nd][md]);
    }
    fesModel.link(P);

    FesAggregateResult<T> out(fesModel);
    out.fesNode = fesNode;
    out.deagg.subsetIndices = subsetIndices;
    out.deagg.complementIndices = complementIndices;
    out.deagg.throughputTable = scalingTable;
    out.deagg.cutoffs = cutoffs;
    out.deagg.stochCompSubset = stochCompSubset;
    out.deagg.stochCompComplement = stochCompComplement;
    out.deagg.isolatedDemands = iso.L;
    out.deagg.isolatedVisits = iso.visits;
    out.deagg.isolatedServers = mi;
    out.deagg.isolatedIsDelay = isDelay;
    out.deagg.escape = escape;
    out.deagg.fesNode = fesNode;
    return out;
}

}  // namespace fes
}  // namespace line

#endif  // LINE_API_FES_FES_AGGREGATE_H
