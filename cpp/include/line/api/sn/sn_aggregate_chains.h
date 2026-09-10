/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_SN_SN_AGGREGATE_CHAINS_H
#define LINE_API_SN_SN_AGGREGATE_CHAINS_H

/**
 * Collapse every chain onto one class, port of `ModelAdapter.aggregateChains`.
 *
 * The reference is `matlab/src/io/@@ModelAdapter/aggregateChains.m`, reached
 * from the model as `model.aggregateChains(suffix)`; the JAR twin is
 * `ModelAdapter.aggregateChains` and the python one
 * `ModelAdapter.aggregate_chains`.
 *
 * A CHAIN is a set of classes that can switch into one another, so the jobs of
 * a chain are one circulating population however many class labels they wear on
 * the way round. The transform replaces each chain by a SINGLE class carrying
 * that population (or that total arrival rate, for an open chain), the chain's
 * service demands and the chain-level routing. The result has no class
 * switching at all, which is what makes it worth building: the state space of a
 * multiclass model is the product over classes, and a model whose K classes
 * form C < K chains solves in the state space of C.
 *
 * WHAT IS PRESERVED, and it is chain-level and not class-level. Per chain: the
 * population, the arrival rate, the per-station service DEMAND, and the routing
 * between stations. Per class: nothing -- that is the point, and it is what the
 * `alpha` matrix returned here exists to undo. `alpha(i,r)` is class r's share
 * of its chain's visits at station i, so a chain result at station i splits
 * back over its classes in those proportions (see the reference's
 * `sn_deaggregate_chain_results`, which consumes exactly `Lchain`, `STchain`,
 * `Vchain` and `alpha` from the block returned here).
 *
 * THE AGGREGATE IS EXACT ON A PRODUCT-FORM MODEL and an approximation
 * otherwise, for one reason: the chain's service law at a station is refitted
 * from the chain MEAN and SCV, mixing the classes' own laws. Where the model is
 * product form the station is insensitive to everything past the mean and the
 * refit costs nothing; where it is not, the refit is a two-moment
 * approximation. The fit follows the reference's ladder exactly --
 *
 *   |SCV - 1| < FineTol   Exp(1/ST)              the insensitive case
 *   SCV < FineTol         Det(ST)                the degenerate lower end
 *   SCV < 1               Erlang(round(1/SCV))   the reference rounds, not ceils
 *   SCV > 1               HyperExp fitted to (ST, SCV)
 *
 * -- and the Erlang order is `round(1/SCV)` here because that is what
 * `aggregateChains.m` writes; `Erlang.fitMeanAndSCV` ceils instead, and the two
 * disagree at every SCV that is not the reciprocal of an integer. The first
 * three rungs are exact-arithmetic clean; the HyperExp one needs a square root,
 * so at a non-transcendental T it is FITTED IN DOUBLE and its three parameters
 * lifted back, which costs nothing that was exact to begin with -- the rung is a
 * two-moment approximation of a mixture and carries no exactness claim, while the
 * routing, the demands, the populations and the aggregation around it stay exact.
 *
 * THE ROUTING IS READ OFF `rt`, NOT `rtnodes`, and that is load bearing. `rt`
 * is the stochastic complement of the routing over the STATEFUL nodes, so the
 * nodes the aggregate does not carry -- Routers, and the ClassSwitch nodes that
 * `link` synthesizes for exactly the switching this transform is eliminating --
 * are already folded away, along with the Sink, which is not stateful in this
 * port nor in the reference. That last one is why an open chain still routes:
 * a `Q -> Sink -> Source` path in the original appears in `rt` as `Q -> Source`
 * (the sink closure is inside the complement), so the aggregate's open chain
 * returns to the Source and the refresh re-derives its own closure. Measured on
 * a two-class open chain in MATLAB 2026-08-15: the aggregate routes
 * Source -> Q1 -> Q2 -> Source and reports the exact per-station utilizations.
 *
 * ONE DELIBERATE DEPARTURE FROM MATLAB. The reference `line_warning`s and SKIPS
 * a node kind it does not carry (anything that is not Source, Sink, Queue,
 * Delay, Router or ClassSwitch), which silently returns a model with a node
 * missing and its routing rewired around the hole. This port refuses by name
 * instead, as the rest of the port does: a model the transform cannot carry is
 * a diagnostic, not a quietly different model. A ClassSwitch is still dropped
 * rather than refused, since eliminating class switching is the transform's
 * whole purpose and the reference drops it on both of its branches.
 *
 * ARITHMETIC: rational in the routing and the visits; the HyperExp fit is
 * transcendental (`hyperexp_fit_mean_scv`) and drops to a double fit at an exact
 * T, as `network_reader.h` already does for every T.
 */

#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

#include "line/lang/dist_fitters.h"
#include "line/lang/lang_types.h"
#include "line/lang/qn/network_builder.h"
#include "line/lang/qn/network_struct.h"
#include "line/num/number.h"
#include "line/solvers/mva/sn_chain.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace api {

/**
 * Everything needed to map a chain-level result back onto the classes.
 *
 * The reference's `deaggInfo`, minus its `originalSn` field: a struct is a
 * value here and the caller already holds the one it passed in, so carrying a
 * second copy would double the memory for nothing. Nothing in
 * `sn_deaggregate_chain_results` reads it either -- the documented call is
 * `sn_deaggregate_chain_results(sn, Lchain, [], STchain, Vchain, alpha, ...)`.
 */
template <class T>
struct ChainAggregationDeagg {
    Matrix<T> alpha;     ///< (M x K) class share of its chain's visits at a station
    Matrix<T> Lchain;    ///< (M x C) chain demand
    Matrix<T> STchain;   ///< (M x C) chain mean service time
    Matrix<T> Vchain;    ///< (M x C) chain visits
    Matrix<T> SCVchain;  ///< (M x C)
    std::vector<double> Nchain;      ///< (C) population, infinite for an open chain
    std::vector<T> lambdachain;      ///< (C) aggregate arrival rate, 0 for a closed chain
    std::vector<bool> isopenchain;   ///< (C)
    std::vector<std::vector<std::size_t> > inchain;  ///< (C) 1-based class indices
    std::vector<std::size_t> refstat;                ///< (K) 1-based reference station
    std::vector<std::size_t> refstatchain;           ///< (C) 1-based reference station
    bool isaggregated = false;  ///< false when C == K and the model was merely copied
    std::size_t nclasses = 0;   ///< K of the ORIGINAL model
    std::size_t nchains = 0;    ///< C, the class count of the aggregate
};

/** What `sn_aggregate_chains` returns. */
template <class T>
struct ChainAggregationResult {
    explicit ChainAggregationResult(const qn::Network<T>& m) : model(m) {}
    qn::Network<T> model;  ///< the aggregate, one class per chain
    Matrix<T> alpha;       ///< (M x K), the same matrix as `deagg.alpha`
    /** `chainclass[c]` is the 1-based class index of chain c+1 in `model`. */
    std::vector<std::size_t> chainclass;
    /** `stationnode[i]` is the aggregate's 1-based node for original station i+1, 0 if dropped. */
    std::vector<std::size_t> stationnode;
    ChainAggregationDeagg<T> deagg;
};

/**
 * @param sn     a REFRESHED model (the chains, `rt` and the rates are read)
 * @param suffix appended to each chain class name, as the reference's argument
 */
template <class T>
ChainAggregationResult<T> sn_aggregate_chains(const qn::NetworkStruct<T>& sn,
                                              const std::string& suffix = std::string()) {
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    const double fineTol = lang::GlobalConstants::FineTol;
    const std::size_t M = sn.nstations, K = sn.nclasses, C = sn.nchains;

    if (M == 0 || K == 0)
        throw InputError("sn_aggregate_chains: the model has no stations or no classes");
    if (C == 0 || sn.inchain.size() != C)
        throw InputError(
            "sn_aggregate_chains: the model is not refreshed (it declares no chains); call "
            "refresh_struct() or get_struct() first");

    // ---- C == K: nothing to merge ------------------------------------------
    //
    // The reference returns a COPY rather than running the general path, and
    // that is not an optimisation: the general path refits every service law
    // from its first two moments, so running it here would replace each
    // station's own distribution by a two-moment surrogate of itself and change
    // every non-product-form answer for no gain.
    if (C == K) {
        qn::Network<T> copy(sn.name);
        copy.raw_struct() = sn;
        ChainAggregationResult<T> out(copy);
        // MATLAB returns eye(M,K) here, which is what `alpha` degenerates to
        // only when the station and class axes happen to line up; it is carried
        // verbatim so a caller comparing against the reference sees the same
        // matrix, and `isaggregated` is what says not to trust it.
        out.alpha = Matrix<T>(M, K, zero);
        for (std::size_t i = 0; i < M && i < K; ++i) out.alpha(i, i) = one;
        out.deagg.alpha = out.alpha;
        out.deagg.inchain = sn.inchain;
        out.deagg.isaggregated = false;
        out.deagg.nclasses = K;
        out.deagg.nchains = C;
        out.chainclass.resize(C, 0);
        for (std::size_t c = 0; c < C; ++c) out.chainclass[c] = c + 1;
        out.stationnode.resize(M, 0);
        for (std::size_t i = 0; i < M; ++i) out.stationnode[i] = sn.station_to_node[i];
        return out;
    }

    const mva::ChainDemands<T> dem = mva::sn_get_demands_chain(sn);

    // ---- which chains are open, and at what aggregate rate ------------------
    std::vector<bool> isOpenChain(C, false);
    std::vector<T> lambdaChain(C, zero);
    for (std::size_t c = 0; c < C; ++c) {
        const std::vector<std::size_t>& ic = sn.inchain[c];
        for (std::size_t a = 0; a < ic.size(); ++a)
            if (!std::isfinite(sn.classes[ic[a] - 1].population)) isOpenChain[c] = true;
        if (!isOpenChain[c] || sn.sourceIdx == 0) continue;
        // The reference sums `sn.rates(sourceStation, inchain)` with 'omitnan';
        // a disabled arrival is the NaN it is omitting.
        for (std::size_t a = 0; a < ic.size(); ++a) {
            const std::size_t rr = ic[a];
            if (sn.disabled[sn.sourceIdx - 1][rr - 1]) continue;
            const double v = num_traits<T>::to_double(sn.rates(sn.sourceIdx - 1, rr - 1));
            if (std::isfinite(v)) lambdaChain[c] += sn.rates(sn.sourceIdx - 1, rr - 1);
        }
    }

    // ---- the nodes the aggregate carries ------------------------------------
    qn::Network<T> chainModel(sn.name.empty() ? std::string("aggregated")
                                              : sn.name + "_aggregated");
    std::vector<std::size_t> nodeMap(sn.nodes.size() + 1, 0);  // old 1-based node -> new node
    std::size_t sourceNode = 0;
    for (std::size_t i = 1; i <= sn.nodes.size(); ++i) {
        const qn::NodeDef& nd = sn.nodes[i - 1];
        switch (nd.nodetype) {
            case qn::NodeType::Source:
                nodeMap[i] = chainModel.add_source(nd.name);
                sourceNode = nodeMap[i];
                break;
            case qn::NodeType::Sink: nodeMap[i] = chainModel.add_sink(nd.name); break;
            case qn::NodeType::Delay: nodeMap[i] = chainModel.add_delay(nd.name); break;
            case qn::NodeType::Router: nodeMap[i] = chainModel.add_router(nd.name); break;
            case qn::NodeType::Queue: {
                const qn::Station<T>& st = sn.stations[nd.station - 1];
                nodeMap[i] = chainModel.add_queue(nd.name, st.sched);
                if (!std::isinf(st.nservers))
                    chainModel.set_number_of_servers(nodeMap[i], st.nservers);
                if (st.cap > 0.0 && std::isfinite(st.cap))
                    chainModel.set_capacity(nodeMap[i], st.cap);
                break;
            }
            case qn::NodeType::ClassSwitch:
                // DROPPED, not refused: the aggregate has one class per chain
                // and therefore nothing left to switch. This covers both the
                // user's own ClassSwitch nodes and the ones `link` synthesizes,
                // which is the reference's behaviour on both of its branches.
                break;
            default:
                throw UnsupportedError(
                    std::string("sn_aggregate_chains: node '") + nd.name + "' is a " +
                    lang::node_type_to_text(nd.nodetype) +
                    ", which the chain aggregate cannot carry; the transform keeps Source, Sink, "
                    "Queue, Delay and Router and drops ClassSwitch");
        }
    }

    std::vector<std::size_t> stationNode(M + 1, 0);  // old 1-based station -> new node
    for (std::size_t i = 1; i <= M; ++i) stationNode[i] = nodeMap[sn.station_to_node[i - 1]];

    // ---- one class per chain -------------------------------------------------
    std::vector<std::size_t> chainClass(C, 0);
    for (std::size_t c = 0; c < C; ++c) {
        const std::vector<std::size_t>& ic = sn.inchain[c];
        std::string nm = (ic.size() == 1) ? sn.classes[ic[0] - 1].name
                                          : (std::string("Chain") + std::to_string(c + 1));
        if (!suffix.empty()) nm += suffix;
        if (isOpenChain[c]) {
            chainClass[c] = chainModel.add_open_class(nm);
        } else {
            const std::size_t refst = dem.refstatchain[c];
            if (refst == 0 || refst > M || stationNode[refst] == 0)
                throw InputError("sn_aggregate_chains: the reference station of chain " +
                                 std::to_string(c + 1) + " is not carried by the aggregate");
            chainClass[c] = chainModel.add_closed_class(nm, dem.Nchain[c], stationNode[refst]);
        }
    }

    // ---- the aggregate arrival rates ----------------------------------------
    for (std::size_t c = 0; c < C; ++c) {
        if (!isOpenChain[c]) continue;
        if (sourceNode == 0)
            throw InputError(
                "sn_aggregate_chains: chain " + std::to_string(c + 1) +
                " is open but the model carries no Source to attach its arrivals to");
        if (num_traits<T>::to_double(lambdaChain[c]) > 0.0)
            chainModel.set_arrival(sourceNode, chainClass[c],
                                   lang::Distrib<T>::exp_rate(lambdaChain[c]));
    }

    // ---- the chain service laws ---------------------------------------------
    for (std::size_t i = 1; i <= M; ++i) {
        if (stationNode[i] == 0) continue;
        const qn::NodeType ty = sn.nodes[sn.station_to_node[i - 1] - 1].nodetype;
        if (ty == qn::NodeType::Source || ty == qn::NodeType::Sink) continue;
        for (std::size_t c = 0; c < C; ++c) {
            const double st = num_traits<T>::to_double(dem.STchain(i - 1, c));
            if (!(st > 0.0) || !std::isfinite(st)) {
                // No service for this chain here. Left at the table's own
                // Disabled, which is what the reference sets explicitly.
                chainModel.set_service(stationNode[i], chainClass[c],
                                       lang::Distrib<T>::disabled_dist());
                continue;
            }
            double scv = num_traits<T>::to_double(dem.SCVchain(i - 1, c));
            if (!std::isfinite(scv) || scv <= 0.0) scv = 1.0;  // default to exponential
            const T mean = dem.STchain(i - 1, c);
            lang::Distrib<T> d;
            if (std::fabs(scv - 1.0) < fineTol) {
                d = lang::Distrib<T>::exp_mean(mean);
            } else if (scv < 1.0) {
                if (scv < fineTol) {
                    d = lang::Distrib<T>::det(mean);
                } else {
                    // `Erlang.fitMeanAndOrder(mean, k)` is Erlang(k/mean, k);
                    // the reference ROUNDS 1/SCV here rather than ceiling it.
                    long k = std::lround(1.0 / scv);
                    if (k < 1) k = 1;
                    const T ph = T(num_traits<T>::from_int(k) / mean);
                    d = lang::Distrib<T>::erlang(ph, static_cast<std::size_t>(k));
                }
            } else {
                // THE HYPEREXP RUNG IS FITTED IN DOUBLE AT EXACT ARITHMETIC, and it
                // has to be selected at COMPILE time. `hyperexp_fit_mean_scv`
                // static_asserts on transcendental arithmetic (the moment
                // discriminant is a square root), and a template instantiates every
                // branch of this ladder whatever the runtime SCV is -- so calling it
                // unguarded made the whole of `sn_aggregate_chains`, and every caller
                // of it up to `solver_ctmc_run_analyzer`, uninstantiable at exact arithmetic,
                // including the Exp, Det and Erlang rungs that never reach here.
                //
                // Fitting in double and lifting the three parameters is what
                // `network_reader.h::hyperexp_fit_mean_scv` already does for every T,
                // and it costs nothing that was exact to begin with: this rung IS a
                // two-moment approximation of a chain that mixes several classes'
                // service laws, so its parameters carry no exactness claim. What
                // stays exact is everything around it -- the routing, the demands,
                // the populations and the aggregation itself. The alternative,
                // refusing, would leave a chain with SCV > 1 with no exact-arithmetic
                // path at all.
                if constexpr (!num_traits<T>::has_transcendental) {
                    // params = (p, lambda1, lambda2), the order `Distrib::hyperexp`
                    // pushes them in.
                    const lang::Distrib<double> dd = lang::hyperexp_fit_mean_scv<double>(
                        num_traits<T>::to_double(mean), scv);
                    d = lang::Distrib<T>::hyperexp(num_traits<T>::from_double(dd.params[0]),
                                                   num_traits<T>::from_double(dd.params[1]),
                                                   num_traits<T>::from_double(dd.params[2]));
                } else {
                    d = lang::hyperexp_fit_mean_scv(mean, num_traits<T>::from_double(scv));
                }
            }
            chainModel.set_service(stationNode[i], chainClass[c], d);
        }
    }

    // ---- the chain routing ---------------------------------------------------
    //
    // p_chain(i -> j) = sum_{k,s in chain} alpha(i,k) * rt((i,k) -> (j,s)),
    // i.e. the class-level routing weighted by how the chain's visits at i are
    // shared among its classes. `rt` is indexed (stateful-1)*K + class.
    qn::RoutingMatrix<T> P;
    for (std::size_t c = 0; c < C; ++c) {
        const std::vector<std::size_t>& ic = sn.inchain[c];
        Matrix<T> Pc(chainModel.raw_struct().nodes.size() + 1,
                     chainModel.raw_struct().nodes.size() + 1, zero);
        for (std::size_t i = 1; i <= M; ++i) {
            if (stationNode[i] == 0) continue;
            const qn::NodeType tyi = sn.nodes[sn.station_to_node[i - 1] - 1].nodetype;
            // A closed chain neither leaves through the Sink nor enters at the
            // Source, so those rows are not its routing.
            if (!isOpenChain[c] &&
                (tyi == qn::NodeType::Source || tyi == qn::NodeType::Sink))
                continue;
            const std::size_t isf_i = sn.stateful_of_station(i);
            for (std::size_t j = 1; j <= M; ++j) {
                if (stationNode[j] == 0) continue;
                const std::size_t isf_j = sn.stateful_of_station(j);
                T pij = zero;
                for (std::size_t a = 0; a < ic.size(); ++a) {
                    const std::size_t k = ic[a];
                    if (!(num_traits<T>::to_double(dem.alpha(i - 1, k - 1)) > 0.0)) continue;
                    for (std::size_t b = 0; b < ic.size(); ++b) {
                        const std::size_t s = ic[b];
                        const std::size_t from = (isf_i - 1) * K + k - 1;
                        const std::size_t to = (isf_j - 1) * K + s - 1;
                        if (from >= sn.rt.rows() || to >= sn.rt.cols()) continue;
                        const T p_ks = sn.rt(from, to);
                        if (num_traits<T>::to_double(p_ks) > 0.0)
                            pij += T(dem.alpha(i - 1, k - 1) * p_ks);
                    }
                }
                if (num_traits<T>::to_double(pij) > fineTol)
                    Pc(stationNode[i], stationNode[j]) = pij;
            }
        }
        // Renormalise each row. The weights already sum to one when the chain
        // decomposition and the alpha shares agree; the reference renormalises
        // anyway and warns above 1%, and the correction is kept here so a
        // rounding residue cannot turn into a leaking chain.
        for (std::size_t a = 1; a < Pc.rows(); ++a) {
            T rowSum = zero;
            for (std::size_t b = 1; b < Pc.cols(); ++b) rowSum += Pc(a, b);
            if (!(num_traits<T>::to_double(rowSum) > fineTol)) continue;
            for (std::size_t b = 1; b < Pc.cols(); ++b)
                if (num_traits<T>::to_double(Pc(a, b)) != 0.0) {
                    Pc(a, b) = T(Pc(a, b) / rowSum);
                    P.set(chainClass[c], chainClass[c], a, b, Pc(a, b));
                }
        }
    }
    chainModel.link(P);

    // ---- the result ----------------------------------------------------------
    ChainAggregationResult<T> out(chainModel);
    out.alpha = dem.alpha;
    out.chainclass = chainClass;
    out.stationnode.assign(M, 0);
    for (std::size_t i = 1; i <= M; ++i) out.stationnode[i - 1] = stationNode[i];
    out.deagg.alpha = dem.alpha;
    out.deagg.Lchain = dem.Lchain;
    out.deagg.STchain = dem.STchain;
    out.deagg.Vchain = dem.Vchain;
    out.deagg.SCVchain = dem.SCVchain;
    out.deagg.Nchain = dem.Nchain;
    out.deagg.lambdachain = lambdaChain;
    out.deagg.isopenchain = isOpenChain;
    out.deagg.inchain = sn.inchain;
    out.deagg.refstat.assign(K, 0);
    for (std::size_t r = 0; r < K; ++r) out.deagg.refstat[r] = sn.classes[r].refstat;
    out.deagg.refstatchain = dem.refstatchain;
    out.deagg.isaggregated = true;
    out.deagg.nclasses = K;
    out.deagg.nchains = C;
    return out;
}

}  // namespace api
}  // namespace line

#endif  // LINE_API_SN_SN_AGGREGATE_CHAINS_H
