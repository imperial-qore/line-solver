/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Regression: `set_global_dependence`, the globally state-dependent rate
 * scaling phi(n) over the FULL (nstations x nclasses) population matrix, and
 * the CTMC generator that consumes it.
 *
 * WHY IT IS NOT `set_joint_dependence` WITH A WIDER ARGUMENT. Every existing
 * hook -- lldscaling, cdscaling, jdscaling -- is handed the population of ONE
 * station (`cd_factor`, state_events.h), so none of them can express a rate
 * that reads the whole state. A Whittle network needs exactly that, and so does
 * bandwidth sharing, where one route holds several links at once.
 *
 * The oracles, in increasing strength:
 *   1. phi == 1 must reproduce the model declared without it, bit for bit.
 *   2. a phi reproducing a per-station alpha_i(n_i) must equal the same model
 *      declared through `set_load_dependence`, which pins the fold point.
 *   3. two stations sharing one unit of capacity by phi_s(n) = n_s/|n| is
 *      multiclass processor sharing, whose closed form is exact.
 *   4. the same balanced model with Erlang service of equal mean must give the
 *      same means: INSENSITIVITY, the defining property of a Whittle network,
 *      and the only oracle here that fails if PHASE transitions go unscaled.
 *
 * The MATLAB twin is line-test.git/test_ctmc_global_dependence.m, the JAR twin
 * is SolverCTMCGlobalDependenceTest.java and the python twin is
 * python/tests/test_ctmc_global_dependence.py.
 */
#include <algorithm>
#include <cmath>
#include <vector>

#include "doctest.h"
#include "line/api/sn/sn_gd_balance.h"
#include "line/lang/qn/network_builder.h"
#include "line/solvers/ctmc/solver_ctmc_analyzer.h"
#include "line/io/network_reader.h"
#include "line/io/network_writer.h"
#include "line/solvers/ssa/ssa_dispatch.h"

namespace qn = line::qn;
namespace ctmc = line::ctmc;
using line::lang::SchedStrategy;
using Dist = line::lang::Distrib<double>;

namespace {

enum class Gd { None, Identity, LoadEquivalent, Shared };

/**
 * Closed cyclic PS pair, one class. `erlang` switches the service to an Erlang-3
 * of the same mean, which changes nothing under a balanced phi and everything
 * under a broken one.
 */
qn::Network<double> closed_pair(Gd gd, double njobs, bool erlang = false) {
    qn::Network<double> m("gd_pair");
    const std::size_t q1 = m.add_queue("Q1", SchedStrategy::PS);
    const std::size_t q2 = m.add_queue("Q2", SchedStrategy::PS);
    const std::size_t c = m.add_closed_class("C1", njobs, q1);
    if (erlang) {
        m.set_service(q1, c, Dist::erlang(3.0 / 1.0, 3));
        m.set_service(q2, c, Dist::erlang(3.0 / 0.5, 3));
    } else {
        m.set_service(q1, c, Dist::exp_rate(1.0));
        m.set_service(q2, c, Dist::exp_rate(2.0));
    }
    qn::RoutingMatrix<double> P;
    P.set(q1, q2, 1.0);
    P.set(q2, q1, 1.0);
    m.link(P);
    // The struct must exist before the handles below can index stations by row.
    const std::size_t M = m.get_struct().nstations;
    switch (gd) {
        case Gd::None:
            break;
        case Gd::Identity:
            m.set_global_dependence(
                [](const std::vector<double>&) { return std::vector<double>(1, 1.0); },
                std::vector<double>(1, 1.0));
            break;
        case Gd::LoadEquivalent:
            // alpha(n) = min(n, 2) at station 2 only, i.e. a per-station rule
            // expressed through the global handle
            m.set_global_dependence(
                [M](const std::vector<double>& n) {
                    std::vector<double> v(M, 1.0);
                    if (M >= 2) v[1] = std::min(n[1], 2.0);
                    if (M >= 2 && n[1] == 0.0) v[1] = 1.0;
                    return v;
                },
                std::vector<double>(1, 2.0));
            break;
        case Gd::Shared:
            // phi_s(n) = n_s/|n|: the two stations share one unit of capacity.
            // Balanced, hence reversible, product form and insensitive.
            m.set_global_dependence(
                [M](const std::vector<double>& n) {
                    std::vector<double> v(M, 1.0);
                    double tot = 0;
                    for (std::size_t i = 0; i < M; ++i) tot += n[i];
                    if (tot > 0)
                        for (std::size_t i = 0; i < M; ++i) v[i] = n[i] / tot;
                    return v;
                },
                std::vector<double>(1, 1.0));
            break;
    }
    return m;
}

}  // namespace

TEST_CASE("ctmc global dependence: an identity phi changes nothing") {
    const ctmc::CtmcOptions opt;
    qn::Network<double> plain = closed_pair(Gd::None, 3);
    qn::Network<double> withgd = closed_pair(Gd::Identity, 3);
    const ctmc::CtmcSolution<double> d0 = ctmc::solver_ctmc_analyzer(plain.get_struct(), opt);
    const ctmc::CtmcSolution<double> d1 = ctmc::solver_ctmc_analyzer(withgd.get_struct(), opt);
    for (std::size_t i = 0; i < plain.get_struct().nstations; ++i) {
        CHECK(d1.avg.QN(i, 0) == doctest::Approx(d0.avg.QN(i, 0)).epsilon(1e-12));
        CHECK(d1.avg.TN(i, 0) == doctest::Approx(d0.avg.TN(i, 0)).epsilon(1e-12));
    }
}

TEST_CASE("ctmc global dependence: reproduces a per-station load dependence") {
    const ctmc::CtmcOptions opt;
    qn::Network<double> ld = closed_pair(Gd::None, 3);
    // alpha(n) = min(n, 2) at Q2, the local way
    ld.set_load_dependence(2, std::vector<double>{1.0, 2.0, 2.0, 2.0});
    qn::Network<double> gd = closed_pair(Gd::LoadEquivalent, 3);
    const ctmc::CtmcSolution<double> a = ctmc::solver_ctmc_analyzer(ld.get_struct(), opt);
    const ctmc::CtmcSolution<double> b = ctmc::solver_ctmc_analyzer(gd.get_struct(), opt);
    for (std::size_t i = 0; i < ld.get_struct().nstations; ++i)
        CHECK(b.avg.QN(i, 0) == doctest::Approx(a.avg.QN(i, 0)).epsilon(1e-9));
}

TEST_CASE("ctmc global dependence: a shared link is processor sharing, and is insensitive") {
    const ctmc::CtmcOptions opt;
    const double N = 3;
    qn::Network<double> mexp = closed_pair(Gd::Shared, N);
    qn::Network<double> merl = closed_pair(Gd::Shared, N, true);
    const ctmc::CtmcSolution<double> a = ctmc::solver_ctmc_analyzer(mexp.get_struct(), opt);
    const ctmc::CtmcSolution<double> b = ctmc::solver_ctmc_analyzer(merl.get_struct(), opt);

    double tot = 0;
    for (std::size_t i = 0; i < mexp.get_struct().nstations; ++i) tot += a.avg.QN(i, 0);
    CHECK(tot == doctest::Approx(N).epsilon(1e-9));

    // Insensitivity: an Erlang-3 of the same mean must not move the means. This
    // is what fails if PHASE transitions are left unscaled.
    for (std::size_t i = 0; i < mexp.get_struct().nstations; ++i)
        CHECK(b.avg.QN(i, 0) == doctest::Approx(a.avg.QN(i, 0)).epsilon(1e-6));
}

TEST_CASE("ctmc global dependence: a malformed declaration is refused at declaration time") {
    qn::Network<double> m = closed_pair(Gd::None, 3);
    // no peak
    CHECK_THROWS(m.set_global_dependence(
        [](const std::vector<double>&) { return std::vector<double>(1, 1.0); },
        std::vector<double>()));
    // non-positive peak
    CHECK_THROWS(m.set_global_dependence(
        [](const std::vector<double>&) { return std::vector<double>(1, 1.0); },
        std::vector<double>(1, -1.0)));
    // wrong output shape
    CHECK_THROWS(m.set_global_dependence(
        [](const std::vector<double>&) { return std::vector<double>(7, 1.0); },
        std::vector<double>(1, 1.0)));
}

TEST_CASE("ctmc global dependence: the balance property is checkable") {
    // phi_s(n) = n_s/|n| is balanced, so a Whittle network; n_s/(1+n_{s+1}) is
    // not, and is exactly the kind of plausible-looking rule that yields no
    // product form and no insensitivity.
    const std::function<std::vector<double>(const std::vector<double>&)> balanced =
        [](const std::vector<double>& n) {
            std::vector<double> v(n.size(), 0.0);
            double tot = 0;
            for (std::size_t i = 0; i < n.size(); ++i) tot += n[i];
            if (tot > 0)
                for (std::size_t i = 0; i < n.size(); ++i) v[i] = n[i] / tot;
            return v;
        };
    const std::function<std::vector<double>(const std::vector<double>&)> unbalanced =
        [](const std::vector<double>& n) {
            std::vector<double> v(n.size(), 0.0);
            for (std::size_t i = 0; i < n.size(); ++i)
                v[i] = n[i] / (1.0 + n[(i + 1) % n.size()]);
            return v;
        };
    const std::vector<std::size_t> cut(2, 3);
    CHECK(line::sn::sn_gd_balance<double>(balanced, cut) < 1e-12);
    CHECK(line::sn::sn_gd_balance<double>(unbalanced, cut) > 1e-3);
}


namespace {

std::size_t bf_code(const std::vector<std::size_t>& n, std::size_t base) {
    std::size_t idx = 0, mult = 1;
    for (std::size_t s = 0; s < n.size(); ++s) {
        idx += n[s] * mult;
        mult *= base;
    }
    return idx;
}

/** Phi(n) = max_l (1/C_l) sum_{s in l} Phi(n-e_s), Phi(0) = 1. */
std::vector<double> bf_phi(const std::vector<std::vector<int>>& A, const std::vector<double>& C,
                           std::size_t cutoff) {
    const std::size_t S = A[0].size(), base = cutoff + 1;
    std::size_t ns = 1;
    for (std::size_t s = 0; s < S; ++s) ns *= base;
    std::vector<double> Phi(ns, 0.0);
    Phi[0] = 1.0;
    std::vector<std::pair<std::size_t, std::size_t>> order;
    for (std::size_t k = 0; k < ns; ++k) {
        std::size_t rem = k, tot = 0;
        for (std::size_t s = 0; s < S; ++s) {
            tot += rem % base;
            rem /= base;
        }
        order.push_back(std::make_pair(tot, k));
    }
    std::stable_sort(order.begin(), order.end());
    for (std::size_t oi = 0; oi < order.size(); ++oi) {
        const std::size_t k = order[oi].second;
        std::vector<std::size_t> n(S, 0);
        std::size_t rem = k, tot = 0;
        for (std::size_t s = 0; s < S; ++s) {
            n[s] = rem % base;
            rem /= base;
            tot += n[s];
        }
        if (tot == 0) continue;
        double best = 0;
        for (std::size_t l = 0; l < A.size(); ++l) {
            double acc = 0;
            for (std::size_t s = 0; s < S; ++s)
                if (A[l][s] > 0 && n[s] > 0) {
                    n[s] -= 1;
                    acc += Phi[bf_code(n, base)];
                    n[s] += 1;
                }
            best = std::max(best, acc / C[l]);
        }
        Phi[k] = best;
    }
    return Phi;
}

}  // namespace

TEST_CASE("ctmc global dependence: open bandwidth sharing matches the product form") {
    // The 2-link linear network: route 1 crosses both links, routes 2 and 3 use
    // one link each. No per-station scaling can express that coupling. Under
    // balanced fairness the chain is reversible, so truncating it preserves the
    // conditional law and the CTMC must reproduce pi(n) ~ Phi(n) prod rho^n
    // EXACTLY rather than approximately.
    const std::vector<std::vector<int>> A = {{1, 1, 0}, {1, 0, 1}};
    const std::vector<double> C = {1.0, 1.0};
    const double nu[3] = {0.20, 0.30, 0.30}, mu[3] = {1.0, 1.0, 1.0};
    const std::size_t CUT = 3, S = 3, base = CUT + 1;
    const std::vector<double> Phi = bf_phi(A, C, CUT);

    qn::Network<double> m("gd_bws");
    const std::size_t src = m.add_source("Source");
    std::vector<std::size_t> rq(S), cl(S);
    for (std::size_t s = 0; s < S; ++s)
        rq[s] = m.add_queue("Route" + std::to_string(s + 1), SchedStrategy::PS);
    const std::size_t snk = m.add_sink("Sink");
    for (std::size_t s = 0; s < S; ++s) {
        cl[s] = m.add_open_class("F" + std::to_string(s + 1));
        m.set_arrival(src, cl[s], Dist::exp_rate(nu[s]));
    }
    for (std::size_t s = 0; s < S; ++s)
        for (std::size_t t = 0; t < S; ++t)
            m.set_service(rq[t], cl[s], s == t ? Dist::exp_rate(mu[s]) : Dist::disabled_dist());
    qn::RoutingMatrix<double> P;
    for (std::size_t s = 0; s < S; ++s) {
        P.set(cl[s], cl[s], src, rq[s], 1.0);
        P.set(cl[s], cl[s], rq[s], snk, 1.0);
    }
    m.link(P);

    const qn::NetworkStruct<double>& sn0 = m.get_struct();
    std::vector<std::size_t> idx(S);
    for (std::size_t s = 0; s < S; ++s) idx[s] = sn0.nodes[rq[s] - 1].station - 1;
    const std::size_t M = sn0.nstations, K = sn0.nclasses;

    m.set_global_dependence(
        [Phi, idx, M, K, S, base](const std::vector<double>& n) {
            std::vector<std::size_t> np(S, 0);
            std::size_t tot = 0;
            for (std::size_t s = 0; s < S; ++s) {
                np[s] = static_cast<std::size_t>(n[idx[s] * K + s] + 0.5);
                tot += np[s];
            }
            std::vector<double> v(M * K, 1.0);
            if (tot == 0) return v;
            const double den = Phi[bf_code(np, base)];
            for (std::size_t s = 0; s < S; ++s) {
                if (np[s] > 0) {
                    np[s] -= 1;
                    v[idx[s] * K + s] = Phi[bf_code(np, base)] / den;
                    np[s] += 1;
                } else {
                    v[idx[s] * K + s] = 0.0;
                }
            }
            return v;
        },
        std::vector<double>{1.0});

    ctmc::CtmcOptions opt;
    opt.cutoff = static_cast<double>(CUT);
    const ctmc::CtmcSolution<double> d = ctmc::solver_ctmc_analyzer(m.get_struct(), opt);

    // closed-form means from the same Phi
    std::size_t ns = 1;
    for (std::size_t s = 0; s < S; ++s) ns *= base;
    std::vector<double> w(ns, 0.0), En(S, 0.0);
    double wsum = 0;
    for (std::size_t k = 0; k < ns; ++k) {
        std::vector<std::size_t> n(S, 0);
        std::size_t rem = k;
        double prod = 1.0;
        for (std::size_t s = 0; s < S; ++s) {
            n[s] = rem % base;
            rem /= base;
            prod *= std::pow(nu[s] / mu[s], static_cast<double>(n[s]));
        }
        w[k] = Phi[k] * prod;
        wsum += w[k];
    }
    for (std::size_t k = 0; k < ns; ++k) {
        std::size_t rem = k;
        for (std::size_t s = 0; s < S; ++s) {
            En[s] += (w[k] / wsum) * static_cast<double>(rem % base);
            rem /= base;
        }
    }

    for (std::size_t s = 0; s < S; ++s)
        CHECK(d.avg.QN(idx[s], s) == doctest::Approx(En[s]).epsilon(1e-9));
    // the values every codebase agrees on at this cutoff
    CHECK(En[0] == doctest::Approx(0.468423).epsilon(1e-5));
    CHECK(En[1] == doctest::Approx(0.512029).epsilon(1e-5));
}

TEST_CASE("ssa global dependence: the serial engine reproduces the CTMC means") {
    // SolverSSA carries the same factorization on the sample path: phi(n) is a
    // constant within a state, so it is evaluated once per state and multiplies
    // every station service rate there. The NRM cannot (its propensity closures
    // see one station's population slice), so the dispatcher must route a
    // gd-bearing model to the serial engine rather than simulate it unscaled --
    // which is what this checks, since an unscaled run gives QLen = [1.5, 1.5]
    // against the exact [2, 1].
    qn::Network<double> m = closed_pair(Gd::Shared, 3);
    const ctmc::CtmcSolution<double> exact =
        ctmc::solver_ctmc_analyzer(m.get_struct(), ctmc::CtmcOptions());

    line::ssa::SsaOptions opt;
    opt.samples = 200000;
    opt.seed = 23000;
    const line::ssa::SsaSolution sim = line::ssa::solver_ssa(m.get_struct(), opt);
    CHECK(sim.method.find("serial") != std::string::npos);
    for (std::size_t i = 0; i < m.get_struct().nstations; ++i)
        CHECK(sim.QN(i, 0) == doctest::Approx(exact.avg.QN(i, 0)).epsilon(2e-2));
}

TEST_CASE("ssa global dependence: method='nrm' is refused by name") {
    qn::Network<double> m = closed_pair(Gd::Shared, 3);
    line::ssa::SsaOptions opt;
    opt.method = "nrm";
    opt.samples = 1000;
    CHECK_THROWS_AS(line::ssa::solver_ssa(m.get_struct(), opt), line::UnsupportedError);
}

TEST_CASE("json wire: a global dependence round-trips through model.json") {
    // phi(n) is a handle, so it reaches the wire only as a TABLE: the writer
    // materializes it over the lattice of the whole network state, restricted to
    // the (station,class) slots a class can occupy. The oracle is the answer, not
    // the document -- a table read back at the wrong coordinate still parses.
    qn::Network<double> m = closed_pair(Gd::Shared, 3);
    const line::io::detail::json doc = line::io::network_json_envelope(m.get_struct());
    const line::io::detail::json& blk = doc.at("model").at("globalDependence");
    CHECK(blk.at("type").get<std::string>() == "globalDependent");
    CHECK(blk.at("slots").size() == 2u);
    CHECK(blk.at("scaling").size() == 16u);  // the (3+1)^2 box over the two slots

    qn::Network<double> back = line::io::build_network_from_json<double>(doc);
    CHECK(static_cast<bool>(back.get_struct().gdscaling));
    const ctmc::CtmcSolution<double> a =
        ctmc::solver_ctmc_analyzer(m.get_struct(), ctmc::CtmcOptions());
    const ctmc::CtmcSolution<double> b =
        ctmc::solver_ctmc_analyzer(back.get_struct(), ctmc::CtmcOptions());
    for (std::size_t i = 0; i < a.avg.QN.rows(); ++i)
        CHECK(b.avg.QN(i, 0) == doctest::Approx(a.avg.QN(i, 0)).epsilon(1e-10));
}
