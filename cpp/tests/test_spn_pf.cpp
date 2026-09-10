/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * `spn_pf`: derive the product form of a stochastic Petri net.
 *
 * The MDD-rec paper (Balsamo-Marin-Stojic, FGCS 111 (2020) 475-490) takes the
 * g_l as GIVEN -- deriving them is declared out of scope in its Sec. 3.2 -- so
 * everything the api/spn functions do downstream was, until `spn_pf`,
 * unreachable from a solver. These cases pin the derivation itself.
 *
 * FOUR INDEPENDENT ORACLES, because a product form that is merely plausible is
 * worse than none:
 *
 *  1. GLOBAL BALANCE. pi built from the derived g_l must satisfy pi Q = 0 on the
 *     generator assembled here from the net's own rate law, which shares no code
 *     path with the derivation. This is the definition of stationarity.
 *  2. THE CLOSED FORM. On the cyclic net the Gordon-Newell factors are known, so
 *     y must come out proportional to 1/mu.
 *  3. THE CERTIFICATE. Deficiency, linkage classes, stoichiometric rank and weak
 *     reversibility are structural facts readable by hand on these nets.
 *  4. THE OTHER CODEBASES. y, G and the token counts are pinned at 12 decimals,
 *     and the cyclic tokens are the value the four already agree on for
 *     `spn_rec`.
 *
 * THE MODELS are the 3-place cyclic net at N = 4 with rates {1, 1.5, 2}; the
 * same net under INFINITE-SERVER firing, which is mass action rather than a
 * constant rate and so selects the other psi; and a FORK-JOIN net, whose marking
 * is not a conserved job population at all -- a mode consumes one token and
 * produces two -- which is the case the MDD-rec paper exists to serve.
 *
 * ARITHMETIC. `spn_pf` exponentiates a least-squares solution, so it is refused
 * under exact arithmetic by name; the last case pins that refusal.
 */

#include <cmath>
#include <limits>
#include <map>
#include <string>
#include <vector>

#include "doctest.h"
#include "line/api/spn/spn_metrics.h"
#include "line/api/spn/spn_pf.h"
#include "line/lang/qn/network_builder.h"
#include "line/num/number.h"

using namespace line;
namespace qn = line::qn;
namespace spn = line::spn;
using line::lang::TimingStrategy;
using Dist = line::lang::Distrib<double>;

namespace {

const double RATES[3] = {1.0, 1.5, 2.0};
const int NJOBS = 4;

/** MATLAB, native python and the JAR on these nets, at %.12f. */
const double CYCLIC_Y[3] = {1.442249570307, 0.961499713538, 0.721124785154};
const double CYCLIC_G = 19.934426352559;
const double CYCLIC_TOKENS[3] = {2.249874392899, 1.069837548149, 0.680288058952};
const double CYCLIC_MODEX = 0.869201138838;
const double FJ_Y[4] = {1.028383350947, 1.381974961646, 1.381974961646, 0.703630713806};
const double FJ_G = 20.320471787308;
const double FJ_TOKENS[4] = {0.710262429604, 1.864198441415, 1.864198441415, 0.425539128981};
const double FJ_MODEX = 0.607360882508;

/** One timed mode consuming `from` and producing into every place in `to`. */
qn::TransitionParam<double> mode_of(std::size_t nnodes, const std::vector<std::size_t>& from,
                                    const std::vector<std::size_t>& to, double rate,
                                    double servers) {
    qn::TransitionParam<double> tp;
    tp.nmodes = 1;
    tp.modenames.push_back("fire");
    tp.enabling.assign(1, line::Matrix<double>(nnodes, 1, 0.0));
    tp.inhibiting.assign(1,
                         line::Matrix<double>(nnodes, 1, std::numeric_limits<double>::infinity()));
    tp.firing.assign(1, line::Matrix<double>(nnodes, 1, 0.0));
    for (std::size_t i = 0; i < from.size(); ++i) tp.enabling[0](from[i] - 1, 0) = 1.0;
    for (std::size_t i = 0; i < to.size(); ++i) tp.firing[0](to[i] - 1, 0) = 1.0;
    tp.nmodeservers.push_back(servers);
    tp.firingphases.push_back(1);
    tp.timing.push_back(TimingStrategy::TIMED);
    tp.fireweight.push_back(1.0);
    tp.firingproc.push_back(Dist::exp_rate(rate));
    return tp;
}

/**
 * P0 -> T0 -> P1 -> T1 -> P2 -> T2 -> P0, one token class.
 *
 * `servers` = 1 makes every mode fire at its rate constant, which is psi = 1;
 * an infinite server count makes it fire at rate lambda*m_p, which is mass
 * action.
 */
qn::Network<double> cyclic_spn(int ntokens, double servers) {
    qn::Network<double> m("spn");
    std::vector<std::size_t> pl;
    for (std::size_t i = 0; i < 3; ++i) pl.push_back(m.add_place("P" + std::to_string(i)));
    const std::size_t c = m.add_closed_class("Class1", static_cast<double>(ntokens), pl[0]);
    for (std::size_t i = 0; i < 3; ++i) m.set_service(pl[i], c, Dist::exp_rate(1.0));
    for (std::size_t i = 0; i < 3; ++i)
        m.add_transition("T" + std::to_string(i),
                         mode_of(6, std::vector<std::size_t>(1, pl[i]),
                                 std::vector<std::size_t>(1, pl[(i + 1) % 3]), RATES[i], servers));
    return m;
}

/**
 * P0 -(Tf)-> P1 + P2 -(Tj)-> P3 -(Tb)-> P0.
 *
 * Tf consumes ONE token and produces TWO, Tj the reverse, so the marking is not
 * a conserved population and the net has no queueing-network counterpart. Its
 * place invariant is (2, 1, 1, 2).
 */
qn::Network<double> forkjoin_spn(int ntokens) {
    qn::Network<double> m("fj");
    std::vector<std::size_t> pl;
    for (std::size_t i = 0; i < 4; ++i) pl.push_back(m.add_place("P" + std::to_string(i)));
    const std::size_t c = m.add_closed_class("C", static_cast<double>(ntokens), pl[0]);
    for (std::size_t i = 0; i < 4; ++i) m.set_service(pl[i], c, Dist::exp_rate(1.0));
    std::vector<std::size_t> in1(1, pl[0]);
    std::vector<std::size_t> out2;
    out2.push_back(pl[1]);
    out2.push_back(pl[2]);
    m.add_transition("Tf", mode_of(7, in1, out2, 1.3, 1.0));
    std::vector<std::size_t> in2;
    in2.push_back(pl[1]);
    in2.push_back(pl[2]);
    m.add_transition("Tj", mode_of(7, in2, std::vector<std::size_t>(1, pl[3]), 0.7, 1.0));
    m.add_transition("Tb", mode_of(7, std::vector<std::size_t>(1, pl[3]),
                                   std::vector<std::size_t>(1, pl[0]), 1.9, 1.0));
    return m;
}

/**
 * max |pi Q| over the enumerated reachable set. The generator is assembled here
 * from the net's own rate law, independently of the derivation, so agreement is
 * evidence and not a tautology.
 */
double balance_residual(const spn::SpnPfResult<double>& pf) {
    const spn::SpnInfo<double>& info = pf.spn.info;
    const std::size_t L = info.nplacelevels;
    const std::vector<std::vector<int>> st = info.diagram.enumerate();
    const std::size_t n = st.size();
    std::map<std::string, std::size_t> idx;
    for (std::size_t i = 0; i < n; ++i) {
        std::string k;
        for (std::size_t l = 0; l < L; ++l) k += std::to_string(st[i][l]) + ",";
        idx[k] = i;
    }
    std::vector<std::vector<double>> Q(n, std::vector<double>(n, 0.0));
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t e = 0; e < info.modes.size(); ++e) {
            const spn::SpnMode<double>& mde = info.modes[e];
            bool enabled = true;
            for (std::size_t l = 0; l < L && enabled; ++l)
                if (st[i][l] < mde.enab[l]) enabled = false;
            if (!enabled) continue;
            double rate = mde.D1(0, 0);
            if (pf.kind == "massaction") {
                for (std::size_t l = 0; l < L; ++l)
                    for (int j = 0; j < static_cast<int>(mde.enab[l]); ++j) rate *= (st[i][l] - j);
            } else {
                double deg = std::numeric_limits<double>::infinity();
                for (std::size_t l = 0; l < L; ++l)
                    if (mde.enab[l] > 0) deg = std::min(deg, std::floor(st[i][l] / mde.enab[l]));
                if (!std::isfinite(deg)) deg = 1.0;
                rate *= std::min(deg, mde.srv);
            }
            std::string k;
            for (std::size_t l = 0; l < L; ++l)
                k += std::to_string(
                         static_cast<int>(st[i][l] - mde.enab[l] + mde.fire[l])) + ",";
            REQUIRE(idx.count(k) == 1);
            Q[i][idx[k]] += rate;
        }
    }
    std::vector<double> pi(n, 0.0);
    double tot = 0;
    for (std::size_t i = 0; i < n; ++i) {
        double p = 1;
        for (std::size_t l = 0; l < L; ++l) p *= pf.g[l][static_cast<std::size_t>(st[i][l])];
        pi[i] = p;
        tot += p;
    }
    for (std::size_t i = 0; i < n; ++i) {
        pi[i] /= tot;
        double s = 0;
        for (std::size_t j = 0; j < n; ++j) s += Q[i][j];
        Q[i][i] -= s;
    }
    double res = 0;
    for (std::size_t j = 0; j < n; ++j) {
        double s = 0;
        for (std::size_t i = 0; i < n; ++i) s += pi[i] * Q[i][j];
        res = std::max(res, std::fabs(s));
    }
    return res;
}

}  // namespace

TEST_CASE("the cyclic net's product form is the Gordon-Newell one") {
    qn::Network<double> m = cyclic_spn(NJOBS, 1.0);
    const spn::SpnPfResult<double> pf = spn::spn_pf<double>(m.get_struct());
    CHECK(pf.kind == "geometric");
    // y is proportional to 1/mu, up to the gauge the minimum-norm solution fixes
    for (std::size_t l = 0; l < 3; ++l) {
        CHECK(pf.y[l] / pf.y[0] == doctest::Approx(RATES[0] / RATES[l]).epsilon(1e-11));
        CHECK(pf.y[l] == doctest::Approx(CYCLIC_Y[l]).epsilon(1e-11));
    }
}

TEST_CASE("the cyclic net's certificate is the structural one") {
    qn::Network<double> m = cyclic_spn(NJOBS, 1.0);
    const spn::SpnPfResult<double> pf = spn::spn_pf<double>(m.get_struct());
    // three complexes {e0, e1, e2}, one linkage class (the cycle is strongly
    // connected), stoichiometric rank 2, so deficiency 3 - 1 - 2 = 0
    CHECK(pf.complexes.size() == 3);
    CHECK(pf.linkage == 1);
    CHECK(pf.srank == 2);
    CHECK(pf.deficiency == 0);
    CHECK(pf.weakly_reversible);
    CHECK(pf.residual < 1e-12);
}

TEST_CASE("the derived law satisfies global balance") {
    {
        qn::Network<double> m = cyclic_spn(NJOBS, 1.0);
        CHECK(balance_residual(spn::spn_pf<double>(m.get_struct())) < 1e-12);
    }
    {
        qn::Network<double> m = forkjoin_spn(3);
        CHECK(balance_residual(spn::spn_pf<double>(m.get_struct())) < 1e-12);
    }
    {
        qn::Network<double> m = cyclic_spn(NJOBS, std::numeric_limits<double>::infinity());
        CHECK(balance_residual(spn::spn_pf<double>(m.get_struct())) < 1e-12);
    }
}

TEST_CASE("infinite-server firing selects mass action") {
    qn::Network<double> m = cyclic_spn(NJOBS, std::numeric_limits<double>::infinity());
    const spn::SpnPfResult<double> pf = spn::spn_pf<double>(m.get_struct());
    CHECK(pf.kind == "massaction");
    // the factors carry the 1/k! that psi = prod 1/m! puts there
    for (std::size_t l = 0; l < 3; ++l) {
        double fact = 1;
        for (std::size_t k = 0; k < pf.g[l].size(); ++k) {
            if (k > 0) fact *= static_cast<double>(k);
            CHECK(pf.g[l][k] ==
                  doctest::Approx(std::pow(pf.y[l], static_cast<double>(k)) / fact).epsilon(1e-11));
        }
    }
}

TEST_CASE("the cyclic net's values are pinned across the codebases") {
    qn::Network<double> m = cyclic_spn(NJOBS, 1.0);
    const spn::SpnPfResult<double> pf = spn::spn_pf<double>(m.get_struct());
    const spn::SpnMetrics<double> met =
        spn::spn_metrics<double>(pf.spn.mdds, pf.g, pf.spn.info);
    CHECK(met.G == doctest::Approx(CYCLIC_G).epsilon(1e-10));
    for (std::size_t l = 0; l < 3; ++l) {
        CHECK(met.tokens[l] == doctest::Approx(CYCLIC_TOKENS[l]).epsilon(1e-10));
        CHECK(met.mode_tput[l] == doctest::Approx(CYCLIC_MODEX).epsilon(1e-10));
    }
}

TEST_CASE("the fork-join net's values are pinned across the codebases") {
    qn::Network<double> m = forkjoin_spn(3);
    const spn::SpnPfResult<double> pf = spn::spn_pf<double>(m.get_struct());
    const spn::SpnMetrics<double> met =
        spn::spn_metrics<double>(pf.spn.mdds, pf.g, pf.spn.info);
    CHECK(pf.deficiency == 0);
    CHECK(pf.weakly_reversible);
    CHECK(met.G == doctest::Approx(FJ_G).epsilon(1e-10));
    for (std::size_t l = 0; l < 4; ++l) {
        CHECK(pf.y[l] == doctest::Approx(FJ_Y[l]).epsilon(1e-11));
        CHECK(met.tokens[l] == doctest::Approx(FJ_TOKENS[l]).epsilon(1e-10));
    }
    // every token that forks must later join and return, so on this cycle the
    // three modes share one throughput -- a conservation law the derivation was
    // never told about
    for (std::size_t e = 0; e < 3; ++e)
        CHECK(met.mode_tput[e] == doctest::Approx(FJ_MODEX).epsilon(1e-10));
}

TEST_CASE("the derivation is refused under exact arithmetic") {
    // y is the exponential of a least-squares solution, so there is no rational
    // answer to give; the refusal is by name rather than a silent double detour,
    // and it fires before the struct is read at all. What CONSUMES the g_l
    // afterwards stays rational -- see test_spn_rec.cpp.
    const qn::NetworkStruct<Rational> empty;
    CHECK_THROWS_AS(spn::spn_pf<Rational>(empty), UnsupportedError);
}
