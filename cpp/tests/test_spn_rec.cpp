/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * MDD-rec and the SPN measures built on it: `mdd_rec`, `mdd_rec_marginal`,
 * `spn_rec_enabled`, `spn_metrics`, `spn_sinvariants`, `spn_conv`.
 *
 * S. Balsamo, A. Marin, I. Stojic, "Computation of the normalising constant for
 * product-form models of distributed systems with synchronisation", FGCS 111
 * (2020) 475-490.
 *
 * THREE INDEPENDENT ORACLES, because a normalising constant is a single number
 * that a wrong recursion can still produce plausibly:
 *
 *  1. EXPLICIT SUM. G is also sum over the enumerated reachable set of the
 *     product of the g_l. The diagram walk and the enumeration share no code
 *     path, so agreement pins the recursion itself.
 *  2. THE CONVOLUTION. `spn_conv` decomposes {m : S m = V} instead of walking
 *     the diagram; on this net the two must agree exactly (FGCS Sec. 5.2).
 *  3. THE OTHER CODEBASES. The marginals must reproduce MATLAB's and native
 *     python's queue lengths on the same net, [2.249874392899, 1.069837548149,
 *     0.680288058952], and the mode throughputs must match the aggregation's X,
 *     which comes from an entirely different algorithm (`mdd_mcd`).
 *
 * THE MODEL is the 3-place cyclic net at N = 4 with firing rates {1, 1.5, 2}: a
 * Gordon-Newell network, so its product form is known in closed form,
 * g_l(n) = (1/mu_l)^n with unit visit ratios, and no product-form TEST is needed
 * to obtain the g_l -- which is the part the paper itself declares out of scope.
 */

#include <cmath>
#include <limits>
#include <string>
#include <vector>

#include "doctest.h"
#include "line/api/mdd/mdd.h"
#include "line/api/mdd/mdd_descriptor.h"
#include "line/api/mdd/mdd_mcd.h"
#include "line/api/mdd/mdd_rec.h"
#include "line/api/mdd/mdd_reachset.h"
#include "line/api/spn/spn_conv.h"
#include "line/api/spn/spn_mdd.h"
#include "line/api/spn/spn_metrics.h"
#include "line/api/spn/spn_rec_enabled.h"
#include "line/api/spn/spn_sinvariants.h"
#include "line/lang/qn/network_builder.h"
#include "line/num/number.h"

using namespace line;
namespace qn = line::qn;
namespace mdd = line::mdd;
namespace spn = line::spn;
using line::lang::TimingStrategy;
using Dist = line::lang::Distrib<double>;

namespace {

const double RATES[3] = {1.0, 1.5, 2.0};
const int NJOBS = 4;
/** MATLAB and native-python queue lengths on this net, at %.12f. */
const double QLEN_REF[3] = {2.249874392899, 1.069837548149, 0.680288058952};

Matrix<double> cyclic_routing(std::size_t K) {
    Matrix<double> P(K, K, 0.0);
    for (std::size_t i = 0; i < K; ++i) P(i, (i + 1) % K) = 1.0;
    return P;
}

qn::TransitionParam<double> one_mode(std::size_t nnodes, std::size_t from, std::size_t to,
                                     double rate) {
    qn::TransitionParam<double> tp;
    tp.nmodes = 1;
    tp.modenames.push_back("fire");
    tp.enabling.assign(1, line::Matrix<double>(nnodes, 1, 0.0));
    tp.inhibiting.assign(1, line::Matrix<double>(nnodes, 1, std::numeric_limits<double>::infinity()));
    tp.firing.assign(1, line::Matrix<double>(nnodes, 1, 0.0));
    tp.enabling[0](from - 1, 0) = 1.0;
    tp.firing[0](to - 1, 0) = 1.0;
    tp.nmodeservers.push_back(1.0);
    tp.firingphases.push_back(1);
    tp.timing.push_back(TimingStrategy::TIMED);
    tp.fireweight.push_back(1.0);
    tp.firingproc.push_back(Dist::exp_rate(rate));
    return tp;
}

qn::Network<double> cyclic_spn(int ntokens) {
    qn::Network<double> m("spn");
    std::vector<std::size_t> pl;
    for (std::size_t i = 0; i < 3; ++i) pl.push_back(m.add_place("P" + std::to_string(i)));
    const std::size_t c = m.add_closed_class("Class1", static_cast<double>(ntokens), pl[0]);
    for (std::size_t i = 0; i < 3; ++i) m.set_service(pl[i], c, Dist::exp_rate(1.0));
    for (std::size_t i = 0; i < 3; ++i)
        m.add_transition("T" + std::to_string(i), one_mode(6, pl[i], pl[(i + 1) % 3], RATES[i]));
    return m;
}

/** Gordon-Newell factors: g_l(n) = (v_l/mu_l)^n, unit visit ratios on a cycle. */
std::vector<std::vector<double>> gordon_newell_g(const std::vector<int>& domain) {
    std::vector<std::vector<double>> g(domain.size());
    for (std::size_t l = 0; l < domain.size(); ++l) {
        g[l].assign(static_cast<std::size_t>(domain[l]), 0.0);
        double x = 1.0;
        for (int n = 0; n < domain[l]; ++n) {
            g[l][static_cast<std::size_t>(n)] = x;
            x /= RATES[l];
        }
    }
    return g;
}

/** G by explicit enumeration, the independent oracle for the recursion. */
double explicit_norm_const(const mdd::MDD& d, const std::vector<std::vector<double>>& g) {
    const std::vector<std::vector<int>> states = d.enumerate();
    double G = 0;
    for (std::size_t s = 0; s < states.size(); ++s) {
        double p = 1.0;
        for (std::size_t l = 0; l < states[s].size(); ++l)
            p *= g[l][static_cast<std::size_t>(states[s][l])];
        G += p;
    }
    return G;
}

}  // namespace

TEST_CASE("MDD-rec reproduces the explicit normalising constant and the marginals") {
    const std::size_t K = 3;
    std::vector<double> mu(RATES, RATES + 3);
    std::vector<double> servers(K, 1.0);
    const mdd::MddDescriptor<double> desc =
        mdd::mdd_descriptor(mu, cyclic_routing(K), servers, NJOBS);
    const mdd::MDD diagram = mdd::mdd_reachset(desc.domain, desc.init, desc.nextfun);
    const mdd::MddStruct mdds = diagram.to_struct();
    const std::vector<std::vector<double>> g = gordon_newell_g(desc.domain);

    const double G = mdd::mdd_rec(mdds, g);
    CHECK(G == doctest::Approx(explicit_norm_const(diagram, g)).epsilon(1e-12));

    // The marginals are a probability law per level, and their means are the
    // queue lengths the other codebases report through a different algorithm.
    for (std::size_t l = 0; l < K; ++l) {
        const std::vector<double> mass = mdd::mdd_rec_marginal(mdds, g, l);
        double total = 0, mean = 0;
        for (std::size_t k = 0; k < mass.size(); ++k) {
            total += mass[k] / G;
            mean += static_cast<double>(k) * mass[k] / G;
        }
        CHECK(total == doctest::Approx(1.0).epsilon(1e-12));
        CHECK(mean == doctest::Approx(QLEN_REF[l]).epsilon(1e-9));
    }
}

TEST_CASE("the convolution agrees with MDD-rec on an S-invariant reachable net") {
    qn::Network<double> m = cyclic_spn(NJOBS);
    const spn::SpnResult<double> r = spn::spn_mdd(m.get_struct());
    const spn::SpnInvariants inv = spn::spn_sinvariants(m.get_struct());

    // One conservation law, "the tokens are conserved", with weight one on each
    // place: the net is a marked graph and its only minimal support is the whole
    // place set.
    REQUIRE(inv.S.size() == 1);
    CHECK(inv.S[0].size() == 3);
    for (std::size_t p = 0; p < 3; ++p) CHECK(inv.S[0][p] == 1);
    CHECK(inv.V[0] == NJOBS);

    const std::vector<std::vector<double>> g = gordon_newell_g(r.desc.domain);
    const double G_mdd = mdd::mdd_rec(r.mdds, g);
    const double G_conv = spn::spn_conv(inv, g);
    CHECK(G_conv == doctest::Approx(G_mdd).epsilon(1e-12));
}

TEST_CASE("the SPN measures agree with the level aggregation") {
    qn::Network<double> m = cyclic_spn(NJOBS);
    const spn::SpnResult<double> r = spn::spn_mdd(m.get_struct());
    const std::vector<std::vector<double>> g = gordon_newell_g(r.desc.domain);
    const spn::SpnMetrics<double> met = spn::spn_metrics(r.mdds, g, r.info);

    for (std::size_t l = 0; l < 3; ++l) {
        CHECK(met.tokens[l] == doctest::Approx(QLEN_REF[l]).epsilon(1e-9));
        // Every enabling multiplicity here is one, so the enabling degree of
        // mode l IS the marking of place l and the two utilizations coincide.
        CHECK(met.mode_util[l] == doctest::Approx(met.place_util[l]).epsilon(1e-12));
    }

    // Throughput is conserved around the cycle, and equals the throughput the
    // Miner-Ciardo-Donatelli aggregation computes for the equivalent queueing
    // network -- two algorithms with nothing in common but the model.
    std::vector<double> mu(RATES, RATES + 3);
    std::vector<double> servers(3, 1.0);
    const mdd::MddDescriptor<double> qdesc =
        mdd::mdd_descriptor(mu, cyclic_routing(3), servers, NJOBS);
    const mdd::MDD qdiag = mdd::mdd_reachset(qdesc.domain, qdesc.init, qdesc.nextfun);
    const mdd::MddMcdResult<double> qout = mdd::mdd_mcd(qdiag.to_struct(), qdesc);
    for (std::size_t e = 0; e < 3; ++e)
        CHECK(met.mode_tput[e] == doctest::Approx(qout.X[0]).epsilon(1e-9));
    // A place feeds exactly one single-server transition, so what leaves it per
    // unit time is that transition's throughput.
    for (std::size_t l = 0; l < 3; ++l)
        CHECK(met.place_tput[l] == doctest::Approx(met.mode_tput[l]).epsilon(1e-12));
}

TEST_CASE("the enabling-degree law is a proper distribution over the degrees") {
    qn::Network<double> m = cyclic_spn(NJOBS);
    const spn::SpnResult<double> r = spn::spn_mdd(m.get_struct());
    const std::vector<std::vector<double>> g = gordon_newell_g(r.desc.domain);
    const double G = mdd::mdd_rec(r.mdds, g);
    const spn::SpnEnabling<double> en = spn::spn_rec_enabled(r.mdds, g, r.info.modes[0], 3);

    CHECK(en.max_degree == static_cast<std::size_t>(NJOBS));  // one token per enabling set
    CHECK(en.ge[0] == doctest::Approx(G).epsilon(1e-12));      // every marking has e >= 0
    double total = 0;
    for (std::size_t k = 0; k <= en.max_degree; ++k) total += en.eq[k] / G;
    CHECK(total == doctest::Approx(1.0).epsilon(1e-12));
    // With unit arcs the degree is the marking, so its law IS the place marginal.
    const std::vector<double> mass = mdd::mdd_rec_marginal(r.mdds, g, 0);
    for (std::size_t k = 0; k <= en.max_degree; ++k)
        CHECK(en.eq[k] == doctest::Approx(mass[k]).epsilon(1e-12));
}

TEST_CASE("MDD-rec runs in exact arithmetic") {
    // The recursion is field operations only -- no square root, no tolerance --
    // so the normalising constant of a rational net is itself exactly rational.
    const std::size_t K = 3;
    std::vector<Rational> mu;
    mu.push_back(Rational(1));
    mu.push_back(Rational(3, 2));
    mu.push_back(Rational(2));
    Matrix<Rational> P(K, K, Rational(0));
    for (std::size_t i = 0; i < K; ++i) P(i, (i + 1) % K) = Rational(1);
    std::vector<double> servers(K, 1.0);
    const mdd::MddDescriptor<Rational> desc = mdd::mdd_descriptor(mu, P, servers, NJOBS);
    const mdd::MDD diagram = mdd::mdd_reachset(desc.domain, desc.init, desc.nextfun);

    std::vector<std::vector<Rational>> g(K);
    for (std::size_t l = 0; l < K; ++l) {
        g[l].assign(static_cast<std::size_t>(desc.domain[l]), Rational(0));
        Rational x(1);
        for (int n = 0; n < desc.domain[l]; ++n) {
            g[l][static_cast<std::size_t>(n)] = x;
            x = x / mu[l];
        }
    }
    const Rational G = mdd::mdd_rec(diagram.to_struct(), g);
    std::vector<std::vector<double>> gd(K);
    for (std::size_t l = 0; l < K; ++l) {
        gd[l].assign(g[l].size(), 0.0);
        for (std::size_t n = 0; n < g[l].size(); ++n)
            gd[l][n] = num_traits<Rational>::to_double(g[l][n]);
    }
    CHECK(num_traits<Rational>::to_double(G) ==
          doctest::Approx(explicit_norm_const(diagram, gd)).epsilon(1e-12));
    CHECK(G > Rational(0));
}
