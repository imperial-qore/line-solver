/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * The decision-diagram domain (`mdd`) and the SPN front end (`spn_mdd`).
 *
 * THE ORACLES ARE THE OTHER CODEBASES. `mdd_mcd` is a fixed point of a coupled
 * level system, not a quantity with a closed form, so a structural check ("the
 * marginals sum to one") is satisfied by an iteration with the orientation
 * reversed or an event matrix transposed. The numbers below are MATLAB's and
 * native python's own output, which agree to all twelve digits printed, on the
 * same 3-place cyclic net at N=4 with rates {1, 1.5, 2}:
 *
 *   |S|  = 15 = C(N+K-1, K-1)
 *   QLen = [2.249874392899, 1.069837548149, 0.680288058952]
 *
 * WHAT EACH CASE SEPARATES. The QN and the SPN reach that same answer through
 * different descriptors -- `mdd_descriptor` builds a departure/arrival event
 * pair per routing entry, `spn_mdd` builds one event per (transition, mode)
 * gated by the enabling arcs -- so agreement between them pins the SPN
 * translation and not merely the aggregation they share. The orientation trap is
 * covered by asserting the PER-STATION order of QLen: paper level k is station
 * K-1-k, and reversing it leaves every sum invariant while mislabelling every
 * station.
 */

#include <cmath>
#include <limits>
#include <string>
#include <vector>

#include "doctest.h"
#include "line/api/mdd/mdd.h"
#include "line/api/mdd/mdd_descriptor.h"
#include "line/api/mdd/mdd_mcd.h"
#include "line/api/mdd/mdd_ps.h"
#include "line/api/mdd/mdd_reachset.h"
#include "line/api/spn/spn_mdd.h"
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
/** MATLAB and native-python output on this net, printed at %.12f. */
const double QLEN_REF[3] = {2.249874392899, 1.069837548149, 0.680288058952};

/** Cyclic routing over K stations: station i feeds station i+1 mod K. */
Matrix<double> cyclic_routing(std::size_t K) {
    Matrix<double> P(K, K, 0.0);
    for (std::size_t i = 0; i < K; ++i) P(i, (i + 1) % K) = 1.0;
    return P;
}

/** One mode of a marked graph: consume one from `from`, produce one to `to`. */
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

/** The 3-place cyclic net with `ntokens` tokens, one class. */
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

/** Mean tokens per PLACE level, dropping any phase levels. */
std::vector<double> place_tokens(const spn::SpnResult<double>& r,
                                 const mdd::MddMcdResult<double>& out) {
    std::vector<double> q;
    for (std::size_t l = 0; l < r.info.nplacelevels; ++l) q.push_back(out.QLen[l]);
    return q;
}

}  // namespace

TEST_CASE("the diagram stores a set and indexes it bijectively") {
    // Every 2-tuple over {0,1,2} x {0,1}: cardinality, membership and the
    // lexicographic rank are all checkable by hand.
    std::vector<int> domain;
    domain.push_back(3);
    domain.push_back(2);
    std::vector<std::vector<int>> states;
    for (int a = 0; a < 3; ++a)
        for (int b = 0; b < 2; ++b) {
            std::vector<int> s;
            s.push_back(a);
            s.push_back(b);
            states.push_back(s);
        }
    const mdd::MDD d = mdd::MDD::from_states(domain, states);
    CHECK(d.cardinality() == 6);
    for (std::size_t i = 0; i < states.size(); ++i) {
        CHECK(d.member(states[i]));
        CHECK(d.index(states[i]) == static_cast<long long>(i));  // level 0 most significant
    }
    std::vector<int> absent;
    absent.push_back(2);
    absent.push_back(1);
    mdd::MDD partial(domain);
    partial.insert(states[0]);
    CHECK(!partial.member(absent));
    CHECK(partial.index(absent) == -1);
    CHECK(d.enumerate().size() == 6);
    CHECK(d.stats().num_states == 6);
}

TEST_CASE("the closed cyclic QN reproduces the JAR and python aggregation") {
    const std::size_t K = 3;
    std::vector<double> mu(RATES, RATES + 3);
    std::vector<double> servers(K, 1.0);
    const mdd::MddDescriptor<double> desc =
        mdd::mdd_descriptor(mu, cyclic_routing(K), servers, NJOBS);
    const mdd::MDD diagram = mdd::mdd_reachset(desc.domain, desc.init, desc.nextfun);
    CHECK(diagram.cardinality() == 15);  // C(N+K-1, K-1)

    const mdd::MddMcdResult<double> out = mdd::mdd_mcd(diagram.to_struct(), desc);
    for (std::size_t i = 0; i < K; ++i)
        CHECK(out.QLen[i] == doctest::Approx(QLEN_REF[i]).epsilon(1e-9));
    double total = 0;
    for (std::size_t i = 0; i < K; ++i) total += out.QLen[i];
    CHECK(total == doctest::Approx(static_cast<double>(NJOBS)).epsilon(1e-9));
    // A cyclic single-class net is product form, so every station carries the
    // same throughput and the utilizations are in the ratio of the demands.
    for (std::size_t i = 1; i < K; ++i)
        CHECK(out.X[i] == doctest::Approx(out.X[0]).epsilon(1e-9));
    for (std::size_t i = 0; i < K; ++i)
        CHECK(out.U[i] == doctest::Approx(out.X[i] / RATES[i]).epsilon(1e-9));
}

TEST_CASE("the cyclic SPN reproduces the closed queueing network") {
    // A place with a single-server timed transition IS a station, so the two
    // descriptors must agree to numerical noise on the same rates.
    qn::Network<double> m = cyclic_spn(NJOBS);
    const spn::SpnResult<double> r = spn::spn_mdd(m.get_struct());
    CHECK(r.info.nplacelevels == 3);
    CHECK(r.info.diagram.cardinality() == 15);
    CHECK(r.desc.invariant_value == doctest::Approx(static_cast<double>(NJOBS)).epsilon(1e-12));

    const mdd::MddMcdResult<double> out = mdd::mdd_mcd(r.mdds, r.desc);
    const std::vector<double> q = place_tokens(r, out);
    for (std::size_t i = 0; i < 3; ++i)
        CHECK(q[i] == doctest::Approx(QLEN_REF[i]).epsilon(1e-9));
    // No queueing parameters ride on an SPN descriptor, so X and U are absent
    // rather than filled with a plausible zero.
    CHECK(out.X.empty());
    CHECK(out.U.empty());
    CHECK(r.info.levelname[0] == "P0");
    CHECK(r.info.levelkind[0] == 1);
}

TEST_CASE("the SPN translation refuses what has no Kronecker form") {
    // An IMMEDIATE mode makes vanishing states, which must be eliminated first.
    qn::Network<double> imm("imm");
    const std::size_t p1 = imm.add_place("P1");
    const std::size_t p2 = imm.add_place("P2");
    const std::size_t c = imm.add_closed_class("Tok", 1.0, p1);
    imm.set_service(p1, c, Dist::exp_rate(1.0));
    imm.set_service(p2, c, Dist::exp_rate(1.0));
    qn::TransitionParam<double> tp = one_mode(4, p1, p2, 1.0);
    tp.timing[0] = TimingStrategy::IMMEDIATE;
    imm.add_transition("T1", tp);
    imm.add_transition("T2", one_mode(4, p2, p1, 1.0));
    CHECK_THROWS_AS(spn::spn_mdd(imm.get_struct()), line::UnsupportedError);

    // A marking-dependent firing rate is not a product of per-level terms.
    qn::Network<double> dep("dep");
    const std::size_t q1 = dep.add_place("P1");
    const std::size_t q2 = dep.add_place("P2");
    const std::size_t c2 = dep.add_closed_class("Tok", 1.0, q1);
    dep.set_service(q1, c2, Dist::exp_rate(1.0));
    dep.set_service(q2, c2, Dist::exp_rate(1.0));
    qn::TransitionParam<double> tdep = one_mode(4, q1, q2, 1.0);
    tdep.firingdep.push_back([](const std::vector<double>& mk) { return mk.empty() ? 1.0 : mk[0]; });
    dep.add_transition("T1", tdep);
    dep.add_transition("T2", one_mode(4, q2, q1, 1.0));
    CHECK_THROWS_AS(spn::spn_mdd(dep.get_struct()), line::UnsupportedError);
}

TEST_CASE("a multiclass net is refused rather than class-aggregated") {
    // A NetworkStruct holds the arcs per (mode, node) with no class dimension,
    // so per-(place,class) levels cannot be recovered from it. Aggregating would
    // solve a DIFFERENT net, so the translation refuses by name.
    qn::Network<double> m("multiclass");
    const std::size_t p1 = m.add_place("P1");
    const std::size_t p2 = m.add_place("P2");
    const std::size_t c1 = m.add_closed_class("A", 1.0, p1);
    const std::size_t c2 = m.add_closed_class("B", 1.0, p1);
    m.set_service(p1, c1, Dist::exp_rate(1.0));
    m.set_service(p2, c1, Dist::exp_rate(1.0));
    m.set_service(p1, c2, Dist::exp_rate(1.0));
    m.set_service(p2, c2, Dist::exp_rate(1.0));
    m.add_transition("T1", one_mode(4, p1, p2, 1.0));
    m.add_transition("T2", one_mode(4, p2, p1, 1.0));
    CHECK_THROWS_AS(spn::spn_mdd(m.get_struct()), line::UnsupportedError);
}

TEST_CASE("mdd_ps counts jobs per phase and refuses a multi-server station") {
    // Exponential service: the per-phase-count encoding collapses to the plain
    // count, so mdd_ps must reproduce mdd_descriptor's domain and answer at PS.
    const std::size_t K = 3;
    std::vector<double> mu(RATES, RATES + 3);
    std::vector<double> servers(K, 1.0);
    const mdd::MddDescriptor<double> desc = mdd::mdd_ps(mu, cyclic_routing(K), servers, NJOBS);
    for (std::size_t i = 0; i < K; ++i) CHECK(desc.domain[i] == NJOBS + 1);
    const mdd::MDD diagram = mdd::mdd_reachset(desc.domain, desc.init, desc.nextfun);
    CHECK(diagram.cardinality() == 15);
    const mdd::MddMcdResult<double> out = mdd::mdd_mcd(diagram.to_struct(), desc);
    for (std::size_t i = 0; i < K; ++i)
        CHECK(out.QLen[i] == doctest::Approx(QLEN_REF[i]).epsilon(1e-9));

    std::vector<double> multi(K, 1.0);
    multi[1] = 2.0;
    CHECK_THROWS_AS(mdd::mdd_ps(mu, cyclic_routing(K), multi, NJOBS), line::InputError);
}

TEST_CASE("the level aggregation runs in exact arithmetic") {
    // Rational has no square root, so the Householder QR of the double path
    // cannot serve it; the exact backend solves the same appended-normalisation
    // system through line::lstsq, where there is no condition number to square.
    // The answer must be the double one to double precision, and exactly
    // conservative: the level marginals sum to the population with NO residual.
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
    CHECK(diagram.cardinality() == 15);
    const mdd::MddMcdResult<Rational> out = mdd::mdd_mcd(diagram.to_struct(), desc);
    for (std::size_t i = 0; i < K; ++i)
        CHECK(num_traits<Rational>::to_double(out.QLen[i]) ==
              doctest::Approx(QLEN_REF[i]).epsilon(1e-9));
    Rational total(0);
    for (std::size_t i = 0; i < K; ++i) total += out.QLen[i];
    CHECK(total == Rational(NJOBS));  // exact, not to a tolerance
}
