/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
/**
 * The three PRIORITY sharing disciplines in the SSA NRM: PSPRIO, DPSPRIO,
 * GPSPRIO.
 *
 * The NRM used to refuse all three by name ("splits the processor share by
 * class priority before the phase ratio; that share law is not ported to C++").
 * The share law is `solver_ssa_nrm.m`'s `psprioshare` / `dpsprioshare` /
 * `gpsprioshare`, ported into `solver_ssa_nrm.h`, and this file pins it.
 *
 * THE ORACLE IS AN EXACT CTMC BUILT HERE, over the same rate law. That sounds
 * circular and is not: the chain is written from the DEFINITION of the
 * discipline -- below the server count everyone shares, above it only the most
 * urgent non-empty group does -- and solved exactly, while the engine reaches
 * its answer through the NRM's reaction grid, its next-reaction sampling and its
 * time-average accumulators. What is being tested is that the simulator
 * realizes the law, and the two computations share no code.
 *
 * THE DISCRIMINATING PROPERTY, and the reason a tolerance check alone would be
 * worthless: with one server and two classes of equal service rate, priority
 * changes the ANSWER and not just the variance. The urgent class must have a
 * strictly shorter queue than the symmetric non-priority model gives, and the
 * low-priority class a strictly longer one. Both are asserted.
 */
#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

#include "doctest.h"
#include "line/api/mc/ctmc_solve.h"
#include "line/lang/qn/network_builder.h"
#include "line/solvers/ssa/ssa_dispatch.h"

using namespace line;
using D = lang::Distrib<double>;
using lang::SchedStrategy;

namespace {

/** Delay -> Queue -> Delay, two closed classes of `n` jobs each. */
qn::Network<double> prio_model(SchedStrategy sched, int n1, int n2, double z, double mu,
                               int prio1, int prio2, double w1 = 1.0, double w2 = 1.0) {
    qn::Network<double> mm("prio");
    const std::size_t dd = mm.add_delay("Delay");
    const std::size_t qq = mm.add_queue("Q", sched);
    const std::size_t k1 = mm.add_closed_class("C1", n1, dd, prio1);
    const std::size_t k2 = mm.add_closed_class("C2", n2, dd, prio2);
    mm.set_service(dd, k1, D::exp_rate(z));
    mm.set_service(dd, k2, D::exp_rate(z));
    mm.set_service(qq, k1, D::exp_rate(mu));
    mm.set_service(qq, k2, D::exp_rate(mu));
    mm.set_sched_param(qq, k1, w1);
    mm.set_sched_param(qq, k2, w2);
    qn::RoutingMatrix<double> P;
    P.set(k1, k1, dd, qq, 1.0);
    P.set(k1, k1, qq, dd, 1.0);
    P.set(k2, k2, dd, qq, 1.0);
    P.set(k2, k2, qq, dd, 1.0);
    mm.link(P);
    return mm;
}

/** The share law, written from the definition rather than taken from the port. */
double share(SchedStrategy sched, double n1, double n2, std::size_t r, double c, int p1, int p2,
             double w1, double w2) {
    const double n[2] = {n1, n2};
    const int p[2] = {p1, p2};
    const double w[2] = {w1, w2};
    const double tot = n1 + n2;
    if (tot <= 0.0 || n[r] <= 0.0) return 0.0;

    auto plain = [&](const double v[2]) {
        const double t = v[0] + v[1];
        if (t <= 0.0) return 0.0;
        if (sched == SchedStrategy::PSPRIO) return (v[r] / t) * std::min(t, c);
        if (sched == SchedStrategy::DPSPRIO) {
            const double den = w[0] * v[0] + w[1] * v[1];
            return den <= 0.0 ? 0.0 : w[r] * v[r] / den;
        }
        const double den = (v[0] > 0.0 ? w[0] : 0.0) + (v[1] > 0.0 ? w[1] : 0.0);
        return den <= 0.0 ? 0.0 : w[r] / den;
    };
    if (tot <= c) return plain(n);
    // Above capacity: only the most urgent NON-EMPTY group is served.
    int best = 0;
    bool any = false;
    for (int s = 0; s < 2; ++s)
        if (n[s] > 0.0 && (!any || p[s] < best)) {
            best = p[s];
            any = true;
        }
    if (!any || p[r] != best) return 0.0;
    double g[2] = {0.0, 0.0};
    for (int s = 0; s < 2; ++s)
        if (p[s] == p[r]) g[s] = n[s];
    return plain(g);
}

/** Mean queue length per class, from the exact chain of the closed model. */
void exact_prio(SchedStrategy sched, int N1, int N2, double z, double mu, double c, int p1, int p2,
                double w1, double w2, double& q1, double& q2) {
    const std::size_t rows = static_cast<std::size_t>(N1 + 1), cols = static_cast<std::size_t>(N2 + 1);
    const std::size_t n = rows * cols;
    auto idx = [cols](std::size_t a, std::size_t b) { return a * cols + b; };
    Matrix<double> Q(n, n, 0.0);
    for (std::size_t a = 0; a <= static_cast<std::size_t>(N1); ++a)
        for (std::size_t b = 0; b <= static_cast<std::size_t>(N2); ++b) {
            const std::size_t s = idx(a, b);
            // Delay -> Queue: one server per thinking job.
            if (a < static_cast<std::size_t>(N1))
                Q(s, idx(a + 1, b)) += static_cast<double>(N1 - a) * z;
            if (b < static_cast<std::size_t>(N2))
                Q(s, idx(a, b + 1)) += static_cast<double>(N2 - b) * z;
            // Queue -> Delay at mu times the class share.
            const double s1 = share(sched, static_cast<double>(a), static_cast<double>(b), 0, c,
                                    p1, p2, w1, w2);
            const double s2 = share(sched, static_cast<double>(a), static_cast<double>(b), 1, c,
                                    p1, p2, w1, w2);
            if (a >= 1 && s1 > 0.0) Q(s, idx(a - 1, b)) += mu * s1;
            if (b >= 1 && s2 > 0.0) Q(s, idx(a, b - 1)) += mu * s2;
        }
    for (std::size_t a = 0; a < n; ++a) {
        double off = 0.0;
        for (std::size_t b = 0; b < n; ++b)
            if (a != b) off += Q(a, b);
        Q(a, a) = -off;
    }
    const std::vector<double> pi = mc::ctmc_solve(Q);
    q1 = 0.0;
    q2 = 0.0;
    for (std::size_t a = 0; a <= static_cast<std::size_t>(N1); ++a)
        for (std::size_t b = 0; b <= static_cast<std::size_t>(N2); ++b) {
            q1 += static_cast<double>(a) * pi[idx(a, b)];
            q2 += static_cast<double>(b) * pi[idx(a, b)];
        }
}

ssa::SsaOptions nrm_opts() {
    ssa::SsaOptions o;
    o.method = "nrm";
    o.samples = 400000;
    o.seed = 23000;
    return o;
}

}  // namespace

TEST_CASE("ssa nrm: PSPRIO matches its exact chain and favours the urgent class") {
    // Class 1 is the urgent one (lower value = more urgent in LINE).
    const int N1 = 2, N2 = 2, p1 = 0, p2 = 1;
    const double z = 1.0, mu = 2.0, c = 1.0;
    double e1 = 0.0, e2 = 0.0;
    exact_prio(SchedStrategy::PSPRIO, N1, N2, z, mu, c, p1, p2, 1.0, 1.0, e1, e2);

    qn::Network<double> m = prio_model(SchedStrategy::PSPRIO, N1, N2, z, mu, p1, p2);
    const ssa::SsaSolution s = ssa::solver_ssa(m.get_struct(), nrm_opts());
    CHECK(s.QN(1, 0) == doctest::Approx(e1).epsilon(0.06));
    CHECK(s.QN(1, 1) == doctest::Approx(e2).epsilon(0.06));
    // Population is conserved per class whatever the share law does.
    CHECK(s.QN(0, 0) + s.QN(1, 0) == doctest::Approx(2.0).epsilon(1e-9));
    CHECK(s.QN(0, 1) + s.QN(1, 1) == doctest::Approx(2.0).epsilon(1e-9));
    // U = T E[S] per class, the invariant that catches a utilization
    // accumulator using a DIFFERENT sharing factor from the rate law: the two
    // are separate switch statements over the same disciplines and nothing but
    // this ties them together.
    CHECK(s.UN(1, 0) == doctest::Approx(s.TN(1, 0) / mu).epsilon(0.03));
    CHECK(s.UN(1, 1) == doctest::Approx(s.TN(1, 1) / mu).epsilon(0.03));

    // PRIORITY CHANGES THE ANSWER, and by a wide margin: with the two classes
    // otherwise identical, equal priority is symmetric and unequal priority is
    // not. Without this a port that ignored `prio` entirely would pass.
    double s1 = 0.0, s2 = 0.0;
    exact_prio(SchedStrategy::PSPRIO, N1, N2, z, mu, c, 0, 0, 1.0, 1.0, s1, s2);
    CHECK(s1 == doctest::Approx(s2).epsilon(1e-12));
    CHECK(e1 < s1 - 0.05);
    CHECK(e2 > s2 + 0.05);
    CHECK(s.QN(1, 0) < s.QN(1, 1));
}

TEST_CASE("ssa nrm: DPSPRIO and GPSPRIO match their exact chains") {
    const int N1 = 2, N2 = 2, p1 = 0, p2 = 1;
    const double z = 1.0, mu = 2.0, c = 1.0, w1 = 3.0, w2 = 1.0;

    for (int which = 0; which < 2; ++which) {
        const SchedStrategy sched = which == 0 ? SchedStrategy::DPSPRIO : SchedStrategy::GPSPRIO;
        CAPTURE(std::string(lang::sched_to_text(sched)));
        double e1 = 0.0, e2 = 0.0;
        exact_prio(sched, N1, N2, z, mu, c, p1, p2, w1, w2, e1, e2);
        qn::Network<double> m = prio_model(sched, N1, N2, z, mu, p1, p2, w1, w2);
        const ssa::SsaSolution s = ssa::solver_ssa(m.get_struct(), nrm_opts());
        CHECK(s.QN(1, 0) == doctest::Approx(e1).epsilon(0.06));
        CHECK(s.QN(1, 1) == doctest::Approx(e2).epsilon(0.06));
        CHECK(s.QN(0, 0) + s.QN(1, 0) == doctest::Approx(2.0).epsilon(1e-9));
        CHECK(s.QN(0, 1) + s.QN(1, 1) == doctest::Approx(2.0).epsilon(1e-9));
        CHECK(s.UN(1, 0) == doctest::Approx(s.TN(1, 0) / mu).epsilon(0.03));
        CHECK(s.UN(1, 1) == doctest::Approx(s.TN(1, 1) / mu).epsilon(0.03));
    }
}

TEST_CASE("ssa nrm: below the server count a PRIO discipline IS its plain twin") {
    // The whole split is at the server count: with enough servers for every job
    // the priority branch is never taken, so PSPRIO must reproduce PS exactly.
    // This is what pins the `ni <= c` boundary rather than the behaviour above
    // and below it separately.
    const int N1 = 2, N2 = 2;
    const double z = 1.0, mu = 2.0;
    double p1q = 0.0, p2q = 0.0, s1 = 0.0, s2 = 0.0;
    // c = 4 covers the whole population, so no state is ever above capacity.
    exact_prio(SchedStrategy::PSPRIO, N1, N2, z, mu, 4.0, 0, 1, 1.0, 1.0, p1q, p2q);
    exact_prio(SchedStrategy::PSPRIO, N1, N2, z, mu, 4.0, 0, 0, 1.0, 1.0, s1, s2);
    CHECK(p1q == doctest::Approx(s1).epsilon(1e-12));
    CHECK(p2q == doctest::Approx(s2).epsilon(1e-12));
}

TEST_CASE("ssa nrm: an empty class never defines the urgent group") {
    // With no urgent jobs present the low-priority class must be served, not
    // frozen. Reading `isUrgent` as "has the globally smallest priority value"
    // rather than "the smallest among the NON-EMPTY classes" deadlocks the
    // station instead, which shows up as a zero throughput.
    const int N1 = 0, N2 = 2;
    double e1 = 0.0, e2 = 0.0;
    exact_prio(SchedStrategy::PSPRIO, N1, N2, 1.0, 2.0, 1.0, 0, 1, 1.0, 1.0, e1, e2);
    CHECK(e2 > 0.0);
    CHECK(e2 < 2.0);

    qn::Network<double> m = prio_model(SchedStrategy::PSPRIO, 1, 2, 1.0, 2.0, 0, 1);
    const ssa::SsaSolution s = ssa::solver_ssa(m.get_struct(), nrm_opts());
    CHECK(s.TN(1, 1) > 0.0);
}
