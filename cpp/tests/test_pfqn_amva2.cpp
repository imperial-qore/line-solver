/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Approximate-MVA drivers: pfqn_ab_amva, pfqn_schmidt_ext, pfqn_rd, and the
 * shortest-job-next pair pfqn_mvasjn / pfqn_amvasjn (see the second banner
 * below; both are approximate, the station equation being Kant's whether the
 * population recursion is exact or closed).
 *
 * Oracles:
 *  - MATLAB reference values, produced by running the reference implementations
 *    on the models below and printing at 15 significant digits;
 *  - the operational laws that any consistent solution must obey regardless of
 *    the approximation: population conservation sum_i Q(i,r) = N_r, Little's
 *    law Q(i,r) = X_r W(i,r) (with unit visit ratios), and the utilization law
 *    U(i,r) = X_r D(i,r)/c_i;
 *  - the exact solution from pfqn_mva / pfqn_ncld, on the models where the
 *    approximation is known to be exact (single-server product form) and, where
 *    it is not, at the documented accuracy of the approximation.
 */
#include <cmath>
#include <vector>

#include "doctest.h"
#include "line/api/pfqn/pfqn_ab_amva.h"
#include "line/api/pfqn/pfqn_ca.h"
#include "line/api/pfqn/pfqn_mva.h"
#include "line/api/pfqn/pfqn_ncld.h"
#include "line/api/pfqn/pfqn_rd.h"
#include "line/api/pfqn/pfqn_schmidt.h"
#include "line/api/pfqn/pfqn_schmidt_ext.h"
#include "line/api/pfqn/pfqn_sjn.h"
#include "line/api/pfqn/pfqn_tay.h"

using line::Matrix;
using line::Rational;
using namespace line::pfqn;

namespace {

/** Two stations, two classes: D = [0.5 0.3; 0.4 0.6]. */
Matrix<double> amvaD() {
    Matrix<double> D(2, 2);
    D(0, 0) = 0.5; D(0, 1) = 0.3;
    D(1, 0) = 0.4; D(1, 1) = 0.6;
    return D;
}

Matrix<double> unitVisits(std::size_t M, std::size_t R) { return Matrix<double>(M, R, 1.0); }

}  // namespace

// ---------------------------------------------------------------------------
// pfqn_ab_amva
// ---------------------------------------------------------------------------

TEST_CASE("Akyildiz-Bolch AMVA matches MATLAB on a PS plus two-server FCFS model") {
    const std::vector<int> N{4, 3};
    const std::vector<int> ns{1, 2};
    const std::vector<SchedStrategy> sched{SchedStrategy::PS, SchedStrategy::FCFS};
    const auto r = pfqn_ab_amva(amvaD(), N, unitVisits(2, 2), ns, sched, false,
                                AbMarginalMethod::Ab);
    // MATLAB pfqn_ab_amva(D,[4 3],ones(2,2),[1;2],[PS;FCFS],false,'ab')
    CHECK(r.QN(0, 0) == doctest::Approx(3.39745553070038).epsilon(1e-11));
    CHECK(r.QN(0, 1) == doctest::Approx(2.19420329544517).epsilon(1e-11));
    CHECK(r.QN(1, 0) == doctest::Approx(0.602544469299617).epsilon(1e-11));
    CHECK(r.QN(1, 1) == doctest::Approx(0.805796704554835).epsilon(1e-11));
    CHECK(r.UN(0, 0) == doctest::Approx(0.612977087993678).epsilon(1e-11));
    CHECK(r.UN(0, 1) == doctest::Approx(0.375).epsilon(1e-11));
    CHECK(r.UN(1, 0) == doctest::Approx(0.245190835197471).epsilon(1e-11));
    CHECK(r.UN(1, 1) == doctest::Approx(0.375).epsilon(1e-11));
    CHECK(r.RN(0, 0) == doctest::Approx(2.77127448745313).epsilon(1e-11));
    CHECK(r.RN(1, 1) == doctest::Approx(0.644637363643868).epsilon(1e-11));
    CHECK(r.CN[0] == doctest::Approx(3.26276469247188).epsilon(1e-11));
    CHECK(r.CN[1] == doctest::Approx(2.4).epsilon(1e-11));
    CHECK(r.XN[0] == doctest::Approx(1.22595417598736).epsilon(1e-11));
    CHECK(r.XN[1] == doctest::Approx(1.25).epsilon(1e-11));
    CHECK(r.totiter == 8);
}

TEST_CASE("Akyildiz-Bolch AMVA obeys the operational laws") {
    const std::vector<int> N{4, 3};
    const std::vector<int> ns{1, 2};
    const std::vector<SchedStrategy> sched{SchedStrategy::PS, SchedStrategy::FCFS};
    const Matrix<double> D = amvaD();
    const auto r = pfqn_ab_amva(D, N, unitVisits(2, 2), ns, sched, false, AbMarginalMethod::Ab);

    for (std::size_t c = 0; c < 2; ++c) {
        // Population conservation.
        CHECK(r.QN(0, c) + r.QN(1, c) == doctest::Approx(N[c]).epsilon(1e-10));
        // Little's law at each station, with unit visit ratios so that the
        // station-1 throughput the reference reports IS the system throughput.
        for (std::size_t i = 0; i < 2; ++i) {
            CHECK(r.QN(i, c) == doctest::Approx(r.XN[c] * r.RN(i, c)).epsilon(1e-10));
            // Utilization law.
            CHECK(r.UN(i, c) ==
                  doctest::Approx(r.XN[c] * D(i, c) / static_cast<double>(ns[i])).epsilon(1e-10));
        }
        // Cycle time is the sum of the residence times.
        CHECK(r.CN[c] == doctest::Approx(r.RN(0, c) + r.RN(1, c)).epsilon(1e-10));
    }
}

TEST_CASE("Akyildiz-Bolch AMVA matches MATLAB with a delay station") {
    const std::vector<int> N{2, 2};
    const std::vector<int> ns{1, 1};
    const std::vector<SchedStrategy> sched{SchedStrategy::INF, SchedStrategy::PS};
    const auto r = pfqn_ab_amva(amvaD(), N, unitVisits(2, 2), ns, sched, false,
                                AbMarginalMethod::Ab);
    CHECK(r.QN(0, 0) == doctest::Approx(0.546326121979193).epsilon(1e-11));
    CHECK(r.QN(0, 1) == doctest::Approx(0.275451656971627).epsilon(1e-11));
    CHECK(r.QN(1, 0) == doctest::Approx(1.45367387802081).epsilon(1e-11));
    CHECK(r.QN(1, 1) == doctest::Approx(1.72454834302837).epsilon(1e-11));
    CHECK(r.XN[0] == doctest::Approx(1.09265224395839).epsilon(1e-11));
    CHECK(r.XN[1] == doctest::Approx(0.918172189905423).epsilon(1e-11));
    CHECK(r.totiter == 3);

    // The same model solved exactly, with station 1 treated as the delay: the
    // approximation is within 0.7 percent on both throughputs.
    Matrix<double> Lq(1, 2);
    Lq(0, 0) = 0.4;
    Lq(0, 1) = 0.6;
    Matrix<double> Z(1, 2);
    Z(0, 0) = 0.5;
    Z(0, 1) = 0.3;
    const auto exact = pfqn_mva(Lq, N, Z);
    CHECK(std::fabs(r.XN[0] - exact.XN[0]) / exact.XN[0] < 0.01);
    CHECK(std::fabs(r.XN[1] - exact.XN[1]) / exact.XN[1] < 0.01);
}

TEST_CASE("the scatter marginal rule is a distinct approximation, and matches MATLAB") {
    const std::vector<int> N{4, 3};
    const std::vector<int> ns{1, 2};
    const std::vector<SchedStrategy> sched{SchedStrategy::PS, SchedStrategy::FCFS};
    const auto r = pfqn_ab_amva(amvaD(), N, unitVisits(2, 2), ns, sched, false,
                                AbMarginalMethod::Scat);
    CHECK(r.QN(0, 0) == doctest::Approx(3.23946456835715).epsilon(1e-11));
    CHECK(r.QN(0, 1) == doctest::Approx(1.92234081274822).epsilon(1e-11));
    CHECK(r.QN(1, 0) == doctest::Approx(0.760535431642846).epsilon(1e-11));
    CHECK(r.QN(1, 1) == doctest::Approx(1.07765918725178).epsilon(1e-11));
    CHECK(r.XN[0] == doctest::Approx(1.26755905273808).epsilon(1e-11));
    CHECK(r.XN[1] == doctest::Approx(1.19739909694642).epsilon(1e-11));
    CHECK(r.totiter == 4);
}

TEST_CASE("the Schmidt FCFS wait is a distinct branch, and matches MATLAB") {
    const std::vector<int> N{2, 2};
    const std::vector<int> ns{1, 2};
    const std::vector<SchedStrategy> sched{SchedStrategy::PS, SchedStrategy::FCFS};
    const auto r = pfqn_ab_amva(amvaD(), N, unitVisits(2, 2), ns, sched, true,
                                AbMarginalMethod::Ab);
    CHECK(r.QN(0, 0) == doctest::Approx(1.54241099460735).epsilon(1e-11));
    CHECK(r.QN(0, 1) == doctest::Approx(1.22332409549041).epsilon(1e-11));
    CHECK(r.QN(1, 0) == doctest::Approx(0.457589005392651).epsilon(1e-11));
    CHECK(r.QN(1, 1) == doctest::Approx(0.77667590450959).epsilon(1e-11));
    CHECK(r.RN(1, 0) == doctest::Approx(0.429505934794476).epsilon(1e-11));
    CHECK(r.XN[0] == doctest::Approx(1.06538459267533).epsilon(1e-11));
    CHECK(r.XN[1] == doctest::Approx(1.26989131854211).epsilon(1e-11));
    CHECK(r.totiter == 3);
}

TEST_CASE("Akyildiz-Bolch AMVA approaches exact MVA on a single-class single-server model") {
    Matrix<double> D(2, 1);
    D(0, 0) = 0.5;
    D(1, 0) = 0.4;
    const std::vector<int> N{4};
    const std::vector<int> ns{1, 1};
    const std::vector<SchedStrategy> sched{SchedStrategy::PS, SchedStrategy::PS};
    const auto r = pfqn_ab_amva(D, N, unitVisits(2, 1), ns, sched, false, AbMarginalMethod::Ab);
    // MATLAB pfqn_ab_amva on the same model.
    CHECK(r.QN(0, 0) == doctest::Approx(2.40034830783106).epsilon(1e-11));
    CHECK(r.QN(1, 0) == doctest::Approx(1.59965169216894).epsilon(1e-11));
    CHECK(r.XN[0] == doctest::Approx(1.75996516921689).epsilon(1e-11));
    CHECK(r.totiter == 9);

    // Exact MVA: X = 1.75630652070443, Q = [2.43693479295574; 1.56306520704426].
    const auto exact = pfqn_mva(D, N);
    CHECK(exact.XN[0] == doctest::Approx(1.75630652070443).epsilon(1e-12));
    CHECK(std::fabs(r.XN[0] - exact.XN[0]) / exact.XN[0] < 3e-3);
    CHECK(std::fabs(r.QN(0, 0) - exact.QN(0, 0)) / exact.QN(0, 0) < 2e-2);
}

TEST_CASE("Akyildiz-Bolch AMVA rejects malformed input") {
    const std::vector<SchedStrategy> sched{SchedStrategy::PS, SchedStrategy::PS};
    CHECK_THROWS_AS(pfqn_ab_amva(amvaD(), std::vector<int>{4}, unitVisits(2, 2),
                                 std::vector<int>{1, 1}, sched),
                    line::InputError);
    CHECK_THROWS_AS(pfqn_ab_amva(amvaD(), std::vector<int>{4, 3}, unitVisits(2, 2),
                                 std::vector<int>{1, 0}, sched),
                    line::InputError);
}

// ---------------------------------------------------------------------------
// pfqn_schmidt_ext
// ---------------------------------------------------------------------------

TEST_CASE("extended Schmidt matches MATLAB on a class-dependent two-server FCFS station") {
    const std::vector<int> N{2, 2};
    Matrix<int> S(2, 1);
    S(0, 0) = 2;
    S(1, 0) = 1;
    const std::vector<SchedStrategy> sched{SchedStrategy::FCFS, SchedStrategy::PS};
    const auto r = pfqn_schmidt_ext(amvaD(), N, S, sched);
    // MATLAB pfqn_schmidt_ext(D,[2 2],[2;1],[FCFS;PS]). MATLAB replicates XN
    // over the stations before returning it; the port returns the single row.
    CHECK(r.XN[0] == doctest::Approx(1.09555761590953).epsilon(1e-11));
    CHECK(r.XN[1] == doctest::Approx(0.924330274133215).epsilon(1e-11));
    CHECK(r.QN(0, 0) == doctest::Approx(0.547778807954763).epsilon(1e-11));
    CHECK(r.QN(0, 1) == doctest::Approx(0.278859182028621).epsilon(1e-11));
    CHECK(r.QN(1, 0) == doctest::Approx(1.45222119204524).epsilon(1e-11));
    CHECK(r.QN(1, 1) == doctest::Approx(1.72114081797138).epsilon(1e-11));
    CHECK(r.UN(0, 0) == doctest::Approx(0.273889403977382).epsilon(1e-11));
    CHECK(r.UN(1, 1) == doctest::Approx(0.554598164479929).epsilon(1e-11));
    CHECK(r.CN(0, 0) == doctest::Approx(0.5).epsilon(1e-11));
    CHECK(r.CN(0, 1) == doctest::Approx(0.301687816392382).epsilon(1e-11));
    CHECK(r.CN(1, 0) == doctest::Approx(1.3255543761061).epsilon(1e-11));
    CHECK(r.CN(1, 1) == doctest::Approx(1.86204094590039).epsilon(1e-11));

    // The alpha correction genuinely changes the answer: plain Schmidt on the
    // same model returns X = [1.08099859066699, 0.908860778745597].
    const auto plain = pfqn_schmidt(amvaD(), N, S, sched);
    CHECK(plain.XN[0] == doctest::Approx(1.08099859066699).epsilon(1e-11));
    CHECK(r.XN[0] != doctest::Approx(plain.XN[0]).epsilon(1e-6));
}

TEST_CASE("extended Schmidt obeys the operational laws") {
    const std::vector<int> N{2, 2};
    Matrix<int> S(2, 1);
    S(0, 0) = 2;
    S(1, 0) = 1;
    const std::vector<SchedStrategy> sched{SchedStrategy::FCFS, SchedStrategy::PS};
    const Matrix<double> D = amvaD();
    const auto r = pfqn_schmidt_ext(D, N, S, sched);
    for (std::size_t c = 0; c < 2; ++c) {
        CHECK(r.QN(0, c) + r.QN(1, c) == doctest::Approx(N[c]).epsilon(1e-10));
        for (std::size_t i = 0; i < 2; ++i) {
            CHECK(r.QN(i, c) == doctest::Approx(r.XN[c] * r.CN(i, c)).epsilon(1e-10));
            CHECK(r.UN(i, c) ==
                  doctest::Approx(r.XN[c] * D(i, c) / static_cast<double>(S(i, 0))).epsilon(1e-10));
        }
    }
}

TEST_CASE("extended Schmidt is exact on a single-server product-form model") {
    // Class-dependent demands but single servers everywhere: the FCFS branch
    // reduces to the arrival theorem and the recursion is exact MVA.
    Matrix<double> D(2, 2);
    D(0, 0) = 0.5; D(0, 1) = 0.5;
    D(1, 0) = 0.4; D(1, 1) = 0.6;
    const std::vector<int> N{2, 1};
    Matrix<int> S(2, 1, 1);
    const std::vector<SchedStrategy> sched{SchedStrategy::FCFS, SchedStrategy::PS};
    const auto r = pfqn_schmidt_ext(D, N, S, sched);
    const auto exact = pfqn_mva(D, N);
    CHECK(r.XN[0] == doctest::Approx(exact.XN[0]).epsilon(1e-12));
    CHECK(r.XN[1] == doctest::Approx(exact.XN[1]).epsilon(1e-12));
    CHECK(r.QN(0, 0) == doctest::Approx(exact.QN(0, 0)).epsilon(1e-12));
    CHECK(r.QN(1, 1) == doctest::Approx(exact.QN(1, 1)).epsilon(1e-12));
    // MATLAB pfqn_schmidt_ext on the same model.
    CHECK(r.XN[0] == doctest::Approx(1.11027756939235).epsilon(1e-11));
    CHECK(r.XN[1] == doctest::Approx(0.4576144036009).epsilon(1e-11));
    CHECK(r.QN(0, 0) == doctest::Approx(1.11777944486122).epsilon(1e-11));
    CHECK(r.QN(0, 1) == doctest::Approx(0.491372843210803).epsilon(1e-11));
}

TEST_CASE("extended Schmidt is exact when every station is class independent") {
    Matrix<double> D(2, 2);
    D(0, 0) = 0.5; D(0, 1) = 0.5;
    D(1, 0) = 0.4; D(1, 1) = 0.4;
    const std::vector<int> N{2, 1};
    Matrix<int> S(2, 1, 1);
    const std::vector<SchedStrategy> sched{SchedStrategy::FCFS, SchedStrategy::FCFS};
    const auto r = pfqn_schmidt_ext(D, N, S, sched);
    const auto exact = pfqn_mva(D, N);
    // MATLAB: X = [1.10207768744354, 0.55103884372177].
    CHECK(r.XN[0] == doctest::Approx(1.10207768744354).epsilon(1e-11));
    CHECK(r.XN[1] == doctest::Approx(0.55103884372177).epsilon(1e-11));
    CHECK(r.QN(0, 0) == doctest::Approx(1.18337850045167).epsilon(1e-11));
    CHECK(r.QN(1, 1) == doctest::Approx(0.408310749774164).epsilon(1e-11));
    CHECK(r.XN[0] == doctest::Approx(exact.XN[0]).epsilon(1e-12));
    CHECK(r.QN(0, 0) == doctest::Approx(exact.QN(0, 0)).epsilon(1e-12));
}

TEST_CASE("extended Schmidt says so when the alpha solve cannot carry per-class servers") {
    const std::vector<int> N{2, 2};
    Matrix<int> S(2, 2);
    S(0, 0) = 2; S(0, 1) = 2;
    S(1, 0) = 1; S(1, 1) = 1;
    const std::vector<SchedStrategy> sched{SchedStrategy::FCFS, SchedStrategy::PS};
    CHECK_THROWS_AS(pfqn_schmidt_ext(amvaD(), N, S, sched), line::InputError);
}

// ---------------------------------------------------------------------------
// pfqn_rd
// ---------------------------------------------------------------------------

TEST_CASE("the reduction heuristic matches MATLAB and the exact constant on a two-server model") {
    Matrix<double> L(2, 2);
    L(0, 0) = 1.0; L(0, 1) = 0.5;
    L(1, 0) = 0.3; L(1, 1) = 0.9;
    Matrix<double> mu(2, 3);
    mu(0, 0) = 1; mu(0, 1) = 2; mu(0, 2) = 2;
    mu(1, 0) = 1; mu(1, 1) = 1; mu(1, 2) = 1;
    const std::vector<int> N{2, 1};
    const auto r = pfqn_rd(L, N, Matrix<double>(), mu);
    // MATLAB pfqn_rd: lGN = 0.58945194422118, Cgamma = 1.7624633431085.
    CHECK(r.lGN == doctest::Approx(0.58945194422118).epsilon(1e-12));
    CHECK(r.Cgamma == doctest::Approx(1.7624633431085).epsilon(1e-12));
    // The exact load-dependent constant. MATLAB pfqn_ncld reports
    // 0.589451944690575 and this tree's pfqn_ncld 0.5894519442211803, a
    // relative difference of 8e-10 that belongs to the two convolutions, not to
    // the heuristic; the heuristic agrees with the tree's own exact route to
    // the last bit, which is the statement being made here.
    // 'exact' is NAMED rather than left to `default`, which clears the
    // Choudhury-Leung-Whitt cost gate on a model this small and answers by
    // contour inversion to ~1e-9 instead of by the exact ladder.
    const auto exact = pfqn_ncld(L, N, Matrix<double>(), mu, NcldMethod::Exact, 0.0, NcOptions());
    CHECK(exact.method == "exact/gld");
    CHECK(exact.lG == doctest::Approx(0.589451944690575).epsilon(1e-8));
    CHECK(std::fabs(r.lGN - exact.lG) < 1e-12);
}

TEST_CASE("the reduction heuristic reproduces MATLAB's correction factor on a three-server model") {
    Matrix<double> L(2, 2);
    L(0, 0) = 0.7; L(0, 1) = 0.2;
    L(1, 0) = 0.4; L(1, 1) = 0.4;
    Matrix<double> mu(2, 4);
    mu(0, 0) = 1; mu(0, 1) = 2; mu(0, 2) = 3; mu(0, 3) = 3;
    for (std::size_t k = 0; k < 4; ++k) mu(1, k) = 1;
    const std::vector<int> N{2, 2};
    const Matrix<double> Z = Matrix<double>::row({1.0, 1.0});
    const auto r = pfqn_rd(L, N, Z, mu);
    // MATLAB pfqn_rd: Cgamma = 1.84790669834143, lGN = 1.46788225760505. The
    // correction factor matches to the last digit; lGN differs by 3.2e-6
    // because the reference's method='default' dispatches its load-independent
    // constant to the cub/le family while the port uses the exact convolution
    // (that difference is visible on RD3 below, where the heuristic is inert).
    CHECK(r.Cgamma == doctest::Approx(1.84790669834143).epsilon(1e-12));
    CHECK(r.lGN == doctest::Approx(1.46788225760505).epsilon(1e-5));
    // The heuristic overshoots the exact constant (MATLAB pfqn_ncld
    // 1.36986852641263) by about 7 percent in the log on this model.
    const auto exact = pfqn_ncld(L, N, Z, mu);
    CHECK(exact.lG == doctest::Approx(1.36986852641263).epsilon(1e-8));
    CHECK(std::fabs(r.lGN - exact.lG) / std::fabs(exact.lG) < 0.08);
}

TEST_CASE("an all-load-independent model reduces the heuristic to the plain constant") {
    Matrix<double> L(2, 2);
    L(0, 0) = 0.5; L(0, 1) = 0.3;
    L(1, 0) = 0.4; L(1, 1) = 0.6;
    const Matrix<double> mu(2, 3, 1.0);
    const std::vector<int> N{2, 1};
    const Matrix<double> Z = Matrix<double>::row({1.0, 0.0});
    const auto r = pfqn_rd(L, N, Z, mu);
    // Every beta is infinite, so the residual is empty and the answer is the
    // load-independent constant unchanged. MATLAB reaches the same branch but
    // RETURNS WITHOUT ASSIGNING Cgamma, so [lGN,Cgamma] = pfqn_rd(...) raises
    // "Output argument Cgamma not assigned"; the port sets it to one.
    CHECK(r.Cgamma == 1.0);
    const auto ca = pfqn_ca(L, N, Z);
    CHECK(ca.lG == doctest::Approx(1.00099945980111).epsilon(1e-12));
    CHECK(r.lGN == doctest::Approx(ca.lG).epsilon(1e-14));
}

TEST_CASE("the reduction heuristic handles the degenerate populations") {
    const Matrix<double> L(2, 2, 0.5);
    const Matrix<double> mu(2, 1, 1.0);
    const auto zero = pfqn_rd(L, std::vector<int>{0, 0}, Matrix<double>(), mu);
    CHECK(zero.lGN == 0.0);
    CHECK(zero.Cgamma == 1.0);
    const auto neg = pfqn_rd(L, std::vector<int>{-1, 0}, Matrix<double>(), mu);
    CHECK(neg.lGN == -std::numeric_limits<double>::infinity());
}

// ---------------------------------------------------------------------------
// arithmetic classification
// ---------------------------------------------------------------------------

TEST_CASE("plain Schmidt stays exact-capable while its extension does not") {
    // pfqn_schmidt is a finite rational recursion and instantiates at Rational;
    // pfqn_schmidt_ext, pfqn_ab_amva and pfqn_rd are gated on
    // has_transcendental and would not compile there, which is the point of the
    // gate. This case pins the boundary from the side that must keep working.
    Matrix<Rational> D(2, 2);
    D(0, 0) = Rational(1, 2); D(0, 1) = Rational(3, 10);
    D(1, 0) = Rational(2, 5); D(1, 1) = Rational(3, 5);
    Matrix<int> S(2, 1, 1);
    const std::vector<SchedStrategy> sched{SchedStrategy::FCFS, SchedStrategy::PS};
    // Every station is single server, so Schmidt's recursion IS exact MVA and
    // the two must agree as RATIONALS, not merely to rounding. Both return
    // X = [400/357, 610/1071] at N = [2,1] and 3/4 for both classes at N = [1,1].
    const std::vector<int> N11{1, 1}, N21{2, 1};
    const auto r = pfqn_schmidt(D, N11, S, sched);
    CHECK(r.XN[0] == Rational(3, 4));
    CHECK(r.XN[1] == Rational(3, 4));
    CHECK(r.QN(0, 0) == Rational(1, 2));
    const auto r2 = pfqn_schmidt(D, N21, S, sched);
    const auto e2 = pfqn_mva(D, N21);
    CHECK(r2.XN[0] == Rational(400, 357));
    CHECK(r2.XN[1] == Rational(610, 1071));
    CHECK(r2.XN[0] == e2.XN[0]);
    CHECK(r2.XN[1] == e2.XN[1]);
}

// ---------------------------------------------------------------------------
// pfqn_mvasjn / pfqn_amvasjn: shortest-job-next stations (Kant 1992)
//
// The api layer under its OWN oracles, independently of the SolverMVA analyzer
// that consumes it. Every reference number below comes from calling MATLAB's
// pfqn_mvasjn / pfqn_amvasjn directly on the same arguments, printed at 15
// significant digits.
//
// The MATLAB numbers pin the port; the operational laws asserted alongside them
// pin the answer itself, and they are what catches a port that is
// self-consistently wrong. Two hold at every SJN station whatever the
// discipline does:
//
//   population   sum_m Q(m,r) + X_r Z_r = N_r     (Q covers only the queueing
//                                                  stations, Z holds the rest)
//   utilization  U(m,r) = X_r L(m,r)
//   Little       Q(m,r) = X_r C(m,r)
//
// The population law is the sharp one under the utilization cap: sjn_cap acts
// on the residence time and never on the throughput precisely so that it keeps
// holding, and a port that capped X directly would lose jobs and fail it while
// still reporting a plausible U.
// ---------------------------------------------------------------------------

namespace {

/** Row vector as the 1 x n / n x 1 Matrix the SJN api takes. */
Matrix<double> sjnMat(std::size_t M, std::size_t R, std::initializer_list<double> vals) {
    Matrix<double> m(M, R);
    std::size_t k = 0;
    for (double v : vals) {
        m(k / R, k % R) = v;
        ++k;
    }
    return m;
}

/** The three operational laws, asserted on whatever the recursion returned. */
void sjnLaws(const SjnResult& r, const Matrix<double>& L, const std::vector<double>& N,
             const std::vector<double>& Z) {
    const std::size_t M = L.rows(), R = L.cols();
    for (std::size_t s = 0; s < R; ++s) {
        double tot = 0.0;
        for (std::size_t m = 0; m < M; ++m) {
            tot += r.QN(m, s);
            CHECK(r.UN(m, s) == doctest::Approx(r.XN[s] * L(m, s)).epsilon(1e-12));
            CHECK(r.QN(m, s) == doctest::Approx(r.XN[s] * r.CN(m, s)).epsilon(1e-12));
        }
        CHECK(tot + r.XN[s] * Z[s] == doctest::Approx(N[s]).epsilon(1e-10));
    }
}

}  // namespace

TEST_CASE("pfqn_mvasjn solves the SJN population lattice") {
    SUBCASE("one SJN station, one class, exponential") {
        const Matrix<double> L = sjnMat(1, 1, {0.5}), V = sjnMat(1, 1, {1.0});
        const Matrix<double> scv = sjnMat(1, 1, {1.0});
        const std::vector<double> N{3.0}, Z{1.0};
        const SjnResult r = pfqn_mvasjn(L, N, Z, scv, {0}, V, SjnOptions());
        CHECK(r.XN[0] == doctest::Approx(1.57374156251021).epsilon(1e-11));
        CHECK(r.QN(0, 0) == doctest::Approx(1.42625843748979).epsilon(1e-11));
        CHECK(r.UN(0, 0) == doctest::Approx(0.786870781255104).epsilon(1e-11));
        CHECK(r.CN(0, 0) == doctest::Approx(0.906285041627057).epsilon(1e-11));
        CHECK(r.iter == 1);  // the lattice is not an iteration
        CHECK_FALSE(r.capped);
        sjnLaws(r, L, N, Z);
    }
    SUBCASE("an ordinary queueing station alongside the SJN one") {
        // Row 2 takes the plain single-server MVA equation, which is the other
        // half of the recursion and is silent when wrong.
        const Matrix<double> L = sjnMat(2, 1, {1.0 / 3.0, 0.5}), V = sjnMat(2, 1, {1.0, 1.0});
        const Matrix<double> scv = sjnMat(2, 1, {1.0, 1.0});
        const std::vector<double> N{4.0}, Z{1.0};
        const SjnResult r = pfqn_mvasjn(L, N, Z, scv, {0}, V, SjnOptions());
        CHECK(r.XN[0] == doctest::Approx(1.5646100304345).epsilon(1e-11));
        CHECK(r.QN(0, 0) == doctest::Approx(0.813639699460493).epsilon(1e-11));
        CHECK(r.QN(1, 0) == doctest::Approx(1.62175027010501).epsilon(1e-11));
        CHECK(r.CN(0, 0) == doctest::Approx(0.520027152858366).epsilon(1e-11));
        CHECK(r.CN(1, 0) == doctest::Approx(1.0365204354817).epsilon(1e-11));
        sjnLaws(r, L, N, Z);
    }
    SUBCASE("two SJN stations at once") {
        // Each carries its own grid and its own profile; sharing one by mistake
        // is invisible until the two stations differ in demand, as they do here.
        const Matrix<double> L = sjnMat(2, 1, {0.4, 0.3}), V = sjnMat(2, 1, {1.0, 1.0});
        const Matrix<double> scv = sjnMat(2, 1, {1.0, 1.0});
        const std::vector<double> N{3.0}, Z{1.0};
        const SjnResult r = pfqn_mvasjn(L, N, Z, scv, {0, 1}, V, SjnOptions());
        CHECK(r.XN[0] == doctest::Approx(1.48282602418192).epsilon(1e-11));
        CHECK(r.QN(0, 0) == doctest::Approx(0.906297266160705).epsilon(1e-11));
        CHECK(r.QN(1, 0) == doctest::Approx(0.610876709657371).epsilon(1e-11));
        CHECK(r.WX.size() == 2);
        CHECK(r.WX[0].station == 1);
        CHECK(r.WX[1].station == 2);
        sjnLaws(r, L, N, Z);
    }
    SUBCASE("a non-unit visit ratio, so that S = L/V is exercised") {
        // The job size the discipline compares is one VISIT's service time, not
        // the demand accumulated over all visits. A port that fed L rather than
        // L/V to the size distribution still satisfies every law above.
        const Matrix<double> L = sjnMat(1, 1, {1.0}), V = sjnMat(1, 1, {2.0});
        const Matrix<double> scv = sjnMat(1, 1, {1.0});
        const std::vector<double> N{3.0}, Z{1.0};
        const SjnResult r = pfqn_mvasjn(L, N, Z, scv, {0}, V, SjnOptions());
        CHECK(r.XN[0] == doctest::Approx(0.9142077494212).epsilon(1e-11));
        CHECK(r.QN(0, 0) == doctest::Approx(2.0857922505788).epsilon(1e-11));
        CHECK(r.CN(0, 0) == doctest::Approx(2.28152982940623).epsilon(1e-11));
        sjnLaws(r, L, N, Z);
    }
}

TEST_CASE("the SJN size distribution is reconstructed from two moments") {
    // The density is not an input: sjn_fit rebuilds it from the mean and the
    // SCV, taking a branching Erlang below CV^2 = 1 and a balanced-means
    // hyperexponential above it. Both arms are exercised here on the SAME
    // instance, and the SCV must move the answer in the right direction: a
    // more variable job size lengthens the queue under SJN as under any
    // work-conserving discipline.
    const Matrix<double> L = sjnMat(1, 1, {0.5}), V = sjnMat(1, 1, {1.0});
    const std::vector<double> N{3.0}, Z{1.0};
    const SjnResult lo = pfqn_mvasjn(L, N, Z, sjnMat(1, 1, {0.5}), {0}, V, SjnOptions());
    const SjnResult mid = pfqn_mvasjn(L, N, Z, sjnMat(1, 1, {1.0}), {0}, V, SjnOptions());
    const SjnResult hi = pfqn_mvasjn(L, N, Z, sjnMat(1, 1, {4.0}), {0}, V, SjnOptions());

    CHECK(lo.XN[0] == doctest::Approx(1.63001837090197).epsilon(1e-11));
    CHECK(lo.QN(0, 0) == doctest::Approx(1.36998162909803).epsilon(1e-11));
    CHECK(lo.CN(0, 0) == doctest::Approx(0.840470054542974).epsilon(1e-11));
    CHECK(hi.XN[0] == doctest::Approx(1.32157817139047).epsilon(1e-11));
    CHECK(hi.QN(0, 0) == doctest::Approx(1.67842182860953).epsilon(1e-11));
    CHECK(hi.CN(0, 0) == doctest::Approx(1.27001328029172).epsilon(1e-11));
    // Monotone in the SCV, which no single reference number above states.
    CHECK(lo.CN(0, 0) < mid.CN(0, 0));
    CHECK(mid.CN(0, 0) < hi.CN(0, 0));
    sjnLaws(lo, L, N, Z);
    sjnLaws(hi, L, N, Z);
}

TEST_CASE("the SJN multiclass readings are selected by options.prio") {
    // Same demands both ways. POOLED (prio empty) compares the jobs of every
    // class by size directly; PRIORITY (distinct levels, 1 = highest) is Kant's
    // method A, SJN applying only within a class. The two must NOT agree, which
    // is what proves the branch is reached rather than defaulted.
    const Matrix<double> L = sjnMat(1, 2, {0.5, 0.25}), V = sjnMat(1, 2, {1.0, 1.0});
    const Matrix<double> scv = sjnMat(1, 2, {1.0, 0.5});
    const std::vector<double> N{2.0, 2.0}, Z{1.0, 0.5};

    const SjnResult pooled = pfqn_mvasjn(L, N, Z, scv, {0}, V, SjnOptions());
    CHECK(pooled.XN[0] == doctest::Approx(0.929695526050463).epsilon(1e-11));
    CHECK(pooled.XN[1] == doctest::Approx(1.70651782674172).epsilon(1e-11));
    CHECK(pooled.QN(0, 0) == doctest::Approx(1.07030447394954).epsilon(1e-11));
    CHECK(pooled.QN(0, 1) == doctest::Approx(1.14674108662914).epsilon(1e-11));
    CHECK(pooled.CN(0, 0) == doctest::Approx(1.15124193239523).epsilon(1e-11));
    CHECK(pooled.CN(0, 1) == doctest::Approx(0.671977209179603).epsilon(1e-11));
    sjnLaws(pooled, L, N, Z);

    SjnOptions po;
    po.prio = {1, 2};
    const SjnResult prio = pfqn_mvasjn(L, N, Z, scv, {0}, V, po);
    CHECK(prio.XN[0] == doctest::Approx(1.13133658646935).epsilon(1e-11));
    CHECK(prio.XN[1] == doctest::Approx(0.963074357977113).epsilon(1e-11));
    CHECK(prio.QN(0, 0) == doctest::Approx(0.868663413530649).epsilon(1e-11));
    CHECK(prio.QN(0, 1) == doctest::Approx(1.51846282101144).epsilon(1e-11));
    CHECK(prio.CN(0, 0) == doctest::Approx(0.767820491195776).epsilon(1e-11));
    CHECK(prio.CN(0, 1) == doctest::Approx(1.57668284741886).epsilon(1e-11));
    sjnLaws(prio, L, N, Z);
    // Class 1 has the higher priority, so it gains and class 2 pays.
    CHECK(prio.CN(0, 0) < pooled.CN(0, 0));
    CHECK(prio.CN(0, 1) > pooled.CN(0, 1));
}

TEST_CASE("pfqn_amvasjn is the Schweitzer closure of the same station equation") {
    SUBCASE("it tracks the lattice on the same instance") {
        const Matrix<double> L = sjnMat(1, 1, {0.5}), V = sjnMat(1, 1, {1.0});
        const Matrix<double> scv = sjnMat(1, 1, {1.0});
        const std::vector<double> N{3.0}, Z{1.0};
        const SjnResult r = pfqn_amvasjn(L, N, Z, scv, {0}, V, SjnOptions());
        CHECK(r.XN[0] == doctest::Approx(1.58783888347139).epsilon(1e-9));
        CHECK(r.QN(0, 0) == doctest::Approx(1.41216111652861).epsilon(1e-9));
        CHECK(r.UN(0, 0) == doctest::Approx(0.793919441735695).epsilon(1e-9));
        CHECK(r.CN(0, 0) == doctest::Approx(0.889360457933423).epsilon(1e-9));
        CHECK(r.converged);
        CHECK(r.iter == 11);
        sjnLaws(r, L, N, Z);
    }
    SUBCASE("at light load the closure and the lattice nearly agree") {
        // The closure fixes the SHAPE of W(x) and lets only its level scale, so
        // its error grows with the load. At rho well below one the two routes
        // must agree to several digits, which is a check on the closure that no
        // single reference number provides.
        const Matrix<double> L = sjnMat(1, 1, {0.2}), V = sjnMat(1, 1, {1.0});
        const Matrix<double> scv = sjnMat(1, 1, {1.0});
        const std::vector<double> N{2.0}, Z{5.0};
        const SjnResult lat = pfqn_mvasjn(L, N, Z, scv, {0}, V, SjnOptions());
        const SjnResult fix = pfqn_amvasjn(L, N, Z, scv, {0}, V, SjnOptions());
        CHECK(lat.XN[0] == doctest::Approx(0.384041686852909).epsilon(1e-11));
        CHECK(fix.XN[0] == doctest::Approx(0.384036895550972).epsilon(1e-9));
        CHECK(lat.QN(0, 0) == doctest::Approx(0.0797915657354536).epsilon(1e-11));
        CHECK(fix.QN(0, 0) == doctest::Approx(0.0798155222451383).epsilon(1e-9));
        CHECK(fix.XN[0] == doctest::Approx(lat.XN[0]).epsilon(1e-4));
        sjnLaws(lat, L, N, Z);
        sjnLaws(fix, L, N, Z);
    }
    SUBCASE("the utilization cap binds and the population law survives it") {
        // N = 6 at rho -> 1 puts the station in the starvation regime, where the
        // conditional waiting time equation underestimates the residence time
        // and the implied throughput would exceed the station capacity. sjn_cap
        // bisects the waiting time until U = umax = 0.999 EXACTLY: that number
        // is the cap itself, not a solved value. Capping X instead would lose
        // jobs, so the population law below is the assertion that matters.
        const Matrix<double> L = sjnMat(1, 1, {0.5}), V = sjnMat(1, 1, {1.0});
        const Matrix<double> scv = sjnMat(1, 1, {1.0});
        const std::vector<double> N{6.0}, Z{1.0};
        const SjnResult r = pfqn_amvasjn(L, N, Z, scv, {0}, V, SjnOptions());
        CHECK(r.capped);
        CHECK(r.XN[0] == doctest::Approx(1.998).epsilon(1e-9));
        CHECK(r.QN(0, 0) == doctest::Approx(4.002).epsilon(1e-9));
        CHECK(r.UN(0, 0) == doctest::Approx(0.999).epsilon(1e-12));
        CHECK(r.CN(0, 0) == doctest::Approx(2.003003003003).epsilon(1e-9));
        CHECK(r.iter == 20);
        sjnLaws(r, L, N, Z);
    }
}

TEST_CASE("the SJN api rejects the arguments its equations cannot honour") {
    const Matrix<double> L = sjnMat(1, 1, {0.5}), V = sjnMat(1, 1, {1.0});
    const Matrix<double> scv = sjnMat(1, 1, {1.0});
    const std::vector<double> N{3.0}, Z{1.0};
    SUBCASE("an odd grid, which composite Simpson cannot integrate") {
        SjnOptions o;
        o.ns = 33;
        CHECK_THROWS_AS(pfqn_mvasjn(L, N, Z, scv, {0}, V, o), line::InputError);
    }
    SUBCASE("a utilization cap at or above one, where the equation is singular") {
        SjnOptions o;
        o.umax = 1.0;
        CHECK_THROWS_AS(pfqn_mvasjn(L, N, Z, scv, {0}, V, o), line::InputError);
    }
    SUBCASE("tied priority levels, which the SJN priority equations do not cover") {
        const Matrix<double> L2 = sjnMat(1, 2, {0.5, 0.25}), V2 = sjnMat(1, 2, {1.0, 1.0});
        const Matrix<double> scv2 = sjnMat(1, 2, {1.0, 1.0});
        SjnOptions o;
        o.prio = {1, 1};
        CHECK_THROWS_AS(pfqn_mvasjn(L2, {2.0, 2.0}, {1.0, 1.0}, scv2, {0}, V2, o),
                        line::InputError);
    }
    SUBCASE("a station index outside the demand matrix") {
        CHECK_THROWS_AS(pfqn_mvasjn(L, N, Z, scv, {1}, V, SjnOptions()), line::InputError);
    }
    SUBCASE("a repeated station index") {
        const Matrix<double> L2 = sjnMat(2, 1, {0.4, 0.3}), V2 = sjnMat(2, 1, {1.0, 1.0});
        const Matrix<double> scv2 = sjnMat(2, 1, {1.0, 1.0});
        CHECK_THROWS_AS(pfqn_mvasjn(L2, N, Z, scv2, {0, 0}, V2, SjnOptions()), line::InputError);
    }
    SUBCASE("a negative population") {
        CHECK_THROWS_AS(pfqn_mvasjn(L, {-1.0}, Z, scv, {0}, V, SjnOptions()), line::InputError);
    }
}

// ---------------------------------------------------------------------------
// pfqn_tay: arrival-instant AMVA from the throughput elasticities
//
// Reference values from calling MATLAB's pfqn_tay directly, printed at 15
// significant digits. What distinguishes Tay from every other driver in this
// file is WHERE the arrival-instant queue length comes from: not a population
// shift, but the solution of an R x R linear system in the throughput
// elasticities, one per (station, class) pair per sweep. A port that got the
// system's off-diagonal coupling or its right-hand side wrong still converges
// and still satisfies Little's law, so the reference numbers are doing real
// work here and the laws alone would not catch it.
// ---------------------------------------------------------------------------

TEST_CASE("pfqn_tay matches MATLAB on a two-class two-station model") {
    const std::vector<double> N{4.0, 3.0}, Z{1.0, 2.0};
    const auto r = pfqn_tay(amvaD(), N, Z);
    CHECK(r.XN[0] == doctest::Approx(1.14074445702499).epsilon(1e-9));
    CHECK(r.XN[1] == doctest::Approx(0.647505484723929).epsilon(1e-9));
    CHECK(r.QN(0, 0) == doctest::Approx(1.42592771135494).epsilon(1e-9));
    CHECK(r.QN(0, 1) == doctest::Approx(0.558853983419989).epsilon(1e-9));
    CHECK(r.QN(1, 0) == doctest::Approx(1.43332783162008).epsilon(1e-9));
    CHECK(r.QN(1, 1) == doctest::Approx(1.14613504713215).epsilon(1e-9));
    CHECK(r.UN(0, 0) == doctest::Approx(0.570372228512494).epsilon(1e-9));
    CHECK(r.UN(1, 1) == doctest::Approx(0.388503290834357).epsilon(1e-9));
    CHECK(r.RN(0, 0) == doctest::Approx(1.24999749292992).epsilon(1e-9));
    CHECK(r.RN(1, 1) == doctest::Approx(1.77007774323459).epsilon(1e-9));
    CHECK(r.iterations == 20);

    const Matrix<double> D = amvaD();
    for (std::size_t c = 0; c < 2; ++c) {
        // Population conservation with the think time, Little's law and the
        // utilization law: the elasticities may be wrong and these still hold.
        double tot = 0.0;
        for (std::size_t i = 0; i < 2; ++i) {
            tot += r.QN(i, c);
            CHECK(r.QN(i, c) == doctest::Approx(r.XN[c] * r.RN(i, c)).epsilon(1e-9));
            CHECK(r.UN(i, c) == doctest::Approx(r.XN[c] * D(i, c)).epsilon(1e-9));
        }
        CHECK(tot + r.XN[c] * Z[c] == doctest::Approx(N[c]).epsilon(1e-9));
    }
}

TEST_CASE("pfqn_tay without a think time is a different branch of the elasticity system") {
    // Z enters the elasticity DENOMINATOR only, as the AS-server term Z_j X_j.
    // Dropping it is therefore not a rescaling of the answer, and a port that
    // put Z in the numerator, or omitted it, agrees with the case above and
    // disagrees here.
    const std::vector<double> N{4.0, 3.0}, Z0{0.0, 0.0};
    const auto r = pfqn_tay(amvaD(), N, Z0);
    CHECK(r.XN[0] == doctest::Approx(1.13832431058883).epsilon(1e-9));
    CHECK(r.XN[1] == doctest::Approx(0.812884010194267).epsilon(1e-9));
    CHECK(r.QN(0, 0) == doctest::Approx(1.78681103057868).epsilon(1e-9));
    CHECK(r.QN(0, 1) == doctest::Approx(0.90184324893248).epsilon(1e-9));
    CHECK(r.QN(1, 0) == doctest::Approx(2.21318896942132).epsilon(1e-9));
    CHECK(r.QN(1, 1) == doctest::Approx(2.09815675106752).epsilon(1e-9));
    // With no think time every job is at a station, so the columns sum to N.
    CHECK(r.QN(0, 0) + r.QN(1, 0) == doctest::Approx(4.0).epsilon(1e-9));
    CHECK(r.QN(0, 1) + r.QN(1, 1) == doctest::Approx(3.0).epsilon(1e-9));
}

TEST_CASE("pfqn_tay solves an empty class out rather than through") {
    // An empty class makes its own elasticity denominator identically zero, so
    // the R x R system is SINGULAR, not merely redundant: leaving the class in
    // does not give a large answer, it gives no answer. The reference drops the
    // class, solves, and re-expands.
    //
    // The invariant is exact and needs no second reference: N = [4, 0] on the
    // two-class model must reproduce the single-class solve on column 0 alone,
    // which MATLAB puts at X = 1.48183415477026.
    const std::vector<double> Z{1.0, 2.0};
    const auto r = pfqn_tay(amvaD(), {4.0, 0.0}, Z);
    CHECK(r.XN[0] == doctest::Approx(1.48183415477026).epsilon(1e-9));
    CHECK(r.XN[1] == doctest::Approx(0.0));
    CHECK(r.QN(0, 0) == doctest::Approx(1.48027058405875).epsilon(1e-9));
    CHECK(r.QN(1, 0) == doctest::Approx(1.03789526117099).epsilon(1e-9));
    CHECK(r.QN(0, 1) == doctest::Approx(0.0));
    CHECK(r.QN(1, 1) == doctest::Approx(0.0));

    // The same model with the empty class actually absent, solved directly.
    Matrix<double> D1(2, 1);
    D1(0, 0) = 0.5;
    D1(1, 0) = 0.4;
    const auto s = pfqn_tay(D1, std::vector<double>{4.0}, std::vector<double>{1.0});
    CHECK(s.XN[0] == doctest::Approx(r.XN[0]).epsilon(1e-12));
    CHECK(s.QN(0, 0) == doctest::Approx(r.QN(0, 0)).epsilon(1e-12));
    CHECK(s.QN(1, 0) == doctest::Approx(r.QN(1, 0)).epsilon(1e-12));
    CHECK(s.RN(0, 0) == doctest::Approx(0.998944840955063).epsilon(1e-9));
    CHECK(s.RN(1, 0) == doctest::Approx(0.700412564948541).epsilon(1e-9));
    CHECK(s.iterations == 13);
    // Single class, single server, so Tay is close to exact MVA but not equal:
    // MATLAB's exact X on this model is 1.48328323985994.
    CHECK(s.XN[0] == doctest::Approx(1.48328323985994).epsilon(1e-2));
    CHECK(s.XN[0] < 1.48328323985994);
}

TEST_CASE("pfqn_tay rejects the shapes its recursion cannot take") {
    SUBCASE("a population vector of the wrong length") {
        CHECK_THROWS_AS(pfqn_tay(amvaD(), std::vector<double>{4.0}, std::vector<double>{}),
                        line::InputError);
    }
    SUBCASE("a think-time vector of the wrong length") {
        CHECK_THROWS_AS(
            pfqn_tay(amvaD(), std::vector<double>{4.0, 3.0}, std::vector<double>{1.0}),
            line::InputError);
    }
    SUBCASE("an all-empty population returns zeros rather than dividing by zero") {
        const auto r = pfqn_tay(amvaD(), std::vector<double>{0.0, 0.0}, std::vector<double>{});
        CHECK(r.XN[0] == doctest::Approx(0.0));
        CHECK(r.XN[1] == doctest::Approx(0.0));
        CHECK(r.QN(0, 0) == doctest::Approx(0.0));
    }
}
