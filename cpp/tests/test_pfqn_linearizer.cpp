/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * The Linearizer family of approximate MVA: pfqn_linearizer, pfqn_gflinearizer,
 * pfqn_egflinearizer, pfqn_linearizermx, pfqn_linearizerms, pfqn_conwayms
 * and pfqn_schmidt.
 *
 * Every reference value below was produced by running the MATLAB original on
 * the same model (matlab/src/api/pfqn), printed at %.14g. Unless stated
 * otherwise the assertions sit at the tolerance the method itself stops on --
 * 1e-8 on enorm(dQ) for the Linearizer family -- and NOT tighter: the last digits of a fixed point stopped on a tolerance are
 * an artifact of the stopping test, not a property of the algorithm. The
 * observed agreement is in fact far better -- the worst relative deviation
 * over every quantity of every function below is 4.1e-14, i.e. the references
 * are reproduced to the last digit MATLAB printed -- but asserting that would
 * be asserting the stopping rule rather than the algorithm.
 *
 * pfqn_schmidt is the exception: it is a finite recursion over the population
 * lattice with no tolerance at all, so it is held to 1e-12, which is rounding.
 */
#include <cmath>
#include <vector>

#include "doctest.h"
#include "line/api/pfqn/pfqn_bs.h"
#include "line/api/pfqn/pfqn_conwayms.h"
#include "line/api/pfqn/pfqn_egflinearizer.h"
#include "line/api/pfqn/pfqn_gflinearizer.h"
#include "line/api/pfqn/pfqn_linearizer.h"
#include "line/api/pfqn/pfqn_linearizermx.h"
#include "line/api/pfqn/pfqn_linearizerms.h"
#include "line/api/pfqn/pfqn_mva.h"
#include "line/api/pfqn/pfqn_schmidt.h"

using line::Matrix;
using namespace line::pfqn;

namespace {

// Two stations, two classes, heterogeneous demands: neither class is dominant
// at either station, so the Linearizer correction is actually exercised.
Matrix<double> demands() {
    Matrix<double> L(2, 2);
    L(0, 0) = 0.5; L(0, 1) = 0.3;
    L(1, 0) = 0.4; L(1, 1) = 0.6;
    return L;
}

// Think time, present on both classes: without it the closed subnetwork would
// degenerate and the C = N/X - Z branch would never be tested.
Matrix<double> thinktimes() {
    Matrix<double> Z(1, 2);
    Z(0, 0) = 1.0; Z(0, 1) = 2.0;
    return Z;
}

const std::vector<int> kN{4, 3};
const std::vector<double> kZ{1.0, 2.0};
const double kTol = 1e-8;   // the tolerance every Linearizer variant stops on
const int kMaxIter = 1000;

std::vector<SchedStrategy> allPS() {
    return std::vector<SchedStrategy>{SchedStrategy::PS, SchedStrategy::PS};
}

/**
 * Little's law per class on a closed network: the jobs of class r are either
 * queued somewhere or thinking, N_r = sum_i Q(i,r) + X_r Z_r. This is a
 * conservation identity the algorithm must satisfy exactly at its fixed point,
 * independently of how accurate the approximation is.
 */
void checkLittleClosed(const Matrix<double>& Q, const std::vector<double>& X,
                       const std::vector<double>& Z, const std::vector<int>& N, double eps) {
    for (std::size_t r = 0; r < N.size(); ++r) {
        double q = 0.0;
        for (std::size_t i = 0; i < Q.rows(); ++i) q += Q(i, r);
        CHECK(q + X[r] * Z[r] == doctest::Approx(static_cast<double>(N[r])).epsilon(eps));
    }
}

/** Little's law per station-class: Q(i,r) = X_r W(i,r). */
void checkLittleStation(const Matrix<double>& Q, const Matrix<double>& W,
                        const std::vector<double>& X, double eps) {
    for (std::size_t i = 0; i < Q.rows(); ++i)
        for (std::size_t r = 0; r < Q.cols(); ++r)
            CHECK(Q(i, r) == doctest::Approx(X[r] * W(i, r)).epsilon(eps));
}

/** Utilization law: U(i,r) = X_r L(i,r) / c_i. */
void checkUtilLaw(const Matrix<double>& U, const Matrix<double>& L, const std::vector<double>& X,
                  const std::vector<int>& nservers, double eps) {
    for (std::size_t i = 0; i < U.rows(); ++i)
        for (std::size_t r = 0; r < U.cols(); ++r)
            CHECK(U(i, r) ==
                  doctest::Approx(X[r] * L(i, r) / static_cast<double>(nservers[i])).epsilon(eps));
}

const std::vector<int> kSingle{1, 1};

}  // namespace

// --------------------------------------------------------------------------
// pfqn_linearizer
// --------------------------------------------------------------------------

TEST_CASE("pfqn_linearizer matches the MATLAB reference") {
    // MATLAB: pfqn_linearizer([0.5 0.3;0.4 0.6],[4 3],[1 2],[PS;PS],1e-8,1000)
    const LinearizerResult<double> a =
        pfqn_linearizer(demands(), kN, thinktimes(), allPS(), kTol, kMaxIter, Matrix<double>());
    CHECK(a.X[0] == doctest::Approx(1.1427873365872).epsilon(kTol));
    CHECK(a.X[1] == doctest::Approx(0.64809976902412).epsilon(kTol));
    CHECK(a.Q(0, 0) == doctest::Approx(1.4333313783761).epsilon(kTol));
    CHECK(a.Q(0, 1) == doctest::Approx(0.55907885909637).epsilon(kTol));
    CHECK(a.Q(1, 0) == doctest::Approx(1.4238812850366).epsilon(kTol));
    CHECK(a.Q(1, 1) == doctest::Approx(1.1447216028554).epsilon(kTol));
    CHECK(a.U(0, 0) == doctest::Approx(0.57139366829362).epsilon(kTol));
    CHECK(a.U(1, 1) == doctest::Approx(0.38885986141447).epsilon(kTol));
    CHECK(a.W(0, 0) == doctest::Approx(1.2542415657637).epsilon(kTol));
    CHECK(a.W(1, 1) == doctest::Approx(1.7662737398889).epsilon(kTol));
    CHECK(a.C[0] == doctest::Approx(2.5002137947603).epsilon(kTol));
    CHECK(a.C[1] == doctest::Approx(2.6289169405464).epsilon(kTol));
    // MATLAB reports totiter = 265. This is not a numerical assertion but a
    // structural one: the same inner-loop schedule visited the same sequence
    // of populations, which is what makes the agreement above meaningful
    // rather than coincidental.
    CHECK(a.totiter == 265);
}

TEST_CASE("pfqn_linearizer obeys Little's law and the utilization law") {
    const LinearizerResult<double> a =
        pfqn_linearizer(demands(), kN, thinktimes(), allPS(), kTol, kMaxIter, Matrix<double>());
    checkLittleClosed(a.Q, a.X, kZ, kN, 1e-9);
    checkLittleStation(a.Q, a.W, a.X, 1e-9);
    checkUtilLaw(a.U, demands(), a.X, kSingle, 1e-9);
    // C_r = N_r/X_r - Z_r is the cycle time excluding the delay.
    for (std::size_t r = 0; r < 2; ++r) {
        double w = 0.0;
        for (std::size_t i = 0; i < 2; ++i) w += a.W(i, r);
        CHECK(a.C[r] == doctest::Approx(w).epsilon(1e-9));
    }
    // Well-posed model: the run finished inside the iteration budget rather
    // than being cut off by it, which is what convergence means for a routine
    // whose only signal is the iteration count.
    CHECK(a.totiter < kMaxIter);
}

TEST_CASE("pfqn_linearizer is at least as accurate as Bard-Schweitzer") {
    // No think time, so exact MVA is directly comparable.
    // MATLAB pfqn_mva gives X = [1.1379298638694, 0.81284062035666];
    // pfqn_bs gives [1.1396670114877, 0.7862446110151];
    // pfqn_linearizer gives [1.13757741313, 0.81344542017133].
    const Matrix<double> L = demands();
    const MvaResult<double> exact = pfqn_mva(L, kN);
    const LinearizerResult<double> lin =
        pfqn_linearizer(L, kN, Matrix<double>(), allPS(), kTol, kMaxIter, Matrix<double>());
    const AmvaResult<double> bs =
        pfqn_bs(L, std::vector<double>{4.0, 3.0}, std::vector<double>{0.0, 0.0});
    CHECK(exact.XN[0] == doctest::Approx(1.1379298638694).epsilon(1e-12));
    CHECK(exact.XN[1] == doctest::Approx(0.81284062035666).epsilon(1e-12));
    CHECK(lin.X[0] == doctest::Approx(1.13757741313).epsilon(kTol));
    CHECK(lin.X[1] == doctest::Approx(0.81344542017133).epsilon(kTol));

    for (std::size_t r = 0; r < 2; ++r) {
        const double el = std::fabs(lin.X[r] - exact.XN[r]) / exact.XN[r];
        const double eb = std::fabs(bs.XN[r] - exact.XN[r]) / exact.XN[r];
        INFO("class ", r, ": Linearizer ", el, " vs Bard-Schweitzer ", eb);
        CHECK(el <= eb);
        CHECK(el < 1e-3);
    }
    // The same ordering on the queue lengths, which is what Linearizer's
    // Delta correction actually refines.
    for (std::size_t i = 0; i < 2; ++i)
        for (std::size_t r = 0; r < 2; ++r) {
            const double el = std::fabs(lin.Q(i, r) - exact.QN(i, r)) / exact.QN(i, r);
            const double eb = std::fabs(bs.QN(i, r) - exact.QN(i, r)) / exact.QN(i, r);
            INFO("station ", i, " class ", r, ": Linearizer ", el, " vs BS ", eb);
            CHECK(el <= eb);
        }
}

// --------------------------------------------------------------------------
// pfqn_gflinearizer / pfqn_egflinearizer
// --------------------------------------------------------------------------

TEST_CASE("pfqn_gflinearizer matches the MATLAB reference at alpha = 2") {
    // MATLAB: pfqn_gflinearizer(L,[4 3],[1 2],[PS;PS],1e-8,1000,2.0)
    const LinearizerResult<double> g = pfqn_gflinearizer(demands(), kN, thinktimes(), allPS(),
                                                         kTol, kMaxIter, 2.0, Matrix<double>());
    CHECK(g.X[0] == doctest::Approx(1.147509388812).epsilon(kTol));
    CHECK(g.X[1] == doctest::Approx(0.64760660144857).epsilon(kTol));
    CHECK(g.Q(0, 0) == doctest::Approx(1.4135524106165).epsilon(kTol));
    CHECK(g.Q(0, 1) == doctest::Approx(0.5508780735473).epsilon(kTol));
    CHECK(g.Q(1, 0) == doctest::Approx(1.4389382005715).epsilon(kTol));
    CHECK(g.Q(1, 1) == doctest::Approx(1.1539087235556).epsilon(kTol));
    checkLittleClosed(g.Q, g.X, kZ, kN, 1e-9);
    checkLittleStation(g.Q, g.W, g.X, 1e-9);
    checkUtilLaw(g.U, demands(), g.X, kSingle, 1e-9);
    CHECK(g.totiter < kMaxIter);
}

TEST_CASE("pfqn_egflinearizer matches the MATLAB reference at a per-class alpha") {
    // MATLAB: pfqn_egflinearizer(L,[4 3],[1 2],[PS;PS],1e-8,1000,[0.9 1.1])
    const std::vector<double> alpha{0.9, 1.1};
    const LinearizerResult<double> e = pfqn_egflinearizer(demands(), kN, thinktimes(), allPS(),
                                                          kTol, kMaxIter, alpha, Matrix<double>());
    CHECK(e.X[0] == doctest::Approx(1.1425284834319).epsilon(kTol));
    CHECK(e.X[1] == doctest::Approx(0.64832918672096).epsilon(kTol));
    CHECK(e.Q(0, 0) == doctest::Approx(1.4346134268157).epsilon(kTol));
    CHECK(e.Q(0, 1) == doctest::Approx(0.55992750719538).epsilon(kTol));
    CHECK(e.Q(1, 0) == doctest::Approx(1.4228580897524).epsilon(kTol));
    CHECK(e.Q(1, 1) == doctest::Approx(1.1434141193627).epsilon(kTol));
    checkLittleClosed(e.Q, e.X, kZ, kN, 1e-9);
    checkLittleStation(e.Q, e.W, e.X, 1e-9);
    checkUtilLaw(e.U, demands(), e.X, kSingle, 1e-9);
    CHECK(e.totiter < kMaxIter);
}

TEST_CASE("pfqn_linearizer is pfqn_egflinearizer at alpha = 1") {
    const std::vector<double> ones{1.0, 1.0};
    const LinearizerResult<double> a =
        pfqn_linearizer(demands(), kN, thinktimes(), allPS(), kTol, kMaxIter, Matrix<double>());
    const LinearizerResult<double> e = pfqn_egflinearizer(demands(), kN, thinktimes(), allPS(),
                                                          kTol, kMaxIter, ones, Matrix<double>());
    for (std::size_t r = 0; r < 2; ++r) CHECK(a.X[r] == doctest::Approx(e.X[r]).epsilon(1e-14));
}

// --------------------------------------------------------------------------
// pfqn_linearizerms
// --------------------------------------------------------------------------

TEST_CASE("pfqn_linearizerms matches the MATLAB reference") {
    // MATLAB: pfqn_linearizerms(L,[4 3],[1 2],[2;1],[PS;PS],1e-8,1000)
    //
    // RE-RECORDED 2026-07-31 against the current reference. The previous values
    // came from MATLAB BEFORE its own multiserver-Linearizer fix, which corrected
    // three defects at once: the partially-idle-server term summed every class's
    // demand instead of the arriving job's own over m, the Estimate step left the
    // marginals at population N while reducing the queue lengths, and the marginal
    // recursion was iterated as a Jacobi sweep that amplifies by the number of busy
    // servers. This port now reproduces the fixed reference to 13 digits, at the
    // same iteration count (159); the numbers below were obtained by running MATLAB
    // on these very inputs, not by recording this port's output.
    const std::vector<int> ns{2, 1};
    const LinearizerResult<double> m = pfqn_linearizerms(demands(), kN, thinktimes(), ns, allPS(),
                                                         kTol, kMaxIter, Matrix<double>());
    CHECK(m.X[0] == doctest::Approx(1.33563463439118).epsilon(kTol));
    CHECK(m.X[1] == doctest::Approx(0.66966793510685).epsilon(kTol));
    CHECK(m.Q(0, 0) == doctest::Approx(0.76137067654236).epsilon(kTol));
    CHECK(m.Q(0, 1) == doctest::Approx(0.23907832616540).epsilon(kTol));
    CHECK(m.Q(1, 0) == doctest::Approx(1.90299468906646).epsilon(kTol));
    CHECK(m.Q(1, 1) == doctest::Approx(1.42158580362090).epsilon(kTol));
    CHECK(m.U(0, 0) == doctest::Approx(0.33390865859779).epsilon(kTol));
    CHECK(m.U(0, 1) == doctest::Approx(0.10045019026603).epsilon(kTol));
    CHECK(m.U(1, 0) == doctest::Approx(0.53425385375647).epsilon(kTol));
    CHECK(m.U(1, 1) == doctest::Approx(0.40180076106411).epsilon(kTol));
    checkLittleClosed(m.Q, m.X, kZ, kN, 1e-9);
    checkLittleStation(m.Q, m.W, m.X, 1e-9);
    checkUtilLaw(m.U, demands(), m.X, ns, 1e-9);
    CHECK(m.totiter < kMaxIter);
    // Adding a server can only speed the network up.
    const LinearizerResult<double> one =
        pfqn_linearizerms(demands(), kN, thinktimes(), kSingle, allPS(), kTol, kMaxIter,
                          Matrix<double>());
    CHECK(m.X[0] > one.X[0]);
    CHECK(m.X[1] > one.X[1]);
}

// --------------------------------------------------------------------------
// pfqn_conwayms
// --------------------------------------------------------------------------

TEST_CASE("pfqn_conwayms matches the MATLAB reference") {
    // MATLAB: pfqn_conwayms(L,[4 3],[1 2],[2;1]) -- note that this routine
    // defaults `type` to all-FCFS, not all-PS like the rest of the family, so
    // the reference call deliberately omits it.
    const std::vector<int> ns{2, 1};
    const LinearizerResult<double> c =
        pfqn_conwayms(demands(), kN, thinktimes(), ns, std::vector<SchedStrategy>(), kTol,
                      kMaxIter, Matrix<double>());
    CHECK(c.X[0] == doctest::Approx(1.2348029172801).epsilon(kTol));
    CHECK(c.X[1] == doctest::Approx(0.71890458124158).epsilon(kTol));
    CHECK(c.Q(0, 0) == doctest::Approx(0.74888304630864).epsilon(kTol));
    CHECK(c.Q(0, 1) == doctest::Approx(0.31365809957937).epsilon(kTol));
    CHECK(c.Q(1, 0) == doctest::Approx(2.0163140364112).epsilon(kTol));
    CHECK(c.Q(1, 1) == doctest::Approx(1.2485327379375).epsilon(kTol));
    CHECK(c.U(0, 0) == doctest::Approx(0.30870072932003).epsilon(kTol));
    CHECK(c.U(1, 1) == doctest::Approx(0.43134274874495).epsilon(kTol));
    CHECK(c.W(0, 0) == doctest::Approx(0.60647981619463).epsilon(kTol));
    CHECK(c.W(1, 1) == doctest::Approx(1.7367155120659).epsilon(kTol));
    checkLittleClosed(c.Q, c.X, kZ, kN, 1e-9);
    checkLittleStation(c.Q, c.W, c.X, 1e-9);
    checkUtilLaw(c.U, demands(), c.X, ns, 1e-9);
    CHECK(c.totiter < kMaxIter);
}

// --------------------------------------------------------------------------
// pfqn_schmidt
// --------------------------------------------------------------------------

TEST_CASE("pfqn_schmidt matches the MATLAB reference on a multiserver PS model") {
    // MATLAB: pfqn_schmidt(L,[4 3],[2;1],[PS;PS]). Held to 1e-12, not to a
    // method tolerance: pfqn_schmidt is a finite recursion over the population
    // lattice with no stopping test, so the only error is rounding.
    Matrix<int> S(2, 1);
    S(0, 0) = 2;
    S(1, 0) = 1;
    const std::vector<SchedStrategy> sched = allPS();
    const SchmidtResult<double> s = pfqn_schmidt(demands(), kN, S, sched);
    CHECK(s.XN[0] == doctest::Approx(1.3281255990795).epsilon(1e-12));
    CHECK(s.XN[1] == doctest::Approx(0.77706555316029).epsilon(1e-12));
    CHECK(s.QN(0, 0) == doctest::Approx(0.81731274242371).epsilon(1e-12));
    CHECK(s.QN(0, 1) == doctest::Approx(0.30753533861694).epsilon(1e-12));
    CHECK(s.QN(1, 0) == doctest::Approx(3.1826872575763).epsilon(1e-12));
    CHECK(s.QN(1, 1) == doctest::Approx(2.6924646613831).epsilon(1e-12));
    CHECK(s.UN(0, 0) == doctest::Approx(0.33203139976988).epsilon(1e-12));
    CHECK(s.UN(0, 1) == doctest::Approx(0.11655983297404).epsilon(1e-12));
    CHECK(s.UN(1, 0) == doctest::Approx(0.53125023963181).epsilon(1e-12));
    CHECK(s.UN(1, 1) == doctest::Approx(0.46623933189617).epsilon(1e-12));
    CHECK(s.CN(0, 0) == doctest::Approx(0.61538814024077).epsilon(1e-12));
    CHECK(s.CN(1, 1) == doctest::Approx(3.4649131600712).epsilon(1e-12));

    const std::vector<int> ns{2, 1};
    checkLittleClosed(s.QN, s.XN, std::vector<double>{0.0, 0.0}, kN, 1e-12);
    checkLittleStation(s.QN, s.CN, s.XN, 1e-12);
    checkUtilLaw(s.UN, demands(), s.XN, ns, 1e-12);
}

TEST_CASE("pfqn_schmidt reduces to exact MVA on single-server exponential FCFS") {
    // A single-server FCFS station with exponential service is product-form,
    // so Schmidt's recursion must return the exact MVA solution, not an
    // approximation of it. MATLAB agrees to the printed digits.
    Matrix<int> S(2, 1);
    S(0, 0) = 1;
    S(1, 0) = 1;
    const std::vector<SchedStrategy> sched{SchedStrategy::FCFS, SchedStrategy::FCFS};
    const SchmidtResult<double> s = pfqn_schmidt(demands(), kN, S, sched);
    const MvaResult<double> exact = pfqn_mva(demands(), kN);
    CHECK(s.XN[0] == doctest::Approx(1.1379298638694).epsilon(1e-12));
    CHECK(s.XN[1] == doctest::Approx(0.81284062035666).epsilon(1e-12));
    CHECK(s.QN(0, 0) == doctest::Approx(1.7931221780901).epsilon(1e-12));
    CHECK(s.QN(1, 1) == doctest::Approx(2.098365022288).epsilon(1e-12));
    for (std::size_t r = 0; r < 2; ++r) {
        CHECK(s.XN[r] == doctest::Approx(exact.XN[r]).epsilon(1e-10));
        for (std::size_t i = 0; i < 2; ++i)
            CHECK(s.QN(i, r) == doctest::Approx(exact.QN(i, r)).epsilon(1e-10));
    }
}

// --------------------------------------------------------------------------
// pfqn_linearizermx
// --------------------------------------------------------------------------

namespace {

// Class 1 open at lambda = 0.4, class 2 closed at N = 3 with think time 2.
Matrix<double> mxThink() {
    Matrix<double> Z(1, 2);
    Z(0, 0) = 0.0;
    Z(0, 1) = 2.0;
    return Z;
}
const std::vector<double> kLambda{0.4, 0.0};
const std::vector<int> kNmx{kOpenClass, 3};

LinearizerResult<double> runMx(LinearizerMxMethod m, const std::vector<int>& ns) {
    return pfqn_linearizermx(kLambda, demands(), kNmx, mxThink(), ns, allPS(), kTol, kMaxIter, m,
                             Matrix<double>());
}

}  // namespace

TEST_CASE("pfqn_linearizermx matches the MATLAB reference with method 'lin'") {
    // MATLAB: pfqn_linearizermx([0.4 0],L,[Inf 3],[0 2],[1;1],[PS;PS],1e-8,1000,'lin')
    const LinearizerResult<double> a = runMx(LinearizerMxMethod::Lin, kSingle);
    CHECK(a.X[0] == doctest::Approx(0.4).epsilon(1e-14));
    CHECK(a.X[1] == doctest::Approx(0.84186985221022).epsilon(kTol));
    CHECK(a.Q(0, 0) == doctest::Approx(0.34914456045114).epsilon(kTol));
    CHECK(a.Q(0, 1) == doctest::Approx(0.39657824180458).epsilon(kTol));
    CHECK(a.Q(1, 0) == doctest::Approx(0.36565372452857).epsilon(kTol));
    CHECK(a.Q(1, 1) == doctest::Approx(0.91968205377498).epsilon(kTol));
    CHECK(a.U(0, 0) == doctest::Approx(0.2).epsilon(1e-12));
    CHECK(a.U(0, 1) == doctest::Approx(0.25256095566307).epsilon(kTol));
    CHECK(a.U(1, 0) == doctest::Approx(0.16).epsilon(1e-12));
    CHECK(a.U(1, 1) == doctest::Approx(0.50512191132613).epsilon(kTol));
    CHECK(a.W(0, 0) == doctest::Approx(0.87286140112786).epsilon(kTol));
    CHECK(a.W(0, 1) == doctest::Approx(0.47106834953575).epsilon(kTol));
    CHECK(a.W(1, 0) == doctest::Approx(0.91413431132142).epsilon(kTol));
    CHECK(a.W(1, 1) == doctest::Approx(1.0924278276036).epsilon(kTol));
    CHECK(a.C[0] == doctest::Approx(1.7869957124493).epsilon(kTol));
    CHECK(a.C[1] == doctest::Approx(1.5634961771394).epsilon(kTol));
    CHECK(a.totiter < kMaxIter);
}

TEST_CASE("pfqn_linearizermx obeys Little's law and the utilization law") {
    const LinearizerResult<double> a = runMx(LinearizerMxMethod::Lin, kSingle);
    // Little's law on every class, station by station.
    checkLittleStation(a.Q, a.W, a.X, 1e-9);
    // The utilization law uses the ORIGINAL demands, not the open-corrected
    // ones the closed subnetwork was solved with.
    checkUtilLaw(a.U, demands(), a.X, kSingle, 1e-9);
    // Closed class: the population identity still holds.
    double q = 0.0;
    for (std::size_t i = 0; i < 2; ++i) q += a.Q(i, 1);
    CHECK(q + a.X[1] * 2.0 == doctest::Approx(3.0).epsilon(1e-9));
    // Open class: throughput is the given arrival rate and the cycle time is
    // the sum of the residence times.
    CHECK(a.X[0] == doctest::Approx(0.4).epsilon(1e-14));
    CHECK(a.C[0] == doctest::Approx(a.W(0, 0) + a.W(1, 0)).epsilon(1e-12));
}

TEST_CASE("pfqn_linearizermx matches the MATLAB reference with method 'gflin'") {
    const LinearizerResult<double> g = runMx(LinearizerMxMethod::Gflin, kSingle);
    CHECK(g.X[1] == doctest::Approx(0.84550240432855).epsilon(kTol));
    CHECK(g.Q(0, 0) == doctest::Approx(0.34916578447467).epsilon(kTol));
    CHECK(g.Q(0, 1) == doctest::Approx(0.39666313789869).epsilon(kTol));
    CHECK(g.Q(1, 0) == doctest::Approx(0.36425372446556).epsilon(kTol));
    CHECK(g.Q(1, 1) == doctest::Approx(0.91233205344421).epsilon(kTol));
    checkLittleStation(g.Q, g.W, g.X, 1e-9);
    checkUtilLaw(g.U, demands(), g.X, kSingle, 1e-9);
}

TEST_CASE("pfqn_linearizermx 'egflin' uses a closed-class-indexed Gompertz alpha") {
    // MATLAB pfqn_linearizermx builds alphaM over the GLOBAL class indices but
    // pfqn_egflinearizer consumes it over the CLOSED ones, so here -- class 1
    // open, class 2 closed -- the closed class silently receives alphaM(1) = 0
    // and the Gompertz scaling is switched off. The JAR
    // (Pfqn_linearizermx.java) indexes alphaM over Nclosed and is correct;
    // this port follows the JAR. See the header of pfqn_linearizermx.h.
    //
    // Reference produced by driving MATLAB's own pfqn_egflinearizer through
    // the pfqn_linearizermx wrapper arithmetic with the intended exponent
    // alpha = 0.6 + 1.4*exp(-8*exp(-0.8*3)) = 1.2775503648739.
    const LinearizerResult<double> e = runMx(LinearizerMxMethod::Egflin, kSingle);
    CHECK(e.X[0] == doctest::Approx(0.4).epsilon(1e-14));
    CHECK(e.X[1] == doctest::Approx(0.8426013436523).epsilon(kTol));
    CHECK(e.Q(0, 0) == doctest::Approx(0.34909741175553).epsilon(kTol));
    CHECK(e.Q(0, 1) == doctest::Approx(0.39638964702213).epsilon(kTol));
    CHECK(e.Q(1, 0) == doctest::Approx(0.36541098393777).epsilon(kTol));
    CHECK(e.Q(1, 1) == doctest::Approx(0.91840766567328).epsilon(kTol));
    CHECK(e.W(0, 1) == doctest::Approx(0.47043557431793).epsilon(kTol));
    CHECK(e.C[1] == doctest::Approx(1.5604025825503).epsilon(kTol));
    // MATLAB as shipped returns 0.85177305755470 here, which is exactly what a
    // direct pfqn_egflinearizer call with alpha = 0 produces; that value is
    // the defect, so the port must NOT reproduce it.
    CHECK(std::fabs(e.X[1] - 0.85177305755470) > 1e-4);
    checkLittleStation(e.Q, e.W, e.X, 1e-9);
    checkUtilLaw(e.U, demands(), e.X, kSingle, 1e-9);
}

TEST_CASE("pfqn_linearizermx routes to the multiserver Linearizer above one server") {
    // MATLAB: pfqn_linearizermx(...,[2;1],...,'lin'). The method argument is
    // ignored on this branch, exactly as in the reference.
    const std::vector<int> ns{2, 1};
    // RE-RECORDED 2026-07-31 with the case above, and for the same reason: this
    // branch delegates to pfqn_linearizerms, so it moved when that reference did.
    const LinearizerResult<double> m = runMx(LinearizerMxMethod::Lin, ns);
    CHECK(m.X[0] == doctest::Approx(0.4).epsilon(1e-14));
    CHECK(m.X[1] == doctest::Approx(0.86205332202134).epsilon(kTol));
    CHECK(m.Q(0, 0) == doctest::Approx(0.33136779471920).epsilon(kTol));
    CHECK(m.Q(0, 1) == doctest::Approx(0.32547117887682).epsilon(kTol));
    CHECK(m.Q(1, 0) == doctest::Approx(0.37150898611057).epsilon(kTol));
    CHECK(m.Q(1, 1) == doctest::Approx(0.95042217708051).epsilon(kTol));
    const LinearizerResult<double> ignored = runMx(LinearizerMxMethod::Egflin, ns);
    CHECK(ignored.X[1] == doctest::Approx(m.X[1]).epsilon(1e-14));
}

TEST_CASE("pfqn_linearizermx refuses ill-posed mixed models") {
    // An arrival rate on a class that also carries a finite population.
    const std::vector<int> Nbad{2, 3};
    CHECK_THROWS(pfqn_linearizermx(kLambda, demands(), Nbad, mxThink(), kSingle, allPS(), kTol,
                                   kMaxIter, LinearizerMxMethod::Lin, Matrix<double>()));
    // Open traffic that saturates a station leaves no capacity for the closed
    // classes: 1/(1 - U_open) is not a usable demand correction. Station 1 has
    // demand 0.5 for the open class, so lambda = 2.5 drives it to U = 1.25.
    const std::vector<double> lamBig{2.5, 0.0};
    CHECK_THROWS(pfqn_linearizermx(lamBig, demands(), kNmx, mxThink(), kSingle, allPS(), kTol,
                                   kMaxIter, LinearizerMxMethod::Lin, Matrix<double>()));
}
