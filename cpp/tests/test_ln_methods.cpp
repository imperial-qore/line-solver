/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * The @@SolverLN entry points beyond the mean fixed point: `moment3`, the
 * response-time CDF, the transient, the sensitivity table, the Majumdar-Woodside
 * box bounds and the alternative layer engines.
 *
 * ONE MODEL throughout, `randomLQN/model_C1_L2_T2_P2_c4_z1_s1_t1_p1_y1_4.lqnx`:
 * a two-layer chain c0 -> t0 -> t1 on three processors, closed with four jobs.
 * It is small enough that a fluid layer transient and a finite-difference
 * sensitivity sweep both run in well under a second, and deep enough that a
 * dropped normalisation shows up (t0 and t1 sit on different hosts, so an entry
 * whose service is charged in the wrong index space cannot hide).
 *
 * Reference numbers are MATLAB SolverLN on the identical file (2026-07-31,
 * defaultOptions plus the method under test). Where the assertion is structural
 * rather than numeric -- a CDF is monotone, a bound brackets the fixed point,
 * two branches of the same derivative agree -- it is written as such, because
 * those are the properties a port can violate while still printing plausible
 * numbers.
 */

#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

#include "doctest.h"
#include "line/lang/lqn/lqn_reader.h"
#include "line/solvers/ln/lqn_analyzers.h"
#include "line/solvers/ln/solver_ln.h"

using namespace line;

namespace {

std::string model_path() {
    return std::string(LINE_MP_REPO_ROOT) +
           "/jar/src/test/resources/lqn/randomLQN/model_C1_L2_T2_P2_c4_z1_s1_t1_p1_y1_4.lqnx";
}

lqn::LqnStruct<double> model() { return lqn::read_lqnx<double>(model_path()); }

std::size_t idx_of(const lqn::LqnStruct<double>& l, const std::string& hashname) {
    for (std::size_t i = 1; i <= l.nidx; ++i)
        if (l.hashnames[i] == hashname) return i;
    return 0;
}

}  // namespace

TEST_CASE("LN methods: the default fixed point still matches MATLAB") {
    // ONE TABLE PER ENCODING, and the default is asserted against ITS OWN.
    // The 'srvn' alias resolves to 'srvn.ph' on this model, and the two
    // encodings are DIFFERENT fixed points -- 0.34144 against 0.32651 at the
    // reference task, some 4.4% apart -- so a single table cannot serve both.
    // The 'srvn.cs' numbers below are the ones this case carried while the
    // alias still resolved there; they are kept, under the method that
    // actually produces them, so the alias change rebased the golden without
    // discarding it.
    //
    // MATLAB SolverLN AvgTable on the identical file, RE-RECORDED 2026-08-12
    // with options.method set explicitly to each encoding.
    struct Row { const char* hn; double util; double tput; double respt; };
    const Row ph_rows[] = {
        {"P:p0", 0.254920484015, 0.0, 0.0},
        {"P:p1", 0.696196664539, 0.0, 0.0},
        {"E:c0", 0.682096310233, 0.326511847659, 12.2507039956},
        {"E:e0", 0.254920484015, 0.418068975043, 2.27502447305},
        {"E:e1", 0.696196664539, 0.308957989384, 2.25337},
    };
    const Row cs_rows[] = {
        {"P:p0", 0.266578956773, 0.0, 0.0},
        {"P:p1", 0.728036353951, 0.0, 0.0},
        {"E:c0", 0.713291045538, 0.34144441731, 11.7149374451},
        {"E:e0", 0.266578956773, 0.437188842068, 2.27502447044},
        {"E:e1", 0.728036353951, 0.323087799141, 2.25337},
    };

    const lqn::LqnStruct<double> l = model();
    for (int enc = 0; enc < 2; ++enc) {
        const bool ph = (enc == 0);
        ln::LnOptions opt;
        // The first pass leaves the method at its default ON PURPOSE: what is
        // under test there is where the alias LANDS, not the encoding named.
        if (!ph) opt.method = "srvn.cs";
        CAPTURE(ph ? "default (srvn.ph)" : "srvn.cs");
        ln::SolverLN<double> s(l, opt);
        const ln::LnSolution<double> sol = s.get_ensemble_avg();
        REQUIRE(sol.converged);
        CHECK(!sol.is_bound);

        const Row* rows = ph ? ph_rows : cs_rows;
        for (std::size_t k = 0; k < 5; ++k) {
            const Row& row = rows[k];
            const std::size_t i = idx_of(l, row.hn);
            REQUIRE(i > 0);
            CAPTURE(std::string(row.hn));
            if (row.util > 1e-7) CHECK(sol.UN[i] == doctest::Approx(row.util).epsilon(1e-3));
            if (row.tput > 1e-7) CHECK(sol.TN[i] == doctest::Approx(row.tput).epsilon(1e-3));
            if (row.respt > 1e-7) CHECK(sol.RN[i] == doctest::Approx(row.respt).epsilon(1e-3));
        }
    }
}

TEST_CASE("LN methods: moment3 fits every entry and keeps the throughputs") {
    const lqn::LqnStruct<double> l = model();
    ln::LnOptions def;
    // moment3 builds the ROUTING layers and adds a distribution pass to them, so
    // the baseline it must not move is that same encoding. The 'srvn' alias
    // resolves to 'srvn.ph' on this model, which is a different fixed point.
    def.method = "srvn.cs";
    ln::SolverLN<double> sd(l, def);
    const ln::LnSolution<double> base = sd.get_ensemble_avg();

    ln::LnOptions opt;
    opt.method = "moment3";
    ln::SolverLN<double> s(l, opt);
    const ln::LnSolution<double> sol = s.get_ensemble_avg();
    REQUIRE(sol.converged);

    // The fixed point is driven by the same mean update in both methods -- the
    // APH fitting only replaces the ENTRY service law -- so the throughputs and
    // the host utilizations must be identical, not merely close. A moment3 run
    // that moved them would mean the distribution pass had leaked into the
    // iteration.
    const char* same[] = {"P:p0", "P:p1", "E:c0", "E:e0", "E:e1"};
    for (const char* hn : same) {
        const std::size_t i = idx_of(l, hn);
        REQUIRE(i > 0);
        CAPTURE(std::string(hn));
        CHECK(sol.TN[i] == doctest::Approx(base.TN[i]).epsilon(1e-9));
        CHECK(sol.UN[i] == doctest::Approx(base.UN[i]).epsilon(1e-9));
    }

    // THE ENTRY SERVICE IS REFITTED, so it moves -- but only by the difference
    // between an exponential and an APH with the same first moment. These are
    // the JAR's numbers for the same model and method (jline.solvers.ln.SolverLN
    // with options.method = "moment3", 2026-07-31); MATLAB does not finish this
    // method on this model in reasonable time, so the JAR is the reference here.
    //
    // THE TOLERANCE IS WHAT IT IS BECAUSE THE GRIDS DIFFER. Both codebases refine
    // the passage-time grid until no CDF jump exceeds 5e-4, but the JAR starts
    // from its ODE solver's own output steps and this port from a uniform grid,
    // so the quadrature that reads the moments off the curve is not the same
    // sum. An entry left at zero, or one per cent off, is a defect; a per-mille
    // difference is the grid.
    struct M3 { const char* hn; double respt; };
    const M3 m3[] = {{"E:c0", 11.6795}, {"E:e0", 2.27477}, {"E:e1", 2.25748}};
    for (const M3& row : m3) {
        const std::size_t i = idx_of(l, row.hn);
        REQUIRE(i > 0);
        CAPTURE(std::string(row.hn));
        CHECK(sol.RN[i] > 0.0);
        CHECK(sol.RN[i] == doctest::Approx(row.respt).epsilon(5e-3));
    }
}

TEST_CASE("LN methods: getCdfRespT refuses the phase-type encoding") {
    // The distribution pass reads the routing encoding, so a solver built for
    // srvn.ph has nothing to run it over and says so, rather than reporting an
    // empty table. Same refusal as @SolverLN/getCdfRespT.
    const lqn::LqnStruct<double> l = model();
    ln::LnOptions ph;
    ph.method = "srvn.ph";
    ln::SolverLN<double> s(l, ph);
    CHECK_THROWS_AS(ln::lqn_cdf_respt(s), UnsupportedError);
}

TEST_CASE("LN methods: getCdfRespT returns a distribution per entry") {
    const lqn::LqnStruct<double> l = model();
    ln::LnOptions opt;
    // the routing encoding, which the distribution pass needs; the method itself
    // is left unset, so getCdfRespT is the one that switches to moment3
    opt.method = "srvn.cs";
    ln::SolverLN<double> s(l, opt);
    const std::vector<ln::LnCdf> cdf = ln::lqn_cdf_respt(s);
    REQUIRE(cdf.size() > l.nentries);

    for (std::size_t e = 1; e <= l.nentries; ++e) {
        CAPTURE(e);
        const ln::LnCdf& c = cdf[e];
        REQUIRE(!c.empty());
        REQUIRE(c.t.size() == c.cdf.size());
        CHECK(c.t.front() == doctest::Approx(0.0));
        // A CDF: nondecreasing in t, inside [0,1], reaching the upper tail.
        for (std::size_t k = 0; k + 1 < c.t.size(); ++k) {
            REQUIRE(c.t[k + 1] >= c.t[k]);
            REQUIRE(c.cdf[k + 1] >= c.cdf[k] - 1e-12);
            REQUIRE(c.cdf[k] >= -1e-12);
            REQUIRE(c.cdf[k] <= 1.0 + 1e-12);
        }
        CHECK(c.cdf.back() > 0.99);
    }

    // The fitted law's mean is the entry service the same run reports, which is
    // the consistency the two passes of moment3 are supposed to have: the
    // distribution pass writes servt from the SAME convolution it evaluates.
    ln::LnOptions m3;
    m3.method = "moment3";
    ln::SolverLN<double> s3(l, m3);
    const ln::LnSolution<double> sol = s3.get_ensemble_avg();
    for (std::size_t e = 1; e <= l.nentries; ++e) {
        const ln::LnCdf& c = cdf[e];
        double mean = 0.0;
        for (std::size_t k = 0; k + 1 < c.t.size(); ++k)
            mean += 0.5 * (c.t[k + 1] + c.t[k]) * (c.cdf[k + 1] - c.cdf[k]);
        CAPTURE(e);
        CHECK(mean == doctest::Approx(sol.RN[l.eshift + e]).epsilon(0.02));
    }
}

TEST_CASE("LN methods: the box bounds bracket the fixed point") {
    const lqn::LqnStruct<double> l = model();
    ln::LnOptions def;
    ln::SolverLN<double> sd(l, def);
    const ln::LnSolution<double> base = sd.get_ensemble_avg();

    ln::LnOptions up;
    up.method = "mwba.upper";
    ln::SolverLN<double> su(l, up);
    const ln::LnSolution<double> upper = su.get_ensemble_avg();

    ln::LnOptions lo;
    lo.method = "mwba.lower";
    ln::SolverLN<double> sl(l, lo);
    const ln::LnSolution<double> lower = sl.get_ensemble_avg();

    CHECK(upper.is_bound);
    CHECK(lower.is_bound);

    // A BOUND DEFINES THROUGHPUT AND PROCESSOR UTILIZATION AND NOTHING ELSE.
    // Reporting a queue length there would be a claim the bound does not make,
    // so the undefined flags are part of the contract, not cosmetics.
    for (std::size_t i = 1; i <= l.nidx; ++i) {
        CHECK(!upper.defined_Q[i]);
        CHECK(!upper.defined_R[i]);
        CHECK(!lower.defined_Q[i]);
        CHECK(!lower.defined_R[i]);
    }

    const char* tput_of[] = {"R:c0", "T:t0", "T:t1", "E:c0", "E:e0", "E:e1"};
    for (const char* hn : tput_of) {
        const std::size_t i = idx_of(l, hn);
        REQUIRE(i > 0);
        CAPTURE(std::string(hn));
        REQUIRE(upper.defined_T[i]);
        REQUIRE(lower.defined_T[i]);
        CHECK(lower.TN[i] <= base.TN[i] + 1e-9);
        CHECK(upper.TN[i] >= base.TN[i] - 1e-9);
    }

    // The processor utilization is bracketed too, and the upper bound saturates
    // the bottleneck host: p1 carries the heaviest demand in this model.
    const std::size_t p1 = idx_of(l, "P:p1");
    REQUIRE(p1 > 0);
    CHECK(upper.UN[p1] == doctest::Approx(1.0).epsilon(1e-9));
    CHECK(lower.UN[p1] <= base.UN[p1] + 1e-9);
}

TEST_CASE("LN methods: the analytic sensitivity agrees with central differences") {
    // THE AGREEMENT IS A PRODUCT-FORM CLAIM, so it is asserted on the encoding
    // whose layers are product-form. 'srvn.cs' charges each call as its own
    // exponential class, so a layer is exactly the recursion pfqn_sens
    // differentiates and the two branches coincide to eight digits. Under
    // 'srvn.ph', which is where the 'srvn' alias lands on this model, the entry
    // service is a COMPOSED PHASE-TYPE law: the layer solve reads its higher
    // moments, the analytic branch differentiates a product-form recursion that
    // reads only the mean, and the two legitimately part company. MATLAB splits
    // exactly the same way (2026-08-12), so the divergence is pinned below
    // rather than tolerated by a wider epsilon.
    const lqn::LqnStruct<double> l = model();
    ln::LnOptions opt;
    opt.method = "srvn.cs";

    ln::SolverLN<double> se(l, opt);
    sens::SensOptions eo;
    eo.method = "exact";
    const ln::LnSensTable<double> exact = ln::lqn_sensitivity_table(se, eo);

    ln::SolverLN<double> sf(l, opt);
    sens::SensOptions fo;
    fo.method = "fd";
    fo.scheme = "central";
    const ln::LnSensTable<double> fd = ln::lqn_sensitivity_table(sf, fo);

    CHECK(exact.method == "exact");
    CHECK(fd.method == "fd");
    REQUIRE(exact.rows.size() == fd.rows.size());
    REQUIRE(exact.rows.size() > 0);

    for (std::size_t r = 0; r < exact.rows.size(); ++r) {
        const auto& a = exact.rows[r];
        const auto& b = fd.rows[r];
        CAPTURE(a.layer);
        CAPTURE(a.station);
        CAPTURE(a.jobclass);
        CHECK(a.station == b.station);
        CHECK(a.jobclass == b.jobclass);
        CHECK(a.dTput == doctest::Approx(b.dTput).epsilon(1e-3));
        CHECK(a.dRespT == doctest::Approx(b.dRespT).epsilon(1e-3));
        CHECK(a.dQLen == doctest::Approx(b.dQLen).epsilon(1e-3));
        CHECK(a.dUtil == doctest::Approx(b.dUtil).epsilon(1e-3));

        // SIGNS ARE THE POINT of the table, and they are not free: a faster
        // server can only raise the throughput and can only lower the response
        // time, the queue length and the utilization it charges. An unsigned
        // clean-up on the way out printed exact zeros for the last three, which
        // is why this is asserted here rather than left to the numbers.
        CHECK(a.dTput >= -1e-9);
        CHECK(a.dRespT <= 1e-9);
        CHECK(a.dQLen <= 1e-9);
        CHECK(a.dUtil <= 1e-9);
    }
}

TEST_CASE("LN methods: under srvn.ph the two sensitivity branches part company") {
    // The companion of the case above, and the reason it names its encoding.
    // Under the default alias the composed phase-type entry law puts the task
    // layer outside the product-form recursion, so `exact` and `fd` answer
    // different questions at T:t0 -- and only there, the host layers being
    // ordinary either way. Both branches are pinned to MATLAB SolverLN
    // getSensitivityTable on the identical file (2026-08-12), so a port that
    // silently moved either one is caught even though they disagree.
    const lqn::LqnStruct<double> l = model();
    ln::LnOptions opt;  // default: the 'srvn' alias resolves to 'srvn.ph' here

    ln::SolverLN<double> se(l, opt);
    sens::SensOptions eo;
    eo.method = "exact";
    const ln::LnSensTable<double> exact = ln::lqn_sensitivity_table(se, eo);

    ln::SolverLN<double> sf(l, opt);
    sens::SensOptions fo;
    fo.method = "fd";
    fo.scheme = "central";
    const ln::LnSensTable<double> fd = ln::lqn_sensitivity_table(sf, fo);

    struct Row { const char* station; const char* jobclass;
                 double dT, dR, dQ, dU; };
    const Row exact_rows[] = {
        {"P:p0", "T:t0", 0.064984453, -0.3718036, -0.11581482, -0.11581482},
        {"P:p1", "T:t1", 0.4846898, -5.0776764, -0.47660323, -0.47660323},
        {"T:t0", "R:c0", 0.97691746, -33.517979, -2.0408197, -0.051553249},
        {"T:t1", "T:t0", 0.4846898, -2.7731158, -0.35221551, -0.35221551},
    };
    const Row fd_rows[] = {
        {"P:p0", "T:t0", 0.064984454, -0.3718036, -0.11581483, -0.11581483},
        {"P:p1", "T:t1", 0.4846898, -5.0776764, -0.47660323, -0.47660323},
        {"T:t0", "R:c0", 0.89818466, -33.69983, -1.8763414, -0.15419044},
        {"T:t1", "T:t0", 0.4846898, -2.7665339, -0.34946383, -0.35221551},
    };

    auto check = [](const ln::LnSensTable<double>& t, const Row* want) {
        REQUIRE(t.rows.size() == 4u);
        for (std::size_t r = 0; r < 4; ++r) {
            const auto& a = t.rows[r];
            CAPTURE(a.station);
            CAPTURE(a.jobclass);
            CHECK(a.station == want[r].station);
            CHECK(a.jobclass == want[r].jobclass);
            CHECK(a.dTput == doctest::Approx(want[r].dT).epsilon(1e-5));
            CHECK(a.dRespT == doctest::Approx(want[r].dR).epsilon(1e-5));
            CHECK(a.dQLen == doctest::Approx(want[r].dQ).epsilon(1e-5));
            CHECK(a.dUtil == doctest::Approx(want[r].dU).epsilon(1e-5));
            // the signs hold on both branches, product form or not
            CHECK(a.dTput >= -1e-9);
            CHECK(a.dRespT <= 1e-9);
            CHECK(a.dQLen <= 1e-9);
            CHECK(a.dUtil <= 1e-9);
        }
    };
    CHECK(exact.method == "exact");
    CHECK(fd.method == "fd");
    check(exact, exact_rows);
    check(fd, fd_rows);

    // and the gap is at the composed layer alone: the two host layers agree
    // there as tightly as they do under 'srvn.cs'.
    for (std::size_t r = 0; r < exact.rows.size(); ++r) {
        if (exact.rows[r].station[0] != 'P') continue;
        CAPTURE(exact.rows[r].station);
        CHECK(exact.rows[r].dTput == doctest::Approx(fd.rows[r].dTput).epsilon(1e-6));
        CHECK(exact.rows[r].dQLen == doctest::Approx(fd.rows[r].dQLen).epsilon(1e-6));
    }
}

TEST_CASE("LN methods: NC layers reproduce the MVA fixed point") {
    const lqn::LqnStruct<double> l = model();
    // The encoding is named because the claim below is an EXACTNESS claim, and
    // it is recorded for the routing layers. Under 'srvn.ph', which is what the
    // 'srvn' alias resolves to on this model, the two layer solvers stop the
    // outer Picard iteration at iterates ~3e-5 apart -- well inside its own
    // iter_tol of 5e-3, so nothing is wrong there, but far outside the 1e-6 that
    // says the two algorithms evaluated the SAME fixed point.
    ln::LnOptions mva;
    mva.method = "srvn.cs";
    ln::SolverLN<double> sm(l, mva);
    const ln::LnSolution<double> a = sm.get_ensemble_avg();

    ln::LnOptions nc;
    nc.method = "srvn.cs";
    nc.layer_solver = "nc";
    ln::SolverLN<double> sn(l, nc);
    const ln::LnSolution<double> b = sn.get_ensemble_avg();
    REQUIRE(b.converged);

    // Every layer here is product-form, so the normalizing constant and MVA are
    // the SAME algorithm evaluated differently: the fixed points must coincide
    // to solver tolerance, not merely be close.
    for (std::size_t i = 1; i <= l.nidx; ++i) {
        CAPTURE(l.hashnames[i]);
        if (a.defined_T[i] && b.defined_T[i])
            CHECK(b.TN[i] == doctest::Approx(a.TN[i]).epsilon(1e-6));
        if (a.defined_U[i] && b.defined_U[i])
            CHECK(b.UN[i] == doctest::Approx(a.UN[i]).epsilon(1e-6));
        if (a.defined_R[i] && b.defined_R[i])
            CHECK(b.RN[i] == doctest::Approx(a.RN[i]).epsilon(1e-6));
    }
}

TEST_CASE("LN methods: the transient converges to the steady state") {
    const lqn::LqnStruct<double> l = model();
    ln::LnOptions opt;
    opt.layer_solver = "fluid";
    opt.timespan_end = 40.0;
    opt.tran_points = 81;

    ln::SolverLN<double> ss(l, opt);
    const ln::LnSolution<double> steady = ss.get_ensemble_avg();
    REQUIRE(steady.converged);

    // The DEFAULT is the coupled relaxation, as in the reference; the decoupled
    // transient has to be asked for by name.
    ln::LnOptions dop = opt;
    dop.ln_transient = "decoupled";
    ln::SolverLN<double> sd(l, dop);
    const ln::LnTranSolution dec = ln::lqn_tran_avg(sd);
    CHECK(dec.mode == "decoupled");
    CHECK(dec.iterations == 0);
    REQUIRE(dec.layers.size() > 0);

    ln::SolverLN<double> sc(l, opt);
    const ln::LnTranSolution cou = ln::lqn_tran_avg(sc);
    CHECK(cou.mode == "coupled");
    CHECK(cou.iterations >= 1);
    REQUIRE(cou.layers.size() == dec.layers.size());

    for (std::size_t e = 0; e < dec.layers.size(); ++e) {
        const ln::LnTranLayer& d = dec.layers[e];
        const ln::LnTranLayer& c = cou.layers[e];
        CAPTURE(e);
        REQUIRE(d.t.size() == std::size_t(opt.tran_points));
        REQUIRE(c.t.size() == d.t.size());
        CHECK(d.t.front() == doctest::Approx(0.0));
        CHECK(d.t.back() == doctest::Approx(opt.timespan_end));
        for (std::size_t k = 0; k + 1 < d.t.size(); ++k) REQUIRE(d.t[k + 1] > d.t[k]);

        // Waveform relaxation reconciles the layers in model time; on a model
        // this loosely coupled it barely moves the trajectories, so the two
        // modes must agree at the horizon. A coupled run that diverged from the
        // decoupled one here would mean the rate schedule was mis-scaled.
        for (std::size_t i = 0; i < d.QN.size(); ++i)
            for (std::size_t r = 0; r < d.QN[i].size(); ++r) {
                const double qd = d.QN[i][r].back(), qc = c.QN[i][r].back();
                REQUIRE(qd >= -1e-9);
                CHECK(qc == doctest::Approx(qd).epsilon(0.02).scale(1.0));
            }
    }
}

TEST_CASE("LN methods: the refusals name the reason") {
    const lqn::LqnStruct<double> l = model();

    // The transient is a FLUID quantity: an MVA layer has no trajectory to
    // report, and answering with the steady state repeated would be a fiction.
    ln::LnOptions mva;
    mva.timespan_end = 10.0;
    ln::SolverLN<double> sm(l, mva);
    CHECK_THROWS_AS(ln::lqn_tran_avg(sm), UnsupportedError);

    // An unknown mode is rejected by name rather than silently defaulted.
    ln::LnOptions bad;
    bad.layer_solver = "fluid";
    bad.timespan_end = 10.0;
    bad.ln_transient = "sequential";
    ln::SolverLN<double> sb(l, bad);
    CHECK_THROWS_AS(ln::lqn_tran_avg(sb), InputError);
}

/**
 * An unset horizon is a REQUEST, not an omission. There is then no common
 * interval to co-evolve over, so the COUPLED default defers to the decoupled
 * path and each layer picks its own horizon by the fluid analyzer's rule --
 * thirty mean events of its slowest transition. That is the reference's rule
 * exactly: `@SolverLN/getTranAvgCoupled.m` returns `getTranAvgDecoupled()` when
 * `options.timespan` is not two finite endpoints. The grids are therefore
 * per-layer and need not share an endpoint, so what is asserted here is that
 * every one of them is a usable time axis rather than that they agree.
 */
TEST_CASE("LN methods: an unset horizon defers to the per-layer choice") {
    const lqn::LqnStruct<double> l = model();
    ln::LnOptions nospan;
    nospan.layer_solver = "fluid";
    nospan.tran_points = 41;  // ln_transient left at the coupled default
    ln::SolverLN<double> sn(l, nospan);

    ln::LnTranSolution tr;
    REQUIRE_NOTHROW(tr = ln::lqn_tran_avg(sn));
    // the coupled default stepped aside rather than co-evolving over nothing
    CHECK(tr.mode == "decoupled");
    CHECK(tr.iterations == 0);
    REQUIRE(tr.layers.size() > 0);
    for (std::size_t e = 0; e < tr.layers.size(); ++e) {
        const ln::LnTranLayer& L = tr.layers[e];
        CAPTURE(e);
        REQUIRE(L.t.size() == std::size_t(nospan.tran_points));
        CHECK(L.t.front() == doctest::Approx(0.0));
        CHECK(std::isfinite(L.t.back()));
        CHECK(L.t.back() > 0.0);
        for (std::size_t k = 0; k + 1 < L.t.size(); ++k) REQUIRE(L.t[k + 1] > L.t[k]);
        // and the layer actually reported on that axis
        REQUIRE(L.QN.size() > 0);
        for (std::size_t i = 0; i < L.QN.size(); ++i)
            for (std::size_t r = 0; r < L.QN[i].size(); ++r) {
                REQUIRE(L.QN[i][r].size() == L.t.size());
                CHECK(L.QN[i][r].back() >= -1e-9);
            }
    }
}
