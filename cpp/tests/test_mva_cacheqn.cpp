/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * Integrated caching-queueing analyzer (WS-E, dispatch branch 9). The model is
 * gallery_cache_routing: Source -> Cache(LRU, 4 items, cap 2) -> HitQueue /
 * MissQueue -> Sink, three open classes (Init reads, Hit/Miss are the switched
 * classes). Reference numbers are MATLAB SolverMVA(model).getAvgTable.
 */

#include "doctest.h"
#include "line/lang/qn/network_builder.h"
#include "line/solvers/mva/solver_mva_runner.h"

using namespace line;
using lang::Distrib;
using lang::ReplacementStrategy;
using lang::SchedStrategy;

namespace {

using D = Distrib<double>;

std::size_t station_of(const qn::NetworkStruct<double>& sn, const std::string& nm) {
    for (std::size_t i = 0; i < sn.nstations; ++i)
        if (sn.stations[i].name == nm) return i;
    FAIL("no station named ", nm);
    return 0;
}

TEST_CASE("cacheqn: an embedded cache routing hit/miss to distinct queues") {
    qn::Network<double> m("Cache-Routing");
    const std::size_t src = m.add_source("Source");
    qn::CacheParam<double> ch;
    ch.nitems = 4;
    ch.itemcap = std::vector<int>{2};
    ch.replacestrat = ReplacementStrategy::LRU;
    // InitClass reads uniformly; Hit/Miss classes do not read (empty rows).
    ch.pread = std::vector<std::vector<double> >{{0.25, 0.25, 0.25, 0.25}, {}, {}};
    ch.hitclass = std::vector<std::size_t>{2, 0, 0};   // InitClass -> HitClass (class 2)
    ch.missclass = std::vector<std::size_t>{3, 0, 0};  //           -> MissClass (class 3)
    const std::size_t cache = m.add_cache("Cache", ch);
    const std::size_t hq = m.add_queue("HitQueue", SchedStrategy::FCFS);
    const std::size_t mq = m.add_queue("MissQueue", SchedStrategy::FCFS);
    const std::size_t snk = m.add_sink("Sink");
    const std::size_t init = m.add_open_class("InitClass");
    const std::size_t hit = m.add_open_class("HitClass");
    const std::size_t miss = m.add_open_class("MissClass");
    m.set_arrival(src, init, D::exp_rate(1.0));
    m.set_service(hq, hit, D::exp_rate(2.0));
    m.set_service(mq, miss, D::exp_rate(1.0));
    qn::RoutingMatrix<double> P;
    P.set(init, init, src, cache, 1.0);
    P.set(hit, hit, cache, hq, 1.0);
    P.set(hit, hit, hq, snk, 1.0);
    P.set(miss, miss, cache, mq, 1.0);
    P.set(miss, miss, mq, snk, 1.0);
    m.link(P);

    mva::MvaOptions opt;
    Matrix<double> init_sol;
    const mva::AvgResult<double> r = mva::solver_mva_run_analyzer(m.get_struct(), opt, init_sol);
    const qn::NetworkStruct<double>& sn = m.get_struct();
    const std::size_t H = station_of(sn, "HitQueue"), M = station_of(sn, "MissQueue");

    // HitQueue serves HitClass (class index 1); MissQueue serves MissClass (2).
    // All six columns of the MATLAB getAvgTable are checked.
    CHECK(r.QN(H, 1) == doctest::Approx(0.14286).epsilon(1e-4));
    CHECK(r.UN(H, 1) == doctest::Approx(0.125).epsilon(1e-4));
    CHECK(r.RN(H, 1) == doctest::Approx(0.57143).epsilon(1e-4));
    CHECK(r.WN(H, 1) == doctest::Approx(0.28571).epsilon(1e-4));  // ResidT
    CHECK(r.AN(H, 1) == doctest::Approx(0.5).epsilon(1e-4));      // ArvR (offered)
    CHECK(r.TN(H, 1) == doctest::Approx(0.25).epsilon(1e-4));
    CHECK(r.QN(M, 2) == doctest::Approx(0.33333).epsilon(1e-4));
    CHECK(r.UN(M, 2) == doctest::Approx(0.25).epsilon(1e-4));
    CHECK(r.RN(M, 2) == doctest::Approx(1.3333).epsilon(1e-4));
    CHECK(r.WN(M, 2) == doctest::Approx(0.66667).epsilon(1e-4));
    CHECK(r.AN(M, 2) == doctest::Approx(0.5).epsilon(1e-4));
    CHECK(r.TN(M, 2) == doctest::Approx(0.25).epsilon(1e-4));
    // The over-routing pass-through: each queue carries the other class at the
    // hit/miss rate with no service (QLen 0). The SERVED metrics above match
    // MATLAB exactly (QLen 0.14286 / 0.33333, served Tput 0.25); this pass-through
    // Tput is where the C++ da_cacheqn diverges: it reports 0.25, while MATLAB's
    // un-normalized visit propagation reports 0.125 (and Source InitClass Tput
    // 0.75 vs the C++ 1.0). Reproducing MATLAB's exact secondary metrics needs
    // its reducible-over-route visit arithmetic, which native Python also fails
    // to match (it concentrates the flow instead). See _kb/07 cacheqn divergence.
    CHECK(r.TN(H, 2) == doctest::Approx(0.25).epsilon(1e-4));  // MissClass through HitQueue
    CHECK(r.TN(M, 1) == doctest::Approx(0.25).epsilon(1e-4));  // HitClass through MissQueue
    CHECK(r.QN(H, 2) == doctest::Approx(0.0).epsilon(1e-9));
    CHECK(r.QN(M, 1) == doctest::Approx(0.0).epsilon(1e-9));

    // The cache split rides with the result on this branch too. THE ORACLE IS
    // SYMMETRY, not another solver: the reads are uniform over 4 items and the
    // single list holds 2, so every item is resident with probability exactly
    // 1/2 and the hit fraction is exactly 1/2, whatever the analyzer does
    // internally. `mva_dispatch` used to drop all of this.
    REQUIRE(r.cache.caches.size() == 1);
    const solvers::CacheNodeMetrics<double>& cm = r.cache.caches[0];
    CHECK(cm.name == "Cache");
    CHECK(cm.nitems == 4);
    CHECK(cm.hitprob[0] == doctest::Approx(0.5).epsilon(1e-6));
    CHECK(cm.missprob[0] == doctest::Approx(0.5).epsilon(1e-6));
    // Per-item occupancy, (nitems x lists+1) with column 0 the miss. LRU takes
    // the characteristic-time branch, so this is a genuine per-item law.
    REQUIRE(cm.itemprob.rows() == 4);
    REQUIRE(cm.itemprob.cols() == 2);
    for (std::size_t i = 0; i < 4; ++i) {
        CHECK(cm.itemprob(i, 0) + cm.itemprob(i, 1) == doctest::Approx(1.0).epsilon(1e-6));
        // symmetry again: no item is preferred
        CHECK(cm.itemprob(i, 1) == doctest::Approx(0.5).epsilon(1e-6));
    }
}

}  // namespace
