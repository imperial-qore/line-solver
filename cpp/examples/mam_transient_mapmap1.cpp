/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * `python/examples/mam_transient_mapmap1.py`: the transient MAP/MAP/1 queue.
 *
 * The reference builds a three-phase MAP arrival and a two-phase MAP service,
 * rescales the service so that rho = 0.6, and asks SolverMAM for the transient
 * mean queue length, utilization and throughput over [0, 40] starting empty.
 *
 * WHICH ENGINE RUNS is not chosen here and is not chosen by the method either:
 * `solver_mam_get_tran_avg` forces `ldqbd` and then lets
 * `mam_transient_qbd_applicable` decide, which sends a correlated MAP arrival
 * with a correlated MAP service to the Laplace-domain transient QBD. That is
 * the same dispatch the Python `getTranAvg` docstring describes, so the two
 * codebases run the same algorithm on the same model rather than agreeing by
 * coincidence.
 *
 * The rho = 0.6 rescaling is done here exactly as the reference does it --
 * elementwise multiplication of BOTH service blocks by a scalar -- and not
 * through `map_scale`, which normalizes to a target MEAN and would land on a
 * different pair of matrices.
 */

#include <cstdio>
#include <vector>

#include "examples_common.h"
#include "line/api/mam/map_moment.h"
#include "line/solvers/mam/solver_mam_runner.h"

namespace line {
namespace examples {

namespace {

/** The reference's `D0`, `D1`, `S0` and `S1`, entry for entry. */
mam::Map<double> arrival_map() {
    mam::Map<double> m;
    m.D0 = Matrix<double>(3, 3, 0.0);
    m.D1 = Matrix<double>(3, 3, 0.0);
    m.D0(0, 0) = -8.0;
    m.D0(0, 1) = 1.0;
    m.D0(0, 2) = 3.0;
    m.D0(1, 1) = -6.0;
    m.D0(1, 2) = 4.0;
    m.D0(2, 0) = 2.0;
    m.D0(2, 2) = -3.0;
    m.D1(0, 0) = 3.0;
    m.D1(0, 1) = 1.0;
    m.D1(1, 1) = 2.0;
    m.D1(2, 2) = 1.0;
    return m;
}

mam::Map<double> service_map() {
    mam::Map<double> m;
    m.D0 = Matrix<double>(2, 2, 0.0);
    m.D1 = Matrix<double>(2, 2, 0.0);
    m.D0(0, 0) = -3.0;
    m.D0(0, 1) = 1.0;
    m.D0(1, 0) = 6.0;
    m.D0(1, 1) = -7.0;
    m.D1(0, 1) = 2.0;
    m.D1(1, 0) = 1.0;
    return m;
}

}  // namespace

void mam_transient_mapmap1() {
    const mam::Map<double> a = arrival_map();
    mam::Map<double> s = service_map();

    // Scale service to rho = 0.6: both blocks by one factor, as the reference.
    const double fac = (mam::map_lambda(a) / 0.6) / mam::map_lambda(s);
    for (std::size_t i = 0; i < 2; ++i)
        for (std::size_t j = 0; j < 2; ++j) {
            s.D0(i, j) *= fac;
            s.D1(i, j) *= fac;
        }

    Net m("MAP/MAP/1 transient");
    Source source(m, "Source");
    Queue queue(m, "Queue", SchedStrategy::FCFS);
    Sink sink(m, "Sink");
    OpenClass oclass(m, "Class1");
    source.set_arrival(oclass, D::map_dist(a.D0, a.D1, lang::ProcessType::MAP));
    queue.set_service(oclass, D::map_dist(s.D0, s.D1, lang::ProcessType::MAP));
    Routing P;
    serial(P, {source, queue, sink});
    m.link(P);

    mam::MamOptions opt;
    opt.timespan_start = 0.0;
    opt.timespan_end = 40.0;
    const mam::TranResult<double> tr = mam::solver_mam_get_tran_avg(m.get_struct(), opt);

    // Station 1 is the Queue; the Source carries no transient curve.
    const mam::TranCurve<double>& q = tr.Qt[1][0];
    const mam::TranCurve<double>& u = tr.Ut[1][0];
    const mam::TranCurve<double>& t = tr.Tt[1][0];
    std::printf("MAP/MAP/1 transient (rho=0.6), start empty:\n");
    std::printf("  t=%5.1f  E[N]=%.5f  U=%.5f  Tput=%.5f\n", q.times.back(), q.values.back(),
                u.values.back(), t.values.back());
    std::printf("  steady-state E[N] approaches 1.5458 as t -> inf.\n");
}

LINE_EXAMPLE(".", mam_transient_mapmap1);

}  // namespace examples
}  // namespace line
