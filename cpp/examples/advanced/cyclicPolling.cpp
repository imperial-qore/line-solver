/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * `python/examples/advanced/cyclicPolling/`: a single server visiting one
 * buffer per class in turn, under the four service disciplines.
 *
 * EACH CLASS IS A QUEUE OF THE POLLING MODEL, which is what makes a
 * Source/Queue/Sink model with a POLLING station a cyclic polling system:
 * `solver_mva_polling_analyzer` reads the per-class arrival rate, the two
 * service moments and the switchover moments off exactly this shape.
 *
 * JMT is refused by name: this port does not carry it.
 */

#include <cstdio>
#include <string>

#include "example_util.h"
#include "examples_common.h"
#include "line/solvers/mva/solver_mva_runner.h"

namespace line {
namespace examples {

namespace {

/** mySource -> myQueue (POLLING) -> mySink, with one buffer per class. */
Net polling_model(const std::string& nm, const D& arrival1, const D& service1, const D& arrival2,
                  const D& service2, std::size_t& queue, std::size_t& c1, std::size_t& c2) {
    Net m(nm);
    Source source(m, "mySource");
    queue = m.add_queue("myQueue", SchedStrategy::POLLING);
    Sink sink(m, "mySink");
    c1 = m.add_open_class("myClass1");
    source.set_arrival(c1, arrival1);
    m.set_service(queue, c1, service1);
    c2 = m.add_open_class("myClass2");
    source.set_arrival(c2, arrival2);
    m.set_service(queue, c2, service2);
    Routing P;
    serial(P, c1, {source, queue, sink});
    serial(P, c2, {source, queue, sink});
    m.link(P);
    return m;
}

/**
 * The MVA table, which is what each reference script prints first.
 *
 * A REFUSAL IS REPORTED, NOT RETHROWN. `solver_mva_polling_analyzer` refuses a
 * polling system whose total switchover time is zero, which is what an absent
 * or Immediate switchover reduces to here; the reference reaches the same
 * formula with a switchover of mean 1/GlobalConstants.Immediate instead and
 * returns finite waiting times, so this is a divergence to read rather than a
 * number to invent.
 */
void polling_mva(Net& m) {
    const Sn& sn = m.get_struct();
    section("MVA");
    const mva::MvaOptions opt;
    const Matrix<double> init;
    try {
        print_avg(sn, mva::solver_mva_run_analyzer(sn, opt, init));
    } catch (const Error& e) {
        std::printf("N/A: %s\n", e.what());
    }
}

/** The JMT cross-check every one of these examples ends on. */
void polling_jmt(Net& m, std::size_t samples) {
    section("JMT");
    print_avg(m.get_struct(), jmt_avg(m, sim_opts(23000, samples)));
}

}  // namespace

/** `polling_exhaustive_exp.py`: exhaustive service, exponential everywhere. */
void polling_exhaustive_exp() {
    std::size_t queue = 0, c1 = 0, c2 = 0;
    Net m = polling_model("M[2]/M[2]/1-Gated", Exp(0.1), Exp(1.0),
                          Exp(0.1), Exp(1.5), queue, c1, c2);
    m.set_polling_type(queue, lang::PollingType::EXHAUSTIVE);
    polling_mva(m);
    // TODO(cpp): avg_table_jmt = JMT(model, seed=23000).getAvgTable(); print(avg_table_jmt)
    polling_jmt(m, 0);
}

/** `polling_exhaustive_det.py`: the same, with deterministic times and immediate walks. */
void polling_exhaustive_det() {
    std::size_t queue = 0, c1 = 0, c2 = 0;
    Net m = polling_model("M[2]/M[2]/1-Gated", Det(1.0), Det(0.001), Det(1.0),
                          Det(0.001), queue, c1, c2);
    m.set_polling_type(queue, lang::PollingType::EXHAUSTIVE);
    m.set_switchover(queue, c1, Immediate());
    m.set_switchover(queue, c2, Immediate());
    polling_mva(m);
    // TODO(cpp): avg_table_jmt = JMT(model, seed=23000, samples=int(1e5)).getAvgTable(); print(avg_table_jmt)
    polling_jmt(m, 100000);
}

/** `polling_gated.py`: gated service with exponential switchover times. */
void polling_gated() {
    std::size_t queue = 0, c1 = 0, c2 = 0;
    Net m = polling_model("M[2]/M[2]/1-Exhaustive", Exp(1.0), Exp(4.0),
                          Exp(0.8), Exp(1.5), queue, c1, c2);
    m.set_polling_type(queue, lang::PollingType::GATED);
    m.set_switchover(queue, c1, Exp(1.0));
    m.set_switchover(queue, c2, Exp(0.5));
    polling_mva(m);
    // TODO(cpp): avg_table_jmt = JMT(model, seed=23000, samples=int(1e6)).getAvgTable(); print(avg_table_jmt)
    polling_jmt(m, 1000000);
}

/** `polling_klimited.py`: at most one job served per visit. */
void polling_klimited() {
    std::size_t queue = 0, c1 = 0, c2 = 0;
    Net m = polling_model("M[2]/M[2]/1-Gated", Exp(0.2), Exp(1.0),
                          Exp(0.3), Exp(1.5), queue, c1, c2);
    m.set_polling_type(queue, lang::PollingType::KLIMITED, 1);
    m.set_switchover(queue, c1, Exp(1.0));
    m.set_switchover(queue, c2, Immediate());
    polling_mva(m);
    // TODO(cpp): avg_table_jmt = JMT(model, samples=int(1e5), seed=23000).getAvgTable(); print(avg_table_jmt)
    polling_jmt(m, 100000);
}

LINE_EXAMPLE("advanced/cyclicPolling", polling_exhaustive_exp);
LINE_EXAMPLE("advanced/cyclicPolling", polling_exhaustive_det);
LINE_EXAMPLE("advanced/cyclicPolling", polling_gated);
LINE_EXAMPLE("advanced/cyclicPolling", polling_klimited);

}  // namespace examples
}  // namespace line
