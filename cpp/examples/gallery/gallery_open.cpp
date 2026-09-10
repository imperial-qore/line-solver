/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * The open-queue half of the gallery: `gallery_mm1` through `gallery_um1`,
 * plus the multiclass, feedback, reentrant and tandem variants.
 *
 * Each factory is its reference file line for line. The one thing that is NOT
 * a transcription is the pair of random MAPs: `MAP.rand(seed=23000)` draws from
 * Python's Mersenne Twister, which no C++ RNG reproduces, so the REALIZED
 * matrices are written out below exactly as that seed produces them. Anything
 * else would give this port a different model under the same name.
 */

#include <cstdio>
#include <string>
#include <vector>

#include "gallery.h"
#include "line/api/mam/map_transform.h"
#include "line/lang/dist_fitters.h"

namespace line {
namespace examples {

namespace {

/**
 * `MAP.rand(seed=23000)`, the two-phase MAP every `gallery_map*` entry uses.
 *
 * `gallery_mapm1.py` already carries these sixteen digits inline; the other
 * three entries call `MAP.rand(seed=23000)` and get the same pair, which was
 * verified by drawing it rather than assumed.
 */
mam::Map<double> map_rand_23000() {
    mam::Map<double> m;
    m.D0 = Matrix<double>(2, 2, 0.0);
    m.D1 = Matrix<double>(2, 2, 0.0);
    m.D0(0, 0) = -0.6984901916396979;
    m.D0(0, 1) = 0.45234650636128054;
    m.D0(1, 0) = 0.34690024319398277;
    m.D0(1, 1) = -0.8194057961021199;
    m.D1(0, 0) = 0.2125067546435463;
    m.D1(0, 1) = 0.033636930634871165;
    m.D1(1, 0) = 0.4441520099524867;
    m.D1(1, 1) = 0.028353542955650513;
    return m;
}

/** `MAP.rand(n=3, seed=23001)`, the second service process of mmap1_multiclass. */
mam::Map<double> map_rand_23001_n3() {
    mam::Map<double> m;
    m.D0 = Matrix<double>(3, 3, 0.0);
    m.D1 = Matrix<double>(3, 3, 0.0);
    m.D0(0, 0) = -1.4381234478114808;
    m.D0(0, 1) = 0.23527923807749257;
    m.D0(0, 2) = 0.12505764264346458;
    m.D0(1, 0) = 0.16432447428913083;
    m.D0(1, 1) = -1.342134785997235;
    m.D0(1, 2) = 0.23494442318500375;
    m.D0(2, 0) = 0.39249255613161854;
    m.D0(2, 1) = 0.03744558449493818;
    m.D0(2, 2) = -1.6697294651684789;
    m.D1(0, 0) = 0.35840271271024693;
    m.D1(0, 1) = 0.25319705703370393;
    m.D1(0, 2) = 0.46618679734657287;
    m.D1(1, 0) = 0.3901108570600689;
    m.D1(1, 1) = 0.2793353004054196;
    m.D1(1, 2) = 0.27341973105761186;
    m.D1(2, 0) = 0.3740965715993369;
    m.D1(2, 1) = 0.46805435432660947;
    m.D1(2, 2) = 0.3976403986159759;
    return m;
}

D map_of(const mam::Map<double>& m) {
    return D::map_dist(m.D0, m.D1, lang::ProcessType::MAP);
}

/**
 * The service means of the `*_linear` family: `linspace(0.1, Umax, n/2)`
 * mirrored, with Umax in the middle when n is odd, and MATLAB's
 * `linspace(a, b, 1) == b` reproduced.
 */
std::vector<double> linear_means(std::size_t n, double umax) {
    const std::size_t half = n / 2;
    std::vector<double> means;
    if (half == 1) {
        means.push_back(umax);
    } else {
        for (std::size_t i = 0; i < half; ++i)
            means.push_back(0.1 + (umax - 0.1) * static_cast<double>(i) /
                                      static_cast<double>(half - 1));
    }
    std::vector<double> out = means;
    if (n % 2 != 0) out.push_back(umax);
    for (std::size_t i = means.size(); i-- > 0;) out.push_back(means[i]);
    return out;
}

/** Source -> Queue -> Sink, the shape most of the gallery is. */
Net one_queue(const std::string& model_name, const std::string& src, const std::string& q,
              const std::string& snk, SchedStrategy sched, const D& arrival, const D& service,
              double nservers = 1.0) {
    Net m(model_name);
    Source source(m, src);
    Queue queue(m, q, sched);
    Sink sink(m, snk);
    OpenClass oclass(m, "myClass");
    source.set_arrival(oclass, arrival);
    queue.set_service(oclass, service);
    if (nservers != 1.0) queue.set_number_of_servers(nservers);
    Routing P;
    serial(P, {source, queue, sink});
    m.link(P);
    return m;
}

/**
 * The reentrant shape: Class1 arrives, is served, switches to Class2 at the
 * same queue and only then leaves. Class2 has no arrivals.
 */
Net reentrant(const std::string& model_name, SchedStrategy sched, const D& arrival,
              const D& service1, const D& service2, double nservers = 1.0) {
    Net m(model_name);
    Source source(m, "Source");
    Queue queue(m, "Queue", sched);
    Sink sink(m, "Sink");
    if (nservers != 1.0) queue.set_number_of_servers(nservers);
    OpenClass c1(m, "Class1");
    OpenClass c2(m, "Class2");
    source.set_arrival(c1, arrival);
    source.set_arrival(c2, Disabled());
    queue.set_service(c1, service1);
    queue.set_service(c2, service2);
    Routing P;
    P.set(c1, c1, source, queue, 1.0);
    P.set(c1, c2, queue, queue, 1.0);
    P.set(c2, c2, queue, sink, 1.0);
    m.link(P);
    return m;
}

}  // namespace

// ---------------------------------------------------------------------------
// Single-class open queues
// ---------------------------------------------------------------------------

Net gallery_mm1() {
    return one_queue("M/M/1", "mySource", "myQueue", "mySink", SchedStrategy::FCFS,
                     Exp(1.0), Exp(2.0));
}

Net gallery_mm1_ps() {
    return one_queue("M/M/1-PS", "mySource", "myQueue", "mySink", SchedStrategy::PS,
                     Exp(1.0), Exp(2.0));
}

Net gallery_mmk(double k) {
    return one_queue("M/M/k", "mySource", "myQueue", "mySink", SchedStrategy::FCFS,
                     Exp(1.0), Exp(2.0), k);
}

Net gallery_mm1k(double K) {
    Net m("M/M/1/K");
    Source source(m, "Source");
    Queue queue(m, "Queue", SchedStrategy::FCFS);
    Sink sink(m, "Sink");
    queue.set_number_of_servers(1.0);
    queue.set_capacity(K);
    OpenClass oclass(m, "Class1", 0);
    source.set_arrival(oclass, Exp(0.8));
    queue.set_service(oclass, Exp(1.0));
    Routing P;
    P.set(oclass, oclass, source, queue, 1.0);
    P.set(oclass, oclass, queue, sink, 1.0);
    m.link(P);
    return m;
}

Net gallery_mdk(double k) {
    return one_queue("M/D/k", "mySource", "myQueue", "mySink", SchedStrategy::FCFS,
                     D::exp_mean(1.0), Det(2.0 / k), k);
}

Net gallery_merl1() {
    return one_queue("M/E/1", "mySource", "myQueue", "mySink", SchedStrategy::FCFS,
                     Exp(1.0), lang::erlang_fit_mean_order<double>(0.5, 2));
}

Net gallery_merlk(double k) {
    return one_queue("M/Erl/" + std::to_string(static_cast<int>(k)), "mySource", "myQueue",
                     "mySink", SchedStrategy::FCFS, Exp(1.0),
                     lang::erlang_fit_mean_order<double>(0.5, 2), k);
}

Net gallery_mhyp1() {
    return one_queue("M/H/1", "mySource", "myQueue", "mySink", SchedStrategy::FCFS,
                     Exp(1.0), lang::hyperexp_fit_mean_scv<double>(0.5, 4.0));
}

Net gallery_mhypk(double k) {
    return one_queue("M/Hyper/" + std::to_string(static_cast<int>(k)), "mySource", "myQueue",
                     "mySink", SchedStrategy::FCFS, Exp(1.0),
                     lang::coxian_fit_mean_scv<double>(0.5, 4.0), k);
}

Net gallery_mpar1() {
    Net m("M/Par/1");
    Source source(m, "Source");
    Queue queue(m, "Queue", SchedStrategy::FCFS);
    Sink sink(m, "Sink");
    OpenClass c1(m, "Class1");
    source.set_arrival(c1, Exp(1.0));
    queue.set_service(c1, lang::pareto_fit_mean_scv<double>(0.5, 64.0));
    Routing P;
    P.set(c1, c1, source, queue, 1.0);
    P.set(c1, c1, queue, sink, 1.0);
    m.link(P);
    return m;
}

Net gallery_erlm1() {
    return one_queue("Er/M/1", "mySource", "myQueue", "mySink", SchedStrategy::FCFS,
                     lang::erlang_fit_mean_order<double>(1.0, 5), Exp(2.0));
}

Net gallery_erlm1_ps() {
    return one_queue("Er/M/1-PS", "mySource", "myQueue", "mySink", SchedStrategy::PS,
                     lang::erlang_fit_mean_order<double>(1.0, 5), Exp(2.0));
}

Net gallery_erlm1ps() { return gallery_erlm1_ps(); }

Net gallery_erldk(double k) {
    return one_queue("Erl/D/k", "mySource", "myQueue", "mySink", SchedStrategy::FCFS,
                     lang::erlang_fit_mean_order<double>(1.0, 5), Det(2.0 / k), k);
}

Net gallery_hyperlk(double k) {
    return one_queue("Hyper/Erl/k", "mySource", "myQueue", "mySink", SchedStrategy::FCFS,
                     lang::hyperexp_fit_mean_scv_balanced<double>(1.0 / 1.8, 4.0),
                     D::erlang_fit(1.0, 0.25), k);
}

Net gallery_hypm1() {
    return one_queue("H2/M/1", "mySource", "myQueue", "mySink", SchedStrategy::FCFS,
                     lang::hyperexp_fit_mean_scv<double>(1.0, 64.0), Exp(2.0));
}

Net gallery_detm1() {
    return one_queue("D/M/1", "mySource", "myQueue", "mySink", SchedStrategy::FCFS, Det(1.0),
                     Exp(2.0));
}

Net gallery_dm1() {
    Net m("D/M/1");
    Source source(m, "Source");
    Queue queue(m, "Queue", SchedStrategy::FCFS);
    Sink sink(m, "Sink");
    OpenClass c1(m, "Class1");
    source.set_arrival(c1, Det(1.0));
    queue.set_service(c1, Exp(2.0));
    Routing P;
    serial(P, {source, queue, sink});
    m.link(P);
    return m;
}

Net gallery_gamm1() {
    Net m("Gam/M/1");
    Source source(m, "Source");
    Queue queue(m, "Queue", SchedStrategy::FCFS);
    Sink sink(m, "Sink");
    OpenClass c1(m, "Class1");
    source.set_arrival(c1, lang::gamma_fit_mean_scv<double>(1.0, 1.0 / 5.0));
    queue.set_service(c1, Exp(2.0));
    Routing P;
    serial(P, {source, queue, sink});
    m.link(P);
    return m;
}

Net gallery_parm1() {
    Net m("Par/M/1");
    Source source(m, "Source");
    Queue queue(m, "Queue", SchedStrategy::FCFS);
    Sink sink(m, "Sink");
    OpenClass c1(m, "Class1");
    source.set_arrival(c1, lang::pareto_fit_mean_scv<double>(1.0, 64.0));
    queue.set_service(c1, Exp(2.0));
    Routing P;
    P.set(c1, c1, source, queue, 1.0);
    P.set(c1, c1, queue, sink, 1.0);
    m.link(P);
    return m;
}

Net gallery_um1() {
    Net m("U/M/1");
    Source source(m, "Source");
    Queue queue(m, "Queue", SchedStrategy::FCFS);
    Sink sink(m, "Sink");
    OpenClass c1(m, "Class1");
    source.set_arrival(c1, Uniform(1.0, 2.0));
    queue.set_service(c1, Exp(2.0));
    Routing P;
    P.set(c1, c1, source, queue, 1.0);
    P.set(c1, c1, queue, sink, 1.0);
    m.link(P);
    return m;
}

Net gallery_aphm1() {
    return one_queue("APH/M/1", "mySource", "myQueue", "mySink", SchedStrategy::FCFS,
                     lang::aph_fit_central<double>(1.0, 0.99, 1.999), Exp(2.0));
}

Net gallery_coxm1() {
    return one_queue("Cox/M/1", "mySource", "myQueue", "mySink", SchedStrategy::FCFS,
                     lang::coxian_fit_central<double>(1.0, 0.99, 1.999), Exp(2.0));
}

Net gallery_mapm1() {
    return one_queue("MAP/M/1", "mySource", "myQueue", "mySink", SchedStrategy::FCFS,
                     map_of(map_rand_23000()), Exp(2.0));
}

Net gallery_mapmk(double k) {
    return one_queue("MAP/M/k", "mySource", "myQueue", "mySink", SchedStrategy::FCFS,
                     map_of(map_rand_23000()), Exp(2.0), k);
}

Net gallery_mmap1() {
    return one_queue("M/MAP/1", "mySource", "myQueue", "mySink", SchedStrategy::FCFS,
                     Exp(1.0), map_of(mam::map_scale(map_rand_23000(), 0.5)));
}

Net gallery_mmapk(double k) {
    return one_queue("M/MAP/k", "mySource", "myQueue", "mySink", SchedStrategy::FCFS,
                     Exp(1.0), map_of(map_rand_23000()), k);
}

/**
 * `Replayer(example_trace.txt)`: the trace-driven arrival process.
 *
 * The file is the reference's own, read from the repository rather than
 * copied, so the two codebases replay the same 10000 samples.
 */
Net gallery_replayerm1() {
    const std::vector<double> trace =
        read_trace(std::string(LINE_EXAMPLES_REPO_ROOT) +
                   "/python/examples/gettingstarted/example_trace.txt");
    const D replayer = D::replayer(trace);
    return one_queue("Trace/M/1", "mySource", "myQueue", "mySink", SchedStrategy::FCFS, replayer,
                     Exp(3.0 / replayer.mean));
}

// ---------------------------------------------------------------------------
// Multiclass, feedback, reentrant, tandem
// ---------------------------------------------------------------------------

Net gallery_mm1_multiclass() {
    Net m("M[2]/M[2]/1");
    Source source(m, "mySource");
    Queue queue(m, "myQueue", SchedStrategy::FCFS);
    Sink sink(m, "mySink");
    OpenClass c1(m, "myClass1");
    source.set_arrival(c1, Exp(1.0));
    queue.set_service(c1, Exp(4.0));
    OpenClass c2(m, "myClass2");
    source.set_arrival(c2, Exp(0.5));
    queue.set_service(c2, Exp(4.0));
    Routing P;
    serial(P, c1, {source, queue, sink});
    serial(P, c2, {source, queue, sink});
    m.link(P);
    return m;
}

Net gallery_mm1_ps_multiclass() {
    Net m("M[2]/M[2]/1-PS");
    Source source(m, "mySource");
    Queue queue(m, "myQueue", SchedStrategy::PS);
    Sink sink(m, "mySink");
    OpenClass c1(m, "myClass1");
    source.set_arrival(c1, Exp(1.0));
    queue.set_service(c1, Exp(4.0));
    OpenClass c2(m, "myClass2");
    source.set_arrival(c2, Exp(0.5));
    queue.set_service(c2, Exp(4.0));
    Routing P;
    serial(P, c1, {source, queue, sink});
    serial(P, c2, {source, queue, sink});
    m.link(P);
    return m;
}

/** HOL at the queue: myClass1 carries priority 1, myClass2 priority 0. */
Net gallery_mm1_prio() {
    Net m("M[2]/M[2]/1");
    Source source(m, "mySource");
    Queue queue(m, "myQueue", SchedStrategy::HOL);
    Sink sink(m, "mySink");
    OpenClass c1(m, "myClass1", 1);
    source.set_arrival(c1, Exp(1.0));
    queue.set_service(c1, Exp(4.0));
    OpenClass c2(m, "myClass2", 0);
    source.set_arrival(c2, Exp(0.5));
    queue.set_service(c2, Exp(4.0));
    Routing P;
    serial(P, c1, {source, queue, sink});
    serial(P, c2, {source, queue, sink});
    m.link(P);
    return m;
}

Net gallery_mm1_feedback() {
    const double p = 1.0 / 3.0;
    Net m("M/M/1-Feedback");
    Source source(m, "Source");
    Queue queue(m, "Queue", SchedStrategy::FCFS);
    Sink sink(m, "Sink");
    OpenClass c1(m, "Class1");
    source.set_arrival(c1, D::exp_mean(1.0));
    queue.set_service(c1, D::exp_mean(0.5));
    Routing P;
    P.set(c1, c1, source, queue, 1.0);
    P.set(c1, c1, queue, queue, p);
    P.set(c1, c1, queue, sink, 1.0 - p);
    m.link(P);
    return m;
}

Net gallery_mm1_ps_feedback() {
    const double p = 1.0 / 3.0;
    Net m("M/M/1-PS-Feedback");
    Source source(m, "Source");
    Queue queue(m, "Queue", SchedStrategy::PS);
    Sink sink(m, "Sink");
    OpenClass c1(m, "Class1");
    source.set_arrival(c1, D::exp_mean(1.0));
    queue.set_service(c1, D::exp_mean(0.5));
    Routing P;
    P.set(c1, c1, source, queue, 1.0);
    P.set(c1, c1, queue, queue, p);
    P.set(c1, c1, queue, sink, 1.0 - p);
    m.link(P);
    return m;
}

Net gallery_mm1_reentrant() {
    return reentrant("M/M/1-Reentrant", SchedStrategy::FCFS, Exp(1.0), Exp(2.0),
                     Exp(3.0));
}

Net gallery_mm1_ps_reentrant() {
    return reentrant("M/M/1", SchedStrategy::PS, Exp(1.0), Exp(2.0),
                     Exp(3.0));
}

Net gallery_hypm1_reentrant() {
    return reentrant("Hyper/M/1-Reentrant", SchedStrategy::FCFS,
                     lang::hyperexp_fit_mean_scv<double>(1.0, 4.0), Exp(2.0),
                     Exp(3.0));
}

Net gallery_merl1_reentrant() {
    return reentrant("M/Erl/1-Reentrant", SchedStrategy::FCFS, Exp(1.0),
                     lang::erlang_fit_mean_order<double>(0.5, 5), Exp(3.0));
}

Net gallery_erlm1_reentrant() {
    return reentrant("Er/M/1-Reentrant", SchedStrategy::FCFS,
                     lang::erlang_fit_mean_order<double>(1.0, 5), Exp(2.0),
                     Exp(3.0));
}

Net gallery_mhyp1_reentrant() {
    return reentrant("M/Hyper/1-Reentrant", SchedStrategy::FCFS, Exp(1.0),
                     lang::coxian_fit_mean_scv<double>(0.5, 4.0), Exp(3.0));
}

Net gallery_hyphyp1_reentrant() {
    return reentrant("Hyper/Hyper/1-Reentrant", SchedStrategy::FCFS,
                     lang::hyperexp_fit_mean_scv<double>(1.0, 64.0),
                     lang::hyperexp_fit_mean_scv<double>(0.5, 4.0), Exp(3.0), 2.0);
}

Net gallery_hyperl1_reentrant() {
    return reentrant("Hyper/Erl/1-Reentrant", SchedStrategy::FCFS,
                     lang::hyperexp_fit_mean_scv<double>(1.0, 64.0),
                     lang::erlang_fit_mean_order<double>(0.5, 5), Exp(3.0));
}

/** Class1 feeds back into itself 90% of the time; there is no class switch. */
Net gallery_hyperl1_feedback() {
    Net m("Hyper/Erl/1-Feedback");
    Source source(m, "Source");
    Queue queue(m, "Queue", SchedStrategy::FCFS);
    Sink sink(m, "Sink");
    OpenClass c1(m, "Class1");
    source.set_arrival(c1, lang::hyperexp_fit_mean_scv<double>(1.0, 64.0));
    queue.set_service(c1, lang::erlang_fit_mean_order<double>(0.05, 5));
    Routing P;
    P.set(c1, c1, source, queue, 1.0);
    P.set(c1, c1, queue, queue, 0.9);
    P.set(c1, c1, queue, sink, 0.1);
    m.link(P);
    return m;
}

/** Erl/Erl/1: Class2 exists only to carry the exit, and never arrives. */
Net gallery_erlerl1(std::size_t n) {
    Net m("Erl/Erl/1");
    Source source(m, "Source");
    Queue queue(m, "Queue", SchedStrategy::FCFS);
    Sink sink(m, "Sink");
    OpenClass c1(m, "Class1");
    OpenClass c2(m, "Class2");
    source.set_arrival(c1, lang::erlang_fit_mean_order<double>(1.0, n));
    source.set_arrival(c2, Disabled());
    queue.set_service(c1, lang::erlang_fit_mean_order<double>(0.5, n));
    queue.set_service(c2, Exp(3.0));
    Routing P;
    P.set(c1, c1, source, queue, 1.0);
    P.set(c2, c2, queue, sink, 1.0);
    m.link(P);
    return m;
}

/** The reentrant variant, whose class switch and exit each carry probability 1/2. */
Net gallery_erlerl1_reentrant() {
    const std::size_t n = 5;
    Net m("Erl/Erl/1-Reentrant");
    Source source(m, "Source");
    Queue queue(m, "Queue", SchedStrategy::FCFS);
    Sink sink(m, "Sink");
    OpenClass c1(m, "Class1");
    OpenClass c2(m, "Class2");
    source.set_arrival(c1, lang::erlang_fit_mean_order<double>(1.0, n));
    source.set_arrival(c2, Disabled());
    queue.set_service(c1, lang::erlang_fit_mean_order<double>(0.1, n));
    queue.set_service(c2, Exp(10.0));
    Routing P;
    P.set(c1, c1, source, queue, 1.0);
    P.set(c1, c2, queue, queue, 0.5);
    P.set(c2, c2, queue, sink, 0.5);
    m.link(P);
    return m;
}

namespace {

/** The `*_linear` chain: Source, n FCFS queues, Sink, one class through all. */
Net linear_chain(const std::string& model_name, std::size_t n, const D& arrival,
                 const std::vector<D>& services) {
    Net m(model_name);
    std::vector<std::size_t> line;
    line.push_back(m.add_source("mySource"));
    for (std::size_t i = 0; i < n; ++i)
        line.push_back(m.add_queue("Queue" + std::to_string(i + 1), SchedStrategy::FCFS));
    line.push_back(m.add_sink("mySink"));
    OpenClass oclass(m, "myClass");
    m.set_arrival(line[0], oclass, arrival);
    for (std::size_t i = 0; i < n; ++i) m.set_service(line[i + 1], oclass, services[i]);
    Routing P;
    serial(P, line);
    m.link(P);
    return m;
}

}  // namespace

Net gallery_mm1_linear(std::size_t n, double umax) {
    const std::vector<double> means = linear_means(n, umax);
    std::vector<D> svc;
    for (std::size_t i = 0; i < n; ++i) svc.push_back(D::exp_mean(means[i]));
    return linear_chain("M/M/1-Linear", n, Exp(1.0), svc);
}

Net gallery_mm1_tandem() { return gallery_mm1_linear(2); }

Net gallery_merl1_linear(std::size_t n, double umax) {
    const std::vector<double> means = linear_means(n, umax);
    std::vector<D> svc;
    for (std::size_t i = 0; i < n; ++i)
        svc.push_back(lang::erlang_fit_mean_order<double>(means[i], i + 1));
    return linear_chain("M/Erl/1-Linear", n, Exp(1.0), svc);
}

Net gallery_merl1_tandem() { return gallery_merl1_linear(2); }

Net gallery_mhyp1_linear(std::size_t n, double umax) {
    const std::vector<double> means = linear_means(n, umax);
    std::vector<D> svc;
    for (std::size_t i = 0; i < n; ++i)
        svc.push_back(lang::hyperexp_fit_mean_scv<double>(means[i], static_cast<double>(n)));
    return linear_chain("M/Hyp/1-Linear", n, Exp(1.0), svc);
}

Net gallery_mhyp1_tandem() { return gallery_mhyp1_linear(2); }

Net gallery_hyphyp1_linear(std::size_t n, double umax) {
    const std::vector<double> means = linear_means(n, umax);
    std::vector<D> svc;
    for (std::size_t i = 0; i < n; ++i)
        svc.push_back(
            lang::hyperexp_fit_mean_scv<double>(means[i], 1.0 + static_cast<double>(i) + 1.0));
    return linear_chain("Hyp/Hyp/1-Linear", n, lang::hyperexp_fit_mean_scv<double>(1.0, 2.0), svc);
}

Net gallery_hyphyp1_tandem() { return gallery_hyphyp1_linear(2); }

Net gallery_mm1_tandem_multiclass() {
    Net m("M[2]/M[2]/1 -> -/M[2]/1");
    Source source(m, "Source");
    Queue q1(m, "Queue1", SchedStrategy::FCFS);
    Queue q2(m, "Queue2", SchedStrategy::FCFS);
    Sink sink(m, "mySink");
    OpenClass c1(m, "myClass1");
    source.set_arrival(c1, Exp(1.0));
    q1.set_service(c1, Exp(4.0));
    q2.set_service(c1, Exp(6.0));
    OpenClass c2(m, "myClass2");
    source.set_arrival(c2, Exp(0.5));
    q1.set_service(c2, Exp(2.0));
    q2.set_service(c2, Exp(6.0));
    Routing P;
    serial(P, c1, {source, q1, q2, sink});
    serial(P, c2, {source, q1, q2, sink});
    m.link(P);
    return m;
}

/** Two MAP-served classes, each arriving so that the offered load is fixed. */
Net gallery_mmap1_multiclass() {
    const mam::Map<double> m1 = mam::map_scale(map_rand_23000(), 0.5);
    const mam::Map<double> m2 = mam::map_scale(map_rand_23001_n3(), 0.5);
    Net m("M/MAP/1");
    Source source(m, "mySource");
    Queue queue(m, "myQueue", SchedStrategy::FCFS);
    Sink sink(m, "mySink");
    OpenClass c1(m, "myClass1");
    source.set_arrival(c1, Exp(0.35 / mam::map_mean(m1)));
    queue.set_service(c1, map_of(m1));
    OpenClass c2(m, "myClass2");
    source.set_arrival(c2, Exp(0.15 / mam::map_mean(m2)));
    queue.set_service(c2, map_of(m2));
    Routing P;
    serial(P, c1, {source, queue, sink});
    serial(P, c2, {source, queue, sink});
    m.link(P);
    return m;
}

/**
 * The Lu-Kumar / Rybko-Stolyar reentrant network: two chains crossing two
 * stations in opposite order, with the FBFS priority that makes it unstable
 * under HOL despite both utilizations sitting at 0.88.
 */
Net gallery_lukumar_reentrant(const std::string& sched) {
    SchedStrategy s = SchedStrategy::FCFS;
    if (sched == "HOL" || sched == "hol") s = SchedStrategy::HOL;
    else if (sched == "PS" || sched == "ps") s = SchedStrategy::PS;

    Net m("Lu-Kumar-Reentrant");
    Source source(m, "Source");
    Queue st1(m, "Station1", s);
    Queue st2(m, "Station2", s);
    Sink sink(m, "Sink");

    OpenClass c1(m, "Class1", 1);
    OpenClass c2(m, "Class2", 0);
    OpenClass c3(m, "Class3", 1);
    OpenClass c4(m, "Class4", 0);

    const double lambda = 0.08;
    source.set_arrival(c1, Exp(lambda));
    source.set_arrival(c2, Disabled());
    source.set_arrival(c3, Exp(lambda));
    source.set_arrival(c4, Disabled());

    // Slow first visits (mean 10), fast second visits (mean 1): the asymmetry
    // that creates the virtual bottleneck.
    st1.set_service(c1, Exp(1.0 / 10.0));
    st1.set_service(c2, Disabled());
    st1.set_service(c3, Disabled());
    st1.set_service(c4, Exp(1.0 / 1.0));
    st2.set_service(c1, Disabled());
    st2.set_service(c2, Exp(1.0 / 1.0));
    st2.set_service(c3, Exp(1.0 / 10.0));
    st2.set_service(c4, Disabled());

    Routing P;
    P.set(c1, c1, source, st1, 1.0);
    P.set(c1, c2, st1, st2, 1.0);
    P.set(c2, c2, st2, sink, 1.0);
    P.set(c3, c3, source, st2, 1.0);
    P.set(c3, c4, st2, st1, 1.0);
    P.set(c4, c4, st1, sink, 1.0);
    m.link(P);
    return m;
}

}  // namespace examples
}  // namespace line
