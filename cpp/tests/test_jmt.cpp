/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * SolverJMT: the JSIM and JMVA writers, and an end-to-end solve when a JVM and
 * common/JMT.jar are both present.
 *
 * THE WRITER TESTS DO NOT NEED JMT. They assert on the document, which is where
 * a port defect actually lives: a wrong `classPath`, a dropped `refClass`, a
 * priority written without the LINE-to-JMT inversion. An end-to-end run would
 * catch those only through a wrong number, and only on a host with the jar.
 *
 * THE END-TO-END CASES SKIP, they do not fail, when JMT is absent: the jar is a
 * 32 MB artifact that is not in the repository, and a red suite on a machine
 * that simply has not downloaded it says nothing about the port.
 */

#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "doctest.h"
#include "line/io/jmva_writer.h"
#include "line/lang/qn/network_builder.h"
#include "line/solvers/wrappers/jmt/jmt_logs.h"
#include "line/solvers/wrappers/jmt/solver_jmt.h"
#include "line/util/tempdir.h"

using namespace line;
using lang::Distrib;
using lang::SchedStrategy;

namespace {

using D = Distrib<double>;

/** Is `needle` a substring of `hay`? */
bool has(const std::string& hay, const std::string& needle) {
    return hay.find(needle) != std::string::npos;
}

/** How many times `needle` occurs in `hay`. */
std::size_t count(const std::string& hay, const std::string& needle) {
    std::size_t n = 0, p = 0;
    while ((p = hay.find(needle, p)) != std::string::npos) {
        ++n;
        p += needle.size();
    }
    return n;
}

/** M/M/1: Source -> Queue -> Sink, arrival 0.5, service 1.0 (rho = 0.5). */
qn::Network<double> mm1(double arrival_rate = 0.5, double service_rate = 1.0) {
    qn::Network<double> net("mm1");
    const std::size_t src = net.add_source("Source");
    const std::size_t q = net.add_queue("Queue1", SchedStrategy::FCFS);
    const std::size_t snk = net.add_sink("Sink");
    net.add_open_class("Class1");
    net.set_arrival(src, 1, D::exp_rate(arrival_rate));
    net.set_service(q, 1, D::exp_rate(service_rate));
    qn::RoutingMatrix<double> P;
    P.set(src, q, 1.0);
    P.set(q, snk, 1.0);
    net.link(P);
    return net;
}

/** A closed two-station cycle: a Delay think station and an FCFS queue. */
qn::Network<double> cqn(double n = 3.0) {
    qn::Network<double> net("cqn");
    const std::size_t d = net.add_delay("Delay");
    const std::size_t q = net.add_queue("Queue1", SchedStrategy::PS);
    net.add_closed_class("Class1", n, d);
    net.set_service(d, 1, D::exp_mean(1.0));
    net.set_service(q, 1, D::exp_mean(0.5));
    qn::RoutingMatrix<double> P;
    P.set(d, q, 1.0);
    P.set(q, d, 1.0);
    net.link(P);
    return net;
}

/** Write a JMVA document to a scratch file and read it back as a string. */
std::string jmva_text(const qn::NetworkStruct<double>& sn, const std::string& method) {
    line::util::TempDir dir("jmvatest");
    const std::string path = dir.path() + "/model.jmva";
    io::write_jmva(sn, path, method, 10000);
    std::ifstream f(path.c_str());
    std::ostringstream ss;
    ss << f.rdbuf();
    return ss.str();
}

}  // namespace

TEST_CASE("jmt writer: the JSIM header and the open-class declaration") {
    qn::Network<double> net = mm1();
    io::JmtWriteOptions opt;
    opt.seed = 4242;
    opt.max_samples = 20000;
    const std::string xml = io::jmt_write_jsim(net.get_struct(), opt);

    CHECK(has(xml, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>"));
    CHECK(has(xml, "name=\"model.jsimg\""));
    CHECK(has(xml, "xsi:noNamespaceSchemaLocation=\"SIMmodeldefinition.xsd\""));
    CHECK(has(xml, "maxSamples=\"20000\""));
    CHECK(has(xml, "seed=\"4242\""));
    CHECK(has(xml, "maxEvents=\"-1\""));
    // An unbounded horizon must NOT emit maxSimulated at all; emitting it as
    // "inf" would make JMT stop at parse time.
    CHECK_FALSE(has(xml, "maxSimulated"));
    CHECK(has(xml, "<userClass name=\"Class1\" type=\"open\""));
    CHECK(has(xml, "referenceSource=\"Source\""));
    // A closed model's `customers` attribute must be absent on an open class.
    CHECK_FALSE(has(xml, "customers="));
}

TEST_CASE("jmt writer: the sections of a Source, a Queue and a Sink") {
    const std::string xml = io::jmt_write_jsim(mm1().get_struct(), io::JmtWriteOptions());

    // Source: RandomSource / ServiceTunnel / Router, in that order.
    CHECK(has(xml, "<section className=\"RandomSource\">"));
    CHECK(has(xml, "<section className=\"ServiceTunnel\"/>"));
    CHECK(has(xml, "<section className=\"JobSink\"/>"));
    // A LINE Buffer is JMT's Queue and a LINE Dispatcher is JMT's Router.
    CHECK(has(xml, "<section className=\"Queue\">"));
    CHECK(has(xml, "<section className=\"Router\">"));
    CHECK(has(xml, "<section className=\"Server\">"));
    // The exponential service of the queue, and the arrival at the source.
    CHECK(has(xml, "classPath=\"jmt.engine.random.Exponential\" name=\"Exponential\""));
    CHECK(count(xml, "name=\"lambda\"") == 2);
    // FCFS is a tail insertion; there is no scheduling field in JMT.
    CHECK(has(xml, "QueuePutStrategies.TailStrategy"));
    CHECK(has(xml, "QueueGetStrategies.FCFSstrategy"));
    // The topology.
    CHECK(has(xml, "<connection source=\"Source\" target=\"Queue1\"/>"));
    CHECK(has(xml, "<connection source=\"Queue1\" target=\"Sink\"/>"));
}

TEST_CASE("jmt writer: a closed model declares its population and its preload") {
    const std::string xml = io::jmt_write_jsim(cqn(3.0).get_struct(), io::JmtWriteOptions());
    CHECK(has(xml, "type=\"closed\""));
    CHECK(has(xml, "customers=\"3\""));
    CHECK(has(xml, "referenceSource=\"Delay\""));
    // The whole population starts at the reference station, and the other
    // station is declared present with zero.
    CHECK(has(xml, "<preload>"));
    CHECK(has(xml, "<stationPopulations stationName=\"Delay\">"));
    CHECK(has(xml, "<classPopulation population=\"3\" refClass=\"Class1\"/>"));
    CHECK(has(xml, "<stationPopulations stationName=\"Queue1\">"));
    CHECK(has(xml, "<classPopulation population=\"0\" refClass=\"Class1\"/>"));
    // A PS station is JMT's PSServer with an EPS sharing strategy, and an
    // infinite server is JMT's Delay.
    CHECK(has(xml, "<section className=\"PSServer\">"));
    CHECK(has(xml, "PSStrategies.EPSStrategy"));
    CHECK(has(xml, "<section className=\"Delay\">"));
}

TEST_CASE("jmt writer: the measures requested for each station and class") {
    const std::string xml = io::jmt_write_jsim(mm1().get_struct(), io::JmtWriteOptions());
    // The Source reports throughput and arrival rate but no queue length,
    // response time or utilization; the Queue reports all five.
    CHECK(count(xml, "type=\"Number of Customers\"") == 1);
    CHECK(count(xml, "type=\"Utilization\"") == 1);
    CHECK(count(xml, "type=\"Response Time\"") == 1);
    CHECK(count(xml, "type=\"Throughput\"") == 2);
    CHECK(count(xml, "type=\"Arrival Rate\"") == 2);
    // Residence time is never requested: JMT's definition disagrees with LINE's
    // on class switching, and it is recomputed instead.
    CHECK_FALSE(has(xml, "type=\"Residence Time\""));
    CHECK(has(xml, "alpha=\"0.01\""));
    CHECK(has(xml, "precision=\"0.03\""));
}

TEST_CASE("jmt writer: LINE priorities are inverted for JMT") {
    qn::Network<double> net("prio");
    const std::size_t src = net.add_source("Source");
    const std::size_t q = net.add_queue("Queue1", SchedStrategy::HOL);
    const std::size_t snk = net.add_sink("Sink");
    // LINE orders priorities with the SMALLER value more urgent.
    net.add_open_class("Hi", 0);
    net.add_open_class("Lo", 4);
    net.set_arrival(src, 1, D::exp_rate(0.2));
    net.set_arrival(src, 2, D::exp_rate(0.2));
    net.set_service(q, 1, D::exp_rate(1.0));
    net.set_service(q, 2, D::exp_rate(1.0));
    qn::RoutingMatrix<double> P;
    P.set(1, 1, src, q, 1.0);
    P.set(1, 1, q, snk, 1.0);
    P.set(2, 2, src, q, 1.0);
    P.set(2, 2, q, snk, 1.0);
    net.link(P);

    const std::string xml = io::jmt_write_jsim(net.get_struct(), io::JmtWriteOptions());
    // max(prio) - prio: the LINE-urgent class 'Hi' must come out as the JMT
    // HIGH number, or the service order is reversed with no diagnostic.
    CHECK(has(xml, "<userClass name=\"Hi\" type=\"open\" softDeadline=\"0.0\" priority=\"4\""));
    CHECK(has(xml, "<userClass name=\"Lo\" type=\"open\" softDeadline=\"0.0\" priority=\"0\""));
    CHECK(has(xml, "QueuePutStrategies.TailStrategyPriority"));
}

TEST_CASE("jmt writer: a multiserver station and its buffer") {
    qn::Network<double> net = mm1();
    net.set_number_of_servers(2, 3);  // node 2 is Queue1
    const std::string xml = io::jmt_write_jsim(net.get_struct(), io::JmtWriteOptions());
    CHECK(has(xml, "name=\"maxJobs\""));
    CHECK(has(xml, "<value>3</value>"));
    // An unbounded buffer is JMT's -1, never a large finite number.
    CHECK(has(xml, "name=\"size\""));
    CHECK(has(xml, "<value>-1</value>"));
}

TEST_CASE("jmt writer: a binding closed-class station capacity is refused") {
    // BUG-81's JMT half. LINE blocks a closed job that finds no room -- the
    // upstream departure is disabled and the job stays where it is -- and no JMT
    // drop strategy reproduces that. `waiting queue`, which is what a WAITQ
    // station maps onto, does not enforce <size> at all: on this fixture JSIM
    // returned the UNCONSTRAINED [2.03 1.99 1.98], X = 0.750, against the exact
    // [3.6090 0.9711 1.4199], X = 0.6522. `BAS blocking` does enforce it, but
    // completes the service BEFORE blocking, so the blocked job moves the instant
    // room frees -- a different queueing model, [2.871 1.373 1.756], X = 0.7126.
    // Refusing by name is the same call `assert_class_cap_exportable` already
    // makes for the per-class capacity.
    qn::Network<double> net("tandem");
    const std::size_t q1 = net.add_queue("Q1", SchedStrategy::FCFS);
    const std::size_t q2 = net.add_queue("Q2", SchedStrategy::FCFS);
    const std::size_t q3 = net.add_queue("Q3", SchedStrategy::FCFS);
    const std::size_t c = net.add_closed_class("C1", 6, q1);
    net.set_service(q1, c, D::exp_rate(1.0));
    net.set_service(q2, c, D::exp_rate(1.0));
    net.set_service(q3, c, D::exp_rate(1.0));
    net.set_capacity(q2, 2.0);
    qn::RoutingMatrix<double> P;
    P.set(c, c, q1, q2, 1.0);
    P.set(c, c, q2, q3, 1.0);
    P.set(c, c, q3, q1, 1.0);
    net.link(P);
    CHECK_THROWS_AS(io::jmt_write_jsim(net.get_struct(), io::JmtWriteOptions()),
                    UnsupportedError);

    // The OPEN half must still export: JMT loses a refused open arrival exactly
    // as LINE does, and M/M/1/K is the model that says so.
    qn::Network<double> open = mm1();
    open.set_capacity(2, 3.0);  // node 2 is Queue1
    CHECK_NOTHROW(io::jmt_write_jsim(open.get_struct(), io::JmtWriteOptions()));
}

TEST_CASE("jmt writer: a DERIVED multi-class capacity is not a buffer") {
    // `refresh_capacity` derives `sn.cap` for a station nobody capped, as
    // min(sum_c chaincap, sum_r classcap) -- and the classcap row is the CHAIN
    // population repeated per class, so a station serving two classes of one
    // 4-job chain gets 8, which those 4 jobs can never reach. The reachability
    // test was `cap != total`, an EQUALITY that only the single-class case
    // satisfies, so every multi-class station fell through and was refused as a
    // "binding" buffer the model never declared. It is `<` now, and such a
    // station exports JMT's unbounded -1.
    //
    // The two classes must share ONE chain: two INDEPENDENT closed classes are
    // two chains, both chaincap columns are then set and the sums agree at the
    // population, so the inflation does not arise.
    qn::Network<double> net("twoclass");
    const std::size_t d = net.add_delay("Delay");
    const std::size_t q = net.add_queue("Queue1", SchedStrategy::PS);
    const std::size_t c1 = net.add_closed_class("Class1", 2, d);
    const std::size_t c2 = net.add_closed_class("Class2", 2, d);
    net.set_service(d, c1, D::exp_rate(1.0));
    net.set_service(d, c2, D::exp_rate(1.0));
    net.set_service(q, c1, D::exp_rate(2.0));
    net.set_service(q, c2, D::exp_rate(2.0));
    // Every (class, node) row sums to 1, and each class switches into the other
    // with probability 1/2 at both stations, which is what puts the two of them
    // in ONE chain.
    qn::RoutingMatrix<double> P;
    P.set(c1, c1, d, q, 0.5);
    P.set(c1, c2, d, q, 0.5);
    P.set(c2, c1, d, q, 0.5);
    P.set(c2, c2, d, q, 0.5);
    P.set(c1, c1, q, d, 0.5);
    P.set(c1, c2, q, d, 0.5);
    P.set(c2, c1, q, d, 0.5);
    P.set(c2, c2, q, d, 0.5);
    net.link(P);

    const qn::NetworkStruct<double>& sn = net.get_struct();
    const std::size_t iq = sn.nodes[q - 1].station;
    CHECK(sn.nchains == 1);
    CHECK(sn.cap[iq - 1] == 8.0);   // (2 classes served) x (chain population 4)

    std::string xml;
    CHECK_NOTHROW(xml = io::jmt_write_jsim(sn, io::JmtWriteOptions()));
    CHECK(has(xml, "name=\"size\""));
    CHECK_FALSE(has(xml, "<value>8</value>"));
}

TEST_CASE("jmt writer: a DECLARED capacity above the population is not a buffer") {
    // The same equality's other false positive: setCapacity(100) on a 4-job
    // closed model was refused, though 4 jobs cannot fill 100 either. Only
    // `cap < total` is a buffer.
    qn::Network<double> net("roomy");
    const std::size_t d = net.add_delay("Delay");
    const std::size_t q = net.add_queue("Queue1", SchedStrategy::FCFS);
    const std::size_t c = net.add_closed_class("C", 4, d);
    net.set_service(d, c, D::exp_rate(1.0));
    net.set_service(q, c, D::exp_rate(2.0));
    net.set_capacity(q, 100.0);
    qn::RoutingMatrix<double> P;
    P.set(c, c, d, q, 1.0);
    P.set(c, c, q, d, 1.0);
    net.link(P);
    CHECK_NOTHROW(io::jmt_write_jsim(net.get_struct(), io::JmtWriteOptions()));
}

TEST_CASE("jmt writer: the pool block follows the service strategies") {
    // JMT's SimLoader picks the Server constructor by the POSITIONAL types of the
    // parameters, so classParallelism and the four pool parameters must come after
    // ServiceStrategy: written before it, the engine dies with "Constructor of
    // Section not found" and every heterogeneous model is unsolvable.
    qn::Network<double> net = mm1();
    typename qn::Station<double>::ServerType fast;
    fast.name = "Fast";
    fast.count = 2;
    typename qn::Station<double>::ServerType slow;
    slow.name = "Slow";
    slow.count = 3;
    net.add_server_type(2, fast);
    net.add_server_type(2, slow);
    net.set_hetero_sched_policy(2, lang::HeteroSchedPolicy::FSF);
    const std::string xml = io::jmt_write_jsim(net.get_struct(), io::JmtWriteOptions());

    CHECK(xml.find("name=\"ServiceStrategy\"") < xml.find("name=\"classParallelism\""));
    CHECK(xml.find("name=\"classParallelism\"") < xml.find("name=\"serverNames\""));
    CHECK(xml.find("name=\"serverNames\"") < xml.find("name=\"serversPerServerType\""));
    CHECK(xml.find("name=\"serversPerServerType\"") < xml.find("name=\"serverCompatibilities\""));
    CHECK(xml.find("name=\"serverCompatibilities\"") < xml.find("name=\"schedulingPolicy\""));
    // The pools size the station, so its server count is their sum
    CHECK(has(xml, "name=\"maxJobs\""));
    CHECK(has(xml, "<value>5</value>"));
    CHECK(has(xml, "FSF (Fastest Servers First)"));
    // Parallelism is not declared here, so every class seizes one server
    CHECK(count(xml, "name=\"serverParallelism\"") == 1);
}

TEST_CASE("jmt writer: parallelism alone still needs the whole pool block") {
    // There is no Server constructor taking classParallelism without the four pool
    // parameters, so a station declaring only parallelism is given one synthetic
    // pool. Its count must be the station's own, since the pools, not maxJobs,
    // size the server pool once any pool is declared.
    qn::Network<double> net = mm1();
    net.set_number_of_servers(2, 4);
    net.set_server_parallelism(2, 1, 2);
    const std::string xml = io::jmt_write_jsim(net.get_struct(), io::JmtWriteOptions());

    CHECK(has(xml, "name=\"classParallelism\""));
    CHECK(has(xml, "<value>Queue1 - Server Type 1</value>"));
    CHECK(has(xml, "name=\"serverTypesNumOfServers\""));
    CHECK(has(xml, "Order (Assign according to order below)"));
    CHECK(count(xml, "name=\"serverTypesNames\"") == 1);
}

TEST_CASE("jmt writer: a station declaring neither writes no pool block") {
    const std::string xml = io::jmt_write_jsim(mm1().get_struct(), io::JmtWriteOptions());
    CHECK(!has(xml, "name=\"classParallelism\""));
    CHECK(!has(xml, "name=\"serverNames\""));
    CHECK(!has(xml, "name=\"schedulingPolicy\""));
}

TEST_CASE("jmt writer: parallelism above the server count is refused") {
    qn::Network<double> net = mm1();
    net.set_number_of_servers(2, 2);
    CHECK_THROWS_AS(net.set_server_parallelism(2, 1, 3), InputError);
}

TEST_CASE("jmva writer: the model is per chain and the algorithm is named") {
    qn::Network<double> net = cqn(3.0);
    const std::string xml = jmva_text(net.get_struct(), "jmva");
    CHECK(has(xml, "xsi:noNamespaceSchemaLocation=\"JMTmodel.xsd\""));
    CHECK(has(xml, "<closedclass population=\"3\" name=\"Chain01\"/>"));
    CHECK(has(xml, "<delaystation name=\"Delay\">"));
    CHECK(has(xml, "<listation "));
    CHECK(has(xml, "name=\"Queue1\""));
    CHECK(has(xml, "<algType name=\"MVA\" tolerance=\"1.0E-7\""));
    CHECK(has(xml, "<compareAlgs value=\"false\"/>"));
    CHECK(has(xml, "<Class name=\"Chain01\" refStation=\"Delay\"/>"));
    // A Source is not a JMVA station; the closed model has none, and the count
    // must equal the number of exported stations.
    CHECK(has(xml, "<stations number=\"2\">"));
}

TEST_CASE("jmva writer: a single-server algorithm refuses a multiserver station") {
    qn::Network<double> net = cqn(3.0);
    net.set_number_of_servers(2, 2);
    CHECK_THROWS_AS(jmva_text(net.get_struct(), "jmva.lin"), UnsupportedError);
    // The exact algorithm has no such restriction.
    CHECK_NOTHROW(jmva_text(net.get_struct(), "jmva"));
}

TEST_CASE("jmt linkAndLog: a Logger on each side of the logged node") {
    // The Network must OUTLIVE the reference: binding `mm1().get_struct()` to a
    // const reference does not extend the temporary's lifetime through the
    // member call, and the struct is destroyed at the end of the statement.
    qn::Network<double> net = mm1();
    const qn::NetworkStruct<double>& sn = net.get_struct();
    std::vector<bool> logged(sn.nodes.size(), false);
    logged[1] = true;  // Queue1
    const qn::NetworkStruct<double> inst = jmt::jmt_link_and_log(sn, logged, "/tmp/logs");

    REQUIRE(inst.nodes.size() == sn.nodes.size() + 2);
    CHECK(inst.nodes[3].name == "Arv_Queue1");
    CHECK(inst.nodes[4].name == "Dep_Queue1");
    CHECK(inst.log_path == "/tmp/logs");
    CHECK(inst.nodes[3].logger.file_name == "Queue1-Arv.csv");
    CHECK(inst.nodes[4].logger.file_name == "Queue1-Dep.csv");
    // The queue keeps its station identity and its service process.
    CHECK(inst.nstations == sn.nstations);
    CHECK(inst.rates(sn.nodes[1].station - 1, 0) == doctest::Approx(1.0));
    // Source -> Arv_Queue1 -> Queue1 -> Dep_Queue1 -> Sink.
    const std::vector<std::vector<bool>> conn = io::jmt_conn_matrix(inst);
    CHECK(conn[0][3]);
    CHECK(conn[3][1]);
    CHECK(conn[1][4]);
    CHECK(conn[4][2]);
    CHECK_FALSE(conn[0][1]);
}

TEST_CASE("jmt end-to-end: an M/M/1 solved by JSIM") {
    if (!jmt::jmt_available()) return;  // no jar or no JVM on this host
    jmt::JmtOptions opt;
    opt.method = "jsim";
    opt.samples = 20000;
    opt.seed = 23000;
    qn::Network<double> net = mm1(0.5, 1.0);
    const qn::NetworkStruct<double>& sn = net.get_struct();
    const jmt::JmtResult<double> res = jmt::solver_jmt_run_analyzer(sn, opt);

    const std::size_t q = sn.nodes[1].station;
    // M/M/1 at rho = 0.5: QLen 1, RespT 2, Util 0.5, Tput 0.5. The tolerance is
    // the simulation's, not the port's: a 20k-sample run of a rho = 0.5 queue
    // has a few percent of standard error on the queue length.
    CHECK(res.avg.QN(q - 1, 0) == doctest::Approx(1.0).epsilon(0.15));
    CHECK(res.avg.RN(q - 1, 0) == doctest::Approx(2.0).epsilon(0.15));
    CHECK(res.avg.UN(q - 1, 0) == doctest::Approx(0.5).epsilon(0.10));
    CHECK(res.avg.TN(q - 1, 0) == doctest::Approx(0.5).epsilon(0.10));
}

TEST_CASE("jmt end-to-end: a closed network solved by JMVA") {
    if (!jmt::jmt_available()) return;
    jmt::JmtOptions opt;
    opt.method = "jmva";
    qn::Network<double> net = cqn(3.0);
    const qn::NetworkStruct<double>& sn = net.get_struct();
    const jmt::JmtResult<double> res = jmt::solver_jmt_run_analyzer(sn, opt);

    // Exact MVA on Delay(1.0) + PS(0.5) with N = 3, by the recursion itself:
    //   n=1  R = 1 + 0.5                 X = 1/1.5     Q = (0.6667, 0.3333)
    //   n=2  R = 1 + 0.5(1 + 0.3333)     X = 2/1.6667  Q = (1.2, 0.8)
    //   n=3  R = 1 + 0.5(1 + 0.8) = 1.9  X = 3/1.9     Q = (1.5789, 1.4211)
    // The delay holds MORE than the queue: its demand is twice as large.
    const std::size_t d = sn.nodes[0].station, q = sn.nodes[1].station;
    CHECK(res.avg.TN(d - 1, 0) == doctest::Approx(3.0 / 1.9).epsilon(1e-4));
    CHECK(res.avg.QN(d - 1, 0) == doctest::Approx(3.0 / 1.9).epsilon(1e-4));
    CHECK(res.avg.QN(q - 1, 0) == doctest::Approx(3.0 - 3.0 / 1.9).epsilon(1e-4));
}

TEST_CASE("jmt end-to-end: the response-time CDF from the JMT logs") {
    if (!jmt::jmt_available()) return;
    jmt::JmtOptions opt;
    opt.method = "jsim";
    opt.samples = 10000;
    qn::Network<double> net = mm1(0.5, 1.0);
    const qn::NetworkStruct<double>& sn = net.get_struct();
    const std::map<std::pair<std::size_t, std::size_t>,
                   std::vector<std::pair<double, double>>>
        rd = jmt::jmt_get_cdf_resp_t(sn, opt);

    const std::pair<std::size_t, std::size_t> key(sn.nodes[1].station, 1);
    REQUIRE(rd.count(key) == 1);
    const std::vector<std::pair<double, double>>& fx = rd.at(key);
    REQUIRE(fx.size() > 10);
    CHECK(fx.front().first == doctest::Approx(0.0));
    CHECK(fx.back().first == doctest::Approx(1.0));
    // The mean of the empirical law must agree with the M/M/1 response time.
    double mean = 0.0;
    for (std::size_t i = 1; i < fx.size(); ++i)
        mean += (fx[i].first - fx[i - 1].first) * fx[i].second;
    CHECK(mean == doctest::Approx(2.0).epsilon(0.20));
}
