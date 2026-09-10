/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * `sn_print`, the full NetworkStruct debug dump.
 *
 * THE ORACLE IS THE MODEL ITSELF: every asserted line restates a struct value
 * the builder tests already pin against MATLAB (rates, njobs, chains, visits),
 * so the check here is that the dump carries those values verbatim in the
 * reference's `field: value` notation, that integer-valued entries print as
 * integers, and that every advertised field name appears exactly once.
 */

#include <string>
#include <vector>

#include "doctest.h"
#include "line/api/sn/sn_print.h"
#include "line/lang/qn/network_builder.h"

using namespace line;
using lang::Distrib;
using lang::SchedStrategy;
using D = Distrib<double>;

namespace {
std::size_t count_of(const std::string& hay, const std::string& needle) {
    std::size_t n = 0;
    for (std::size_t p = hay.find(needle); p != std::string::npos;
         p = hay.find(needle, p + needle.size()))
        ++n;
    return n;
}
}  // namespace

TEST_CASE("sn_print dumps the closed two-station struct in the reference notation") {
    qn::Network<double> m("cqn");
    const std::size_t d = m.add_delay("Delay");
    const std::size_t q = m.add_queue("Queue", SchedStrategy::FCFS);
    const std::size_t c = m.add_closed_class("C1", 3.0, d);
    m.set_service(d, c, D::exp_rate(1.0));
    m.set_service(q, c, D::exp_rate(2.0));
    qn::RoutingMatrix<double> P;
    P.set(d, q, 1.0);
    P.set(q, d, 1.0);
    m.link(P);

    const std::string s = api::sn_print(m.get_struct());

    CHECK(s.find("nstations: 2\n") != std::string::npos);
    CHECK(s.find("nnodes: 2\n") != std::string::npos);
    CHECK(s.find("nclasses: 1\n") != std::string::npos);
    CHECK(s.find("nclosedjobs: 3\n") != std::string::npos);
    CHECK(s.find("nchains: 1\n") != std::string::npos);
    CHECK(s.find("njobs: [3]") != std::string::npos);     // integer-rendered
    CHECK(s.find("rates: [1; 2]") != std::string::npos);  // column per class
    CHECK(s.find("chains: [1]") != std::string::npos);
    CHECK(s.find("classnames: [\"C1\"]") != std::string::npos);
    CHECK(s.find("nodenames: [\"Delay\", \"Queue\"]") != std::string::npos);
    CHECK(s.find("nodetype: [Delay, Queue]") != std::string::npos);
    CHECK(s.find("\"Queue\": \"fcfs\"") != std::string::npos);
    CHECK(s.find("\"type\": \"Exp\"") != std::string::npos);
    // every advertised field appears exactly once
    const char* fields[] = {"nstations: ", "refstat: ",  "nservers: ",  "scv: ",
                            "phases: ",    "rt: ",       "rtnodes: ",   "cap: ",
                            "classcap: ",  "refclass: ", "sched: ",     "proc: ",
                            "inchain: ",   "visits: ",   "droprule: ",  "csmatrix: "};
    for (const char* f : fields) CHECK_MESSAGE(count_of(s, std::string("\n") + f) +
                                                       (s.rfind(f, 0) == 0 ? 1 : 0) ==
                                                   1,
                                               f);
    CHECK(s.back() == '\n');
}
