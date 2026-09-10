/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * The .lqnx `pri` discipline is preemptive priority resume, not head-of-line.
 *
 * lqns writes SCHEDULE_PPR as "pri" (LQIO::SCHEDULE::PPR in lqiolib
 * labels.cpp), and PPR preempts: a more urgent arrival takes the server and the
 * displaced job resumes its residual work (lqns/processor.cc builds a
 * PR_FCFS_Server). This reader used to fold "pri" onto HOL, which is
 * NON-preemptive, so a model wrote out one discipline and was answered as
 * another with nothing in the result saying so. It is FCFSPRPRIO.
 *
 * "pp" is the spelling of lqn-core.xsd, which is stale: it appears nowhere in
 * the lqns 6.2.31 sources, so no file lqns produces carries it. It is accepted
 * on the read side and never written.
 *
 * The twins are jar LqnxPriSchedTest.java, python test_lqnx_pri_sched.py and
 * MATLAB test_lqnx_pri_sched.m.
 */

#include <fstream>
#include <string>

#include "doctest.h"
#include "line/lang/lqn/lqn_reader.h"
#include "line/lang/lqn/lqn_writer.h"
#include "line/util/tempdir.h"

using namespace line;
using line::lang::SchedStrategy;

namespace {

/** A reference task on an inf host calling one server task on host PS. */
std::string document(const std::string& proc_sched, const std::string& task_sched) {
    return std::string(
               "<?xml version=\"1.0\"?>\n"
               "<lqn-model name=\"pri\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\">\n"
               "  <processor name=\"PC\" scheduling=\"inf\">\n"
               "    <task name=\"TC\" scheduling=\"ref\" multiplicity=\"2\">\n"
               "      <entry name=\"EC\" type=\"PH1PH2\">\n"
               "        <entry-phase-activities>\n"
               "          <activity name=\"EC_ph1\" phase=\"1\" host-demand-mean=\"1\">\n"
               "            <synch-call dest=\"ES\" calls-mean=\"1\"/>\n"
               "          </activity>\n"
               "        </entry-phase-activities>\n"
               "      </entry>\n"
               "    </task>\n"
               "  </processor>\n"
               "  <processor name=\"PS\" scheduling=\"") +
           proc_sched +
           "\">\n"
           "    <task name=\"TS\" scheduling=\"" +
           task_sched +
           "\" priority=\"3\">\n"
           "      <entry name=\"ES\" type=\"PH1PH2\">\n"
           "        <entry-phase-activities>\n"
           "          <activity name=\"ES_ph1\" phase=\"1\" host-demand-mean=\"2\"/>\n"
           "        </entry-phase-activities>\n"
           "      </entry>\n"
           "    </task>\n"
           "  </processor>\n"
           "</lqn-model>\n";
}

std::string write_model(const util::TempDir& dir, const std::string& name,
                        const std::string& body) {
    const std::string path = dir.file(name + ".lqnx");
    std::ofstream out(path.c_str());
    out << body;
    out.close();
    return path;
}

std::string slurp(const std::string& path) {
    std::ifstream in(path.c_str());
    return std::string((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
}

}  // namespace

TEST_CASE("lqnx pri is preemptive priority resume on both a processor and a task") {
    util::TempDir dir("lqnx_pri");

    const lqn::LqnModel<double> m =
        lqn::read_lqnx_model<double>(write_model(dir, "pri", document("pri", "pri")));
    REQUIRE(m.procs.size() == 2u);
    REQUIRE(m.tasks.size() == 2u);
    CHECK(m.procs[1].sched == SchedStrategy::FCFSPRPRIO);
    CHECK(m.tasks[1].sched == SchedStrategy::FCFSPRPRIO);
    CHECK(m.procs[1].sched != SchedStrategy::HOL);

    // the stale schema spelling reads the same way
    const lqn::LqnModel<double> mpp =
        lqn::read_lqnx_model<double>(write_model(dir, "pp", document("pp", "pp")));
    CHECK(mpp.procs[1].sched == SchedStrategy::FCFSPRPRIO);
    CHECK(mpp.tasks[1].sched == SchedStrategy::FCFSPRPRIO);
}

TEST_CASE("lqnx writer emits pri, the spelling lqns reads back") {
    util::TempDir dir("lqnx_pri_rt");

    const lqn::LqnModel<double> m =
        lqn::read_lqnx_model<double>(write_model(dir, "in", document("pri", "pri")));
    const std::string out = dir.file("out.lqnx");
    lqn::write_lqnx(m, out);

    const std::string txt = slurp(out);
    CHECK(txt.find("scheduling=\"pri\"") != std::string::npos);
    // fcfsprprio is LINE's own name for the discipline and is not valid LQN
    CHECK(txt.find("fcfsprprio") == std::string::npos);

    const lqn::LqnModel<double> back = lqn::read_lqnx_model<double>(out);
    CHECK(back.procs[1].sched == SchedStrategy::FCFSPRPRIO);
    CHECK(back.tasks[1].sched == SchedStrategy::FCFSPRPRIO);
}
