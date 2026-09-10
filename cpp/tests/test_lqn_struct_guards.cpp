/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * Two structural refusals that belong to lqn_finalize (getStruct), not to the
 * parser.
 *
 * Both documents are valid against lqn-core.xsd and both are accepted by lqns,
 * so nothing in the reader stops them; what they describe is a model LINE's
 * LqnStruct cannot represent, and the refusal is by name.
 *
 *   1. An entry whose <entry-phase-activities> is EMPTY reaches no activity, so
 *      it has no service and no reply. The JAR and python used to build the
 *      struct anyway and report a row of NaN for it; MATLAB and this port
 *      refused, and this port's message is now getStruct.m's own.
 *   2. An activity that REPLIES and then continues into a phase-1 successor
 *      ends phase 1 of its entry, so the tail is a second phase the struct has
 *      no slot for. MATLAB, the JAR and python refused; this port had NO guard,
 *      which is the silent case: it would have served the tail as though the
 *      reply had not happened.
 *
 * The wording is getStruct.m's, verbatim, in all four codebases. The twins are
 * jar LqnStructGuardsTest.java, python/tests/test_lqn_struct_guards.py and
 * MATLAB test_lqn_struct_guards.m. `user-models/blsr.lqnx` in the LQNS corpus
 * carries BOTH defects, and refuses on the first, as MATLAB does.
 */

#include <fstream>
#include <string>

#include "doctest.h"
#include "line/lang/lqn/lqn_reader.h"
#include "line/util/tempdir.h"

using namespace line;

namespace {

const char* HEAD =
    "<?xml version=\"1.0\"?>\n"
    "<lqn-model name=\"g\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\">\n";
const char* TAIL = "</lqn-model>\n";

const char* SERVER =
    "    <task name=\"TS\" scheduling=\"fcfs\">\n"
    "      <entry name=\"ES\" type=\"PH1PH2\">\n"
    "        <entry-phase-activities>\n"
    "          <activity name=\"ES_ph1\" phase=\"1\" host-demand-mean=\"2\"/>\n"
    "        </entry-phase-activities>\n"
    "      </entry>\n"
    "    </task>\n";

/** An entry that names no activity at all. */
const char* EMPTY_ENTRY =
    "    <task name=\"TN\" scheduling=\"inf\">\n"
    "      <entry name=\"EN\" type=\"PH1PH2\">\n"
    "        <entry-phase-activities>\n"
    "        </entry-phase-activities>\n"
    "      </entry>\n"
    "    </task>\n";

/** a1 replies, a2 follows it and carries no phase attribute, so it is phase 1. */
const char* REPLY_THEN_CONTINUE =
    "    <task name=\"TR\" scheduling=\"fcfs\">\n"
    "      <entry name=\"ER\" type=\"NONE\"/>\n"
    "      <task-activities>\n"
    "        <activity name=\"a1\" bound-to-entry=\"ER\" host-demand-mean=\"1\"/>\n"
    "        <activity name=\"a2\" host-demand-mean=\"1\"/>\n"
    "        <precedence>\n"
    "          <pre><activity name=\"a1\"/></pre>\n"
    "          <post><activity name=\"a2\"/></post>\n"
    "        </precedence>\n"
    "        <reply-entry name=\"ER\"><reply-activity name=\"a1\"/></reply-entry>\n"
    "      </task-activities>\n"
    "    </task>\n";

std::string client(const std::string& dest) {
    return std::string(
               "  <processor name=\"PC\" scheduling=\"inf\">\n"
               "    <task name=\"TC\" scheduling=\"ref\" multiplicity=\"2\">\n"
               "      <entry name=\"EC\" type=\"PH1PH2\">\n"
               "        <entry-phase-activities>\n"
               "          <activity name=\"EC_ph1\" phase=\"1\" host-demand-mean=\"1\">\n"
               "            <synch-call dest=\"") +
           dest +
           "\" calls-mean=\"1\"/>\n"
           "          </activity>\n"
           "        </entry-phase-activities>\n"
           "      </entry>\n"
           "    </task>\n"
           "  </processor>\n";
}

std::string build(const util::TempDir& dir, const std::string& name, const std::string& extra_task,
                  const std::string& dest) {
    const std::string path = dir.file(name + ".lqnx");
    std::ofstream out(path.c_str());
    out << HEAD << client(dest) << "  <processor name=\"PS\" scheduling=\"fcfs\">\n"
        << SERVER << extra_task << "  </processor>\n"
        << TAIL;
    out.close();
    return path;
}

/** The message lqn_finalize throws, or the empty string. */
std::string refusal(const std::string& path) {
    try {
        lqn::read_lqnx<double>(path);
    } catch (const InputError& e) {
        return std::string(e.what());
    }
    return std::string();
}

}  // namespace

TEST_CASE("lqn_finalize refuses an entry with no bound activity") {
    util::TempDir dir("lqn_guards_entry");
    CHECK(refusal(build(dir, "empty_entry", EMPTY_ENTRY, "ES")) ==
          "An entry does not have any boundTo activity.");
}

TEST_CASE("lqn_finalize refuses a reply that continues into phase 1") {
    util::TempDir dir("lqn_guards_reply");
    CHECK(refusal(build(dir, "reply_tail", REPLY_THEN_CONTINUE, "ER")) ==
          "Unsupported replyTo in non-terminal activity.");
}

TEST_CASE("lqn_finalize still builds a sound model") {
    util::TempDir dir("lqn_guards_ok");
    const lqn::LqnStruct<double> l = lqn::read_lqnx<double>(build(dir, "ok", "", "ES"));
    CHECK(l.nentries == 2u);
    CHECK(l.ntasks == 2u);
}
