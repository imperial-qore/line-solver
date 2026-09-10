/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * The .lqnx reader must refuse a structurally inconsistent document.
 *
 * A defective input used to be accepted in silence and to surface much later,
 * as a NaN metric, an aliased element or a model with no customers at all. The
 * reader now names the defect at its source. The same fourteen documents below
 * are refused with the SAME message by the MATLAB, JAR and Python readers, so
 * this file also pins the wording that keeps the four codebases
 * interchangeable; the twin is python/tests/test_lqn_input_validation.py.
 */

#include <fstream>
#include <string>
#include <vector>

#include "doctest.h"
#include "line/lang/lqn/lqn_reader.h"
#include "line/util/tempdir.h"

using namespace line;

namespace {

const char* HEAD =
    "<?xml version=\"1.0\"?>\n"
    "<lqn-model name=\"t\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\">\n";
const char* TAIL = "</lqn-model>\n";

std::string entry(const std::string& name,
                  const std::vector<std::string>& calls = std::vector<std::string>(),
                  const std::vector<std::pair<std::string, std::string>>& fwd =
                      std::vector<std::pair<std::string, std::string>>(),
                  const std::string& arrival = std::string()) {
    std::string s = "      <entry name=\"" + name + "\" type=\"PH1PH2\"";
    if (!arrival.empty()) s += " open-arrival-rate=\"" + arrival + "\"";
    s += ">\n";
    for (const auto& f : fwd)
        s += "        <forwarding dest=\"" + f.first + "\" prob=\"" + f.second + "\"/>\n";
    s += "        <entry-phase-activities>\n";
    s += "          <activity name=\"" + name + "_ph1\" phase=\"1\" host-demand-mean=\"1\">\n";
    for (const std::string& d : calls)
        s += "            <synch-call dest=\"" + d + "\" calls-mean=\"1\"/>\n";
    s += "          </activity>\n        </entry-phase-activities>\n      </entry>\n";
    return s;
}

std::string task(const std::string& name, const std::string& sched, const std::string& entries,
                 const std::string& extra = std::string()) {
    return "    <task name=\"" + name + "\" scheduling=\"" + sched + "\">\n" + entries + extra +
           "    </task>\n";
}

std::string proc(const std::string& name, const std::string& tasks,
                 const std::string& sched = "fcfs") {
    return "  <processor name=\"" + name + "\" scheduling=\"" + sched + "\">\n" + tasks +
           "  </processor>\n";
}

std::string or_fork(const std::string& p2, const std::string& p3) {
    return std::string(
               "      <task-activities>\n"
               "        <activity name=\"a1\" bound-to-entry=\"e1\" host-demand-mean=\"1\"/>\n"
               "        <activity name=\"a2\" host-demand-mean=\"1\"/>\n"
               "        <activity name=\"a3\" host-demand-mean=\"1\"/>\n"
               "        <precedence>\n"
               "          <pre><activity name=\"a1\"/></pre>\n"
               "          <post-OR><activity name=\"a2\" prob=\"") +
           p2 + "\"/><activity name=\"a3\" prob=\"" + p3 +
           "\"/></post-OR>\n"
           "        </precedence>\n"
           "        <reply-entry name=\"e1\"><reply-activity name=\"a2\"/></reply-entry>\n"
           "      </task-activities>\n";
}

const std::string REF = task("t0", "ref", entry("e0", {"e1"}));
const std::string SRV = task("t1", "fcfs", entry("e1"));

/** Writes the document to the temporary directory and returns its path. */
std::string write_model(const util::TempDir& dir, const std::string& name,
                        const std::string& body) {
    const std::string path = dir.file(name + ".lqnx");
    std::ofstream out(path.c_str());
    out << HEAD << body << TAIL;
    out.close();
    return path;
}

/** The message of the InputError the reader throws, or the empty string. */
std::string refusal(const std::string& path) {
    try {
        lqn::read_lqnx_model<double>(path);
    } catch (const InputError& e) {
        return std::string(e.what());
    }
    return std::string();
}

}  // namespace

TEST_CASE("lqnx reader refuses a structurally inconsistent document") {
    util::TempDir dir("lqn_validation");

    CHECK(refusal(write_model(dir, "dup_processor", proc("p0", REF, "inf") + proc("p0", SRV))) ==
          "Duplicate processor name \"p0\".");

    CHECK(refusal(write_model(dir, "dup_task", proc("p0", REF, "inf") +
                                                   proc("p1", task("t0", "fcfs", entry("e1"))))) ==
          "Duplicate task name \"t0\".");

    CHECK(refusal(write_model(dir, "dup_entry", proc("p0", REF, "inf") +
                                                    proc("p1", task("t1", "fcfs", entry("e0"))))) ==
          "Duplicate entry name \"e0\".");

    CHECK(refusal(write_model(
              dir, "dup_activity",
              proc("p0", REF, "inf") +
                  proc("p1", task("t1", "fcfs", entry("e1"),
                                  "      <task-activities>\n"
                                  "        <activity name=\"e1_ph1\" host-demand-mean=\"1\"/>\n"
                                  "      </task-activities>\n")))) ==
          "Duplicate activity name \"e1_ph1\" in task \"t1\".");

    CHECK(refusal(write_model(dir, "no_entries", proc("p0", REF, "inf") + proc("p1", SRV) +
                                                     proc("p2", task("t2", "fcfs", "")))) ==
          "Task \"t2\" has no entries.");

    CHECK(refusal(write_model(dir, "no_reference_task",
                              proc("p0", task("t0", "fcfs", entry("e0", {"e1"}))) +
                                  proc("p1", SRV))) ==
          "The model has no reference task and no open arrivals.");

    CHECK(refusal(write_model(dir, "ref_receiver",
                              proc("p0", task("t0", "ref", entry("e0", {"e1"})), "inf") +
                                  proc("p1", task("t1", "fcfs", entry("e1", {"e0"}))))) ==
          "Entry \"e0\" belongs to reference task \"t0\" and cannot receive requests.");

    CHECK(refusal(write_model(
              dir, "ref_replies",
              proc("p0",
                   task("t0", "ref", entry("e0", {"e1"}),
                        "      <task-activities>\n"
                        "        <reply-entry name=\"e0\">\n"
                        "          <reply-activity name=\"e0_ph1\"/>\n"
                        "        </reply-entry>\n"
                        "      </task-activities>\n"),
                   "inf") +
                  proc("p1", SRV))) ==
          "Entry \"e0\" belongs to reference task \"t0\" and cannot be replied to.");

    CHECK(refusal(write_model(
              dir, "ref_forwarding",
              proc("p0", task("t0", "ref", entry("e0", {"e1"}, {{"e1", "0.5"}})), "inf") +
                  proc("p1", SRV))) ==
          "Entry \"e0\" belongs to reference task \"t0\" and cannot forward requests.");

    CHECK(refusal(write_model(
              dir, "ref_open_arrivals",
              proc("p0", task("t0", "ref", entry("e0", {"e1"}, {}, "0.5")), "inf") +
                  proc("p1", SRV))) ==
          "Entry \"e0\" belongs to reference task \"t0\" and cannot have open arrivals.");

    CHECK(refusal(write_model(
              dir, "forwarding_probability_negative",
              proc("p0", REF, "inf") +
                  proc("p1", task("t1", "fcfs", entry("e1", {}, {{"e2", "-0.5"}}))) +
                  proc("p2", task("t2", "fcfs", entry("e2"))))) ==
          "Forwarding from entry \"e1\" to entry \"e2\" has an invalid probability of -0.5.");

    CHECK(refusal(write_model(
              dir, "forwarding_probability_total",
              proc("p0", REF, "inf") +
                  proc("p1", task("t1", "fcfs",
                                  entry("e1", {}, {{"e2", "0.7"}, {"e3", "0.7"}}))) +
                  proc("p2", task("t2", "fcfs", entry("e2") + entry("e3"))))) ==
          "Entry \"e1\" has a total forwarding probability of 1.4.");

    CHECK(refusal(write_model(dir, "or_branch_probability_invalid",
                              proc("p0", REF, "inf") +
                                  proc("p1", task("t1", "fcfs",
                                                  "      <entry name=\"e1\" type=\"NONE\"/>\n",
                                                  or_fork("1.4", "-0.4"))))) ==
          "Activity \"a2\" in task \"t1\" has an invalid branch probability of 1.4.");

    CHECK(refusal(write_model(dir, "or_branch_probabilities_sum",
                              proc("p0", REF, "inf") +
                                  proc("p1", task("t1", "fcfs",
                                                  "      <entry name=\"e1\" type=\"NONE\"/>\n",
                                                  or_fork("0.4", "0.4"))))) ==
          "Branch probabilities of an OR-fork in task \"t1\" sum to 0.8 instead of 1.");
}

TEST_CASE("lqnx reader accepts a consistent document") {
    util::TempDir dir("lqn_validation_ok");
    const std::string path = write_model(dir, "good", proc("p0", REF, "inf") + proc("p1", SRV));
    CHECK(refusal(path) == "");
    const lqn::LqnModel<double> m = lqn::read_lqnx_model<double>(path);
    CHECK(m.tasks.size() == 2u);
}
