/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.layered;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assertions.fail;

/**
 * The .lqnx reader must refuse a structurally inconsistent document.
 *
 * <p>A defective input used to be accepted in silence and to surface much later, as a NaN
 * metric, an aliased element or a model with no customers at all. The reader now names the
 * defect at its source. The same fourteen documents below are refused with the SAME message by
 * the MATLAB, Python and C++ readers, so this file also pins the wording that keeps the four
 * codebases interchangeable; the twin is python/tests/test_lqn_input_validation.py.</p>
 */
public class LqnInputValidationTest {

    private static final String HEAD =
            "<?xml version=\"1.0\"?>\n"
            + "<lqn-model name=\"t\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\">\n";
    private static final String TAIL = "</lqn-model>\n";

    private static String entry(String name, List<String> calls, String[][] fwd, String arrival) {
        StringBuilder s = new StringBuilder("      <entry name=\"" + name + "\" type=\"PH1PH2\"");
        if (arrival != null) {
            s.append(" open-arrival-rate=\"").append(arrival).append("\"");
        }
        s.append(">\n");
        for (String[] f : fwd) {
            s.append("        <forwarding dest=\"").append(f[0]).append("\" prob=\"").append(f[1]).append("\"/>\n");
        }
        s.append("        <entry-phase-activities>\n");
        s.append("          <activity name=\"").append(name).append("_ph1\" phase=\"1\" host-demand-mean=\"1\">\n");
        for (String d : calls) {
            s.append("            <synch-call dest=\"").append(d).append("\" calls-mean=\"1\"/>\n");
        }
        s.append("          </activity>\n        </entry-phase-activities>\n      </entry>\n");
        return s.toString();
    }

    private static String entry(String name) {
        return entry(name, Arrays.<String>asList(), new String[0][], null);
    }

    private static String entry(String name, String call) {
        return entry(name, Arrays.asList(call), new String[0][], null);
    }

    private static String task(String name, String sched, String entries, String extra) {
        return "    <task name=\"" + name + "\" scheduling=\"" + sched + "\">\n" + entries + extra + "    </task>\n";
    }

    private static String task(String name, String sched, String entries) {
        return task(name, sched, entries, "");
    }

    private static String proc(String name, String tasks, String sched) {
        return "  <processor name=\"" + name + "\" scheduling=\"" + sched + "\">\n" + tasks + "  </processor>\n";
    }

    private static String proc(String name, String tasks) {
        return proc(name, tasks, "fcfs");
    }

    private static String orFork(String p2, String p3) {
        return "      <task-activities>\n"
                + "        <activity name=\"a1\" bound-to-entry=\"e1\" host-demand-mean=\"1\"/>\n"
                + "        <activity name=\"a2\" host-demand-mean=\"1\"/>\n"
                + "        <activity name=\"a3\" host-demand-mean=\"1\"/>\n"
                + "        <precedence>\n"
                + "          <pre><activity name=\"a1\"/></pre>\n"
                + "          <post-OR><activity name=\"a2\" prob=\"" + p2 + "\"/>"
                + "<activity name=\"a3\" prob=\"" + p3 + "\"/></post-OR>\n"
                + "        </precedence>\n"
                + "        <reply-entry name=\"e1\"><reply-activity name=\"a2\"/></reply-entry>\n"
                + "      </task-activities>\n";
    }

    private static final String REF = task("t0", "ref", entry("e0", "e1"));
    private static final String SRV = task("t1", "fcfs", entry("e1"));

    @TempDir
    Path tmp;

    private String write(String name, String body) throws IOException {
        File f = new File(tmp.toFile(), name + ".lqnx");
        FileWriter w = new FileWriter(f);
        w.write(HEAD + body + TAIL);
        w.close();
        return f.getAbsolutePath();
    }

    /** The message of the error the reader raises, or the empty string when it accepts. */
    private static String refusal(String path) {
        try {
            LayeredNetwork.parseXML(path, false);
        } catch (Throwable e) {
            String msg = e.getMessage() == null ? "" : e.getMessage();
            int at = msg.indexOf("] ");
            return at >= 0 ? msg.substring(at + 2) : msg;
        }
        return "";
    }

    @Test
    public void duplicateNamesAreRefused() throws IOException {
        assertEquals("Duplicate processor name \"p0\".",
                refusal(write("dup_processor", proc("p0", REF, "inf") + proc("p0", SRV))));
        assertEquals("Duplicate task name \"t0\".",
                refusal(write("dup_task", proc("p0", REF, "inf") + proc("p1", task("t0", "fcfs", entry("e1"))))));
        assertEquals("Duplicate entry name \"e0\".",
                refusal(write("dup_entry", proc("p0", REF, "inf") + proc("p1", task("t1", "fcfs", entry("e0"))))));
        assertEquals("Duplicate activity name \"e1_ph1\" in task \"t1\".",
                refusal(write("dup_activity", proc("p0", REF, "inf")
                        + proc("p1", task("t1", "fcfs", entry("e1"),
                        "      <task-activities>\n"
                                + "        <activity name=\"e1_ph1\" host-demand-mean=\"1\"/>\n"
                                + "      </task-activities>\n")))));
    }

    @Test
    public void taskWithoutEntriesIsRefused() throws IOException {
        assertEquals("Task \"t2\" has no entries.",
                refusal(write("no_entries", proc("p0", REF, "inf") + proc("p1", SRV)
                        + proc("p2", task("t2", "fcfs", "")))));
    }

    @Test
    public void modelWithoutCustomersIsRefused() throws IOException {
        assertEquals("The model has no reference task and no open arrivals.",
                refusal(write("no_reference_task", proc("p0", task("t0", "fcfs", entry("e0", "e1"))) + proc("p1", SRV))));
    }

    @Test
    public void referenceTaskEntryIsNeitherServerNorForwarder() throws IOException {
        assertEquals("Entry \"e0\" belongs to reference task \"t0\" and cannot receive requests.",
                refusal(write("ref_receiver", proc("p0", task("t0", "ref", entry("e0", "e1")), "inf")
                        + proc("p1", task("t1", "fcfs", entry("e1", "e0"))))));
        assertEquals("Entry \"e0\" belongs to reference task \"t0\" and cannot be replied to.",
                refusal(write("ref_replies", proc("p0", task("t0", "ref", entry("e0", "e1"),
                        "      <task-activities>\n"
                                + "        <reply-entry name=\"e0\">\n"
                                + "          <reply-activity name=\"e0_ph1\"/>\n"
                                + "        </reply-entry>\n"
                                + "      </task-activities>\n"), "inf") + proc("p1", SRV))));
        assertEquals("Entry \"e0\" belongs to reference task \"t0\" and cannot forward requests.",
                refusal(write("ref_forwarding", proc("p0", task("t0", "ref",
                        entry("e0", Arrays.asList("e1"), new String[][]{{"e1", "0.5"}}, null)), "inf")
                        + proc("p1", SRV))));
        assertEquals("Entry \"e0\" belongs to reference task \"t0\" and cannot have open arrivals.",
                refusal(write("ref_open_arrivals", proc("p0", task("t0", "ref",
                        entry("e0", Arrays.asList("e1"), new String[0][], "0.5")), "inf")
                        + proc("p1", SRV))));
    }

    @Test
    public void invalidProbabilitiesAreRefused() throws IOException {
        assertEquals("Forwarding from entry \"e1\" to entry \"e2\" has an invalid probability of -0.5.",
                refusal(write("forwarding_probability_negative", proc("p0", REF, "inf")
                        + proc("p1", task("t1", "fcfs", entry("e1", Arrays.<String>asList(), new String[][]{{"e2", "-0.5"}}, null)))
                        + proc("p2", task("t2", "fcfs", entry("e2"))))));
        assertEquals("Entry \"e1\" has a total forwarding probability of 1.4.",
                refusal(write("forwarding_probability_total", proc("p0", REF, "inf")
                        + proc("p1", task("t1", "fcfs", entry("e1", Arrays.<String>asList(),
                        new String[][]{{"e2", "0.7"}, {"e3", "0.7"}}, null)))
                        + proc("p2", task("t2", "fcfs", entry("e2") + entry("e3"))))));
        assertEquals("Activity \"a2\" in task \"t1\" has an invalid branch probability of 1.4.",
                refusal(write("or_branch_probability_invalid", proc("p0", REF, "inf")
                        + proc("p1", task("t1", "fcfs", "      <entry name=\"e1\" type=\"NONE\"/>\n",
                        orFork("1.4", "-0.4"))))));
        assertEquals("Branch probabilities of an OR-fork in task \"t1\" sum to 0.8 instead of 1.",
                refusal(write("or_branch_probabilities_sum", proc("p0", REF, "inf")
                        + proc("p1", task("t1", "fcfs", "      <entry name=\"e1\" type=\"NONE\"/>\n",
                        orFork("0.4", "0.4"))))));
    }

    @Test
    public void consistentDocumentIsAccepted() throws IOException {
        String path = write("good", proc("p0", REF, "inf") + proc("p1", SRV));
        assertEquals("", refusal(path));
        LayeredNetwork model = LayeredNetwork.parseXML(path, false);
        if (model == null) {
            fail("the reader returned no model");
        }
        assertTrue(model.getTasks().size() == 2);
    }
}
