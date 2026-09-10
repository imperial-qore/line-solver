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

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Two structural refusals that belong to getStruct, not to the parser.
 *
 * <p>Both documents are valid against {@code lqn-core.xsd} and both are accepted by lqns, so
 * nothing in the reader stops them; what they describe is a model LINE's
 * {@code LayeredNetworkStruct} cannot represent, and the refusal is by name.</p>
 *
 * <ol>
 * <li>An entry whose {@code <entry-phase-activities>} is EMPTY reaches no activity, so it has
 * no service and no reply. This codebase and python used to build the struct anyway and report
 * a row of NaN for it; MATLAB and C++ refused.</li>
 * <li>An activity that REPLIES and then continues into a phase-1 successor ends phase 1 of its
 * entry, so the tail is a second phase the struct has no slot for. MATLAB, this codebase and
 * python refused; C++ had no guard.</li>
 * </ol>
 *
 * <p>The wording is {@code getStruct.m}'s, verbatim, in all four codebases. The twins are
 * cpp/tests/test_lqn_struct_guards.cpp, python/tests/test_lqn_struct_guards.py and MATLAB
 * test_lqn_struct_guards.m. {@code user-models/blsr.lqnx} in the LQNS corpus carries BOTH
 * defects, and refuses on the first, as MATLAB does.</p>
 */
public class LqnStructGuardsTest {

    private static final String HEAD =
            "<?xml version=\"1.0\"?>\n"
            + "<lqn-model name=\"g\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\">\n";
    private static final String TAIL = "</lqn-model>\n";

    private static final String SERVER =
            "    <task name=\"TS\" scheduling=\"fcfs\">\n"
            + "      <entry name=\"ES\" type=\"PH1PH2\">\n"
            + "        <entry-phase-activities>\n"
            + "          <activity name=\"ES_ph1\" phase=\"1\" host-demand-mean=\"2\"/>\n"
            + "        </entry-phase-activities>\n"
            + "      </entry>\n"
            + "    </task>\n";

    /** An entry that names no activity at all. */
    private static final String EMPTY_ENTRY =
            "    <task name=\"TN\" scheduling=\"inf\">\n"
            + "      <entry name=\"EN\" type=\"PH1PH2\">\n"
            + "        <entry-phase-activities>\n"
            + "        </entry-phase-activities>\n"
            + "      </entry>\n"
            + "    </task>\n";

    /** a1 replies, a2 follows it and carries no phase attribute, so it is phase 1. */
    private static final String REPLY_THEN_CONTINUE =
            "    <task name=\"TR\" scheduling=\"fcfs\">\n"
            + "      <entry name=\"ER\" type=\"NONE\"/>\n"
            + "      <task-activities>\n"
            + "        <activity name=\"a1\" bound-to-entry=\"ER\" host-demand-mean=\"1\"/>\n"
            + "        <activity name=\"a2\" host-demand-mean=\"1\"/>\n"
            + "        <precedence>\n"
            + "          <pre><activity name=\"a1\"/></pre>\n"
            + "          <post><activity name=\"a2\"/></post>\n"
            + "        </precedence>\n"
            + "        <reply-entry name=\"ER\"><reply-activity name=\"a1\"/></reply-entry>\n"
            + "      </task-activities>\n"
            + "    </task>\n";

    private static String client(String dest) {
        return "  <processor name=\"PC\" scheduling=\"inf\">\n"
                + "    <task name=\"TC\" scheduling=\"ref\" multiplicity=\"2\">\n"
                + "      <entry name=\"EC\" type=\"PH1PH2\">\n"
                + "        <entry-phase-activities>\n"
                + "          <activity name=\"EC_ph1\" phase=\"1\" host-demand-mean=\"1\">\n"
                + "            <synch-call dest=\"" + dest + "\" calls-mean=\"1\"/>\n"
                + "          </activity>\n"
                + "        </entry-phase-activities>\n"
                + "      </entry>\n"
                + "    </task>\n"
                + "  </processor>\n";
    }

    private static String build(Path dir, String name, String extraTask, String dest)
            throws IOException {
        String body = HEAD + client(dest) + "  <processor name=\"PS\" scheduling=\"fcfs\">\n"
                + SERVER + extraTask + "  </processor>\n" + TAIL;
        File f = new File(dir.toFile(), name + ".lqnx");
        FileWriter w = new FileWriter(f);
        w.write(body);
        w.close();
        return f.getAbsolutePath();
    }

    /** The message getStruct throws, or the empty string. */
    private static String refusal(String path) {
        try {
            LayeredNetwork.parseXML(path).getStruct();
        } catch (Exception e) {
            return e.getMessage() == null ? "" : e.getMessage();
        }
        return "";
    }

    @Test
    public void entryWithNoBoundActivityIsRefused(@TempDir Path dir) throws IOException {
        String msg = refusal(build(dir, "empty_entry", EMPTY_ENTRY, "ES"));
        assertTrue(msg.contains("An entry does not have any boundTo activity."), msg);
    }

    @Test
    public void replyThenContinueIsRefused(@TempDir Path dir) throws IOException {
        String msg = refusal(build(dir, "reply_tail", REPLY_THEN_CONTINUE, "ER"));
        assertTrue(msg.contains("Unsupported replyTo in non-terminal activity."), msg);
    }

    @Test
    public void aSoundModelStillBuilds(@TempDir Path dir) throws IOException {
        LayeredNetworkStruct lsn = LayeredNetwork.parseXML(build(dir, "ok", "", "ES")).getStruct();
        assertEquals(2, lsn.nentries);
        assertEquals(2, lsn.ntasks);
    }
}
