/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.layered;

import jline.lang.constant.SchedStrategy;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * The .lqnx "pri" discipline is preemptive priority resume, not head-of-line.
 *
 * <p>lqns writes SCHEDULE_PPR as "pri" (LQIO::SCHEDULE::PPR in lqiolib labels.cpp), and PPR
 * preempts: a more urgent arrival takes the server and the displaced job resumes its residual
 * work (lqns/processor.cc builds a PR_FCFS_Server). LINE holds that as FCFSPRPRIO. The readers
 * used to disagree on the token -- MATLAB and this one refused it outright, Python read it as
 * FCFS and C++ as the NON-preemptive HOL -- so the same file was three different models.</p>
 *
 * <p>"pp" is the spelling of lqn-core.xsd, which is stale: it appears nowhere in the lqns
 * 6.2.31 sources, so no file lqns produces carries it. It is accepted on the read side and never
 * written. The twins are cpp/tests/test_lqnx_pri_sched.cpp, python/tests/test_lqnx_pri_sched.py
 * and MATLAB test_lqnx_pri_sched.m.</p>
 */
public class LqnxPriSchedTest {

    /** A reference task on an inf host calling one server task, both spelled SCHED. */
    private static String document(String procSched, String taskSched) {
        return "<?xml version=\"1.0\"?>\n"
                + "<lqn-model name=\"pri\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\">\n"
                + "  <processor name=\"PC\" scheduling=\"inf\">\n"
                + "    <task name=\"TC\" scheduling=\"ref\" multiplicity=\"2\">\n"
                + "      <entry name=\"EC\" type=\"PH1PH2\">\n"
                + "        <entry-phase-activities>\n"
                + "          <activity name=\"EC_ph1\" phase=\"1\" host-demand-mean=\"1\">\n"
                + "            <synch-call dest=\"ES\" calls-mean=\"1\"/>\n"
                + "          </activity>\n"
                + "        </entry-phase-activities>\n"
                + "      </entry>\n"
                + "    </task>\n"
                + "  </processor>\n"
                + "  <processor name=\"PS\" scheduling=\"" + procSched + "\">\n"
                + "    <task name=\"TS\" scheduling=\"" + taskSched + "\" priority=\"3\">\n"
                + "      <entry name=\"ES\" type=\"PH1PH2\">\n"
                + "        <entry-phase-activities>\n"
                + "          <activity name=\"ES_ph1\" phase=\"1\" host-demand-mean=\"2\"/>\n"
                + "        </entry-phase-activities>\n"
                + "      </entry>\n"
                + "    </task>\n"
                + "  </processor>\n"
                + "</lqn-model>\n";
    }

    private static String write(Path dir, String name, String body) throws IOException {
        File f = new File(dir.toFile(), name + ".lqnx");
        FileWriter w = new FileWriter(f);
        w.write(body);
        w.close();
        return f.getAbsolutePath();
    }

    /** The scheduling strategy of the host named NAME. */
    private static SchedStrategy hostSched(LayeredNetwork model, String name) {
        for (Host h : model.getHosts().values()) {
            if (h.getName().equals(name)) {
                return h.getScheduling();
            }
        }
        throw new IllegalArgumentException("no host " + name);
    }

    /** The scheduling strategy of the task named NAME. */
    private static SchedStrategy taskSched(LayeredNetwork model, String name) {
        for (Task t : model.getTasks().values()) {
            if (t.getName().equals(name)) {
                return t.getScheduling();
            }
        }
        throw new IllegalArgumentException("no task " + name);
    }

    @Test
    public void priIsPreemptivePriorityResume() {
        assertEquals(SchedStrategy.FCFSPRPRIO, SchedStrategy.fromText("pri"));
        assertEquals(SchedStrategy.FCFSPRPRIO, SchedStrategy.fromText("pp"));
    }

    @Test
    public void readerGivesPriToAProcessorAndATask(@TempDir Path dir) throws IOException {
        LayeredNetwork model = LayeredNetwork.parseXML(write(dir, "pri", document("pri", "pri")));
        assertEquals(SchedStrategy.FCFSPRPRIO, hostSched(model, "PS"));
        assertEquals(SchedStrategy.FCFSPRPRIO, taskSched(model, "TS"));

        LayeredNetwork stale = LayeredNetwork.parseXML(write(dir, "pp", document("pp", "pp")));
        assertEquals(SchedStrategy.FCFSPRPRIO, hostSched(stale, "PS"));
        assertEquals(SchedStrategy.FCFSPRPRIO, taskSched(stale, "TS"));
    }

    @Test
    public void writerEmitsPriTheSpellingLqnsReadsBack(@TempDir Path dir) throws IOException {
        LayeredNetwork model = LayeredNetwork.parseXML(write(dir, "in", document("pri", "pri")));
        File out = new File(dir.toFile(), "out.lqnx");
        model.writeXML(out.getAbsolutePath());

        String txt = new String(Files.readAllBytes(out.toPath()));
        assertTrue(txt.contains("scheduling=\"pri\""), "writer must emit the lqns spelling");
        // fcfsprprio is LINE's own name for the discipline and is not valid LQN
        assertFalse(txt.contains("fcfsprprio"), "LINE's own name is not a valid LQN token");

        LayeredNetwork back = LayeredNetwork.parseXML(out.getAbsolutePath());
        assertEquals(SchedStrategy.FCFSPRPRIO, hostSched(back, "PS"));
        assertEquals(SchedStrategy.FCFSPRPRIO, taskSched(back, "TS"));
    }
}
