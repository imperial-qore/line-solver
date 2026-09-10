/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.ag;

import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.processes.Exp;
import jline.solvers.NetworkAvgTable;
import jline.solvers.mva.SolverMVA;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * The RCAT method arm -- 'inap', 'inapplus' and the vestigial 'exact' -- on the
 * two shapes the reference exercises: a closed Delay+Queue and an open M/M/1.
 *
 * These cases used to live in SolverMAMCoverageTest and address SolverMAM,
 * because RCAT was a MAM method arm. It is SolverAG's since the move, and
 * SolverMAM refuses the names outright, so the cases address SolverAG here.
 * The models, the cross-validation against MVA and the tolerances are the ones
 * they carried: what changed is which solver answers, not what is expected of
 * the answer.
 */
public class SolverAGMethodsTest {

    /** A closed Delay+Queue cycle, one class, exponential everywhere. */
    private Network buildClosedDelayQueue(int N, double thinkRate, double serviceRate) {
        Network model = new Network("ClosedDQ");
        Delay delay = new Delay(model, "Delay");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        ClosedClass cc = new ClosedClass(model, "Class1", N, delay);
        delay.setService(cc, new Exp(thinkRate));
        queue.setService(cc, new Exp(serviceRate));
        model.link(model.serialRouting(delay, queue));
        return model;
    }

    /** An open M/M/1. */
    private Network buildOpenMM1(double lambda, double mu) {
        Network model = new Network("OpenMM1");
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");
        OpenClass oc = new OpenClass(model, "Class1", 0);
        source.setArrival(oc, new Exp(lambda));
        queue.setService(oc, new Exp(mu));
        model.link(model.serialRouting(source, queue, sink));
        return model;
    }

    private double getMetric(NetworkAvgTable table, String stationName, String className, String metric) {
        for (int i = 0; i < table.getStationNames().size(); i++) {
            if (table.getStationNames().get(i).equals(stationName)
                    && table.getClassNames().get(i).equals(className)) {
                switch (metric) {
                    case "QLen": return table.getQLen().get(i);
                    case "Util": return table.getUtil().get(i);
                    case "RespT": return table.getRespT().get(i);
                    case "Tput": return table.getTput().get(i);
                    case "ArvR": return table.getArvR().get(i);
                    default: return Double.NaN;
                }
            }
        }
        return Double.NaN;
    }

    @Test
    public void testINAP_closedNetwork() {
        Network model = buildClosedDelayQueue(5, 1.0, 2.0);

        SolverAG solver = new SolverAG(model, "inap");
        NetworkAvgTable table = solver.getAvgTable();
        assertNotNull(table, "INAP should produce results");

        // Cross-validate against MVA
        Network model2 = buildClosedDelayQueue(5, 1.0, 2.0);
        SolverMVA mvaSolver = new SolverMVA(model2);
        NetworkAvgTable mvaTable = mvaSolver.getAvgTable();

        double agTput = getMetric(table, "Queue", "Class1", "Tput");
        double mvaTput = getMetric(mvaTable, "Queue", "Class1", "Tput");

        assertTrue(agTput > 0, "INAP throughput should be positive");
        // INAP is an iterative approximation; wider tolerance needed
        assertEquals(mvaTput, agTput, 0.35 * mvaTput,
            "INAP should approximate MVA within 35% for Exp service");
    }

    @Test
    public void testINAPPlus_closedNetwork() {
        Network model = buildClosedDelayQueue(5, 1.0, 2.0);

        SolverAG solver = new SolverAG(model, "inapplus");
        NetworkAvgTable table = solver.getAvgTable();
        assertNotNull(table, "INAPPLUS should produce results");

        double tput = getMetric(table, "Queue", "Class1", "Tput");
        assertTrue(tput > 0, "INAPPLUS throughput should be positive");
    }

    @Test
    public void testExact_closedNetwork() {
        // "exact" falls back to INAP in JAR
        Network model = buildClosedDelayQueue(5, 1.0, 2.0);

        SolverAG solver = new SolverAG(model, "exact");
        NetworkAvgTable table = solver.getAvgTable();
        assertNotNull(table, "exact (fallback to INAP) should produce results");
    }

    @Test
    public void testINAP_openNetwork() {
        Network model = buildOpenMM1(0.5, 1.0);

        SolverAG solver = new SolverAG(model, "inap");
        NetworkAvgTable table = solver.getAvgTable();
        assertNotNull(table, "INAP should produce results for open network");
    }

    /**
     * The four RCAT names SolverMAM handed over are the ones SolverAG lists, so
     * a caller that followed the redirect finds them where it was sent.
     */
    @Test
    public void testListValidMethods() {
        SolverAG solver = new SolverAG(buildOpenMM1(0.5, 1.0));
        java.util.List<String> methods = solver.listValidMethods();

        assertTrue(methods.contains("default"));
        assertTrue(methods.contains("inap"));
        assertTrue(methods.contains("inapplus"));
        assertTrue(methods.contains("inapinf"));
        assertTrue(methods.contains("exact"));
    }
}
