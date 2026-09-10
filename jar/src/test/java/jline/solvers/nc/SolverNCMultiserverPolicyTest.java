/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.nc;

import jline.VerboseLevel;
import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.RoutingMatrix;
import jline.lang.constant.SchedStrategy;
import jline.lang.constant.SolverType;
import jline.lang.nodes.Queue;
import jline.lang.processes.Exp;
import jline.solvers.SolverOptions;
import jline.solvers.mva.SolverMVA;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

/**
 * SolverNC's <code>config.multiserver</code> handling.
 *
 * <p>SolverNC represents a finite multiserver station in one of two ways:
 * Seidmann's approximation (demand L/c plus a delay L(c-1)/c, applied in
 * Solver_nc) or the exact load-dependent lattice mu(n)=min(n,c), which routes
 * the model to Solver_ncld. Which one it used was decided by the method name
 * alone -- Seidmann on <code>default</code>, the lattice on <code>exact</code>
 * -- with no way to ask for either, even though <code>config.multiserver</code>
 * exists on the shared SolverOptions and SolverMVA honours it.
 *
 * <p>The load-bearing property is that THE SHIPPED DEFAULT MOVES NO RESULT:
 * <code>default</code>, and an options object that never mentions the field,
 * must both reproduce the historical dispatch exactly; only an explicit
 * <code>seidmann</code> or <code>lld</code> changes anything.
 *
 * <p>The model is Bolch, Greiner, de Meer &amp; Trivedi Example 8.1: three FCFS
 * stations with 2, 3 and 1 servers and one closed class of 3 jobs, whose exact
 * load-dependent answer is known in closed form and is what
 * SolverMVA('exact') returns. Twin of
 * <code>python/tests/test_nc_multiserver_policy.py</code>.
 */
public class SolverNCMultiserverPolicyTest {

    /** Seidmann's approximation on Example 8.1, i.e. what the default path returns. */
    private static final double[] SEIDMANN_QLEN = {1.35037, 1.05217, 0.59745};
    /** The exact load-dependent product form, mu(n)=min(n,c). */
    private static final double[] EXACT_QLEN = {1.28981, 1.04988, 0.66031};

    private static final double TOL = 1e-4;

    /** Bolch Example 8.1: three FCFS stations, 2/3/1 servers, K=3. */
    private static Network multiserverModel() {
        Network model = new Network("Bolch_08_01");
        Queue node1 = new Queue(model, "Node1", SchedStrategy.FCFS);
        node1.setNumberOfServers(2);
        Queue node2 = new Queue(model, "Node2", SchedStrategy.FCFS);
        node2.setNumberOfServers(3);
        Queue node3 = new Queue(model, "Node3", SchedStrategy.FCFS);
        ClosedClass jobclass = new ClosedClass(model, "Class1", 3, node1);
        node1.setService(jobclass, new Exp(0.8));
        node2.setService(jobclass, new Exp(0.6));
        node3.setService(jobclass, new Exp(0.4));
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(jobclass, jobclass, node1, node1, 0.133);
        P.set(jobclass, jobclass, node1, node2, 0.667);
        P.set(jobclass, jobclass, node1, node3, 0.2);
        P.set(jobclass, jobclass, node2, node1, 1.0);
        P.set(jobclass, jobclass, node3, node1, 1.0);
        model.link(P);
        return model;
    }

    /** Bolch Example 7.5: the same shape with one server everywhere. */
    private static Network singleServerModel() {
        Network model = new Network("Bolch_07_05");
        Queue node1 = new Queue(model, "Node1", SchedStrategy.FCFS);
        Queue node2 = new Queue(model, "Node2", SchedStrategy.FCFS);
        Queue node3 = new Queue(model, "Node3", SchedStrategy.FCFS);
        ClosedClass jobclass = new ClosedClass(model, "Class1", 3, node1);
        node1.setService(jobclass, new Exp(0.8));
        node2.setService(jobclass, new Exp(0.6));
        node3.setService(jobclass, new Exp(0.4));
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(jobclass, jobclass, node1, node1, 0.6);
        P.set(jobclass, jobclass, node1, node2, 0.3);
        P.set(jobclass, jobclass, node1, node3, 0.1);
        P.set(jobclass, jobclass, node2, node1, 0.2);
        P.set(jobclass, jobclass, node2, node2, 0.3);
        P.set(jobclass, jobclass, node2, node3, 0.5);
        P.set(jobclass, jobclass, node3, node1, 0.4);
        P.set(jobclass, jobclass, node3, node3, 0.6);
        model.link(P);
        return model;
    }

    /** @param multiserver null leaves config.multiserver untouched. */
    private static Matrix qlen(Network model, String method, String multiserver) {
        SolverOptions o = new SolverOptions(SolverType.NC);
        o.verbose = VerboseLevel.SILENT;
        if (method != null) {
            o.method = method;
        }
        if (multiserver != null) {
            o.config.multiserver = multiserver;
        }
        return new SolverNC(model, o).getAvgQLen();
    }

    private static void assertQLen(double[] expected, Matrix got, String what) {
        for (int i = 0; i < expected.length; i++) {
            assertEquals(expected[i], got.get(i, 0), TOL, what + ", station " + (i + 1));
        }
    }

    /** The shipped default must not move: still Seidmann on 'default'. */
    @Test
    public void testDefaultPolicyIsUnchanged() {
        assertQLen(SEIDMANN_QLEN, qlen(multiserverModel(), null, null), "default, field unset");
        assertQLen(SEIDMANN_QLEN, qlen(multiserverModel(), null, "default"), "default, field 'default'");
    }

    /** 'exact' keeps converting to the load-dependent lattice. */
    @Test
    public void testExactMethodIsUnchanged() {
        assertQLen(EXACT_QLEN, qlen(multiserverModel(), "exact", null), "exact, field unset");
        assertQLen(EXACT_QLEN, qlen(multiserverModel(), "exact", "default"), "exact, field 'default'");
    }

    /** config.multiserver='lld' reaches the exact answer from the default method. */
    @Test
    public void testLldMakesTheDefaultExact() {
        assertQLen(EXACT_QLEN, qlen(multiserverModel(), null, "lld"), "default + lld");
    }

    /** config.multiserver='seidmann' keeps Seidmann even when 'exact' is asked for. */
    @Test
    public void testSeidmannOptsTheExactMethodOut() {
        assertQLen(SEIDMANN_QLEN, qlen(multiserverModel(), "exact", "seidmann"), "exact + seidmann");
    }

    /** The opt-in answer is SolverMVA's exact multiserver answer, not a near miss. */
    @Test
    public void testLldMatchesMvaExact() {
        Matrix nc = qlen(multiserverModel(), null, "lld");
        SolverOptions mvaOpt = new SolverOptions(SolverType.MVA);
        mvaOpt.method = "exact";
        mvaOpt.verbose = VerboseLevel.SILENT;
        Matrix mva = new SolverMVA(multiserverModel(), mvaOpt).getAvgQLen();
        for (int i = 0; i < 3; i++) {
            assertEquals(mva.get(i, 0), nc.get(i, 0), 1e-6, "NC lld vs MVA exact, station " + (i + 1));
        }
    }

    /** A SolverMVA-only rule warns and is not silently honoured as something else. */
    @Test
    public void testUnimplementedMvaValueFallsBackToDefault() {
        assertQLen(SEIDMANN_QLEN, qlen(multiserverModel(), null, "softmin"), "default + softmin");
        assertQLen(SEIDMANN_QLEN, qlen(multiserverModel(), null, "conway"), "default + conway");
    }

    /** With one server everywhere the policy has nothing to select, under any method. */
    @Test
    public void testSingleServerModelIsPolicyIndependent() {
        Matrix baseline = qlen(singleServerModel(), null, null);
        for (String method : new String[] {null, "exact"}) {
            for (String ms : new String[] {null, "default", "lld", "seidmann"}) {
                Matrix got = qlen(singleServerModel(), method, ms);
                for (int i = 0; i < 3; i++) {
                    assertEquals(baseline.get(i, 0), got.get(i, 0), 1e-9,
                            "single-server, method=" + method + " multiserver=" + ms
                                    + ", station " + (i + 1));
                }
            }
        }
    }
}
