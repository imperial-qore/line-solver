/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.examples.java.advanced;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import jline.api.pfqn.Pfqn_sdr;
import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.StateDepRouting;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Node;
import jline.lang.nodes.Queue;
import jline.lang.processes.Exp;
import jline.solvers.ctmc.SolverCTMC;
import jline.util.matrix.Matrix;

/**
 * Multi-center branches under Krzesinski (1987) state-dependent routing.
 *
 * <p>Checks that the product form of eq. (16) extends to a branch holding
 * several centers, and that the coefficients xi are the branch traffic
 * equations rather than the Section 3.2 shorthand xi = xi_e, which is exact
 * only when the branch departure center is visited once.</p>
 *
 * <p>Case A is a plain series branch, where the two readings coincide. Case B
 * feeds the branch departure center back onto its entry center, where they do
 * not: the traffic-equation xi reproduces the exact CTMC to machine precision
 * while the literal xi = 1 is out by 2.98e-1 in Q. See
 * _kb/16-state-dependent-routing.md</p>
 */
public class SdrMultiBranchExample {

    public static void main(String[] args) {
        System.out.println();
        System.out.println("==== SDR with multi-center branches ====");
        runCase("A: branch 2 = 2a -> 2b (series)", 0.0);
        runCase("B: branch 2 = 2a -> 2b, 2b -> 2a w.p. 0.5 (feedback onto the branch departure)", 0.5);
    }

    /**
     * Builds and reports one case.
     *
     * @param label  what the case is
     * @param pback  probability that the branch departure center feeds back onto its entry center
     */
    private static void runCase(String label, double pback) {
        double[] mu = {1.0, 0.9, 0.7, 0.5};   // centers 1, 2a, 2b, 3
        int N = 3;

        Network model = new Network("sdr_multi");
        Queue cpu = new Queue(model, "CPU", SchedStrategy.FCFS);
        Queue b2a = new Queue(model, "B2a", SchedStrategy.FCFS);
        Queue b2b = new Queue(model, "B2b", SchedStrategy.FCFS);
        Queue b3 = new Queue(model, "B3", SchedStrategy.FCFS);
        ClosedClass jobclass = new ClosedClass(model, "Class1", N, cpu, 0);
        cpu.setService(jobclass, new Exp(mu[0]));
        b2a.setService(jobclass, new Exp(mu[1]));
        b2b.setService(jobclass, new Exp(mu[2]));
        b3.setService(jobclass, new Exp(mu[3]));

        model.addLink(cpu, cpu);      // denied entry: the busy form of waiting
        model.addLink(cpu, b2a);
        model.addLink(cpu, b3);
        model.addLink(b2a, b2b);
        model.addLink(b2b, cpu);
        model.addLink(b3, cpu);
        b2a.setProbRouting(jobclass, b2b, 1.0);
        if (pback > 0) {
            model.addLink(b2b, b2a);
            b2b.setProbRouting(jobclass, b2a, pback);
            b2b.setProbRouting(jobclass, cpu, 1.0 - pback);
        } else {
            b2b.setProbRouting(jobclass, cpu, 1.0);
        }
        b3.setProbRouting(jobclass, cpu, 1.0);

        List<List<Node>> branches = new ArrayList<List<Node>>();
        branches.add(new ArrayList<Node>());
        branches.add(new ArrayList<Node>(Arrays.asList((Node) b2a, (Node) b2b)));
        branches.add(new ArrayList<Node>(Arrays.asList((Node) b3)));
        double[][] d = {{0, 2, 2}, {0, 0, 2}};
        cpu.setStateDepRouting(jobclass, cpu, branches, new int[]{0, 1, 2}, new double[]{-1, -1}, d);

        StateDepRouting sdr = model.getStruct(true).sdr;

        // The state-INDEPENDENT part of the routing: the complement M-V and the
        // arcs inside a branch. The state-dependent arcs out of the entry center
        // are not part of it.
        Matrix P = new Matrix(4, 4);
        P.set(0, 0, 1.0);                     // complement M-V = {1}
        P.set(1, 2, 1.0);                     // inside branch 2
        if (pback > 0) {
            P.set(2, 1, pback);
        }
        List<Matrix> Pchain = new ArrayList<Matrix>();
        Pchain.add(P);
        Matrix xi = Pfqn_sdr.pfqn_sdrvisits(sdr, Pchain);

        Matrix S = new Matrix(4, 1);
        for (int i = 0; i < 4; i++) {
            S.set(i, 0, 1.0 / mu[i]);
        }
        Matrix pop = new Matrix(1, 1);
        pop.set(0, 0, N);
        Pfqn_sdr.Result res = Pfqn_sdr.pfqn_sdr(S, xi, pop, sdr, null);

        SolverCTMC ctmc = new SolverCTMC(model);
        Matrix Qc = ctmc.getAvgQLen();
        Matrix Xc = ctmc.getAvgTput();

        System.out.println();
        System.out.printf("  %s%n", label);
        System.out.printf("    xi          = %s%n", row(xi));
        System.out.printf("    pfqn_sdr Q  = %s%n", row(res.Q));
        System.out.printf("    CTMC     Q  = %s%n", row(Qc));
        System.out.printf("    pfqn_sdr X  = %s%n", row(res.X));
        System.out.printf("    CTMC     X  = %s%n", row(Xc));
        System.out.printf("    max|dQ| = %.3e   max|dX| = %.3e%n",
                maxAbsDiff(Qc, res.Q), maxAbsDiff(Xc, res.X));

        // the reading asserted verbatim by the paper for E+D: xi = xi_e everywhere
        Matrix xiflat = new Matrix(4, 1);
        for (int i = 0; i < 4; i++) {
            xiflat.set(i, 0, 1.0);
        }
        Pfqn_sdr.Result flat = Pfqn_sdr.pfqn_sdr(S, xiflat, pop, sdr, null);
        System.out.printf("    xi==1 everywhere: max|dQ| = %.3e   max|dX| = %.3e%n",
                maxAbsDiff(Qc, flat.Q), maxAbsDiff(Xc, flat.X));
    }

    /** One matrix printed as a row of six-figure numbers. */
    private static String row(Matrix m) {
        StringBuilder sb = new StringBuilder("[");
        for (int i = 0; i < m.getNumRows(); i++) {
            for (int j = 0; j < m.getNumCols(); j++) {
                if (sb.length() > 1) {
                    sb.append(" ");
                }
                sb.append(String.format("%.6g", m.get(i, j)));
            }
        }
        return sb.append("]").toString();
    }

    /** Largest elementwise absolute difference of two equally shaped matrices. */
    private static double maxAbsDiff(Matrix a, Matrix b) {
        double worst = 0.0;
        for (int i = 0; i < a.getNumRows(); i++) {
            for (int j = 0; j < a.getNumCols(); j++) {
                worst = Math.max(worst, Math.abs(a.get(i, j) - b.get(i, j)));
            }
        }
        return worst;
    }
}
