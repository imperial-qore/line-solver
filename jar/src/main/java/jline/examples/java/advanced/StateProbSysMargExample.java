/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.examples.java.advanced;

import java.util.ArrayList;
import java.util.List;

import jline.io.Ret.ProbabilityResult;
import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.RoutingMatrix;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.processes.Exp;
import jline.solvers.ctmc.SolverCTMC;
import jline.solvers.nc.SolverNC;
import jline.util.matrix.Matrix;

/**
 * Joint probability of the per-station TOTAL queue lengths, all classes summed
 * out, from {@link SolverNC#getProbSysMarg(Matrix)}.
 *
 * <p>Compare with {@code getProbSysAggr}, which fixes the PER-CLASS population
 * of every station: each probability here is the sum of that one over every
 * per-class table with these row sums. The fibre grows combinatorially, so the
 * quantity is evaluated as a permanent of the demand matrix replicated once per
 * job (H. J. Ryser, "Combinatorial Mathematics", MAA 1963) rather than by
 * enumerating it.</p>
 *
 * <p>The law is exact, so its first moments are the mean queue lengths and the
 * example checks them against {@link SolverCTMC}. The approximate permanent
 * engines are also exercised: they trade accuracy for cost on models whose class
 * count makes the exact expansion dear, need a demand matrix with full support,
 * and refuse a structural zero rather than flooring it.</p>
 */
public class StateProbSysMargExample {

    public static void main(String[] args) {
        int[] N = {2, 1};

        Network model = new Network("model");
        Delay delay = new Delay(model, "Delay");
        Queue queue1 = new Queue(model, "Queue1", SchedStrategy.PS);
        Queue queue2 = new Queue(model, "Queue2", SchedStrategy.PS);

        ClosedClass class1 = new ClosedClass(model, "Class1", N[0], delay, 0);
        ClosedClass class2 = new ClosedClass(model, "Class2", N[1], delay, 0);

        delay.setService(class1, Exp.fitMean(1.5));
        delay.setService(class2, Exp.fitMean(2.0));
        queue1.setService(class1, Exp.fitMean(0.7));
        queue1.setService(class2, Exp.fitMean(0.4));
        queue2.setService(class1, Exp.fitMean(0.3));
        queue2.setService(class2, Exp.fitMean(0.9));

        RoutingMatrix P = model.initRoutingMatrix();
        P.set(class1, class1, delay, queue1, 1.0);
        P.set(class1, class1, queue1, queue2, 1.0);
        P.set(class1, class1, queue2, delay, 1.0);
        P.set(class2, class2, delay, queue1, 1.0);
        P.set(class2, class2, queue1, queue2, 1.0);
        P.set(class2, class2, queue2, delay, 1.0);
        model.link(P);

        SolverNC solver = new SolverNC(model);

        // Every way of splitting the closed population across the stations
        int M = 3, total = N[0] + N[1];
        List<int[]> states = compositions(M, total);

        System.out.printf("%n  n(Delay) n(Queue1) n(Queue2)        P(n)%n");
        double[] p = new double[states.size()];
        double sum = 0;
        double[] mean = new double[M];
        for (int j = 0; j < states.size(); j++) {
            int[] n = states.get(j);
            ProbabilityResult pr = solver.getProbSysMarg(rowVector(n));
            p[j] = pr.probability.get(0);
            sum += p[j];
            for (int i = 0; i < M; i++) {
                mean[i] += n[i] * p[j];
            }
            System.out.printf("  %8d %9d %9d  %10.6f%n", n[0], n[1], n[2], p[j]);
        }
        System.out.printf("  ------------------------------------------%n");
        System.out.printf("  sum %39.6f%n", sum);

        // The law is exact, so its first moments are the queue lengths
        Matrix qlen = new SolverCTMC(model).getAvgQLen();
        System.out.printf("%n  E[n] from the joint law : ");
        for (int i = 0; i < M; i++) {
            System.out.printf("%12.8f", mean[i]);
        }
        System.out.printf("%n  QLen from SolverCTMC    : ");
        for (int i = 0; i < M; i++) {
            double s = 0;
            for (int r = 0; r < qlen.getNumCols(); r++) {
                s += qlen.get(i, r);
            }
            System.out.printf("%12.8f", s);
        }
        System.out.println();

        double[] pb = new double[states.size()];
        double sb = 0;
        for (int j = 0; j < states.size(); j++) {
            pb[j] = solver.getProbSysMarg(rowVector(states.get(j)), "bethe").probability.get(0);
            sb += pb[j];
        }
        double err = 0;
        for (int j = 0; j < states.size(); j++) {
            err += Math.abs(pb[j] / sb - p[j]) / p[j];
        }
        System.out.printf("  Bethe engine, mean relative error : %.2f%%%n",
                100.0 * err / states.size());

        System.out.println(solver.citations());
    }

    /** Every composition of {@code total} into {@code M} nonnegative parts. */
    private static List<int[]> compositions(int M, int total) {
        List<int[]> out = new ArrayList<int[]>();
        compositionsRec(M, total, new int[M], 0, out);
        return out;
    }

    private static void compositionsRec(int M, int left, int[] acc, int pos, List<int[]> out) {
        if (pos == M - 1) {
            acc[pos] = left;
            out.add(acc.clone());
            return;
        }
        for (int k = 0; k <= left; k++) {
            acc[pos] = k;
            compositionsRec(M, left - k, acc, pos + 1, out);
        }
    }

    private static Matrix rowVector(int[] v) {
        Matrix m = new Matrix(1, v.length);
        for (int i = 0; i < v.length; i++) {
            m.set(0, i, v[i]);
        }
        return m;
    }
}
