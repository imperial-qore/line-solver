/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.examples.java.advanced;

import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.RoutingMatrix;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.processes.Exp;
import jline.solvers.NetworkSolver;
import jline.solvers.ctmc.SolverCTMC;
import jline.solvers.ldes.LDESOptions;
import jline.solvers.ldes.SolverLDES;
import jline.solvers.mva.SolverMVA;
import jline.util.matrix.Matrix;

/**
 * Demo of the LDES warm-start feature: an auxiliary solver is passed to SolverLDES
 * as an argument, its steady-state distribution is computed, and that distribution
 * decides the initial state of the simulation. Because the simulation then starts in
 * (approximately) steady state, the initialization bias vanishes and no warmup
 * samples are discarded, so a target accuracy is reached with fewer simulated events
 * than a cold-started run.
 *
 * Auxiliary-solver dispatch inside SolverLDES.initFromSolver:
 *  - SolverCTMC: the exact stationary distribution over the aggregate state space is
 *    computed and the initial state is its mode (most probable aggregate state);
 *  - any other solver (e.g. SolverMVA): the steady-state mean queue lengths are
 *    rounded to an integer placement that conserves the closed populations.
 *
 * The benchmark model is a closed near-balanced tandem, Think(Exp,1) -> Queue1(Exp,1.0)
 * -> Queue2(Exp,0.98). Near-balanced closed networks mix slowly: the split of jobs
 * between the two queues drifts on a long time scale, so the bias of the default
 * cold start (all jobs at the reference station) persists beyond what the MSER-5
 * transient filter can remove.
 *
 * Part A (sample efficiency, N=40): model small enough for SolverCTMC, comparing the
 * cold start against warm starts from the CTMC stationary-distribution mode and from
 * the rounded MVA mean queue lengths.
 *
 * Part B (wall-clock speedup, N=100): the auxiliary solver is exact MVA
 * (milliseconds), so the reduction in required samples translates directly into a
 * wall-clock speedup, auxiliary solver time included.
 */
public class LDESWarmStartExample {

    private static final double TOL = 0.10;   // 10% target relative error
    private static final int[] SEEDS = {23000, 23001, 23002, 23003, 23004};

    /** Closed near-balanced tandem: Think(Exp,1) -> Queue1(Exp,1.0) -> Queue2(Exp,0.98). */
    public static Network buildModel(int njobs) {
        Network model = new Network("ldesWarmStart");
        Delay think = new Delay(model, "Think");
        Queue queue1 = new Queue(model, "Queue1", SchedStrategy.FCFS);
        Queue queue2 = new Queue(model, "Queue2", SchedStrategy.FCFS);

        ClosedClass jobs = new ClosedClass(model, "Jobs", njobs, think);
        think.setService(jobs, new Exp(1.0));
        queue1.setService(jobs, new Exp(1.0));
        queue2.setService(jobs, new Exp(0.98));  // near-balanced bottleneck

        RoutingMatrix routing = model.initRoutingMatrix();
        routing.set(jobs, jobs, think, queue1, 1.0);
        routing.set(jobs, jobs, queue1, queue2, 1.0);
        routing.set(jobs, jobs, queue2, think, 1.0);
        model.link(routing);
        return model;
    }

    /** L1 relative error of the simulated station mean queue lengths vs the exact ones. */
    private static double relError(Matrix exact, Matrix sim) {
        double num = 0.0;
        double den = 0.0;
        for (int i = 0; i < exact.getNumRows(); i++) {
            num += Math.abs(sim.get(i, 0) - exact.get(i, 0));
            den += exact.get(i, 0);
        }
        return num / den;
    }

    /**
     * Run LDES over the seed set at a given sample budget and return {avgRelError,
     * totalRuntimeSec}. When initSol is non-null the simulation is warm-started from
     * that placement (steady-state initial state, no warmup discarded), reusing the
     * placement computed once by the auxiliary solver.
     */
    private static double[] runAtBudget(int njobs, Matrix exactQLen, int samples, Matrix initSol) {
        double errSum = 0.0;
        double tSum = 0.0;
        for (int s = 0; s < SEEDS.length; s++) {
            Network model = buildModel(njobs);
            long t0 = System.nanoTime();
            SolverLDES ldes = new SolverLDES(model, "samples", samples, "seed", SEEDS[s], "verbose", false);
            if (initSol != null) {
                LDESOptions ldesOptions = (LDESOptions) ldes.getOptions();
                ldesOptions.init_sol = initSol;
                ldesOptions.tranfilter = "fixed";
                ldesOptions.warmupfrac = 0.0;
            }
            Matrix simQLen = ldes.getAvgQLen();
            tSum += (System.nanoTime() - t0) / 1.0e9;
            errSum += relError(exactQLen, simQLen);
        }
        return new double[]{errSum / SEEDS.length, tSum};
    }

    /**
     * Increase the sample budget along the grid until TOL is met; return {samples,
     * avgError, cumulativeRuntimeSec}. initTime is the auxiliary solver cost, paid
     * once upfront (zero for the cold start).
     */
    private static double[] timeToAccuracy(int njobs, Matrix exactQLen, int[] grid,
                                           Matrix initSol, double initTime) {
        double cumTime = initTime;
        for (int g = 0; g < grid.length; g++) {
            double[] res = runAtBudget(njobs, exactQLen, grid[g], initSol);
            cumTime += res[1];
            System.out.printf("  samples=%8d  avg rel.err=%6.2f%%  batch time=%6.2fs%n",
                    grid[g], 100.0 * res[0], res[1]);
            if (res[0] < TOL) {
                return new double[]{grid[g], res[0], cumTime};
            }
        }
        return new double[]{-1, Double.NaN, cumTime};
    }

    /**
     * Compute a warm-start placement by passing the auxiliary solver to SolverLDES,
     * returning {init_sol, initTimeSec}. This exercises the new
     * SolverLDES(model, initSolver, args) constructor.
     */
    private static Object[] warmPlacement(Network model, NetworkSolver initSolver) {
        long t0 = System.nanoTime();
        SolverLDES proto = new SolverLDES(model, initSolver, "verbose", false);
        Matrix initSol = proto.getOptions().init_sol;
        double initTime = (System.nanoTime() - t0) / 1.0e9;
        return new Object[]{initSol, Double.valueOf(initTime)};
    }

    private static void report(String label, double[] res) {
        if (res[0] < 0) {
            System.out.printf("  %s: target not reached within grid (%.2fs spent)%n", label, res[2]);
        } else {
            System.out.printf("  %s: samples=%d, time=%.2fs%n", label, (int) res[0], res[2]);
        }
    }

    private static String rowToString(Matrix row) {
        StringBuilder sb = new StringBuilder("[");
        for (int j = 0; j < row.getNumCols(); j++) {
            if (j > 0) sb.append(' ');
            sb.append((int) Math.round(row.get(0, j)));
        }
        return sb.append(']').toString();
    }

    public static void main(String[] args) {
        // ================= Part A: sample efficiency (N=40, CTMC feasible) =================
        int nA = 40;
        System.out.println("PART A: sample efficiency, near-balanced tandem, N=" + nA);
        System.out.println();

        Matrix exactA = new SolverCTMC(buildModel(nA), "verbose", false).getAvgQLen();
        System.out.println("Exact CTMC mean queue lengths: " + rowToString(exactA.transpose()));

        Network m1 = buildModel(nA);
        Object[] wMva = warmPlacement(m1, new SolverMVA(m1, "exact", "verbose", false));
        Network m2 = buildModel(nA);
        Object[] wCtmc = warmPlacement(m2, new SolverCTMC(m2, "verbose", false));
        System.out.printf("Warm placement from SolverMVA  (%.3fs, rounded mean qlen): %s%n",
                (Double) wMva[1], rowToString((Matrix) wMva[0]));
        System.out.printf("Warm placement from SolverCTMC (%.3fs, distribution mode): %s%n",
                (Double) wCtmc[1], rowToString((Matrix) wCtmc[0]));
        System.out.println();

        int[] gridA = {5000, 10000, 20000, 50000, 100000, 200000};
        System.out.println("COLD start (default init, MSER-5 warmup removal):");
        double[] coldA = timeToAccuracy(nA, exactA, gridA, null, 0.0);
        System.out.println("WARM-MVA start (init from MVA mean queue lengths, no warmup):");
        double[] warmMvaA = timeToAccuracy(nA, exactA, gridA, (Matrix) wMva[0], (Double) wMva[1]);
        System.out.println("WARM-CTMC start (init from mode of exact stationary distribution, no warmup):");
        double[] warmCtmcA = timeToAccuracy(nA, exactA, gridA, (Matrix) wCtmc[0], (Double) wCtmc[1]);
        System.out.println();
        System.out.println("Samples needed for " + (int) (100 * TOL) + "% relative error (Part A):");
        report("COLD     ", coldA);
        report("WARM-MVA ", warmMvaA);
        report("WARM-CTMC", warmCtmcA);
        if (coldA[0] > 0 && warmCtmcA[0] > 0) {
            System.out.printf("  Sample reduction WARM-CTMC vs COLD: %.1fx%n", coldA[0] / warmCtmcA[0]);
        }
        System.out.println();

        // ================= Part B: wall-clock speedup (N=100, MVA warm start) =================
        int nB = 100;
        System.out.println("PART B: wall-clock speedup, near-balanced tandem, N=" + nB);
        System.out.println();

        Matrix exactB = new SolverMVA(buildModel(nB), "exact", "verbose", false).getAvgQLen();
        System.out.println("Exact MVA mean queue lengths: " + rowToString(exactB.transpose()));

        Network m3 = buildModel(nB);
        Object[] wMvaB = warmPlacement(m3, new SolverMVA(m3, "exact", "verbose", false));
        System.out.printf("Warm placement from SolverMVA (%.3fs): %s%n",
                (Double) wMvaB[1], rowToString((Matrix) wMvaB[0]));
        System.out.println();

        int[] gridB = {10000, 20000, 50000, 100000, 200000, 500000, 1000000};
        System.out.println("COLD start (default init, MSER-5 warmup removal):");
        double[] coldB = timeToAccuracy(nB, exactB, gridB, null, 0.0);
        System.out.println("WARM-MVA start (init from MVA mean queue lengths, no warmup):");
        double[] warmB = timeToAccuracy(nB, exactB, gridB, (Matrix) wMvaB[0], (Double) wMvaB[1]);
        System.out.println();
        System.out.println("Time to reach " + (int) (100 * TOL) + "% relative error (Part B):");
        report("COLD    ", coldB);
        report("WARM-MVA", warmB);
        if (coldB[2] > 0 && warmB[2] > 0 && coldB[0] > 0 && warmB[0] > 0) {
            System.out.printf("  Speedup WARM-MVA vs COLD: %.2fx (%.0fx fewer samples)%n",
                    coldB[2] / warmB[2], coldB[0] / warmB[0]);
        }
    }
}
