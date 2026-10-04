/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.tr;

import java.util.ArrayList;
import java.util.List;

import jline.lang.Network;
import jline.lang.NetworkStruct;
import jline.solvers.SolverOptions;
import jline.util.matrix.Matrix;

/**
 * Solver-agnostic driver of a model TRANSFORMATION, the sibling of
 * {@link jline.solvers.fj.FJFixedPoint}.
 *
 * <p>Where the fork-join fixed point drives ONE transformation (MMT/HT), this
 * drives whichever one {@code options.config.transform} names: the strategy
 * rewrites the model into subproblems, each subproblem is solved by a REAL
 * solver, and the strategy maps the metrics back onto the original classes and
 * stations.
 *
 * <pre>
 * expand -&gt; [ for e in subproblems: solve(e); couple(e) ] -&gt; converged -&gt; lift
 * </pre>
 *
 * <p>SINGLE PASS BY DEFAULT: unless the expand phase declares itself iterated
 * the loop runs one sweep and lifts, so the common case never pays for a
 * convergence test it does not need.
 *
 * <p>THE INNER SOLVER IS THE OUTER SOLVER. The caller supplies {@link InnerSolve},
 * which constructs an instance of its own class, so a transformation written
 * once serves every solver rather than the one it was first written for.
 * {@code SolverCTMC.runChainAggregationAnalyzer} used to hard-wire
 * {@code new SolverCTMC(...)}, which is what kept a transform with nothing
 * CTMC-specific in it out of reach of MVA, NC and FLD.
 *
 * <p>COUPLING IS GAUSS-SEIDEL BY CONSTRUCTION: the couple step runs immediately
 * after each subproblem's solve, inside the sweep, so subproblem e+1 sees the
 * updated state of 1..e.
 *
 * <p>Mirrors MATLAB {@code @NetworkSolver/transformSolve.m} and python
 * {@code line_solver/solvers/transform_driver.py}.
 */
public final class TransformSolve {

    private TransformSolve() {
    }

    /** The inner solve, which the caller implements with its own solver class. */
    public interface InnerSolve {
        /**
         * Solves one subproblem.
         *
         * @param submodel the transformed network
         * @param options  the inner options, with the transform method name cleared
         * @return the metrics of that solve
         */
        Inner solve(Network submodel, SolverOptions options);
    }

    /**
     * One inner solve's answer, on all FOUR channels.
     *
     * <p>{@code getAvg} alone is not enough: its sixth output is the RESIDENCE
     * time rather than the system throughput, and the normalizing constant and
     * the reported method name live on the result object.
     */
    public static final class Inner {
        public final Matrix Q;
        public final Matrix U;
        public final Matrix R;
        public final Matrix T;
        /** The per-chain SYSTEM throughput, from getAvgSysTput. */
        public final Matrix X;
        /** log G, NaN when the inner solver does not produce one. */
        public final double lG;
        /** The method the inner solve reported, for compound method names. */
        public final String method;

        public Inner(Matrix Q, Matrix U, Matrix R, Matrix T, Matrix X, double lG, String method) {
            this.Q = Q;
            this.U = U;
            this.R = R;
            this.T = T;
            this.X = X;
            this.lG = lG;
            this.method = method;
        }
    }

    /**
     * What a strategy's expand phase produces: the subproblems, plus whatever
     * state its lift needs. A strategy subclasses this to carry its own context.
     */
    public static class Expanded {
        public final List<Network> submodels;
        /** False for a single-pass transformation, which is the default. */
        public final boolean iterated;

        public Expanded(List<Network> submodels, boolean iterated) {
            this.submodels = submodels;
            this.iterated = iterated;
        }
    }

    /** The lifted metrics, in ORIGINAL class and station coordinates. */
    public static final class Lifted {
        public final Matrix Q;
        public final Matrix U;
        public final Matrix R;
        public final Matrix T;
        public final Matrix C;
        public final Matrix X;

        public Lifted(Matrix Q, Matrix U, Matrix R, Matrix T, Matrix C, Matrix X) {
            this.Q = Q;
            this.U = U;
            this.R = R;
            this.T = T;
            this.C = C;
            this.X = X;
        }
    }

    /** A transformation, as a pair of phases. */
    public interface Strategy {
        /** Builds the subproblems and the context the lift will need. */
        Expanded expand(Network model, NetworkStruct sn, SolverOptions options);

        /** Maps the subproblem metrics back onto the original model. */
        Lifted lift(Expanded ctx, List<Inner> res);

        /**
         * Updates the subproblems from the sweep so far. Called only by an
         * iterated strategy, immediately after each subproblem's solve, which is
         * what makes the coupling Gauss-Seidel.
         */
        Expanded couple(Expanded ctx, List<Inner> res, int e);

        /** True when the sweep may stop. Called only by an iterated strategy. */
        boolean converged(Expanded ctx, List<Inner> res, int it);
    }

    /** What the driver answers with. */
    public static final class Result {
        public final Matrix Q;
        public final Matrix U;
        public final Matrix R;
        public final Matrix T;
        public final Matrix C;
        public final Matrix X;
        public final double lG;
        public final double runtime;
        public final int iter;
        public final String method;

        Result(Lifted l, double lG, double runtime, int iter, String method) {
            this.Q = l.Q;
            this.U = l.U;
            this.R = l.R;
            this.T = l.T;
            this.C = l.C;
            this.X = l.X;
            this.lG = lG;
            this.runtime = runtime;
            this.iter = iter;
            this.method = method;
        }
    }

    /** The canonical method name asked for, or {@link TransformMethod#NONE}. */
    public static String requested(SolverOptions options) {
        return TransformMethod.canonical(options.config.get("transform"));
    }

    /**
     * Whether a transformation was asked for at all.
     *
     * <p>Every solver that hosts a transformation calls this at the top of its
     * {@code runAnalyzer}, the counterpart of the branch MATLAB puts in its
     * shared {@code runAnalyzerPreamble}.
     */
    public static boolean isRequested(SolverOptions options) {
        return !TransformMethod.NONE.equals(requested(options));
    }

    /**
     * Runs the transformation named by {@code options.config.transform}.
     *
     * @param model       the original model
     * @param sn          its struct
     * @param options     the caller's options; the transform method name is read here
     * @param innerSolve  how a subproblem is solved, normally with the caller's
     *                    own solver class
     * @return the metrics in original coordinates
     */
    public static Result run(Network model, NetworkStruct sn, SolverOptions options,
                             InnerSolve innerSolve) {
        long t0 = System.nanoTime();
        String token = TransformMethod.canonical(options.config.get("transform"));
        if (TransformMethod.NONE.equals(token)) {
            throw new RuntimeException("TransformSolve.run called with transform='none'.");
        }
        int depth = TransformMethod.depth(options);
        if (depth > 0) {
            throw new RuntimeException("a model transformation ('" + token + "') cannot be nested "
                    + "inside another one; options.config.transform_depth is " + depth + ".");
        }
        Strategy strategy = TransformMethod.strategy(token);

        Expanded ctx = strategy.expand(model, sn, options);

        // The inner solve must not re-enter this driver. A KERNEL selection
        // inside the inner solver is not a transform and is correctly not cut.
        SolverOptions inner = options.copy();
        inner.config.put("transform", TransformMethod.NONE);
        inner.config.put("transform_depth", Integer.valueOf(depth + 1));

        List<Inner> res = new ArrayList<Inner>();
        for (int e = 0; e < ctx.submodels.size(); e++) {
            res.add(null);
        }
        int iters = 0;
        int iterMax = Math.max(1, options.iter_max);
        for (int it = 1; it <= iterMax; it++) {
            iters = it;
            for (int e = 0; e < ctx.submodels.size(); e++) {
                res.set(e, innerSolve.solve(ctx.submodels.get(e), inner));
                if (ctx.iterated) {
                    ctx = strategy.couple(ctx, res, e);
                }
            }
            if (!ctx.iterated) {
                break;
            }
            if (strategy.converged(ctx, res, it)) {
                break;
            }
        }

        Lifted lifted = strategy.lift(ctx, res);
        double runtime = (System.nanoTime() - t0) / 1000000000.0;
        return new Result(lifted, res.get(0).lG, runtime, iters, token);
    }
}
