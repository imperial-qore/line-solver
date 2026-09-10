/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.npfqn;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.optim.MaxIter;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.linear.LinearConstraint;
import org.apache.commons.math3.optim.linear.LinearConstraintSet;
import org.apache.commons.math3.optim.linear.LinearObjectiveFunction;
import org.apache.commons.math3.optim.linear.NonNegativeConstraint;
import org.apache.commons.math3.optim.linear.Relationship;
import org.apache.commons.math3.optim.linear.SimplexSolver;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;

import jline.util.matrix.Matrix;

/**
 * First-order linear-programming relaxation of the achievable region of a
 * multiclass open Markovian queueing network.
 *
 * <p>Port of matlab/src/api/npfqn/npfqn_bnd_bpt.m. Returns a LOWER bound on
 * {@code sum_r c_r x_r}, where {@code x_r} is the mean sojourn time of class
 * {@code r}, valid for EVERY non-idling scheduling policy.</p>
 *
 * <p>A "class" here is a buffer with its own exponential service rate and its
 * own Markovian routing, so a station serving several customer types owns one
 * class per type. The network is open: class {@code r} receives external
 * Poisson arrivals at rate {@code lambda0[r]} and, on completing service,
 * becomes class {@code r'} with probability {@code P(r,r')} or leaves with the
 * row deficit.</p>
 *
 * <p>METHOD. Uniformize the chain and let {@code R(t) = sum_r f(r) n_r(t)} for
 * an arbitrary vector {@code f}. The steady-state balance of {@code E[R^2]} is
 * an identity quadratic in {@code f}; since it holds for every {@code f}, the
 * coefficient matrices of the two sides agree entrywise. Diagonal entries give
 * one equation per class, off-diagonal entries one per unordered pair, in the
 * variables {@code x_r}, {@code I(r,l) = E[1{sigma(r) busy with r} n_l]} and
 * {@code N(i,l) = E[1{station i idle} n_l]}. A third block states that the
 * events "station i serves class r" and "station i idle" are mutually
 * exclusive and exhaustive, so their terms sum to {@code E[n_l] = lambda_l
 * x_l}. Minimizing over this polyhedron is a relaxation, hence a lower
 * bound.</p>
 *
 * <p>EXACT ON M/M/1: the LP reduces to {@code mu*I11 - lambda^2*x = lambda} and
 * {@code I11 + N11 = lambda*x} with {@code N11 >= 0}, whence
 * {@code x >= 1/(mu-lambda)} with equality.</p>
 *
 * <p>NOT INCLUDED, DELIBERATELY: the valid inequality {@code I(r,r) >= rho_r}
 * would tighten the relaxation but is not part of the reference's
 * characterization, and reproducing the reference's published bounds is the
 * acceptance test here.</p>
 *
 * <p>Reference: D. Bertsimas, I. Paschalidis, J. Tsitsiklis (1994).
 * Optimization of multiclass queueing networks: polyhedral and nonlinear
 * characterizations of achievable performance. Annals of Applied Probability
 * 4(1), 43-75. See also D. Bertsimas (1995), Queueing Systems 21, 337-389,
 * Theorem 9, which restates the same characterization.</p>
 */
public final class Npfqn_bnd_bpt {

    private Npfqn_bnd_bpt() {}

    /** Outcome of one solve. */
    public static class Result {
        /** Lower bound on {@code sum_r c_r x_r}. */
        public double zlb;
        /** The {@code x} block of the LP optimizer, length K. */
        public double[] x;
        /** Effective arrival rate of each class. */
        public double[] lambda;
        /** Per-class utilization {@code lambda_r/mu_r}. */
        public double[] rho;
        /** Per-station utilization. */
        public double[] rhoStation;
        /** Number of LP variables. */
        public int nvars;
        /** Number of LP rows. */
        public int nrows;
    }

    /**
     * @param lambda0   external Poisson arrival rate into each class (0 if none)
     * @param mu        exponential service rate of each class
     * @param P         K x K routing, P(r,r') = P(class r becomes r' after service)
     * @param stationOf zero-based station index of each class
     * @param c         objective weights; null means all ones
     */
    public static Result npfqn_bnd_bpt(double[] lambda0, double[] mu, Matrix P,
                                       int[] stationOf, double[] c) {
        final int K = lambda0.length;
        if (mu.length != K || stationOf.length != K) {
            throw new IllegalArgumentException("lambda0, mu and stationOf must all have " + K + " entries.");
        }
        if (P.getNumRows() != K || P.getNumCols() != K) {
            throw new IllegalArgumentException("P must be " + K + "x" + K + ".");
        }
        double[] cost = c;
        if (cost == null) {
            cost = new double[K];
            for (int r = 0; r < K; r++) {
                cost[r] = 1.0;
            }
        }
        if (cost.length != K) {
            throw new IllegalArgumentException("c must have " + K + " entries.");
        }
        for (int r = 0; r < K; r++) {
            if (!(mu[r] > 0)) {
                throw new IllegalArgumentException("Every class needs a strictly positive service rate.");
            }
            double rowSum = 0;
            for (int s = 0; s < K; s++) {
                rowSum += P.get(r, s);
            }
            if (rowSum > 1 + 1e-9) {
                throw new IllegalArgumentException("The routing matrix has a row summing above one.");
            }
        }
        int M = 0;
        for (int r = 0; r < K; r++) {
            M = Math.max(M, stationOf[r] + 1);
        }

        // ---- traffic equations, lambda = lambda0 + P' lambda ----
        Matrix ImPt = new Matrix(K, K);
        for (int i = 0; i < K; i++) {
            for (int j = 0; j < K; j++) {
                ImPt.set(i, j, (i == j ? 1.0 : 0.0) - P.get(j, i));
            }
        }
        Matrix rhsL = new Matrix(K, 1);
        for (int r = 0; r < K; r++) {
            rhsL.set(r, 0, lambda0[r]);
        }
        Matrix lamM = ImPt.inv().mult(rhsL);
        double[] lambda = new double[K];
        double[] rho = new double[K];
        double[] rhoStation = new double[M];
        for (int r = 0; r < K; r++) {
            lambda[r] = Math.max(lamM.get(r, 0), 0.0);
            rho[r] = lambda[r] / mu[r];
            rhoStation[stationOf[r]] += rho[r];
        }
        for (int i = 0; i < M; i++) {
            if (rhoStation[i] >= 1 - 1e-12) {
                throw new RuntimeException("Station " + (i + 1) + " is saturated (rho=" + rhoStation[i]
                        + "): no policy stabilizes the network.");
            }
        }

        // ---- variable layout: x(r) | I(r,l) | N(i,l) ----
        final int oI = K;
        final int oN = K + K * K;
        final int nv = K + K * K + M * K;

        List<LinearConstraint> rows = new ArrayList<LinearConstraint>();

        // (a) diagonal equations, test function n_r^2
        for (int r = 0; r < K; r++) {
            double[] a = new double[nv];
            a[oI + r * K + r] += 2 * mu[r];
            for (int w = 0; w < K; w++) {
                if (P.get(w, r) != 0) {
                    a[oI + w * K + r] += -2 * mu[w] * P.get(w, r);
                }
            }
            a[r] += -2 * lambda0[r] * lambda[r];
            rows.add(new LinearConstraint(a, Relationship.EQ, 2 * lambda[r] * (1 - P.get(r, r))));
        }

        // (b) off-diagonal equations, test function n_r*n_s
        for (int r = 1; r < K; r++) {
            for (int s = 0; s < r; s++) {
                double[] a = new double[nv];
                a[oI + r * K + s] += mu[r];
                a[oI + s * K + r] += mu[s];
                for (int w = 0; w < K; w++) {
                    if (P.get(w, r) != 0) {
                        a[oI + w * K + s] += -mu[w] * P.get(w, r);
                    }
                    if (P.get(w, s) != 0) {
                        a[oI + w * K + r] += -mu[w] * P.get(w, s);
                    }
                }
                a[s] += -lambda0[r] * lambda[s];
                a[r] += -lambda0[s] * lambda[r];
                rows.add(new LinearConstraint(a, Relationship.EQ,
                        -lambda[r] * P.get(r, s) - lambda[s] * P.get(s, r)));
            }
        }

        // (c) exhaustiveness at each station
        for (int i = 0; i < M; i++) {
            for (int l = 0; l < K; l++) {
                double[] a = new double[nv];
                for (int r = 0; r < K; r++) {
                    if (stationOf[r] == i) {
                        a[oI + r * K + l] += 1.0;
                    }
                }
                a[oN + i * K + l] += 1.0;
                a[l] += -lambda[l];
                rows.add(new LinearConstraint(a, Relationship.EQ, 0.0));
            }
        }

        double[] obj = new double[nv];
        for (int r = 0; r < K; r++) {
            obj[r] = cost[r];
        }
        LinearObjectiveFunction f = new LinearObjectiveFunction(obj, 0.0);
        SimplexSolver solver = new SimplexSolver();
        PointValuePair sol = solver.optimize(new MaxIter(100000), f,
                new LinearConstraintSet(rows), GoalType.MINIMIZE, new NonNegativeConstraint(true));
        if (sol == null) {
            throw new RuntimeException("The achievable-region LP did not solve to optimality.");
        }

        Result out = new Result();
        out.zlb = sol.getValue();
        out.x = new double[K];
        double[] pt = sol.getPoint();
        for (int r = 0; r < K; r++) {
            out.x[r] = pt[r];
        }
        out.lambda = lambda;
        out.rho = rho;
        out.rhoStation = rhoStation;
        out.nvars = nv;
        out.nrows = rows.size();
        return out;
    }
}
