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

/**
 * Piecewise-linear Lyapunov UPPER bound on the steady-state queue lengths of a
 * multitype (deterministic-routing) multiclass Markovian queueing network,
 * valid for EVERY work-conserving Markovian policy.
 *
 * <p>Port of matlab/src/api/npfqn/npfqn_bnd_bgt.m.</p>
 *
 * <p>MODEL. J single-server stations; I customer types; type i arrives as a
 * Poisson stream of rate {@code lambda[i]} and passes through stages
 * k = 0..mu[i].length-1, stage k being served at station {@code sigma[i][k]}
 * at exponential rate {@code mu[i][k]}. Class (i,k) is the buffer of type i at
 * stage k; N = sum_i mu[i].length is the number of classes.</p>
 *
 * <p>METHOD. Solve the Down-Meyn global-stability linear program GLP[dm], eq.
 * (25)-(28) of the reference, in the piecewise-linear Lyapunov function
 * {@code Phi(x) = max_j L^j'x}:</p>
 *
 * <pre>
 *   L^j(i,1) lambda_i + mu(i,k) (L^j(i,k+1) - L^j(i,k)) + V_j &lt;= -gamma
 *                                                    for (i,k) in station j
 *   mu(i,k) (L^j(i,k+1) - L^j(i,k)) &lt;= V_j           for (i,k) not in j
 *   (1/(J-1)) sum_{j' != j} L^j'(i,k) &gt;= L^j(i,k)    for (i,k) not in j
 *   L, V, gamma &gt;= 0
 * </pre>
 *
 * <p>with {@code L^j(i,Ji+1) = 0}. A feasible solution with gamma &gt; 0
 * certifies that EVERY work-conserving policy is stable, and a smoothed Phi is
 * then a Lyapunov function with drift gamma/4 and an explicit exception
 * parameter, giving the reference's Theorem 4 bound</p>
 *
 * <pre>
 *   E[L^j'Q] &lt;= 16 N J^2 (J-1) (Lmax+gamma)^3 / gamma^2
 *               + 8 (Lmax + gamma/2)^2 / gamma  =: U
 * </pre>
 *
 * <p>for every j, whence {@code E[Q(i,k)] &lt;= U / max_j L^j(i,k)}.</p>
 *
 * <p>THE RATES ARE RESCALED so that {@code sum_i lambda_i + sum_{i,k} mu(i,k)
 * = 1}, the uniformization the reference imposes before Theorem 4. Queue
 * lengths are counts and are unaffected by the time scale.</p>
 *
 * <p>NORMALIZATION, WHICH THE REFERENCE LEAVES OPEN. GLP[dm] is homogeneous and
 * so is the bound, so this routine fixes {@code L^j(i,k) &lt;= 1} and MAXIMIZES
 * gamma, then breaks ties among gamma-optimal solutions by maximizing sum L: a
 * degenerate optimum can otherwise zero some {@code L^j(i,k)} and report an
 * infinite bound for a class for no reason.</p>
 *
 * <p>THE BOUND IS LOOSE, and knowingly so: the exception parameter carries
 * {@code (Lmax+gamma)^3/gamma^2} and dominates as soon as J &gt; 1. What is
 * sharp is the STABILITY CERTIFICATE gamma &gt; 0 and the geometric tail
 * RATE.</p>
 *
 * <p>Reference: D. Bertsimas, D. Gamarnik, J. N. Tsitsiklis (2001). Performance
 * of multiclass Markovian queueing networks via piecewise linear Lyapunov
 * functions. Annals of Applied Probability 11(4), 1384-1428, Section 5.1
 * (GLP[dm] of Down and Meyn 1997, and Theorem 4).</p>
 */
public final class Npfqn_bnd_bgt {

    private Npfqn_bnd_bgt() {}

    /** Outcome of one solve. */
    public static class Result {
        /** Per type and stage, the upper bound on E[Q(i,k)]; +Inf where max_j L is 0. */
        public double[][] Qub;
        /** The drift certificate; strictly positive on success. */
        public double gamma;
        /** max over j and (i,k) of L. */
        public double Lmax;
        /** The Lyapunov coefficients, [J][N]. */
        public double[][] L;
        /** The per-station slack V_j. */
        public double[] V;
        /** The exception parameter of the smoothed Lyapunov function. */
        public double B;
        /** The Theorem 4 bound on E[L^j'Q], the same for every j. */
        public double U;
        /** Geometric decay ratio of the tail bound. */
        public double tailRatio;
        /** Step of the tail bound, 2(Lmax+gamma/2). */
        public double tailStep;
        /** Per-class nominal load. */
        public double[] rho;
        /** Per-station nominal load. */
        public double[] rhoStation;
        /** The uniformization divisor applied to lambda and mu. */
        public double scale;
        /** Per-class type index. */
        public int[] classType;
        /** Per-class stage index, zero based. */
        public int[] classStage;
        /** Per-class station index, zero based. */
        public int[] classStation;
    }

    /**
     * @param lambda Poisson arrival rate of each type
     * @param mu     mu[i][k] = service rate of stage k of type i
     * @param sigma  sigma[i][k] = zero-based station of stage k of type i
     * @param J      number of stations
     */
    public static Result npfqn_bnd_bgt(double[] lambda, double[][] mu, int[][] sigma, int J) {
        final int I = lambda.length;
        if (mu.length != I || sigma.length != I) {
            throw new IllegalArgumentException("mu and sigma must both have " + I + " entries.");
        }

        // ---- flatten (i,k) into a class index ----
        List<Integer> ctype = new ArrayList<Integer>();
        List<Integer> cstage = new ArrayList<Integer>();
        List<Integer> cstation = new ArrayList<Integer>();
        List<Double> cmu = new ArrayList<Double>();
        int[] firstOf = new int[I];
        for (int i = 0; i < I; i++) {
            if (mu[i].length != sigma[i].length) {
                throw new IllegalArgumentException("mu[" + i + "] and sigma[" + i
                        + "] have different lengths.");
            }
            if (mu[i].length == 0) {
                throw new IllegalArgumentException("Type " + (i + 1) + " has no stage.");
            }
            firstOf[i] = ctype.size();
            for (int k = 0; k < mu[i].length; k++) {
                ctype.add(i);
                cstage.add(k);
                cstation.add(sigma[i][k]);
                cmu.add(mu[i][k]);
            }
            if (!(lambda[i] > 0)) {
                throw new IllegalArgumentException("Every type needs a strictly positive arrival rate.");
            }
        }
        final int N = cmu.size();
        int[] classType = new int[N], classStage = new int[N], classStation = new int[N];
        double[] muc = new double[N];
        for (int c = 0; c < N; c++) {
            classType[c] = ctype.get(c);
            classStage[c] = cstage.get(c);
            classStation[c] = cstation.get(c);
            muc[c] = cmu.get(c);
            if (!(muc[c] > 0)) {
                throw new IllegalArgumentException("Every stage needs a strictly positive service rate.");
            }
        }
        int[] nextOf = new int[N];
        for (int c = 0; c < N; c++) {
            nextOf[c] = (c + 1 < N && classType[c + 1] == classType[c]) ? (c + 1) : -1;
        }

        // ---- loads ----
        double[] rho = new double[N];
        double[] rhoStation = new double[J];
        for (int c = 0; c < N; c++) {
            rho[c] = lambda[classType[c]] / muc[c];
            rhoStation[classStation[c]] += rho[c];
        }
        for (int j = 0; j < J; j++) {
            if (rhoStation[j] >= 1) {
                throw new RuntimeException("Station " + (j + 1) + " is saturated (rho="
                        + rhoStation[j] + "): the load condition of the reference fails.");
            }
        }

        // ---- uniformization ----
        double scale = 0;
        for (int i = 0; i < I; i++) {
            scale += lambda[i];
        }
        for (int c = 0; c < N; c++) {
            scale += muc[c];
        }
        double[] lam = new double[I];
        for (int i = 0; i < I; i++) {
            lam[i] = lambda[i] / scale;
        }
        double[] mus = new double[N];
        for (int c = 0; c < N; c++) {
            mus[c] = muc[c] / scale;
        }

        // ---- LP layout: L(j,c) -> j*N + c ; V(j) -> J*N + j ; gamma -> J*N+J ----
        final int oV = J * N;
        final int ig = J * N + J;
        final int nv = J * N + J + 1;

        List<LinearConstraint> rows = new ArrayList<LinearConstraint>();
        for (int j = 0; j < J; j++) {
            for (int c = 0; c < N; c++) {
                double[] a = new double[nv];
                if (classStation[c] == j) {
                    a[j * N + firstOf[classType[c]]] += lam[classType[c]];
                    a[j * N + c] += -mus[c];
                    if (nextOf[c] >= 0) {
                        a[j * N + nextOf[c]] += mus[c];
                    }
                    a[oV + j] += 1.0;
                    a[ig] += 1.0;
                    rows.add(new LinearConstraint(a, Relationship.LEQ, 0.0));
                } else {
                    a[j * N + c] += -mus[c];
                    if (nextOf[c] >= 0) {
                        a[j * N + nextOf[c]] += mus[c];
                    }
                    a[oV + j] += -1.0;
                    rows.add(new LinearConstraint(a, Relationship.LEQ, 0.0));
                    if (J > 1) {
                        double[] b = new double[nv];
                        b[j * N + c] += 1.0;
                        for (int jp = 0; jp < J; jp++) {
                            if (jp != j) {
                                b[jp * N + c] += -1.0 / (J - 1);
                            }
                        }
                        rows.add(new LinearConstraint(b, Relationship.LEQ, 0.0));
                    }
                }
            }
        }
        // the homogeneous normalization Lmax <= 1
        for (int v = 0; v < J * N; v++) {
            double[] a = new double[nv];
            a[v] = 1.0;
            rows.add(new LinearConstraint(a, Relationship.LEQ, 1.0));
        }

        double[] obj = new double[nv];
        obj[ig] = 1.0;
        SimplexSolver solver = new SimplexSolver();
        PointValuePair sol = solver.optimize(new MaxIter(200000),
                new LinearObjectiveFunction(obj, 0.0), new LinearConstraintSet(rows),
                GoalType.MAXIMIZE, new NonNegativeConstraint(true));
        if (sol == null) {
            throw new RuntimeException("GLP[dm] did not solve to optimality.");
        }
        double gamma = sol.getValue();
        if (!(gamma > 0)) {
            throw new RuntimeException("GLP[dm] has no solution with gamma > 0: this network is "
                    + "not certified globally stable, so no finite piecewise-linear Lyapunov "
                    + "bound exists.");
        }

        // Tie-break among gamma-optimal solutions: maximize sum L.
        double[] pt = sol.getPoint();
        List<LinearConstraint> rows2 = new ArrayList<LinearConstraint>(rows);
        double[] pin = new double[nv];
        pin[ig] = 1.0;
        rows2.add(new LinearConstraint(pin, Relationship.GEQ, gamma));
        double[] obj2 = new double[nv];
        for (int v = 0; v < J * N; v++) {
            obj2[v] = 1.0;
        }
        try {
            PointValuePair sol2 = solver.optimize(new MaxIter(200000),
                    new LinearObjectiveFunction(obj2, 0.0), new LinearConstraintSet(rows2),
                    GoalType.MAXIMIZE, new NonNegativeConstraint(true));
            if (sol2 != null) {
                pt = sol2.getPoint();
                gamma = pt[ig];
            }
        } catch (RuntimeException e) {
            // the tie-break is an improvement, not a requirement: keep stage one
        }

        Result out = new Result();
        out.L = new double[J][N];
        out.Lmax = 0;
        for (int j = 0; j < J; j++) {
            for (int c = 0; c < N; c++) {
                out.L[j][c] = pt[j * N + c];
                out.Lmax = Math.max(out.Lmax, out.L[j][c]);
            }
        }
        out.V = new double[J];
        for (int j = 0; j < J; j++) {
            out.V[j] = pt[oV + j];
        }
        out.gamma = gamma;
        out.B = 16.0 * N * J * J * (J - 1) * Math.pow(out.Lmax + gamma, 3) / (gamma * gamma);
        out.U = out.B + 8.0 * Math.pow(out.Lmax + gamma / 2, 2) / gamma;
        out.tailStep = 2 * (out.Lmax + gamma / 2);
        out.tailRatio = (out.Lmax + gamma / 2) / (out.Lmax + 0.75 * gamma);

        out.Qub = new double[I][];
        for (int i = 0; i < I; i++) {
            out.Qub[i] = new double[mu[i].length];
        }
        for (int c = 0; c < N; c++) {
            double best = 0;
            for (int j = 0; j < J; j++) {
                best = Math.max(best, out.L[j][c]);
            }
            out.Qub[classType[c]][classStage[c]] =
                    best > 0 ? out.U / best : Double.POSITIVE_INFINITY;
        }
        out.rho = rho;
        out.rhoStation = rhoStation;
        out.scale = scale;
        out.classType = classType;
        out.classStage = classStage;
        out.classStation = classStation;
        return out;
    }
}
