/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.fluid.analyzers;

import jline.GlobalConstants;
import jline.lang.NetworkStruct;
import jline.lang.constant.SchedStrategy;
import jline.solvers.SolverOptions;
import jline.solvers.SolverResult;
import jline.util.matrix.Matrix;
import org.apache.commons.math3.util.FastMath;

import java.util.Random;

import static jline.io.InputOutput.line_error;
import static jline.io.InputOutput.mfilename;

/**
 * Diffusion approximation for closed multiclass BCMP networks.
 *
 * <p>Queue lengths are modelled as continuous variables driven by Brownian noise
 * and integrated with the Euler-Maruyama scheme, which is the stochastic
 * extension of the deterministic fluid limit. Each step takes the flow drift
 *
 * <pre>
 *   d(i,r) = sum_(j,s) x(j,s)/mu_inv(j,s) * P[(j,s)-&gt;(i,r)]  -  x(i,r)/mu_inv(i,r)
 * </pre>
 *
 * adds sqrt(dt)*N(0,1) per coordinate, reflects at zero and projects each class
 * back onto its own population, which is the closed-network constraint the
 * unconstrained SDE does not preserve.
 *
 * <p>Supported models: closed multiclass networks, scheduling PS / FCFS / INF /
 * SIRO, single-server or infinite-server stations only.
 *
 * <p>IT IS A STOCHASTIC METHOD. Like SSA and LDES it draws from a seeded stream,
 * so a run is reproducible for a fixed {@code options.seed} but is NOT expected
 * to agree digit for digit with MATLAB, whose {@code randn} is a different
 * generator. Everything else -- the drift, the reflection, the projection and
 * the metric read-back -- is the reference's own.
 *
 * <p>Port of {@code matlab/src/solvers/FLD/solver_fluid_diffusion.m}.
 */
public class DiffusionAnalyzer implements FluidAnalyzer {

    private Matrix xvecIt;

    /** Default number of Euler-Maruyama steps, as in the reference. */
    private static final int DEFAULT_STEPS = 10000;
    /** Default integration step, as in the reference. */
    private static final double DEFAULT_DT = 0.01;

    @Override
    public void analyze(NetworkStruct sn, SolverOptions options, SolverResult result) {
        int M = sn.nstations;
        int K = sn.nclasses;

        // -- model validation: diffusion supports closed networks only ---------
        for (int r = 0; r < K; r++) {
            if (Double.isInfinite(sn.njobs.get(0, r))) {
                line_error(mfilename(new Object() {
                }), "Diffusion method only supports closed queueing networks (no open classes).");
            }
        }
        for (int i = 0; i < M; i++) {
            SchedStrategy sched = sn.sched.get(sn.stations.get(i));
            if (sched == SchedStrategy.EXT) {
                line_error(mfilename(new Object() {
                }), "Diffusion method does not support Source nodes. Use closed networks only.");
            }
            if (sched != SchedStrategy.PS && sched != SchedStrategy.FCFS
                    && sched != SchedStrategy.INF && sched != SchedStrategy.SIRO) {
                line_error(mfilename(new Object() {
                }), String.format("Diffusion method does not support scheduling strategy %s at station %d.",
                        sched, i + 1));
            }
            double c = sn.nservers.get(i, 0);
            if (c > 1 && !Double.isInfinite(c)) {
                line_error(mfilename(new Object() {
                }), String.format("Diffusion method only supports single-server (c=1) or infinite-server "
                        + "stations. Station %d has %g servers.", i + 1, c));
            }
        }

        int steps = (options.iter_max > 0) ? options.iter_max : DEFAULT_STEPS;
        double dt = DEFAULT_DT;

        // -- inverse service rates: Inf where the class is not served here -----
        double[][] muInv = new double[M][K];
        for (int i = 0; i < M; i++) {
            for (int r = 0; r < K; r++) {
                double rate = sn.rates.get(i, r);
                muInv[i][r] = (rate > 0 && !Double.isNaN(rate))
                        ? 1.0 / rate : Double.POSITIVE_INFINITY;
            }
        }

        // -- station-level routing, class-blocked ------------------------------
        // The reference re-indexes this into MATLAB's column-major sub2ind layout
        // and then decodes it back with ind2sub inside the drift, so the same
        // permutation is applied to rows and columns and cancels exactly; the
        // row-major (station-major) block used here is the same operator.
        double[][] P = new double[M * K][M * K];
        for (int i = 0; i < M; i++) {
            int isf = (int) sn.stationToStateful.get(i);
            for (int r = 0; r < K; r++) {
                int src = isf * K + r;
                for (int j = 0; j < M; j++) {
                    int jsf = (int) sn.stationToStateful.get(j);
                    for (int s = 0; s < K; s++) {
                        int dst = jsf * K + s;
                        P[i * K + r][j * K + s] = sn.rt.get(src, dst);
                    }
                }
            }
        }

        double[] Nvec = new double[K];
        for (int r = 0; r < K; r++) {
            Nvec[r] = sn.njobs.get(0, r);
        }

        // -- Euler-Maruyama ----------------------------------------------------
        Random rng = new Random(options.seed);
        double sqdt = FastMath.sqrt(dt);
        double[][] x = new double[M][K];
        double[][] xavg = new double[M][K];
        double[][][] xt = new double[steps][M][K];
        for (int r = 0; r < K; r++) {
            for (int i = 0; i < M; i++) {
                x[i][r] = Nvec[r] / M;
            }
        }
        for (int i = 0; i < M; i++) {
            for (int r = 0; r < K; r++) {
                xt[0][i][r] = x[i][r];
                xavg[i][r] = x[i][r] / steps;
            }
        }
        for (int step = 1; step < steps; step++) {
            double[][] drift = drift(x, muInv, P, M, K);
            double[][] xnew = new double[M][K];
            for (int i = 0; i < M; i++) {
                for (int r = 0; r < K; r++) {
                    double v = x[i][r] + drift[i][r] * dt + sqdt * rng.nextGaussian();
                    xnew[i][r] = FastMath.max(v, 0);
                }
            }
            // project each class back onto its own population: the closed-network
            // constraint the unconstrained SDE does not preserve
            for (int r = 0; r < K; r++) {
                double tot = 0;
                for (int i = 0; i < M; i++) {
                    tot += xnew[i][r];
                }
                if (tot > 0) {
                    for (int i = 0; i < M; i++) {
                        xnew[i][r] = xnew[i][r] * Nvec[r] / tot;
                    }
                } else {
                    for (int i = 0; i < M; i++) {
                        xnew[i][r] = Nvec[r] / M;
                    }
                }
            }
            x = xnew;
            for (int i = 0; i < M; i++) {
                for (int r = 0; r < K; r++) {
                    xt[step][i][r] = x[i][r];
                    xavg[i][r] += x[i][r] / steps;
                }
            }
        }

        // -- metrics -----------------------------------------------------------
        result.QN = new Matrix(M, K);
        result.UN = new Matrix(M, K);
        result.RN = new Matrix(M, K);
        result.TN = new Matrix(M, K);
        result.WN = new Matrix(0, 0);
        result.AN = new Matrix(0, 0);
        result.QNt = new Matrix[M][K];
        result.UNt = new Matrix[M][K];
        result.TNt = new Matrix[M][K];
        result.t = new Matrix(steps, 1);
        for (int step = 0; step < steps; step++) {
            result.t.set(step, 0, step * dt);
        }

        for (int i = 0; i < M; i++) {
            double c = sn.nservers.get(i, 0);
            for (int r = 0; r < K; r++) {
                double q = xavg[i][r];
                result.QN.set(i, r, q);

                double tn = 0;
                boolean served = muInv[i][r] > 0 && !Double.isInfinite(muInv[i][r]);
                if (served) {
                    // an infinite server clears q jobs at once; a single server
                    // clears at most one, which is Little's law read at the
                    // saturation bound
                    tn = (Double.isInfinite(c) ? q : FastMath.min(q, 1.0)) / muInv[i][r];
                }
                result.TN.set(i, r, tn);
                result.UN.set(i, r, Double.isInfinite(c) ? q : FastMath.min(q / c, 1.0));
                // See SolverFluid: TN is zero only to the integrator's accuracy.
                result.RN.set(i, r, tn > GlobalConstants.Zero ? q / tn : 0.0);

                result.QNt[i][r] = new Matrix(steps, 1);
                result.UNt[i][r] = new Matrix(steps, 1);
                result.TNt[i][r] = new Matrix(steps, 1);
                for (int step = 0; step < steps; step++) {
                    double qs = xt[step][i][r];
                    result.QNt[i][r].set(step, 0, qs);
                    result.UNt[i][r].set(step, 0,
                            Double.isInfinite(c) ? qs : FastMath.min(qs / c, 1.0));
                    result.TNt[i][r].set(step, 0, served
                            ? (Double.isInfinite(c) ? qs : FastMath.min(qs, 1.0)) / muInv[i][r]
                            : 0.0);
                }
            }
        }

        xvecIt = new Matrix(1, M * K);
        for (int i = 0; i < M; i++) {
            for (int r = 0; r < K; r++) {
                xvecIt.set(0, i * K + r, xavg[i][r]);
            }
        }
        result.iter = 1;
    }

    @Override
    public Matrix getXVecIt() {
        return xvecIt;
    }

    /** Net flow into each (station, class) pair. */
    private static double[][] drift(double[][] x, double[][] muInv, double[][] P, int M, int K) {
        double[] outflow = new double[M * K];
        for (int j = 0; j < M; j++) {
            for (int s = 0; s < K; s++) {
                double mi = muInv[j][s];
                outflow[j * K + s] = (mi > 0 && !Double.isInfinite(mi)) ? x[j][s] / mi : 0.0;
            }
        }
        double[][] d = new double[M][K];
        for (int i = 0; i < M; i++) {
            for (int r = 0; r < K; r++) {
                double inflow = 0;
                for (int j = 0; j < M * K; j++) {
                    if (outflow[j] != 0) {
                        inflow += outflow[j] * P[j][i * K + r];
                    }
                }
                d[i][r] = inflow - outflow[i * K + r];
            }
        }
        return d;
    }
}
