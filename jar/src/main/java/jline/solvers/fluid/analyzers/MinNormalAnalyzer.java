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
import jline.solvers.fluid.moments.FluidLyapunov;
import jline.solvers.fluid.moments.FluidMomentTerms;
import jline.solvers.fluid.moments.FluidRefineMeanfield;
import jline.solvers.fluid.moments.FluidNonHyperbolicException;
import jline.util.matrix.Matrix;
import org.apache.commons.math3.util.FastMath;

import static jline.io.InputOutput.line_error;
import static jline.io.InputOutput.mfilename;

/**
 * Second-order moment-closure fluid analysis, backing
 * {@code options.method='minnormal'}. Java twin of the MATLAB
 * {@code solver_fluid_moments}.
 *
 * <p>The default fluid methods close the moment hierarchy at first order: the
 * drift of the mean depends on {@code E[min(X_i,c_i)]}, which they replace by
 * {@code min(E[X_i],c_i)}. No second moment ever enters, so no variance is
 * produced and the mean itself is biased wherever min() is not locally linear.
 * This analyzer reinstates the second moment with the min-normal closure of
 * Guenther, Stefanek and Bradley: the drift uses {@code E[min(X_i,c_i)]} under a
 * normal marginal whose variance is produced by the covariance (Lyapunov)
 * equation, so mean and covariance are solved self-consistently by fixed-point
 * iteration. This corrects the mean, most visibly near rho = 1 where the
 * first-order closure is worst, and it is the only fluid method that can
 * represent GPS at all.</p>
 *
 * <p>All performance measures are read back from the same event representation
 * that defines the drift, so throughputs balance flow at the fixed point under
 * whichever closure was used. That is also why this analyzer does not reuse the
 * closing metric reader, which accepts SIRO as FCFS and has no GPS branch.</p>
 *
 * @see FluidMomentTerms
 * @see FluidLyapunov
 */
public class MinNormalAnalyzer extends ClosingAndStateDepMethodsAnalyzer {

    /** State-level stationary covariance at the converged fixed point. */
    public Matrix sigmaMatrix;
    /** Per-(station,class) queue-length variance. */
    public Matrix qVar;
    /** State coordinates of each (station,class), flattened as {@code i*K+k}. */
    public int[][] classBlock;
    /** State coordinates of each station, i.e. every class and phase it serves. */
    public int[][] stationBlock;
    /** Per-station population variance. */
    public double[] sigma2;
    /**
     * The same variance as it enters the DRIFT: zero at the delay stations, which
     * have no min() to close. This, not {@link #sigma2}, is what any later solve on
     * the same fixed point (the passage-time ODE) must close its capacity term at,
     * or it evaluates a first-order drift at a second-order fixed point.
     */
    public double[] sigma2Drift;
    /** Number of outer (mean, covariance) iterations performed. */
    public int outerIters;

    /**
     * Relative offset at which the two sides of a saturation kink are probed.
     * Well outside the sqrt(eps) band that defines the kink, and small enough
     * that the Jacobian is the one-sided limit rather than a nearby point's.
     */
    private static final double KINK_PROBE = 1e-6;

    @Override
    public void analyze(NetworkStruct sn, SolverOptions options, SolverResult result) {
        long t0 = System.nanoTime();
        int M = sn.nstations;
        int K = sn.nclasses;

        FluidMomentTerms terms = new FluidMomentTerms(sn, options);
        terms.checkSupported();

        // the covariance is a dense nstate-by-nstate object and the Lyapunov solve is
        // cubic in it, so refuse rather than silently crawl
        int maxstate = (options.config == null) ? 200 : options.config.moment_maxstate;
        if (maxstate <= 0) {
            maxstate = 200;
        }
        if (terms.nstate > maxstate) {
            line_error(mfilename(new Object() {
            }), String.format("The moment-closure methods solve a %dx%d Lyapunov equation, above the limit of %d "
                    + "set by options.config.moment_maxstate. Raise that limit or use options.method=\"closing\".",
                    terms.nstate, terms.nstate, maxstate));
        }

        boolean[] shareSched = new boolean[M];
        for (int i = 0; i < M; i++) {
            shareSched[i] = terms.sched[i] == SchedStrategy.PS || terms.sched[i] == SchedStrategy.FCFS
                    || terms.sched[i] == SchedStrategy.DPS || terms.sched[i] == SchedStrategy.GPS;
        }

        int outerMax = 20;
        if (options.iter_max > 0) {
            outerMax = FastMath.min(20, FastMath.max(2, options.iter_max));
        }
        // THE CLOSURE IS JUDGED FAR TIGHTER THAN CoarseTol, so it must not stop
        // there. Converged only to 1e-3, this alternation is not a fixed point to
        // two machines: on mqn_singleserver_ps the MATLAB twin answered 42.3962 on
        // two hosts and 42.8207 on a third, 1e-2 relative apart, because the
        // transient iterate below fell on opposite sides. 1e-6 is the loosest that
        // reproduces; MIN, not assignment, so a caller may still ask tighter.
        //
        // It governs the INNER mean solve too, through meanopt. That is not
        // incidental: the transient non-hyperbolic iterate this method used to
        // abort on was an artefact of a loosely converged inner solve, so
        // tightening the outer loop alone would leave the abort in place.
        double momTol = 1e-6;
        if (options.iter_tol > 0) {
            momTol = FastMath.min(momTol, options.iter_tol);
        }
        double outerTol = momTol;

        SolverOptions meanopt = options.copy();
        meanopt.method = "closing"; // the closure enters through config.moment_sigma2
        meanopt.iter_tol = momTol;

        double[] sigma2cur = new double[M];
        // the capacity share is a RATIO of populations, so closing it needs the covariance
        // BETWEEN the station coordinates, not only the station total
        Matrix[] covblk = new Matrix[M];
        double[] sigma2solve = sigma2cur.clone();
        Matrix[] covsolve = covblk.clone();
        // The last variance whose Lyapunov solve SUCCEEDED, and the floor on the
        // step taken toward the next one. See the damping in the loop below.
        double[] sigma2ok = new double[M];
        Matrix[] covok = new Matrix[M];
        final double DAMP_MIN = 1.0 / 64.0;
        Matrix sigma = new Matrix(terms.nstate, terms.nstate);
        double[] x = null;
        int outer = 0;

        // A DELAY STATION HAS NO min() TO CLOSE, so its variance must never reach the
        // drift -- only the report. This mask used to be applied to sigma2drift after
        // the loop and nowhere inside it, so every mean solve of the fixed point ran
        // with the delay variance switched on. The rate factor there is mu*n, which
        // the Gaussian correction turns into something that does not vanish with n:
        // the coordinate is driven NEGATIVE, the drift is conservative so another
        // coordinate grows to match, and the trajectory leaves the simplex for good.
        // On CQN_Cox_CS_9 (Delay + PS + PS(c=5), N=6) the first window past
        // sigma2 = 0 moved 8.7e3 of mass and the drift norm reached 5.4e9.
        boolean[] noDriftVar = new boolean[M];
        for (int i = 0; i < M; i++) {
            noDriftVar[i] = terms.sched[i] == SchedStrategy.INF || terms.sched[i] == SchedStrategy.EXT;
        }

        for (outer = 1; outer <= outerMax; outer++) {
            // A TRANSIENT ITERATE MUST NOT VETO THE METHOD. The Lyapunov gate asks
            // whether the linear noise approximation has a stationary covariance at
            // the point THIS iterate landed on; a fixed point that fails it is a
            // model the closure cannot answer, but an intermediate iterate that
            // fails it is only a variance step that overshot. Letting one abort the
            // solve discarded answers the method reaches: on mqn_singleserver_ps
            // iterate 1 was stable at -4.93e-03, iterate 2 declined at +1.07e+01,
            // and the fixed point the fallback then found was stable at -4.92e-03.
            // So a failing iterate RETREATS toward the last variance that
            // succeeded, halving until the LNA is defined again; only a step below
            // DAMP_MIN, or a failure at the seed where there is nothing to retreat
            // toward, is the model's own non-hyperbolicity and still throws.
            double step = 1.0;
            double[] sigma2try;
            Matrix[] covtry;
            while (true) {
            sigma2try = blendSigma(sigma2ok, sigma2cur, step);
            covtry = blendCov(covok, covblk, step);
            sigma2solve = sigma2try.clone();
            covsolve = covtry.clone();
            for (int i = 0; i < M; i++) {
                if (noDriftVar[i]) {
                    sigma2solve[i] = 0;
                    covsolve[i] = null;
                }
            }
            // A station whose occupancy cannot reach its server count has min(n,c) = n
            // on the whole support, so the closure there must stay first order: see
            // FluidMomentTerms.buildMinExact, which decides it from the chain
            // populations. Its share closure follows, because mu_r*(n_r/n)*min(n,c)
            // collapses to mu_r*n_r once the min is the identity. The covariance is
            // still solved for these stations and still reported, it just does not
            // enter the drift, exactly as at the delay stations below.
            for (int i = 0; i < terms.M; i++) {
                if (terms.minExact[i]) {
                    sigma2solve[i] = 0;
                    covsolve[i] = null;
                }
            }
            meanopt.config.moment_sigma2 = sigma2solve;
            meanopt.config.moment_cov = covsolve;
            solver_fluid(sn, meanopt, result);

            x = new double[terms.nstate];
            for (int j = 0; j < terms.nstate; j++) {
                x[j] = this.xvec_it.get(0, j);
            }
            Matrix r = terms.rates(x, sigma2solve, covsolve);
            // A point ON a saturation kink has no Jacobian: the rate factor
            // min(n_i,c_i) has slope 1 below c_i and 0 above, and jacobian() would
            // hand back the left slope on the strict n_i > c_i tie-break. Which side
            // the integrator stops on is a rounding residue, so a verdict read off
            // that side depends on arithmetic noise rather than on the model.
            //
            // The verdict, not the point, is what has to be side-independent. Both
            // one-sided Jacobians exist and are ordinary matrices, so ASK BOTH: only
            // if they disagree on hyperbolicity is the answer being decided by the
            // rounding residue, and only then is declining right. Refusing at every
            // kink instead threw away models the reference solves, because the first
            // outer iterate runs at sigma2 = 0 and a saturated model's first-order
            // fixed point lands on the kink by construction -- Delay(1) -> PS(0.5),
            // N=3 is one, where MATLAB and the two other ports return the closure's
            // 1.5751 while this analyzer declined to the first-order 2.0. Later
            // iterates carry a positive sigma2 and are smooth, so the check costs two
            // Jacobians on the seed and nothing after it.
            int kink = terms.driftKinkStation(x, sigma2solve);
            if (kink >= 0) {
                for (int side = 0; side < 2; side++) {
                    double rel = (side == 0) ? -KINK_PROBE : KINK_PROBE;
                    double[] xs = terms.nudgedOffKink(x, sigma2solve, rel);
                    try {
                        localLyapunov(terms.jacobian(xs, sigma2solve, covsolve),
                                terms.rates(xs, sigma2solve, covsolve), terms);
                    } catch (FluidNonHyperbolicException e) {
                        throw new FluidNonHyperbolicException(String.format(
                                "The fluid fixed point sits on the saturation kink of station %d (population "
                                + "equals its %g servers) and the two one-sided drift Jacobians there disagree "
                                + "on hyperbolicity, so which of them the linear noise approximation would use "
                                + "is decided by the integrator's rounding residue rather than by the model. "
                                + "This is the saturated boundary of a continuum of equilibria; use "
                                + "options.method=\"closing\" for the mean only. Underlying: %s",
                                kink + 1, sn.nservers.get(kink, 0), e.getMessage()));
                    }
                }
            }
            Matrix A = terms.jacobian(x, sigma2solve, covsolve);
            try {
                sigma = localLyapunov(A, r, terms);
                break;
            } catch (FluidNonHyperbolicException e) {
                boolean atSeed = norm1Diff(sigma2cur, sigma2ok) == 0 && allNull(covblk);
                if (step <= DAMP_MIN || atSeed) {
                    throw e;
                }
                step = step / 2;
            }
            }
            sigma2ok = sigma2try;
            covok = covtry;

            double[] sigma2new = new double[M];
            Matrix[] covnew = new Matrix[M];
            for (int i = 0; i < M; i++) {
                int[] blk = terms.stationBlock[i];
                if (blk.length == 0) {
                    continue;
                }
                sigma2new[i] = FastMath.max(0, blockSum(sigma, blk));
                if (shareSched[i]) {
                    covnew[i] = subMatrix(sigma, blk);
                }
            }

            double delta = norm1Diff(sigma2new, sigma2try) / FastMath.max(1, norm1(sigma2new));
            for (int i = 0; i < M; i++) {
                // sigma2 is the SUM of a block, so it can converge while the
                // off-diagonals the share closure reads are still moving
                if (covnew[i] == null) {
                    continue;
                }
                double num = (covtry[i] == null) ? matNorm1(covnew[i]) : matNorm1Diff(covnew[i], covtry[i]);
                delta = FastMath.max(delta, num / FastMath.max(1, matNorm1(covnew[i])));
            }
            sigma2cur = sigma2new;
            covblk = covnew;
            if (delta < outerTol) {
                break;
            }
        }
        this.outerIters = FastMath.min(outer, outerMax);

        // Metrics must be read at the SAME variance the mean solve used, not at the
        // variance that solve produced: the latter evaluates the rate functions at a
        // point that is not their fixed point, and throughput then fails to balance.
        // The delay stations have no min() to close, so their variance must not enter
        // the drift; keep it for reporting only. noDriftVar already masked sigma2solve
        // on the way in, so this IS the variance the mean solve used, which is what the
        // paragraph above requires.
        double[] sigma2drift = sigma2solve.clone();
        Matrix[] covdrift = covsolve.clone();

        // THE 'refined' METHOD: the O(1/N) correction of Gast (POMACS 2017),
        // taken about the MEAN-FIELD fixed point and not the Gaussian one.
        // Adding it to the 'minnormal' fixed point would count the same O(1/N)
        // term TWICE, since the Gaussian closure already resums it: expanding
        // E[F(X)] to second order and setting it to zero reproduces exactly the
        // Gast correction equation. So the base point is recomputed with the
        // first-order closure, while the Hessian and the Jacobian are taken from
        // the smooth Gaussian drift, the hard min being only piecewise linear
        // and, at saturation, kinked exactly at the fixed point.
        // Port of solver_fluid_moments.m's own 'refined' branch.
        boolean isRefined = "refined".equals(options.method) || "fluid.refined".equals(options.method);
        if (isRefined) {
            SolverOptions mfopt = options.copy();
            mfopt.method = "closing";
            mfopt.config.moment_sigma2 = new double[M];
            mfopt.config.moment_cov = new Matrix[M];
            solver_fluid(sn, mfopt, result);
            double[] xmf = new double[terms.nstate];
            for (int j = 0; j < terms.nstate; j++) {
                xmf[j] = this.xvec_it.get(0, j);
            }
            Matrix Amf = terms.jacobian(xmf, sigma2drift, covdrift);
            Matrix sigmaMf = localLyapunov(Amf, terms.rates(xmf, sigma2drift, covdrift), terms);

            // A LINEAR DRIFT NEEDS NO REFINEMENT, and that is not the degenerate
            // call FluidRefineMeanfield refuses. When every station is either an
            // infinite server or minExact -- min(n,c) is the identity on the
            // reachable set, the population bound never reaching c -- the drift
            // is exactly affine there, its Hessian vanishes and the O(1/N)
            // correction is identically zero. The mask above then zeroes all of
            // sigma2drift, which the refinement reads as "the caller handed me
            // the first-order closure" and rejects. Settle it here, where the
            // reason for the zero is known: a null correction, not an error.
            // Delay + PS(c=2) at N=2 is the smallest case.
            boolean allLinear = true;
            for (int i = 0; i < M; i++) {
                if (!(noDriftVar[i] || terms.minExact[i])) {
                    allLinear = false;
                    break;
                }
            }
            double[] refinement = new double[terms.nstate];
            if (!allLinear) {
                refinement = FluidRefineMeanfield.refine(xmf, sigma2drift, sigmaMf, terms, covdrift).V;
            }
            x = new double[terms.nstate];
            for (int j = 0; j < terms.nstate; j++) {
                x[j] = FastMath.max(0, xmf[j] + refinement[j]);
            }
            // the corrected point is a correction OF the mean-field fixed point,
            // so its rates are read with the mean-field (zero) variance
            sigma2drift = new double[M];
            covdrift = new Matrix[M];
            sigma = sigmaMf;
            double[] sigma2ref = new double[M];
            for (int i = 0; i < M; i++) {
                int[] blk = terms.stationBlock[i];
                if (blk.length == 0) {
                    continue;
                }
                sigma2ref[i] = FastMath.max(0, blockSum(sigmaMf, blk));
            }
            sigma2cur = sigma2ref;
        }

        Matrix r = terms.rates(x, sigma2drift, covdrift);
        Matrix gfac = terms.factors(x, sigma2drift, covdrift);

        // a load-dependent station clears alpha(n) times the nominal work, so its
        // utilization normalises by the peak scaling (T*S/peak, as in the CTMC)
        Matrix Seff = terms.S.copy();
        if (terms.lldscaling != null) {
            for (int i = 0; i < FastMath.min(M, terms.lldscaling.getNumRows()); i++) {
                double peak = terms.S.get(i, 0);
                for (int j = 0; j < terms.lldscaling.getNumCols(); j++) {
                    peak = FastMath.max(peak, terms.lldscaling.get(i, j));
                }
                Seff.set(i, 0, peak);
            }
        }

        result.QN = new Matrix(M, K);
        result.UN = new Matrix(M, K);
        result.RN = new Matrix(M, K);
        result.TN = new Matrix(M, K);
        for (int i = 0; i < M; i++) {
            for (int k = 0; k < K; k++) {
                int[] blk = terms.classBlock[i * K + k];
                if (blk.length == 0) {
                    continue;
                }
                double q = 0;
                double g = 0;
                for (int j = 0; j < blk.length; j++) {
                    q += x[blk[j]];
                    g += gfac.get(blk[j], 0);
                }
                // A Source holds no jobs: the coordinate is the arrival
                // process's normalisation constant, not a population (same
                // exclusion as in the covariance solve below), and the
                // station-level routing folds departures of an open class back
                // onto it. See ClosingAndStateDepMethodsAnalyzer.
                result.QN.set(i, k, terms.sched[i] == SchedStrategy.EXT ? 0.0 : q);
                result.UN.set(i, k, terms.sched[i] == SchedStrategy.INF ? q : g / Seff.get(i, 0));
                // Summed over ORIGINAL events through Emap: a reduced event folded through an
                // immediate coordinate is a completion at more than one (station,class), and its
                // rate has to reach every one of them. Emap is the identity without a reduction.
                double tn = 0;
                for (int e = 0; e < terms.Emap.getNumRows(); e++) {
                    double w = departureWeight(terms, e, i, k);
                    if (w != 0) {
                        tn += r.get(e, 0) * w;
                    }
                }
                result.TN.set(i, k, tn);
                // See SolverFluid: TN is zero only to the integrator's accuracy.
                if (tn > GlobalConstants.Zero) {
                    result.RN.set(i, k, q / tn);
                }
            }
        }

        // transients, evaluated on the trajectory the ODE actually integrated
        int nt = this.xvec_t.getNumRows();
        result.QNt = new Matrix[M][K];
        result.UNt = new Matrix[M][K];
        result.TNt = new Matrix[M][K];
        Matrix Gt = new Matrix(nt, terms.nstate);
        Matrix Rt = new Matrix(nt, terms.eventIdx.length);
        double[] xs = new double[terms.nstate];
        for (int s = 0; s < nt; s++) {
            for (int j = 0; j < terms.nstate; j++) {
                xs[j] = this.xvec_t.get(s, j);
            }
            Matrix gs = terms.factors(xs, sigma2drift, covdrift);
            Matrix rs = terms.rates(xs, sigma2drift, covdrift);
            for (int j = 0; j < terms.nstate; j++) {
                Gt.set(s, j, gs.get(j, 0));
            }
            for (int e = 0; e < terms.eventIdx.length; e++) {
                Rt.set(s, e, rs.get(e, 0));
            }
        }
        for (int i = 0; i < M; i++) {
            for (int k = 0; k < K; k++) {
                result.QNt[i][k] = new Matrix(nt, 1);
                result.UNt[i][k] = new Matrix(nt, 1);
                result.TNt[i][k] = new Matrix(nt, 1);
                int[] blk = terms.classBlock[i * K + k];
                if (blk.length == 0) {
                    continue;
                }
                for (int s = 0; s < nt; s++) {
                    double q = 0;
                    double g = 0;
                    for (int j = 0; j < blk.length; j++) {
                        q += this.xvec_t.get(s, blk[j]);
                        g += Gt.get(s, blk[j]);
                    }
                    result.QNt[i][k].set(s, 0, terms.sched[i] == SchedStrategy.EXT ? 0.0 : q);
                    result.UNt[i][k].set(s, 0,
                            terms.sched[i] == SchedStrategy.INF ? q : g / Seff.get(i, 0));
                    double tn = 0;
                    for (int e = 0; e < terms.Emap.getNumRows(); e++) {
                        double w = departureWeight(terms, e, i, k);
                        if (w != 0) {
                            tn += Rt.get(s, e) * w;
                        }
                    }
                    result.TNt[i][k].set(s, 0, tn);
                }
            }
        }

        this.sigmaMatrix = sigma;
        this.sigma2 = sigma2cur;
        this.sigma2Drift = sigma2drift;
        this.classBlock = terms.classBlock;
        this.stationBlock = terms.stationBlock;
        this.qVar = new Matrix(M, K);
        for (int i = 0; i < M; i++) {
            for (int k = 0; k < K; k++) {
                int[] blk = terms.classBlock[i * K + k];
                if (blk.length > 0) {
                    this.qVar.set(i, k, FastMath.max(0, blockSum(sigma, blk)));
                }
            }
        }

        result.method = "minnormal";
        result.runtime = (System.nanoTime() - t0) / 1.0e9;
    }

    /**
     * Stationary covariance on the coordinates that carry a real population.
     *
     * <p>For a closed model covIdx is every coordinate and this is the plain
     * solve. For an open or mixed model it drops the EXT source pool, whose
     * coordinate is a normalisation constant rather than a job count. The result
     * is scattered back to full size with zeros on the dropped rows so that the
     * station and class blocks index unchanged downstream.</p>
     */
    private static Matrix localLyapunov(Matrix A, Matrix r, FluidMomentTerms terms) {
        int[] idx = terms.covIdx;
        int nc = idx.length;
        int nevents = terms.eventIdx.length;
        Matrix Dc = new Matrix(nc, nevents);
        for (int i = 0; i < nc; i++) {
            for (int e = 0; e < nevents; e++) {
                Dc.set(i, e, terms.D.get(idx[i], e));
            }
        }
        Matrix Qc = new Matrix(nc, nc);
        for (int i = 0; i < nc; i++) {
            for (int j = 0; j < nc; j++) {
                double acc = 0;
                for (int e = 0; e < nevents; e++) {
                    acc += Dc.get(i, e) * r.get(e, 0) * Dc.get(j, e);
                }
                Qc.set(i, j, acc);
            }
        }
        Matrix Ac = new Matrix(nc, nc);
        for (int i = 0; i < nc; i++) {
            for (int j = 0; j < nc; j++) {
                Ac.set(i, j, A.get(idx[i], idx[j]));
            }
        }
        Matrix Sc = FluidLyapunov.solve(Ac, Qc, Dc, 0).sigma;
        Matrix sigma = new Matrix(terms.nstate, terms.nstate);
        for (int i = 0; i < nc; i++) {
            for (int j = 0; j < nc; j++) {
                sigma.set(idx[i], idx[j], Sc.get(i, j));
            }
        }
        return sigma;
    }

    private static double blockSum(Matrix m, int[] blk) {
        double acc = 0;
        for (int i = 0; i < blk.length; i++) {
            for (int j = 0; j < blk.length; j++) {
                acc += m.get(blk[i], blk[j]);
            }
        }
        return acc;
    }

    private static Matrix subMatrix(Matrix m, int[] blk) {
        Matrix s = new Matrix(blk.length, blk.length);
        for (int i = 0; i < blk.length; i++) {
            for (int j = 0; j < blk.length; j++) {
                s.set(i, j, m.get(blk[i], blk[j]));
            }
        }
        return s;
    }

    private static double norm1(double[] v) {
        double acc = 0;
        for (int i = 0; i < v.length; i++) {
            acc += FastMath.abs(v[i]);
        }
        return acc;
    }

    /**
     * a + step*(b - a), entry by entry. The scalar half of the damped variance
     * step; the twin of {@link #blendCov}, so the two stay consistent.
     */
    private static double[] blendSigma(double[] a, double[] b, double step) {
        double[] out = new double[b.length];
        for (int i = 0; i < b.length; i++) {
            double ai = (a != null && i < a.length) ? a[i] : 0.0;
            out[i] = ai + step * (b[i] - ai);
        }
        return out;
    }

    /**
     * a + step*(b - a) for the covariance blocks, with a null read as the zero
     * matrix and an entry left null when both sides are null.
     */
    private static Matrix[] blendCov(Matrix[] a, Matrix[] b, double step) {
        Matrix[] out = new Matrix[b.length];
        for (int i = 0; i < b.length; i++) {
            Matrix ai = (a != null && i < a.length) ? a[i] : null;
            Matrix bi = b[i];
            if (ai == null && bi == null) {
                out[i] = null;
            } else if (ai == null) {
                out[i] = bi.scale(step);
            } else if (bi == null) {
                out[i] = ai.scale(1.0 - step);
            } else {
                out[i] = ai.add(step, bi.sub(1.0, ai));
            }
        }
        return out;
    }

    /** True when every covariance block is still unset, i.e. the seed iterate. */
    private static boolean allNull(Matrix[] a) {
        if (a == null) {
            return true;
        }
        for (int i = 0; i < a.length; i++) {
            if (a[i] != null) {
                return false;
            }
        }
        return true;
    }

    private static double norm1Diff(double[] a, double[] b) {
        double acc = 0;
        for (int i = 0; i < a.length; i++) {
            acc += FastMath.abs(a[i] - b[i]);
        }
        return acc;
    }

    /** Maximum absolute column sum, matching the MATLAB matrix 1-norm. */
    private static double matNorm1(Matrix m) {
        double best = 0;
        for (int j = 0; j < m.getNumCols(); j++) {
            double acc = 0;
            for (int i = 0; i < m.getNumRows(); i++) {
                acc += FastMath.abs(m.get(i, j));
            }
            best = FastMath.max(best, acc);
        }
        return best;
    }

    private static double matNorm1Diff(Matrix a, Matrix b) {
        double best = 0;
        for (int j = 0; j < a.getNumCols(); j++) {
            double acc = 0;
            for (int i = 0; i < a.getNumRows(); i++) {
                acc += FastMath.abs(a.get(i, j) - b.get(i, j));
            }
            best = FastMath.max(best, acc);
        }
        return best;
    }

    /**
     * Weight with which the reduced event {@code e} counts as a class-{@code k} completion at
     * station {@code i}: the total, over the ORIGINAL events it stands for, of Emap times the
     * departure indicator. It is 0 or 1 when nothing was eliminated, and can exceed the single flag
     * a composed event would otherwise carry, because a job that passes instantly through a station
     * completes at more than one of them at once.
     *
     * @param terms the moment terms carrying Emap and the original-event classification
     * @param e     reduced event index
     * @param i     station index
     * @param k     class index
     * @return the departure weight of that event for that (station, class)
     */
    private static double departureWeight(FluidMomentTerms terms, int e, int i, int k) {
        double w = 0;
        for (int o = 0; o < terms.evIsDeparture.length; o++) {
            if (terms.evIsDeparture[o] && terms.evStation[o] == i && terms.evClass[o] == k) {
                w += terms.Emap.get(e, o);
            }
        }
        return w;
    }

}
