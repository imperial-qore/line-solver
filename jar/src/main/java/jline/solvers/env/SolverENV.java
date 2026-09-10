/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Integration note: Enhanced with changes from fyp25-yiran-generalized-blending.git
 * integrated from commit 83efd675f3c291737f199674f3e982427b0c0212 onwards.
 */

package jline.solvers.env;

import jline.GlobalConstants;
import static jline.GlobalConstants.Inf;
import static jline.io.InputOutput.line_debug;
import static jline.io.InputOutput.line_warning;

import jline.io.Ret;
import jline.lang.*;
import jline.lang.constant.SchedStrategy;
import jline.lang.constant.SolverType;
import jline.lang.nodes.StatefulNode;
import jline.lang.nodes.Station;
import jline.lang.nodes.ServiceStation;
import jline.lang.nodes.Source;
import jline.lang.processes.Markovian;
import jline.lang.processes.ContinuousDistribution;
import jline.lang.processes.Exp;
import jline.solvers.*;
import jline.solvers.ctmc.CTMCResult;
import jline.solvers.ctmc.SolverCTMC;
import jline.solvers.ldes.SolverLDES;
import jline.lang.Model;
import jline.lang.layered.LayeredNetwork;
import jline.solvers.Solver;
import jline.solvers.ln.SolverLN;
import jline.solvers.fluid.SolverFluid;
import jline.solvers.fluid.analyzers.FluidCacheTran;
import jline.solvers.fluid.analyzers.CacheTranResult;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;
import jline.util.Pair;
import jline.VerboseLevel;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.integration.SimpsonIntegrator;
import org.apache.commons.math3.analysis.integration.UnivariateIntegrator;
import org.apache.commons.math3.analysis.interpolation.SplineInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;
import org.apache.commons.math3.util.FastMath;

import java.text.NumberFormat;
import java.util.*;
import java.util.List;
import java.util.Set;
import java.util.HashSet;
import java.util.Comparator;
import java.util.Arrays;
import java.util.function.Function;

import static java.lang.Math.*;
import static jline.api.mam.Map_cdf.map_cdf;
import static jline.api.mam.Map_mean.map_mean;
import static jline.api.mam.Map_normalize.map_normalize;
import static jline.api.mam.Map_pie.map_pie;
import static jline.api.mc.Ctmc_makeinfgen.ctmc_makeinfgen;
import static jline.api.mc.Ctmc_solve.ctmc_solve;
import static jline.api.mc.Ctmc_kms.ctmc_kms;
import static jline.api.mc.Ctmc_takahashi.ctmc_takahashi;
import static jline.api.mc.Ctmc_multi.ctmc_multi;
import static jline.api.mc.Dtmc_solve.dtmc_solve;
import static jline.api.mam.Map_mean.map_mean;
import static jline.util.Utils.isInf;

import jline.lang.processes.Exp;
import jline.solvers.ctmc.SolverCTMC;
import jline.solvers.mam.SolverMAM;
import jline.solvers.ctmc.ResultCTMC;
import jline.solvers.ctmc.handlers.Solver_ctmc;
import jline.solvers.ctmc.handlers.Ctmc_avg_from_pi;
import jline.solvers.mam.handlers.Solver_mam_ldqbd_statevec;
import jline.api.mc.Ctmc_timeaverage;
import jline.api.mc.Ctmc_transient;
import jline.api.mc.Ctmc_solve_reducible;
import jline.util.Pair;

/**
 * ENV - Ensemble environment solver for models immersed in a random environment.
 */
public class SolverENV extends EnsembleSolver {

    // User-supplied representation of each stage transition
    private final Environment envObj;
    // see _kb/06-solver-catalog.md (JAR-only implementation notes: stage-model/stage-solver adapters)
    private Model[] envModels;
    private Solver[] envSolvers;
    private NetworkStruct[] sn;
    private final Environment.ResetQueueLengthsFunction[][] resetFromMarginal;
    private final Environment.ResetEnvRatesFunction[][] resetEnvRates;
    private MatrixCell ServerNum;
    private MatrixCell SRates;
    private String stateDepMethod;
    private Matrix E0;
    private Matrix Eutil;
    private Matrix pi;
    private Function<Double, Double>[] sojournCdfs;
    private Function<Double, Double>[] sojournCdfsUtil;
    private Function<Double, Double>[][] transitionCdfs;
    private Matrix dtmcP;
    private Matrix holdTime;
    private long startTime;
    private int ref = 0;
    private boolean SMPMethod = false;
    private boolean compression = false;
    private Compression_result compressionResult;
    private int Ecompress;
    private List<Matrix[][]> UNtStages;
    private MatrixCell tStages;

    // ---- State-vector analyzer (options.method='statevec') working data ----
    private boolean statevecMethod;
    private Matrix E0rate;            // raw env rate matrix (zero diagonal); cf. MATLAB self.E0
    private SvStage[] svStages;       // per-stage generator + backend metadata
    private Matrix[] piEnter;         // per-stage entry distribution
    private Matrix[] piEnterPrev;     // previous-iteration entry distribution
    private Matrix[][] piExitDest;    // per-stage per-destination exit distribution
    private Matrix[] piTimeAvg;       // per-stage sojourn-end distribution

    /** Per-stage state-vector working data (CTMC or MAM/LDQBD backend). */
    private static final class SvStage {
        String backend;               // "ctmc" or "mam"
        Matrix Q;                     // per-stage generator
        double[] timespan;
        // CTMC backend
        Matrix SS;
        Matrix SSaggr;
        double[][][] arvRates;
        double[][][] depRates;
        NetworkStruct snE;
        // MAM/LDQBD backend
        Solver_mam_ldqbd_statevec.Ld ld;
        int[] levelOf;
    }

    public SolverENV(Environment renv, Solver[] solvers) {
        super(renv, "SolverENV", new SolverOptions(SolverType.ENV));
        this.envObj = renv;
        this.envModels = renv.getStageModels();
        this.envSolvers = solvers;
        int E = getNumberOfModels();
        this.sn = new NetworkStruct[E];
        this.resetFromMarginal = new Environment.ResetQueueLengthsFunction[E][E];
        this.resetEnvRates = new Environment.ResetEnvRatesFunction[E][E];
        this.result = new SolverResult();
        line_debug(options.verbose, String.format("ENV solver starting: %d stages", E));

        for (int e = 0; e < E; e++) {
            this.sn[e] = envStructOf(this.envModels[e]);
             if (!solverSupportsStage(envSolvers[e], envModels[e])) {
               throw new RuntimeException("Model is not supported by the solver.");
             }
            System.arraycopy(renv.resetQLFun[e], 0, resetFromMarginal[e], 0, E);
            System.arraycopy(renv.resetEnvRatesFun[e], 0, resetEnvRates[e], 0, E);
        }

        // Auto-detect state-dependent environment
        boolean hasStateDependentRates = false;
        for (int e = 0; e < E && !hasStateDependentRates; e++) {
            for (int h = 0; h < E; h++) {
                if (resetEnvRates[e][h] != null) {
                    // Non-null resetEnvRates indicates state-dependent transitions
                    hasStateDependentRates = true;
                    break;
                }
            }
        }

        if (hasStateDependentRates) {
            this.stateDepMethod = "statedep";
            line_debug(options.verbose, "ENV: auto-detected state-dependent environment rates");
        }

        // Validate incompatible method combinations
        if (this.SMPMethod && "statedep".equalsIgnoreCase(this.stateDepMethod)) {
            throw new IllegalArgumentException(
                "SMP method (method='smp') is incompatible with state-dependent environments.\n" +
                "SMP method computes environment probabilities once at initialization,\n" +
                "but state-dependent environments modify transition rates during iterations.\n" +
                "Please use either method='smp' OR state-dependent rates, but not both.");
        }

        // Validate distributions are Markovian unless SMP method is enabled
        for (int e = 0; e < E; e++) {
            for (int h = 0; h < E; h++) {
                ContinuousDistribution dist = envObj.env[e][h];
                if (dist != null && !(dist instanceof Markovian) && !this.SMPMethod) {
                    throw new IllegalArgumentException(
                        String.format("The distribution of the environment transition from stage %d to %d is not supported by the SolverENV solver. Use method='smp' for non-Markovian distributions.", e, h));
                }
            }
        }
    }

    public SolverENV(Environment renv, Solver[] solvers, SolverOptions options) {
        super(renv, "SolverENV", options);

        // Enable SMP method if specified in options
        if (options != null && options.method != null && options.method.equalsIgnoreCase("smp")) {
            this.SMPMethod = true;
            if (options.verbose == VerboseLevel.DEBUG) {
                System.out.println("ENV solver: SMP method enabled via options.method='smp'");
            }
        }

        this.envObj = renv;
        this.envModels = renv.getStageModels();
        this.envSolvers = solvers;
        int E = getNumberOfModels();
        this.sn = new NetworkStruct[E];
        this.resetFromMarginal = new Environment.ResetQueueLengthsFunction[E][E];
        this.resetEnvRates = new Environment.ResetEnvRatesFunction[E][E];
        this.result = new SolverResult();

        for (int e = 0; e < E; e++) {
            this.sn[e] = envStructOf(this.envModels[e]);
             if (!solverSupportsStage(envSolvers[e], envModels[e])) {
               throw new RuntimeException("Model is not supported by the solver.");
             }
            System.arraycopy(renv.resetQLFun[e], 0, resetFromMarginal[e], 0, E);
            System.arraycopy(renv.resetEnvRatesFun[e], 0, resetEnvRates[e], 0, E);
        }

        // Auto-detect state-dependent environment
        boolean hasStateDependentRates = false;
        for (int e = 0; e < E && !hasStateDependentRates; e++) {
            for (int h = 0; h < E; h++) {
                if (resetEnvRates[e][h] != null) {
                    // Non-null resetEnvRates indicates state-dependent transitions
                    hasStateDependentRates = true;
                    break;
                }
            }
        }

        if (hasStateDependentRates) {
            this.stateDepMethod = "statedep";
            if (options != null && options.verbose == VerboseLevel.DEBUG) {
                System.out.println("ENV solver: Auto-detected state-dependent environment rates");
            }
        }

        // Validate incompatible method combinations
        if (this.SMPMethod && "statedep".equalsIgnoreCase(this.stateDepMethod)) {
            throw new IllegalArgumentException(
                "SMP method (method='smp') is incompatible with state-dependent environments.\n" +
                "SMP method computes environment probabilities once at initialization,\n" +
                "but state-dependent environments modify transition rates during iterations.\n" +
                "Please use either method='smp' OR state-dependent rates, but not both.");
        }

        // Validate distributions are Markovian unless SMP method is enabled
        for (int e = 0; e < E; e++) {
            for (int h = 0; h < E; h++) {
                ContinuousDistribution dist = envObj.env[e][h];
                if (dist != null && !(dist instanceof Markovian) && !this.SMPMethod) {
                    throw new IllegalArgumentException(
                        String.format("The distribution of the environment transition from stage %d to %d is not supported by the SolverENV solver. Use method='smp' for non-Markovian distributions.", e, h));
                }
            }
        }
    }

    public static FeatureSet getFeatureSet() {
        FeatureSet featSupported = new FeatureSet();
        featSupported.setTrue(
                new String[]{
                        // Nodes
                        "ClassSwitch",
                        "Delay",
                        "DelayStation",
                        "Queue",
                        "Sink",
                        "JobSink",
                        "Source",
                        // Distributions
                        "Coxian",
                        "Cox2",
                        "Erlang",
                        "Exp",
                        "HyperExp",
                        // Sections
                        "StatelessClassSwitcher",
                        "InfiniteServer",
                        "SharedServer",
                        "Buffer",
                        "Dispatcher",
                        "Server",
                        "RandomSource",
                        "ServiceTunnel",
                        // Scheduling strategies
                        "SchedStrategy_INF",
                        "SchedStrategy_PS",
                        "SchedStrategy_FCFS",
                        "RoutingStrategy_PROB",
                        "RoutingStrategy_RAND",
                        "RoutingStrategy_RROBIN",
                        // Customer Classes
                        "ClosedClass",
                        "OpenClass"
                });
        return featSupported;
    }

    @Override
    public boolean supports(Network model) {
        FeatureSet featUsed = model.getUsedLangFeatures();
        FeatureSet featSupported = SolverENV.getFeatureSet();
        return FeatureSet.supports(featSupported, featUsed);
    }

    public static SolverOptions defaultOptions() {
        return new SolverOptions(SolverType.ENV);
    }

    public void setStateDepMethod(String method) {
        if (method == null || method.isEmpty()) {
            throw new IllegalArgumentException("State-dependent method cannot be null or empty.");
        }
        this.stateDepMethod = method;
    }

    public void setSMPMethod(boolean SMPMethod) {
        this.SMPMethod = SMPMethod;
    }

    public void setCompression(boolean compression) {
        this.compression = compression;
    }

    /**
     * Evaluates the per-stage sojourn-time CDF Pr{S_e &lt;= t} for stage {@code e}.
     * Only available when the semi-Markov method (method='smp') has been used,
     * which is when the per-stage sojourn CDFs are computed.
     *
     * @param e the stage index (0-based)
     * @param t the time at which to evaluate the CDF
     * @return Pr{S_e &lt;= t}
     */
    public double computeSojournCdf(int e, double t) {
        if (sojournCdfs == null || sojournCdfs[e] == null) {
            throw new IllegalStateException(
                "Sojourn-time CDFs are only available with the semi-Markov method (options.method='smp').");
        }
        return sojournCdfs[e].apply(t);
    }

    // Convergence test at iteration it
    @Override
    protected boolean converged(int it) {

        int M = sn[0].nstations;
        int K = sn[0].nclasses;
        int E;
        if (compression) {
            E = compressionResult.pMacro.getNumCols();
        }
        else {
            E = getNumberOfModels();
        }

        boolean converged = true;
        if (it <= 1) {
            return false;
        }
        for (int k = 0; k < K; k++) {
            Matrix QEntry = new Matrix(M, E);
            Matrix QExit = new Matrix(M, E);
            for (int e = 0; e < E; e++) {
                SolverResult resPrev = results.get(it - 1).get(e);
                SolverResult resCurr = results.get(it).get(e);
                if (resPrev != null && resPrev.QN != null) {
                    QEntry.setColumn(e, resPrev.QN.getColumn(k));
                }
                if (resCurr != null && resCurr.QN != null) {
                    QExit.setColumn(e, resCurr.QN.getColumn(k));
                }
            }
            double tol = options.iter_tol;
            // Match MATLAB maxpe(): max(abs(1 - approx./exact)) for exact > 0
            double maxDiff = 0;
            for (int idx = 0; idx < QEntry.getNumElements(); idx++) {
                double exact = QEntry.get(idx);
                double approx = QExit.get(idx);
                if (exact > 0) {
                    double relDiff = Math.abs(1.0 - approx / exact);
                    maxDiff = Math.max(maxDiff, relDiff);
                }
            }
            if (Double.isNaN(maxDiff) || Double.isInfinite(maxDiff) || maxDiff >= tol) {
                converged = false;
            }
        }
        return converged;

//        Matrix mapes = new Matrix(1, E);
//        for (int e = 0; e < E; e++) {
//            for (int i = 0; i < M; i++) {
//                for (int j = 0; j < K; j++) {
//                    // Error is calculated only on entry value (t = 0)
//                    mapes.set(
//                            0,
//                            e,
//                            FastMath.max(
//                                    mapes.get(0, e),
//                                    Utils.mape(
//                                            Matrix.extractRows(results.get(it).get(e).QNt[i][j], 0, 1, null),
//                                            Matrix.extractRows(results.get(it - 1).get(e).QNt[i][j], 0, 1, null))));
//                }
//            }
//        }
//        return mapes.elementMax() < options.iter_tol;

        // TODO: refactor converge test from blending method
//        for (int e = 0; e < E; e++) {
//            for (int i = 0; i < M; i++) {
//                for (int j = 0; j < K; j++) {
//                    // Error is calculated only on entry value (t = 0)
//                    double diff = (results.get(it).get(e).QN.get(i,j) - results.get(it - 1).get(e).QN.get(i,j)) /
//                            results.get(it - 1).get(e).QN.get(i,j);;
//                    if (FastMath.abs(diff) > options.iter_tol) {
//                        return false;
//                    }
//                }
//            }
//        }
//        return true;

    }

    @Override
    public int getNumberOfModels() {
        return envModels != null ? envModels.length : 0;
    }

    // see _kb/06-solver-catalog.md (JAR-only implementation notes: stage-model/stage-solver adapters)

    private NetworkStruct envStructOf(Model m) {
        if (m instanceof LayeredNetwork) return ((LayeredNetwork) m).getEnvStruct();
        return ((Network) m).getStruct(true);
    }

    private void initFromMarginalOf(Model m, Matrix qn) {
        if (m instanceof LayeredNetwork) ((LayeredNetwork) m).initFromMarginal(qn);
        else ((Network) m).initFromMarginal(qn);
    }

    private java.util.List<Station> stationsOf(Model m) {
        if (m instanceof LayeredNetwork) return ((LayeredNetwork) m).getStations();
        return ((Network) m).getStations();
    }

    private int nodesOf(Model m) {
        if (m instanceof LayeredNetwork) return ((LayeredNetwork) m).getNumberOfNodes();
        return ((Network) m).getNumberOfNodes();
    }

    private void solverTranAvgOf(Solver s) {
        if (s instanceof SolverLN) ((SolverLN) s).getTranAvg();
        else if (s instanceof NetworkSolver) ((NetworkSolver) s).getTranAvg();
    }

    private void solverAvgOf(Solver s) {
        if (s instanceof SolverLN) ((SolverLN) s).getEnsembleAvg();
        else if (s instanceof NetworkSolver) ((NetworkSolver) s).getAvg();
    }

    private boolean solverSupportsStage(Solver s, Model m) {
        // Flat stages use the existing feature check; LQN stages defer to the
        // per-layer support check inside SolverLN.getEnsembleAvg.
        if (s instanceof NetworkSolver && m instanceof Network) {
            return ((NetworkSolver) s).supports((Network) m);
        }
        return true;
    }

    @Override
    protected void init() {
        envObj.init();

        int M = sn[0].nstations;
        int K = sn[0].nclasses;
        int E = getNumberOfModels();
        line_debug(options.verbose, String.format("ENV init: %d stages, %d stations, %d classes, SMPMethod=%b, stateDepMethod=%s",
            E, M, K, SMPMethod, stateDepMethod != null ? stateDepMethod : "none"));

        // see _kb/06-solver-catalog.md (JAR-only implementation notes: ODE MaxStep tuning for ENV convergence)
        for (int e = 0; e < E; e++) {
            if (Double.isInfinite(this.envSolvers[e].options.odesolvers.odemaxstep)) {
                double tEnd = this.envSolvers[e].options.timespan[1];
                if (!Double.isInfinite(tEnd) && tEnd > 0) {
                    this.envSolvers[e].options.setODEMaxStep(tEnd / 100.0);
                }
            }
        }

        // Initialize UNtStages and tStages before any analyze() calls
        UNtStages = new ArrayList<>();
        tStages = new MatrixCell(E);

        ServerNum = new MatrixCell(K);
        SRates = new MatrixCell(K);
        for (int k = 0; k < K; k++) {
            Matrix serverNumPerClass = new Matrix(M, E);
            Matrix ratesPerClass = new Matrix(M, E);
            for (int e = 0; e < E; e++) {
                for (int m = 0; m < M; m++) {
                    serverNumPerClass.set(m, e, sn[e].nservers.get(m, 0));
                    ratesPerClass.set(m, e, sn[e].rates.get(m, k));
                }
            }
            ServerNum.set(k, serverNumPerClass);
            SRates.set(k, ratesPerClass);
        }

        final ContinuousDistribution[][] E1 = envObj.env;
        final Matrix ELocal = new Matrix(E, E);
        for (int i = 0; i < E; i++) {
            for (int j = 0; j < E; j++) {
                if (E1[i][j] != null && E1[i][j] instanceof Markovian) {
                    ELocal.set(i, j, ((Markovian) E1[i][j]).getRate());
                } else {
                    ELocal.set(i, j, 0.0);
                }
            }
        }
        Eutil = ctmc_makeinfgen(ELocal);
        E0 = new Matrix(Eutil);
        // see _kb/06-solver-catalog.md (JAR-only implementation notes: E0 vs E0rate distinction)
        this.E0rate = new Matrix(ELocal);
        this.pi = envObj.probEnv;

        // see _kb/06-solver-catalog.md (JAR-only implementation notes: embweight CTMC-generator gate)
        Matrix embweight = new Matrix(E, E);
        double[] piArray = pi.toArray1D();
        for (int e = 0; e < E; e++) {
            double sum = 0.0;
            for (int h = 0; h < E; h++) {
                if (h != e) {
                    sum += piArray[h] * E0.get(h, e);
                }
            }
            for (int k = 0; k < E; k++) {
                if (k == e) {
                    embweight.set(k, e, 0);
                } else {
                    embweight.set(k, e, piArray[k] * E0.get(k, e) / sum);
                }
            }
        }
        // Only overwrite probOrig for state-dependent or SMP methods
        // (matching MATLAB where this overwrite only occurs inside SMPMethod block)
        if ("statedep".equalsIgnoreCase(stateDepMethod) || SMPMethod) {
            this.envObj.probOrig = embweight;
        }


        @SuppressWarnings("unchecked")
        Function<Double,Double>[][] transitionCdfs = (Function<Double,Double>[][]) new Function[E][E];
        this.transitionCdfs = transitionCdfs;
        // calculate transition CDFs
        for (int k = 0; k < E; k++) {
            final int kk = k;
            for (int h = 0; h < E; h++) {
                final int hh = h;
                this.transitionCdfs[k][h] = (Double t) -> {
                    ContinuousDistribution dist = envObj.env[kk][hh];
                    if (dist == null || !(dist instanceof Markovian)) {
                        return 0.0;
                    }
                    Markovian m = (Markovian) dist;
                    double Fkh = m.evalCDF(t);
                    this.transitionCdfs[kk][hh] = m::evalCDF;
                    return Fkh;
                };
            }
        }

        if (!compression) {
            @SuppressWarnings("unchecked")
            Function<Double,Double>[] sojournCdfs = (Function<Double,Double>[]) new Function[E];
            for (int k = 0; k < E; k++) {
                final int kk = k;
                final int finalE = E;
                sojournCdfs[k] = (Double t) -> {
                    double surv = 1.0;
                    for (int i = 0; i < finalE; i++) {
                        if (i == kk) continue;
                        surv *= (1.0 - this.transitionCdfs[kk][i].apply(t));
                    }
                    return 1.0 - surv;
                };
            }
            this.sojournCdfs = sojournCdfs;
            @SuppressWarnings("unchecked")
            Function<Double,Double>[] sojournCdfsUtil = (Function<Double,Double>[]) new Function[E];
            this.sojournCdfsUtil = sojournCdfsUtil;
            System.arraycopy(sojournCdfs, 0, sojournCdfsUtil, 0, E);

        // SMPMethod: Use DTMC-based computation for Semi-Markov Processes
        // Verified numerical integration for Semi-Markov Process DTMC transition probabilities
        if (SMPMethod) {
            dtmcP = new Matrix(E, E);
            for (int k = 0; k < E; k++) {
                for (int e = 0; e < E; e++) {
                    if (k == e || envObj.env[k][e] == null) {
                        dtmcP.set(k, e, 0.0);
                    } else {
                        // compute the upper limit of the sojourn time
                        double epsilon = 1e-8;
                        double T = 1;
                        while (transitionCdfs[k][e].apply(T) < 1.0 - epsilon) {
                            T *= 2;
                        }
                        // Adaptive number of integration intervals based on T
                        int N = Math.max(1000, (int)(T * 100));
                        double dt = T / N;
                        double sum = 0;
                        for (int i = 0; i < N; i++) {
                            double t0 = i * dt;
                            double t1 = t0 + dt;
                            double deltaF = transitionCdfs[k][e].apply(t1) - transitionCdfs[k][e].apply(t0);
                            double survival = 1;
                            for (int h = 0; h < E; h++) {
                                if (h != k && h != e && envObj.env[k][h] != null) {
                                    // Use midpoint for better accuracy in survival probability calculation
                                    double tmid = (t0 + t1) / 2.0;
                                    survival *= (1.0 - envObj.env[k][h].evalCDF(tmid));
                                }
                            }
                            sum += deltaF * survival;
                        }
                        dtmcP.set(k, e, sum);
                    }
                }
            }
            Matrix dtmcPie = dtmc_solve(dtmcP);
            // Calculate the survival function for the sojourn time in the environment h_k
            Matrix holdTime = getHoldTime(E, sojournCdfs);
            this.holdTime = holdTime;

            Matrix pi = new Matrix(1, E);
            for (int k = 0; k < E; k++) {
                double sum = 0;
                for (int e = 0; e < E; e++) {
                    sum += dtmcPie.get(e) * holdTime.get(0, e);
                }
                pi.set(0, k, dtmcPie.get(k) * holdTime.get(0, k) / sum);
            }
            this.pi = pi;
            this.envObj.probEnv = pi;

            // update embweight and store it in envObj.probOrig
            Matrix newEmbweight = new Matrix(E, E);
            double[] newPiArray = pi.toArray1D();
            for (int e = 0; e < E; e++) {
                double sum = 0.0;
                for (int h = 0; h < E; h++) {
                    if (h != e) {
                        sum += newPiArray[h] * E0.get(h, e);
                    }
                }
                for (int k = 0; k < E; k++) {
                    if (k == e) {
                        newEmbweight.set(k, e, 0);
                    } else {
                        newEmbweight.set(k, e, newPiArray[k] * E0.get(k, e) / sum);
                    }
                }
            }
            this.envObj.probOrig = newEmbweight;
        }
        } else {
            // Compression block
            MatrixCell MS;
            if (E <= 10) {
                MS = findBestPartition(E0);
                if (MS == null) {
                    MS = new MatrixCell(E);
                    for (int i = 0; i < E; i++) {
                        MS.set(i, new Matrix(new int[]{i}));
                    }
                    Ecompress = E;
                }
                E = Ecompress;
            } else {
                int B = 3;
                // Alpha parameter for Courtois decomposition - controls the coupling threshold
                // Smaller values allow weaker coupling between groups (default: 0.01)
                double alpha = (options.config != null && options.config.env_alpha != null) ?
                               options.config.env_alpha : 0.01;
                List<MatrixCell> beam = new ArrayList<>();
                {
                    MatrixCell singletons = new MatrixCell(E);
                    for (int i = 0; i < E; i++) {
                        singletons.set(i, new Matrix(new int[]{i}));
                    }
                    beam.add(singletons);
                }
                MatrixCell bestSeen = beam.get(0);
                double bestEps = ctmc_decompose(E0, bestSeen, options).eps;

                for (int depth = 1; depth < E; depth++) {
                    List<Pair<MatrixCell, Double>> candidates = new ArrayList<>();
                    for (MatrixCell ms : beam) {
                        List<Set<Integer>> blocks = new ArrayList<>();
                        for (int b = 0; b < ms.size(); b++) {
                            double[][] rows = ms.get(b).toArray2D();
                            Set<Integer> s = new HashSet<>();
                            for (double[] row : rows) {
                                s.add((int) row[0]);
                            }
                            blocks.add(s);
                        }
                        for (int i = 0; i < blocks.size(); i++) {
                            for (int j = i + 1; j < blocks.size(); j++) {
                                List<Set<Integer>> trial = new ArrayList<>();
                                for (int k = 0; k < blocks.size(); k++) {
                                    if (k == i) {
                                        Set<Integer> merged = new HashSet<>(blocks.get(i));
                                        merged.addAll(blocks.get(j));
                                        trial.add(merged);
                                    } else if (k != j) {
                                        trial.add(new HashSet<>(blocks.get(k)));
                                    }
                                }
                                MatrixCell child = new MatrixCell(trial.size());
                                for (int k = 0; k < trial.size(); k++) {
                                    int[] idx = trial.get(k).stream().mapToInt(x -> x).toArray();
                                    child.set(k, new Matrix(idx));
                                }
                                Compression_result cr = ctmc_decompose(E0, child, options);
                                double epsChild = cr.eps;
                                double epsMaxChild = cr.epsMax;

                                double cost = epsChild - epsMaxChild + alpha * depth;
                                if (epsChild == 0) {
                                    continue;
                                }
                                candidates.add(new Pair<>(child, cost));
                                if (cost < bestEps) {
                                    bestEps  = cost;
                                    bestSeen = child;
                                }
                            }
                        }
                    }
                    candidates.sort(Comparator.comparingDouble(Pair::getRight));
                    beam.clear();
                    for (int i = 0; i < Math.min(B, candidates.size()); i++) {
                        beam.add(candidates.get(i).getLeft());
                    }
                }
                MS = bestSeen;
                Ecompress = MS.size();
            }

            Compression_result compressionResult = ctmc_decompose(E0, MS, options);

            if (compressionResult.eps > compressionResult.epsMax) {
                System.out.println("This model cannot be compressed, its eps is larger than epsMax");
            }

            this.pi = new Matrix(compressionResult.pMacro);
            envObj.probEnv = pi;
            this.compressionResult = compressionResult;
            Matrix Pmacro = compressionResult.G;
            double q      = compressionResult.q;
            Matrix Qmacro = new Matrix(Pmacro);
            Qmacro.subEq(Matrix.eye(Ecompress));
            Qmacro.scaleEq(q);
            Eutil = ctmc_makeinfgen(Qmacro);
            E0 = new Matrix(Eutil);
            System.out.println("eps: "+ compressionResult.eps + ", epsMax: " + compressionResult.epsMax);

            // embedding weight matrix
            Matrix embweight2 = new Matrix(Ecompress, Ecompress);
            double[] piArrayCompressed = pi.toArray1D();
            for (int e = 0; e < Ecompress; e++) {
                double sum = 0.0;
                for (int h = 0; h < Ecompress; h++) {
                    if (h != e) {
                        sum += piArrayCompressed[h] * E0.get(h, e);
                    }
                }
                for (int k = 0; k < Ecompress; k++) {
                    if (k == e) {
                        embweight2.set(k, e, 0);
                    } else {
                        embweight2.set(k, e, piArrayCompressed[k] * E0.get(k, e) / sum);
                    }
                }
            }
            this.envObj.probOrig = embweight2;

            //sojourn CDFs: weighted mixture over micro-states in each macro-block
            @SuppressWarnings("unchecked")
            Function<Double, Double>[] sojournCdfs = (Function<Double, Double>[]) new Function[Ecompress];
            int procSize = 0;
            for (int i = 0; i < Ecompress; i++) {
                final int ii = i;
                int subSize = MS.get(i).getNumRows();
                final int[] micros = MS.get(i).toIntArray1D();
                final double[] weights = new double[subSize];
                for (int r = 0; r < subSize; r++) {
                    weights[r] = compressionResult.pmicro.get(procSize + r, 0);
                }
                final int finalEcompress = Ecompress;
                final MatrixCell finalMS = MS;
                sojournCdfs[ii] = (Double t) -> {
                    if (t < 0) return 0.0;
                    double cdfSum = 0.0;
                    for (int r = 0; r < micros.length; r++) {
                        int mi = micros[r];
                        double surv = 1.0;
                        for (int j = 0; j < finalEcompress; j++) {
                            if (j == ii) continue;
                            for (int idx = 0; idx < finalMS.get(j).getNumRows(); idx++) {
                                int ds = (int) finalMS.get(j).get(idx, 0);
                                Function<Double, Double> tcdf = transitionCdfs[mi][ds];
                                if (tcdf != null) {
                                    surv *= 1.0 - tcdf.apply(t);
                                }
                            }
                        }
                        double Fr = 1.0 - surv;
                        cdfSum += weights[r] * Fr;
                    }
                    return cdfSum;
                };
                procSize += subSize;
            }
            this.sojournCdfs     = sojournCdfs;
            this.sojournCdfsUtil = Arrays.copyOf(sojournCdfs, sojournCdfs.length);

            // Update transition CDFs computation for compressed environment
            @SuppressWarnings("unchecked")
            Function<Double, Double>[][] macroTransitionCdfs = new Function[Ecompress][Ecompress];
            procSize = 0;
            for (int i = 0; i < Ecompress; i++) {
                final int subSize = MS.get(i).getNumRows();
                final int[] micros = new int[subSize];
                for (int r = 0; r < subSize; r++) {
                    micros[r] = (int) MS.get(i).get(r, 0);
                }

                final double[] weights = new double[subSize];
                for (int r = 0; r < subSize; r++) {
                    weights[r] = compressionResult.pmicro.get(procSize + r, 0);
                }

                for (int j = 0; j < Ecompress; j++) {
                    if (i == j) continue;

                    final int macroI = i;
                    final int macroJ = j;

                    MatrixCell finalMS1 = MS;
                    macroTransitionCdfs[macroI][macroJ] = (Double t) -> {
                        if (t < 0) return 0.0;
                        double total = 0.0;

                        for (int r = 0; r < micros.length; r++) {
                            int m = micros[r];
                            double wm = weights[r];

                            double survival = 1.0;
                            for (int h = 0; h < Ecompress; h++) {
                                if (h == macroJ) continue;
                                for (double dstMicro : Matrix.columnMatrixToDoubleArray(finalMS1.get(h))) {
                                    int f = (int) dstMicro;
                                    Function<Double, Double> F_mf = envObj.env[m][f] == null ? null : envObj.env[m][f]::evalCDF;
                                    if (F_mf != null) {
                                        survival *= (1.0 - F_mf.apply(t));
                                    }
                                }
                            }

                            double toJ = 0.0;
                            for (double dstMicro : Matrix.columnMatrixToDoubleArray(finalMS1.get(macroJ))) {
                                int e = (int) dstMicro;
                                Function<Double, Double> F_me = envObj.env[m][e] == null ? null : envObj.env[m][e]::evalCDF;
                                if (F_me != null) {
                                    toJ += F_me.apply(t);
                                }
                            }

                            total += wm * survival * toJ;
                        }

                        return total;
                    };
                }

                procSize += subSize;
            }
            this.transitionCdfs = macroTransitionCdfs;

            // Build macro‐state networks for the compressed environment
            Network[] macroEnsemble = new Network[Ecompress];
            NetworkSolver[] macroSolvers = new NetworkSolver[Ecompress];
            NetworkStruct[] macroSn = new NetworkStruct[Ecompress];
            procSize = 0;
            for (int i = 0; i < Ecompress; i++) {
                Network CompressedNetwork = ((Network) envModels[i]).copy();
                for (int m = 0; m < sn[0].nstations; m++) {
                    for (int k = 0; k < sn[0].nclasses; k++) {
                        double rateSum = 0.0;
                        // Sum rates weighted by micro‐state probabilities
                        for (int r = 0; r < MS.get(i).getNumRows(); r++) {
                            int microIdx = (int) MS.get(i).get(r, 0);
                            double w = compressionResult.pmicro.get(procSize + r, 0);
                            rateSum += w * sn[microIdx].rates.get(m, k);
                        }
                        JobClass jobclass = CompressedNetwork.getClasses().get(k);
                        Station st = CompressedNetwork.getStations().get(m);
                        if (st instanceof jline.lang.nodes.Queue) {
                            ((jline.lang.nodes.Queue) st).setService(jobclass, new Exp(rateSum));
                        } else if (st instanceof ServiceStation) {
                            ((ServiceStation) st).setService(jobclass, new Exp(rateSum));
                        }
                        CompressedNetwork.refreshRates(null, null);
                    }
                }
                macroEnsemble[i] = CompressedNetwork;
                macroSolvers[i]  = new SolverFluid(CompressedNetwork, envSolvers[i].options);
                macroSn[i]       = CompressedNetwork.getStruct(true);
                procSize += MS.get(i).getNumRows();
            }
            this.ensemble = macroEnsemble;
            this.solvers  = macroSolvers;
            this.sn       = macroSn;
            dtmcP = new Matrix(Ecompress, Ecompress);

            ServerNum = new MatrixCell(K);
            SRates = new MatrixCell(K);
            E = getNumberOfModels();
            for (int k = 0; k < K; k++) {
                Matrix serverNumPerClass = new Matrix(M, E);
                Matrix ratesPerClass = new Matrix(M, E);
                for (int e = 0; e < E; e++) {
                    for (int m = 0; m < M; m++) {
                        serverNumPerClass.set(m, e, sn[e].nservers.get(m, 0));
                        ratesPerClass.set(m, e, sn[e].rates.get(m, k));
                    }
                }
                ServerNum.set(k, serverNumPerClass);
                SRates.set(k, ratesPerClass);
            }
        }
    }

    private static Matrix getHoldTime(int E, Function<Double, Double>[] sojournCdfs) {
        Matrix holdTime = new Matrix(1, E);
        for (int k = 0; k < E; k++) {
            int kk = k;
            UnivariateFunction surv = t -> 1 - sojournCdfs[kk].apply(t);
            UnivariateIntegrator integrator = new SimpsonIntegrator(1e-8, 1e-8, 3, 64);
            // compute the upper limit of the sojourn time
            double upperLimit = 10;
            while (surv.value(upperLimit) > 1e-8) {
                upperLimit *= 2;
            }
            double h_k = integrator.integrate(10000, surv, 0, upperLimit);
            holdTime.set(0, k, h_k);
        }
        return holdTime;
    }

    // Numerical integration for sojourn time using discrete time points
    private void getHoldTime(int E, MatrixCell tVectors) {
        Matrix holdTime = new Matrix(1, E);
        double epsilon = 1e-6;

        for (int k = 0; k < E; k++) {
            int kk = k;
            double[] tVec = tVectors.get(kk).toArray1D();
            Function<Double, Double> survFunc = t -> 1.0 - sojournCdfsUtil[kk].apply(t);

            double integral = 0.0;
            for (int i = 0; i < tVec.length - 1; i++) {
                double t1 = tVec[i];
                double t2 = tVec[i + 1];
                double mid = 0.5 * (t1 + t2);

                double survMid = survFunc.apply(mid);
                if (survMid < epsilon) break;

                double dt = t2 - t1;
                integral += survMid * dt;
            }

            holdTime.set(0, kk, integral);
        }
        this.holdTime = holdTime;
    }


    @Override
    protected void pre(int it) {

        int E = getNumberOfModels();

        if (it == 1) {
            for (int e = 0; e < E; e++) {
                try {
                    if (isInf(this.envSolvers[e].options.timespan[1])) {
                        solverAvgOf(this.envSolvers[e]);
                    } else {
                        solverTranAvgOf(this.envSolvers[e]);
                    }
                    // Match MATLAB: extract last time point from transient for initial state
                    // MATLAB pre(): QN(i,k) = QNt{i,k}.metric(end)
                    SolverResult preResult = this.envSolvers[e].result;
                    int M = sn[0].nstations;
                    int K = sn[0].nclasses;
                    Matrix QN = new Matrix(M, K);
                    if (preResult.QNt != null) {
                        for (int i = 0; i < M; i++) {
                            for (int k = 0; k < K; k++) {
                                if (preResult.QNt[i] != null && preResult.QNt[i][k] != null) {
                                    int lastRow = preResult.QNt[i][k].getNumRows() - 1;
                                    QN.set(i, k, preResult.QNt[i][k].get(lastRow, 0));
                                }
                            }
                        }
                    } else if (preResult.QN != null) {
                        QN = preResult.QN;
                    }
                    if (!(envSolvers[e] instanceof SolverFluid) && !(envModels[e] instanceof LayeredNetwork)) {
                        roundMarginalForDiscreteSolver(QN, sn[0]);
                    }
                    initFromMarginalOf(this.envModels[e], QN);
                } catch (Exception ex) {
                    // Skip pre-initialization for stages where getAvg fails
                }
            }
        }
    }

    @Override
    protected void post(int it) {
        line_debug(options.verbose, String.format("ENV post: iteration %d", it));

        int M = sn[0].nstations;
        int K = sn[0].nclasses;
        int E = getNumberOfModels();
        Matrix[][] QExit = new Matrix[E][E];
        Matrix[][] UExit = new Matrix[E][E];
        Matrix[][] TExit = new Matrix[E][E];
        Matrix[][] w = new Matrix[E][E];
        Matrix[] QEntry = new Matrix[E]; // Average entry queue-length

        // see _kb/06-solver-catalog.md (JAR-only implementation notes: EXT/Source skip in QExit)
        boolean[] isExtStation = new boolean[M];
        {
            java.util.List<jline.lang.nodes.Station> stationList = stationsOf(envModels[0]);
            for (int i = 0; i < M; i++) {
                SchedStrategy s = sn[0].sched.get(stationList.get(i));
                isExtStation[i] = (s == SchedStrategy.EXT);
            }
        }

        for (int e = 0; e < E; e++) {
            // Skip stages where analysis failed (empty result)
            SolverResult resultE = results.get(it).get(e);
            if (resultE == null || resultE.t == null || resultE.QNt == null) {
                for (int h = 0; h < E; h++) {
                    QExit[e][h] = new Matrix(M, K);
                    UExit[e][h] = new Matrix(M, K);
                    TExit[e][h] = new Matrix(M, K);
                }
                continue;
            }
            for (int h = 0; h < E; h++) {
                QExit[e][h] = new Matrix(M, K);
                UExit[e][h] = new Matrix(M, K);
                TExit[e][h] = new Matrix(M, K);
                // THE WEIGHT IS THE STAGE SOJOURN, NOT THE e -> h CLOCK: competing
                // exponentials leave the exit TIME independent of which destination
                // won, so the exit average does not depend on h -- h enters only
                // through the reset applied on the way in. Weighting by
                // `proc[e][h]` read the transient over the mean of ONE risk rather
                // than of their minimum: on renv_twostages_repairmen, whose Stage2
                // competes a 0.5 self arc with a 0.5 arc back, that is a mean of 2
                // against the sojourn's 1, and it reported Queue1 QLen 0.55879
                // against MATLAB's 0.55550.
                // COPIED because map_normalize rescales IN PLACE, and holdTime is
                // the environment's own, read again by every later iteration.
                Matrix D0 = envObj.holdTime[e].get(0).copy();
                Matrix D1 = envObj.holdTime[e].get(1).copy();
                map_normalize(D0, D1);
                // The weight lives on the sojourn scale, so the grid summed over
                // has to as well; see refineForCdf. The stage carries ONE grid
                // for every (i,r), so it is refined once here.
                CdfGrid grid = refineForCdf(results.get(it).get(e).t, D0, D1);
                for (int i = 0; i < M; i++) {
                    if (isExtStation[i]) continue; // Skip Source stations (disabled in MATLAB)
                    for (int r = 0; r < K; r++) {
                        w[e][h] = new Matrix(1, 1);
                        Matrix cdf1 =
                                map_cdf(
                                        D0, D1,
                                        Matrix.extractRows(
                                                grid.t, 1, grid.t.getNumRows(), null));

                        Matrix cdf2 =
                                map_cdf(
                                        D0, D1,
                                        Matrix.extractRows(
                                                grid.t,
                                                0,
                                                grid.t.getNumRows() - 1,
                                                null));

                        Matrix cdfDiff = cdf1.sub(1, cdf2); // probability of leaving stage e to h, which is also \pi
                        cdfDiff = cdfDiff.transpose();
                        w[e][h] = Matrix.concatRows(w[e][h], cdfDiff, null);
                        // TODO: refactor from blending method
//                        int tR = results.get(it).get(e).t.getNumRows();
//
//                        // weight is CDF of the exponential distribution
//                        double lambda = tmE.get(e, e);
//                        Matrix weight = new Matrix(tR, 1);
//                        double[] expVals = new double[tR];
//                        for (int j = 0; j < tR; j++) {
//                            expVals[j] = 1 - exp(lambda * results.get(it).get(e).t.get(j, 0));
//                        }
//
//                        for (int j = 0; j < tR-1; j++) {
//                            weight.set(j, 0, expVals[j+1] - expVals[j]);
//                        }
//                        weight.set(tR - 1,0, 0.0);

                        if (!w[e][h].hasNaN()) {
                            QExit[e][h].set(
                                    i,
                                    r,
                                    grid.on(results.get(it).get(e).QNt[i][r]).transpose().mult(w[e][h], null).value()
                                            / w[e][h].elementSum());
                            UExit[e][h].set(
                                    i,
                                    r,
                                    grid.on(results.get(it).get(e).UNt[i][r]).transpose().mult(w[e][h], null).value()
                                            / w[e][h].elementSum());
                            TExit[e][h].set(
                                    i,
                                    r,
                                    grid.on(results.get(it).get(e).TNt[i][r]).transpose().mult(w[e][h], null).value()
                                            / w[e][h].elementSum());
                        }
                    }
                }
            }
        }

        for (int e = 0; e < E; e++) {
            // Skip stages where analysis failed (no valid results)
            SolverResult resultEForEntry = results.get(it).get(e);
            if (resultEForEntry == null || resultEForEntry.t == null || resultEForEntry.QNt == null) {
                continue;
            }
            QEntry[e] = new Matrix(M, K);
            for (int h = 0; h < E; h++) {
                // Probability of coming from h to e \times resetFun(Qexit from h to e
                if (envObj.probOrig.get(h, e) > 0) {
                    Matrix partialQEntry = new Matrix(0, 0);
                    resetFromMarginal[h][e]
                            .reset(QExit[h][e])
                            .scaleEq(envObj.probOrig.get(h, e), partialQEntry);
                    QEntry[e] = QEntry[e].add(1, partialQEntry);
                }
            }

            // see _kb/06-solver-catalog.md (JAR-only implementation notes: QEntry normalization)
            NetworkStruct snRef = sn[0];
            for (int c = 0; c < snRef.nchains; c++) {
                double njobs_chain = 0;
                double state_chain = 0;
                java.util.List<Integer> chainClasses = new java.util.ArrayList<>();
                for (int k = 0; k < K; k++) {
                    if (snRef.chains.get(c, k) > 0) {
                        chainClasses.add(k);
                        njobs_chain += snRef.njobs.get(0, k);
                    }
                }
                if (Double.isInfinite(njobs_chain)) continue; // open chain
                for (int i = 0; i < M; i++) {
                    for (int k : chainClasses) {
                        state_chain += QEntry[e].get(i, k);
                    }
                }
                if (state_chain > 0 && Math.abs(state_chain - njobs_chain) > 1e-10) {
                    double scale = njobs_chain / state_chain;
                    for (int i = 0; i < M; i++) {
                        for (int k : chainClasses) {
                            QEntry[e].set(i, k, QEntry[e].get(i, k) * scale);
                        }
                    }
                }
            }

            if (!(envSolvers[e] instanceof SolverFluid) && !(envModels[e] instanceof LayeredNetwork)) {
                roundMarginalForDiscreteSolver(QEntry[e], sn[0]);
            }
            envSolvers[e].reset();
            initFromMarginalOf(envModels[e], QEntry[e]);
        }

        // Update transition rates between stages if State Dependent
        // Auto-detected or manually specified via options.method='statedep'
        if ("statedep".equalsIgnoreCase(stateDepMethod) || Objects.equals(options.method, "statedep")) {
            line_debug(options.verbose, "ENV post: updating state-dependent transition rates");
            boolean anyRateUpdated = false;
            for (int e = 0; e < E; e++) {
                for (int h = 0; h < E; h++) {
                    if (envObj.env[e][h] != null && envObj.env[e][h] instanceof Markovian) {
                        // If not defined, rates are left unchanged
                        if (resetEnvRates[e][h] != null) {
                            envObj.env[e][h] =
                                    resetEnvRates[e][h].reset(
                                            (Markovian) envObj.env[e][h], QExit[e][h], UExit[e][h], TExit[e][h]);
                            anyRateUpdated = true;
                        }
                    }
                }
            }
            // see _kb/06-solver-catalog.md (JAR-only implementation notes: conditional re-init on resetEnvRates callbacks)
            if (anyRateUpdated) {
                envObj.init();
            }
        }
    }

    /**
     * Round fractional queue lengths to integers using the largest remainder method,
     * preserving closed chain populations exactly. Modifies Q in-place.
     */
    private void roundMarginalForDiscreteSolver(Matrix Q, NetworkStruct snRef) {
        jline.lang.state.State.roundMarginalPreservingChains(Q, snRef);
    }

    @Override
    protected void finish() {
        line_debug(options.verbose, "ENV finish: computing CDF-weighted steady-state metrics");

        // Use last iteration — matching MATLAB SolverENV.finish()
        // CDF-weight transient trajectories using holdTime{e} (total sojourn time distribution)
        int it = results.size();
        int M = sn[0].nstations;
        int K = sn[0].nclasses;
        int E = getNumberOfModels();

        Matrix[] QExit = new Matrix[E];
        Matrix[] UExit = new Matrix[E];
        Matrix[] TExit = new Matrix[E];

        for (int e = 0; e < E; e++) {
            QExit[e] = new Matrix(M, K);
            UExit[e] = new Matrix(M, K);
            TExit[e] = new Matrix(M, K);

            SolverResult resE = results.get(it).get(e);
            if (resE == null || resE.QNt == null || resE.t == null) {
                continue;
            }

            // Get holdTime MAP for stage e: D0 = holdTime[e].get(0), D1 = holdTime[e].get(1)
            Matrix D0 = envObj.holdTime[e].get(0);
            Matrix D1 = envObj.holdTime[e].get(1);
            map_normalize(D0, D1);
            // Same refinement as post(): the sum runs over the grid the weight
            // can see, not over the one the integrator chose for the horizon.
            CdfGrid grid = refineForCdf(resE.t, D0, D1);

            for (int i = 0; i < M; i++) {
                for (int r = 0; r < K; r++) {
                    if (resE.QNt[i] == null || resE.QNt[i][r] == null) {
                        continue;
                    }
                    int tR = grid.t.getNumRows();
                    if (tR < 2) {
                        continue;
                    }

                    // Compute CDF weights: w = [0, cdf(t2)-cdf(t1), cdf(t3)-cdf(t2), ...]
                    // matching MATLAB: w{e} = [0, map_cdf(holdTime{e}, t(2:end)) - map_cdf(holdTime{e}, t(1:end-1))]'
                    Matrix tLater = Matrix.extractRows(grid.t, 1, tR, null);
                    Matrix tEarlier = Matrix.extractRows(grid.t, 0, tR - 1, null);
                    Matrix cdfLater = map_cdf(D0, D1, tLater);
                    Matrix cdfEarlier = map_cdf(D0, D1, tEarlier);
                    Matrix cdfDiff = cdfLater.sub(1, cdfEarlier);
                    cdfDiff = cdfDiff.transpose();

                    // Prepend 0 for t=0
                    Matrix w = new Matrix(tR, 1);
                    w.set(0, 0, 0.0);
                    for (int j = 0; j < tR - 1; j++) {
                        w.set(j + 1, 0, cdfDiff.get(j, 0));
                    }

                    double wSum = w.elementSum();
                    if (wSum > 0 && !w.hasNaN()) {
                        // QExit{e}(i,r) = QNt' * w / sum(w)
                        QExit[e].set(i, r, grid.on(resE.QNt[i][r]).transpose().mult(w, null).value() / wSum);

                        if (resE.UNt != null && resE.UNt[i] != null && resE.UNt[i][r] != null) {
                            UExit[e].set(i, r, grid.on(resE.UNt[i][r]).transpose().mult(w, null).value() / wSum);
                        }
                        if (resE.TNt != null && resE.TNt[i] != null && resE.TNt[i][r] != null) {
                            TExit[e].set(i, r, grid.on(resE.TNt[i][r]).transpose().mult(w, null).value() / wSum);
                        }
                    }
                }
            }
        }

        // Aggregate: Qval = sum_e probEnv(e) * QExit{e}
        this.result.QN = new Matrix(M, K);
        this.result.UN = new Matrix(M, K);
        this.result.TN = new Matrix(M, K);

        for (int e = 0; e < E; e++) {
            double p = envObj.probEnv.get(0, e);
            for (int i = 0; i < M; i++) {
                for (int k = 0; k < K; k++) {
                    this.result.QN.set(i, k, this.result.QN.get(i, k) + QExit[e].get(i, k) * p);
                    this.result.UN.set(i, k, this.result.UN.get(i, k) + UExit[e].get(i, k) * p);
                    this.result.TN.set(i, k, this.result.TN.get(i, k) + TExit[e].get(i, k) * p);
                }
            }
        }

        this.result.XN = new Matrix(1, K);
        for (int k = 0; k < K; k++) {
            this.result.XN.set(0, k, this.result.TN.get(ref, k));
        }

        // Cache-hit aggregation across the environment for fluid inner solvers.
        aggregateCacheMeanfield();

        result.runtime = (System.nanoTime() - startTime) / 1000000000.0;
        if (options.verbose != VerboseLevel.SILENT) {
            // THE COMPLETION LINE NAMES THE SOLVER, in the shape every other JAR
            // solver uses (NetworkSolver's banner). The former wording named the
            // coupling ("blending completed in N iterations") and no solver at
            // all, so any reader that pairs a table with the solver its banner
            // names -- the parity comparator among them -- saw an UNLABELLED
            // table and adopted it onto whatever golden key was left over. On
            // renv_basic that was the example's standalone per-stage MVA table,
            // i.e. the coupled answer scored against a single stage's.
            jline.io.LineConsole.deferPrint(
                    "ENV analysis [method: %s; type: approximate, deterministic; lang: java; "
                            + "env: %s] completed in %fs. Iterations: %d.\n",
                    (options.method == null || options.method.isEmpty()) ? "default" : options.method,
                    System.getProperty("java.version"), result.runtime, it);
            System.out.flush();
        }
    }

    /**
     * Cache-hit aggregation for fluid inner solvers. Runs a self-contained
     * mean-field fixed point that carries each cache's mean occupancy across
     * environment switches (the cache analog of the queue-length handoff): each
     * stage is integrated over its sojourn from an entry occupancy that mixes
     * the exit occupancy of its predecessors by probOrig. At convergence the
     * per-class hit throughput is the probEnv-weighted, sojourn-averaged
     * arrival x hit-prob, and the reported hit ratio is hit/(hit+miss), written
     * onto the reference model's Cache nodes. No-op unless every stage uses a
     * fluid inner solver and the reference struct has a Cache node. Mirrors the
     * MATLAB aggregateCacheMeanfield_.
     */
    private void aggregateCacheMeanfield() {
        int E = getNumberOfModels();
        for (int e = 0; e < E; e++) {
            if (!(envSolvers[e] instanceof SolverFluid)) {
                return;
            }
        }
        List<Integer> caches = new ArrayList<Integer>();
        for (int i = 0; i < sn[0].nodetype.size(); i++) {
            if (sn[0].nodetype.get(i) == jline.lang.constant.NodeType.Cache) {
                caches.add(i);
            }
        }
        if (caches.isEmpty()) {
            return;
        }
        int ncaches = caches.size();
        int K = sn[0].nclasses;
        int nPoints = 200;

        // Finite integration window per stage (fall back to a few mean holding
        // times when the inner solver left the timespan open).
        double[] time = new double[E];
        for (int e = 0; e < E; e++) {
            double[] ts = envSolvers[e].options.timespan;
            if (ts != null && ts.length >= 2 && Double.isFinite(ts[1]) && ts[1] > 0) {
                time[e] = ts[1];
            } else {
                time[e] = 20.0 * map_mean(envObj.holdTime[e].get(0), envObj.holdTime[e].get(1));
            }
        }

        // Fixed point over per-stage entry occupancy.
        double[][][] entryOcc = new double[E][ncaches][];
        CacheTranResult[] res = new CacheTranResult[E];
        double[][][] wmass = new double[E][ncaches][];
        int maxSweep = Math.max(1, options.iter_max);
        double tol = options.iter_tol;
        double[] prevFlat = null;
        for (int sweep = 0; sweep < maxSweep; sweep++) {
            double[][][] exitOcc = new double[E][ncaches][];
            for (int e = 0; e < E; e++) {
                res[e] = FluidCacheTran.tran(sn[e], time[e], nPoints, entryOcc[e]);
                for (int c = 0; c < ncaches; c++) {
                    double Lam = 0.0;
                    for (int k = 0; k < K; k++) {
                        Lam += res[e].arate[c][k];
                    }
                    wmass[e][c] = sojournWeights(res[e].t, e, Lam);
                    double sw = 0.0;
                    for (double w : wmass[e][c]) {
                        sw += w;
                    }
                    double[][] xo = res[e].xocc[c];
                    if (xo == null || !(sw > 0)) {
                        exitOcc[e][c] = null;
                        continue;
                    }
                    int dim = xo.length;
                    int nt = res[e].t.length;
                    double[] ex = new double[dim];
                    for (int d = 0; d < dim; d++) {
                        double s = 0.0;
                        for (int t = 0; t < nt; t++) {
                            s += xo[d][t] * wmass[e][c][t];
                        }
                        ex[d] = s / sw;
                    }
                    exitOcc[e][c] = ex;
                }
            }
            double[][][] newEntry = new double[E][ncaches][];
            for (int e = 0; e < E; e++) {
                for (int c = 0; c < ncaches; c++) {
                    double[] acc = null;
                    for (int hh = 0; hh < E; hh++) {
                        double po = envObj.probOrig.get(hh, e);
                        if (po > 0 && exitOcc[hh][c] != null) {
                            if (acc == null) {
                                acc = new double[exitOcc[hh][c].length];
                            }
                            for (int d = 0; d < acc.length; d++) {
                                acc[d] += po * exitOcc[hh][c][d];
                            }
                        }
                    }
                    newEntry[e][c] = acc;
                }
            }
            List<Double> flatList = new ArrayList<Double>();
            for (int e = 0; e < E; e++) {
                for (int c = 0; c < ncaches; c++) {
                    if (newEntry[e][c] != null) {
                        for (double v : newEntry[e][c]) {
                            flatList.add(v);
                        }
                    }
                }
            }
            double[] flat = new double[flatList.size()];
            for (int i = 0; i < flat.length; i++) {
                flat[i] = flatList.get(i);
            }
            entryOcc = newEntry;
            if (prevFlat != null && prevFlat.length == flat.length) {
                double mx = 0.0;
                for (int i = 0; i < flat.length; i++) {
                    mx = Math.max(mx, Math.abs(flat[i] - prevFlat[i]));
                }
                if (mx < tol) {
                    prevFlat = flat;
                    break;
                }
            }
            prevFlat = flat;
        }

        // Aggregate converged hit/miss throughputs and write onto the reference
        // model's cache nodes.
        Network refModel = (Network) envModels[0];
        for (int c = 0; c < ncaches; c++) {
            double[] hitT = new double[K];
            double[] missT = new double[K];
            for (int e = 0; e < E; e++) {
                double[] w_ec = wmass[e][c];
                if (w_ec == null) {
                    continue;
                }
                double sw = 0.0;
                boolean bad = false;
                for (double w : w_ec) {
                    sw += w;
                    if (Double.isNaN(w)) {
                        bad = true;
                    }
                }
                if (!(sw > 0) || bad) {
                    continue;
                }
                double pe = envObj.probEnv.get(0, e);
                int nt = res[e].t.length;
                for (int k = 0; k < K; k++) {
                    double a = res[e].arate[c][k];
                    if (a <= 0) {
                        continue;
                    }
                    double hbar = 0.0;
                    double mbar = 0.0;
                    for (int t = 0; t < nt; t++) {
                        hbar += res[e].hitprob[c][k][t] * w_ec[t];
                        mbar += res[e].missprob[c][k][t] * w_ec[t];
                    }
                    hbar /= sw;
                    mbar /= sw;
                    hitT[k] += pe * a * hbar;
                    missT[k] += pe * a * mbar;
                }
            }
            Matrix hitVec = new Matrix(1, K);
            Matrix missVec = new Matrix(1, K);
            for (int k = 0; k < K; k++) {
                double tot = hitT[k] + missT[k];
                if (tot > 0) {
                    hitVec.set(0, k, hitT[k] / tot);
                    missVec.set(0, k, missT[k] / tot);
                } else {
                    hitVec.set(0, k, Double.NaN);
                    missVec.set(0, k, Double.NaN);
                }
            }
            jline.lang.nodes.Node node = refModel.getNodes().get(caches.get(c));
            if (node instanceof jline.lang.nodes.Cache) {
                ((jline.lang.nodes.Cache) node).setResultHitProb(hitVec);
                ((jline.lang.nodes.Cache) node).setResultMissProb(missVec);
            }
        }
    }

    /**
     * Sojourn (residence-time) weights over a trajectory grid for stage e. For
     * a Markovian environment the holding-time-mass weighting equals the
     * time-average of the metric over the stage. Mirrors the CDF-mass weights
     * used in finish().
     */
    private double[] sojournWeights(double[] t, int e, double scale) {
        int nt = (t == null) ? 0 : t.length;
        double[] w = new double[nt];
        if (nt < 2) {
            return w;
        }
        // The RMF drift evolves in normalized (per-request) time; dividing by the
        // cache's total request rate maps it to real time so the holding-time
        // weighting is consistent across phases with different arrival rates.
        double s = (scale > 0) ? scale : 1.0;
        Matrix D0 = envObj.holdTime[e].get(0).copy();
        Matrix D1 = envObj.holdTime[e].get(1).copy();
        map_normalize(D0, D1);
        Matrix tLater = new Matrix(nt - 1, 1);
        Matrix tEarlier = new Matrix(nt - 1, 1);
        for (int i = 0; i < nt - 1; i++) {
            tEarlier.set(i, 0, t[i] / s);
            tLater.set(i, 0, t[i + 1] / s);
        }
        Matrix cdfLater = map_cdf(D0, D1, tLater);
        Matrix cdfEarlier = map_cdf(D0, D1, tEarlier);
        Matrix cdfDiff = cdfLater.sub(1, cdfEarlier).transpose();
        w[0] = 0.0;
        for (int i = 0; i < nt - 1; i++) {
            w[i + 1] = cdfDiff.get(i, 0);
        }
        return w;
    }

    public String getName() {
        return "SolverENV";
    }

    public EnvGeneratorResult getGenerator() {
        int E = getNumberOfModels();
        Matrix[] stageInfGen =  new Matrix[E];
        MatrixCell[] stageEventFilt = new MatrixCell[E];
        @SuppressWarnings("unchecked")
        Map<Integer, Sync>[] stageEvents = (Map<Integer, Sync>[]) new Map[E];

        for (int e = 0; e < E; e++) {
            if (envSolvers[e] instanceof SolverCTMC) {
                stageInfGen[e] = ((SolverCTMC) envSolvers[e]).getGenerator().infGen;
                stageEventFilt[e] = ((SolverCTMC) envSolvers[e]).getGenerator().eventFilt;
                stageEvents[e] = ((SolverCTMC) envSolvers[e]).getGenerator().ev;
            }
            else {
                throw new RuntimeException(
                        "This method requires SolverENV to be instantiated with the CTMC solver.");
            }
        }


        int[] nstates = new int[E];
        for (int e = 0; e < E; e++) {
            nstates[e] = stageInfGen[e].getNumRows();
        }

        Matrix nphases = new Matrix(E, E);
        for (int i = 0; i < E; i++) {
            for (int j = 0; j < E; j++) {
                if (envObj.env[i][j] != null && envObj.env[i][j] instanceof Markovian) {
                    nphases.set(i, j, ((Markovian) envObj.env[i][j]).getNumberOfPhases());
                } else {
                    nphases.set(i, j, 1);
                }
            }
        }
        for (int i = 0; i < E; i++) {
            nphases.set(i, i, nphases.get(i, i) - 1);
        }


        MatrixCell[] renvInfGen = new MatrixCell[E];
        for (int e = 0; e < E; e++) {
            renvInfGen[e] = new MatrixCell(E);
            renvInfGen[e].set(e, stageInfGen[e].copy());
            for (int h = 0; h < E; h++) {
                if (h != e) {
                    Matrix resetMatrixEH = new Matrix(nstates[e], nstates[h]);
                    resetMatrixEH.fill(0.0);
                    int minStates = FastMath.min(nstates[e], nstates[h]);
                    for (int i = 0; i < minStates; i++) {
                        resetMatrixEH.set(i, i, 1.0);
                    }
                    renvInfGen[e].set(h, resetMatrixEH.copy());
                }
            }
        }


        List<RenvEvent> renvEvents = new ArrayList<>();
        for (int e = 0; e < E; e++) {
            for (int h = 0; h < E; h++) {
                if (h != e) {
                    Matrix D0;
                    if (envObj.env[e][h] == null) {
                        D0 = new Matrix(new int[] {0});
                    }
                    else {
                        D0 = envObj.env[e][h].getProcess().get(0);
                    }
                    renvInfGen[e].set(e, renvInfGen[e].get(e).krons(D0));

                    Matrix pie, D1;
                    if (envObj.env[h][e] == null || map_pie(envObj.env[h][e].getProcess()).hasNaN()) {
                        pie = new Matrix((int) nphases.get(h, e), 1);
                        pie.fill(1.0);
                    }
                    else {
                        pie = map_pie(envObj.env[h][e].getProcess());
                    }
                    if (envObj.env[e][h] == null) {
                        D1 = new Matrix(new int[] {0});
                    }
                    else {
                        D1 = envObj.env[e][h].getProcess().get(1);
                    }
                    Matrix onePhase = Matrix.ones((int) nphases.get(e, h), 1);
                    onePhase.multEq(pie);
                    Matrix kronArg = D1.mult(onePhase, null);
                    renvInfGen[e].set(h, renvInfGen[e].get(h).kron(kronArg));

                    for (int i = 0; i < nodesOf(envModels[e]); i++) {
                        RenvEvent ev = new RenvEvent(i, -1, Double.NaN, new Matrix(0, 0), Double.NaN, Double.NaN, new Pair<>(e, h));
                        renvEvents.add(ev);
                    }

                    for (int f = 0; f < E; f++) {
                        if (f != h && f != e) {
                            Matrix pie_fh;
                            if (envObj.env[f][h] == null || map_pie(envObj.env[f][h].getProcess()).hasNaN()) {
                                pie_fh = new Matrix((int) nphases.get(f, h), 1);
                                pie_fh.fill(1.0);
                            }
                            else {
                                pie_fh = map_pie(envObj.env[f][h].getProcess());
                            }
                            Matrix oneVec = new Matrix((int) nphases.get(e, h), 1);
                            oneVec.fill(1.0);
                            renvInfGen[e].set(f, renvInfGen[e].get(f).kron(oneVec.mult(pie_fh, null)));

                        }
                    }
                }
            }
        }

        Matrix[][] renvEventFilt = new Matrix[E][E];
        for (int e = 0; e < E; e++) {
            for (int h = 0; h < E; h++) {
                MatrixCell[] tmpCell = new MatrixCell[E];
                for (int e1 = 0; e1 < E; e1++) {
                    tmpCell[e1] = new MatrixCell(renvInfGen[e1]);
                    tmpCell[e1].set(e1, tmpCell[e1].get(e1).fill(0));
                    for (int h1 = 0; h1 < E; h1++) {
                        if (e != e1 && h != h1) {
                            tmpCell[e1].set(h1, tmpCell[e1].get(h1).fill(0));
                        }
                    }
                }
                Matrix flatTmpCell = flatMatrix(tmpCell);
                renvEventFilt[e][h] = flatTmpCell;
            }
        }

        Matrix newRenvInfGen = flatMatrix(renvInfGen);
        newRenvInfGen = ctmc_makeinfgen(newRenvInfGen);

        return new EnvGeneratorResult(stageInfGen, newRenvInfGen, stageEventFilt, renvEventFilt, stageEvents, renvEvents);
    }

    public void setRef(int i) {
        this.ref = i;
    }

    /**
     * Container class holding the generator matrices and related data structures
     * for both stage-specific and random environment transitions.
     * Used to encapsulate results from environment-aware model generation.
     */
    public static class EnvGeneratorResult {
        public Matrix[] stageInfGen;
        public Matrix renvInfGen;
        public MatrixCell[] stageEventFilt;
        public Matrix[][] renvEventFilt;
        public Map<Integer, Sync>[] stageEvents;
        public List<RenvEvent> renvEvents;

        public EnvGeneratorResult(Matrix[] stageInfGen, Matrix renvInfGen, MatrixCell[] stageEventFilt,
                                  Matrix[][] renvEventFilt, Map<Integer, Sync>[] stageEvents,
                                  List<RenvEvent> renvEvents) {
            this.stageInfGen = stageInfGen;
            this.renvInfGen = renvInfGen;
            this.stageEventFilt = stageEventFilt;
            this.renvEventFilt = renvEventFilt;
            this.stageEvents = stageEvents;
            this.renvEvents = renvEvents;
        }
    }

    /**
     * Result container for sample path analysis, containing metrics for each segment.
     */
    public static class SamplePathResult {
        public List<SamplePathSegment> segments;

        public SamplePathResult() {
            this.segments = new ArrayList<>();
        }

        /**
         * Data for a single segment in the sample path.
         */
        public static class SamplePathSegment {
            public int segmentIndex;
            public int stageIndex;
            public String stageName;
            public double duration;
            public Matrix initialQ, initialU, initialT;
            public Matrix finalQ, finalU, finalT;
            public Matrix[][] QNt, UNt, TNt;
            public Matrix t;
        }
    }


    /**
     * Flattens a row of MatrixCell objects into a single large Matrix.
     * Each MatrixCell represents a mapping from column indices to Matrix blocks.
     * The resulting Matrix is constructed by concatenating the blocks in the correct positions.
     *
     * @param blocks A row (array) of MatrixCell objects.
     * @return The flattened Matrix.
     */
    private Matrix flatMatrix(MatrixCell[] blocks) {
        int blockRows = blocks.length;
        int blockCols = blocks[0].size();

        int[] rowHeights = new int[blockRows];
        int[] colWidths = new int[blockCols];
        int totalRows = 0;
        int totalCols = 0;

        for (int i = 0; i < blockRows; i++) {
            rowHeights[i] = blocks[i].get(0).getNumRows();
            totalRows += rowHeights[i];
        }

        for (int j = 0; j < blockCols; j++) {
            colWidths[j] = blocks[0].get(j).getNumCols();
            totalCols += colWidths[j];
        }

        Matrix result = new Matrix(totalRows, totalCols);

        int rowOffset = 0;
        for (int i = 0; i < blockRows; i++) {
            int colOffset = 0;
            for (int j = 0; j < blockCols; j++) {
                Matrix block = blocks[i].get(j);
                for (int r = 0; r < block.getNumRows(); r++) {
                    for (int c = 0; c < block.getNumCols(); c++) {
                        result.set(rowOffset + r, colOffset + c, block.get(r, c));
                    }
                }
                colOffset += colWidths[j];
            }
            rowOffset += rowHeights[i];
        }

        return result;
    }


    public void getAvg() {
        getEnsembleAvg();
    }

    @Override
    public AvgTable getEnsembleAvg() {
        if (this.result == null || this.result.QN == null || this.result.QN.isEmpty() || this.options.force) {
            // Solver console: ENV reaches its analysis through blending(), not
            // through EnsembleSolver.iterate, so the narrated run is opened
            // here -- the point every accessor of the averages goes through.
            jline.io.LineConsole.beginRun(this, this.options);
            try {
                blending();
            } finally {
                jline.io.LineConsole.closeRun(this);
            }
        }

        List<Double> Qval = this.result.QN.toList1D();
        List<Double> Uval = this.result.UN.toList1D();
        List<Double> Tval = this.result.TN.toList1D();
        List<Double> Rval = new ArrayList<>();
        List<Double> Residval = new ArrayList<>();
        List<Double> Aval = new ArrayList<>();
        for (int i = 0; i < Qval.size(); i++) {
            Residval.add(0.0);
            Aval.add(0.0);
            if (Tval.get(i) == 0 && Qval.get(i) == 0) {
                Rval.add(0.0);
            } else if (Tval.get(i) == 0) {
                Rval.add(Inf);
            } else {
                Rval.add(Qval.get(i) / Tval.get(i));
            }
        }

        AvgTable avgTable = new NetworkAvgTable(Qval, Uval, Rval, Residval, Aval, Tval);
        avgTable.setOptions(this.options);
        return avgTable;
    }

    /**
     * Computes transient performance metrics for a sample path through environment states.
     * The method runs transient analysis for each segment and extracts initial and final metric values.
     *
     * @param samplePath List of entries where each entry is an Object[] containing:
     *                   - entry[0]: stage identifier (Integer for 0-based index, or String for stage name)
     *                   - entry[1]: duration (Double, time spent in that stage)
     * @return SamplePathResult containing metrics for each segment
     * @throws IllegalArgumentException if sample path is empty, stage not found, or duration non-positive
     *
     * Example:
     * <pre>
     * List<Object[]> path = new ArrayList<>();
     * path.add(new Object[]{"Fast", 5.0});
     * path.add(new Object[]{"Slow", 10.0});
     * path.add(new Object[]{"Fast", 3.0});
     * SamplePathResult result = solver.getSamplePathTable(path);
     * </pre>
     */
    public SamplePathResult getSamplePathTable(List<Object[]> samplePath) {
        if (samplePath == null || samplePath.isEmpty()) {
            throw new IllegalArgumentException("Sample path cannot be empty.");
        }

        // Initialize if needed
        if (envObj.probEnv == null || envObj.probEnv.isEmpty()) {
            init();
        }

        int E = getNumberOfModels();
        int M = sn[0].nstations;
        int K = sn[0].nclasses;

        SamplePathResult result = new SamplePathResult();

        // Initialize queue lengths (uniform distribution for closed classes)
        Matrix Q_current = new Matrix(M, K);
        for (int k = 0; k < K; k++) {
            double njobs = sn[0].njobs.get(k);
            // Only initialize from njobs for closed classes (finite population)
            if (Double.isFinite(njobs) && njobs > 0) {
                for (int i = 0; i < M; i++) {
                    Q_current.set(i, k, njobs / M);
                }
            }
        }

        int segIdx = 0;
        for (Object[] entry : samplePath) {
            Object stageSpec = entry[0];
            double duration = (Double) entry[1];

            // Resolve stage index
            int e;
            String stageName;
            if (stageSpec instanceof String) {
                e = envObj.findStageByName((String) stageSpec);
                if (e == -1) {
                    throw new IllegalArgumentException("Stage not found: " + stageSpec);
                }
                stageName = (String) stageSpec;
            } else {
                e = (Integer) stageSpec;
                if (e < 0 || e >= E) {
                    throw new IllegalArgumentException("Stage index out of range [0, " + (E - 1) + "]: " + e);
                }
                stageName = envObj.getStageName(e);
            }

            if (duration <= 0) {
                throw new IllegalArgumentException("Duration must be positive.");
            }

            // Initialize from current queue lengths
            initFromMarginalOf(envModels[e], Q_current);

            // Set solver timespan and run transient analysis
            envSolvers[e].options.timespan = new double[]{0, duration};
            envSolvers[e].reset();

            // Build initial state vector from Q_current
            Matrix initial = new Matrix(1, M * K);
            int idx = 0;
            for (int i = 0; i < M; i++) {
                for (int k = 0; k < K; k++) {
                    initial.set(idx++, Q_current.get(i, k));
                }
            }

            // Run the analysis based on solver type
            SolverResult stageResult;
            if (envSolvers[e] instanceof SolverFluid) {
                SolverFluid solverFluid = (SolverFluid) envSolvers[e];
                solverFluid.options.init_sol = initial;
                // Sync sn.state from init_sol so SolverFluid sees the real initial queue lengths
                NetworkStruct sn_fl = solverFluid.model.getStruct(false);
                List<StatefulNode> statefulNodes = sn_fl.stateful;
                int offset = 0;
                for (StatefulNode node : statefulNodes) {
                    int cols = sn_fl.nclasses;
                    Matrix nodeState = Matrix.extract(initial, 0, 1, offset, offset + cols);
                    sn_fl.state.put(node, nodeState);
                    int nodeIdx = node.getNodeIndex();
                    int stationIdx = (int) sn_fl.nodeToStation.get(0, nodeIdx);
                    solverFluid.model.getStations().get(stationIdx).setState(nodeState);
                    offset += cols;
                }
                solverFluid.sn = sn_fl;
                stageResult = solverFluid.runMethodSpecificAnalyzer();
            } else if (envSolvers[e] instanceof SolverLDES) {
                SolverLDES solverLDES = (SolverLDES) envSolvers[e];
                solverLDES.options.init_sol = initial;
                // Sync sn.state from init_sol so SolverLDES sees the real initial queue lengths
                NetworkStruct sn_des = solverLDES.model.getStruct(false);
                List<StatefulNode> statefulNodes = sn_des.stateful;
                int offset = 0;
                for (StatefulNode node : statefulNodes) {
                    int cols = sn_des.nclasses;
                    Matrix nodeState = Matrix.extract(initial, 0, 1, offset, offset + cols);
                    sn_des.state.put(node, nodeState);
                    int nodeIdx = node.getNodeIndex();
                    int stationIdx = (int) sn_des.nodeToStation.get(0, nodeIdx);
                    solverLDES.model.getStations().get(stationIdx).setState(nodeState);
                    offset += cols;
                }
                solverLDES.sn = sn_des;
                stageResult = solverLDES.runMethodSpecificAnalyzer();
            } else {
                throw new UnsupportedOperationException(
                    "getSamplePathTable requires SolverFluid or SolverLDES as stage solver, got: " +
                    envSolvers[e].getClass().getSimpleName());
            }

            // Extract metrics
            SamplePathResult.SamplePathSegment seg = new SamplePathResult.SamplePathSegment();
            seg.segmentIndex = segIdx;
            seg.stageIndex = e;
            seg.stageName = stageName;
            seg.duration = duration;
            seg.QNt = stageResult.QNt;
            seg.UNt = stageResult.UNt;
            seg.TNt = stageResult.TNt;
            seg.t = stageResult.t;

            seg.initialQ = new Matrix(M, K);
            seg.initialU = new Matrix(M, K);
            seg.initialT = new Matrix(M, K);
            seg.finalQ = new Matrix(M, K);
            seg.finalU = new Matrix(M, K);
            seg.finalT = new Matrix(M, K);

            for (int i = 0; i < M; i++) {
                for (int k = 0; k < K; k++) {
                    if (stageResult.QNt != null && stageResult.QNt[i][k] != null && !stageResult.QNt[i][k].isEmpty()) {
                        int nRows = stageResult.QNt[i][k].getNumRows();
                        seg.initialQ.set(i, k, stageResult.QNt[i][k].get(0, 0));
                        seg.finalQ.set(i, k, stageResult.QNt[i][k].get(nRows - 1, 0));
                    }
                    if (stageResult.UNt != null && stageResult.UNt[i][k] != null && !stageResult.UNt[i][k].isEmpty()) {
                        int nRows = stageResult.UNt[i][k].getNumRows();
                        seg.initialU.set(i, k, stageResult.UNt[i][k].get(0, 0));
                        seg.finalU.set(i, k, stageResult.UNt[i][k].get(nRows - 1, 0));
                    }
                    if (stageResult.TNt != null && stageResult.TNt[i][k] != null && !stageResult.TNt[i][k].isEmpty()) {
                        int nRows = stageResult.TNt[i][k].getNumRows();
                        seg.initialT.set(i, k, stageResult.TNt[i][k].get(0, 0));
                        seg.finalT.set(i, k, stageResult.TNt[i][k].get(nRows - 1, 0));
                    }
                }
            }

            result.segments.add(seg);
            Q_current = seg.finalQ.copy();
            segIdx++;
        }

        return result;
    }


    @SuppressWarnings("unchecked")
    private void blending() {
        if (options != null && options.method != null
                && (options.method.equalsIgnoreCase("avg") || options.method.equalsIgnoreCase("dec"))) {
            solveEnvLimit();
            return;
        }
        this.statevecMethod = options != null && options.method != null
                && options.method.equalsIgnoreCase("statevec");
        if (this.statevecMethod) {
            // The state-vector analyzer needs a single per-stage CTMC generator;
            // an LQN decomposes into multiple layer submodels and exposes no
            // single generator, so it is supported only by the mean-field path.
            for (int e = 0; e < getNumberOfModels(); e++) {
                if (envModels[e] instanceof LayeredNetwork) {
                    throw new RuntimeException(String.format(
                        "The state-vector (statevec) analyzer does not support LayeredNetwork stages (stage %d): "
                        + "an LQN has no single stage generator. Use the default mean-field analyzer.", e));
                }
            }
            blendingStatevec();
            return;
        }
        init();
        // Match MATLAB iterate(): call pre(1) to run initial ODE for each stage,
        // extract steady-state, and initialize models via initFromMarginal().
        // Without this, the iteration loop starts from an arbitrary initial state
        // instead of the ODE steady-state, causing slow/failed convergence.
        pre(1);
        int max_iter = options.iter_max;
        startTime = System.nanoTime();

       int E = getNumberOfModels();
       int M = sn[0].nstations;
       int K = sn[0].nclasses;

        this.results = new HashMap<>();

        for (int it = 1; it <= max_iter; it++) {
            UNtStages = new ArrayList<>();
            tStages = new MatrixCell(E);
            for (int e = 0; e < E; e++) {
                SolverResult results_it_e = analyze(it, e);
                if (results.containsKey(it)) {
                    results.get(it).put(e, results_it_e);
                }
                else {
                    Map<Integer, SolverResult> map = new HashMap<>();
                    map.put(e, results_it_e);
                    results.put(it, map);
                }
            }

            // Call post() to compute CDF-weighted exit metrics, entry states,
            // and re-initialize models for next iteration (matching MATLAB iterate loop)
            post(it);

            // State-dependent updates (resetEnvRates + envObj.init()) are handled inside post(),
            // matching MATLAB's EnsembleSolver iterate loop. No additional Eutil scaling needed here.


            if (converged(it)) {
                break;
            }

            if (it == max_iter) {
                int itLast = (int) round(max_iter * 0.9);
                for (int e = 0; e < E; e++) {
                    Matrix tmpQN = new Matrix(M, K);
                    for (int itTmp = itLast; itTmp <= it; itTmp++) {
                        tmpQN = tmpQN.add(1, results.get(itTmp).get(e).QN);
                    }
                    tmpQN.scaleEq(1.0 / (it - itLast + 1));
                    this.results.get(it).get(e).QN = tmpQN;
                }
            }
        }

        finish();
    }


    /**
     * Closed-form fast/slow random-environment limits (options.method='avg' or
     * 'dec'), which require no inter-stage coupling iteration and therefore no
     * transient analysis from the stage solvers.
     *
     * <p>'avg' is the fast-environment limit: the base model sees the
     * stationary-probability-weighted average of the modulated rates, so a
     * single rate-averaged model is built and solved once. Exact as the
     * stage-switching rate tends to infinity.
     *
     * <p>'dec' is the slow-environment (quasi-stationary) decomposition: each
     * stage is solved independently in steady state and the per-stage metrics
     * are averaged with weights probEnv(e). Exact as the stage-switching rate
     * tends to zero.
     *
     * <p>Mirrors matlab/src/solvers/@SolverENV/SolverENV.m solveEnvLimit().
     */
    private void solveEnvLimit() {
        startTime = System.nanoTime();
        init();
        int E = getNumberOfModels();
        int M = sn[0].nstations;
        int K = sn[0].nclasses;
        Matrix Qval = new Matrix(M, K);
        Matrix Uval = new Matrix(M, K);
        Matrix Tval = new Matrix(M, K);

        if (options.method.equalsIgnoreCase("dec")) {
            for (int e = 0; e < E; e++) {
                double p = envObj.probEnv.get(0, e);
                Solver se = envSolvers[e];
                se.reset();
                SolverResult re = ((NetworkSolver) se).getAvg();
                accumulate(Qval, re.QN, p);
                accumulate(Uval, re.UN, p);
                accumulate(Tval, re.TN, p);
            }
        } else {
            Network avgModel = buildRateAveragedModel();
            NetworkSolver innerSolver = newStageSolver(envSolvers[0], avgModel);
            SolverResult re = innerSolver.getAvg();
            accumulate(Qval, re.QN, 1.0);
            accumulate(Uval, re.UN, 1.0);
            accumulate(Tval, re.TN, 1.0);
        }

        this.result = new SolverResult();
        this.result.QN = Qval;
        this.result.UN = Uval;
        this.result.TN = Tval;
        this.result.XN = new Matrix(1, K);
        for (int k = 0; k < K; k++) {
            this.result.XN.set(0, k, Tval.get((int) sn[0].refstat.get(k, 0), k));
        }
        this.result.runtime = (System.nanoTime() - startTime) / 1000000000.0;
    }

    /** dest += weight * src, tolerating a null or mis-shaped source. */
    private void accumulate(Matrix dest, Matrix src, double weight) {
        if (src == null || src.isEmpty()) {
            return;
        }
        for (int i = 0; i < dest.getNumRows(); i++) {
            for (int k = 0; k < dest.getNumCols(); k++) {
                dest.set(i, k, dest.get(i, k) + weight * src.get(i, k));
            }
        }
    }

    /**
     * Fast-environment model: replace every environment-modulated (i.e.
     * stage-varying) station rate by its probEnv-weighted average, represented
     * as an exponential. Non-modulated parameters keep their original
     * distribution, so the base model is preserved exactly outside the
     * modulated rates.
     */
    private Network buildRateAveragedModel() {
        if (!(envModels[0] instanceof Network)) {
            throw new RuntimeException("The rate-averaged (avg) environment limit requires flat Network stages.");
        }
        int E = getNumberOfModels();
        Network avgModel = ((Network) envModels[0]).copy();
        int M = avgModel.getNumberOfStations();
        int K = avgModel.getNumberOfClasses();
        java.util.List<Station> stations = avgModel.getStations();
        java.util.List<JobClass> classes = avgModel.getClasses();
        for (int i = 0; i < M; i++) {
            Station node = stations.get(i);
            // Stateful/absorbing nodes (Cache, Sink) are not Stations in the JAR
            // hierarchy, so they cannot appear here; a station whose rate is NaN
            // in some stage is skipped below.
            for (int k = 0; k < K; k++) {
                double[] r = new double[E];
                boolean usable = true;
                double rmin = Inf;
                double rmax = -Inf;
                for (int e = 0; e < E; e++) {
                    r[e] = sn[e].rates.get(i, k);
                    if (Double.isNaN(r[e]) || r[e] <= 0) {
                        usable = false;
                        break;
                    }
                    rmin = Math.min(rmin, r[e]);
                    rmax = Math.max(rmax, r[e]);
                }
                if (!usable) {
                    continue; // disabled for some stage: leave as configured
                }
                if ((rmax - rmin) <= 1e-12 * Math.max(1.0, rmax)) {
                    continue; // not modulated: keep the original distribution
                }
                double ravg = 0;
                for (int e = 0; e < E; e++) {
                    ravg += envObj.probEnv.get(0, e) * r[e];
                }
                if (node instanceof Source) {
                    ((Source) node).setArrival(classes.get(k), new Exp(ravg));
                } else if (node instanceof ServiceStation) {
                    ((ServiceStation) node).setService(classes.get(k), new Exp(ravg));
                }
            }
        }
        // The copy inherited a cached NetworkStruct from the stage model; force a
        // hard rebuild so the averaged rates take effect.
        avgModel.refreshStruct(true);
        return avgModel;
    }

    /** A stage solver of the same class and options as the template, on another model. */
    private NetworkSolver newStageSolver(Solver template, Network model) {
        try {
            return (NetworkSolver) template.getClass()
                    .getConstructor(Network.class, SolverOptions.class)
                    .newInstance(model, template.options);
        } catch (Exception e) {
            throw new RuntimeException("Cannot instantiate the stage solver " + template.getClass().getSimpleName()
                    + " on the rate-averaged model: " + e.getMessage(), e);
        }
    }

    // =====================================================================
    // State-vector analyzer (options.method='statevec').
    //
    // Carries the full per-stage state distribution across environment
    // switches instead of collapsing it to marginal mean queue lengths.
    // Mirrors matlab/src/solvers/ENV/solver_env_statevec_analyzer.m.
    // =====================================================================
    private void blendingStatevec() {
        init();
        int E = getNumberOfModels();
        svPre();
        int max_iter = options.iter_max;
        startTime = System.nanoTime();
        boolean hasConverged = false;
        for (int it = 1; it <= max_iter; it++) {
            for (int e = 0; e < E; e++) {
                svAnalyze(e);
            }
            svPost();
            if (svConverged()) {
                hasConverged = true;
                break;
            }
        }
        if (!hasConverged) {
            // Exiting on iter_max is a NON-convergence: the last iterate can be
            // far from the fixed point. Returning it silently is what makes a
            // wrong number indistinguishable from a right one.
            line_warning("SolverENV", "The statevec fixed point did not converge in options.iter_max=%d "
                    + "iterations; the returned solution is the last iterate and may be far from the fixed "
                    + "point. Raise options.iter_max, or loosen options.iter_tol only if the residual is "
                    + "already small.", max_iter);
        }
        svFinish();
    }

    // Build the per-stage generator, state space and metric-mapping data once,
    // and initialise the per-stage entry distributions.
    private void svPre() {
        int E = getNumberOfModels();
        svStages = new SvStage[E];
        piEnter = new Matrix[E];
        piExitDest = new Matrix[E][E];
        piTimeAvg = new Matrix[E];
        for (int e = 0; e < E; e++) {
            NetworkSolver solver_e = (NetworkSolver) envSolvers[e];
            SolverOptions opts_e = solver_e.options;
            if (opts_e.timespan == null || opts_e.timespan.length < 2 || !Double.isFinite(opts_e.timespan[1])) {
                throw new RuntimeException(String.format("The statevec analyzer requires a finite inner-solver "
                        + "timespan for stage %d, e.g. CTMC(model,'timespan',[0,T]).", e));
            }
            SvStage st = new SvStage();
            st.timespan = opts_e.timespan;
            if (solver_e instanceof SolverCTMC) {
                ResultCTMC r = Solver_ctmc.solver_ctmc(sn[e], opts_e);
                st.backend = "ctmc";
                st.Q = r.getQ();
                st.SS = r.getStateSpace();
                st.SSaggr = r.getStateSpaceAggr();
                st.arvRates = r.getArvRates();
                st.depRates = r.getDepRates();
                st.snE = r.getSn();
            } else if (solver_e instanceof SolverMAM) {
                Solver_mam_ldqbd_statevec.Ld ld = Solver_mam_ldqbd_statevec.solver_mam_ldqbd_ld(sn[e], opts_e);
                Solver_mam_ldqbd_statevec.Flat f = Solver_mam_ldqbd_statevec.solver_mam_ldqbd_flatten(ld);
                st.backend = "mam";
                st.Q = f.Q;
                st.ld = ld;
                st.levelOf = f.levelOf;
            } else {
                throw new RuntimeException(String.format("The statevec analyzer requires a SolverCTMC or SolverMAM "
                        + "inner solver, but environment stage %d uses %s.", e, solver_e.getClass().getSimpleName()));
            }
            svStages[e] = st;
        }
        svSeedEntryDistributions();
        piEnterPrev = piEnter.clone();
    }

    // Warm start. A stage's OWN stationary distribution is not a usable seed
    // here: a stage that is individually unstable or critical (arrival rate >=
    // its own service rate) has no stationary law at all, and
    // ctmc_solve_reducible then returns the stationary law of the TRUNCATED
    // generator, which piles mass against the truncation wall and whose mean
    // grows linearly with the cutoff (for a critical M/M/1 truncated at N it is
    // uniform, with mean N/2).
    //
    // The fixed point chained in svPost is exact -- it is the stationary
    // equation of the joint (queue,stage) chain,
    // phi_e = (sum_h phi_h q_he) (s_e I - Q_e)^-1 -- and it does contract to
    // the right answer from that seed, but the number of sweeps needed grows
    // with the cutoff. At a finite iter_max the reported result therefore
    // drifts further from the truth as the cutoff is RAISED, i.e. the natural
    // user response to a suspect number makes it worse.
    //
    // Seed instead from the environment-averaged generator sum_e probEnv(e)*Q_e,
    // which is positive recurrent exactly when the model is stable on average --
    // the regime in which the answer exists -- so its stationary law is
    // cutoff-independent. Averaging needs one common state space; when the
    // stages differ in size (resetStateFun is what bridges them) fall back to
    // the per-stage law, which is no worse than before.
    private void svSeedEntryDistributions() {
        int E = getNumberOfModels();
        boolean sameSpace = E > 1;
        for (int e = 1; e < E && sameSpace; e++) {
            sameSpace = svStages[e].Q.getNumRows() == svStages[0].Q.getNumRows();
        }

        Matrix shared = null;
        if (sameSpace) {
            double[] w = new double[E];
            double wsum = 0.0;
            boolean usable = envObj.probEnv != null && envObj.probEnv.getNumCols() >= E;
            for (int e = 0; e < E && usable; e++) {
                double pe = envObj.probEnv.get(0, e);
                if (!Double.isFinite(pe) || pe < 0.0) {
                    usable = false;
                } else {
                    wsum += pe;
                }
            }
            if (usable && wsum > 0.0) {
                for (int e = 0; e < E; e++) {
                    w[e] = envObj.probEnv.get(0, e) / wsum;
                }
            } else {
                // Stage probabilities unavailable or degenerate: weight equally.
                for (int e = 0; e < E; e++) {
                    w[e] = 1.0 / E;
                }
            }
            Matrix Qbar = svStages[0].Q.copy();
            Qbar.zero();
            for (int e = 0; e < E; e++) {
                Qbar = Qbar.add(w[e], svStages[e].Q);
            }
            shared = rowNormalizeNonneg(Ctmc_solve_reducible.ctmc_solve_reducible(Qbar).getLeft());
        }

        for (int e = 0; e < E; e++) {
            if (shared == null) {
                piEnter[e] = rowNormalizeNonneg(
                        Ctmc_solve_reducible.ctmc_solve_reducible(svStages[e].Q).getLeft());
            } else {
                piEnter[e] = shared.copy();
            }
        }
    }

    // Propagate the entry distribution of stage e through its sojourn and store
    // the per-destination exit distributions and the sojourn-end distribution.
    private void svAnalyze(int e) {
        int E = getNumberOfModels();
        SvStage st = svStages[e];
        Matrix Q = st.Q;
        Matrix pi0 = piEnter[e];
        double t0 = st.timespan[0];
        double t1 = st.timespan[1];

        // Deterministic sojourn option (off by default): the stage lasts exactly
        // its mean duration d_e, so the exit equals pi0*exp(Q*d_e) toward every
        // destination, and the blend uses the time average over [0,d_e].
        if (options.sojourn != null && options.sojourn.equalsIgnoreCase("deterministic")) {
            double d_e = Math.max(map_mean(envObj.holdTime[e]), 2.220446049250313e-16);
            Pair<Matrix, Matrix> ta = Ctmc_timeaverage.ctmc_timeaverage(pi0, Q, d_e);
            Matrix piAvg = ta.getLeft();
            Matrix piEx = ta.getRight();
            for (int h = 0; h < E; h++) {
                piExitDest[e][h] = (E0rate.get(e, h) > 0) ? piEx : null;
            }
            piTimeAvg[e] = piAvg;
            return;
        }

        // Exponential environment sojourn (every enabled transition from e is
        // Exp): the exit equals the time average given exactly by the resolvent
        // s*pi*(sI - Q)^{-1}, s = total exit rate, the same toward every dest.
        boolean expSojourn = true;
        for (int h = 0; h < E; h++) {
            if (E0rate.get(e, h) > 0 && !(envObj.env[e][h] instanceof Exp)) {
                expSojourn = false;
                break;
            }
        }
        if (expSojourn) {
            double s_e = 0.0;
            for (int h = 0; h < E; h++) {
                s_e += E0rate.get(e, h);
            }
            int d = Q.getNumRows();
            Matrix A = Matrix.eye(d).scale(s_e).sub(Q);
            Matrix piRes = pi0.mult(A.inv()).scale(s_e);
            for (int h = 0; h < E; h++) {
                piExitDest[e][h] = (E0rate.get(e, h) > 0) ? piRes : null;
            }
            piTimeAvg[e] = piRes;
            return;
        }

        // General Markovian (PH/Erlang) sojourn: transient pi(t) = pi0*exp(Q t),
        // averaged over the random holding-time / transition CDFs.
        Pair<double[], List<double[]>> tr = Ctmc_transient.ctmc_transient(Q, pi0, t0, t1);
        double[] t = tr.getLeft();
        Matrix pit = listToMatrix(tr.getRight());

        for (int h = 0; h < E; h++) {
            if (E0rate.get(e, h) <= 0) {
                piExitDest[e][h] = null;
                continue;
            }
            Matrix w = cdfIncrementWeights(envObj.proc[e][h], t);
            piExitDest[e][h] = weightedAvgRows(pit, w);
        }

        Matrix wh = cdfIncrementWeights(envObj.holdTime[e], t);
        Matrix avg = weightedAvgRows(pit, wh);
        piTimeAvg[e] = (avg != null) ? avg : extractRowVec(pit, pit.getNumRows() - 1);
    }

    // Chain the entry distributions: carry each stage's exit distributions into
    // the stages they feed, weighted by the origin probabilities probOrig.
    private void svPost() {
        int E = getNumberOfModels();
        piEnterPrev = new Matrix[E];
        for (int e = 0; e < E; e++) {
            piEnterPrev[e] = (piEnter[e] != null) ? piEnter[e].copy() : null;
        }
        Matrix[] piEnterNew = new Matrix[E];
        for (int e = 0; e < E; e++) {
            int ns = svStages[e].Q.getNumRows();
            Matrix acc = new Matrix(1, ns);
            acc.zero();
            double wsum = 0.0;
            for (int h = 0; h < E; h++) {
                double po = envObj.probOrig.get(h, e);
                if (po > 0 && piExitDest[h][e] != null) {
                    Matrix pex = asRowVec(envObj.resetStateFun[h][e].reset(piExitDest[h][e]));
                    if (pex.length() != ns) {
                        throw new RuntimeException(String.format("resetStateFun[%d][%d] returned a %d-element vector "
                                + "but stage %d has %d states. Supply a resetStateFun[%d][%d] that maps the state "
                                + "space of stage %d onto that of stage %d.", h, e, pex.length(), e, ns, h, e, h, e));
                    }
                    acc = acc.add(po, pex);
                    wsum += po;
                }
            }
            if (wsum > 0) {
                acc.scaleEq(1.0 / wsum);
            } else {
                acc = piEnter[e].copy();
            }
            piEnterNew[e] = rowNormalizeNonneg(acc);
        }
        piEnter = piEnterNew;
    }

    // Converged when the max L1 change across all entry distributions over a
    // full cycle falls below iter_tol.
    private boolean svConverged() {
        int E = getNumberOfModels();
        if (piEnterPrev == null || piEnter == null) {
            return false;
        }
        double l1 = 0.0;
        for (int e = 0; e < E; e++) {
            Matrix a = piEnter[e];
            Matrix b = piEnterPrev[e];
            if (a == null || b == null || a.length() != b.length()) {
                return false;
            }
            double s = 0.0;
            for (int i = 0; i < a.length(); i++) {
                s += Math.abs(a.get(i) - b.get(i));
            }
            l1 = Math.max(l1, s);
        }
        if (Double.isNaN(l1) || Double.isInfinite(l1)) {
            return false;
        }
        return l1 < options.iter_tol;
    }

    // Environment-averaged blend: map each stage's sojourn-end distribution to
    // marginal metrics and weight by the stage probability probEnv.
    private void svFinish() {
        int E = getNumberOfModels();
        int M = sn[0].nstations;
        int K = sn[0].nclasses;
        Matrix Qval = new Matrix(M, K); Qval.zero();
        Matrix Uval = new Matrix(M, K); Uval.zero();
        Matrix Tval = new Matrix(M, K); Tval.zero();
        for (int e = 0; e < E; e++) {
            Matrix piF = piTimeAvg[e];
            if (piF == null) {
                continue;
            }
            SvStage st = svStages[e];
            Matrix QN;
            Matrix UN;
            Matrix TN;
            if ("mam".equals(st.backend)) {
                Solver_mam_ldqbd_statevec.Avg a = Solver_mam_ldqbd_statevec.solver_mam_ldqbd_avg(st.ld, piF, st.levelOf);
                QN = a.QN; UN = a.UN; TN = a.TN;
            } else {
                Ctmc_avg_from_pi.Result r = Ctmc_avg_from_pi.ctmc_avg_from_pi(st.snE, piF, st.SS, st.SSaggr, st.arvRates, st.depRates);
                QN = r.QN; UN = r.UN; TN = r.TN;
            }
            double pe = envObj.probEnv.get(0, e);
            Qval = Qval.add(pe, QN);
            Uval = Uval.add(pe, UN);
            Tval = Tval.add(pe, TN);
        }
        this.result = new SolverResult();
        this.result.QN = Qval;
        this.result.UN = Uval;
        this.result.TN = Tval;
        Matrix RN = new Matrix(M, K); RN.zero();
        for (int i = 0; i < M; i++) {
            for (int k = 0; k < K; k++) {
                double tq = Tval.get(i, k);
                RN.set(i, k, tq > 0 ? Qval.get(i, k) / tq : 0.0);
            }
        }
        this.result.RN = RN;
    }

    // ---- state-vector helpers ----
    private static Matrix rowNormalizeNonneg(Matrix v) {
        Matrix r = asRowVec(v);
        for (int i = 0; i < r.length(); i++) {
            if (r.get(i) < 0) {
                r.set(i, 0.0);
            }
        }
        double s = r.elementSum();
        if (s > 0) {
            r.divideEq(s);
        }
        return r;
    }

    private static Matrix asRowVec(Matrix v) {
        if (v.getNumRows() == 1) {
            return v.copy();
        }
        if (v.getNumCols() == 1) {
            return v.transpose();
        }
        // General: flatten in row-major order into a 1xN row.
        int n = v.length();
        Matrix r = new Matrix(1, n);
        for (int i = 0; i < n; i++) {
            r.set(0, i, v.get(i));
        }
        return r;
    }

    private static Matrix extractRowVec(Matrix m, int row) {
        int d = m.getNumCols();
        Matrix r = new Matrix(1, d);
        for (int j = 0; j < d; j++) {
            r.set(0, j, m.get(row, j));
        }
        return r;
    }

    private static Matrix listToMatrix(List<double[]> rows) {
        int nt = rows.size();
        int d = (nt > 0) ? rows.get(0).length : 0;
        Matrix m = new Matrix(nt, d);
        for (int i = 0; i < nt; i++) {
            double[] r = rows.get(i);
            for (int j = 0; j < d; j++) {
                m.set(i, j, r[j]);
            }
        }
        return m;
    }

    // Increment weights w[i] = F(t_i) - F(t_{i-1}) of a MAP/PH CDF over the grid
    // t, with w[0] = 0 (matches the meanfield post()/finish() weighting).
    private static Matrix cdfIncrementWeights(MatrixCell map, double[] t) {
        int nt = t.length;
        Matrix w = new Matrix(nt, 1);
        w.zero();
        if (nt < 2 || map == null || map.size() < 2) {
            return w;
        }
        Matrix D0 = map.get(0);
        Matrix D1 = map.get(1);
        Matrix tLater = new Matrix(nt - 1, 1);
        Matrix tEarlier = new Matrix(nt - 1, 1);
        for (int i = 1; i < nt; i++) {
            tLater.set(i - 1, 0, t[i]);
            tEarlier.set(i - 1, 0, t[i - 1]);
        }
        Matrix cdfL = map_cdf(D0, D1, tLater);
        Matrix cdfE = map_cdf(D0, D1, tEarlier);
        for (int i = 1; i < nt; i++) {
            w.set(i, 0, cdfL.get(i - 1, 0) - cdfE.get(i - 1, 0));
        }
        return w;
    }

    // Weighted average of the rows of pit (nt x d) by w (nt x 1): (w'*pit)/sum(w).
    // Returns null when the weight mass is non-positive or contains NaN (matching
    // the MATLAB guard that falls back to the terminal distribution).
    private static Matrix weightedAvgRows(Matrix pit, Matrix w) {
        int nt = pit.getNumRows();
        int d = pit.getNumCols();
        double sw = 0.0;
        for (int i = 0; i < nt; i++) {
            double wi = w.get(i, 0);
            if (Double.isNaN(wi)) {
                return null;
            }
            sw += wi;
        }
        if (!(sw > 0)) {
            return null;
        }
        Matrix r = new Matrix(1, d);
        r.zero();
        for (int j = 0; j < d; j++) {
            double acc = 0.0;
            for (int i = 0; i < nt; i++) {
                acc += w.get(i, 0) * pit.get(i, j);
            }
            r.set(0, j, acc / sw);
        }
        return r;
    }

    protected SolverResult analyze(int it, int e) {
        // Matching MATLAB SolverENV.analyze(): simply reset solver and get transient averages.
        // Initial state is set by pre() (iteration 1) or post() (subsequent iterations)
        // via initFromMarginal() on the ensemble model.
        try {
            this.envSolvers[e].reset();
            // Clear stale init_sol from previous FLD run. In MATLAB, solver_fluid_analyzer
            // receives options by value, so init_sol modifications are local. In JAR, the FCFS
            // convergence loop modifies this.options.init_sol directly. Without clearing,
            // initSol() would be skipped and the updated model state (from initFromMarginal)
            // would not be read.
            this.envSolvers[e].options.init_sol = new Matrix(0, 0);

            // Note: timespan is NOT adjusted here - getTranAvg() uses 30/minrate
            // when timespan is unset, matching MATLAB's getTranAvg.m behavior.
            //
            // ASK FOR THE POINTS THE SOJOURN WEIGHT NEEDS, as MATLAB's analyze_
            // and the native-python twin do. refineForCdf below still resamples
            // whatever grid comes back, but reading the quadrature grid off a
            // LINEAR interpolation of the integrator's steps is a first-order
            // error its own dense output does not have: it left this engine
            // 3.7e-4 from MATLAB on renv_fourstages_repairmen's Queue1 QLen
            // while every other codebase agreed to 1e-4. Set for THIS stage
            // solve and put back, since a stage solver's own getters must not
            // inherit it.
            double[] savedTranpoints = this.envSolvers[e].options.tranpoints;
            this.envSolvers[e].options.tranpoints = stageCdfGrid(e);
            try {
                solverTranAvgOf(this.envSolvers[e]);
            } finally {
                this.envSolvers[e].options.tranpoints = savedTranpoints;
            }

            // see _kb/06-solver-catalog.md (JAR-only implementation notes: per-iteration result deep copy)
            SolverResult stageResult = this.envSolvers[e].result.deepCopy();

            // Set QN/UN/TN to entry state (first time point) for convergence checking.
            // Matching MATLAB converged(): compares Qik.metric(1) which is the t=0 value.
            int M = sn[0].nstations;
            int K = sn[0].nclasses;
            if (stageResult.QNt != null) {
                stageResult.QN = new Matrix(M, K);
                stageResult.UN = new Matrix(M, K);
                stageResult.TN = new Matrix(M, K);
                for (int i = 0; i < M; i++) {
                    for (int k = 0; k < K; k++) {
                        if (stageResult.QNt[i] != null && stageResult.QNt[i][k] != null) {
                            stageResult.QN.set(i, k, stageResult.QNt[i][k].get(0, 0));
                        }
                        if (stageResult.UNt != null && stageResult.UNt[i] != null && stageResult.UNt[i][k] != null) {
                            stageResult.UN.set(i, k, stageResult.UNt[i][k].get(0, 0));
                        }
                        if (stageResult.TNt != null && stageResult.TNt[i] != null && stageResult.TNt[i][k] != null) {
                            stageResult.TN.set(i, k, stageResult.TNt[i][k].get(0, 0));
                        }
                    }
                }
            }
            return stageResult;
        } catch (Exception ex) {
            // Skip analysis for stages where transient analysis fails
            SolverResult emptyResult = new SolverResult();
            return emptyResult;
        }
    }

    private void sanityCheck(SolverResult iterativeResult, int e) {
        int M = iterativeResult.QN.getNumRows();
        int K = iterativeResult.QN.getNumCols();
        double N = sn[0].nclosedjobs;
        for (int k = 0; k < K; k++) {
            boolean invalid = false;
            for (int m = 0; m < M; m++) {
                double u = iterativeResult.UN.get(m, k);
                if (u < 0 || u > ServerNum.get(k).get(m, e)) {
                    invalid = true;
                    break;
                }
            }
            if (invalid) {
                for (int m = 0; m < M; m++ ) {
                    double u = iterativeResult.UN.get(m, k);
                    iterativeResult.UN.set(m, k, max(min(u, 1), 0));
                }
            }

            invalid = false;
            for (int m = 0; m < M; m++) {
                    double l = iterativeResult.QN.get(m, k);
                    if (l < 0) {
                        invalid = true;
                        break;
                    }
            }
            Matrix QSameClass = iterativeResult.QN.getColumn(k);
            if (QSameClass.elementSum() > sn[0].njobs.get(0, k)) {
                invalid = true;
            }
            if (invalid) {
                for (int m = 0; m < M; m++ ) {
                    double l = iterativeResult.QN.get(m, k);
                    double Nk = sn[e].njobs.get(0, k);
                    iterativeResult.QN.set(m, k, max(min(l, Nk), 0));
                }
            }
        }


        for (int k = 0; k < K; k++) {
            double Nk = sn[0].njobs.get(0, k);
            // Only scale for closed classes (finite population)
            if (!Double.isFinite(Nk)) {
                continue;
            }
            double colSum = iterativeResult.QN.getColumn(k).elementSum();
            if (colSum > 0 && colSum != Nk) {
                double scale = Nk / colSum;
                for (int m = 0; m < M; m++) {
                    iterativeResult.QN.set(m, k, iterativeResult.QN.get(m, k) * scale);
                }
            }
        }
    }

    private Matrix getXiMin(Matrix xi, double si) {
        Matrix result = new Matrix(1, xi.getNumCols());
        for (int i = 0; i < xi.getNumCols(); i++) {
            result.set(0,  i, min(xi.get(0, i), si));
        }
        return result;
    }

    public CTMCResult runAnalyzerByCTMC() {
        line_debug(options.verbose, "ENV: running CTMC-based analyzer for random environment");
        init();
        int E = getNumberOfModels();
        int M = sn[0].nstations;
        int K = sn[0].nclasses;
        Matrix[] InfgenMatrices = new Matrix[E];
        MatrixCell stateSpace = new MatrixCell(E);
        for (int e = 0; e < E; e++) {
            Network network = (Network) envModels[e];
            SolverCTMC solverCTMC = new SolverCTMC(network, this.options);
            SolverCTMC.generatorResult generatorResult = solverCTMC.getGenerator();
            CTMCResult ctmcResult = (CTMCResult) solverCTMC.result;
            InfgenMatrices[e] = ctmcResult.infGen;
            stateSpace.set(e, ctmcResult.spaceAggr);
        }
        int states = InfgenMatrices[0].getNumRows();
        Matrix Egen = new Matrix(E0);
        Matrix Q = new Matrix(E*states, E*states);
        for (int e = 0; e < E; e++) {
           for (int h = 0; h < E; h++) {
               if (e == h) {
                   Matrix block = InfgenMatrices[e];
                   Matrix identity = Matrix.eye(states);
                   identity.scaleEq(Egen.get(e, e));
                   block.subEq(identity);
                     for (int i = 0; i < states; i++) {
                          for (int j = 0; j < states; j++) {
                            Q.set(e * states + i, h * states + j, block.get(i, j));
                          }
                     }
               }
               else {
                   // Environment transitions are state-independent: rate Egen[e,h]
                   // applied uniformly to all system states (diagonal block = Egen[e,h] * I)
                   for (int i = 0; i < states; i++) {
                       Q.set(e * states + i, h * states + i, Egen.get(e, h));
                   }
               }

           }
        }

        Matrix combinedInfGen = ctmc_makeinfgen(Q);
        Matrix pi = ctmc_solve(combinedInfGen);
        Matrix QN = new Matrix(M, K);
        Matrix UN = new Matrix(M, K);
        Matrix TN = new Matrix(M, K);

        for (int e = 0; e < E; e++) {
            for (int state = 0; state < states; state++) {
                double p = pi.get(0, e * states + state);
                for (int m = 0; m < M; m++) {
                    // Compute total jobs at station m in this state
                    double nTotal = 0;
                    for (int r = 0; r < K; r++) {
                        nTotal += stateSpace.get(e).get(state, m * K + r);
                    }
                    double c = ServerNum.get(0).get(m, e); // nservers at station m

                    for (int k = 0; k < K; k++) {
                        double prob = stateSpace.get(e).get(state, m * K + k);
                        QN.set(m, k, QN.get(m, k) + p * prob);

                        if (Double.isInfinite(c)) {
                            // Infinite-server (Delay): UN = QN, TN = n_k * mu_k
                            UN.set(m, k, UN.get(m, k) + p * prob);
                            TN.set(m, k, TN.get(m, k) + p * prob * SRates.get(k).get(m, e));
                        } else {
                            // Finite-server (PS): apply capacity sharing factor
                            double scaling = (nTotal > 0) ? min(nTotal, c) / nTotal : 0;
                            TN.set(m, k, TN.get(m, k) + p * prob * SRates.get(k).get(m, e) * scaling);
                            double uContrib = (nTotal > 0) ? prob * min(nTotal, c) / (nTotal * c) : 0;
                            UN.set(m, k, UN.get(m, k) + p * uContrib);
                        }
                    }
                }
            }
        }
        CTMCResult result = new CTMCResult();
        result.QN = QN;
        result.UN = UN;
        result.TN = TN;
        result.infGen = Q;
        return result;
    }

    public final NetworkAvgTable getAvgTable() {
        return getAvgTable(this.options);
    }

    public final NetworkAvgTable getAvgTable(SolverOptions options, boolean keepDisabled) {
        return jline.io.LineResultRecorder.around(this, "avg", () -> getAvgTableImpl(options, keepDisabled));
    }

    /**
     * Body of {@link #getAvgTable(SolverOptions, boolean)}, split out so {@link jline.io.LineResultRecorder}
     * sees what the getter RETURNED. The JAVA cross-codebase parity row is
     * measured from that rather than from what an example printed.
     */
    protected final NetworkAvgTable getAvgTableImpl(SolverOptions options, boolean keepDisabled) {

        AvgTable avgTable = this.getEnsembleAvg();
        Matrix QN = this.result.QN;
        Matrix UN = this.result.UN;
        Matrix TN = this.result.TN;

        int M = QN.getNumRows();
        int K = QN.getNumCols();

        NetworkAvgTable networkAvgTable;

        if (QN.isEmpty()) {
            throw new RuntimeException(
                    "Unable to compute results and therefore unable to print AvgTable.");
        }

        List<Double> Qval = new ArrayList<>();
        List<Double> Uval = new ArrayList<>();
        List<Double> Residval = new ArrayList<>();
        List<Double> ArvR = new ArrayList<>();
        List<Double> Tval = new ArrayList<>();
        List<Double> respTVal = new ArrayList<>();
        List<String> className = new ArrayList<>();
        List<String> stationName = new ArrayList<>();
        if (!keepDisabled) {
            for (int i = 0; i < M; i++) {
                for (int k = 0; k < K; k++) {
                    if (QN.get(i, k) + UN.get(i, k) + TN.get(i, k) > 0) {
                        Qval.add(QN.get(i, k));
                        Uval.add(UN.get(i, k));
                        Residval.add(0.0);
                        ArvR.add(0.0);
                        Tval.add(TN.get(i, k));
                        if (TN.get(i, k) == 0 && QN.get(i, k) == 0) {
                            respTVal.add(0.0);
                        }
                        else if (TN.get(i, k) == 0) {
                            respTVal.add(Inf);
                        } else {
                            respTVal.add(QN.get(i, k) / TN.get(i, k));
                        }
                        className.add(sn[0].jobclasses.get(k).getName());
                        stationName.add(sn[0].stations.get(i).getName());
                    }
                }
            }
            networkAvgTable = new NetworkAvgTable(Qval, Uval, respTVal, Residval, ArvR, Tval);
        } else {
            networkAvgTable = (NetworkAvgTable) avgTable;
            for (int i = 0; i < M; i++) {
                for (int k = 0; k < K; k++) {
                    className.add(sn[0].jobclasses.get(k).getName());
                    stationName.add(sn[0].stations.get(i).getName());
                }
            }
        }

        // ENV defaults verbose to SILENT to mute the ensemble progress log; that must not
        // also mute an explicitly requested print(), which MATLAB/Python always emit.
        SolverOptions tableOptions = options.copy();
        if (tableOptions.verbose == VerboseLevel.SILENT) {
            tableOptions.verbose = VerboseLevel.STD;
        }
        networkAvgTable.setOptions(tableOptions);
        networkAvgTable.setClassNames(className);
        networkAvgTable.setStationNames(stationName);
        return networkAvgTable;
    }

    public final NetworkAvgTable getAvgTable(SolverOptions options) {
        boolean keepDisabled = false;
        return getAvgTable(options, keepDisabled);
    }

    // Return table of average station metrics
    // TODO: this method should be based on getAvgTable and return an appropriate table class
    /**
     * Prints the average metrics table for all stations and job classes.
     *
     * @param keepDisabled If true, includes disabled/inactive entries (0-valued metrics).
     */
    public void printAvgTable(boolean keepDisabled) {

        NetworkAvgTable networkAvgTable = this.getAvgTable(this.options, keepDisabled);
        Matrix QN = this.result.QN;

        if (QN.isEmpty()) {
            throw new RuntimeException(
                    "Unable to compute results and therefore unable to print AvgTable.");
        }

        List<Double> Qval = networkAvgTable.get(0);
        List<Double> Uval = networkAvgTable.get(1);
        List<Double> Tval = networkAvgTable.get(5);
        List<Double> respTVal = networkAvgTable.get(2);
        List<String> className = networkAvgTable.getClassNames();
        List<String> stationName = networkAvgTable.getStationNames();

            System.out.printf(
                    "\n%-12s\t %-12s\t %-10s\t %-10s\t %-10s\t %-10s",
                    "Station", "JobClass", "QLen", "Util", "RespT", "Tput");
            System.out.println(
                    "\n------------------------------------------------------------------------------");
            NumberFormat nf = NumberFormat.getNumberInstance();
            nf.setMinimumFractionDigits(5);
            for (int i = 0; i < stationName.size(); i++) {
                System.out.format(
                        "%-12s\t %-12s\t %-10s\t %-10s\t %-10s\t %-10s\n",
                        stationName.get(i),
                        className.get(i),
                        nf.format(Qval.get(i)),
                        nf.format(Uval.get(i)),
                        nf.format(respTVal.get(i)),
                        nf.format(Tval.get(i)));
            }
            System.out.println(
                    "------------------------------------------------------------------------------");
    }

    public void printAvgTable() {
        boolean keepDisabled = false;
        this.printAvgTable(keepDisabled);
    }

    @Override
    public void runAnalyzer() {
        // Propagate solver verbose level to global
        if (this.options != null) {
            GlobalConstants.Verbose = options.verbose;
        }
        line_debug(options.verbose, String.format("ENV runAnalyzer: method=%s, iter_max=%d",
            options.method != null ? options.method : "default", options.iter_max));
        iterate();
    }

    /**
     * Returns the network structures for all stages in the environment.
     * Each element in the returned array corresponds to a stage's network structure.
     *
     * @return Array of NetworkStruct objects, one for each stage
     */
    public NetworkStruct[] getStruct() {
        int E = getNumberOfModels();
        NetworkStruct[] envsn = new NetworkStruct[E];
        for (int e = 0; e < E; e++) {
            envsn[e] = envStructOf(envModels[e]);
        }
        return envsn;
    }

    /**
     * Returns the list of valid solution methods supported by this solver.
     * Currently only supports the "default" method.
     *
     * @return Array of method names
     */
    public String[] listValidMethods() {
        // SolverENV.m verbatim. Each name selects a COUPLING -- what crosses an
        // environment switch -- and each is dispatched in this class or its
        // analyzers: 'default'/'meanfield' carry the marginal means, 'statevec'/
        // 'blend' the whole joint distribution, 'smp' lifts the Markovian-arc
        // check, 'statedep' makes the transition depend on the state it leaves,
        // and 'avg'/'dec' are the closed-form fast/slow environment limits. The
        // list used to name only "default", so six dispatched methods were
        // reachable only by reading the source.
        return new String[]{"default", "meanfield", "smp", "statedep", "statevec", "blend", "avg", "dec"};
    }

    private MatrixCell findBestPartition(Matrix E0) {
        double thetaMin = 0.01 * E0.elementMax();
        int E = getNumberOfModels();
        double bestEps = Inf;
        MatrixCell bestMS = null;
        for (int e = 2; e < E; e++) {
            List<int[]> rgs = new ArrayList<>();
            generatePartitionsRG(0, E, e, new int[E], -1, rgs);
            for (int[] rg : rgs) {
                boolean valid = true;
                MatrixCell MS = msFromRG(rg, e);
                for (int i = 0; i < MS.size(); i++) {
                    Matrix mi = MS.get(i);
                    Matrix subE0 = E0.getSubMatrix(mi, mi);
                    subE0.absEq();
                    if (subE0.elementMin() < thetaMin) {
                        valid = false;
                        break;
                    }
                }
                if (!valid) continue;
                Compression_result cr = ctmc_decompose(E0, MS, options);
                if (cr.eps < bestEps) {
                    bestEps = cr.eps;
                    bestMS = MS;
                    this.Ecompress = e;
                }
            }
        }
        return bestMS;
    }

    private MatrixCell msFromRG(int[] rg, int k) {
        MatrixCell MS = new MatrixCell(k);
        for (int block = 0; block < k; block++) {
            List<Integer> members = new ArrayList<>();
            for (int i = 0; i < rg.length; i++) {
                if (rg[i] == block) {
                    members.add(i);
                }
            }
            // turn that list into a int[] and wrap in a Matrix
            int[] idx = members.stream().mapToInt(Integer::intValue).toArray();
            MS.set(block, new Matrix(idx));
        }
        return MS;
    }

    private void generatePartitionsRG(int pos, int n, int k, int[] a, int currentMax, List<int[]> out) {
        if (pos == n) {
            if (currentMax == k - 1) {
                out.add(a.clone());
            }
            return;
        }
        int limit = min(currentMax + 1, k - 1);
        for (int label = 0; label <= limit; label++) {
            a[pos] = label;
            generatePartitionsRG(pos + 1, n, k, a, max(currentMax, label), out);
        }
    }

    public static Compression_result ctmc_courtois(Matrix Q, MatrixCell MS, double q) {
        Matrix v = new Matrix(Q.getNumRows(), 1);
        // fill v with the ordered macro-states
        int index = 0;
        for (int i = 0; i < MS.size(); i++) {
            Matrix macroState = MS.get(i);
            for (int j = 0; j < macroState.getNumRows(); j++) {
                v.set(index++, macroState.get(j, 0));
            }
        }
        Matrix Qperm = Q.getSubMatrix(v, v);
        Matrix Qdec = new Matrix(Qperm);
        int procRows = 0;
        for (int i = 0; i < MS.size(); i++) {
            int subSize = MS.get(i).getNumRows();
            if (procRows > 0) {
                // fill the left part of the diagonal with zeros
                for (int subrow = procRows; subrow < procRows + subSize; subrow++) {
                    for (int subcol = 0; subcol < procRows; subcol++) {
                        Qdec.set(subrow, subcol, 0.0);
                    }
                }
            }
            // fill the right part of the diagonal with zeros
            for (int subrow = procRows; subrow < procRows + subSize; subrow++) {
                for (int subcol = procRows + subSize; subcol < Q.getNumCols(); subcol++) {
                    Qdec.set(subrow, subcol, 0.0);
                }
            }
            procRows += subSize;
        }
        Qdec = ctmc_makeinfgen(Qdec);

        // Compute NCD Error Index
        Matrix epsC = Qperm.sub(Qdec);
        Matrix C = new Matrix(epsC);
        double eps = 0;

        // apply randomization coefficient
        Matrix P = new Matrix(Qperm.getNumRows(), Qperm.getNumCols());
        Qperm.divide(q, P, true);
        P.addEq(Matrix.eye(Qperm.getNumRows()));

        Matrix A = new Matrix(P);
        procRows = 0;
        for (int i = 0; i < MS.size(); i++) {
            int subSize = MS.get(i).getNumRows();
            if (procRows > 0) {
                for (int subrow = procRows; subrow < procRows + subSize; subrow++) {
                    for (int subcol = 0; subcol < procRows; subcol++) {
                        A.set(subrow, subcol, 0.0);
                    }
                }
            }
            for (int subrow = procRows; subrow < procRows + subSize; subrow++) {
                for (int subcol = procRows + subSize; subcol < Q.getNumCols(); subcol++) {
                    A.set(subrow, subcol, 0.0);
                }
            }
            procRows+= subSize;
        }
        Matrix B = P.sub(A);
        for (int i = 0; i < B.getNumRows(); i++) {
            double eleSum = B.getRow(i).elementSum();
            if (eleSum > eps) {
                eps = eleSum;
            }
        }

        // Compute epsMAX
        // We normalize A by changing the diagonal elements
        procRows = 0;
        for (int i = 0; i < MS.size(); i++) {
            int subSize = MS.get(i).getNumRows();
            for (int subrow = procRows; subrow < procRows + subSize; subrow++) {
                double diagSum = 0.0;
                for (int subcol = procRows; subcol < procRows + subSize; subcol++) {
                    if (subrow != subcol) {
                        diagSum += A.get(subrow, subcol);
                    }
                }
                A.set(subrow, subrow, 1.0 - diagSum);
            }
            procRows += subSize;
        }

        // Compute epsMAX: second-largest eigenvalue per macro-block
        procRows = 0;
        Matrix eigMS = Matrix.zeros(MS.size(), 1);
        for (int i = 0; i < MS.size(); i++) {
            int subSize = MS.get(i).getNumRows();
            int[] idxArr = new int[subSize];
            for (int j = 0; j < subSize; j++) {
                idxArr[j] = procRows + j;
            }
            Matrix idxMat = new Matrix(idxArr);
            Matrix subP = A.getSubMatrix(idxMat, idxMat);
            Ret.Eigs eigenvalue = subP.eigval();
            if (eigenvalue.values.getNumCols() > 1) {
                eigenvalue.values.absEq();
                Matrix sortedValues = eigenvalue.values.sort();
                double secondLargest = sortedValues.get( 0, sortedValues.getNumCols() - 2);
                eigMS.set(i, 0, secondLargest);
            } else {
                eigMS.set(i, 0, 0.0);
            }
            procRows += subSize;
        }
        double epsMax = (1 - eigMS.elementMax()) / 2;

        // Compute Microprobabilities
        Matrix pmicro = Matrix.zeros(Q.getNumRows(), 1);
        procRows = 0;
        for (int i = 0; i < MS.size(); i++) {
            int subSize = MS.get(i).getNumRows();
            int[] idxArr = new int[subSize];
            for (int j = 0; j < subSize; j++) {
                idxArr[j] = procRows + j;
            }
            Matrix idxMat = new Matrix(idxArr);
            Matrix subP = A.getSubMatrix(idxMat, idxMat);
            Matrix subPmicro = dtmc_solve(subP);
            for (int j = 0; j < subPmicro.getNumCols(); j++) {
                pmicro.set(procRows + j, 0, subPmicro.get(0, j));
            }
            procRows += subSize;
        }

        // Compute Macroprobabilities
        Matrix G = Matrix.zeros(MS.size(), MS.size());
        procRows = 0;
        for (int i = 0; i < MS.size(); i++) {
            int procCols = 0;
            int subSize = MS.get(i).getNumRows();
            for (int j = 0; j < MS.size(); j++) {
                if (i != j) {
                    double GammaIJ = G.get(i, j);
                    for (int iState = 0; iState < subSize; iState++) {
                        for (int jState = 0; jState < MS.get(j).getNumRows(); jState++) {
                            GammaIJ += pmicro.get(procRows + iState, 0) * P.get(procRows + iState, procCols + jState);
                        }
                    }
                    G.set(i, j, GammaIJ);
                }
                procCols += MS.get(j).getNumRows();
            }
            procRows += subSize;
        }

        // now deal with diagonal elements
        for (int i = 0; i < MS.size(); i++) {
            double GammaSum = 0.0;
            for (int j = 0; j < MS.size(); j++) {
                if (i != j) {
                    GammaSum += G.get(i, j);
                }
            }
            G.set(i, i, 1.0 - GammaSum);
        }

        // Calculate the approximate steady-state probability vector
        Matrix pMacro = dtmc_solve(G);
        Matrix p = new Matrix(Q.getNumRows(), 1);
        procRows = 0;
        for (int i = 0; i < MS.size(); i++) {
            int subSize = MS.get(i).getNumRows();
            for (int j = 0; j < subSize; j++) {
                p.set(procRows + j, 0, pMacro.get(0, i) * pmicro.get(procRows + j, 0));
            }
            procRows += subSize;
        }

        Matrix pOut = new Matrix(p.getNumRows(), 1);
        for (int i = 0; i < v.getNumRows(); i++) {
            pOut.set((int) v.get(i, 0), 0, p.get(i, 0));
        }
        p = pOut;
        return new Compression_result(p, Qperm, Qdec, eps, epsMax, P, B, C, q, pMacro, G, pmicro);
    }

    public static Compression_result ctmc_courtois(Matrix Q, MatrixCell MS) {
        Matrix v = new Matrix(Q.getNumRows(), 1);
        int index = 0;
        for (int i = 0; i < MS.size(); i++) {
            Matrix macroState = MS.get(i);
            for (int j = 0; j < macroState.getNumRows(); j++) {
                v.set(index++, macroState.get(j, 0));
            }
        }
        Matrix Qperm = Q.getSubMatrix(v, v);
        Matrix AbsQperm = new Matrix(Qperm);
        AbsQperm.absEq();
        double q = 1.05 * AbsQperm.elementMax();
        return ctmc_courtois(Q, MS, q);
    }

    /**
     * Perform CTMC decomposition using the configured method.
     * Uses options.config.da to select the decomposition algorithm:
     *   'courtois'  - Courtois decomposition (default)
     *   'kms'       - Koury-McAllister-Stewart method
     *   'takahashi' - Takahashi's method
     *   'multi'     - Multigrid method
     *
     * @param Q Infinitesimal generator matrix
     * @param MS Macro-state partition
     * @param options Solver options containing config.da and config.da_iter
     * @return Compression_result with decomposition results
     */
    public static Compression_result ctmc_decompose(Matrix Q, MatrixCell MS, SolverOptions options) {
        // Get decomposition method from options
        String method = "courtois";
        int numSteps = 10;
        if (options != null && options.config != null) {
            if (options.config.da != null) {
                method = options.config.da.toLowerCase();
            }
            numSteps = options.config.da_iter;
        }

        // Always start with Courtois to get all required fields
        Compression_result cr = ctmc_courtois(Q, MS);

        if (method.equals("courtois")) {
            return cr;
        }

        // Convert MatrixCell to List<List<Int>> for API functions
        List<List<Integer>> msList = new ArrayList<>();
        for (int i = 0; i < MS.size(); i++) {
            List<Integer> block = new ArrayList<>();
            Matrix mi = MS.get(i);
            for (int j = 0; j < mi.getNumRows(); j++) {
                block.add((int) mi.get(j, 0));
            }
            msList.add(block);
        }

        // Call the appropriate decomposition method to refine p
        switch (method) {
            case "kms": {
                jline.util.Triple<Matrix, Double, Double> result = ctmc_kms(Q, msList, numSteps);
                cr.p = result.getFirst();
                cr.eps = result.getSecond();
                cr.epsMax = result.getThird();
                break;
            }
            case "takahashi": {
                jline.util.Triple<Matrix, Double, Double> result = ctmc_takahashi(Q, msList, numSteps);
                cr.p = result.getFirst();
                cr.eps = result.getSecond();
                cr.epsMax = result.getThird();
                break;
            }
            case "multi": {
                // Multi requires MSS (macro-macro-states), default to singletons
                List<List<Integer>> mss = new ArrayList<>();
                for (int i = 0; i < MS.size(); i++) {
                    List<Integer> singleton = new ArrayList<>();
                    singleton.add(i);
                    mss.add(singleton);
                }
                jline.api.mc.Ctmc_multi.CtmcMultiResult result = ctmc_multi(Q, msList, mss);
                cr.p = result.p;
                cr.eps = result.eps;
                cr.epsMax = result.epsMAX;
                break;
            }
            default:
                throw new RuntimeException("Unknown decomposition method: " + method);
        }

        // Recompute pMacro and pmicro from the refined p (consistent with MATLAB)
        int nMacro = MS.size();
        cr.pMacro = new Matrix(1, nMacro);
        cr.pmicro = new Matrix(cr.p.getNumRows(), 1);

        for (int i = 0; i < nMacro; i++) {
            Matrix msBlock = MS.get(i);
            int blockSize = msBlock.getNumRows();

            // Compute sum of probabilities in this macro-state
            double blockSum = 0.0;
            for (int j = 0; j < blockSize; j++) {
                int idx = (int) msBlock.get(j, 0);
                blockSum += cr.p.get(idx, 0);
            }
            cr.pMacro.set(0, i, blockSum);

            // Compute normalized micro probabilities within block
            for (int j = 0; j < blockSize; j++) {
                int idx = (int) msBlock.get(j, 0);
                if (blockSum > 0) {
                    cr.pmicro.set(idx, 0, cr.p.get(idx, 0) / blockSum);
                } else {
                    cr.pmicro.set(idx, 0, 0.0);
                }
            }
        }

        return cr;
    }

    static class Compression_result {
        Matrix p;
        Matrix Qperm;
        Matrix Qdec;
        double eps;
        double epsMax;
        Matrix P;
        Matrix B;
        Matrix C;
        double q;
        Matrix pMacro;
        Matrix G;
        Matrix pmicro;

        Compression_result(Matrix p, Matrix Qperm, Matrix Qdec, double eps, double epsMax,
                           Matrix P, Matrix B, Matrix C, double q, Matrix pMacro, Matrix G, Matrix pmicro) {
            this.p = p;
            this.Qperm = Qperm;
            this.Qdec = Qdec;
            this.eps = eps;
            this.epsMax = epsMax;
            this.P = P;
            this.B = B;
            this.C = C;
            this.q = q;
            this.pMacro = pMacro;
            this.G = G;
            this.pmicro = pmicro;
        }
    }

    /**
     * The transient grid an exit average may be summed on, refined where the
     * sojourn weight cannot see the solver's own grid.
     *
     * <p>The exit metric is the Stieltjes sum {@code sum_k m(t_k) * [F(t_k) -
     * F(t_{k-1})]} evaluated on the ODE solver's OWN output grid. That grid is
     * chosen for the horizon, not for the sojourn, so a stage integrated over
     * {@code [0,1e3]} and read through an {@code Exp(1)} clock puts almost every
     * point where the weight is zero, and the answer becomes an artifact of step
     * placement rather than of the model.
     *
     * <p>The grid is rebuilt UNCONDITIONALLY -- 90% of the points under
     * {@code 5*E[S]} and the rest across the tail -- rather than only when the
     * solver's own grid looks too coarse. A "50 points inside the support is
     * enough" escape (native python's, before this) stops wherever the
     * integrator's steps happened to fall and does not converge: on
     * renv_node_breakdown the sum runs 0.462260, 0.460580, 0.459704, 0.459272,
     * 0.459138, 0.459122 as the point count goes 500 to 5e4. Rebuilding always
     * is also what makes the four codebases sum the SAME points, which is the
     * property parity needs.
     */
    private static final class CdfGrid {
        static final int NINTERP = 5000;

        final Matrix t;
        private final Matrix src;   // null when the original grid was kept

        private CdfGrid(Matrix t, Matrix src) {
            this.t = t;
            this.src = src;
        }

        /** A metric sampled on the original grid, re-sampled onto {@link #t}. */
        Matrix on(Matrix metric) {
            if (src == null || metric == null || metric.length() != src.length()) {
                return metric;
            }
            Matrix out = new Matrix(t.getNumRows(), 1);
            int n = src.length();
            int j = 0;
            for (int k = 0; k < t.getNumRows(); k++) {
                double tk = t.get(k, 0);
                while (j < n - 2 && src.get(j + 1) < tk) {
                    j++;
                }
                double x0 = src.get(j);
                double x1 = src.get(j + 1);
                double y0 = metric.get(j);
                double y1 = metric.get(j + 1);
                double v;
                if (tk <= x0) {
                    v = y0;
                } else if (tk >= x1) {
                    v = y1;
                } else {
                    v = y0 + (y1 - y0) * (tk - x0) / (x1 - x0);
                }
                out.set(k, 0, v);
            }
            return out;
        }
    }

    /**
     * The instants stage {@code e}'s exit average is summed over, as an array
     * the fluid integrator can be asked for, or null when the horizon is not
     * finite (there is then no grid to ask for).
     *
     * <p>Same construction as {@link #refineForCdf}, built from the stage
     * HORIZON rather than from a trajectory, so it can be requested before the
     * solve.</p>
     *
     * @param e stage index
     * @return increasing output instants, or null
     */
    private double[] stageCdfGrid(int e) {
        double[] ts = this.envSolvers[e].options.timespan;
        if (ts == null || ts.length < 2 || !Double.isFinite(ts[1]) || !(ts[1] > ts[0])) {
            return null;
        }
        Matrix span = new Matrix(2, 1);
        span.set(0, 0, ts[0]);
        span.set(1, 0, ts[1]);
        Matrix D0 = envObj.holdTime[e].get(0).copy();
        Matrix D1 = envObj.holdTime[e].get(1).copy();
        map_normalize(D0, D1);
        Matrix fine = refineForCdf(span, D0, D1).t;
        if (fine == null || fine.getNumRows() < 2) {
            return null;
        }
        double[] out = new double[fine.getNumRows()];
        for (int k = 0; k < out.length; k++) {
            out[k] = fine.get(k, 0);
        }
        return out;
    }

    /** The grid {@code t} refined for a sojourn distributed as the MAP (D0, D1). */
    private static CdfGrid refineForCdf(Matrix t, Matrix D0, Matrix D1) {
        int n = (t == null) ? 0 : t.getNumRows();
        if (n < 2) {
            return new CdfGrid(t, null);
        }
        double t0 = t.get(0, 0);
        double tEnd = t.get(n - 1, 0);
        double meanSojourn = map_mean(D0, D1);
        if (!(meanSojourn > 0) || Double.isInfinite(meanSojourn) || Double.isNaN(meanSojourn)) {
            meanSojourn = (tEnd - t0) / 10.0;
        }
        double tCdfEnd = Math.min(tEnd, 5.0 * meanSojourn);
        if (tCdfEnd <= t0) {
            tCdfEnd = tEnd;
        }
        int nDense = (int) (0.9 * CdfGrid.NINTERP);
        int nTail = CdfGrid.NINTERP - nDense;
        boolean withTail = tCdfEnd < tEnd && nTail > 1;
        Matrix fine = new Matrix(withTail ? nDense + nTail : nDense, 1);
        for (int k = 0; k < nDense; k++) {
            fine.set(k, 0, t0 + (tCdfEnd - t0) * k / (double) (nDense - 1));
        }
        if (withTail) {
            for (int k = 1; k <= nTail; k++) {
                fine.set(nDense + k - 1, 0, tCdfEnd + (tEnd - tCdfEnd) * k / (double) (nTail + 1 - 1));
            }
        }
        return new CdfGrid(fine, t);
    }
}
