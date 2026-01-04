/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Integration note: Enhanced with changes from fyp25-yiran-generalized-blending.git
 * integrated from commit 83efd675f3c291737f199674f3e982427b0c0212 onwards.
 */

package jline.solvers.env;

import static jline.GlobalConstants.Inf;
import static jline.io.InputOutputKt.line_warning;

import jline.io.Ret;
import jline.lang.*;
import jline.lang.constant.SolverType;
import jline.lang.nodes.StatefulNode;
import jline.lang.nodes.Station;
import jline.lang.nodes.ServiceStation;
import jline.lang.processes.Markovian;
import jline.lang.processes.ContinuousDistribution;
import jline.lang.processes.Exp;
import jline.solvers.*;
import jline.solvers.ctmc.CTMCResult;
import jline.solvers.ctmc.SolverCTMC;
import jline.solvers.des.SolverDES;
import jline.solvers.fluid.SolverFluid;
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
import static jline.api.mam.Map_cdfKt.map_cdf;
import static jline.api.mam.Map_normalizeKt.map_normalize;
import static jline.api.mam.Map_pieKt.map_pie;
import static jline.api.mc.Ctmc_makeinfgenKt.ctmc_makeinfgen;
import static jline.api.mc.Ctmc_solveKt.ctmc_solve;
import static jline.api.mc.Ctmc_kmsKt.ctmc_kms;
import static jline.api.mc.Ctmc_takahashiKt.ctmc_takahashi;
import static jline.api.mc.Ctmc_multiKt.ctmc_multi;
import static jline.api.mc.Dtmc_solveKt.dtmc_solve;
import static jline.util.Utils.isInf;

/**
 * ENV - Ensemble environment solver for models immersed in a random environment.
 */
public class SolverENV extends EnsembleSolver {

    // User-supplied representation of each stage transition
    private final Environment envObj;
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

    public SolverENV(Environment renv, NetworkSolver[] solvers) {
        super(renv, "SolverENV", new SolverOptions(SolverType.ENV));
        int E = getNumberOfModels();
        this.envObj = renv;
        this.ensemble = renv.getEnsemble().toArray(new Network[0]);
        this.solvers = solvers;
        this.sn = new NetworkStruct[E];
        this.resetFromMarginal = new Environment.ResetQueueLengthsFunction[E][E];
        this.resetEnvRates = new Environment.ResetEnvRatesFunction[E][E];
        this.result = new SolverResult();

        for (int e = 0; e < E; e++) {
            this.sn[e] = this.ensemble[e].getStruct(true);
             if (!solvers[e].supports(ensemble[e])) {
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

    public SolverENV(Environment renv, NetworkSolver[] solvers, SolverOptions options) {
        super(renv, "SolverENV", options);

        // Enable SMP method if specified in options
        if (options != null && options.method != null && options.method.equalsIgnoreCase("smp")) {
            this.SMPMethod = true;
            if (options.verbose == VerboseLevel.DEBUG) {
                System.out.println("ENV solver: SMP method enabled via options.method='smp'");
            }
        }

        int E = getNumberOfModels();
        this.envObj = renv;
        this.ensemble = renv.getEnsemble().toArray(new Network[0]);
        this.solvers = solvers;
        this.sn = new NetworkStruct[E];
        this.resetFromMarginal = new Environment.ResetQueueLengthsFunction[E][E];
        this.resetEnvRates = new Environment.ResetEnvRatesFunction[E][E];
        this.result = new SolverResult();

        for (int e = 0; e < E; e++) {
            this.sn[e] = this.ensemble[e].getStruct(true);
             if (!solvers[e].supports(ensemble[e])) {
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
                QEntry.setColumn(e, results.get(it - 1).get(e).QN.getColumn(k));
                QExit.setColumn(e, results.get(it).get(e).QN.getColumn(k));
            }
            double tol = options.iter_tol;
            double maxDiff = Matrix.maxAbsDiff(QEntry, QExit);
            if (maxDiff >= tol) {
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
    protected void init() {
        envObj.init();

        int M = sn[0].nstations;
        int K = sn[0].nclasses;
        int E = getNumberOfModels();
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
        this.pi = envObj.probEnv;

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
        this.envObj.probOrig = embweight;


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
                Network CompressedNetwork = ensemble[i].copy();
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
                macroSolvers[i]  = new SolverFluid(CompressedNetwork, solvers[i].options);
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


    // Solves model in stage e
//  @Override
//  protected SolverResult analyze(int it, int e) {
//    // TODO: [Qt,Ut,Tt] = self.ensemble{e}.getTranHandles;
//    this.solvers[e].reset();
//    // If pStar values exist, implement p-norm smoothing
//    if (this.solvers[e].options.config.pstar.size() != 0) {
//      PStarSearcher searcher = new PStarSearcher();

    // /      Matrix targetQueueLengths = searcher.generateTargetQueueLengths(this.solvers[e].model);
    // /      PointValuePair pStarValues =
    // /          searcher.findPStarValues(this.solvers[e].model, targetQueueLengths);
    // /      solvers[e].options.config.pstar.clear();
    // /      for (int i = 0; i < this.solvers[e].model.getNumberOfNodes(); i++) {
    // /        solvers[e].options.config.pstar.add(i, pStarValues.getPoint()[i]);
    // /      }
    // /    }
    // /    this.solvers[e].getTranAvg();
//    return this.solvers[e].result;
//  }
    @Override
    protected void pre(int it) {

        int E = getNumberOfModels();

        if (it == 1) {
            for (int e = 0; e < E; e++) {
                if (isInf(this.solvers[e].options.timespan[1])) {
                    this.solvers[e].getAvg();
                } else {
                    this.solvers[e].getTranAvg();
                }
                this.ensemble[e].initFromMarginal(this.solvers[e].result.QN);
            }
        }
    }

    // Solves model in stage e
//    @Override
//    protected SolverResult analyze(int it, int e) {
//        // TODO: [Qt,Ut,Tt] = self.ensemble{e}.getTranHandles;
//        this.solvers[e].reset();
//        this.solvers[e].getTranAvg();
//        return this.solvers[e].result;
//    }

    @Override
    protected void post(int it) {

        int M = sn[0].nstations;
        int K = sn[0].nclasses;
        int E = getNumberOfModels();
        Matrix[][] QExit = new Matrix[E][E];
        Matrix[][] UExit = new Matrix[E][E];
        Matrix[][] TExit = new Matrix[E][E];
        Matrix[][] w = new Matrix[E][E];
        Matrix[] QEntry = new Matrix[E]; // Average entry queue-length

        for (int e = 0; e < E; e++) {
            for (int h = 0; h < E; h++) {
                QExit[e][h] = new Matrix(M, K);
                UExit[e][h] = new Matrix(M, K);
                TExit[e][h] = new Matrix(M, K);
                Matrix D0 = envObj.proc[e][h].get(0);
                Matrix D1 = envObj.proc[e][h].get(1);
                map_normalize(D0, D1);
                for (int i = 0; i < M; i++) {
                    for (int r = 0; r < K; r++) {
                        w[e][h] = new Matrix(1, 1);
                        Matrix cdf1 =
                                map_cdf(
                                        D0, D1,
                                        Matrix.extractRows(
                                                results.get(it).get(e).t, 1, results.get(it).get(e).t.getNumRows(), null));

                        Matrix cdf2 =
                                map_cdf(
                                        D0, D1,
                                        Matrix.extractRows(
                                                results.get(it).get(e).t,
                                                0,
                                                results.get(it).get(e).t.getNumRows() - 1,
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
                                    results.get(it).get(e).QNt[i][r].transpose().mult(w[e][h], null).value()
                                            / w[e][h].elementSum());
                            UExit[e][h].set(
                                    i,
                                    r,
                                    results.get(it).get(e).UNt[i][r].transpose().mult(w[e][h], null).value()
                                            / w[e][h].elementSum());
                            TExit[e][h].set(
                                    i,
                                    r,
                                    results.get(it).get(e).TNt[i][r].transpose().mult(w[e][h], null).value()
                                            / w[e][h].elementSum());
                        }
                    }
                }
            }
        }

        for (int e = 0; e < E; e++) {
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
            solvers[e].reset();

            ensemble[e].initFromMarginal(QEntry[e]);
        }

        // Update transition rates between stages if State Dependent
        // Auto-detected or manually specified via options.method='statedep'
        if ("statedep".equalsIgnoreCase(stateDepMethod) || Objects.equals(options.method, "statedep")) {
            for (int e = 0; e < E; e++) {
                for (int h = 0; h < E; h++) {
                    if (envObj.env[e][h] != null && envObj.env[e][h] instanceof Markovian) {
                        // If not defined, rates are left unchanged
                        if (resetEnvRates[e][h] != null) {
                            envObj.env[e][h] =
                                    resetEnvRates[e][h].reset(
                                            (Markovian) envObj.env[e][h], QExit[e][h], UExit[e][h], TExit[e][h]);
                        }
                    }
                }
            }
            // Reinitialise
            envObj.init();
        }
    }

    @Override
    protected void finish() {

        // Use last iteration
        int it = results.size();
        int M = sn[0].nstations;
        int K = sn[0].nclasses;
        int E = getNumberOfModels();

        this.result.UN = new Matrix(M, K);
        this.result.XN = new Matrix(1, K);
        this.result.TN = new Matrix(M, K);
        this.result.QN = new Matrix(M, K);

        for (int e = 0; e < E; e++) {
            double p = pi.get(0, e);
            SolverResult resE = results.get(it).get(e);
            for (int i = 0; i < M; i++) {
                for (int k = 0; k < K; k++) {
                    this.result.QN.set(i, k, this.result.QN.get(i, k) + resE.QN.get(i, k) * p);
                    this.result.UN.set(i, k, this.result.UN.get(i, k) + resE.UN.get(i, k) * p);
                    this.result.TN.set(i, k, this.result.TN.get(i, k) + resE.TN.get(i, k) * p);
                }
            }
        }

        this.result.XN = new Matrix(1, K);
        for (int k = 0; k < K; k++) {
            this.result.XN.set(0, k, this.result.TN.get(ref, k));
        }
        result.runtime = (System.nanoTime() - startTime) / 1000000000.0;
        if (options.verbose != VerboseLevel.SILENT) {
            System.out.printf("blending completed in %d iterations in %.3f seconds%n", it, result.runtime);
        }
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
            if (solvers[e] instanceof SolverCTMC) {
                stageInfGen[e] = ((SolverCTMC) solvers[e]).getGenerator().infGen;
                stageEventFilt[e] = ((SolverCTMC) solvers[e]).getGenerator().eventFilt;
                stageEvents[e] = ((SolverCTMC) solvers[e]).getGenerator().ev;
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

                    for (int i = 0; i < ensemble[e].getNumberOfNodes(); i++) {
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
            blending();
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
            if (njobs > 0) {
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
            ensemble[e].initFromMarginal(Q_current);

            // Set solver timespan and run transient analysis
            solvers[e].options.timespan = new double[]{0, duration};
            solvers[e].reset();

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
            if (solvers[e] instanceof SolverFluid) {
                SolverFluid solverFluid = (SolverFluid) solvers[e];
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
            } else if (solvers[e] instanceof SolverDES) {
                SolverDES solverDES = (SolverDES) solvers[e];
                solverDES.options.init_sol = initial;
                // Sync sn.state from init_sol so SolverDES sees the real initial queue lengths
                NetworkStruct sn_des = solverDES.model.getStruct(false);
                List<StatefulNode> statefulNodes = sn_des.stateful;
                int offset = 0;
                for (StatefulNode node : statefulNodes) {
                    int cols = sn_des.nclasses;
                    Matrix nodeState = Matrix.extract(initial, 0, 1, offset, offset + cols);
                    sn_des.state.put(node, nodeState);
                    int nodeIdx = node.getNodeIndex();
                    int stationIdx = (int) sn_des.nodeToStation.get(0, nodeIdx);
                    solverDES.model.getStations().get(stationIdx).setState(nodeState);
                    offset += cols;
                }
                solverDES.sn = sn_des;
                stageResult = solverDES.runMethodSpecificAnalyzer();
            } else {
                throw new UnsupportedOperationException(
                    "getSamplePathTable requires SolverFluid or SolverDES as stage solver, got: " +
                    solvers[e].getClass().getSimpleName());
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
        init();
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

            if (stateDepMethod == null || stateDepMethod.isEmpty() || stateDepMethod.equalsIgnoreCase("stateindep")) {

            } else if (stateDepMethod.equalsIgnoreCase("statedep")) {
                if (!SMPMethod) {
                    for (int a = 0; a < E; a++) {
                        for (int b = 0; b < E; b++) {
                            double old_e = E0.get(a, b);
                            Eutil.set(a, b, old_e * results.get(it).get(a).UN.get(ref, 0));
                        }
                    }

                    Eutil = ctmc_makeinfgen(Eutil);
                    pi = ctmc_solve(Eutil);
                    envObj.probEnv = pi;

                    Matrix embweight = new Matrix(E, E);
                    double[] piArray = pi.toArray1D();
                    for (int i = 0; i < E; i++) {
                        double sum = 0.0;
                        for (int h = 0; h < E; h++) {
                            if (h != i) {
                                sum += piArray[h] * Eutil.get(h, i);
                            }
                        }
                        for (int k = 0; k < E; k++) {
                            if (k == i) {
                                embweight.set(k, i, 0);
                            } else {
                                embweight.set(k, i, piArray[k] * Eutil.get(k, i) / sum);
                            }
                        }
                    }
                    this.envObj.probOrig = embweight;

                    for (int k = 0; k < E; k++) {
                        Function<Double, Double> originalCdf = sojournCdfs[k];
                        double uRef = results.get(it).get(k).UN.get(ref, 0);
                        sojournCdfsUtil[k] = t -> originalCdf.apply(uRef * t);
                    }
                }

                else {
                    // New method using spline interpolation for gke
                    Map<Integer, PolynomialSplineFunction> gkeFunctions = new java.util.HashMap<>();
                    for (int k = 0; k < E; k++) {
                        if (!gkeFunctions.containsKey(k)) {
                            Matrix UNt = UNtStages.get(k)[ref][0];
                            Matrix tVector = tStages.get(k);
                            double[] gkeTime = new double[tVector.length()];
                            double[] gkeValues = new double[tVector.length()];
                            double integral = 0.0;
                            for (int i = 0; i < tVector.length(); i++) {
                                gkeTime[i] = tVector.get(i, 0);
                                double dt = (i == 0) ? gkeTime[0] : gkeTime[i] - gkeTime[i - 1];
                                integral += UNt.get(i, 0) * dt;
                                gkeValues[i] = gkeTime[i] > 0 ? integral / gkeTime[i] : 0.0;
                            }
                            PolynomialSplineFunction gkeFunc =
                                    new SplineInterpolator().interpolate(gkeTime, gkeValues);
                            gkeFunctions.put(k, gkeFunc);
                        }
                        for (int h = 0; h < E; h++) {
                            if (envObj.env[k][h] != null) {
                                // Use precomputed spline for transitionCdfs
                                final int kk = k;
                                final int hh = h;
                                transitionCdfs[k][h] = t -> {
                                    PolynomialSplineFunction gkeFunc = gkeFunctions.get(kk);
                                    double tEval = t;
                                    double minT = gkeFunc.getKnots()[0];
                                    double maxT = gkeFunc.getKnots()[gkeFunc.getN()];
                                    if (tEval < minT) tEval = minT;
                                    if (tEval > maxT) tEval = maxT;
                                    double gke = gkeFunc.value(tEval);
                                    return envObj.env[kk][hh].evalCDF(gke * t);
                                };
                            }
                        }
                    }

                    // update sojourn cdfs
                    Function<Double, Double>[] sojournCdfsUtil = new Function[E];
                    for (int k = 0; k < E; k++) {
                        final int kk = k;
                        PolynomialSplineFunction gkeFunc = gkeFunctions.get(kk);
                        sojournCdfsUtil[kk] = (Double t) -> {
                            double minT = gkeFunc.getKnots()[0];
                            double maxT = gkeFunc.getKnots()[gkeFunc.getN()];
                            double tEval = t;
                            if (tEval < minT) tEval = minT;
                            if (tEval > maxT) tEval = maxT;
                            double gke = gkeFunc.value(tEval);

                            double surv = 1.0;
                            for (int j = 0; j < E; j++) {
                                if (j == kk) continue;
                                if (envObj.env[kk][j] == null) continue;
                                double Fkj = envObj.env[kk][j].evalCDF(gke * t);
                                surv *= (1.0 - Fkj);
                            }
                            return 1.0 - surv;
                        };
                    }
                    this.sojournCdfsUtil = sojournCdfsUtil;

                    // update firing rates
                    double[][] tVectors = new double[E][];
                    for (int k = 0; k < E; k++) {
                        Matrix tMat = tStages.get(k);
                        tVectors[k] = tMat.toArray1D();
                    }
                    for (int k = 0; k < E; k++) {
                        for (int e = 0; e < E; e++) {
                            int kk = k;
                            int ee = e;
                            if (transitionCdfs[kk][ee] == null) {
                                Eutil.set(kk, ee, 0.0);
                                continue;
                            }

                            double[] tVector = tVectors[kk];

                            UnivariateFunction transRate = t -> {
                                double val = transitionCdfs[kk][ee].apply(t);
                                if (Double.isNaN(val) || val > 1.0 || val < 0.0) return 0.0;
                                return 1.0 - val;
                            };

                            double integral = 0.0;
                            for (int i = 1; i < tVector.length; i++) {
                                double t1 = tVector[i - 1];
                                double t2 = tVector[i];
                                double y1 = transRate.value(t1);
                                double y2 = transRate.value(t2);
                                integral += 0.5 * (t2 - t1) * (y1 + y2);
                            }

                            if (Double.isFinite(integral) && integral > 1e-8) {
                                Eutil.set(kk, ee, 1.0 / integral);
                            } else {
                                line_warning("SolverENV", "Invalid or zero sojourn time for transition (%d,%d): %f",
                                        kk, ee, integral);
                                Eutil.set(kk, ee, 0.0);
                            }
                        }
                    }
                    Eutil = ctmc_makeinfgen(Eutil);

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
                                        if (h != k && h != e && transitionCdfs[k][h] != null) {
                                            // Use midpoint for better accuracy in survival probability calculation
                                            double tmid = (t0 + t1) / 2.0;
                                            survival *= (1.0 - transitionCdfs[k][h].apply(tmid));
                                        }
                                    }
                                    sum += deltaF * survival;
                                }
                                dtmcP.set(k, e, sum);
                            }
                        }
                    }
                    Matrix dtmcPie = dtmc_solve(dtmcP);
                    getHoldTime(E, tStages);
                    // update pi with hold times
                    for (int k = 0; k < E; k++) {
                        double sum = 0;
                        for (int e = 0; e < E; e++) {
                            sum += dtmcPie.get(e) * holdTime.get(0, e);
                        }
                        this.pi.set(0, k, dtmcPie.get(k) * holdTime.get(0, k) / sum);
                    }
                    this.envObj.probEnv = pi;

                    // update embweight
                    Matrix newEmbweight = new Matrix(E, E);
                    double[] newPiArray = envObj.probEnv.toArray1D();
                    for (int e = 0; e < E; e++) {
                        double sum = 0.0;
                        for (int h = 0; h < E; h++) {
                            if (h != e) {
                                sum += newPiArray[h] * Eutil.get(h, e);
                            }
                        }
                        for (int k = 0; k < E; k++) {
                            if (k == e) {
                                newEmbweight.set(k, e, 0);
                            } else {
                                newEmbweight.set(k, e, newPiArray[k] * Eutil.get(k, e) / sum);
                            }
                        }
                    }
                    this.envObj.probOrig = newEmbweight;
                }
            } else {
                throw new IllegalArgumentException("Unknown state-dependent method: " + this.stateDepMethod);
            }


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

    protected SolverResult analyze(int it, int e) {
        int M = sn[0].nstations;
        int E = getNumberOfModels();
        double epsilon = 1E-30;
        double N = sn[0].nclosedjobs;
        int K = sn[0].nclasses;
        SolverResult iterativeResult = new SolverResult();
        iterativeResult.reset();
        if (it == 1) {
            iterativeResult.QN = new Matrix(M, K);
            for (int k = 0; k < K; k++) {
                double Nk = sn[0].njobs.get(0, k);
                double initVal = Nk / (double) M;
                for (int m = 0; m < M; m++) {
                    iterativeResult.QN.set(m, k, initVal);
                }
            }
        } else {
            iterativeResult.QN = new Matrix(M, K);
        }
        iterativeResult.UN = Matrix.ones(M, K);
        iterativeResult.TN = Matrix.ones(M, K);
        for (int j = 0; j < E; j++) {
            double weight = envObj.probOrig.get(j, e);
            if (it > 1) {
                iterativeResult.QN = iterativeResult.QN.add(weight, results.get(it - 1).get(j).QN);
            }
        }

        // solve the ode
        Function<Double, Double> sojournCDF = sojournCdfsUtil[e];
        double upperLimit = 10;
        double utilRef;
        if (it == 1) {
            utilRef = 1;
        }
        else {
            utilRef = results.get(it-1).get(e).UN.get(ref,0);
        }

        double convergeTol;
        if (abs(utilRef) > epsilon) {
            convergeTol = max(epsilon * abs(utilRef), epsilon);
        }
        else {
            convergeTol = epsilon;
            line_warning("SolverENV", "utilRef is very small (%f), using absolute tolerance for upper limit calculation", utilRef);
        }

        int maxIteration = 1000;
        int iterCount = 0;
        while (iterCount < maxIteration) {
            double residual = abs(sojournCDF.apply(upperLimit) - 1);
            if (residual <= convergeTol) {
                break;
            }
            upperLimit *= 1.2;
            iterCount++;
            // Prevent upperLimit from growing too large
            if (upperLimit > 1e8) {
                break;
            }
        }
        if (iterCount >= maxIteration) {
            line_warning("SolverENV", "Upper limit calculation did not converge after %d iterations. Final upperLimit=%f",
                    maxIteration, upperLimit);
        }

        Matrix initialState = iterativeResult.QN.copy();

        Solver solver = solvers[e];
        Matrix initial = new Matrix(1, initialState.getNumElements());
        int index = 0;
        for (int m = 0; m < M; m++) {
            for (int k = 0; k < K; k++) {
                initial.set(index++, initialState.get(m, k));
            }
        }
        if (solver instanceof SolverFluid) {
            SolverFluid solverFluid = (SolverFluid) solver;
            solverFluid.options.timespan[1] = upperLimit;
            solverFluid.options.init_sol = initial;
            // sync sn.state from init_sol so SolverFluid sees the real initial queue lengths
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
            SolverResult stage_result = solverFluid.runMethodSpecificAnalyzer();

            iterativeResult.QNt = stage_result.QNt;
            iterativeResult.UNt = stage_result.UNt;
            iterativeResult.TNt = stage_result.TNt;
            iterativeResult.t = stage_result.t;
            iterativeResult.iter = it;
            UNtStages.add(stage_result.UNt);
            tStages.set(e, stage_result.t);
        } else if (solver instanceof SolverDES) {
            SolverDES solverDES = (SolverDES) solver;
            solverDES.options.timespan[1] = upperLimit;
            solverDES.options.init_sol = initial;
            // sync sn.state from init_sol so SolverDES sees the real initial queue lengths
            NetworkStruct sn_des = solverDES.model.getStruct(false);
            List<StatefulNode> statefulNodes = sn_des.stateful;
            int offset = 0;
            for (StatefulNode node : statefulNodes) {
                int cols = sn_des.nclasses;
                Matrix nodeState = Matrix.extract(initial, 0, 1, offset, offset + cols);
                sn_des.state.put(node, nodeState);
                int nodeIdx = node.getNodeIndex();
                int stationIdx = (int) sn_des.nodeToStation.get(0, nodeIdx);
                solverDES.model.getStations().get(stationIdx).setState(nodeState);
                offset += cols;
            }
            solverDES.sn = sn_des;
            SolverResult stage_result = solverDES.runMethodSpecificAnalyzer();

            iterativeResult.QNt = stage_result.QNt;
            iterativeResult.UNt = stage_result.UNt;
            iterativeResult.TNt = stage_result.TNt;
            iterativeResult.t = stage_result.t;
            iterativeResult.iter = it;
            UNtStages.add(stage_result.UNt);
            tStages.set(e, stage_result.t);
        }
        int tR = iterativeResult.t.getNumRows();

        Matrix weight = new Matrix(tR, 1);
        double prevCdf = 0.0;
        for (int i = 0; i < tR - 1; i++) {
            double ti   = iterativeResult.t.get(i, 0);
            double cdfI = sojournCDF.apply(ti);
            weight.set(i, 0, cdfI - prevCdf);
            prevCdf = cdfI;
        }
        weight.set(tR - 1, 0, 1.0 - prevCdf);


        for (int i = 0; i < M; i++) {
            for (int k = 0; k < K; k++) {
                iterativeResult.QN.set(i, k, iterativeResult.QNt[i][k].transpose().mult(weight, null).value());
                if (isInf(ServerNum.get(k).get(i, e))) {
                    iterativeResult.UN.set(i, k, iterativeResult.QNt[i][k].transpose().mult(weight, null).value());
                } else {
                    Matrix xi = iterativeResult.QNt[i][k].transpose();
                    double si = ServerNum.get(k).get(i, e);
                    Matrix xiMin = getXiMin(xi, si);
                    iterativeResult.UN.set(i, k, xiMin.mult(weight, null).value());
                }
                Matrix xi = iterativeResult.QNt[i][k].transpose();
                double si = ServerNum.get(k).get(i, e);
                Matrix xiMin = getXiMin(xi, si);
                iterativeResult.TN.set(i, k, SRates.get(k).get(i, e) * xiMin.mult(weight, null).value());
            }
        }
        sanityCheck(iterativeResult, e);
        return iterativeResult;
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
        init();
        int E = getNumberOfModels();
        int M = sn[0].nstations;
        int K = sn[0].nclasses;
        // Dynamically find which station the environment transition is based on
        MatrixCell EnvAtStation = new MatrixCell(K);
        for (int k = 0; k < K; k++) {
            Matrix EnvAtStationK = new Matrix(E, E);
            for (int e = 0; e < E; e++) {
                for (int h = 0; h < E; h++) {
                    boolean found = false;
                    for (int m = 0; m < M; m++) {
                        if (SRates.get(k).get(m, e) != SRates.get(k).get(m, h)) {
                            EnvAtStationK.set(e, h, m);
                            found = true;
                            break;
                        }
                    }
                    if (!found) {
                        EnvAtStationK.set(e, h, -1);
                    }
                }
            }
            EnvAtStation.set(k, EnvAtStationK);
        }

        Matrix[] InfgenMatrices = new Matrix[E];
        MatrixCell stateSpace = new MatrixCell(E);
        for (int e = 0; e < E; e++) {
            Network network = ensemble[e];
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
                   int classIndex = -1;
                   for (int k = 0; k < K; k++) {
                       boolean difference = false;
                       for (int m = 0; m < M; m++) {
                           if (SRates.get(k).get(m, e) != SRates.get(k).get(m, h)) {
                               difference = true;
                               break;
                           }
                       }
                       if (difference) {
                           classIndex = k;
                           break;
                       }
                   }
                   int stationIndex = (int) EnvAtStation.get(classIndex).get(e, h);
                     if (stationIndex < 0) {
                          stationIndex = 0;
                     }

                   Matrix gate = new Matrix(states, states);
                   Matrix occupancy = stateSpace.get(e);
                   for (int s = 0; s < states; s++) {
                      double qlen = occupancy.get(s, stationIndex * K + classIndex);
                      if (qlen >= 1) gate.set(s, s, Egen.get(e, h));
                   }
                   for (int i = 0; i < states; i++) {
                       for (int j = 0; j < states; j++) {
                           Q.set(e * states + i, h * states + j, gate.get(i, j));
                       }
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
                    for (int k = 0; k < K; k++) {
                        double prob = stateSpace.get(e).get(state, m * K + k);
                        QN.set(m, k, QN.get(m, k) + p * prob);
                        int u = (int) min(prob, ServerNum.get(k).get(m, e));
                        UN.set(m, k, UN.get(m, k) + p * u);
                        TN.set(m, k, TN.get(m, k) + p * prob * SRates.get(k).get(m, e));
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

        networkAvgTable.setOptions(options);
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
            envsn[e] = ensemble[e].getStruct(true);
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
        return new String[]{"default"};
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
                kotlin.Triple<Matrix, Double, Double> result = ctmc_kms(Q, msList, numSteps);
                cr.p = result.getFirst();
                cr.eps = result.getSecond();
                cr.epsMax = result.getThird();
                break;
            }
            case "takahashi": {
                kotlin.Triple<Matrix, Double, Double> result = ctmc_takahashi(Q, msList, numSteps);
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
                kotlin.Triple<Matrix, Double, Double> result = ctmc_multi(Q, msList, mss);
                cr.p = result.getFirst();
                cr.eps = result.getSecond();
                cr.epsMax = result.getThird();
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
}
