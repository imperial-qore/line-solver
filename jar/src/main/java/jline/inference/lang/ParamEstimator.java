/*
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 */
package jline.inference.lang;

import jline.inference.api.Infer_compute_ql_at_arrival;
import jline.inference.api.Infer_fmlps;
import jline.inference.api.Infer_gibbs;
import jline.inference.api.Infer_mlps;
import jline.inference.api.Infer_qmle;
import jline.inference.api.Sn_set_service_coc;
import jline.inference.util.NnlsSolver;
import jline.inference.util.OptimUtils;
import jline.inference.util.SplineInterpolator;
import jline.lang.ClosedClass;
import jline.lang.JobClass;
import jline.lang.Network;
import jline.lang.NetworkStruct;
import jline.lang.RoutingMatrix;
import jline.lang.constant.MetricType;
import jline.lang.constant.EventType;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Node;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.nodes.Station;
import jline.lang.processes.Det;
import jline.lang.processes.Distribution;
import jline.lang.processes.Exp;
import jline.lang.processes.Markovian;
import jline.solvers.NetworkSolver;
import jline.solvers.SolverResult;
import jline.solvers.mva.SolverMVA;
import jline.solvers.mva.analyzers.Solver_mva_analyzer;
import jline.util.Pair;
import jline.util.matrix.Matrix;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.TreeSet;
import java.util.function.Function;

/**
 * Service demand estimator for queueing network models.
 */
public class ParamEstimator {

    public final Network model;
    public final EstimatorOptions options;

    public final List<List<List<SampledMetric>>> samples;
    public final List<List<SampledMetric>> samplesAggr;

    public ParamEstimator(Network model) {
        this(model, EstimatorOptions.defaultOptions());
    }

    public ParamEstimator(Network model, EstimatorOptions options) {
        this.model = model;
        this.options = options;
        int nNodes = model.getNumberOfNodes();
        int nClasses = model.getNumberOfClasses();
        this.samples = new ArrayList<List<List<SampledMetric>>>(nNodes);
        for (int i = 0; i < nNodes; i++) {
            List<List<SampledMetric>> row = new ArrayList<List<SampledMetric>>(nClasses);
            for (int r = 0; r < nClasses; r++) {
                row.add(new ArrayList<SampledMetric>());
            }
            this.samples.add(row);
        }
        this.samplesAggr = new ArrayList<List<SampledMetric>>(nNodes);
        for (int i = 0; i < nNodes; i++) {
            this.samplesAggr.add(new ArrayList<SampledMetric>());
        }
    }

    public void addSamples(SampledMetric sampleData) {
        int i = model.getNodeIndex(sampleData.node);
        if (sampleData.isAggregate()) {
            samplesAggr.get(i).add(sampleData);
        } else {
            int r = model.getClassIndex(sampleData.jobClass);
            samples.get(i).get(r).add(sampleData);
        }
    }

    public List<List<List<SampledMetric>>> getData() { return samples; }
    public List<List<SampledMetric>> getDataAggr() { return samplesAggr; }

    public SampledMetric getArvR(Node node, JobClass jobClass) {
        int i = model.getNodeIndex(node);
        int r = model.getClassIndex(jobClass);
        for (SampledMetric sm : samples.get(i).get(r)) {
            if (sm.type == MetricType.ArvR) return sm;
        }
        return null;
    }

    public SampledMetric getUtil(Node node, JobClass jobClass) {
        int i = model.getNodeIndex(node);
        int r = model.getClassIndex(jobClass);
        for (SampledMetric sm : samples.get(i).get(r)) {
            if (sm.type == MetricType.Util) return sm;
        }
        return null;
    }

    public SampledMetric getRespT(Node node, JobClass jobClass) {
        int i = model.getNodeIndex(node);
        int r = model.getClassIndex(jobClass);
        for (SampledMetric sm : samples.get(i).get(r)) {
            if (sm.type == MetricType.RespT) return sm;
        }
        return null;
    }

    public SampledMetric getAggrUtil(Node node) {
        int i = model.getNodeIndex(node);
        for (SampledMetric sm : samplesAggr.get(i)) {
            if (sm.type == MetricType.Util) return sm;
        }
        return null;
    }

    public List<SampledMetric> getQLen(Node node, JobClass jobClass) {
        return getQLen(node, jobClass, null);
    }

    public List<SampledMetric> getQLen(Node node, JobClass jobClass, ConditionEvent ev) {
        int i = model.getNodeIndex(node);
        int r = model.getClassIndex(jobClass);
        List<SampledMetric> nodeData = samples.get(i).get(r);
        List<SampledMetric> result = new ArrayList<SampledMetric>();
        if (ev == null) {
            for (SampledMetric sm : nodeData) {
                if (sm.type == MetricType.QLen) result.add(sm);
            }
        } else {
            for (SampledMetric sm : nodeData) {
                if (sm.type == MetricType.QLen && sm.cond != null
                        && sm.cond.node == ev.node && sm.cond.jobClass == ev.jobClass
                        && sm.cond.event == ev.event) {
                    result.add(sm);
                }
            }
        }
        return result;
    }

    public SampledMetric getAggrQLen(Node node) { return getAggrQLen(node, null); }

    public SampledMetric getAggrQLen(Node node, ConditionEvent ev) {
        int i = model.getNodeIndex(node);
        List<SampledMetric> nodeData = samplesAggr.get(i);
        if (ev == null) {
            for (SampledMetric sm : nodeData) {
                if (sm.type == MetricType.QLen) return sm;
            }
        } else {
            for (SampledMetric sm : nodeData) {
                if (sm.type == MetricType.QLen && sm.cond != null
                        && sm.cond.node == ev.node && sm.cond.jobClass == ev.jobClass
                        && sm.cond.event == ev.event) {
                    return sm;
                }
            }
        }
        return null;
    }

    public SampledMetric getTput(Node node, JobClass jobClass) {
        int i = model.getNodeIndex(node);
        int r = model.getClassIndex(jobClass);
        for (SampledMetric sm : samples.get(i).get(r)) {
            if (sm.type == MetricType.Tput) return sm;
        }
        return null;
    }

    public String autoMethod() {
        boolean hasArvR = false, hasRespT = false, hasUtil = false;
        boolean hasQLen = false, hasTput = false, hasTrace = false;
        boolean hasAggrUtil = false, hasAggrQLen = false;

        for (int i = 0; i < samples.size(); i++) {
            for (int r = 0; r < samples.get(i).size(); r++) {
                for (SampledMetric sm : samples.get(i).get(r)) {
                    if (sm.type == MetricType.ArvR) hasArvR = true;
                    else if (sm.type == MetricType.RespT) hasRespT = true;
                    else if (sm.type == MetricType.Util) hasUtil = true;
                    else if (sm.type == MetricType.QLen) hasQLen = true;
                    else if (sm.type == MetricType.Tput) hasTput = true;
                    if (sm.isTrace()) hasTrace = true;
                }
            }
        }

        for (int i = 0; i < samplesAggr.size(); i++) {
            for (SampledMetric sm : samplesAggr.get(i)) {
                if (sm.type == MetricType.Util) hasAggrUtil = true;
                else if (sm.type == MetricType.QLen) hasAggrQLen = true;
            }
        }

        String method;
        if (hasTrace && hasRespT && hasAggrQLen) method = "erps";
        else if (hasTrace && hasRespT && hasArvR) method = "mlps";
        else if (hasArvR && hasRespT && (hasUtil || hasAggrUtil)) method = "ubo";
        else if (hasArvR && (hasUtil || hasAggrUtil)) method = "ubr";
        else if (hasQLen) method = "qmle";
        else throw new IllegalStateException(
                "Insufficient data to automatically select an estimation method. Please set options.method manually.");

        options.method = method;
        return method;
    }

    public void interpolate() {
        TreeSet<Double> tUnion = new TreeSet<Double>();

        for (int i = 0; i < samples.size(); i++) {
            for (int r = 0; r < samples.get(i).size(); r++) {
                for (SampledMetric sm : samples.get(i).get(r)) {
                    for (double t : sm.t) tUnion.add(t);
                }
            }
        }
        for (int i = 0; i < samplesAggr.size(); i++) {
            for (SampledMetric sm : samplesAggr.get(i)) {
                for (double t : sm.t) tUnion.add(t);
            }
        }

        double[] tNew = new double[tUnion.size()];
        int idx = 0;
        for (Double t : tUnion) tNew[idx++] = t;
        if (tNew.length < 2) return;

        for (int i = 0; i < samples.size(); i++) {
            for (int r = 0; r < samples.get(i).size(); r++) {
                for (SampledMetric sm : samples.get(i).get(r)) {
                    if (sm.t.length >= 2) {
                        sm.data = SplineInterpolator.interpolate(sm.t, sm.data, tNew);
                        sm.t = Arrays.copyOf(tNew, tNew.length);
                    }
                }
            }
        }

        for (int i = 0; i < samplesAggr.size(); i++) {
            for (SampledMetric sm : samplesAggr.get(i)) {
                if (sm.t.length >= 2) {
                    sm.data = SplineInterpolator.interpolate(sm.t, sm.data, tNew);
                    sm.t = Arrays.copyOf(tNew, tNew.length);
                }
            }
        }
    }

    public Matrix estimateAt(List<Station> nodes) {
        NetworkStruct sn = model.getStruct(false);

        Matrix estVal;
        String m = options.method;
        if ("ubr".equals(m)) estVal = estimatorUbr(nodes);
        else if ("ubo".equals(m)) estVal = estimatorUbo(nodes);
        else if ("erps".equals(m)) estVal = estimatorErps(nodes);
        else if ("ekf".equals(m)) estVal = estimatorEkf(nodes);
        else if ("mcmc".equals(m)) estVal = estimatorMcmc(nodes);
        else if ("mle".equals(m)) estVal = estimatorMle(nodes);
        else if ("mlps".equals(m)) estVal = estimatorMlps(nodes);
        else if ("fmlps".equals(m)) estVal = estimatorFmlps(nodes);
        else if ("qmle".equals(m)) estVal = estimatorQmle(nodes);
        else if ("gibbs".equals(m)) estVal = estimatorGibbs(nodes);
        else throw new IllegalArgumentException("Unknown inference method: " + m);

        List<JobClass> jobClasses = model.getJobClasses();
        for (int n = 0; n < nodes.size(); n++) {
            Station nd = nodes.get(n);
            if (nd instanceof Source) continue;
            for (int r = 0; r < sn.nclasses; r++) {
                double v = estVal.get(n, r);
                if (Double.isNaN(v) || v <= 0) continue;
                Object existingDist = (nd instanceof Queue) ? ((Queue) nd).getService(jobClasses.get(r)) : null;
                if (existingDist instanceof Markovian) {
                    ((Markovian) existingDist).setMean(v);
                } else if (existingDist instanceof Det) {
                    ((Det) existingDist).setMean(v);
                } else {
                    nd.setService(jobClasses.get(r), Exp.fitMean(v));
                }
            }
        }
        model.reset();

        return estVal;
    }

    private Matrix estimatorUbr(List<Station> nodes) {
        Station node = nodes.get(0);
        NetworkStruct sn = model.getStruct(false);
        int nClasses = sn.nclasses;
        List<JobClass> jobClasses = model.getJobClasses();

        double[] avgAggrUtil = null;
        if (node instanceof Queue) {
            int ns = ((Queue) node).getNumberOfServers();
            if (ns > 0 && ns < Integer.MAX_VALUE) {
                SampledMetric U = getAggrUtil(node);
                if (U != null) {
                    avgAggrUtil = new double[U.data.length];
                    for (int i = 0; i < U.data.length; i++) avgAggrUtil[i] = U.data[i] * ns;
                }
            }
        }

        boolean[] isUtilKnown = new boolean[nClasses];
        double[][] avgUtil = new double[nClasses][];
        double[][] avgArvR = new double[nClasses][];

        for (int r = 0; r < nClasses; r++) {
            SampledMetric Ur = getUtil(node, jobClasses.get(r));
            if (Ur != null) {
                int nServers = (node instanceof Queue) ? ((Queue) node).getNumberOfServers() : 1;
                avgUtil[r] = new double[Ur.data.length];
                for (int i = 0; i < Ur.data.length; i++) avgUtil[r][i] = Ur.data[i] * nServers;
                isUtilKnown[r] = true;
            }
            SampledMetric ar = getArvR(node, jobClasses.get(r));
            if (ar == null) throw new IllegalStateException("Arrival rate data missing for class " + r);
            avgArvR[r] = ar.data;
        }

        int nSamples = avgArvR[0].length;
        Matrix estVal = new Matrix(1, nClasses);

        for (int r = 0; r < nClasses; r++) {
            if (isUtilKnown[r]) {
                double[][] A = new double[nSamples][1];
                for (int i = 0; i < nSamples; i++) A[i][0] = avgArvR[r][i];
                double[] result = NnlsSolver.lsqnonneg(A, avgUtil[r]);
                estVal.set(0, r, result[0]);
            }
        }

        List<Integer> unknownClasses = new ArrayList<Integer>();
        for (int r = 0; r < nClasses; r++) if (!isUtilKnown[r]) unknownClasses.add(r);
        if (!unknownClasses.isEmpty() && avgAggrUtil != null) {
            double[] residualUtil = Arrays.copyOf(avgAggrUtil, avgAggrUtil.length);
            for (int r = 0; r < nClasses; r++) {
                if (isUtilKnown[r]) {
                    double[] sumUr = avgUtil[r];
                    for (int i = 0; i < residualUtil.length; i++) residualUtil[i] -= sumUr[i];
                }
            }
            double[][] A = new double[nSamples][unknownClasses.size()];
            for (int i = 0; i < nSamples; i++) {
                for (int j = 0; j < unknownClasses.size(); j++) {
                    A[i][j] = avgArvR[unknownClasses.get(j)][i];
                }
            }
            double[] result = NnlsSolver.lsqnonneg(A, residualUtil);
            for (int j = 0; j < unknownClasses.size(); j++) {
                estVal.set(0, unknownClasses.get(j), result[j]);
            }
        }

        return estVal;
    }

    private Matrix estimatorUbo(List<Station> nodes) {
        NetworkStruct sn = model.getStruct(false);
        int R = sn.nclasses;
        int M = nodes.size();
        int MR = M * R;
        List<JobClass> jobClasses = model.getJobClasses();

        double[][] avgU = new double[M][];
        double[][][] avgArvRData = new double[M][R][];
        double[][][] avgRespTData = new double[M][R][];

        for (int n = 0; n < M; n++) {
            Station node = nodes.get(n);
            int nServers = (node instanceof Queue) ? ((Queue) node).getNumberOfServers() : 1;
            SampledMetric U = getAggrUtil(node);
            if (U != null) {
                avgU[n] = new double[U.data.length];
                for (int i = 0; i < U.data.length; i++) avgU[n][i] = U.data[i] * nServers;
            }
            for (int r = 0; r < R; r++) {
                SampledMetric ar = getArvR(node, jobClasses.get(r));
                if (ar == null) throw new IllegalStateException("Arrival rate data missing");
                avgArvRData[n][r] = ar.data;
                SampledMetric rt = getRespT(node, jobClasses.get(r));
                if (rt == null) throw new IllegalStateException("Response time data missing");
                avgRespTData[n][r] = rt.data;
            }
        }

        int N = avgArvRData[0][0].length;

        Matrix H = new Matrix(MR, MR);
        Matrix h = new Matrix(MR, 1);

        for (int n = 0; n < N; n++) {
            double[] rhoN = new double[M];
            double[] betaN = new double[M];
            for (int i = 0; i < M; i++) {
                rhoN[i] = avgU[i][n];
                betaN[i] = 1.0 / (1.0 - rhoN[i]);
            }

            double[][] lambdaN = new double[M][R];
            double[][] RNarr = new double[M][R];
            for (int i = 0; i < M; i++) {
                for (int r = 0; r < R; r++) {
                    lambdaN[i][r] = avgArvRData[i][r][n];
                    RNarr[i][r] = avgRespTData[i][r][n];
                }
            }

            double[] lambdaR = new double[R];
            for (int r = 0; r < R; r++) {
                double s = 0.0;
                for (int i = 0; i < M; i++) s += lambdaN[i][r];
                lambdaR[r] = s;
            }
            double totalLambda = 0.0;
            for (int r = 0; r < R; r++) totalLambda += lambdaR[r];
            double[] wN = new double[R];
            for (int r = 0; r < R; r++) wN[r] = totalLambda > 0 ? lambdaR[r] / totalLambda : 1.0 / R;

            double[] EN = new double[R];
            for (int r = 0; r < R; r++) {
                double s = 0.0;
                for (int i = 0; i < M; i++) s += RNarr[i][r];
                EN[r] = s;
            }

            for (int r1 = 0; r1 < R; r1++) {
                for (int r2 = 0; r2 < R; r2++) {
                    for (int i1 = 0; i1 < M; i1++) {
                        for (int i2 = 0; i2 < M; i2++) {
                            int idx1 = r1 * M + i1;
                            int idx2 = r2 * M + i2;
                            if (r1 == r2) {
                                H.set(idx1, idx2, H.get(idx1, idx2) + 2.0 * wN[r1] * betaN[i1] * betaN[i2]);
                            }
                            if (i1 == i2) {
                                H.set(idx1, idx2, H.get(idx1, idx2) + 2.0 * lambdaN[i1][r1] * lambdaN[i2][r2]);
                            }
                        }
                    }
                }
            }

            for (int r = 0; r < R; r++) {
                for (int i = 0; i < M; i++) {
                    int idx = r * M + i;
                    h.set(idx, 0, h.get(idx, 0) - 2.0 * wN[r] * betaN[i] * EN[r]);
                    h.set(idx, 0, h.get(idx, 0) - 2.0 * lambdaN[i][r] * rhoN[i]);
                }
            }
        }

        Matrix lb = new Matrix(MR, 1);
        Pair<Matrix, Double> qpResult = OptimUtils.quadprog(H, h, lb);
        Matrix sVec = qpResult.getLeft();

        Matrix estVal = new Matrix(M, R);
        for (int r = 0; r < R; r++) {
            for (int i = 0; i < M; i++) {
                estVal.set(i, r, sVec.get(r * M + i, 0));
            }
        }
        return estVal;
    }

    private Matrix estimatorErps(List<Station> nodes) {
        Station node = nodes.get(0);
        NetworkStruct sn = model.getStruct(false);
        int R = sn.nclasses;
        List<JobClass> jobClasses = model.getJobClasses();

        if (!(node instanceof Queue) || sn.sched.get(sn.stations.get(((Queue) node).getStationIdx())) != SchedStrategy.PS) {
            throw new IllegalStateException("ERPS method is only available for PS stations.");
        }

        double[][] avgRespT = new double[R][];
        Matrix[] avgAQLen = new Matrix[R];

        for (int r = 0; r < R; r++) {
            SampledMetric rtData = getRespT(node, jobClasses.get(r));
            if (rtData == null) throw new IllegalStateException("Response time data missing for class " + r);
            avgRespT[r] = rtData.data;

            SampledMetric qlData = getAggrQLen(node, new ConditionEvent(node, jobClasses.get(r), EventType.ARV));
            if (qlData == null) throw new IllegalStateException("Arrival queue-length data missing for class " + r);
            Matrix m = new Matrix(qlData.data.length, 1);
            for (int i = 0; i < qlData.data.length; i++) m.set(i, 0, qlData.data[i]);
            avgAQLen[r] = m;
        }

        List<Double> busyCoresList = new ArrayList<Double>();
        for (int r = 0; r < R; r++) {
            for (int i = 0; i < avgAQLen[r].getNumRows(); i++) {
                busyCoresList.add(avgAQLen[r].get(i, 0));
            }
        }
        int nServers = (node instanceof Queue) ? ((Queue) node).getNumberOfServers() : 1;
        double sumBC = 0.0;
        for (Double d : busyCoresList) sumBC += d;
        double avgBC = busyCoresList.isEmpty() ? 0 : sumBC / busyCoresList.size();
        double avgBusyCores = Math.min(avgBC, (double) nServers);

        Matrix estVal = new Matrix(1, R);
        for (int r = 0; r < R; r++) {
            double[] respTimes = avgRespT[r];
            double[] totalQL = new double[respTimes.length];
            for (int i = 0; i < respTimes.length; i++) {
                totalQL[i] = avgAQLen[r].get(i, 0) / avgBusyCores;
            }

            double[][] A = new double[respTimes.length][1];
            for (int i = 0; i < respTimes.length; i++) A[i][0] = totalQL[i];
            double[] result = NnlsSolver.lsqnonneg(A, respTimes);
            estVal.set(0, r, result[0]);
        }
        return estVal;
    }

    private Matrix estimatorEkf(List<Station> nodes) {
        Station node = nodes.get(0);
        NetworkStruct sn = model.getStruct(false);
        int R = sn.nclasses;
        List<JobClass> jobClasses = model.getJobClasses();

        int nServers = (node instanceof Queue) ? ((Queue) node).getNumberOfServers() : 1;
        SampledMetric U = getAggrUtil(node);
        if (U == null) throw new IllegalStateException("Aggregate utilization data missing");
        double[] avgU = new double[U.data.length];
        for (int i = 0; i < U.data.length; i++) avgU[i] = U.data[i] * nServers;

        double[][] avgA = new double[R][];
        double[][] avgR = new double[R][];
        for (int r = 0; r < R; r++) {
            SampledMetric ar = getArvR(node, jobClasses.get(r));
            if (ar == null) throw new IllegalStateException("Arrival rate data missing for class " + r);
            avgA[r] = ar.data;
            SampledMetric rt = getRespT(node, jobClasses.get(r));
            if (rt == null) throw new IllegalStateException("Response time data missing for class " + r);
            avgR[r] = rt.data;
        }

        int N = avgU.length;

        final Function<NetworkStruct, SolverResult> solverAnalyzer = getSolverAnalyzer(model, options);

        double[] x;
        if (options.x0 != null) {
            x = Arrays.copyOf(options.x0, options.x0.length);
        } else {
            x = new double[R];
            for (int r = 0; r < R; r++) {
                double mx = 0;
                for (double v : avgR[r]) if (v > mx) mx = v;
                x[r] = Math.random() * mx;
            }
        }
        double[][] p = new double[R][R];
        for (int i = 0; i < R; i++) p[i][i] = x[i] * x[i];

        double[][] mCovNoise = new double[R + 1][R + 1];
        for (int i = 0; i < R + 1; i++) mCovNoise[i][i] = 0.01;
        double[][] pCovNoise = new double[R][R];
        for (int i = 0; i < R; i++) pCovNoise[i][i] = 0.001;

        int stIdx = (node instanceof Station) ? ((Station) node).getStationIdx() : 0;
        double stepBound = 0.6;

        int iters = Math.min(N, options.iterMax);
        for (int n = 0; n < iters; n++) {
            double[] xN = Arrays.copyOf(x, x.length);
            double[][] pN = new double[R][R];
            for (int i = 0; i < R; i++) for (int j = 0; j < R; j++) pN[i][j] = p[i][j] + pCovNoise[i][j];

            double stepUtil = avgU[n];
            double[] stepResponse = new double[R];
            for (int r = 0; r < R; r++) stepResponse[r] = avgR[r][n];

            double[] zN = getPredictedMeasurement(xN, R, sn, stIdx, solverAnalyzer);
            double[][] HN = getJacobian(xN, R, sn, stIdx, solverAnalyzer, zN);
            double[] z = new double[R + 1];
            for (int i = 0; i < R + 1; i++) z[i] = (i < R) ? stepResponse[i] : stepUtil;
            double[] yN = new double[R + 1];
            for (int i = 0; i < R + 1; i++) yN[i] = z[i] - zN[i];

            double[][] S = new double[R + 1][R + 1];
            for (int i = 0; i < R + 1; i++) {
                for (int j = 0; j < R + 1; j++) {
                    double sum = mCovNoise[i][j];
                    for (int k = 0; k < R; k++) {
                        for (int l = 0; l < R; l++) {
                            sum += HN[i][k] * pN[k][l] * HN[j][l];
                        }
                    }
                    S[i][j] = sum;
                }
            }

            double[][] SInv = invertMatrix(S);
            double[][] K = new double[R][R + 1];
            for (int i = 0; i < R; i++) {
                for (int j = 0; j < R + 1; j++) {
                    double sum = 0.0;
                    for (int k = 0; k < R; k++) {
                        for (int l = 0; l < R + 1; l++) {
                            sum += pN[i][k] * HN[l][k] * SInv[l][j];
                        }
                    }
                    K[i][j] = sum;
                }
            }

            x = new double[R];
            for (int i = 0; i < R; i++) {
                double sum = xN[i];
                for (int j = 0; j < R + 1; j++) sum += K[i][j] * yN[j];
                double lower = stepBound * 0.0 + (1.0 - stepBound) * sum;
                double upper = stepBound * Double.MAX_VALUE + (1.0 - stepBound) * sum;
                x[i] = Math.min(upper, Math.max(lower, sum));
            }

            double xs = 0.0;
            for (double v : x) xs += v;
            if (xs < 0) {
                for (int i = 0; i < x.length; i++) x[i] = -x[i];
            }

            double[][] IKH = new double[R][R];
            for (int i = 0; i < R; i++) {
                for (int j = 0; j < R; j++) {
                    double kh = 0.0;
                    for (int k = 0; k < R + 1; k++) kh += K[i][k] * HN[k][j];
                    IKH[i][j] = ((i == j) ? 1.0 : 0.0) - kh;
                }
            }
            double[][] pNew = new double[R][R];
            for (int i = 0; i < R; i++) {
                for (int j = 0; j < R; j++) {
                    double sum = 0.0;
                    for (int k = 0; k < R; k++) sum += IKH[i][k] * pN[k][j];
                    pNew[i][j] = sum;
                }
            }
            p = pNew;
        }

        Matrix estVal = new Matrix(1, R);
        for (int r = 0; r < R; r++) estVal.set(0, r, x[r]);
        return estVal;
    }

    private Matrix estimatorMcmc(List<Station> nodes) {
        NetworkStruct sn = model.getStruct(false);
        int R = sn.nclasses;
        List<JobClass> jobClasses = model.getJobClasses();

        int Nopen = options.openPopulation;

        double[] P = new double[R];
        for (int r = 0; r < R; r++) {
            double njobs = sn.njobs.get(r);
            P[r] = (njobs < Double.MAX_VALUE / 2) ? njobs : (double) Nopen;
        }

        double[] Z = new double[R];
        List<Node> allNodes = model.nodes;
        for (Node nd : allNodes) {
            if (nd instanceof Delay) {
                for (int r = 0; r < R; r++) {
                    if (sn.njobs.get(r) < Double.MAX_VALUE / 2) {
                        Object dist = ((Delay) nd).getService(jobClasses.get(r));
                        if (dist != null) {
                            try {
                                Z[r] += (Double) dist.getClass().getMethod("getMean").invoke(dist);
                            } catch (Exception e) {}
                        }
                    }
                }
            } else if (nd instanceof Source) {
                for (int r = 0; r < R; r++) {
                    if (sn.njobs.get(r) >= Double.MAX_VALUE / 2) {
                        Object dist = ((Source) nd).getArrivalProcess(jobClasses.get(r));
                        if (dist != null) {
                            try {
                                double mean = (Double) dist.getClass().getMethod("getMean").invoke(dist);
                                double lambdaR = 1.0 / mean;
                                Z[r] = P[r] / lambdaR;
                            } catch (Exception e) {}
                        }
                    }
                }
            }
        }

        int M = nodes.size();
        List<double[]> avgQLList = new ArrayList<double[]>();
        for (int n = 0; n < M; n++) {
            SampledMetric qlData = getAggrQLen(nodes.get(n));
            if (qlData == null) throw new IllegalStateException("Queue-length data missing for node " + n);
            avgQLList.add(qlData.data);
        }

        int experiments = avgQLList.get(0).length;
        Matrix avgQL = new Matrix(M, R);
        for (int i = 0; i < M; i++) {
            double[] arr = avgQLList.get(i);
            double avg = 0.0;
            for (double v : arr) avg += v;
            if (arr.length > 0) avg /= arr.length;
            for (int r = 0; r < R; r++) {
                avgQL.set(i, r, avg);
            }
        }

        Map<Integer, Matrix> visits = sn.visits;
        return mcmcData(avgQL, visits, experiments, options.iterMax, P, Z);
    }

    private Matrix estimatorMle(List<Station> nodes) {
        Station node = nodes.get(0);
        final NetworkStruct sn = model.getStruct(false);
        final int R = sn.nclasses;
        List<JobClass> jobClasses = model.getJobClasses();

        int nServers = (node instanceof Queue) ? ((Queue) node).getNumberOfServers() : 1;
        SampledMetric U = getAggrUtil(node);
        if (U == null) throw new IllegalStateException("Aggregate utilization data missing");
        final double[] avgU = new double[U.data.length];
        for (int i = 0; i < U.data.length; i++) avgU[i] = U.data[i] * nServers;

        final double[][] avgA = new double[R][];
        final double[][] avgRT = new double[R][];
        for (int r = 0; r < R; r++) {
            SampledMetric ar = getArvR(node, jobClasses.get(r));
            if (ar == null) throw new IllegalStateException("Arrival rate data missing for class " + r);
            avgA[r] = ar.data;
            SampledMetric rt = getRespT(node, jobClasses.get(r));
            if (rt == null) throw new IllegalStateException("Response time data missing for class " + r);
            avgRT[r] = rt.data;
        }

        final int N = avgU.length;
        final Function<NetworkStruct, SolverResult> solverAnalyzer = getSolverAnalyzer(model, options);
        final int stIdx = (node instanceof Station) ? ((Station) node).getStationIdx() : 0;

        double[] maxRT = new double[R];
        for (int r = 0; r < R; r++) {
            double mx = 0;
            for (double v : avgRT[r]) if (v > mx) mx = v;
            maxRT[r] = mx;
        }
        double[] x0 = new double[R];
        for (int r = 0; r < R; r++) x0[r] = Math.random() * maxRT[r];
        double[] xLB = new double[R];
        for (int r = 0; r < R; r++) xLB[r] = 1e-8;
        double[] xUB = maxRT;

        final double[][] w = new double[N][R];
        for (int i = 0; i < N; i++) {
            double totalArvR = 0.0;
            for (int r = 0; r < R; r++) totalArvR += avgA[r][i];
            for (int r = 0; r < R; r++) {
                w[i][r] = totalArvR > 0 ? avgA[r][i] / totalArvR : 1.0 / R;
            }
        }

        Function<double[], Double> objFun = new Function<double[], Double>() {
            @Override
            public Double apply(double[] x) {
                for (int c = 0; c < R; c++) {
                    Sn_set_service_coc.sn_set_service_coc(sn, stIdx, c, 1.0 / x[c]);
                }
                SolverResult result = solverAnalyzer.apply(sn);
                double[] predR = new double[R];
                for (int r = 0; r < R; r++) predR[r] = result.RN.get(stIdx, r);
                double predU = 0.0;
                for (int r = 0; r < R; r++) predU += result.UN.get(stIdx, r);

                double f = 0.0;
                for (int i = 0; i < N; i++) {
                    for (int r = 0; r < R; r++) {
                        double delta = predR[r] - avgRT[r][i];
                        f += w[i][r] * delta * delta;
                    }
                    double eps = predU - avgU[i];
                    f += eps * eps;
                }
                return f;
            }
        };

        Pair<double[], Double> fminResult = OptimUtils.fmincon(objFun, x0, xLB, xUB, options.iterMax);
        double[] demEst = fminResult.getLeft();

        Matrix estVal = new Matrix(1, R);
        for (int r = 0; r < R; r++) estVal.set(0, r, demEst[r]);
        return estVal;
    }

    private Matrix estimatorMlps(List<Station> nodes) {
        Queue node = (Queue) nodes.get(0);
        NetworkStruct sn = model.getStruct(false);
        int R = sn.nclasses;
        List<JobClass> jobClasses = model.getJobClasses();

        if (sn.sched.get(sn.stations.get(node.getStationIdx())) != SchedStrategy.PS) {
            throw new IllegalStateException("MLPS method is only available for PS stations.");
        }

        boolean hasOpen = false;
        for (int r = 0; r < R; r++) if (sn.njobs.get(r) >= Double.MAX_VALUE / 2) { hasOpen = true; break; }
        Network eqModel;
        Queue eqNode;
        if (hasOpen) {
            Pair<Network, Queue> eq = buildClosedEquivalentForPS(node);
            eqModel = eq.getLeft();
            eqNode = eq.getRight();
        } else {
            eqModel = model;
            eqNode = node;
        }

        List<Double> rtAll = new ArrayList<Double>();
        List<Integer> classAll = new ArrayList<Integer>();
        List<Double> atAll = new ArrayList<Double>();

        for (int r = 0; r < R; r++) {
            SampledMetric arvData = getArvR(node, jobClasses.get(r));
            if (arvData == null) throw new IllegalStateException("Arrival data missing for class " + r);
            if (!arvData.isTrace()) throw new IllegalStateException("MLPS requires trace-format data.");

            SampledMetric rtData = getRespT(node, jobClasses.get(r));
            if (rtData == null) throw new IllegalStateException("Response time data missing for class " + r);
            if (!rtData.isTrace()) throw new IllegalStateException("MLPS requires trace-format data.");

            for (int i = 0; i < rtData.data.length; i++) {
                rtAll.add(rtData.data[i]);
                atAll.add(arvData.data[i]);
                classAll.add(r);
            }
        }

        double[] rt = new double[rtAll.size()];
        double[] at = new double[atAll.size()];
        int[] classVec = new int[classAll.size()];
        for (int i = 0; i < rt.length; i++) {
            rt[i] = rtAll.get(i);
            at[i] = atAll.get(i);
            classVec[i] = classAll.get(i);
        }
        int n = at.length;
        int[] jobid = new int[n];
        for (int i = 0; i < n; i++) jobid[i] = i;

        Matrix ql = Infer_compute_ql_at_arrival.infer_compute_ql_at_arrival(at, jobid, rt, jobid, classVec, R);

        Integer[] sortIdxBoxed = new Integer[n];
        for (int i = 0; i < n; i++) sortIdxBoxed[i] = i;
        final double[] atFinal = at;
        Arrays.sort(sortIdxBoxed, new Comparator<Integer>() {
            @Override
            public int compare(Integer a, Integer b) { return Double.compare(atFinal[a], atFinal[b]); }
        });
        int[] sortIdx = new int[n];
        for (int i = 0; i < n; i++) sortIdx[i] = sortIdxBoxed[i];

        double[] rtSorted = new double[n];
        int[] classSorted = new int[n];
        for (int i = 0; i < n; i++) {
            rtSorted[i] = rt[sortIdx[i]];
            classSorted[i] = classVec[sortIdx[i]];
        }
        Matrix qlSorted = new Matrix(n, R);
        for (int i = 0; i < n; i++) {
            for (int r = 0; r < R; r++) qlSorted.set(i, r, ql.get(sortIdx[i], r));
        }

        List<Integer> valid = new ArrayList<Integer>();
        for (int i = 0; i < n; i++) if (rtSorted[i] > 0) valid.add(i);
        double[] rtValid = new double[valid.size()];
        int[] classValid = new int[valid.size()];
        for (int i = 0; i < valid.size(); i++) {
            rtValid[i] = rtSorted[valid.get(i)];
            classValid[i] = classSorted[valid.get(i)];
        }
        Matrix qlValid = new Matrix(valid.size(), R);
        for (int i = 0; i < valid.size(); i++) {
            for (int r = 0; r < R; r++) qlValid.set(i, r, qlSorted.get(valid.get(i), r));
        }

        double[] demandEst = Infer_mlps.infer_mlps(eqModel, eqNode, rtValid, classValid, qlValid);
        Matrix estVal = new Matrix(1, R);
        for (int r = 0; r < R; r++) estVal.set(0, r, demandEst[r]);
        return estVal;
    }

    private Matrix estimatorFmlps(List<Station> nodes) {
        Queue node = (Queue) nodes.get(0);
        NetworkStruct sn = model.getStruct(false);
        int R = sn.nclasses;
        List<JobClass> jobClasses = model.getJobClasses();

        if (sn.sched.get(sn.stations.get(node.getStationIdx())) != SchedStrategy.PS) {
            throw new IllegalStateException("FMLPS method is only available for PS stations.");
        }

        boolean hasOpen = false;
        for (int r = 0; r < R; r++) if (sn.njobs.get(r) >= Double.MAX_VALUE / 2) { hasOpen = true; break; }
        Network eqModel;
        Queue eqNode;
        int W;
        if (hasOpen) {
            Pair<Network, Queue> eq = buildClosedEquivalentForPS(node);
            eqModel = eq.getLeft();
            eqNode = eq.getRight();
            NetworkStruct eqSn = eqModel.getStruct(false);
            W = (int) eqSn.njobs.elementSum();
        } else {
            eqModel = model;
            eqNode = node;
            W = (int) sn.njobs.elementSum();
        }

        List<Double> rtAll = new ArrayList<Double>();
        List<Integer> classAll = new ArrayList<Integer>();
        List<Double> atAll = new ArrayList<Double>();

        for (int r = 0; r < R; r++) {
            SampledMetric arvData = getArvR(node, jobClasses.get(r));
            if (arvData == null) throw new IllegalStateException("Arrival data missing for class " + r);
            SampledMetric rtData = getRespT(node, jobClasses.get(r));
            if (rtData == null) throw new IllegalStateException("Response time data missing for class " + r);

            for (int i = 0; i < rtData.data.length; i++) {
                rtAll.add(rtData.data[i]);
                atAll.add(arvData.data[i]);
                classAll.add(r);
            }
        }

        double[] rt = new double[rtAll.size()];
        double[] at = new double[atAll.size()];
        int[] classVec = new int[classAll.size()];
        for (int i = 0; i < rt.length; i++) {
            rt[i] = rtAll.get(i);
            at[i] = atAll.get(i);
            classVec[i] = classAll.get(i);
        }
        int n = at.length;
        int[] jobid = new int[n];
        for (int i = 0; i < n; i++) jobid[i] = i;

        Matrix ql = Infer_compute_ql_at_arrival.infer_compute_ql_at_arrival(at, jobid, rt, jobid, classVec, R);

        Integer[] sortIdxBoxed = new Integer[n];
        for (int i = 0; i < n; i++) sortIdxBoxed[i] = i;
        final double[] atFinal = at;
        Arrays.sort(sortIdxBoxed, new Comparator<Integer>() {
            @Override
            public int compare(Integer a, Integer b) { return Double.compare(atFinal[a], atFinal[b]); }
        });
        int[] sortIdx = new int[n];
        for (int i = 0; i < n; i++) sortIdx[i] = sortIdxBoxed[i];

        double[] rtSorted = new double[n];
        int[] classSorted = new int[n];
        for (int i = 0; i < n; i++) {
            rtSorted[i] = rt[sortIdx[i]];
            classSorted[i] = classVec[sortIdx[i]];
        }
        Matrix qlSorted = new Matrix(n, R);
        for (int i = 0; i < n; i++) {
            for (int r = 0; r < R; r++) qlSorted.set(i, r, ql.get(sortIdx[i], r));
        }

        List<Integer> valid = new ArrayList<Integer>();
        for (int i = 0; i < n; i++) if (rtSorted[i] > 0) valid.add(i);
        double[] rtValid = new double[valid.size()];
        int[] classValid = new int[valid.size()];
        for (int i = 0; i < valid.size(); i++) {
            rtValid[i] = rtSorted[valid.get(i)];
            classValid[i] = classSorted[valid.get(i)];
        }
        Matrix qlValid = new Matrix(valid.size(), R);
        for (int i = 0; i < valid.size(); i++) {
            for (int r = 0; r < R; r++) qlValid.set(i, r, qlSorted.get(valid.get(i), r));
        }

        double[] demandEst = Infer_fmlps.infer_fmlps(eqModel, eqNode, rtValid, classValid, qlValid, W);
        Matrix estVal = new Matrix(1, R);
        for (int r = 0; r < R; r++) estVal.set(0, r, demandEst[r]);
        return estVal;
    }

    private Matrix estimatorQmle(List<Station> nodes) {
        NetworkStruct sn = model.getStruct(false);
        int R = sn.nclasses;
        int M = nodes.size();
        List<JobClass> jobClasses = model.getJobClasses();

        int Nopen = options.openPopulation;
        double[] N = new double[R];
        for (int r = 0; r < R; r++) {
            double njobs = sn.njobs.get(r);
            N[r] = (njobs < Double.MAX_VALUE / 2) ? njobs : (double) Nopen;
        }

        double[] Z = new double[R];
        List<Node> allNodes = model.nodes;
        for (Node nd : allNodes) {
            if (nd instanceof Delay) {
                for (int r = 0; r < R; r++) {
                    if (sn.njobs.get(r) < Double.MAX_VALUE / 2) {
                        Object dist = ((Delay) nd).getService(jobClasses.get(r));
                        if (dist != null) {
                            try {
                                Z[r] += (Double) dist.getClass().getMethod("getMean").invoke(dist);
                            } catch (Exception e) {}
                        }
                    }
                }
            } else if (nd instanceof Source) {
                for (int r = 0; r < R; r++) {
                    if (sn.njobs.get(r) >= Double.MAX_VALUE / 2) {
                        Object dist = ((Source) nd).getArrivalProcess(jobClasses.get(r));
                        if (dist != null) {
                            try {
                                double mean = (Double) dist.getClass().getMethod("getMean").invoke(dist);
                                double lambdaR = 1.0 / mean;
                                Z[r] = N[r] / lambdaR;
                            } catch (Exception e) {}
                        }
                    }
                }
            }
        }

        Matrix Q = new Matrix(M, R);
        for (int n = 0; n < M; n++) {
            for (int r = 0; r < R; r++) {
                List<SampledMetric> qlData = getQLen(nodes.get(n), jobClasses.get(r));
                if (qlData.isEmpty()) throw new IllegalStateException("Queue-length data missing");
                double avg = 0.0;
                double[] arr = qlData.get(0).data;
                for (double v : arr) avg += v;
                if (arr.length > 0) avg /= arr.length;
                Q.set(n, r, avg);
            }
        }

        return Infer_qmle.infer_qmle(Q, N, Z);
    }

    private Matrix estimatorGibbs(List<Station> nodes) {
        Station node = nodes.get(0);
        NetworkStruct sn = model.getStruct(false);
        int R = sn.nclasses;
        int nServers = (node instanceof Queue) ? ((Queue) node).getNumberOfServers() : 1;
        List<JobClass> jobClasses = model.getJobClasses();

        double[][][] data = new double[7][R + 1][];

        for (int r = 0; r < R; r++) {
            SampledMetric arvData = getArvR(node, jobClasses.get(r));
            if (arvData == null) throw new IllegalStateException("Arrival data missing for class " + r);
            if (!arvData.isTrace()) throw new IllegalStateException("Gibbs requires trace-format data.");
            double[] arvScaled = new double[arvData.data.length];
            for (int i = 0; i < arvData.data.length; i++) arvScaled[i] = arvData.data[i] * 1000;
            data[3][r] = arvScaled;

            SampledMetric rtData = getRespT(node, jobClasses.get(r));
            if (rtData == null) throw new IllegalStateException("Response time data missing for class " + r);
            if (!rtData.isTrace()) throw new IllegalStateException("Gibbs requires trace-format data.");
            data[4][r] = rtData.data;

            SampledMetric tputData = getTput(node, jobClasses.get(r));
            if (tputData == null) throw new IllegalStateException("Throughput data missing for class " + r);
            data[6][r] = tputData.data;
        }

        double[] demandEst = Infer_gibbs.infer_gibbs(data, nServers, options.tol);
        Matrix estVal = new Matrix(1, R);
        int lim = Math.min(R, demandEst.length);
        for (int r = 0; r < lim; r++) {
            estVal.set(0, r, demandEst[r]);
        }
        return estVal;
    }

    public Pair<Network, Queue> buildClosedEquivalentForPS(Queue node) {
        NetworkStruct sn = model.getStruct(false);
        int R = sn.nclasses;
        List<JobClass> jobClasses = model.getJobClasses();
        int Nopen = options.openPopulation;

        int[] N = new int[R];
        double[] Z = new double[R];

        List<Node> allNodes = model.nodes;
        for (Node nd : allNodes) {
            if (nd instanceof Delay) {
                for (int r = 0; r < R; r++) {
                    if (sn.njobs.get(r) < Double.MAX_VALUE / 2) {
                        Object dist = ((Delay) nd).getService(jobClasses.get(r));
                        if (dist != null) {
                            try {
                                Z[r] += (Double) dist.getClass().getMethod("getMean").invoke(dist);
                            } catch (Exception e) {}
                        }
                    }
                }
            } else if (nd instanceof Source) {
                for (int r = 0; r < R; r++) {
                    if (sn.njobs.get(r) >= Double.MAX_VALUE / 2) {
                        Object dist = ((Source) nd).getArrivalProcess(jobClasses.get(r));
                        if (dist != null) {
                            try {
                                double mean = (Double) dist.getClass().getMethod("getMean").invoke(dist);
                                double lambdaR = 1.0 / mean;
                                Z[r] = (double) Nopen / lambdaR;
                            } catch (Exception e) {}
                        }
                    }
                }
            }
        }

        for (int r = 0; r < R; r++) {
            if (sn.njobs.get(r) < Double.MAX_VALUE / 2) {
                N[r] = (int) sn.njobs.get(r);
            } else {
                N[r] = Nopen;
            }
        }

        double[] delayRate = new double[R];
        for (int r = 0; r < R; r++) delayRate[r] = 1.0 / Z[r];

        Network eqModel = new Network("closed_equiv");
        Delay eqDelay = new Delay(eqModel, "Think");
        Queue eqQueue = new Queue(eqModel, "Queue1", SchedStrategy.PS);
        eqQueue.setNumberOfServers(node.getNumberOfServers());

        ClosedClass[] eqClasses = new ClosedClass[R];
        for (int r = 0; r < R; r++) {
            eqClasses[r] = new ClosedClass(eqModel, "Class" + (r + 1), N[r], eqDelay, 0);
            eqDelay.setService(eqClasses[r], new Exp(delayRate[r]));
            Object origDist = node.getService(jobClasses.get(r));
            if (origDist != null) {
                eqQueue.setService(eqClasses[r], (Distribution) origDist);
            } else {
                eqQueue.setService(eqClasses[r], new Exp(1.0));
            }
        }

        Object P = eqModel.initRoutingMatrix();
        List<Node> nodeList = new ArrayList<Node>();
        nodeList.add(eqDelay);
        nodeList.add(eqQueue);
        for (int r = 0; r < R; r++) {
            try {
                java.lang.reflect.Method setM = P.getClass().getMethod("set", JobClass.class, Object.class);
                setM.invoke(P, eqClasses[r], Network.serialRouting(nodeList));
            } catch (Exception e) {}
        }
        eqModel.link((RoutingMatrix) P);

        return new Pair<Network, Queue>(eqModel, eqQueue);
    }

    private Function<NetworkStruct, SolverResult> getSolverAnalyzer(final Network model, final EstimatorOptions options) {
        Function<Network, NetworkSolver> solverFactory = options.solverFactory;
        if (solverFactory != null) {
            final NetworkSolver solver = solverFactory.apply(model);
            final Object solverOpts = solver.options;

            if (solver instanceof SolverMVA) {
                return new Function<NetworkStruct, SolverResult>() {
                    @Override public SolverResult apply(NetworkStruct sn) {
                        return Solver_mva_analyzer.solver_mva_analyzer(sn, (jline.solvers.SolverOptions) solverOpts);
                    }
                };
            } else {
                return new Function<NetworkStruct, SolverResult>() {
                    @Override public SolverResult apply(NetworkStruct sn) {
                        Object mvaOpts = SolverMVA.defaultOptions();
                        return Solver_mva_analyzer.solver_mva_analyzer(sn, (jline.solvers.SolverOptions) mvaOpts);
                    }
                };
            }
        }

        final Object mvaOpts = SolverMVA.defaultOptions();
        return new Function<NetworkStruct, SolverResult>() {
            @Override public SolverResult apply(NetworkStruct sn) {
                return Solver_mva_analyzer.solver_mva_analyzer(sn, (jline.solvers.SolverOptions) mvaOpts);
            }
        };
    }

    private double[] getPredictedMeasurement(double[] x, int R, NetworkStruct sn, int stIdx,
                                              Function<NetworkStruct, SolverResult> solverAnalyzer) {
        for (int c = 0; c < R; c++) {
            Sn_set_service_coc.sn_set_service_coc(sn, stIdx, c, 1.0 / x[c]);
        }

        SolverResult result = solverAnalyzer.apply(sn);
        double[] h = new double[R + 1];
        for (int c = 0; c < R; c++) h[c] = result.RN.get(stIdx, c);
        double sum = 0.0;
        for (int r = 0; r < R; r++) sum += result.UN.get(stIdx, r);
        h[R] = sum;
        return h;
    }

    private double[][] getJacobian(double[] x, int R, NetworkStruct sn, int stIdx,
                                    Function<NetworkStruct, SolverResult> solverAnalyzer, double[] h0) {
        double delta = 1e-6;
        double[][] Hx = new double[R + 1][R];

        for (int c = 0; c < R; c++) {
            double[] xPert = Arrays.copyOf(x, x.length);
            xPert[c] += delta;
            NetworkStruct snPert = model.getStruct(true);
            double[] hPert = getPredictedMeasurement(xPert, R, snPert, stIdx, solverAnalyzer);
            for (int i = 0; i <= R; i++) {
                Hx[i][c] = (hPert[i] - h0[i]) / delta;
            }
        }
        return Hx;
    }

    private Matrix mcmcData(Matrix avgQL, Map<Integer, Matrix> visits, int experiments, int iterMax,
                             double[] P, double[] Z) {
        int M = avgQL.getNumRows();
        int R = avgQL.getNumCols();
        int S = 100;

        String mciVariant = "imci";
        double maxQL = 0.0;
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < R; j++) {
                if (avgQL.get(i, j) > maxQL) maxQL = avgQL.get(i, j);
            }
        }
        double thetaStep = maxQL / 400;
        List<Double> stepsList = new ArrayList<Double>();
        for (double t = 0.0; t <= maxQL; t += thetaStep) stepsList.add(t);
        double[] steps = new double[stepsList.size()];
        for (int i = 0; i < steps.length; i++) steps[i] = stepsList.get(i);

        double[][][] theta = new double[S + 1][M][R];

        Random rng = new Random();
        for (int s = 0; s < S; s++) {
            double[][] sampleTheta = new double[M][R];
            for (int i = 0; i < M; i++) sampleTheta[i] = Arrays.copyOf(theta[s][i], R);
            for (int i = 0; i < M; i++) {
                for (int c = 0; c < R; c++) {
                    Matrix thetaMatrix = new Matrix(M, R);
                    for (int ii = 0; ii < M; ii++) {
                        for (int cc = 0; cc < R; cc++) {
                            thetaMatrix.set(ii, cc, sampleTheta[ii][cc]);
                        }
                    }
                    Matrix popMatrix = new Matrix(1, R);
                    for (int cc = 0; cc < R; cc++) popMatrix.set(0, cc, P[cc]);
                    Matrix thinkMatrix = new Matrix(1, R);
                    for (int cc = 0; cc < R; cc++) thinkMatrix.set(0, cc, Z[cc]);

                    Object gResult = jline.api.pfqn.nc.Pfqn_mci.pfqn_mci(thetaMatrix, popMatrix, thinkMatrix, experiments, mciVariant);
                    double logG;
                    try {
                        java.lang.reflect.Field lGf = gResult.getClass().getField("lG");
                        java.lang.reflect.Field Gf = gResult.getClass().getField("G");
                        Object lGv = lGf.get(gResult);
                        if (lGv != null && Double.isFinite(((Number) lGv).doubleValue())) {
                            logG = ((Number) lGv).doubleValue();
                        } else {
                            logG = Math.log(((Number) Gf.get(gResult)).doubleValue() + 1e-15);
                        }
                    } catch (Exception ex) {
                        logG = 0.0;
                    }
                    double logPrior = Math.log(thetaStep / maxQL);

                    double[] logPosteriors = new double[steps.length];
                    for (int st = 0; st < steps.length; st++) {
                        double stepTheta = steps[st];
                        logPosteriors[st] = experiments * avgQL.get(i, c) * Math.log(stepTheta + 1e-15) -
                                experiments * logG + logPrior;
                    }

                    double maxLogPost = Double.NEGATIVE_INFINITY;
                    for (double v : logPosteriors) if (v > maxLogPost) maxLogPost = v;
                    double[] probs = new double[steps.length];
                    double probSum = 0.0;
                    for (int j = 0; j < steps.length; j++) {
                        probs[j] = Math.exp(logPosteriors[j] - maxLogPost);
                        probSum += probs[j];
                    }
                    for (int j = 0; j < probs.length; j++) probs[j] /= probSum;

                    double[] cumProb = new double[probs.length];
                    cumProb[0] = probs[0];
                    for (int j = 1; j < probs.length; j++) cumProb[j] = cumProb[j - 1] + probs[j];

                    double u = rng.nextDouble();
                    int index = -1;
                    for (int j = 0; j < cumProb.length; j++) if (cumProb[j] > u) { index = j; break; }
                    if (index < 0) index = cumProb.length - 1;
                    sampleTheta[i][c] = steps[index];
                }
            }
            for (int i = 0; i < M; i++) {
                theta[s + 1][i] = Arrays.copyOf(sampleTheta[i], R);
            }
        }

        if (visits != null) {
            for (int s = 0; s <= S; s++) {
                for (int i = 0; i < M; i++) {
                    for (int c = 0; c < R; c++) {
                        Matrix chainVisits = visits.get(c);
                        if (chainVisits != null) {
                            double visitVal = chainVisits.get(i + 1, c);
                            if (visitVal != 0.0) {
                                theta[s][i][c] /= visitVal;
                            }
                        }
                    }
                }
            }
        }

        int cutoff = S / 2;
        Matrix thetaAvg = new Matrix(M, R);
        for (int i = 0; i < M; i++) {
            for (int c = 0; c < R; c++) {
                double sum = 0.0;
                for (int s = cutoff; s <= S; s++) {
                    sum += theta[s][i][c];
                }
                thetaAvg.set(i, c, sum / (S - cutoff + 1));
            }
        }

        return thetaAvg;
    }

    private double[][] invertMatrix(double[][] A) {
        int n = A.length;
        double[][] aug = new double[n][2 * n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < 2 * n; j++) {
                if (j < n) aug[i][j] = A[i][j];
                else if (j - n == i) aug[i][j] = 1.0;
                else aug[i][j] = 0.0;
            }
        }

        for (int col = 0; col < n; col++) {
            int maxRow = col;
            double maxVal = Math.abs(aug[col][col]);
            for (int row = col + 1; row < n; row++) {
                if (Math.abs(aug[row][col]) > maxVal) {
                    maxVal = Math.abs(aug[row][col]);
                    maxRow = row;
                }
            }
            double[] temp = aug[col]; aug[col] = aug[maxRow]; aug[maxRow] = temp;

            double pivot = aug[col][col];
            if (Math.abs(pivot) < 1e-14) continue;
            for (int j = 0; j < 2 * n; j++) aug[col][j] /= pivot;

            for (int row = 0; row < n; row++) {
                if (row == col) continue;
                double factor = aug[row][col];
                for (int j = 0; j < 2 * n; j++) aug[row][j] -= factor * aug[col][j];
            }
        }

        double[][] inv = new double[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) inv[i][j] = aug[i][j + n];
        }
        return inv;
    }

    public static EstimatorOptions defaultOptions() {
        return EstimatorOptions.defaultOptions();
    }

    public static String getRequiredMetrics(String method) {
        if ("ubr".equals(method)) return "ArvR (per-class) + Util (per-class or aggregate)";
        if ("ubo".equals(method)) return "ArvR (per-class) + RespT (per-class) + Util (aggregate)";
        if ("erps".equals(method)) return "RespT (per-class) + QLen (aggregate, conditional on class arrivals). PS stations only.";
        if ("ekf".equals(method)) return "RespT (per-class) + Util (aggregate). Sequential/recursive estimation.";
        if ("mcmc".equals(method)) return "QLen (aggregate). Gibbs sampling with MCMC. Open/mixed via closed equivalence.";
        if ("mle".equals(method)) return "ArvR (per-class) + RespT (per-class) + Util (aggregate)";
        if ("mlps".equals(method)) return "ArvR (per-class, trace) + RespT (per-class, trace). PS stations only.";
        if ("fmlps".equals(method)) return "ArvR (per-class, trace) + RespT (per-class, trace). PS stations only.";
        if ("qmle".equals(method)) return "QLen (per-class). Open/mixed via closed equivalence.";
        if ("gibbs".equals(method)) return "ArvR (per-class, trace) + RespT (per-class, trace) + Tput (per-class). Gibbs sampling.";
        return "Unknown method: " + method;
    }
}
