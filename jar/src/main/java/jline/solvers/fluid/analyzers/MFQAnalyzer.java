package jline.solvers.fluid.analyzers;

import jline.GlobalConstants;
import jline.lang.NetworkStruct;
import jline.lang.constant.NodeType;
import jline.lib.butools.queues.FluFluQueue;
import jline.lib.butools.queues.FluFluResult;
import jline.solvers.SolverOptions;
import jline.solvers.SolverResult;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

/**
 * MFQ (Markovian Fluid Queue) analyzer for single-queue open systems.
 */
public class MFQAnalyzer implements FluidAnalyzer {
    private Matrix xvecIt;

    @Override
    public void analyze(NetworkStruct sn, SolverOptions options, SolverResult result) {
        int M = sn.nstations;
        int K = sn.nclasses;
        result.QN = Matrix.zeros(M, K);
        result.UN = Matrix.zeros(M, K);
        result.RN = Matrix.zeros(M, K);
        result.TN = Matrix.zeros(M, K);
        result.WN = new Matrix(0, 0);
        result.AN = new Matrix(0, 0);

        TopologyInfo topologyInfo = validateTopology(sn);
        if (!topologyInfo.isValid) {
            throw new RuntimeException("MFQ requires single-queue topology: " + topologyInfo.errorMsg);
        }
        int sourceStation = topologyInfo.sourceStation;
        int queueStation = topologyInfo.queueStation;

        for (int k = 0; k < K; k++) {
            if (sn.njobs.get(0, k) == Double.POSITIVE_INFINITY) {
                ClassResult classResult = analyzeOpenClass(sn, k, sourceStation, queueStation, options);
                result.QN.set(queueStation, k, classResult.meanQueueLength);
                result.RN.set(queueStation, k, classResult.meanResponseTime);
                result.TN.set(queueStation, k, classResult.throughput);
                result.UN.set(queueStation, k, classResult.utilization);
                result.TN.set(sourceStation, k, classResult.throughput);
            }
        }
        xvecIt = Matrix.zeros(1, 1);
    }

    @Override
    public Matrix getXVecIt() { return xvecIt; }

    /**
     * The reason MFQ cannot run on this model, or {@code null} when it can.
     *
     * MFQ IS A SINGLE-QUEUE METHOD AND FALLS BACK: `solver_fluid_analyzer.m`
     * warns "MFQ not applicable: ... Falling back to matrix method" and
     * re-enters `solver_fluid_matrix`. The decision has to be taken by the
     * CALLER, before the mfq branch of runAnalyzer diverts around the state-space
     * preparation the matrix method reads (options.init_sol above all), which is
     * why this is exposed rather than handled inside {@link #analyze}.
     *
     * @param sn the model
     * @return the topology error, or null when MFQ applies
     */
    public static String mfqNotApplicableReason(NetworkStruct sn) {
        TopologyInfo info = new MFQAnalyzer().validateTopology(sn);
        return info.isValid ? null : info.errorMsg;
    }

    private TopologyInfo validateTopology(NetworkStruct sn) {
        TopologyInfo info = new TopologyInfo();
        boolean hasOpenClass = false;
        for (int k = 0; k < sn.nclasses; k++) {
            if (sn.njobs.get(0, k) == Double.POSITIVE_INFINITY) { hasOpenClass = true; break; }
        }
        if (!hasOpenClass) { info.errorMsg = "Not an open model - all classes are closed"; return info; }

        int sourceNode = -1, queueNode = -1, sinkNode = -1;
        int sourceCount = 0, queueCount = 0, sinkCount = 0;
        for (int i = 0; i < sn.nnodes; i++) {
            NodeType nt = sn.nodetype.get(i);
            if (nt == NodeType.Source) { sourceNode = i; sourceCount++; }
            else if (nt == NodeType.Queue) { queueNode = i; queueCount++; }
            else if (nt == NodeType.Sink) { sinkNode = i; sinkCount++; }
        }

        if (sourceCount == 0) { info.errorMsg = "No source node found"; return info; }
        if (sourceCount > 1) { info.errorMsg = "Multiple source nodes found (" + sourceCount + ")"; return info; }
        if (sinkCount == 0) { info.errorMsg = "No sink node found"; return info; }
        if (sinkCount > 1) { info.errorMsg = "Multiple sink nodes found (" + sinkCount + ")"; return info; }
        if (queueCount == 0) { info.errorMsg = "No queue node found"; return info; }
        if (queueCount > 1) { info.errorMsg = "Multiple queue nodes found (" + queueCount + ") - MFQ supports single queue only"; return info; }

        int sourceStation = (int) sn.nodeToStation.get(sourceNode);
        int queueStation = (int) sn.nodeToStation.get(queueNode);
        double c = sn.nservers.get(queueStation, 0);
        if (c > 1 && c < Double.POSITIVE_INFINITY) {
            info.errorMsg = "Multi-server queue (c=" + c + ") not supported - MFQ requires c=1 or c=Inf";
            return info;
        }
        info.isValid = true;
        info.sourceNode = sourceNode;
        info.queueNode = queueNode;
        info.sinkNode = sinkNode;
        info.sourceStation = sourceStation;
        info.queueStation = queueStation;
        return info;
    }

    private ClassResult analyzeOpenClass(NetworkStruct sn, int classIdx, int sourceStation, int queueStation, SolverOptions options) {
        double lambda = sn.rates.get(sourceStation, classIdx);
        if (Double.isNaN(lambda) || lambda <= 0) {
            throw new RuntimeException("No valid arrival rate for class " + classIdx + " at source");
        }
        double mu = sn.rates.get(queueStation, classIdx);
        if (Double.isNaN(mu) || mu <= 0) {
            throw new RuntimeException("No valid service rate for class " + classIdx + " at queue");
        }
        double rho = lambda / mu;
        if (rho >= 1.0) {
            return new ClassResult(Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY, lambda, 1.0);
        }
        ProcessParams arrivalParams = extractProcess(sn, sourceStation, classIdx, lambda);
        ProcessParams serviceParams = extractProcess(sn, queueStation, classIdx, mu);
        boolean isSimpleExponential = arrivalParams.isSimple && serviceParams.isSimple;

        double[] flMoms;
        double[] stMoms;
        if (isSimpleExponential) {
            double meanL = rho / (1 - rho);
            double varL = rho / ((1 - rho) * (1 - rho));
            double meanW = 1.0 / (mu - lambda);
            double varW = 1.0 / ((mu - lambda) * (mu - lambda));
            flMoms = new double[]{meanL, meanL * meanL + varL};
            stMoms = new double[]{meanW, meanW * meanW + varW};
        } else {
            double prec = Math.max(options.tol, 1e-14);
            FluFluResult fluResult = FluFluQueue.fluFluQueue(arrivalParams.Q, arrivalParams.R,
                    serviceParams.Q, serviceParams.R, true, 2, 2, prec);
            flMoms = fluResult.getFluidMoments();
            stMoms = fluResult.getSojournMoments();
            if (flMoms == null) throw new RuntimeException("FluFluQueue did not return fluid moments");
            if (stMoms == null) throw new RuntimeException("FluFluQueue did not return sojourn moments");
        }
        double meanQueueLength = flMoms[0];
        double meanResponseTime = stMoms[0];
        double throughput = (meanResponseTime > GlobalConstants.FineTol) ? meanQueueLength / meanResponseTime : lambda;
        double utilization = Math.min(1.0, lambda / mu);
        return new ClassResult(meanQueueLength, meanResponseTime, throughput, utilization);
    }

    private ProcessParams extractProcess(NetworkStruct sn, int stationIdx, int classIdx, double rate) {
        jline.lang.nodes.Station station = sn.stations.get(stationIdx);
        jline.lang.JobClass jobClass = sn.jobclasses.get(classIdx);
        java.util.Map<jline.lang.JobClass, MatrixCell> stationProc = sn.proc.get(station);
        MatrixCell proc = (stationProc != null) ? stationProc.get(jobClass) : null;
        if (proc == null || proc.size() == 0) {
            return new ProcessParams(Matrix.zeros(1, 1), Matrix.singleton(rate), true);
        }
        int nPhases = proc.get(0).getNumRows();
        if (nPhases == 1) {
            return new ProcessParams(Matrix.zeros(1, 1), Matrix.singleton(rate), true);
        }
        Matrix D0 = proc.get(0);
        Matrix D1 = proc.get(1);
        Matrix Q = D0.add(1.0, D1);
        Matrix R = Matrix.zeros(nPhases, nPhases);
        for (int i = 0; i < nPhases; i++) {
            double rowSum = 0.0;
            for (int j = 0; j < D1.getNumCols(); j++) rowSum += D1.get(i, j);
            R.set(i, i, rowSum);
        }
        return new ProcessParams(Q, R, false);
    }

    private static class TopologyInfo {
        boolean isValid = false;
        String errorMsg = "";
        int sourceNode = -1;
        int queueNode = -1;
        int sinkNode = -1;
        int sourceStation = -1;
        int queueStation = -1;
    }

    private static class ProcessParams {
        final Matrix Q;
        final Matrix R;
        final boolean isSimple;
        ProcessParams(Matrix Q, Matrix R, boolean isSimple) {
            this.Q = Q; this.R = R; this.isSimple = isSimple;
        }
    }

    private static class ClassResult {
        final double meanQueueLength;
        final double meanResponseTime;
        final double throughput;
        final double utilization;
        ClassResult(double meanQueueLength, double meanResponseTime, double throughput, double utilization) {
            this.meanQueueLength = meanQueueLength;
            this.meanResponseTime = meanResponseTime;
            this.throughput = throughput;
            this.utilization = utilization;
        }
    }
}
