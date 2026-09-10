package jline.api.fes;

import jline.api.mc.Dtmc_stochcomp;
import jline.api.pfqn.ld.Ljd;
import jline.lang.ClosedClass;
import jline.lang.JobClass;
import jline.lang.Network;
import jline.lang.NetworkStruct;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Station;
import jline.lang.processes.Disabled;
import jline.lang.processes.Exp;
import jline.solvers.mva.SolverMVA;
import jline.util.Utils;
import jline.util.matrix.Matrix;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;

/**
 * Flow-Equivalent Server (FES) Aggregation (Internal Implementation).
 *
 * @deprecated Use jline.lang.ModelAdapter.aggregateFES() for public API
 */
@Deprecated
public final class FESAggregator {

    private FESAggregator() {
    }

    public static FESResult aggregateFES(Network model, List<Station> stationSubset) {
        return aggregateFES(model, stationSubset, FESOptions.defaults());
    }

    public static FESResult aggregateFES(Network model, List<Station> stationSubset, FESOptions options) {
        validateInputs(model, stationSubset);

        NetworkStruct sn = model.getStruct(true);
        int M = sn.nstations;
        int K = sn.nclasses;

        List<Station> modelStations = model.getStations();
        int[] subsetIndices = new int[stationSubset.size()];
        for (int i = 0; i < stationSubset.size(); i++) {
            for (int j = 0; j < M; j++) {
                if (stationSubset.get(i) == modelStations.get(j)) {
                    subsetIndices[i] = j;
                    break;
                }
            }
        }

        List<Integer> subsetList = new ArrayList<Integer>();
        for (int x : subsetIndices) subsetList.add(x);
        List<Integer> complementList = new ArrayList<Integer>();
        for (int i = 0; i < M; i++) {
            if (!subsetList.contains(i)) complementList.add(i);
        }
        int[] complementIndices = new int[complementList.size()];
        for (int i = 0; i < complementList.size(); i++) complementIndices[i] = complementList.get(i);

        Matrix cutoffs;
        if (options.getCutoffs() != null) {
            cutoffs = options.getCutoffs();
        } else {
            cutoffs = new Matrix(1, K);
            for (int k = 0; k < K; k++) {
                cutoffs.set(0, k, sn.njobs.get(k));
            }
        }

        Matrix rt = sn.rt;

        List<Integer> subsetRtIndices = new ArrayList<Integer>();
        for (int i : subsetIndices) {
            int isf = (int) sn.stationToStateful.get(i);
            for (int k = 0; k < K; k++) {
                subsetRtIndices.add(isf * K + k);
            }
        }

        List<Integer> complementRtIndices = new ArrayList<Integer>();
        for (int i : complementIndices) {
            int isf = (int) sn.stationToStateful.get(i);
            for (int k = 0; k < K; k++) {
                complementRtIndices.add(isf * K + k);
            }
        }

        Matrix stochCompSubset = Dtmc_stochcomp.dtmc_stochcomp(rt, new ArrayList<Integer>(subsetRtIndices));
        Matrix stochCompComplement = Dtmc_stochcomp.dtmc_stochcomp(rt, new ArrayList<Integer>(complementRtIndices));

        Network isolatedModel = buildIsolatedModel(model, stationSubset, stochCompSubset, sn);

        // see _kb/03-api-layer.md for rationale
        double[][] pexit = new double[K][subsetIndices.length];
        for (int k = 0; k < K; k++) {
            for (int a = 0; a < subsetIndices.length; a++) {
                int isf_j = (int) sn.stationToStateful.get(subsetIndices[a]);
                double s = 0.0;
                for (int i : complementIndices) {
                    int isf_i = (int) sn.stationToStateful.get(i);
                    int ri = isf_j * K + k, ci = isf_i * K + k;
                    if (ri < rt.getNumRows() && ci < rt.getNumCols()) s += rt.get(ri, ci);
                }
                pexit[k][a] = s;
            }
        }

        List<Matrix> throughputTable = computeThroughputs(isolatedModel, cutoffs, options, pexit);

        Network fesModel = new Network(model.getName() + "_FES");

        HashMap<Integer, Station> stationMap = new HashMap<Integer, Station>();
        for (int i : complementIndices) {
            Station origStation = modelStations.get(i);
            Station newStation;
            if (origStation instanceof Queue) {
                Queue q = new Queue(fesModel, origStation.getName(), ((Queue) origStation).getSchedStrategy());
                if (!Utils.isInf((double) ((Queue) origStation).getNumberOfServers())) {
                    q.setNumberOfServers(((Queue) origStation).getNumberOfServers());
                }
                newStation = q;
            } else if (origStation instanceof Delay) {
                newStation = new Delay(fesModel, origStation.getName());
            } else {
                throw new IllegalArgumentException("Unsupported station type: " + origStation.getClass());
            }
            stationMap.put(i, newStation);
        }

        Queue fesStation = new Queue(fesModel, "FES", SchedStrategy.PS);
        fesStation.setNumberOfServers(1);

        Station refStation;
        if (complementIndices.length > 0) {
            refStation = stationMap.get(complementIndices[0]);
        } else {
            refStation = fesStation;
        }

        ArrayList<ClosedClass> newClasses = new ArrayList<ClosedClass>();
        for (int k = 0; k < K; k++) {
            JobClass origClass = model.getClasses().get(k);
            int population = (int) sn.njobs.get(k);
            ClosedClass newClass = new ClosedClass(fesModel, origClass.getName(), population, refStation, 0);
            newClasses.add(newClass);
        }

        for (int i : complementIndices) {
            Station newStation = stationMap.get(i);
            for (int k = 0; k < K; k++) {
                double rate = sn.rates.get(i, k);
                if (rate > 0 && !Double.isNaN(rate)) {
                    newStation.setService(newClasses.get(k), new Exp(rate));
                } else {
                    newStation.setService(newClasses.get(k), new Disabled());
                }
            }
        }

        for (int k = 0; k < K; k++) {
            fesStation.setService(newClasses.get(k), new Exp(1.0));
        }

        int tableSize = computeTableSize(cutoffs);
        List<Matrix> scalingTables = new ArrayList<Matrix>();
        for (int k = 0; k < K; k++) {
            Matrix stTable = new Matrix(1, tableSize);
            for (int idx = 0; idx < tableSize; idx++) {
                double tput = throughputTable.get(k).get(0, idx);
                // see _kb/03-api-layer.md for rationale
                stTable.set(0, idx, (tput > 1e-10) ? tput : 1e-10);
            }
            scalingTables.add(stTable);
        }

        // see _kb/03-api-layer.md for rationale
        FesBetaFunction fesBeta = new FesBetaFunction(scalingTables, cutoffs);
        // see _kb/03-api-layer.md for rationale
        int Kfes = (int) cutoffs.getNumElements();
        int[] fesNK = new int[Kfes];
        for (int r = 0; r < Kfes; r++) fesNK[r] = (int) Math.round(cutoffs.get(r));
        double fesPeak = jline.api.pfqn.ld.CdPeakScaling.cd_peak_scaling(fesBeta, fesNK, Kfes);
        Matrix fesPeakVec = new Matrix(1, 1);
        fesPeakVec.set(0, 0, fesPeak);
        fesStation.setClassDependence(fesBeta, fesPeakVec);

        List<jline.lang.nodes.Node> nodes = fesModel.getNodes();
        int nNodes = nodes.size();

        int fesNodeIdx = 0;
        for (int n = 0; n < nNodes; n++) {
            if (nodes.get(n) == fesStation) {
                fesNodeIdx = n;
                break;
            }
        }

        HashMap<Integer, Integer> complementNodeMap = new HashMap<Integer, Integer>();
        for (int i : complementIndices) {
            Station station = stationMap.get(i);
            for (int n = 0; n < nNodes; n++) {
                if (nodes.get(n) == station) {
                    complementNodeMap.put(i, n);
                    break;
                }
            }
        }

        Object P = fesModel.initRoutingMatrix();

        for (int k = 0; k < K; k++) {
            Matrix Pk = new Matrix(nNodes, nNodes);

            for (int i : complementIndices) {
                if (!complementNodeMap.containsKey(i)) continue;
                int iNode = complementNodeMap.get(i);
                int isf_i = (int) sn.stationToStateful.get(i);

                for (int j : complementIndices) {
                    if (!complementNodeMap.containsKey(j)) continue;
                    int jNode = complementNodeMap.get(j);

                    // see _kb/03-api-layer.md for rationale
                    int isf_jc = (int) sn.stationToStateful.get(j);
                    int rtIdx_i = isf_i * K + k;
                    int rtIdx_j = isf_jc * K + k;

                    if (rtIdx_i < rt.getNumRows() && rtIdx_j < rt.getNumCols()) {
                        double prob = rt.get(rtIdx_i, rtIdx_j);
                        if (prob > 1e-10) {
                            Pk.set(iNode, jNode, prob);
                        }
                    }
                }

                for (int j : subsetIndices) {
                    int isf_j = (int) sn.stationToStateful.get(j);
                    int rtIdx_i = isf_i * K + k;
                    int rtIdx_j = isf_j * K + k;

                    if (rtIdx_i < rt.getNumRows() && rtIdx_j < rt.getNumCols()) {
                        double prob = rt.get(rtIdx_i, rtIdx_j);
                        if (prob > 1e-10) {
                            Pk.set(iNode, fesNodeIdx, Pk.get(iNode, fesNodeIdx) + prob);
                        }
                    }
                }
            }

            // see _kb/03-api-layer.md for rationale
            double[] visitRatios = classVisitRatios(stochCompSubset, subsetIndices.length, K, k);

            for (int j : complementIndices) {
                if (!complementNodeMap.containsKey(j)) continue;
                int jNode = complementNodeMap.get(j);
                int isf_j = (int) sn.stationToStateful.get(j);

                double probSum = 0.0;
                for (int p = 0; p < subsetIndices.length; p++) {
                    int isf_i = (int) sn.stationToStateful.get(subsetIndices[p]);
                    int rtIdx_i = isf_i * K + k;
                    int rtIdx_j = isf_j * K + k;

                    if (rtIdx_i < rt.getNumRows() && rtIdx_j < rt.getNumCols()) {
                        probSum += visitRatios[p] * rt.get(rtIdx_i, rtIdx_j);
                    }
                }

                if (probSum > 1e-10) {
                    Pk.set(fesNodeIdx, jNode, probSum);
                }
            }

            for (int n = 0; n < nNodes; n++) {
                double rowSum = 0.0;
                for (int m = 0; m < nNodes; m++) {
                    rowSum += Pk.get(n, m);
                }
                if (rowSum > 1e-10) {
                    for (int m = 0; m < nNodes; m++) {
                        Pk.set(n, m, Pk.get(n, m) / rowSum);
                    }
                }
            }

            ((jline.lang.RoutingMatrix) P).set(newClasses.get(k), Pk);
        }

        fesModel.link((jline.lang.RoutingMatrix) P);

        FESDeaggInfo deaggInfo = new FESDeaggInfo(
                model, stationSubset, subsetIndices, complementIndices,
                throughputTable, cutoffs, stochCompSubset, stochCompComplement,
                isolatedModel, fesNodeIdx);

        return new FESResult(fesModel, fesStation, deaggInfo);
    }

    /**
     * Per-class visit ratios of the subset stations: the stationary vector
     * (pi = pi * P) of the class-k routing among subset stations extracted from
     * the stochastic complement. Mirrors fes_build_isolated.m. Scale is
     * irrelevant since the FES out-routing row is normalized by the caller.
     */
    private static double[] classVisitRatios(Matrix stochCompS, int M_sub, int K, int k) {
        double[][] Pk = new double[M_sub][M_sub];
        for (int i = 0; i < M_sub; i++) {
            double rowSum = 0.0;
            for (int j = 0; j < M_sub; j++) {
                int r = i * K + k, c = j * K + k;
                double p = (r < stochCompS.getNumRows() && c < stochCompS.getNumCols())
                        ? stochCompS.get(r, c) : 0.0;
                if (p < 0) p = 0.0;
                Pk[i][j] = p;
                rowSum += p;
            }
            if (rowSum > 1e-10) {
                for (int j = 0; j < M_sub; j++) Pk[i][j] /= rowSum;
            } else {
                Pk[i][i] = 1.0; // no outgoing routing: route to self
            }
        }
        // see _kb/03-api-layer.md for rationale
        double[] pi = new double[M_sub];
        java.util.Arrays.fill(pi, 1.0 / M_sub);
        double[] next = new double[M_sub];
        for (int iter = 0; iter < 10000; iter++) {
            double sum = 0.0, diff = 0.0;
            for (int j = 0; j < M_sub; j++) {
                double acc = 0.0;
                for (int i = 0; i < M_sub; i++) acc += pi[i] * Pk[i][j];
                next[j] = 0.5 * pi[j] + 0.5 * acc;
                sum += next[j];
            }
            for (int j = 0; j < M_sub; j++) {
                double v = (sum > 0) ? next[j] / sum : 1.0 / M_sub;
                diff += Math.abs(v - pi[j]);
                pi[j] = v;
            }
            if (diff < 1e-15) break;
        }
        return pi;
    }

    private static void validateInputs(Network model, List<Station> stationSubset) {
        if (stationSubset.isEmpty()) {
            throw new IllegalArgumentException("Station subset cannot be empty");
        }

        NetworkStruct sn = model.getStruct(true);
        for (int k = 0; k < sn.nclasses; k++) {
            if (Utils.isInf(sn.njobs.get(k))) {
                throw new IllegalArgumentException("FES aggregation only applies to closed queueing networks");
            }
        }

        List<Station> modelStations = model.getStations();
        if (stationSubset.size() >= modelStations.size()) {
            throw new IllegalArgumentException("Cannot aggregate all stations");
        }

        for (Station station : stationSubset) {
            if (!(station instanceof Queue) && !(station instanceof Delay)) {
                throw new IllegalArgumentException("FES aggregation only supports Queue and Delay stations");
            }

            boolean found = false;
            for (Station ms : modelStations) {
                if (station == ms) {
                    found = true;
                    break;
                }
            }
            if (!found) {
                throw new IllegalArgumentException("Station " + station.getName() + " does not belong to the model");
            }
        }
    }

    private static Network buildIsolatedModel(Network model, List<Station> stationSubset,
                                                Matrix stochCompS, NetworkStruct sn) {
        Network isolatedModel = new Network(model.getName() + "_isolated");
        int K = sn.nclasses;
        int M_sub = stationSubset.size();

        List<Station> modelStations = model.getStations();
        int[] subsetIndices = new int[M_sub];
        for (int i = 0; i < M_sub; i++) {
            for (int j = 0; j < sn.nstations; j++) {
                if (stationSubset.get(i) == modelStations.get(j)) {
                    subsetIndices[i] = j;
                    break;
                }
            }
        }

        ArrayList<Station> newStations = new ArrayList<Station>();
        for (int i = 0; i < M_sub; i++) {
            Station origStation = stationSubset.get(i);
            Station newStation;
            if (origStation instanceof Queue) {
                Queue q = new Queue(isolatedModel, origStation.getName(), ((Queue) origStation).getSchedStrategy());
                if (!Utils.isInf((double) ((Queue) origStation).getNumberOfServers())) {
                    q.setNumberOfServers(((Queue) origStation).getNumberOfServers());
                }
                newStation = q;
            } else if (origStation instanceof Delay) {
                newStation = new Delay(isolatedModel, origStation.getName());
            } else {
                throw new IllegalArgumentException("Unsupported station type");
            }
            newStations.add(newStation);
        }

        Station refStation = newStations.get(0);
        ArrayList<ClosedClass> newClasses = new ArrayList<ClosedClass>();
        for (int k = 0; k < K; k++) {
            JobClass origClass = model.getClasses().get(k);
            ClosedClass newClass = new ClosedClass(isolatedModel, origClass.getName(), 1, refStation, 0);
            newClasses.add(newClass);
        }

        for (int i = 0; i < M_sub; i++) {
            int origIdx = subsetIndices[i];
            Station newStation = newStations.get(i);

            for (int k = 0; k < K; k++) {
                double rate = sn.rates.get(origIdx, k);
                if (rate > 0 && !Double.isNaN(rate)) {
                    newStation.setService(newClasses.get(k), new Exp(rate));
                } else {
                    newStation.setService(newClasses.get(k), new Disabled());
                }
            }
        }

        Object P = isolatedModel.initRoutingMatrix();
        List<jline.lang.nodes.Node> nodes = isolatedModel.getNodes();
        int nNodes = nodes.size();

        int[] stationToNode = new int[M_sub];
        for (int i = 0; i < M_sub; i++) {
            for (int n = 0; n < nNodes; n++) {
                if (nodes.get(n) == newStations.get(i)) {
                    stationToNode[i] = n;
                    break;
                }
            }
        }

        for (int k = 0; k < K; k++) {
            Matrix Pk = new Matrix(nNodes, nNodes);

            for (int i = 0; i < M_sub; i++) {
                int rowIdx = i * K + k;
                int iNode = stationToNode[i];

                for (int j = 0; j < M_sub; j++) {
                    int colIdx = j * K + k;
                    int jNode = stationToNode[j];

                    if (rowIdx < stochCompS.getNumRows() && colIdx < stochCompS.getNumCols()) {
                        double prob = stochCompS.get(rowIdx, colIdx);
                        if (prob > 1e-10) {
                            Pk.set(iNode, jNode, prob);
                        }
                    }
                }

                double rowSum = 0.0;
                for (int n = 0; n < nNodes; n++) {
                    rowSum += Pk.get(iNode, n);
                }
                if (rowSum > 1e-10) {
                    for (int n = 0; n < nNodes; n++) {
                        Pk.set(iNode, n, Pk.get(iNode, n) / rowSum);
                    }
                } else {
                    Pk.set(iNode, iNode, 1.0);
                }
            }

            ((jline.lang.RoutingMatrix) P).set(newClasses.get(k), Pk);
        }

        isolatedModel.link((jline.lang.RoutingMatrix) P);
        return isolatedModel;
    }

    private static List<Matrix> computeThroughputs(Network isolatedModel, Matrix cutoffs, FESOptions options, double[][] pexit) {
        int K = isolatedModel.getClasses().size();
        int tableSize = computeTableSize(cutoffs);

        ArrayList<Matrix> throughputTable = new ArrayList<Matrix>();
        for (int k = 0; k < K; k++) {
            throughputTable.add(new Matrix(1, tableSize));
        }

        HashSet<String> visited = new HashSet<String>();
        LinkedList<int[]> queue = new LinkedList<int[]>();
        queue.add(new int[K]);

        while (!queue.isEmpty()) {
            int[] nvec = queue.poll();
            String stateKey = Arrays.toString(nvec);

            if (visited.contains(stateKey)) continue;
            visited.add(stateKey);

            int totalPop = 0;
            for (int v : nvec) totalPop += v;

            Matrix nvecMatrix = new Matrix(1, K);
            for (int i = 0; i < K; i++) nvecMatrix.set(0, i, (double) nvec[i]);
            int idx = Ljd.ljd_linearize(nvecMatrix, cutoffs);

            if (totalPop == 0) {
                for (int k = 0; k < K; k++) {
                    throughputTable.get(k).set(0, idx, 0.0);
                }
            } else {
                try {
                    for (int k = 0; k < K; k++) {
                        ((ClosedClass) isolatedModel.getClasses().get(k)).setPopulation((double) nvec[k]);
                    }
                    // see _kb/03-api-layer.md for rationale
                    isolatedModel.resetStruct();

                    SolverMVA solver = new SolverMVA(isolatedModel);
                    solver.runAnalyzer();
                    Object result = solver.result;

                    // Departure rate per class = sum_j T_{j,r} * P_r(j->complement),
                    // using per-station throughputs TN of the isolated subnetwork.
                    Matrix TN = (Matrix) result.getClass().getField("TN").get(result);
                    for (int k = 0; k < K; k++) {
                        if (nvec[k] > 0 && TN != null) {
                            double dep = 0.0;
                            for (int a = 0; a < pexit[k].length; a++) {
                                dep += TN.get(a, k) * pexit[k][a];
                            }
                            throughputTable.get(k).set(0, idx, dep);
                        } else {
                            throughputTable.get(k).set(0, idx, 0.0);
                        }
                    }
                } catch (Exception e) {
                    for (int k = 0; k < K; k++) {
                        throughputTable.get(k).set(0, idx, 0.0);
                    }
                }
            }

            for (int k = 0; k < K; k++) {
                if (nvec[k] < (int) cutoffs.get(0, k)) {
                    int[] newState = nvec.clone();
                    newState[k] = newState[k] + 1;
                    String newKey = Arrays.toString(newState);
                    if (!visited.contains(newKey)) {
                        queue.add(newState);
                    }
                }
            }
        }

        return throughputTable;
    }

    private static int computeTableSize(Matrix cutoffs) {
        int size = 1;
        for (int k = 0; k < cutoffs.getNumCols(); k++) {
            size *= ((int) cutoffs.get(0, k) + 1);
        }
        return size;
    }

    /** Per-population queue lengths and utilizations of the isolated subnetwork. */
    public static final class ConditionalMetrics {
        /** Indexed by the linearized population state; each (M_sub x K). */
        public final List<Matrix> QN;
        /** Indexed by the linearized population state; each (M_sub x K). */
        public final List<Matrix> UN;

        ConditionalMetrics(List<Matrix> QN, List<Matrix> UN) {
            this.QN = QN;
            this.UN = UN;
        }
    }

    /**
     * Per-station metrics of the ISOLATED subnetwork at every population state.
     *
     * <p>The companion of {@code computeThroughputs}. That returns the aggregate
     * throughput X(n) that becomes the flow-equivalent server's rate, which is
     * all the REDUCED model needs; it is not enough to report the COLLAPSED
     * stations' own metrics. Those are recovered by conditioning on the FES
     * population, E[Q_i] = sum_n P(N_fes = n) * Q_i(n) -- the Chandy-Herzog-Woo
     * hierarchical decomposition, exact when the subnetwork is product-form.
     * This supplies the Q_i(n) and U_i(n) that sum is taken over, indexed by
     * {@code Ljd.ljd_linearize} exactly as the throughput table is.</p>
     *
     * @param isolatedModel the isolated subnetwork built by the transform
     * @param cutoffs per-class population cutoffs
     * @return the per-population tables
     */
    public static ConditionalMetrics computeConditionalMetrics(Network isolatedModel,
                                                               Matrix cutoffs) {
        int K = isolatedModel.getClasses().size();
        int M = isolatedModel.getStations().size();
        int tableSize = computeTableSize(cutoffs);
        ArrayList<Matrix> QN = new ArrayList<Matrix>();
        ArrayList<Matrix> UN = new ArrayList<Matrix>();
        for (int idx = 0; idx < tableSize; idx++) {
            QN.add(new Matrix(M, K));
            UN.add(new Matrix(M, K));
        }
        for (int idx = 0; idx < tableSize; idx++) {
            Matrix nvec = jline.api.pfqn.ld.Ljd.ljd_delinearize(idx, cutoffs);
            double total = 0;
            for (int k = 0; k < K; k++) total += nvec.get(0, k);
            if (total == 0) {
                continue; // an empty subnetwork holds nothing and serves nothing
            }
            try {
                for (int k = 0; k < K; k++) {
                    ((ClosedClass) isolatedModel.getClasses().get(k))
                            .setPopulation(nvec.get(0, k));
                }
                isolatedModel.resetStruct();
                SolverMVA solver = new SolverMVA(isolatedModel);
                solver.runAnalyzer();
                Object result = solver.result;
                Matrix Q = (Matrix) result.getClass().getField("QN").get(result);
                Matrix U = (Matrix) result.getClass().getField("UN").get(result);
                if (Q != null) QN.set(idx, Q.copy());
                if (U != null) UN.set(idx, U.copy());
            } catch (Exception e) {
                // A population the isolated subnetwork cannot hold contributes
                // nothing; the reduced chain gives it probability zero anyway.
            }
        }
        return new ConditionalMetrics(QN, UN);
    }
}
