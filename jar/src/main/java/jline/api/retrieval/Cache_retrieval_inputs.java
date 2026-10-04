/**
 * @file Cache_retrieval_inputs.java
 * @brief Extract delayed-hit retrieval algorithm inputs from a NetworkStruct.
 *
 * Port of matlab/src/api/retrieval/cache_retrieval_inputs.m. Builds the inputs of
 * the retrieval_* algorithms (m, lambda, gamma, eta, alpha, T, R, station_type) for a
 * cache equipped with a retrieval system (Cache.setRetrievalSystem). Single read
 * class (IRM). Supported station types: IS, PS, SIRO, FCFS, LCFSPR; SIRO/FCFS
 * require exponential service with identical per-class rates (LCFSPR is exempt).
 *
 * @since LINE 3.0
 */
package jline.api.retrieval;

import java.util.List;
import java.util.Set;

import jline.api.cache.Cache_gamma_lp;
import jline.io.Ret;
import jline.lang.NetworkStruct;
import jline.lang.JobClass;
import jline.lang.nodes.Station;
import jline.lang.constant.NodeType;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodeparam.CacheNodeParam;
import jline.lang.nodes.Cache;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Cache_retrieval_inputs {
    private Cache_retrieval_inputs() {}

    public static final class Inputs {
        public double[] m;
        public double[] lambda;
        public double[][] gamma;        // n x h
        public double[][] eta;          // n x (r+1)
        public double[][][] alpha;      // [s][i][phase]
        public double[][][][] T;        // [s][i][phase][phase]
        public double[][][] R;          // [i][S+1][S+1]
        public String[] stationType;    // [S]
        public int jobinClass;          // 0-indexed read class
        public int[] queueNode;         // [S] 0-indexed node indices
        public int[] queueStation;      // [S] 0-indexed station indices
        public double sourceRate;       // arrival rate of the read class
    }

    public static Inputs cache_retrieval_inputs(NetworkStruct sn) {
        return cache_retrieval_inputs(sn, Double.NaN);
    }

    /**
     * Builds the retrieval-algorithm inputs, optionally overriding the read-class
     * arrival rate with lambdaOverride (used for a CLOSED integrated cache-queueing
     * sublayer where the read rate comes from the network solution rather than a
     * Source throughput). Pass Double.NaN to use the Source throughput as before.
     */
    public static Inputs cache_retrieval_inputs(NetworkStruct sn, double lambdaOverride) {
        int K = sn.nclasses;

        int ci = -1;
        for (int i = 0; i < sn.nodetype.size(); i++) if (sn.nodetype.get(i) == NodeType.Cache) { ci = i; break; }
        if (ci < 0) throw new RuntimeException("Retrieval analysis requires a Cache node.");
        Cache cache = (Cache) sn.nodes.get(ci);
        CacheNodeParam ch = (CacheNodeParam) sn.nodeparam.get(cache);
        if (ch.retrievalSystemCapacity <= 0) throw new RuntimeException("The Cache node has no retrieval system.");

        Inputs in = new Inputs();
        int h = ch.itemcap.length();
        in.m = new double[h];
        for (int j = 0; j < h; j++) in.m[j] = ch.itemcap.get(j);
        int n = ch.nitems;

        // read class (single-class IRM)
        Set<Integer> keys = ch.retrievalSystemQueueIndices.keySet();
        if (keys.size() != 1) throw new RuntimeException("Retrieval analysis supports a single read class.");
        int jobinClass = keys.iterator().next();
        in.jobinClass = jobinClass;
        List<Integer> qn = ch.retrievalSystemQueueIndices.get(jobinClass);
        int S = qn.size();
        in.queueNode = new int[S];
        in.queueStation = new int[S];
        for (int s = 0; s < S; s++) {
            in.queueNode[s] = qn.get(s);
            in.queueStation[s] = (int) sn.nodeToStation.get(in.queueNode[s]);
        }

        // see _kb/03-api-layer.md for rationale (### retrieval)
        double sourceRate;
        if (!Double.isNaN(lambdaOverride)) {
            sourceRate = lambdaOverride;
        } else {
            int source_ist = -1;
            for (int i = 0; i < sn.nodetype.size(); i++) if (sn.nodetype.get(i) == NodeType.Source) { source_ist = (int) sn.nodeToStation.get(i); break; }
            if (source_ist < 0) throw new RuntimeException("Retrieval analysis of a closed model requires an explicit read rate (cache_retrieval_inputs(sn, lambdaOverride)); no Source node found.");
            sourceRate = sn.rates.get(source_ist, jobinClass);
        }
        if (Double.isNaN(sourceRate)) sourceRate = 0;
        in.sourceRate = sourceRate;
        List<Double> pread = ch.pread.get(jobinClass);
        in.lambda = new double[n];
        for (int i = 0; i < n; i++) in.lambda[i] = sourceRate * pread.get(i);

        // gamma via Cache_gamma_lp (single read class)
        Matrix[] lambda3d = new Matrix[1];
        lambda3d[0] = new Matrix(n, h + 1);
        for (int i = 0; i < n; i++) for (int l = 0; l < h + 1; l++) lambda3d[0].set(i, l, in.lambda[i]);
        Matrix[][] Rcost = ch.accost;
        if (Rcost == null || Rcost.length == 0) {
            Rcost = new Matrix[1][n];
            for (int k = 0; k < n; k++) {
                Matrix Rmat = new Matrix(h + 1, h + 1);
                for (int j = 0; j < h; j++) Rmat.set(j, j + 1, 1);
                Rmat.set(h, h, 1);
                Rcost[0][k] = Rmat;
            }
        }
        Ret.cacheGamma cg = Cache_gamma_lp.cache_gamma_lp(lambda3d, Rcost);
        in.gamma = new double[n][h];
        for (int i = 0; i < n; i++) for (int j = 0; j < h; j++) in.gamma[i][j] = cg.gamma.get(i, j);

        // station types
        in.stationType = new String[S];
        for (int s = 0; s < S; s++) {
            SchedStrategy sd = sn.sched.get(sn.stations.get(in.queueStation[s]));
            if (sd == SchedStrategy.INF) in.stationType[s] = "IS";
            else if (sd == SchedStrategy.PS) in.stationType[s] = "PS";
            else if (sd == SchedStrategy.SIRO) in.stationType[s] = "SIRO";
            else if (sd == SchedStrategy.FCFS) in.stationType[s] = "FCFS";
            else if (sd == SchedStrategy.LCFSPR) in.stationType[s] = "LCFSPR";
            else throw new RuntimeException("Retrieval analysis supports only IS, PS, SIRO, FCFS, LCFSPR stations; node " + in.queueNode[s]);
        }

        // per-item PH service (alpha, T) and routing R
        in.alpha = new double[S][n][];
        in.T = new double[S][n][][];
        int[] fsz = new int[S];
        for (int s = 0; s < S; s++) {
            Station st = sn.stations.get(in.queueStation[s]);
            int rcls0 = (int) ch.retrievalClasses.get(0, jobinClass);
            MatrixCell pc0 = sn.proc.get(st).get(sn.jobclasses.get(rcls0));
            fsz[s] = pc0.get(0).getNumRows();
        }
        in.R = new double[n][S + 1][S + 1];
        for (int i = 0; i < n; i++) {
            int rcls = (int) ch.retrievalClasses.get(i, jobinClass);
            JobClass rc = sn.jobclasses.get(rcls);
            for (int s = 0; s < S; s++) {
                Station st = sn.stations.get(in.queueStation[s]);
                Matrix pieRow = sn.pie.get(st).get(rc);
                Matrix D0 = sn.proc.get(st).get(rc).get(0);
                in.alpha[s][i] = new double[fsz[s]];
                for (int b = 0; b < fsz[s]; b++) in.alpha[s][i][b] = pieRow.get(b);
                in.T[s][i] = new double[fsz[s]][fsz[s]];
                for (int a = 0; a < fsz[s]; a++) for (int b = 0; b < fsz[s]; b++) in.T[s][i][a][b] = D0.get(a, b);
            }
            // routing (index 0 = cache/outside, 1..S = queues); rtnodes flat = node*K + class
            for (int s = 0; s < S; s++) {
                in.R[i][0][s + 1] = sn.rtnodes.get(ci * K + rcls, in.queueNode[s] * K + rcls);
                in.R[i][s + 1][0] = sn.rtnodes.get(in.queueNode[s] * K + rcls, ci * K + rcls);
                for (int sp = 0; sp < S; sp++)
                    in.R[i][s + 1][sp + 1] = sn.rtnodes.get(in.queueNode[s] * K + rcls, in.queueNode[sp] * K + rcls);
            }
        }

        // SIRO/FCFS: exponential single-phase + identical rates
        for (int s = 0; s < S; s++) {
            if (in.stationType[s].equals("SIRO") || in.stationType[s].equals("FCFS")) {
                if (fsz[s] > 1) throw new RuntimeException("SIRO/FCFS retrieval stations require exponential (single-phase) service; node " + in.queueNode[s]);
                double tau0 = phMean(in.alpha[s][0], in.T[s][0]);
                for (int i = 1; i < n; i++)
                    if (Math.abs(phMean(in.alpha[s][i], in.T[s][i]) - tau0) > 1e-9 * Math.max(tau0, 1e-300))
                        throw new RuntimeException("SIRO/FCFS retrieval stations require identical per-class rates; node " + in.queueNode[s]);
            }
        }

        // eta(i,0) = sum IS visits*mean; eta(i,1+p) = PS-type station p
        int r = 0; int[] psIdx = new int[S]; boolean[] isIS = new boolean[S];
        for (int s = 0; s < S; s++) { isIS[s] = in.stationType[s].equals("IS"); if (!isIS[s]) psIdx[r++] = s; }
        in.eta = new double[n][r + 1];
        for (int i = 0; i < n; i++) {
            double[] visits = visits(in.R[i], S);
            double[] tau = new double[S];
            for (int s = 0; s < S; s++) tau[s] = phMean(in.alpha[s][i], in.T[s][i]);
            for (int s = 0; s < S; s++) if (isIS[s]) in.eta[i][0] += visits[s] * tau[s];
            for (int p = 0; p < r; p++) in.eta[i][1 + p] = visits[psIdx[p]] * tau[psIdx[p]];
        }
        return in;
    }

    private static double phMean(double[] al, double[][] Tm) {
        int k = al.length;
        double[][] negT = new double[k][k];
        for (int a = 0; a < k; a++) for (int b = 0; b < k; b++) negT[a][b] = -Tm[a][b];
        double[] e = new double[k]; java.util.Arrays.fill(e, 1.0);
        double[] z = solve(negT, e);
        double mu = 0; for (int a = 0; a < k; a++) mu += al[a] * z[a];
        return mu;
    }

    private static double[] visits(double[][] R, int S) {
        double[] a = new double[S];
        double[][] At = new double[S][S];
        for (int s = 0; s < S; s++) {
            a[s] = R[0][s + 1];
            for (int sp = 0; sp < S; sp++) At[sp][s] = (s == sp ? 1.0 : 0.0) - R[s + 1][sp + 1];
        }
        return solve(At, a);
    }

    private static double[] solve(double[][] Ain, double[] bin) {
        int N = bin.length;
        double[][] A = new double[N][N];
        double[] b = bin.clone();
        for (int i = 0; i < N; i++) A[i] = Ain[i].clone();
        for (int col = 0; col < N; col++) {
            int piv = col; double best = Math.abs(A[col][col]);
            for (int rrow = col + 1; rrow < N; rrow++) if (Math.abs(A[rrow][col]) > best) { best = Math.abs(A[rrow][col]); piv = rrow; }
            if (piv != col) { double[] tr = A[piv]; A[piv] = A[col]; A[col] = tr; double tb = b[piv]; b[piv] = b[col]; b[col] = tb; }
            double d = A[col][col];
            for (int rrow = col + 1; rrow < N; rrow++) {
                double f = A[rrow][col] / d; if (f == 0) continue;
                for (int c = col; c < N; c++) A[rrow][c] -= f * A[col][c];
                b[rrow] -= f * b[col];
            }
        }
        double[] x = new double[N];
        for (int rrow = N - 1; rrow >= 0; rrow--) {
            double sum = b[rrow];
            for (int c = rrow + 1; c < N; c++) sum -= A[rrow][c] * x[c];
            x[rrow] = sum / A[rrow][rrow];
        }
        return x;
    }
}
