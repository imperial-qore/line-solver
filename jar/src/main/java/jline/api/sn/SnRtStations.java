/**
 * @file Station-space projection of the routing table and of the visits
 *
 * sn.rt and sn.visits are indexed by stateful node. Solvers that write traffic
 * equations over stations need them indexed by station, which this class
 * produces by absorbing the stateful nodes that are not stations.
 *
 * @since LINE 3.0
 */
package jline.api.sn;

import jline.lang.NetworkStruct;
import jline.util.matrix.Matrix;
import jline.util.Pair;

public final class SnRtStations {
    private SnRtStations() {}

    /**
     * Station-to-station routing probabilities and per-station visits.
     *
     * <p>Indexing sn.rt or sn.visits by station index silently reads the wrong
     * rows as soon as the model owns a stateful node that is not a station
     * (Router, Cache, stateful class switch). The routing matrix returned here
     * absorbs those nodes,
     *
     * <pre>Pst = P_AA + P_AB * (I - P_BB)^-1 * P_BA,</pre>
     *
     * with A the station rows in station order and B the remaining stateful
     * rows, which is exact because a non-station stateful node holds no jobs:
     * it passes every arrival on instantaneously. When every stateful node is a
     * station the result is sn.rt unchanged.</p>
     *
     * @param sn network structure
     * @return the (M*K,M*K) routing matrix and the (M,K) visits
     */
    public static Pair<Matrix, Matrix> snRtStations(NetworkStruct sn) {
        int K = sn.nclasses;
        int M = sn.nstations;
        int S = sn.nstateful;

        boolean[] isSt = new boolean[S];
        for (int ist = 0; ist < M; ist++) {
            isSt[(int) sn.stationToStateful.get(ist)] = true;
        }
        int[] A = new int[M * K];
        for (int ist = 0; ist < M; ist++) {
            int isf = (int) sn.stationToStateful.get(ist);
            for (int r = 0; r < K; r++) {
                A[ist * K + r] = isf * K + r;
            }
        }
        int nB = (S - M) * K;
        int[] B = new int[nB];
        int pos = 0;
        for (int isf = 0; isf < S; isf++) {
            if (!isSt[isf]) {
                for (int r = 0; r < K; r++) {
                    B[pos++] = isf * K + r;
                }
            }
        }

        Matrix P = sn.rt;
        Matrix rtst = submatrix(P, A, A);
        if (nB > 0) {
            Matrix Pbb = submatrix(P, B, B);
            Matrix ImPbb = Matrix.eye(nB).add(-1.0, Pbb);
            Matrix X = ImPbb.inv().mult(submatrix(P, B, A));
            rtst = rtst.add(1.0, submatrix(P, A, B).mult(X));
        }

        Matrix Vall = Matrix.cellsum(sn.visits);
        Matrix Vst = new Matrix(M, K);
        for (int ist = 0; ist < M; ist++) {
            int isf = (int) sn.stationToStateful.get(ist);
            for (int r = 0; r < K; r++) {
                Vst.set(ist, r, Vall.get(isf, r));
            }
        }
        return new Pair<Matrix, Matrix>(rtst, Vst);
    }

    private static Matrix submatrix(Matrix P, int[] rows, int[] cols) {
        Matrix out = new Matrix(rows.length, cols.length);
        for (int i = 0; i < rows.length; i++) {
            for (int j = 0; j < cols.length; j++) {
                out.set(i, j, P.get(rows[i], cols[j]));
            }
        }
        return out;
    }
}
