/**
 * Multichain convolution algorithm for closed queueing networks with
 * class-dependent service rates.
 */
package jline.api.pfqn.ld;

import java.util.List;

import org.apache.commons.math3.special.Gamma;

import jline.util.SerializableFunction;
import jline.util.matrix.Matrix;

public final class Pfqn_conv {
    private Pfqn_conv() {}

    /**
     * Multichain convolution for networks with class-dependent service rates.
     *
     * <p>Implements the convolution algorithm of Sauer (1983), Section 5.2,
     * "Computational Algorithms for State-Dependent Queueing Networks",
     * ACM TOCS, Vol. 1, No. 1, pp. 67-92.</p>
     *
     * <p>For class-dependent stations, X_m(n) is computed recursively via Sauer
     * eq. (40): X_m(n) = (u_km / mu_km(n)) * X_m(n - e_k), where the
     * class-dependence function supplies mu_km(n) = (n_k/|n|) * beta_{m,k}(n),
     * beta being the DIMENSIONLESS scaling of the service demand (L/beta).
     * cdscaling.get(m) is a function of the per-class population vector n at
     * station m, returning either a 1x1 matrix (chain-independent) or a length-R
     * row vector of per-class rates. Any saturation/cutoff is applied inside the
     * function.</p>
     *
     * @param L         Service demand matrix (M x R)
     * @param N         Population vector, must be finite (closed network)
     * @param Z         Think time vector
     * @param cdscaling class-dependence functions beta_m(n) indexed by station;
     *                  null entries denote load-independent stations
     * @return {G, log(G)}
     */
    public static double[] pfqn_conv(Matrix L, int[] N, double[] Z,
                                     List<SerializableFunction<Matrix, Matrix>> cdscaling) {
        int M = L.getNumRows();
        int R = L.getNumCols();

        int stateSpaceSize = 1;
        for (int n : N) stateSpaceSize *= (n + 1);

        boolean[] isCd = new boolean[M];
        for (int ist = 0; ist < M; ist++) {
            isCd[ist] = cdscaling != null && ist < cdscaling.size() && cdscaling.get(ist) != null;
        }

        double[][] Xm = new double[M][];
        for (int ist = 0; ist < M; ist++) {
            if (isCd[ist]) {
                Xm[ist] = new double[stateSpaceSize];
                Xm[ist][0] = 1.0;

                int[] n = new int[R];
                while (true) {
                    int idx = hashpop(n, N);
                    int sumN = 0;
                    for (int v : n) sumN += v;
                    if (sumN == 0) {
                        Xm[ist][idx] = 1.0;
                    } else {
                        for (int r = 0; r < R; r++) {
                            if (n[r] > 0) {
                                // see _kb/03-api-layer.md for rationale
                                Matrix nvec = new Matrix(1, R);
                                int tot = 0;
                                for (int k = 0; k < R; k++) {
                                    nvec.set(0, k, (double) n[k]);
                                    tot += n[k];
                                }
                                Matrix bval = cdscaling.get(ist).apply(nvec);
                                double beta = (bval.length() > 1) ? bval.get(r) : bval.get(0);

                                // see _kb/03-api-layer.md for rationale
                                double nr = (double) n[r];
                                n[r]--;
                                int idxPrev = hashpop(n, N);
                                n[r]++;

                                if (beta > 0) {
                                    Xm[ist][idx] = ((double) tot / nr) * (L.get(ist, r) / beta) * Xm[ist][idxPrev];
                                }
                                break;
                            }
                        }
                    }
                    if (!pprodNext(n, N)) break;
                }
            }
        }

        double[] gCurr = new double[stateSpaceSize];

        int[] n = new int[R];
        while (true) {
            int idx = hashpop(n, N);
            gCurr[idx] = fz(Z, n);
            if (!pprodNext(n, N)) break;
        }

        for (int ist = 0; ist < M; ist++) {
            if (isCd[ist]) {
                double[] gOld = gCurr.clone();
                gCurr = new double[stateSpaceSize];

                int[] nOuter = new int[R];
                while (true) {
                    int idxN = hashpop(nOuter, N);
                    double convSum = 0.0;

                    int[] i = new int[R];
                    while (true) {
                        int idxI = hashpop(i, N);
                        int[] nmi = new int[R];
                        for (int k = 0; k < R; k++) nmi[k] = nOuter[k] - i[k];
                        int idxNmi = hashpop(nmi, N);
                        convSum += Xm[ist][idxI] * gOld[idxNmi];
                        if (!pprodNextBounded(i, nOuter)) break;
                    }

                    gCurr[idxN] = convSum;
                    if (!pprodNext(nOuter, N)) break;
                }
            } else {
                int[] nInner = new int[R];
                while (true) {
                    int idxN = hashpop(nInner, N);
                    for (int r = 0; r < R; r++) {
                        if (nInner[r] >= 1) {
                            nInner[r]--;
                            int idxN1r = hashpop(nInner, N);
                            nInner[r]++;
                            gCurr[idxN] += L.get(ist, r) * gCurr[idxN1r];
                        }
                    }
                    if (!pprodNext(nInner, N)) break;
                }
            }
        }

        double G = gCurr[stateSpaceSize - 1];
        double lG = (G > 0) ? Math.log(G) : Double.NEGATIVE_INFINITY;
        return new double[] { G, lG };
    }

    private static int hashpop(int[] n, int[] N) {
        int idx = 0;
        int stride = 1;
        for (int r = 0; r < n.length; r++) {
            idx += stride * n[r];
            stride *= (N[r] + 1);
        }
        return idx;
    }

    private static boolean pprodNext(int[] n, int[] N) {
        int R = n.length;
        boolean atMax = true;
        for (int i = 0; i < R; i++) {
            if (n[i] != N[i]) { atMax = false; break; }
        }
        if (atMax) return false;

        int s = R - 1;
        while (s >= 0 && n[s] == N[s]) {
            n[s] = 0;
            s--;
        }
        if (s >= 0) n[s]++;
        return true;
    }

    private static boolean pprodNextBounded(int[] i, int[] upper) {
        int R = i.length;
        boolean atMax = true;
        for (int k = 0; k < R; k++) {
            if (i[k] != upper[k]) { atMax = false; break; }
        }
        if (atMax) return false;

        int s = R - 1;
        while (s >= 0 && i[s] == upper[s]) {
            i[s] = 0;
            s--;
        }
        if (s >= 0) i[s]++;
        return true;
    }

    private static double fz(double[] Z, int[] n) {
        int sumN = 0;
        for (int v : n) sumN += v;
        if (sumN == 0) return 1.0;
        double logF = 0.0;
        for (int r = 0; r < n.length; r++) {
            if (Z[r] > 0) {
                logF += Math.log(Z[r]) * n[r];
                logF -= Gamma.logGamma((double) (1 + n[r]));
            } else if (n[r] > 0) {
                return 0.0;
            }
        }
        return Math.exp(logF);
    }
}
