package jline.api.pfqn.mva;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import jline.GlobalConstants;
import jline.api.pfqn.mva.Pfqn_bs;
import jline.io.Ret;
import jline.lang.constant.SchedStrategy;
import jline.util.matrix.Matrix;

/**
 * Multi-server Krzesinski linearizer approximate MVA.
 */
public final class Pfqn_linearizerms {
    private Pfqn_linearizerms() {}

    public static Ret.pfqnAMVAMS pfqn_linearizerms(Matrix L, Matrix N, Matrix nservers) {
        return pfqn_linearizerms(L, N, new Matrix(1, L.getNumCols()), nservers);
    }

    public static Ret.pfqnAMVAMS pfqn_linearizerms(Matrix L, Matrix N, Matrix Z, Matrix nservers) {
        // one entry PER STATION: the forward MVA step indexes type by station
        List<SchedStrategy> type = new ArrayList<SchedStrategy>();
        for (int i = 0; i < L.getNumRows(); i++) {
            type.add(SchedStrategy.PS);
        }
        return pfqn_linearizerms(L, N, Z, nservers, type, GlobalConstants.FineTol, 1000);
    }

    public static Ret.pfqnAMVAMS pfqn_linearizerms(Matrix L, Matrix N, Matrix Z, Matrix nservers, List<SchedStrategy> type) {
        return pfqn_linearizerms(L, N, Z, nservers, type, GlobalConstants.FineTol, 1000);
    }

    public static Ret.pfqnAMVAMS pfqn_linearizerms(Matrix L, Matrix N, Matrix Z, Matrix nservers,
                                                   List<SchedStrategy> type, double tol, int maxiter) {
        return pfqn_linearizerms(L, N, Z, nservers, type, tol, maxiter, null);
    }

    public static Ret.pfqnAMVAMS pfqn_linearizerms(Matrix L, Matrix N, Matrix Z, Matrix nservers,
                                                   List<SchedStrategy> type, double tol, int maxiter, Matrix QN0) {
        int M = L.getNumRows();
        int R = L.getNumCols();
        if (Z.isEmpty()) Z = new Matrix(1, R);

        Matrix[] Q = new Matrix[M];
        Matrix PB = new Matrix(M, 1 + R);
        Matrix[] P = new Matrix[M];
        Matrix[] Delta = new Matrix[M];
        // Linearizer corrections for the queue-length marginals. Without them the
        // marginals stay at population N while the queue lengths are reduced to
        // N-e_s, which breaks Q + sum_j (m-1-j) p_j >= m-1 and lets W fall below the
        // mean service time.
        Matrix[] DeltaP = new Matrix[M];
        Matrix DeltaPB = new Matrix(M, R);
        for (int i = 0; i < M; i++) {
            Q[i] = new Matrix(R, 1 + R);
            P[i] = new Matrix((int) nservers.elementMax(), 1 + R);
            Delta[i] = new Matrix(R, R);
            DeltaP[i] = new Matrix((int) nservers.elementMax(), R);
        }

        for (int s = -1; s < R; s++) {
            Matrix N_1 = Matrix.oner(N, new ArrayList<Integer>(Collections.singletonList(s)));
            Ret.pfqnAMVA init0 = (QN0 == null || QN0.isEmpty())
                    ? Pfqn_bs.pfqn_bs(L, N_1, Z)
                    : Pfqn_bs.pfqn_bs(L, N_1, Z, tol, maxiter, QN0);
            for (int i = 0; i < M; i++) {
                for (int r = 0; r < R; r++) Q[i].set(r, 1 + s, init0.Q.get(i, r));
            }
        }

        for (int i = 0; i < M; i++) {
            for (int r = 0; r < R; r++) {
                for (int s = -1; s < R; s++) {
                    Matrix N_1 = Matrix.oner(N, new ArrayList<Integer>(Collections.singletonList(s)));
                    double pop = N_1.elementSum();
                    if (nservers.get(i) > 1) {
                        if (pop == 0) {
                            // empty network: the station is idle with probability one
                            for (int j = 0; j < nservers.get(i) - 1; j++) P[i].set(1 + j, 1 + s, 0);
                            PB.set(i, 1 + s, 0);
                            P[i].set(0, 1 + s, 1);
                            continue;
                        }
                        double sumQ = 0.0;
                        for (int k = 0; k < R; k++) sumQ += Q[i].get(k, 1 + s);
                        for (int j = 0; j < nservers.get(i) - 1; j++) {
                            P[i].set(1 + j, 1 + s, 2 * sumQ / (pop * (pop + 1)));
                        }
                        if (pop > nservers.get(i) - 1) {
                            PB.set(i, 1 + s, 2 * sumQ / (pop + 1 - nservers.get(i)) / (pop * (pop + 1)));
                        } else { // fewer jobs than servers: they cannot all be busy
                            PB.set(i, 1 + s, 0);
                        }
                        double sumP = 0.0;
                        for (int k = 0; k < nservers.get(i) - 1; k++) sumP += P[i].get(k + 1, 1 + s);
                        P[i].set(0, 1 + s, 1 - PB.get(i, 1 + s) - sumP);
                    }
                }
            }
        }

        int totiter = 0;
        for (int I = 0; I <= 1; I++) {
            for (int s = -1; s < R; s++) {
                Matrix N_1 = Matrix.oner(N, new ArrayList<Integer>(Collections.singletonList(s)));
                Matrix Q1 = new Matrix(M, R);
                Matrix P1 = new Matrix(M, (int) nservers.elementMax());
                Matrix PB1 = new Matrix(M, 1);
                for (int i = 0; i < M; i++) {
                    for (int j = 0; j < R; j++) Q1.set(i, j, Q[i].get(j, 1 + s));
                    for (int j = 0; j < P1.getNumCols(); j++) P1.set(i, j, P[i].get(j, 1 + s));
                    PB1.set(i, 0, PB.get(i, 1 + s));
                }
                Ret.LinearizerResult ret1 = linearizerms_core(L, M, R, N_1, Z, nservers, Q1, P1, PB1, Delta, DeltaP, DeltaPB, type, tol, maxiter - totiter);
                for (int i = 0; i < M; i++) {
                    for (int j = 0; j < R; j++) Q[i].set(j, 1 + s, ret1.Q.get(i, j));
                    for (int j = 0; j < (int) nservers.elementMax(); j++) P[i].set(j, 1 + s, ret1.P.get(i, j));
                    PB.set(i, 1 + s, ret1.PB.get(i));
                }
                totiter += ret1.iter;
            }
            for (int i = 0; i < M; i++) {
                for (int r = 0; r < R; r++) {
                    for (int s = 0; s < R; s++) {
                        Matrix Ns = Matrix.oner(N, new ArrayList<Integer>(Collections.singletonList(s)));
                        if (Ns.get(r) > 0) {
                            Delta[i].set(r, s, Q[i].get(r, 1 + s) / Ns.get(r) - Q[i].get(r, 0) / N.get(r));
                        } else { // Chandy-Neuse 0/0 convention: F_ir(N-e_s) = 0
                            Delta[i].set(r, s, -Q[i].get(r, 0) / N.get(r));
                        }
                    }
                }
            }
            // Probabilities do not scale with the population, so the analogue of
            // Delta is a plain difference
            for (int i = 0; i < M; i++) {
                if (nservers.get(i) > 1) {
                    for (int s = 0; s < R; s++) {
                        for (int j = 0; j < nservers.get(i); j++) {
                            DeltaP[i].set(j, s, P[i].get(j, 1 + s) - P[i].get(j, 0));
                        }
                        DeltaPB.set(i, s, PB.get(i, 1 + s) - PB.get(i, 0));
                    }
                }
            }
        }

        Matrix Q1 = new Matrix(M, R);
        Matrix P1 = new Matrix(M, (int) nservers.elementMax());
        Matrix PB1 = new Matrix(M, 1);
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < R; j++) Q1.set(i, j, Q[i].get(j, 0));
            for (int j = 0; j < P1.getNumCols(); j++) P1.set(i, j, P[i].get(j, 0));
            PB1.set(i, 0, PB.get(i, 0));
        }
        Ret.LinearizerResult ret1 = linearizerms_core(L, M, R, N, Z, nservers, Q1, P1, PB1, Delta, DeltaP, DeltaPB, type, tol, maxiter - totiter);
        totiter += ret1.iter;
        Matrix newQ = ret1.Q;
        Matrix W = ret1.W;
        Matrix X = ret1.T;
        Matrix U = new Matrix(M, R);
        for (int i = 0; i < M; i++) {
            for (int r = 0; r < R; r++) {
                if (nservers.get(i) == 1.0) U.set(i, r, X.get(r) * L.get(i, r));
                else U.set(i, r, X.get(r) * L.get(i, r) / nservers.get(i));
            }
        }
        Matrix C = new Matrix(1, R);
        for (int i = 0; i < R; i++) C.set(0, i, N.get(i) / X.get(i) - Z.get(i));
        return new Ret.pfqnAMVAMS(newQ, U, W, C, X, totiter);
    }

    static Ret.LinearizerResult linearizerms_core(Matrix L, int M, int R, Matrix N_1, Matrix Z, Matrix nservers,
                                                  Matrix Q, Matrix P, Matrix PB, Matrix[] Delta,
                                                  Matrix[] DeltaP, Matrix DeltaPB,
                                                  List<SchedStrategy> type, double tol, int maxiter) {
        int iter = 0;
        boolean hasConverged = false;
        Matrix W = null;
        Matrix T = null;
        while (!hasConverged) {
            iter++;
            Matrix Qlast = new Matrix(Q);
            Ret.pfqnLinearizerMSEstimate ret1 = linearizerms_estimate(M, R, N_1, nservers, Q, P, PB, Delta, DeltaP, DeltaPB);
            Ret.LinearizerResult ret2 = linearizerms_forwardMVA(L, M, R, N_1, Z, nservers, type, ret1.Q_1, ret1.P_1, ret1.PB_1);
            Q = ret2.Q;
            W = ret2.W;
            T = ret2.T;
            P = ret2.P;
            PB = ret2.PB;
            if (Q.sub(Qlast).norm() < tol || iter > maxiter) hasConverged = true;
        }
        return new Ret.LinearizerResult(Q, W, T, P, PB, iter);
    }

    static Ret.pfqnLinearizerMSEstimate linearizerms_estimate(int M, int R, Matrix N_1, Matrix nservers,
                                                              Matrix Q, Matrix P, Matrix PB, Matrix[] Delta,
                                                              Matrix[] DeltaP, Matrix DeltaPB) {
        Matrix[] P_1 = new Matrix[M];
        Matrix[] Q_1 = new Matrix[M];
        for (int i = 0; i < M; i++) {
            P_1[i] = new Matrix((int) nservers.elementMax(), 1 + R);
            Q_1[i] = new Matrix(R, 1 + R);
        }
        Matrix PB_1 = new Matrix(M, 1 + R);
        for (int i = 0; i < M; i++) {
            if (nservers.get(i) > 1) {
                for (int j = -1; j < nservers.get(i) - 1; j++) {
                    P_1[i].set(1 + j, 0, P.get(i, 1 + j));
                    for (int s = 0; s < R; s++) P_1[i].set(1 + j, 1 + s, P.get(i, 1 + j) + DeltaP[i].get(1 + j, s));
                }
                PB_1.set(i, 0, PB.get(i, 0));
                for (int s = 0; s < R; s++) PB_1.set(i, 1 + s, PB.get(i, 0) + DeltaPB.get(i, s));
            }
            for (int r = 0; r < R; r++) {
                for (int s = 0; s < R; s++) {
                    Matrix Ns = Matrix.oner(N_1, new ArrayList<Integer>(Collections.singletonList(s)));
                    if (N_1.get(r) > 0) {
                        Q_1[i].set(r, 1 + s, Ns.get(r) * (Q.get(i, r) / N_1.get(r) + Delta[i].get(r, s)));
                    } else { // a class with no jobs left has an empty queue everywhere
                        Q_1[i].set(r, 1 + s, 0);
                    }
                }
            }
        }
        return new Ret.pfqnLinearizerMSEstimate(Q_1, P_1, PB_1);
    }

    static Ret.LinearizerResult linearizerms_forwardMVA(Matrix L, int M, int R, Matrix N_1, Matrix Z, Matrix nservers,
                                                        List<SchedStrategy> type, Matrix[] Q_1, Matrix[] P_1, Matrix PB_1) {
        Matrix W = new Matrix(M, R);
        Matrix T = new Matrix(1, R);
        Matrix Q = new Matrix(M, R);
        Matrix P = new Matrix(M, (int) nservers.elementMax());
        Matrix PB = new Matrix(M, 1);
        for (int i = 0; i < M; i++) {
            for (int r = 0; r < R; r++) {
                W.set(i, r, L.get(i, r) / nservers.get(i));
                if (L.get(i, r) == 0.0) continue;
                // The FCFS arm (the other classes' demands, L(i,s)) is taken when EVERY
                // station is FCFS, which is what MATLAB's `if type == SchedStrategy.FCFS`
                // evaluates to on a vector. This test used to be open-coded as "true
                // unless ANY station is FCFS", i.e. exactly inverted: an all-PS model took
                // the FCFS formula and an FCFS model took the PS one. Single-class models
                // cannot see it -- L(i,s) IS L(i,r) there -- which is why it survived.
                boolean allFCFS = M > 0;
                for (int k = 0; k < M; k++) {
                    if (type.get(k) != SchedStrategy.FCFS) { allFCFS = false; break; }
                }
                if (allFCFS) {
                    for (int s = 0; s < R; s++) W.set(i, r, W.get(i, r) + (L.get(i, s) / nservers.get(i)) * Q_1[i].get(s, 1 + r));
                } else {
                    for (int s = 0; s < R; s++) W.set(i, r, W.get(i, r) + (L.get(i, r) / nservers.get(i)) * Q_1[i].get(s, 1 + r));
                }
                // Partially-idle-server correction. It compensates the 1/m scaling of
                // the arriving job's OWN service, so it carries L(i,r)/m and no sum over
                // the other classes: at N=e_r the terms must collapse to W = L(i,r).
                if (nservers.get(i) > 1) {
                    for (int j = 0; j <= nservers.get(i) - 2; j++) {
                        W.set(i, r, W.get(i, r) + (L.get(i, r) / nservers.get(i)) * (nservers.get(i) - 1 - j) * P_1[i].get(j, 1 + r));
                    }
                }
            }
        }
        for (int r = 0; r < R; r++) {
            T.set(r, N_1.get(r) / (Z.get(r) + Matrix.extractColumn(W, r, null).elementSum()));
            for (int i = 0; i < M; i++) Q.set(i, r, T.get(r) * W.get(i, r));
        }
        // Queue-length marginals. The relations
        //   p_j = (A p_{j-1} + d_{j-1}) / j,  pB = (A (pB + p_{ms-1}) + dB) / ms,
        //   p_0 = 1 - pB - sum_j p_j
        // with A = sum_s X_s L_is the mean number of busy servers and d the population
        // corrections, are solved in closed form rather than iterated: as a Jacobi
        // iteration they amplify by A per sweep and diverge once A approaches ms.
        for (int i = 0; i < M; i++) {
            int ms = (int) nservers.get(i);
            if (ms > 1 && !Double.isInfinite(nservers.get(i))) {
                double A = 0.0, dB = 0.0;
                double[] d = new double[ms];
                for (int s = 0; s < R; s++) {
                    double a_s = L.get(i, s) * T.get(s);
                    A += a_s;
                    for (int j = 0; j < ms; j++) d[j] += a_s * (P_1[i].get(j, 1 + s) - P_1[i].get(j, 0));
                    dB += a_s * (PB_1.get(i, 1 + s) - PB_1.get(i, 0));
                }
                if (A >= ms) {
                    throw new RuntimeException(String.format(
                            "pfqn_linearizerms: station %d offers %g busy servers out of %d; the model is "
                            + "saturated and its queue-length marginals do not exist.", i, A, ms));
                }
                double[] alpha = new double[ms];
                double[] beta = new double[ms];
                alpha[0] = 1.0;
                for (int j = 1; j < ms; j++) {
                    alpha[j] = A * alpha[j - 1] / j;
                    beta[j] = (A * beta[j - 1] + d[j - 1]) / j;
                }
                double alphaB = A * alpha[ms - 1] / (ms - A);
                double betaB = (A * beta[ms - 1] + dB + d[ms - 1]) / (ms - A);
                double sumAlpha = 0.0, sumBeta = 0.0;
                for (int j = 1; j < ms; j++) { sumAlpha += alpha[j]; sumBeta += beta[j]; }
                double p0 = (1 - sumBeta - betaB) / (1 + sumAlpha + alphaB);
                P.set(i, 0, p0);
                for (int j = 1; j < ms; j++) P.set(i, j, alpha[j] * p0 + beta[j]);
                PB.set(i, 0, alphaB * p0 + betaB);
            }
        }
        return new Ret.LinearizerResult(Q, W, T, P, PB);
    }
}
