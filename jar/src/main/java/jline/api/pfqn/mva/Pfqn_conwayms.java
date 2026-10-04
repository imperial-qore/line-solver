package jline.api.pfqn.mva;

import java.util.Arrays;

import org.apache.commons.math3.util.FastMath;

import jline.GlobalConstants;
import jline.io.Ret;
import jline.lang.constant.SchedStrategy;
import jline.util.Maths;
import jline.util.PopulationLattice;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

/**
 * Conway-Maxwell approximate MVA for multi-server queueing networks.
 */
public final class Pfqn_conwayms {
    private Pfqn_conwayms() {}

    public static Ret.pfqnAMVAMS pfqn_conwayms(Matrix L, Matrix N, Matrix Z, int[] nservers) {
        int M = L.getNumRows();
        SchedStrategy[] sched = new SchedStrategy[M];
        Arrays.fill(sched, SchedStrategy.FCFS);
        return pfqn_conwayms(L, N, Z, nservers, sched, GlobalConstants.FineTol, 1000);
    }

    public static Ret.pfqnAMVAMS pfqn_conwayms(Matrix L, Matrix N, Matrix Z) {
        int M = L.getNumRows();
        int[] nservers = Matrix.ones(1, M).toIntArray1D();
        SchedStrategy[] sched = new SchedStrategy[M];
        Arrays.fill(sched, SchedStrategy.FCFS);
        return pfqn_conwayms(L, N, Z, nservers, sched, GlobalConstants.FineTol, 1000);
    }

    public static Ret.pfqnAMVAMS pfqn_conwayms(Matrix L, Matrix N, Matrix Z, int[] nservers,
                                                SchedStrategy[] type, double tol, int maxiter) {
        return pfqn_conwayms(L, N, Z, nservers, type, tol, maxiter, null);
    }

    public static Ret.pfqnAMVAMS pfqn_conwayms(Matrix L, Matrix N, Matrix Z, int[] nservers,
                                                SchedStrategy[] type, double tol, int maxiter, Matrix QN0) {
        int M = L.getNumRows();
        int R = L.getNumCols();
        Z = Z.sumCols();

        MatrixCell Q = new MatrixCell(1 + R);
        MatrixCell PB = new MatrixCell(1 + R);
        MatrixCell P = new MatrixCell(1 + R);
        MatrixCell Delta = new MatrixCell(1 + R);
        int maxNumServers = Arrays.stream(nservers).max().orElse(1);
        for (int s = 0; s < 1 + R; s++) {
            Q.set(s, new Matrix(M, R));
            P.set(s, new Matrix(M, maxNumServers));
            PB.set(s, new Matrix(M, 1 + R));
            Delta.set(s, new Matrix(M, R)); // Delta(i,r,s): indexed by STATION and class, not class twice
        }
        for (int i = 0; i < M; i++) {
            for (int r = 0; r < R; r++) {
                for (int s = 0; s <= R; s++) {
                    Matrix N_1 = (s == 0) ? N : Matrix.oner(N, s - 1);
                    if (QN0 != null && !QN0.isEmpty()) {
                        Q.get(s).set(i, r, QN0.get(i, r));   // warm start
                    } else {
                        Q.get(s).set(i, r, N_1.get(r) / M);
                    }
                }
            }
        }
        for (int i = 0; i < M; i++) {
            for (int r = 0; r < R; r++) {
                for (int s = 0; s <= R; s++) {
                    Matrix N_1 = (s == 0) ? N : Matrix.oner(N, s - 1);
                    double pop = N_1.elementSum();
                    if (nservers[i] > 1) {
                        if (pop == 0) {
                            // empty network: the station is idle with probability one
                            for (int j = 1; j < nservers[i]; j++) P.get(s).set(i, j, 0);
                            PB.get(s).set(i, 0);
                            P.get(s).set(i, 0, 1);
                            continue;
                        }
                        for (int j = 1; j < nservers[i]; j++) {
                            P.get(s).set(i, j, 2 * Q.get(s).getRow(i).elementSum() / (pop * (pop + 1)));
                        }
                        if (pop > nservers[i] - 1) {
                            PB.get(s).set(i, 2 * Q.get(s).getRow(i).elementSum() / (pop + 1 - nservers[i]) / (pop * (pop + 1)));
                        } else { // fewer jobs than servers: they cannot all be busy
                            PB.get(s).set(i, 0);
                        }
                        P.get(s).set(i, 0, 1 - PB.get(s).get(i) - P.get(s).sumSubMatrix(i, i + 1, 1, nservers[i]));
                    }
                }
            }
        }

        int totiter = 0;
        for (int I = 0; I <= 1; I++) {
            for (int s = 0; s <= R; s++) {
                Matrix N_1 = (s == 0) ? N : Matrix.oner(N, s - 1);
                Ret.LinearizerResult coreResult = pfqn_conwayms_core(L, M, R, N_1, Z, nservers,
                        Q.get(s), P.get(s), PB.get(s), Delta, type, tol, maxiter - totiter);
                Q.get(s).setTo(coreResult.Q);
                P.get(s).setTo(coreResult.P);
                PB.get(s).setTo(coreResult.PB);
                totiter += coreResult.iter;
            }
            for (int i = 0; i < M; i++) {
                for (int r = 0; r < R; r++) {
                    for (int s = 0; s < R; s++) {
                        Matrix Ns = Matrix.oner(N, s);
                        if (N.get(s) > 2 && N.get(r) > 0) { // an empty class has no F_ir to correct
                            Delta.get(s).set(i, r, Q.get(1 + s).get(i, r) / Ns.get(r) - Q.get(0).get(i, r) / N.get(r));
                        }
                    }
                }
            }
        }

        Ret.LinearizerResult finalCoreResult = pfqn_conwayms_core(L, M, R, N, Z, nservers, Q.get(0), P.get(0), PB.get(0), Delta, type, tol, maxiter);
        Matrix retQ = finalCoreResult.Q;
        Matrix retW = finalCoreResult.W;
        Matrix retX = finalCoreResult.X;
        totiter += finalCoreResult.iter;

        Matrix retU = new Matrix(M, R);
        for (int i = 0; i < M; i++) {
            for (int r = 0; r < R; r++) {
                if (nservers[i] == 1) retU.set(i, r, retX.get(r) * L.get(i, r));
                else retU.set(i, r, retX.get(r) * L.get(i, r) / nservers[i]);
            }
        }
        Matrix retC = N.copy().elementDiv(retX).sub(Z);
        return new Ret.pfqnAMVAMS(retQ, retU, retW, retC, retX, totiter);
    }

    public static Ret.LinearizerResult pfqn_conwayms_core(Matrix L, int M, int R, Matrix N_1, Matrix Z, int[] nservers,
                                                          Matrix Q, Matrix P, Matrix PB, MatrixCell Delta,
                                                          SchedStrategy[] type, double tol, int maxiter) {
        boolean hasConverged = false;
        Matrix W = L.copy();
        Matrix Wlast = null;
        Matrix T = Matrix.createLike(L);
        int iter = 1;
        while (!hasConverged) {
            Matrix Qlast = Q.copy();
            Ret.pfqnEstimate estimateResult = pfqn_conwayms_estimate(M, R, N_1, nservers, Q, P, PB, Delta, W);
            Ret.LinearizerResult forwardMVAResult = pfqn_conwayms_forwardmva(L, M, R, N_1, Z, nservers, type,
                    estimateResult.Q_1, estimateResult.P_1, estimateResult.PB_1, estimateResult.T_1);
            Q = forwardMVAResult.Q;
            W = forwardMVAResult.W;
            T = forwardMVAResult.T;
            P = forwardMVAResult.P;
            PB = forwardMVAResult.PB;
            // W must enter the test: Q alone is satisfied on the FIRST sweep whenever Q
            // cannot move (M=1 seeds Q at its own fixed point), and the residence times
            // returned then are still the seed W=L, so T_1=Q/W is unbounded and the
            // throughput exceeds the station's own service capacity.
            double moved = Double.POSITIVE_INFINITY;
            if (Wlast != null) {
                moved = Math.max(Q.copy().sub(Qlast).norm(), W.copy().sub(Wlast).norm());
            }
            Wlast = W.copy();
            if (moved < tol || iter > maxiter) hasConverged = true;
            iter++;
        }
        return Ret.LinearizerResult.withX(Q, W, T, P, PB, iter);
    }

    public static Ret.pfqnEstimate pfqn_conwayms_estimate(int M, int R, Matrix N_1, int[] nservers,
                                                          Matrix Q, Matrix P, Matrix PB, MatrixCell Delta, Matrix W) {
        int maxNumServers = Arrays.stream(nservers).max().orElse(1);
        MatrixCell P_1 = new MatrixCell(1 + R);
        MatrixCell Q_1 = new MatrixCell(1 + R);
        for (int r = 0; r < R + 1; r++) {
            Q_1.set(r, new Matrix(M, R));
            P_1.set(r, new Matrix(M, maxNumServers));
        }
        Matrix PB_1 = new Matrix(M, 1 + R);
        Matrix T_1 = new Matrix(R, 1 + R);
        for (int i = 0; i < M; i++) {
            if (nservers[i] > 1) {
                for (int j = 0; j < nservers[i]; j++) {
                    for (int s = 0; s <= R; s++) P_1.get(s).set(i, j, P.get(i, j));
                }
                for (int s = 0; s <= R; s++) PB_1.set(i, s, PB.get(i, 0));
            }
            for (int r = 0; r < R; r++) {
                for (int s = 0; s < R; s++) {
                    Matrix Ns = Matrix.oner(N_1, s);
                    // a class with no jobs left has an empty queue everywhere
                    Q_1.get(1 + s).set(i, r, N_1.get(r) > 0
                            ? Ns.get(r) * (Q.get(i, r) / N_1.get(r) + Delta.get(s).get(i, r)) : 0);
                }
            }
        }
        // T_1 is Little's law over the queueing part of the cycle, sum_i Q_1 / sum_i W,
        // and not the ratio at the FIRST station with a positive residence time: the
        // per-station estimates disagree, so picking one made the answer depend on the
        // station order. The demand matrix carries no order, so a model symmetric under
        // permuting classes and stations together must return equal class throughputs,
        // and with the single-station pick it did not.
        for (int r = 0; r < R; r++) {
            for (int s = 0; s < R; s++) {
                Matrix Nr = Matrix.oner(N_1, r);
                double num = 0;
                double den = 0;
                for (int i = 0; i < M; i++) {
                    if (W.get(i, s) > 0 && N_1.get(s) > 0) { // a class with no jobs left has no throughput
                        // Delta is indexed (station, queued class, removed class), as the
                        // Q_1 loop above uses it: here class s queues and class r is removed
                        num += Nr.get(s) * (Q.get(i, s) / N_1.get(s) + Delta.get(r).get(i, s));
                        den += W.get(i, s);
                    }
                }
                // a reduced-population throughput cannot be negative; a negative one
                // makes log(F) complex in the XR/XE sums of the forward step
                if (den > 0) T_1.set(s, 1 + r, FastMath.max(0.0, num / den));
            }
        }
        return new Ret.pfqnEstimate(Q_1, P_1, PB_1, T_1);
    }

    public static Ret.LinearizerResult pfqn_conwayms_forwardmva(Matrix L, int M, int R, Matrix N_1, Matrix Z, int[] nservers,
                                                                SchedStrategy[] type, MatrixCell Q_1, MatrixCell P_1, Matrix PB_1, Matrix T_1) {
        int maxNumServers = Arrays.stream(nservers).max().orElse(1);
        Matrix W = new Matrix(M, R);
        Matrix T = new Matrix(1, R);
        Matrix Q = new Matrix(M, R);
        Matrix P = new Matrix(M, maxNumServers);
        Matrix PB = new Matrix(M, 1);
        Matrix XR = new Matrix(M, R);
        MatrixCell XE = new MatrixCell(R);
        for (int r = 0; r < R; r++) XE.set(r, new Matrix(M, R));

        MatrixCell F = new MatrixCell(R);
        for (int r = 0; r < R; r++) {
            F.set(r, new Matrix(M, R));
            for (int i = 0; i < M; i++) {
                double den = L.getRow(i).mult(T_1.getColumn(1 + r)).value();
                for (int s = 0; s < R; s++) F.get(r).set(i, s, L.get(i, s) * T_1.get(s, 1 + r) / den);
            }
        }

        Matrix mu = L.reciprocal();
        Matrix C = new Matrix(M, R + 1);
        for (int i = 0; i < M; i++) {
            for (int r = 0; r < R; r++) {
                if (nservers[i] > 1) {
                    XR.set(i, r, 0);
                    C.set(i, 1 + r, 0);
                    PopulationLattice.sprodResult sprodRes = PopulationLattice.sprod(R, Matrix.singleton((double) nservers[i]));
                    while (sprodRes.s.value() >= 0) {
                        if (Matrix.compare(sprodRes.n.transpose(), Matrix.oner(N_1, r), "lte")) {
                            Matrix n = sprodRes.n.copy().transpose();
                            double Ai = FastMath.exp(Maths.multinomialln(n) + n.mult(F.get(r).getRow(i).log().transpose()).value());
                            C.set(i, 1 + r, C.get(i, 1 + r) + Ai);
                            XR.set(i, r, XR.get(i, r) + Ai / mu.getRow(i).mult(n.transpose()).value());
                        }
                        sprodRes = PopulationLattice.sprod(sprodRes.s, sprodRes.S, sprodRes.D);
                    }
                    // Br empty: fewer jobs than servers, so all of them can never be busy
                    XR.set(i, r, C.get(i, 1 + r) > 0 ? XR.get(i, r) / C.get(i, 1 + r) : 0);
                }
            }
        }

        Matrix Cx = new Matrix(M, 1 + R);
        for (int i = 0; i < M; i++) {
            for (int r = 0; r < R; r++) {
                if (nservers[i] > 1) {
                    for (int c = 0; c < R; c++) {
                        XE.get(c).set(i, r, 0);
                        Cx.set(i, 1 + r, 0);
                        PopulationLattice.sprodResult sprodResult = PopulationLattice.sprod(R, Matrix.singleton((double) nservers[i]));
                        while (sprodResult.s.value() >= 0) {
                            if (Matrix.compare(sprodResult.n.transpose(), Matrix.oner(N_1, r), "lte") && sprodResult.n.get(c) >= 1) {
                                Matrix n = sprodResult.n.copy().transpose();
                                double Aix = FastMath.exp(Maths.multinomialln(n) + n.mult(F.get(r).getRow(i).log().transpose()).value());
                                Cx.set(i, 1 + r, Cx.get(i, 1 + r) + Aix);
                                XE.get(c).set(i, r, XE.get(c).get(i, r) + Aix / mu.getRow(i).mult(n.transpose()).value());
                            }
                            sprodResult = PopulationLattice.sprod(sprodResult.s, sprodResult.S, sprodResult.D);
                        }
                        // Axr empty: class c has no job left in N_1-e_r, so its term is 0
                        XE.get(c).set(i, r, Cx.get(i, 1 + r) > 0 ? XE.get(c).get(i, r) / Cx.get(i, 1 + r) : 0);
                    }
                }
            }
        }

        for (int i = 0; i < M; i++) {
            for (int r = 0; r < R; r++) {
                if (nservers[i] == 1) {
                    if (type[i] == SchedStrategy.FCFS) {
                        W.set(i, r, L.get(i, r));
                        for (int c = 0; c < R; c++) W.set(i, r, W.get(i, r) + L.get(i, c) * Q_1.get(1 + r).get(i, c));
                    } else {
                        W.set(i, r, L.get(i, r));
                        for (int c = 0; c < R; c++) W.set(i, r, W.get(i, r) + L.get(i, r) * Q_1.get(1 + r).get(i, c));
                    }
                } else {
                    W.set(i, r, L.get(i, r) + PB_1.get(i, 1 + r) * XR.get(i, r));
                    for (int c = 0; c < R; c++) {
                        W.set(i, r, W.get(i, r) + XE.get(c).get(i, r) * (Q_1.get(1 + r).get(i, c) - L.get(i, c) * T_1.get(c, 1 + r)));
                    }
                }
            }
        }
        for (int r = 0; r < R; r++) {
            T.set(r, N_1.get(r) / (Z.get(r) + W.getColumn(r).elementSum()));
            for (int i = 0; i < M; i++) Q.set(i, r, T.get(r) * W.get(i, r));
        }
        // Queue-length marginals. The relations
        //   p_j = A*p_{j-1}/j,  pB = A*(pB + p_{ms-1})/ms,  p_0 = 1 - pB - sum_j p_j
        // with A = sum_s X_s*L_is the mean number of busy servers are solved in closed
        // form rather than iterated. As a Jacobi iteration they amplify by A per sweep,
        // and since the convergence test watches Q and W but not P the routine returned
        // marginals whose mass had run to 334 behind the p_0 = max(0,1-...) floor.
        // The estimate step hands the same marginals to every reduced population, so the
        // population corrections that pfqn_linearizerms carries here are all zero.
        for (int i = 0; i < M; i++) {
            int ms = nservers[i];
            if (ms > 1) {
                double A = 0;
                for (int s = 0; s < R; s++) A += L.get(i, s) * T.get(s);
                for (int j = 0; j < ms; j++) P.set(i, j, 0.0);
                if (A >= ms) {
                    // Saturated: the closed form is singular and its limit is the degenerate
                    // marginal, every server busy with probability one. N = m with Z = 0
                    // reaches it exactly, so this is a legal input.
                    PB.set(i, 1.0);
                } else {
                    double[] alpha = new double[ms];
                    alpha[0] = 1.0;
                    double sumAlpha = 0;
                    for (int j = 1; j < ms; j++) {
                        alpha[j] = A * alpha[j - 1] / j;
                        sumAlpha += alpha[j];
                    }
                    double alphaB = A * alpha[ms - 1] / (ms - A);
                    double p0 = 1.0 / (1.0 + sumAlpha + alphaB);
                    P.set(i, 0, p0);
                    for (int j = 1; j < ms; j++) P.set(i, j, alpha[j] * p0);
                    PB.set(i, alphaB * p0);
                }
            }
        }
        return new Ret.LinearizerResult(Q, W, T, P, PB);
    }
}
