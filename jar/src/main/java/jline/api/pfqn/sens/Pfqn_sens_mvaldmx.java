/**
 * Exact queue-length variances and covariances for mixed open/closed
 * product-form queueing networks with limited load dependence.
 *
 * <p>Mirrors the MATLAB reference {@code pfqn_sens_mvaldmx.m}.</p>
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.sens;

import java.util.ArrayList;
import java.util.List;

import jline.io.Ret;
import jline.util.Maths;
import jline.util.PopulationLattice;
import jline.util.matrix.Matrix;

public final class Pfqn_sens_mvaldmx {
    private Pfqn_sens_mvaldmx() {}

    /**
     * Exact second moments (variances and covariances) of the queue lengths of a
     * mixed open/closed product-form queueing network with limited load-dependent
     * service rates. This is the load-dependent and mixed counterpart of
     * {@link Pfqn_sens_mva}, which is restricted to closed load-independent
     * models.
     *
     * <p>The method is the moment analysis of Akyildiz and Strelen. Their Theorem
     * 1, equation (11), states that multiplying a queue-length moment by one
     * further factor Q_jT costs one derivative with respect to a parameter y_j
     * that scales the service demands s_ir of the classes r in T at station j:</p>
     *
     * <pre>
     *   E[Q_jT^k ...] = d/dy_j E[Q_jT^(k-1) ...]|_{y_j=1}
     *                   + nbar_jT E[Q_jT^(k-1) ...]
     * </pre>
     *
     * <p>Taking k=2 and T={s} gives the second moment, hence</p>
     *
     * <pre>
     *   Cov[n(i,r),n(j,s)] = d nbar(i,r) / dy_(j,s) |_{y=1}
     * </pre>
     *
     * <p>which is evaluated here by forward-mode differentiation of the mixed
     * load-dependent MVA of Bruell-Balbo-Afshari, i.e. of exactly the recursion
     * implemented by {@link jline.api.pfqn.ld.Pfqn_mvaldmx}. The differentiated
     * equations are (13) for the residence times, (15)-(17) for the conditional
     * marginal probabilities, (18) for the throughputs, (19) and (24)-(31) for
     * the effective capacities (delegated to {@link Pfqn_sens_ldmx_ec}), (32) for
     * the closed-class queue lengths and (33) for the open-class ones.</p>
     *
     * <p>Because a demand-scaling parameter perturbs the whole network through
     * the closed-class throughputs, the derivatives must be propagated for every
     * parameter, so the cross-station covariances come out at no extra cost and
     * are returned in {@code QCovFull}. This is unlike {@link Pfqn_sens_mva},
     * whose cheaper same-station recursion cannot reach them.</p>
     *
     * <p>Reference: I. F. Akyildiz and J. C. Strelen, "Moment Analysis for
     * Load-Dependent Mixed Product Form Queueing Networks", IEEE Trans.
     * Communications 39(6):828-832, 1991. The closed load-independent case
     * reduces to E. de Souza e Silva and R. R. Muntz, IEEE Trans. Computers
     * 37(9):1125-1129, 1988, which {@link Pfqn_sens_mva} implements directly.</p>
     *
     * <p>Open classes are supported: an open class contributes to the load Lo(i)
     * that drives the effective capacities, and equation (21) supplies the
     * corresponding dLo/dy. The moments of an open class are those of its queue
     * length at a station, which is finite even though its population is
     * infinite.</p>
     *
     * @param lambda arrival rate vector (1 x R); must be zero on closed classes
     * @param D      service demand matrix (M x R)
     * @param N      population vector (1 x R); Inf entries denote open classes
     * @param Z      think time vector (1 x R)
     * @param mu     load-dependent rate matrix (M x sum(N)), limited load dependence
     * @param S      number of servers per station (M x 1); accepted for signature
     *               compatibility with {@code pfqn_mvaldmx}, which likewise does not
     *               read it: the multiserver behaviour is carried entirely by mu
     * @return the base measures and their exact second moments
     */
    public static Ret.pfqnSensMvaldmx pfqn_sens_mvaldmx(Matrix lambda, Matrix D, Matrix N,
                                                      Matrix Z, Matrix mu, Matrix S) {
        int M = D.getNumRows();
        int R = D.getNumCols();

        double NfiniteSum = 0.0;
        for (int r = 0; r < N.length(); r++) {
            if (Double.isFinite(N.get(r))) {
                NfiniteSum += N.get(r);
            }
        }
        if (mu == null) {
            mu = Matrix.ones(M, (int) NfiniteSum);
            S = Matrix.ones(M, 1);
        }
        if (S == null) {
            S = Matrix.ones(M, 1);
        }
        if (mu.getNumCols() < NfiniteSum) {
            throw new RuntimeException("pfqn_sens_mvaldmx requires to specify the "
                    + "load-dependent rates with one job more than the maximum closed population.");
        }
        for (int r = 0; r < R; r++) {
            if (lambda.get(r) != 0.0 && N.get(r) > 0 && Double.isFinite(N.get(r))) {
                throw new RuntimeException("pfqn_sens_mvaldmx: Arrival rate cannot be specified "
                        + "on closed classes.");
            }
        }

        List<Integer> openClasses = new ArrayList<Integer>();
        List<Integer> closedClasses = new ArrayList<Integer>();
        for (int r = 0; r < R; r++) {
            if (Double.isInfinite(N.get(r))) {
                openClasses.add(Integer.valueOf(r));
            } else {
                closedClasses.add(Integer.valueOf(r));
            }
        }
        int C = closedClasses.size();
        if (C == 0) {
            throw new RuntimeException("pfqn_sens_mvaldmx: requires at least one closed class; "
                    + "use the open-class formulas directly otherwise.");
        }

        // up to sum(N)+1, limited load dependence
        Matrix muExt = new Matrix(mu.getNumRows(), mu.getNumCols() + 1);
        for (int i = 0; i < mu.getNumRows(); i++) {
            for (int j = 0; j < mu.getNumCols(); j++) {
                muExt.set(i, j, mu.get(i, j));
            }
            muExt.set(i, mu.getNumCols(), mu.get(i, mu.getNumCols() - 1));
        }
        Ret.pfqnSensLdmxEc ec = Pfqn_sens_ldmx_ec.pfqn_sens_ldmx_ec(lambda, D, muExt);
        Matrix EC = ec.EC;
        Matrix E = ec.E;
        Matrix Eprime = ec.Eprime;
        Matrix dECdLo = ec.dEC;

        Matrix Dc = new Matrix(M, C);
        double[] Nc = new double[C];
        double[] Zc = new double[C];
        for (int c = 0; c < C; c++) {
            int cls = closedClasses.get(c).intValue();
            for (int ist = 0; ist < M; ist++) {
                Dc.set(ist, c, D.get(ist, cls));
            }
            Nc[c] = N.get(cls);
            Zc[c] = Z.get(cls);
        }
        int NCtot = 0;
        for (int c = 0; c < C; c++) {
            NCtot += (int) Nc[c];
        }

        // ---- parameter list: y(j,r) multiplies the demand D(j,r) ------------
        int P = M * R;
        int[][] pidx = new int[M][R];
        int[] pj = new int[P];
        int[] pr = new int[P];
        int p = 0;
        for (int j = 0; j < M; j++) {
            for (int r = 0; r < R; r++) {
                pidx[j][r] = p;
                pj[p] = j;
                pr[p] = r;
                p++;
            }
        }
        // see _kb/03-api-layer.md for rationale
        double[][] dLo = new double[M][P];
        for (int q = 0; q < P; q++) {
            if (Double.isInfinite(N.get(pr[q]))) {
                dLo[pj[q]][q] = lambda.get(pr[q]) * D.get(pj[q], pr[q]);
            }
        }

        // ---- population recursion -------------------------------------------
        Matrix NcM = new Matrix(1, C);
        for (int c = 0; c < C; c++) {
            NcM.set(0, c, Nc[c]);
        }
        Matrix prods = new Matrix(1, C);
        for (int r = 0; r < C; r++) {
            double pd = 1.0;
            for (int i = 0; i < r; i++) {
                pd *= (Nc[i] + 1);
            }
            prods.set(0, r, pd);
        }
        int NT = 1;
        for (int c = 0; c < C; c++) {
            NT *= (int) (1 + Nc[c]);
        }
        double[][][] Pc = new double[M][1 + NCtot][NT];
        double[][][][] dPc = new double[M][1 + NCtot][NT][P];
        double[][] x = new double[C][NT];
        double[][][] dx = new double[C][NT][P];
        double[][][] w = new double[M][C][NT];
        double[][][][] dw = new double[M][C][NT][P];

        Matrix nvec = PopulationLattice.pprod(NcM);
        for (int ist = 0; ist < M; ist++) {
            Pc[ist][0][PopulationLattice.hashpop(nvec, NcM, C, prods)] = 1.0;   // eq. (16)
        }

        while (!(nvec.getNumRows() == 1 && nvec.getNumCols() == 1 && nvec.value() == -1.0)) {
            int hnvec = PopulationLattice.hashpop(nvec, NcM, C, prods);
            int nc = (int) nvec.elementSum();

            // ---- residence times, eq. (12) and its derivative eq. (13) ------
            for (int ist = 0; ist < M; ist++) {
                for (int c = 0; c < C; c++) {
                    if (nvec.get(c) > 0) {
                        int hnvec_c = PopulationLattice.hashpop(
                                Matrix.oner(nvec, Integer.valueOf(c)), NcM, C, prods);
                        int cls = closedClasses.get(c).intValue();
                        double acc = 0.0;
                        double[] dacc = new double[P];
                        for (int n = 1; n <= nc; n++) {
                            double Pprev = Pc[ist][n - 1][hnvec_c];
                            acc += n * EC.get(ist, n - 1) * Pprev;
                            for (int q = 0; q < P; q++) {
                                dacc[q] += n * (dECdLo.get(ist, n - 1) * dLo[ist][q] * Pprev
                                        + EC.get(ist, n - 1) * dPc[ist][n - 1][hnvec_c][q]);
                            }
                        }
                        w[ist][c][hnvec] = Dc.get(ist, c) * acc;
                        for (int q = 0; q < P; q++) {
                            double dwq = Dc.get(ist, c) * dacc[q];
                            if (pj[q] == ist && pr[q] == cls) {
                                dwq += Dc.get(ist, c) * acc;   // d(D*y)/dy = D
                            }
                            dw[ist][c][hnvec][q] = dwq;
                        }
                    }
                }
            }

            // ---- throughputs, eq. (18) --------------------------------------
            for (int c = 0; c < C; c++) {
                double sumw = 0.0;
                for (int ist = 0; ist < M; ist++) {
                    sumw += w[ist][c][hnvec];
                }
                double den = Zc[c] + sumw;
                x[c][hnvec] = nvec.get(c) / den;
                if (nvec.get(c) > 0) {
                    for (int q = 0; q < P; q++) {
                        double sdw = 0.0;
                        for (int ist = 0; ist < M; ist++) {
                            sdw += dw[ist][c][hnvec][q];
                        }
                        dx[c][hnvec][q] = -nvec.get(c) / (den * den) * sdw;
                    }
                }
            }

            // ---- conditional marginal probabilities, eq. (14)-(15) ----------
            for (int ist = 0; ist < M; ist++) {
                for (int n = 1; n <= nc; n++) {
                    for (int c = 0; c < C; c++) {
                        if (nvec.get(c) > 0) {
                            int hnvec_c = PopulationLattice.hashpop(
                                    Matrix.oner(nvec, Integer.valueOf(c)), NcM, C, prods);
                            int cls = closedClasses.get(c).intValue();
                            double Pprev = Pc[ist][n - 1][hnvec_c];
                            Pc[ist][n][hnvec] += Dc.get(ist, c) * EC.get(ist, n - 1)
                                    * x[c][hnvec] * Pprev;
                            for (int q = 0; q < P; q++) {
                                double dt = Dc.get(ist, c)
                                        * (dECdLo.get(ist, n - 1) * dLo[ist][q] * x[c][hnvec] * Pprev
                                        + EC.get(ist, n - 1) * dx[c][hnvec][q] * Pprev
                                        + EC.get(ist, n - 1) * x[c][hnvec] * dPc[ist][n - 1][hnvec_c][q]);
                                if (pj[q] == ist && pr[q] == cls) {
                                    dt += Dc.get(ist, c) * EC.get(ist, n - 1) * x[c][hnvec] * Pprev;
                                }
                                dPc[ist][n][hnvec][q] += dt;
                            }
                        }
                    }
                }
                // see _kb/03-api-layer.md for rationale
                double sumpc = 0.0;
                for (int k = 1; k <= nc; k++) {
                    sumpc += Pc[ist][k][hnvec];
                }
                Pc[ist][0][hnvec] = Maths.max(Math.ulp(1.0), 1 - sumpc);
                for (int q = 0; q < P; q++) {
                    double sd = 0.0;
                    for (int k = 1; k <= nc; k++) {
                        sd += dPc[ist][k][hnvec][q];
                    }
                    dPc[ist][0][hnvec][q] = -sd;
                }
            }

            nvec = PopulationLattice.pprod(nvec, NcM);
        }

        // ---- measures and their derivatives at the full population ----------
        int hnvec = PopulationLattice.hashpop(NcM, NcM, C, prods);
        Matrix XN = new Matrix(1, R);
        Matrix QN = new Matrix(M, R);
        Matrix UN = new Matrix(M, R);
        Matrix CN = new Matrix(M, R);
        double[][][] dQN = new double[M][R][P];

        // closed classes, eq. (32)
        for (int c = 0; c < C; c++) {
            int cls = closedClasses.get(c).intValue();
            XN.set(0, cls, x[c][hnvec]);
            int hnvec_c = PopulationLattice.hashpop(
                    Matrix.oner(NcM, Integer.valueOf(c)), NcM, C, prods);
            for (int ist = 0; ist < M; ist++) {
                CN.set(ist, cls, w[ist][c][hnvec]);
                QN.set(ist, cls, XN.get(0, cls) * CN.get(ist, cls));
                for (int q = 0; q < P; q++) {
                    dQN[ist][cls][q] = dx[c][hnvec][q] * w[ist][c][hnvec]
                            + x[c][hnvec] * dw[ist][c][hnvec][q];
                }
                double uacc = 0.0;
                for (int n = 1; n <= NCtot; n++) {
                    uacc += Dc.get(ist, c) * x[c][hnvec] * Eprime.get(ist, n - 1)
                            / E.get(ist, n - 1) * Pc[ist][n - 1][hnvec_c];
                }
                UN.set(ist, cls, uacc);
            }
        }

        // open classes, eq. (33)
        for (int ridx = 0; ridx < openClasses.size(); ridx++) {
            int r = openClasses.get(ridx).intValue();
            XN.set(0, r, lambda.get(r));
            for (int ist = 0; ist < M; ist++) {
                double acc = 0.0;
                double[] dacc = new double[P];
                for (int n = 0; n <= NCtot; n++) {
                    double Pn = Pc[ist][n][hnvec];
                    acc += (n + 1) * EC.get(ist, n) * Pn;
                    for (int q = 0; q < P; q++) {
                        dacc[q] += (n + 1) * (dECdLo.get(ist, n) * dLo[ist][q] * Pn
                                + EC.get(ist, n) * dPc[ist][n][hnvec][q]);
                    }
                }
                QN.set(ist, r, lambda.get(r) * D.get(ist, r) * acc);
                CN.set(ist, r, QN.get(ist, r) / lambda.get(r));
                for (int q = 0; q < P; q++) {
                    double dq = lambda.get(r) * D.get(ist, r) * dacc[q];
                    if (pj[q] == ist && pr[q] == r) {
                        dq += lambda.get(r) * D.get(ist, r) * acc;
                    }
                    dQN[ist][r][q] = dq;
                }
                double uacc = 0.0;
                for (int n = 0; n <= NCtot; n++) {
                    uacc += lambda.get(r) * Eprime.get(ist, n + 1) / E.get(ist, n + 1)
                            * Pc[ist][n][hnvec];
                }
                UN.set(ist, r, uacc);
            }
        }

        // ---- moments ---------------------------------------------------------
        // Cov[n(i,r),n(j,s)] = d nbar(i,r) / dy_(j,s)
        double asym = 0.0;
        Matrix[][] QCovFull = new Matrix[M][R];
        for (int i = 0; i < M; i++) {
            for (int r = 0; r < R; r++) {
                Matrix blk = new Matrix(M, R);
                for (int j = 0; j < M; j++) {
                    for (int s = 0; s < R; s++) {
                        double a = dQN[i][r][pidx[j][s]];
                        double bsym = dQN[j][s][pidx[i][r]];
                        asym = Math.max(asym, Math.abs(a - bsym));
                        blk.set(j, s, (a + bsym) / 2.0);
                    }
                }
                QCovFull[i][r] = blk;
            }
        }

        Matrix[] QCov = new Matrix[M];
        Matrix QVar = new Matrix(M, R);
        Matrix QTotVar = new Matrix(M, 1);
        for (int i = 0; i < M; i++) {
            Matrix cov = new Matrix(R, R);
            double tot = 0.0;
            for (int r = 0; r < R; r++) {
                for (int s = 0; s < R; s++) {
                    double v = QCovFull[i][r].get(i, s);
                    cov.set(r, s, v);
                    tot += v;
                }
                QVar.set(i, r, cov.get(r, r));
            }
            QCov[i] = cov;
            QTotVar.set(i, 0, tot);
        }

        return new Ret.pfqnSensMvaldmx(XN, QN, UN, CN, QCov, QCovFull, QVar, QTotVar, asym);
    }

    public static Ret.pfqnSensMvaldmx pfqn_sens_mvaldmx(Matrix lambda, Matrix D, Matrix N, Matrix Z) {
        return pfqn_sens_mvaldmx(lambda, D, N, Z, null, null);
    }
}
