/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.junit.jupiter.api.Test;

import jline.api.pfqn.ld.Pfqn_ldmx_ec;
import jline.api.pfqn.ld.Pfqn_mvaldmx;
import jline.api.pfqn.sens.Pfqn_sens_ldmx_ec;
import jline.api.pfqn.sens.Pfqn_sens_mva;
import jline.api.pfqn.sens.Pfqn_sens_mvaldmx;
import jline.io.Ret;
import jline.util.matrix.Matrix;

import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Validates {@link Pfqn_sens_mvaldmx#pfqn_sens_mvaldmx}, the mixed load-dependent
 * moment analysis of Akyildiz and Strelen (1991), mirroring the MATLAB harness
 * {@code pfqn_sens_mvaldmx_validate.m}. Five independent references, chosen so
 * that every channel of the derivation is exercised by something that does not
 * share its code:
 *
 * <ul>
 *   <li>A. {@link Pfqn_mvaldmx} for the base measures X, Q, U, R. The primal must
 *       be reproduced entry by entry, otherwise the derivative is of the wrong
 *       function.</li>
 *   <li>B. central finite differences of {@link Pfqn_mvaldmx} with respect to the
 *       demand-scaling parameter y(j,s). This checks the differentiated
 *       recursion itself, including the load-dependent channel dEC/dLo and the
 *       open-class channel dLo/dy of eq. (21), but does not check the identity
 *       Cov = d nbar / dy.</li>
 *   <li>C. {@link Pfqn_sens_mva} in the closed load-independent limit. This checks
 *       the identity against the independently validated de Souza e Silva and
 *       Muntz recursion.</li>
 *   <li>D. brute-force enumeration of the product-form equilibrium distribution.
 *       Closed load-dependent models are enumerated exactly; mixed models are
 *       enumerated with the open populations truncated, which converges
 *       geometrically and is therefore checked at a looser tolerance. This is the
 *       only check that closes the loop on the identity in the mixed
 *       load-dependent case.</li>
 *   <li>E. symmetry of QCovFull. Cov[n(i,r),n(j,s)] and Cov[n(j,s),n(i,r)] are
 *       computed by differentiating two different classes' equations, so their
 *       agreement is a nontrivial structural check.</li>
 * </ul>
 *
 * <p>{@link Pfqn_sens_ldmx_ec} is additionally pinned against
 * {@link Pfqn_ldmx_ec} on its primal terms and against finite differences on its
 * derivatives.</p>
 */
public class PfqnSensMvaldmxTest {

    private static final double TOL_MVA = 1e-12;
    private static final double TOL_FD = 1e-6;
    private static final double TOL_MOM = 1e-9;
    private static final double TOL_BRT = 1e-8;     // exact enumeration, closed load-dependent
    private static final double TOL_BRT_T = 5e-5;   // truncated enumeration, mixed
    private static final double TOL_SYM = 1e-8;

    private static Matrix row(double... v) {
        Matrix m = new Matrix(1, v.length);
        for (int i = 0; i < v.length; i++) {
            m.set(0, i, v[i]);
        }
        return m;
    }

    private static double[] flat(Matrix m) {
        double[] v = new double[m.getNumRows() * m.getNumCols()];
        int k = 0;
        for (int i = 0; i < m.getNumRows(); i++) {
            for (int j = 0; j < m.getNumCols(); j++) {
                v[k++] = m.get(i, j);
            }
        }
        return v;
    }

    private static double[] flat(Matrix[] arr) {
        double[] v = new double[arr.length * arr[0].getNumRows() * arr[0].getNumCols()];
        int k = 0;
        for (int i = 0; i < arr.length; i++) {
            for (int r = 0; r < arr[i].getNumRows(); r++) {
                for (int s = 0; s < arr[i].getNumCols(); s++) {
                    v[k++] = arr[i].get(r, s);
                }
            }
        }
        return v;
    }

    private static double relerr(double[] a, double[] b) {
        double e = 0.0;
        for (int k = 0; k < a.length; k++) {
            double scale = Math.max(1.0, Math.max(Math.abs(a[k]), Math.abs(b[k])));
            e = Math.max(e, Math.abs(a[k] - b[k]) / scale);
        }
        return e;
    }

    private static double relerr(Matrix a, Matrix b) {
        return relerr(flat(a), flat(b));
    }

    // =====================================================================
    // A/B/E on random mixed load-dependent models
    // =====================================================================
    @Test
    public void baseMeasuresFiniteDifferencesAndSymmetry() {
        Random rng = new Random(1);
        double errMva = 0.0;
        double errFd = 0.0;
        double errSym = 0.0;
        int nFd = 0;
        int nModels = 0;

        for (int trial = 1; trial <= 12; trial++) {
            int M = 1 + rng.nextInt(2);
            int Ropen = rng.nextInt(2);
            int R = 1 + Ropen;
            Matrix D = new Matrix(M, R);
            for (int i = 0; i < M; i++) {
                for (int r = 0; r < R; r++) {
                    D.set(i, r, 0.2 + 0.6 * rng.nextDouble());
                }
            }
            Matrix N = new Matrix(1, R);
            N.set(0, 0, 1 + rng.nextInt(3));              // closed class
            Matrix lambda = new Matrix(1, R);
            if (Ropen == 1) {
                N.set(0, 1, Double.POSITIVE_INFINITY);    // open class
                lambda.set(0, 1, 0.05 + 0.15 * rng.nextDouble());
            }
            Matrix Z = new Matrix(1, R);
            Z.set(0, 0, 0.5 * rng.nextDouble());
            int NCtot = (int) N.get(0, 0);
            // limited load dependence: rates grow up to level b then saturate
            int b = 1 + rng.nextInt(3);
            Matrix mu = new Matrix(M, Math.max(NCtot, 1));
            for (int i = 0; i < M; i++) {
                for (int n = 1; n <= mu.getNumCols(); n++) {
                    mu.set(i, n - 1, Math.min(n, b) * (0.8 + 0.4 * rng.nextDouble()));
                }
            }
            // keep the geometric tail of the limited load dependence stable
            boolean skip = false;
            for (int i = 0; i < M; i++) {
                double Lo = 0.0;
                for (int r = 0; r < R; r++) {
                    Lo += lambda.get(0, r) * D.get(i, r);
                }
                if (Lo / mu.get(i, mu.getNumCols() - 1) > 0.6) {
                    skip = true;
                }
            }
            if (skip) {
                continue;
            }
            nModels++;

            Matrix S = Matrix.ones(M, 1);
            Ret.pfqnSensMvaldmx mom = Pfqn_sens_mvaldmx.pfqn_sens_mvaldmx(lambda, D, N, Z, mu, S);

            // ---- A. base measures ----------------------------------------
            Ret.pfqnMVALDMX ref = Pfqn_mvaldmx.pfqn_mvaldmx(lambda, D, N, Z, mu, S);
            errMva = Math.max(errMva, relerr(mom.X, ref.X));
            errMva = Math.max(errMva, relerr(mom.Q, ref.Q));
            errMva = Math.max(errMva, relerr(mom.U, ref.U));
            errMva = Math.max(errMva, relerr(mom.R, ref.R));

            // ---- E. symmetry ---------------------------------------------
            errSym = Math.max(errSym, mom.QCovAsym);

            // ---- B. finite differences -----------------------------------
            double h = 1e-6;
            for (int j = 0; j < M; j++) {
                for (int s = 0; s < R; s++) {
                    if (D.get(j, s) <= 0) {
                        continue;
                    }
                    Matrix Dp = D.copy();
                    Matrix Dm = D.copy();
                    Dp.set(j, s, D.get(j, s) * (1 + h));
                    Dm.set(j, s, D.get(j, s) * (1 - h));
                    Ret.pfqnMVALDMX rp = Pfqn_mvaldmx.pfqn_mvaldmx(lambda, Dp, N, Z, mu, S);
                    Ret.pfqnMVALDMX rm = Pfqn_mvaldmx.pfqn_mvaldmx(lambda, Dm, N, Z, mu, S);
                    Matrix fd = new Matrix(M, R);
                    Matrix an = new Matrix(M, R);
                    for (int i = 0; i < M; i++) {
                        for (int r = 0; r < R; r++) {
                            fd.set(i, r, (rp.Q.get(i, r) - rm.Q.get(i, r)) / (2 * h));
                            an.set(i, r, mom.QCovFull[i][r].get(j, s));
                        }
                    }
                    errFd = Math.max(errFd, relerr(an, fd));
                    nFd++;
                }
            }
        }

        assertTrue(nModels > 0, "no mixed load-dependent model was generated");
        assertTrue(nFd > 0, "no finite-difference parameter was checked");
        assertTrue(errMva <= TOL_MVA, "pfqn_mvaldmx base measures: " + errMva);
        assertTrue(errFd <= TOL_FD, "finite differences (" + nFd + " params): " + errFd);
        assertTrue(errSym <= TOL_SYM, "QCovFull symmetry (raw): " + errSym);
    }

    // =====================================================================
    // C. closed load-independent limit against pfqn_sens_mva
    // =====================================================================
    @Test
    public void closedLoadIndependentLimitMatchesSensMva() {
        Random rng = new Random(2);
        double errMom = 0.0;
        int nMom = 0;

        for (int trial = 1; trial <= 12; trial++) {
            int M = 1 + rng.nextInt(3);
            int R = 1 + rng.nextInt(2);
            Matrix D = new Matrix(M, R);
            for (int i = 0; i < M; i++) {
                for (int r = 0; r < R; r++) {
                    D.set(i, r, 0.2 + rng.nextDouble());
                }
            }
            Matrix N = new Matrix(1, R);
            int Ntot = 0;
            for (int r = 0; r < R; r++) {
                int nr = 1 + rng.nextInt(3);
                N.set(0, r, nr);
                Ntot += nr;
            }
            Matrix Z = new Matrix(1, R);
            for (int r = 0; r < R; r++) {
                Z.set(0, r, 0.4 * rng.nextDouble());
            }
            Matrix lambda = new Matrix(1, R);
            Matrix mu = Matrix.ones(M, Ntot);

            Ret.pfqnSensMvaldmx mom =
                    Pfqn_sens_mvaldmx.pfqn_sens_mvaldmx(lambda, D, N, Z, mu, Matrix.ones(M, 1));
            Ret.pfqnSensMva ref = Pfqn_sens_mva.pfqn_sens_mva(D, N, Z);
            errMom = Math.max(errMom, relerr(mom.Q, ref.Q));
            errMom = Math.max(errMom, relerr(flat(mom.QCov), flat(ref.QCov)));
            errMom = Math.max(errMom, relerr(mom.QVar, ref.QVar));
            errMom = Math.max(errMom, relerr(mom.QTotVar, ref.QTotVar));
            nMom++;
        }

        assertTrue(nMom == 12, "unexpected number of closed load-independent models: " + nMom);
        assertTrue(errMom <= TOL_MOM, "pfqn_sens_mva closed LI (" + nMom + " models): " + errMom);
    }

    // =====================================================================
    // D1. brute force, closed load-dependent, exact enumeration
    // =====================================================================
    @Test
    public void closedLoadDependentMatchesBruteForce() {
        Random rng = new Random(3);
        double errBrt = 0.0;
        int nBrt = 0;

        for (int trial = 1; trial <= 10; trial++) {
            int M = 2;
            int R = 1;
            Matrix D = new Matrix(M, R);
            for (int i = 0; i < M; i++) {
                D.set(i, 0, 0.3 + 0.5 * rng.nextDouble());
            }
            Matrix N = row(2 + rng.nextInt(3));
            Matrix Z = row(0.3 * rng.nextDouble());
            Matrix lambda = row(0.0);
            int b = 2 + rng.nextInt(2);
            Matrix mu = new Matrix(M, (int) N.get(0, 0));
            for (int i = 0; i < M; i++) {
                for (int n = 1; n <= mu.getNumCols(); n++) {
                    mu.set(i, n - 1, Math.min(n, b) * (0.8 + 0.4 * rng.nextDouble()));
                }
            }
            Ret.pfqnSensMvaldmx mom =
                    Pfqn_sens_mvaldmx.pfqn_sens_mvaldmx(lambda, D, N, Z, mu, Matrix.ones(M, 1));
            BruteLdmx bm = bruteLdmx(lambda, D, N, Z, mu, 0);
            errBrt = Math.max(errBrt, relerr(mom.Q, bm.Q));
            errBrt = Math.max(errBrt, relerr(flat(mom.QCov), flat(bm.QCov)));
            nBrt++;
        }

        assertTrue(nBrt == 10, "unexpected number of closed load-dependent models: " + nBrt);
        assertTrue(errBrt <= TOL_BRT, "brute force closed LD (" + nBrt + " models): " + errBrt);
    }

    // =====================================================================
    // D2. brute force, mixed load-dependent, truncated enumeration
    // =====================================================================
    @Test
    public void mixedLoadDependentMatchesTruncatedBruteForce() {
        Random rng = new Random(4);
        double errBrtT = 0.0;
        int nBrtT = 0;

        for (int trial = 1; trial <= 6; trial++) {
            int M = 2;
            int R = 2;
            Matrix D = new Matrix(M, R);
            for (int i = 0; i < M; i++) {
                for (int r = 0; r < R; r++) {
                    D.set(i, r, 0.3 + 0.4 * rng.nextDouble());
                }
            }
            Matrix N = row(1 + rng.nextInt(2), Double.POSITIVE_INFINITY);
            Matrix Z = row(0.3 * rng.nextDouble(), 0.0);
            Matrix lambda = row(0.0, 0.05 + 0.1 * rng.nextDouble());
            int b = 1 + rng.nextInt(2);
            Matrix mu = new Matrix(M, (int) N.get(0, 0));
            for (int i = 0; i < M; i++) {
                for (int n = 1; n <= mu.getNumCols(); n++) {
                    mu.set(i, n - 1, Math.min(n, b) * (1.0 + 0.3 * rng.nextDouble()));
                }
            }
            boolean skip = false;
            for (int i = 0; i < M; i++) {
                double Lo = 0.0;
                for (int r = 0; r < R; r++) {
                    Lo += lambda.get(0, r) * D.get(i, r);
                }
                if (Lo / mu.get(i, mu.getNumCols() - 1) > 0.4) {
                    skip = true;
                }
            }
            if (skip) {
                continue;
            }
            Ret.pfqnSensMvaldmx mom =
                    Pfqn_sens_mvaldmx.pfqn_sens_mvaldmx(lambda, D, N, Z, mu, Matrix.ones(M, 1));
            BruteLdmx bm = bruteLdmx(lambda, D, N, Z, mu, 60);
            errBrtT = Math.max(errBrtT, relerr(mom.Q, bm.Q));
            errBrtT = Math.max(errBrtT, relerr(flat(mom.QCov), flat(bm.QCov)));
            nBrtT++;
        }

        assertTrue(nBrtT > 0, "no mixed load-dependent model survived the stability filter");
        assertTrue(errBrtT <= TOL_BRT_T,
                "brute force mixed LD (" + nBrtT + " models): " + errBrtT);
    }

    // =====================================================================
    // Pfqn_sens_ldmx_ec: primal against Pfqn_ldmx_ec, derivatives against
    // finite differences in the open-class load.
    // =====================================================================
    @Test
    public void ldmxEcPrimalMatchesEcAndDerivativesMatchFiniteDifferences() {
        Random rng = new Random(5);
        double errPrimal = 0.0;
        double errDeriv = 0.0;
        int nCases = 0;

        for (int trial = 1; trial <= 10; trial++) {
            int M = 1 + rng.nextInt(2);
            int R = 2;
            Matrix D = new Matrix(M, R);
            for (int i = 0; i < M; i++) {
                for (int r = 0; r < R; r++) {
                    D.set(i, r, 0.3 + 0.4 * rng.nextDouble());
                }
            }
            Matrix lambda = row(0.0, 0.05 + 0.1 * rng.nextDouble());
            int Nt = 2 + rng.nextInt(3);
            int b = 1 + rng.nextInt(3);
            Matrix mu = new Matrix(M, Nt);
            for (int i = 0; i < M; i++) {
                for (int n = 1; n <= Nt; n++) {
                    mu.set(i, n - 1, Math.min(n, b) * (1.0 + 0.3 * rng.nextDouble()));
                }
            }
            boolean skip = false;
            for (int i = 0; i < M; i++) {
                double Lo = 0.0;
                for (int r = 0; r < R; r++) {
                    Lo += lambda.get(0, r) * D.get(i, r);
                }
                if (Lo / mu.get(i, Nt - 1) > 0.4) {
                    skip = true;
                }
            }
            if (skip) {
                continue;
            }
            nCases++;

            Ret.pfqnSensLdmxEc s = Pfqn_sens_ldmx_ec.pfqn_sens_ldmx_ec(lambda, D, mu);
            // Pfqn_ldmx_ec mutates its mu argument in place, so hand it a copy.
            Ret.pfqnLDMXEC ref = Pfqn_ldmx_ec.pfqn_ldmx_ec(lambda, D, new Matrix(mu));
            errPrimal = Math.max(errPrimal, relerr(s.EC, ref.EC));
            errPrimal = Math.max(errPrimal, relerr(s.E, ref.E));
            errPrimal = Math.max(errPrimal, relerr(s.Eprime, ref.Eprime));
            errPrimal = Math.max(errPrimal, relerr(s.Lo, ref.Lo));

            // dEC/dLo etc. by central differences in Lo, perturbed through the
            // open-class arrival rate of the single station under test.
            for (int i = 0; i < M; i++) {
                double h = 1e-7;
                double dLambda = h / D.get(i, 1);
                Matrix lp = lambda.copy();
                Matrix lm = lambda.copy();
                // perturb Lo(i) only: use a demand matrix that isolates station i
                Matrix Di = D.copy();
                for (int j = 0; j < M; j++) {
                    if (j != i) {
                        Di.set(j, 1, 0.0);
                    }
                }
                lp.set(0, 1, lambda.get(0, 1) + dLambda);
                lm.set(0, 1, lambda.get(0, 1) - dLambda);
                Ret.pfqnSensLdmxEc base = Pfqn_sens_ldmx_ec.pfqn_sens_ldmx_ec(lambda, Di, mu);
                Ret.pfqnSensLdmxEc sp = Pfqn_sens_ldmx_ec.pfqn_sens_ldmx_ec(lp, Di, mu);
                Ret.pfqnSensLdmxEc sm = Pfqn_sens_ldmx_ec.pfqn_sens_ldmx_ec(lm, Di, mu);
                double step = sp.Lo.get(i, 0) - sm.Lo.get(i, 0);
                for (int n = 0; n < base.EC.getNumCols(); n++) {
                    double fd = (sp.EC.get(i, n) - sm.EC.get(i, n)) / step;
                    errDeriv = Math.max(errDeriv, relerr(new double[]{base.dEC.get(i, n)},
                            new double[]{fd}));
                }
                for (int n = 0; n < base.E.getNumCols(); n++) {
                    double fdE = (sp.E.get(i, n) - sm.E.get(i, n)) / step;
                    double fdEp = (sp.Eprime.get(i, n) - sm.Eprime.get(i, n)) / step;
                    errDeriv = Math.max(errDeriv, relerr(new double[]{base.dE.get(i, n)},
                            new double[]{fdE}));
                    errDeriv = Math.max(errDeriv, relerr(new double[]{base.dEprime.get(i, n)},
                            new double[]{fdEp}));
                }
            }
        }

        assertTrue(nCases > 0, "no EC case survived the stability filter");
        assertTrue(errPrimal <= 1e-12, "pfqn_sens_ldmx_ec primal vs pfqn_ldmx_ec: " + errPrimal);
        assertTrue(errDeriv <= TOL_FD, "pfqn_sens_ldmx_ec derivatives vs finite differences: "
                + errDeriv);
    }

    // =====================================================================
    private static final class BruteLdmx {
        Matrix Q;
        Matrix[] QCov;
    }

    /**
     * Exact moments by enumerating the mixed load-dependent product form
     * p(n) ~ prod_i [ n_i! prod_r a(i,r)^n(i,r)/n(i,r)! prod_{j=1}^{n_i} 1/mu(i,j) ]
     *        * prod_{closed c} Z(c)^n(0,c)/n(0,c)!
     * with a(i,r) = D(i,r) for a closed class and a(i,r) = lambda(r)*D(i,r) for
     * an open class, n_i the total population at station i, and the closed
     * classes constrained to sum to N. Open classes are truncated at Kopen jobs
     * per station.
     */
    private static BruteLdmx bruteLdmx(Matrix lambda, Matrix D, Matrix N, Matrix Z,
                                       Matrix mu, int Kopen) {
        int M = D.getNumRows();
        int R = D.getNumCols();
        double[][] a = new double[M][R];
        for (int r = 0; r < R; r++) {
            for (int i = 0; i < M; i++) {
                a[i][r] = Double.isInfinite(N.get(0, r))
                        ? lambda.get(0, r) * D.get(i, r)
                        : D.get(i, r);
            }
        }
        List<List<int[]>> alloc = new ArrayList<List<int[]>>();
        for (int r = 0; r < R; r++) {
            if (Double.isInfinite(N.get(0, r))) {
                alloc.add(PfqnSensMvaTest.compositionsLeq(Kopen, M));
            } else {
                alloc.add(PfqnSensMvaTest.compositionsLeq((int) N.get(0, r), M));
            }
        }

        List<int[][]> states = new ArrayList<int[][]>();
        List<Double> weights = new ArrayList<Double>();
        int[] idx = new int[R];
        while (true) {
            int[][] nir = new int[M][R];
            for (int r = 0; r < R; r++) {
                int[] al = alloc.get(r).get(idx[r]);
                for (int i = 0; i < M; i++) {
                    nir[i][r] = al[i];
                }
            }
            double lw = 0.0;
            boolean ok = true;
            for (int i = 0; i < M && ok; i++) {
                int ni = 0;
                for (int r = 0; r < R; r++) {
                    ni += nir[i][r];
                }
                lw += PfqnSensMvaTest.logFactorial(ni);
                for (int j = 1; j <= ni; j++) {
                    lw -= Math.log(muAt(mu, i, j));
                }
                for (int r = 0; r < R; r++) {
                    if (nir[i][r] > 0) {
                        if (a[i][r] <= 0) {
                            ok = false;
                            break;
                        }
                        lw += nir[i][r] * Math.log(a[i][r]) - PfqnSensMvaTest.logFactorial(nir[i][r]);
                    }
                }
            }
            if (ok) {
                for (int c = 0; c < R; c++) {
                    if (Double.isInfinite(N.get(0, c))) {
                        continue;
                    }
                    int sum = 0;
                    for (int i = 0; i < M; i++) {
                        sum += nir[i][c];
                    }
                    int n0c = (int) N.get(0, c) - sum;
                    if (n0c > 0) {
                        if (Z.get(0, c) <= 0) {
                            ok = false;
                            break;
                        }
                        lw += n0c * Math.log(Z.get(0, c)) - PfqnSensMvaTest.logFactorial(n0c);
                    }
                }
            }
            states.add(nir);
            weights.add(Double.valueOf(ok ? Math.exp(lw) : 0.0));

            int r = R - 1;
            while (r >= 0) {
                idx[r]++;
                if (idx[r] < alloc.get(r).size()) {
                    break;
                }
                idx[r] = 0;
                r--;
            }
            if (r < 0) {
                break;
            }
        }

        int K = states.size();
        double[] w = new double[K];
        double tot = 0.0;
        for (int k = 0; k < K; k++) {
            w[k] = weights.get(k).doubleValue();
            tot += w[k];
        }
        for (int k = 0; k < K; k++) {
            w[k] /= tot;
        }

        BruteLdmx out = new BruteLdmx();
        out.Q = new Matrix(M, R);
        for (int k = 0; k < K; k++) {
            int[][] nir = states.get(k);
            for (int i = 0; i < M; i++) {
                for (int r = 0; r < R; r++) {
                    out.Q.set(i, r, out.Q.get(i, r) + w[k] * nir[i][r]);
                }
            }
        }
        out.QCov = new Matrix[M];
        for (int i = 0; i < M; i++) {
            out.QCov[i] = new Matrix(R, R);
            for (int r = 0; r < R; r++) {
                for (int s = 0; s < R; s++) {
                    double m2 = 0.0;
                    for (int k = 0; k < K; k++) {
                        int[][] nir = states.get(k);
                        m2 += w[k] * nir[i][r] * nir[i][s];
                    }
                    out.QCov[i].set(r, s, m2 - out.Q.get(i, r) * out.Q.get(i, s));
                }
            }
        }
        return out;
    }

    /** Limited load dependence: the rate saturates at its last tabulated value. */
    private static double muAt(Matrix mu, int i, int j) {
        return j <= mu.getNumCols() ? mu.get(i, j - 1) : mu.get(i, mu.getNumCols() - 1);
    }
}
