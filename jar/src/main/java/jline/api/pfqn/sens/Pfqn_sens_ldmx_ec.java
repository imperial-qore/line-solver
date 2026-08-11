/**
 * Effective capacity terms of the mixed load-dependent MVA together with their
 * exact derivatives with respect to the open-class load.
 *
 * <p>Mirrors the MATLAB reference {@code pfqn_sens_ldmx_ec.m}.</p>
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.sens;

import org.apache.commons.math3.util.FastMath;

import jline.io.Ret;
import jline.util.matrix.Matrix;

public final class Pfqn_sens_ldmx_ec {
    private Pfqn_sens_ldmx_ec() {}

    /**
     * Computes the effective capacity terms EC, E and Eprime of the mixed
     * load-dependent MVA of Bruell-Balbo-Afshari, exactly as
     * {@link jline.api.pfqn.ld.Pfqn_ldmx_ec} does, and additionally their exact
     * analytic derivatives with respect to the open-class load Lo(i) of each
     * station.
     *
     * <p>Lo(i) = sum_r lambda(r)*D(i,r) is the only channel through which a
     * service demand enters E, Eprime and EC: the load-dependent rates mu are
     * independent of the demands. Station i's terms depend on Lo(i) alone, so a
     * single derivative per station is enough, and the chain rule then yields the
     * derivative with respect to any demand-scaling parameter. This is the
     * factorization behind equations (19), (21) and (24)-(31) of the reference,
     * which are reproduced here term by term.</p>
     *
     * <p>Reference: I. F. Akyildiz and J. C. Strelen, "Moment Analysis for
     * Load-Dependent Mixed Product Form Queueing Networks", IEEE Trans.
     * Communications 39(6):828-832, 1991.</p>
     *
     * <p>The column index n of E, Eprime, dE and dEprime carries the value of the
     * corresponding term at population n, i.e. it matches the MATLAB column 1+n;
     * the column index n of EC and dEC carries the term at population n+1, i.e.
     * the MATLAB column 1+n. This is the same convention used by
     * {@link jline.api.pfqn.ld.Pfqn_ldmx_ec}, so the two are interchangeable.</p>
     *
     * @param lambda arrival rate vector (1 x R); zero for closed classes
     * @param D      service demand matrix (M x R)
     * @param mu     load-dependent rate matrix (M x Nt), limited load dependence
     * @return EC, E, Eprime, Lo and their derivatives with respect to Lo
     */
    public static Ret.pfqnSensLdmxEc pfqn_sens_ldmx_ec(Matrix lambda, Matrix D, Matrix mu) {
        int M = mu.getNumRows();
        int Nt = mu.getNumCols();
        int R = D.getNumCols();

        Matrix Lo = new Matrix(M, 1);
        for (int ist = 0; ist < M; ist++) {
            double acc = 0.0;
            for (int r = 0; r < R; r++) {
                acc += lambda.get(r) * D.get(ist, r);
            }
            Lo.set(ist, 0, acc);
        }

        // limited load dependence level, kept in the 1-based convention of the
        // MATLAB reference: b(i) is the first level at which mu saturates
        int[] b = new int[M];
        int bmax = 0;
        for (int ist = 0; ist < M; ist++) {
            int idx = 0;
            while (idx < Nt && mu.get(ist, idx) != mu.get(ist, Nt - 1)) {
                idx++;
            }
            b[ist] = idx + 1;
            bmax = Math.max(bmax, b[ist]);
        }

        // C(i,j) = 1/mu(i,j) for j = 1..Nt+1+max(b), the rate saturating past Nt.
        // Cv[i][j] holds the MATLAB entry C(i,1+j).
        int Cn = Nt + 1 + bmax;
        double[][] Cv = new double[M][Cn];
        for (int ist = 0; ist < M; ist++) {
            for (int j = 0; j < Cn; j++) {
                Cv[ist][j] = 1.0 / (j < Nt ? mu.get(ist, j) : mu.get(ist, Nt - 1));
            }
        }

        Matrix EC = new Matrix(M, Nt);
        Matrix E = new Matrix(M, 1 + Nt);
        Matrix Eprime = new Matrix(M, 1 + Nt);
        Matrix dEC = new Matrix(M, Nt);
        Matrix dE = new Matrix(M, 1 + Nt);
        Matrix dEprime = new Matrix(M, 1 + Nt);

        for (int ist = 0; ist < M; ist++) {
            int bi = b[ist];
            double Cb = Cv[ist][bi - 1];
            double Loi = Lo.get(ist, 0);
            double den = 1.0 - Loi * Cb;   // geometric tail factor of the limited load dependence

            double[] E1 = new double[1 + Nt];
            double[] dE1 = new double[1 + Nt];
            double[] E2 = new double[1 + Nt];
            double[] dE2 = new double[1 + Nt];
            double[] E3 = new double[1 + Nt];
            double[] dE3 = new double[1 + Nt];
            double[] E2prime = new double[1 + Nt];
            double[] dE2prime = new double[1 + Nt];
            int W = 1 + Math.max(0, bi - 2);
            double[][] F2 = new double[1 + Nt][W];
            double[][] dF2 = new double[1 + Nt][W];
            double[][] F3 = new double[1 + Nt][W];
            double[][] dF3 = new double[1 + Nt][W];
            double[][] F2prime = new double[1 + Nt][W];
            double[][] dF2prime = new double[1 + Nt][W];

            for (int n = 0; n <= Nt; n++) {
                if (n >= bi) {
                    // E(n) = 1/den^(n+1)  =>  dE/dLo = (n+1)*Cb/den^(n+2)
                    E.set(ist, n, 1.0 / FastMath.pow(den, n + 1));
                    dE.set(ist, n, (n + 1) * Cb / FastMath.pow(den, n + 2));
                    Eprime.set(ist, n, Cb * E.get(ist, n));
                    dEprime.set(ist, n, Cb * dE.get(ist, n));
                } else { // n <= bi-1
                    // E1 and its derivative, eq. (25)-(26)
                    if (n == 0) {
                        E1[0] = 1.0 / den;
                        dE1[0] = Cb / (den * den);
                        for (int j = 1; j <= bi - 1; j++) {
                            E1[0] = E1[0] * Cv[ist][j - 1] / Cb;
                            dE1[0] = dE1[0] * Cv[ist][j - 1] / Cb;
                        }
                    } else {
                        double fac = Cb / Cv[ist][n - 1];
                        E1[n] = (1.0 / den) * fac * E1[n - 1];
                        dE1[n] = (Cb / (den * den)) * fac * E1[n - 1] + (1.0 / den) * fac * dE1[n - 1];
                    }

                    // F2 and its derivative, eq. (27)-(28)
                    for (int n0 = 0; n0 <= bi - 2; n0++) {
                        if (n0 == 0) {
                            F2[n][0] = 1.0;
                            dF2[n][0] = 0.0;
                        } else {
                            double coef = ((double) (n + n0)) / n0 * Cv[ist][n + n0 - 1];
                            F2[n][n0] = coef * Loi * F2[n][n0 - 1];
                            dF2[n][n0] = coef * (F2[n][n0 - 1] + Loi * dF2[n][n0 - 1]);
                        }
                    }
                    double s2 = 0.0;
                    double ds2 = 0.0;
                    for (int n0 = 0; n0 <= bi - 2; n0++) {
                        s2 += F2[n][n0];
                        ds2 += dF2[n][n0];
                    }
                    E2[n] = s2;
                    dE2[n] = ds2;

                    // F3 and its derivative, eq. (29)-(30)
                    for (int n0 = 0; n0 <= bi - 2; n0++) {
                        if (n == 0 && n0 == 0) {
                            F3[0][0] = 1.0;
                            for (int j = 1; j <= bi - 1; j++) {
                                F3[0][0] = F3[0][0] * Cv[ist][j - 1] / Cb;
                            }
                            dF3[0][0] = 0.0;
                        } else if (n > 0 && n0 == 0) {
                            double fac = Cb / Cv[ist][n - 1];
                            F3[n][0] = fac * F3[n - 1][0];
                            dF3[n][0] = fac * dF3[n - 1][0];
                        } else {
                            double coef = ((double) (n + n0)) / n0 * Cb;
                            F3[n][n0] = coef * Loi * F3[n][n0 - 1];
                            dF3[n][n0] = coef * (F3[n][n0 - 1] + Loi * dF3[n][n0 - 1]);
                        }
                    }
                    double s3 = 0.0;
                    double ds3 = 0.0;
                    for (int n0 = 0; n0 <= bi - 2; n0++) {
                        s3 += F3[n][n0];
                        ds3 += dF3[n][n0];
                    }
                    E3[n] = s3;
                    dE3[n] = ds3;

                    // F2prime and its derivative
                    for (int n0 = 0; n0 <= bi - 2; n0++) {
                        if (n0 == 0) {
                            F2prime[n][0] = Cv[ist][n];
                            dF2prime[n][0] = 0.0;
                        } else {
                            double coef = ((double) (n + n0)) / n0 * Cv[ist][n + n0];
                            F2prime[n][n0] = coef * Loi * F2prime[n][n0 - 1];
                            dF2prime[n][n0] = coef * (F2prime[n][n0 - 1] + Loi * dF2prime[n][n0 - 1]);
                        }
                    }
                    double s2p = 0.0;
                    double ds2p = 0.0;
                    for (int n0 = 0; n0 <= bi - 2; n0++) {
                        s2p += F2prime[n][n0];
                        ds2p += dF2prime[n][n0];
                    }
                    E2prime[n] = s2p;
                    dE2prime[n] = ds2p;

                    // E = E1 + E2 - E3, eq. (23)-(24)
                    E.set(ist, n, E1[n] + E2[n] - E3[n]);
                    dE.set(ist, n, dE1[n] + dE2[n] - dE3[n]);
                    if (n < bi - 1) {
                        Eprime.set(ist, n, Cb * E1[n] + E2prime[n] - Cb * E3[n]);
                        dEprime.set(ist, n, Cb * dE1[n] + dE2prime[n] - Cb * dE3[n]);
                    } else { // n >= bi-1
                        Eprime.set(ist, n, Cb * E.get(ist, n));
                        dEprime.set(ist, n, Cb * dE.get(ist, n));
                    }
                }
            }

            // EC(n) = C(n)*E(n)/E(n-1), eq. (19); quotient rule for the derivative
            for (int n = 1; n <= Nt; n++) {
                double En = E.get(ist, n);
                double Enm = E.get(ist, n - 1);
                EC.set(ist, n - 1, Cv[ist][n - 1] * En / Enm);
                dEC.set(ist, n - 1, Cv[ist][n - 1]
                        * (dE.get(ist, n) * Enm - En * dE.get(ist, n - 1)) / (Enm * Enm));
            }
        }

        return new Ret.pfqnSensLdmxEc(EC, E, Eprime, Lo, dEC, dE, dEprime);
    }
}
