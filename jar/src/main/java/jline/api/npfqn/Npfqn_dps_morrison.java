/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api.npfqn;

import org.apache.commons.math3.special.Erf;
import org.apache.commons.math3.util.FastMath;

import jline.util.matrix.Matrix;

/**
 * Two-term heavy-usage asymptotic approximation for a closed queueing network with one
 * infinite-server (think) station and one discriminatory processor-sharing (DPS) station, after
 * J.A. Morrison, "Asymptotic analysis of a large closed queueing network with discriminatory
 * processor sharing", Queueing Systems 9 (1991) 191-214.
 *
 * <p>The network is NOT product-form, so nothing here computes a normalizing constant: the method
 * expands the GENERATING FUNCTION of the balance equations. The substitution
 * {@code P(n) = <w,n> f(n)} clears the DPS denominator and turns the balance recursion into a
 * linear PDE with affine coefficients (Morrison eq. 2.5); rescaling {@code z = 1 - xi/sqrt(N)} and
 * expanding in powers of {@code N^(-1/2)} leaves a degenerate leading operator whose kernel is the
 * functions of the similarity variable eta, and the solvability condition along its characteristic
 * gives an ODE for the amplitude (eq. 2.20). RESULT 1 (eq. 4.11) and RESULT 2 (eq. 4.17) are the
 * two-term approximations returned here.</p>
 *
 * <p>Scaling. Morrison writes {@code K_j = N b_j} and {@code lambda_j = N r_j g_j} with N large and
 * usage {@code rho = sum_j b_j/g_j = 1 - a/sqrt(N)}. N is bookkeeping only: b, g and a all move
 * with it and the approximation is invariant, so this routine fixes N = 1, i.e. {@code b = N_pop},
 * {@code g = Z./S}, {@code r = 1./Z} and {@code a = 1 - rho}. Accuracy is governed by the PHYSICAL
 * regime -- large populations with rho near 1 -- and not by any choice made here. {@code rho > 1}
 * is admissible: it is the saturated regime of Morrison's appendix A.</p>
 *
 * <p>This is the Java port of the MATLAB {@code npfqn_dps_morrison}; the two agree to double
 * precision on the parity tests.</p>
 */
public class Npfqn_dps_morrison {

    /** Mean queue lengths, sojourn times and throughputs, with the intermediate constants. */
    public static class Result {
        /** Mean number of class-k jobs at the DPS station, 1 x K. */
        public Matrix Q;
        /** Mean class-k sojourn time per visit to the DPS station, 1 x K. */
        public Matrix R;
        /** Per-class throughput, 1 x K. */
        public Matrix X;
        /** Leading-order (one-term) queue lengths, 1 x K. */
        public Matrix Qlead;
        /** Leading-order (one-term) sojourn times, 1 x K. */
        public Matrix Rlead;
        /** Usage sum_j b_j/g_j. */
        public double rho;
        /** Heavy-usage parameter a = 1 - rho at the N = 1 scale. */
        public double a;
        /** Morrison's constants of eqs. (2.18), (2.19), (3.11), (3.13)-(3.15), (3.19)-(3.21). */
        public double cB, cC, cD, cH, cI, cJ, cK, cL, cM, cQ, delta, cR, cS, cU, cA, cV;
        /** The vector sigma of eq. (4.12), 1 x K. */
        public Matrix sigma;
        /** W_0..W_4 of eq. (3.23). */
        public double[] W;
    }

    /**
     * Evaluates Morrison's two-term approximation.
     *
     * @param N per-class populations, 1 x K, finite and positive
     * @param Z per-class mean think times, 1 x K, finite and positive
     * @param S per-class mean DPS service times, 1 x K, finite and positive
     * @param w per-class DPS weights, 1 x K, finite and positive
     * @return the mean performance measures and the intermediate constants
     */
    public static Result npfqn_dps_morrison(Matrix N, Matrix Z, Matrix S, Matrix w) {
        int p = N.length();
        if (Z.length() != p || S.length() != p || w.length() != p) {
            throw new RuntimeException("N, Z, S and w must have the same number of classes.");
        }
        double[] b = new double[p];
        double[] z = new double[p];
        double[] s = new double[p];
        double[] wt = new double[p];
        for (int i = 0; i < p; i++) {
            b[i] = N.get(i);
            z[i] = Z.get(i);
            s[i] = S.get(i);
            wt[i] = w.get(i);
            if (!Double.isFinite(b[i]) || b[i] <= 0) {
                throw new RuntimeException("The Morrison approximation requires finite positive class "
                        + "populations (closed classes only).");
            }
            if (!Double.isFinite(z[i]) || z[i] <= 0 || !Double.isFinite(s[i]) || s[i] <= 0) {
                throw new RuntimeException("Think times Z and DPS service times S must be finite and positive.");
            }
            if (!Double.isFinite(wt[i]) || wt[i] <= 0) {
                throw new RuntimeException("DPS weights must be finite and positive.");
            }
        }

        // Morrison's parameters at the bookkeeping scale N = 1
        double[] r = new double[p];
        double[] g = new double[p];
        double rho = 0;
        for (int i = 0; i < p; i++) {
            r[i] = 1.0 / z[i];
            g[i] = z[i] / s[i];
            rho += b[i] / g[i];
        }
        double a = 1.0 - rho;

        // constants, eqs. (2.18), (2.19), (3.11), (3.13)-(3.15)
        double cB = 0, cC = 0, cD = 0, cH = 0, cI = 0, cJ = 0, cK = 0, cL = 0, cM = 0, cQ = 0;
        for (int i = 0; i < p; i++) {
            double g2 = g[i] * g[i], g3 = g2 * g[i], w2 = wt[i] * wt[i], r2 = r[i] * r[i];
            cB += b[i] / (r[i] * g2 * wt[i]);
            cC += b[i] / (g2 * wt[i]);
            cD += b[i] / (r[i] * g2);
            cH += b[i] / (r2 * g3 * wt[i]);
            cI += b[i] / (r2 * g3 * w2);
            cJ += b[i] / (r[i] * g3 * w2);
            cK += b[i] / (g3 * w2);
            cL += b[i] / (r[i] * g3 * wt[i]);
            cM += b[i] / (r2 * g3);
            cQ += b[i] / g2;
        }

        // sigma: eq. (4.12) with the normalization (4.13). The p equations have rank p-1 (Morrison
        // p.197), so the last one -- implied by the others -- is REPLACED by (4.13), giving a
        // square nonsingular system. All four codebases use this same scheme so their sigma agree.
        double[][] A = new double[p][p];
        double[] rhs = new double[p];
        for (int i = 0; i < p; i++) {
            A[i][i] += rho;
            for (int j = 0; j < p; j++) {
                double den = r[i] * g[i] * wt[i] + r[j] * g[j] * wt[j];
                A[i][i] -= wt[j] * b[j] * r[j] / den;
                A[i][j] -= wt[j] * b[i] * r[i] / den;
            }
            rhs[i] = rho * (b[i] / g[i]) * (cD / (cB * wt[i]) - 1);
        }
        for (int j = 0; j < p; j++) {
            A[p - 1][j] = 1.0 / (r[j] * g[j]);
        }
        rhs[p - 1] = 0.0;
        double[] sigma = solveSquare(A, rhs);

        // alpha from eq. (4.9), then delta of eq. (3.15)
        double delta = 0;
        double[] alpha = new double[p];
        for (int i = 0; i < p; i++) {
            alpha[i] = sigma[i] - (b[i] / g[i]) * (cD / (cB * wt[i]) - 1);
            delta += alpha[i] / g[i];
        }

        // eqs. (3.19)-(3.21)
        double cR = 3 * (cB * cL - cD * cJ) / (cB * cD);
        double cS = (2 * cB * (cD * cH - cB * cM) - cD * (cD * cI - cB * cH)) / (2 * cB * cB * cD * cD);
        double cU = (cQ - cC * cD / cB - delta) / rho - cD * cR / cB + (a * a - cC * cD / cB) * cS;
        double cA = cS * cC * cC + cR * cC - cK;
        double cV = cR + 2 * cS * cC;

        double[] W = wm(cB, cC, cD, a);

        // RESULT 1 (4.11) and RESULT 2 (4.17), at sqrt(N) = 1. NOTE the numerator bracket carries
        // U*W2: eq. (4.10) of the paper misprints it as U*W1, but (4.7), (4.11), (A6) and (B2) all
        // agree on U*W2, and it is what the derivation from (4.4)-(4.9) gives.
        double eps = cB / cD;
        double num = W[1] - eps * (cA / 3 * W[4] + a / 2 * cV * W[3] + cU * W[2]);
        double den = W[0] - eps * (cA / 3 * W[3] + a / 2 * cV * W[2] + cU * W[1] + cS);
        if (den == 0 || !Double.isFinite(den)) {
            throw new RuntimeException("The Morrison expansion is degenerate for this model (vanishing "
                    + "denominator); the usage is too far from the moderately-heavy regime.");
        }

        Result res = new Result();
        res.Q = new Matrix(1, p);
        res.R = new Matrix(1, p);
        res.X = new Matrix(1, p);
        res.Qlead = new Matrix(1, p);
        res.Rlead = new Matrix(1, p);
        res.sigma = new Matrix(1, p);
        for (int j = 0; j < p; j++) {
            double gw = g[j] * wt[j];
            double qlead = b[j] * W[1] / (gw * W[0]);
            double q = b[j] * num / (gw * den)
                    - b[j] * W[2] / (g[j] * g[j] * wt[j] * wt[j] * W[0])
                    - sigma[j] / rho;
            double rlead = W[1] / (r[j] * gw * W[0]);
            double rr = num / (r[j] * gw * den)
                    + ((W[1] / W[0]) * (W[1] / W[0]) - W[2] / W[0]) / (r[j] * g[j] * g[j] * wt[j] * wt[j])
                    - sigma[j] / (rho * r[j] * b[j]);
            res.Q.set(0, j, q);
            res.R.set(0, j, rr);
            res.X.set(0, j, r[j] * (b[j] - q));
            res.Qlead.set(0, j, qlead);
            res.Rlead.set(0, j, rlead);
            res.sigma.set(0, j, sigma[j]);
        }
        res.rho = rho;
        res.a = a;
        res.cB = cB; res.cC = cC; res.cD = cD; res.cH = cH; res.cI = cI; res.cJ = cJ;
        res.cK = cK; res.cL = cL; res.cM = cM; res.cQ = cQ; res.delta = delta;
        res.cR = cR; res.cS = cS; res.cU = cU; res.cA = cA; res.cV = cV;
        res.W = W;
        return res;
    }

    /**
     * W_m of eq. (3.23), m = 0..4. Substituting z = sigma s with sigma = sqrt(D/(BC)) normalizes
     * the Gaussian to W_m(y) = (B/D)^2 sigma^(m+1) I_m(yh) with yh = y sqrt(B/(CD)) and
     * I_m = int_0^inf s^m exp(-s^2/2 - yh s) ds, so I_0 = sqrt(pi/2) erfcx(yh/sqrt(2)),
     * I_1 = 1 - yh I_0 and I_m = (m-1) I_{m-2} - yh I_{m-1}. The recursion cancels for large yh, so
     * a loss of positivity (the I_m are integrals of positive integrands) falls back to quadrature.
     */
    static double[] wm(double cB, double cC, double cD, double y) {
        final int mmax = 4;
        double sig = FastMath.sqrt(cD / (cB * cC));
        double yh = y * FastMath.sqrt(cB / (cC * cD));

        double[] Iv = new double[mmax + 1];
        Iv[0] = FastMath.sqrt(Math.PI / 2) * erfcx(yh / FastMath.sqrt(2));
        if (!Double.isFinite(Iv[0])) {
            throw new RuntimeException("The usage is so far above saturation (rho = " + (1 - y) + ") that "
                    + "the Morrison expansion overflows. This model is outside the moderately-heavy regime "
                    + "the approximation is derived for; use SolverFLD, SolverMVA or SolverCTMC.");
        }
        Iv[1] = 1 - yh * Iv[0];
        for (int m = 2; m <= mmax; m++) {
            Iv[m] = (m - 1) * Iv[m - 2] - yh * Iv[m - 1];
        }
        boolean positive = true;
        for (int m = 0; m <= mmax; m++) {
            if (!(Iv[m] > 0)) {
                positive = false;
            }
        }
        if (!positive) {
            for (int m = 0; m <= mmax; m++) {
                Iv[m] = quadI(m, yh);
            }
        }
        double[] W = new double[mmax + 1];
        double sp = sig;
        for (int m = 0; m <= mmax; m++) {
            W[m] = (cB / cD) * (cB / cD) * sp * Iv[m];
            sp *= sig;
        }
        return W;
    }

    /** I_m(yh) by composite Simpson on a truncation that follows the integrand's own scale. */
    private static double quadI(int m, double yh) {
        double hi = yh >= 1 ? Math.min(12.0, 40.0 / yh) : (yh >= 0 ? 12.0 : 12.0 - yh);
        int n = 4096;
        double h = hi / n;
        double sum = 0;
        for (int k = 0; k <= n; k++) {
            double s = k * h;
            double f = (m == 0 ? 1.0 : FastMath.pow(s, m)) * FastMath.exp(-s * s / 2 - yh * s);
            sum += (k == 0 || k == n) ? f : ((k % 2 == 1) ? 4 * f : 2 * f);
        }
        return sum * h / 3;
    }

    /**
     * Scaled complementary error function exp(x^2)*erfc(x), valid for either sign. The direct
     * product overflows past x ~ 26, where the asymptotic series is already exact to double
     * precision; for x &lt; 0 the reflection erfcx(x) = 2 exp(x^2) - erfcx(-x) is used, which
     * overflows below x ~ -26 and is reported as such by the caller.
     */
    static double erfcx(double x) {
        if (x < 0) {
            return 2.0 * FastMath.exp(x * x) - erfcx(-x);
        }
        if (x < 25.0) {
            return FastMath.exp(x * x) * Erf.erfc(x);
        }
        double y = 1.0 / (2.0 * x * x);
        double term = 1.0;
        double sum = 1.0;
        for (int k = 1; k <= 12; k++) {
            term *= -(2 * k - 1) * y;
            sum += term;
        }
        return sum / (x * FastMath.sqrt(Math.PI));
    }

    /** Gaussian elimination with partial pivoting on a small dense square system. */
    private static double[] solveSquare(double[][] A, double[] rhs) {
        int n = rhs.length;
        double[][] M = new double[n][n + 1];
        for (int i = 0; i < n; i++) {
            System.arraycopy(A[i], 0, M[i], 0, n);
            M[i][n] = rhs[i];
        }
        for (int c = 0; c < n; c++) {
            int piv = c;
            for (int i = c + 1; i < n; i++) {
                if (Math.abs(M[i][c]) > Math.abs(M[piv][c])) {
                    piv = i;
                }
            }
            if (Math.abs(M[piv][c]) < 1e-300) {
                throw new RuntimeException("Singular system while solving Morrison's sigma equations (4.12).");
            }
            double[] t = M[c]; M[c] = M[piv]; M[piv] = t;
            for (int i = c + 1; i < n; i++) {
                double f = M[i][c] / M[c][c];
                for (int j = c; j <= n; j++) {
                    M[i][j] -= f * M[c][j];
                }
            }
        }
        double[] x = new double[n];
        for (int i = n - 1; i >= 0; i--) {
            double acc = M[i][n];
            for (int j = i + 1; j < n; j++) {
                acc -= M[i][j] * x[j];
            }
            x[i] = acc / M[i][i];
        }
        return x;
    }
}
