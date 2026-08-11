/**
 * @file Srinivasan (1985) Successively Improving Bounds (SIB)
 *
 * Closed-form hierarchy of upper/lower bounds on cycle time and throughput for single-class
 * product-form closed networks of fixed-rate (and delay) stations, based on the power sums
 * S_i = sum_m rho_m^i. Level 1 is Thm 2.1; higher levels use Thms 3.5 (upper) and 3.6 (lower).
 * Only Z=0 is supported (delay needs the Section-3.2 demand substitution). Ported at parity
 * from MATLAB pfqn_sib.m.
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.mva;

import jline.util.matrix.Matrix;
import org.apache.commons.math3.util.FastMath;

public final class Pfqn_sib {
    private Pfqn_sib() {}

    /**
     * Level-`level` successively-improving throughput/cycle-time bounds.
     *
     * @param L     fixed-rate service demand vector (M x 1); delay demand goes in Z
     * @param N     total population (>= 2)
     * @param Z     think time (only Z=0 supported)
     * @param level bound level (>= 1); higher = tighter
     * @return {Xlo, Xhi, Wlo, Whi}
     */
    public static double[] pfqn_sib(Matrix L, double N, double Z, int level) {
        if (level < 1) level = 1;
        if (Z > 0) {
            throw new RuntimeException("pfqn_sib supports Z=0 only (delay needs the "
                + "Section-3.2 demand substitution, not yet implemented).");
        }
        int M = L.length();
        double Lsum = L.elementSum();
        double[] rho = new double[M];
        double rhoU = 0.0;
        for (int m = 0; m < M; m++) { rho[m] = L.get(m) / Lsum; if (rho[m] > rhoU) rhoU = rho[m]; }
        int imax = level + 3;
        double[] S = new double[imax + 1];              // S[i] = sum rho^i, i=1..imax
        for (int i = 1; i <= imax; i++) {
            double s = 0.0;
            for (int m = 0; m < M; m++) s += FastMath.pow(rho[m], i);
            S[i] = s;
        }
        double S2 = S[2];
        double[] alpha = new double[level + 1];         // alpha[i] = alpha_i
        alpha[0] = S2;
        for (int i = 1; i <= level; i++) {
            double acc = 0.0;
            for (int j = 0; j <= i - 1; j++) acc += S[i + 1 - j] * alpha[j];
            alpha[i] = S[i + 2] - acc;
        }

        double NN = N - 1;
        double phiLo = (N - 1) * S2;
        double T1s2 = (N - 1) * rhoU - 1;
        double phiHi = 0.5 * (T1s2 + FastMath.sqrt(T1s2 * T1s2 + 4 * (N - 1) * S2));

        if (N >= 3) {
            double eta = (N - 2) / (N - 1);
            double T1u = (N - 2) * rhoU - 1;
            double su = sigma(NN, level - 1, rhoU, S, S2);
            double phiUn = 0.5 / eta * (T1u + FastMath.sqrt(Math.max(0, T1u * T1u + 4 * (N - 2) * (S2 - su))));
            phiHi = Math.min(phiHi, phiUn);

            double T1l = (N - 2) * S2 - 1;
            double bl = betaL(NN, level - 1, alpha, rhoU, S2);
            double phiLn = (T1l + FastMath.sqrt(Math.max(0, T1l * T1l + 4 * (N - 2) * (S2 + (N - 2) * bl)))) / (2 * level);
            phiLo = Math.max(phiLo, phiLn);
        }

        phiLo = Math.max(0, phiLo);
        if (phiHi < phiLo) phiHi = phiLo;

        double Wlo = Lsum * (1 + phiLo) + Z;
        double Whi = Lsum * (1 + phiHi) + Z;
        double Xlo = N / Whi;
        double Xhi = N / Wlo;
        return new double[]{Xlo, Xhi, Wlo, Whi};
    }

    /** Level-1 upper bound on phi(K) (Thm 3.5 with n=1 / eq. 3.18). */
    private static double phiU1(double K, double rhoU, double S2) {
        if (K <= 0) return 0.0;
        if (K == 1) return S2;
        double eta = (K - 1) / K;
        double T1 = (K - 1) * rhoU - 1;
        return 0.5 / eta * (T1 + FastMath.sqrt(T1 * T1 + 4 * (K - 1) * S2));
    }

    /** eq. (3.22c); NN plays the role of (N-1). */
    private static double sigma(double NN, int i, double rhoU, double[] S, double S2) {
        double s = 0.0;
        if (i <= 0) return 0.0;
        double Dbar = 1 + phiU1(NN - 2, rhoU, S2);
        double pnum = 1.0;
        for (int j = 1; j <= i; j++) {
            pnum = pnum * (NN - 1 - (j - 1));
            s += (rhoU * S[j + 1] - S[j + 2]) * pnum / FastMath.pow(Dbar, j);
        }
        return s;
    }

    /** eq. (3.23c); NN plays the role of (N-1). */
    private static double betaL(double NN, int i, double[] alpha, double rhoU, double S2) {
        double b = 0.0;
        if (i <= 0) return 0.0;
        for (int j = 1; j <= i - 1; j++) {
            double p = 1.0;
            for (int m = 2; m <= j; m++) p = p * (NN - m) / (1 + phiU1(NN - m, rhoU, S2));
            b += alpha[j] * p;
        }
        double p = 1.0;
        for (int m = 2; m <= i; m++) p = p * (NN - m) / (1 + phiU1(NN - m, rhoU, S2));
        double Nim2 = NN - 1 - i - 1;
        double corr = 1 + (alpha[i] / alpha[i - 1]) * Nim2 / (1 + Nim2 * alpha[0]);
        b += alpha[i] * p * corr;
        return b;
    }
}
