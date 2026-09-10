/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api;

import jline.api.pfqn.mva.Pfqn_momlin;
import jline.api.pfqn.mva.Pfqn_mva;
import jline.api.pfqn.sens.Pfqn_sens;
import jline.io.Ret;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Validates the CoMoM-backed dispatch of {@link Pfqn_sens#pfqn_sens} for the
 * M=1 repairman model (base measures vs {@link Pfqn_mva#pfqn_mva}, Jacobian vs
 * central finite differences), the exact queue-length second-moment by-product
 * ({@code QVar}/{@code QCov}), and the self-consistency of the moment
 * linearizer {@link Pfqn_momlin#pfqn_momlin} (its analytic demand-derivatives
 * against finite differences of its own Schweitzer-Bard means).
 */
public class PfqnSensComomTest {

    private static Matrix row(double... v) {
        Matrix m = new Matrix(1, v.length);
        for (int i = 0; i < v.length; i++) {
            m.set(0, i, v[i]);
        }
        return m;
    }

    private static double relerr(double a, double b) {
        return Math.abs(a - b) / Math.max(1.0, Math.max(Math.abs(a), Math.abs(b)));
    }

    @Test
    public void comomPathMatchesMvaAndFiniteDifferences() {
        double[][] Ls = {{0.6, 0.3, 0.5}, {1.0, 0.5}, {0.8, 0.55, 0.4}};
        double[][] Ns = {{3, 2, 2}, {4, 3}, {5, 4, 3}};
        double[][] Zs = {{0.4, 0.8, 0.2}, {0.5, 0.9}, {0.6, 0.3, 0.9}};

        double maxBase = 0.0;
        double maxJac = 0.0;
        for (int t = 0; t < Ls.length; t++) {
            Matrix L = row(Ls[t]);
            Matrix N = row(Ns[t]);
            Matrix Z = row(Zs[t]);
            int R = L.getNumCols();
            Ret.pfqnSens s = Pfqn_sens.pfqn_sens(L, N, Z);

            Ret.pfqnMVA mva = Pfqn_mva.pfqn_mva(L, N, Z);
            for (int r = 0; r < R; r++) {
                maxBase = Math.max(maxBase, relerr(s.X.get(0, r), mva.X.get(0, r)));
                maxBase = Math.max(maxBase, relerr(s.Q.get(0, r), mva.Q.get(0, r)));
                maxBase = Math.max(maxBase, relerr(s.U.get(0, r), mva.U.get(0, r)));
                maxBase = Math.max(maxBase, relerr(s.R.get(0, r), mva.R.get(0, r)));
            }

            // Jacobian vs central finite differences of pfqn_mva
            for (int p = 0; p < s.paramType.length; p++) {
                double h = 1e-6;
                Matrix Lp = L.copy();
                Matrix Lm = L.copy();
                Matrix Zp = Z.copy();
                Matrix Zm = Z.copy();
                double base;
                if (s.paramType[p] == 0) {
                    base = L.get(s.paramStation[p], s.paramClass[p]);
                    double hh = h * Math.max(1.0, Math.abs(base));
                    Lp.set(s.paramStation[p], s.paramClass[p], base + hh);
                    Lm.set(s.paramStation[p], s.paramClass[p], base - hh);
                    h = hh;
                } else {
                    base = Z.get(0, s.paramClass[p]);
                    double hh = h * Math.max(1.0, Math.abs(base));
                    Zp.set(0, s.paramClass[p], base + hh);
                    Zm.set(0, s.paramClass[p], base - hh);
                    h = hh;
                }
                Ret.pfqnMVA mp = Pfqn_mva.pfqn_mva(Lp, N, Zp);
                Ret.pfqnMVA mm = Pfqn_mva.pfqn_mva(Lm, N, Zm);
                for (int r = 0; r < R; r++) {
                    double fdX = (mp.X.get(0, r) - mm.X.get(0, r)) / (2 * h);
                    double fdQ = (mp.Q.get(0, r) - mm.Q.get(0, r)) / (2 * h);
                    double fdU = (mp.U.get(0, r) - mm.U.get(0, r)) / (2 * h);
                    double fdR = (mp.R.get(0, r) - mm.R.get(0, r)) / (2 * h);
                    maxJac = Math.max(maxJac, relerr(s.dX.get(r, p), fdX));
                    maxJac = Math.max(maxJac, relerr(s.dQ[p].get(0, r), fdQ));
                    maxJac = Math.max(maxJac, relerr(s.dU[p].get(0, r), fdU));
                    maxJac = Math.max(maxJac, relerr(s.dR[p].get(0, r), fdR));
                }
                // QCov must equal D_{j,s} dQ_{i,r}/dD_{j,s}
                if (s.paramType[p] == 0) {
                    int js = s.paramClass[p];
                    for (int r = 0; r < R; r++) {
                        double cov = L.get(0, js) * s.dQ[p].get(0, r);
                        maxJac = Math.max(maxJac, relerr(s.QCov[0][r].get(0, js), cov));
                    }
                }
            }
        }
        assertTrue(maxBase < 1e-9, "base measures vs pfqn_mva: " + maxBase);
        assertTrue(maxJac < 1e-4, "Jacobian vs finite differences: " + maxJac);
    }

    @Test
    public void momentLinearizerIsSelfConsistent() {
        Matrix L = new Matrix(2, 3);
        double[][] Ld = {{0.8, 0.5, 0.4}, {0.3, 0.6, 0.5}};
        for (int i = 0; i < 2; i++) {
            for (int r = 0; r < 3; r++) {
                L.set(i, r, Ld[i][r]);
            }
        }
        Matrix N = row(3, 3, 2);
        Matrix Z = row(0.5, 0.3, 0.7);
        Pfqn_momlin.MomlinResult res = Pfqn_momlin.pfqn_momlin(L, N, Z);

        double maxErr = 0.0;
        for (int j = 0; j < 2; j++) {
            for (int s = 0; s < 3; s++) {
                double h = 1e-6 * Math.max(1.0, L.get(j, s));
                Matrix Lp = L.copy();
                Matrix Lm = L.copy();
                Lp.set(j, s, L.get(j, s) + h);
                Lm.set(j, s, L.get(j, s) - h);
                Pfqn_momlin.MomlinResult rp = Pfqn_momlin.pfqn_momlin(Lp, N, Z);
                Pfqn_momlin.MomlinResult rm = Pfqn_momlin.pfqn_momlin(Lm, N, Z);
                for (int i = 0; i < 2; i++) {
                    for (int r = 0; r < 3; r++) {
                        double fd = (rp.Q.get(i, r) - rm.Q.get(i, r)) / (2 * h);
                        double an = res.dQ[j][s].get(i, r);
                        maxErr = Math.max(maxErr, Math.abs(fd - an) / Math.max(1e-8, Math.abs(fd)));
                    }
                }
            }
        }
        assertTrue(maxErr < 1e-4, "momlin analytic dQ vs finite differences: " + maxErr);
    }
}
