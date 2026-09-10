package jline.api;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assertions.fail;

import org.junit.jupiter.api.Test;

import jline.api.pfqn.mva.Pfqn_nintmva;
import jline.api.pfqn.mva.Pfqn_tay;
import jline.api.pfqn.nc.Pfqn_ca;
import jline.api.pfqn.nc.Pfqn_dnc;
import jline.api.pfqn.nc.Pfqn_rgf;
import jline.api.pfqn.nc.Pfqn_rgfmc;
import jline.api.pfqn.sens.Pfqn_hst;
import jline.util.matrix.Matrix;

/**
 * Numerical validation of pfqn_rgf, pfqn_dnc, pfqn_nintmva, pfqn_tay and pfqn_hst against
 * their source papers and against the exact routes already in the API.
 */
public class PfqnNewMethodsTest {

    private static Matrix col(double... v) {
        Matrix m = new Matrix(v.length, 1);
        for (int i = 0; i < v.length; i++) {
            m.set(i, 0, v[i]);
        }
        return m;
    }

    private static Matrix row(double... v) {
        Matrix m = new Matrix(1, v.length);
        for (int i = 0; i < v.length; i++) {
            m.set(0, i, v[i]);
        }
        return m;
    }

    @Test
    public void rgfMatchesConvolution() {
        double[][] demands = {{0.5, 0.3, 0.2}, {1.0, 1.0, 1.0, 0.4}, {0.05, 0.05, 0.05, 0.05, 0.9}};
        double[] thinks = {0.0, 3.0, 12.5};
        for (int d = 0; d < demands.length; d++) {
            Matrix L = col(demands[d]);
            for (int n = 1; n <= 20; n++) {
                for (int z = 0; z < thinks.length; z++) {
                    double lgCa = Pfqn_ca.pfqn_ca(L, row(n), row(thinks[z])).lG;
                    double lgRgf = Pfqn_rgf.pfqn_rgf(L, n, thinks[z]).lG;
                    assertEquals(lgCa, lgRgf, 1e-9 * Math.max(1.0, Math.abs(lgCa)));
                }
            }
        }
    }

    @Test
    public void rgfReproducesCouryHarrisonTables() {
        // Coury-Harrison (1997) Sec. 4: IS terminals of load 1, three groups of m devices
        // with per-device loads 0.001/0.002/0.003, one CPU of load 0.004.
        int[] mv = {2, 10, 20};
        int[] kv = {20, 100, 200, 300};
        // Published exact throughputs; the m=20, k=300 entry is printed 227.52 in both
        // tables and is a dropped-zero typo for 227.052 (see line-gaps.md G37).
        double[][] pub = {{19.67, 97.83, 192.44, 249.67},
                          {18.75, 92.33, 178.49, 242.84},
                          {17.71, 86.444, 164.93, 227.052}};
        double[] p = {0.001, 0.002, 0.003};
        for (int a = 0; a < mv.length; a++) {
            int m = mv[a];
            double[] demands = new double[3 * m + 1];
            for (int g = 0; g < 3; g++) {
                for (int i = 0; i < m; i++) {
                    demands[g * m + i] = p[g];
                }
            }
            demands[3 * m] = 0.004;
            Matrix L = col(demands);
            for (int c = 0; c < kv.length; c++) {
                double lg1 = Pfqn_rgf.pfqn_rgf(L, kv[c], 1.0).lG;
                double lg0 = Pfqn_rgf.pfqn_rgf(L, kv[c] - 1, 1.0).lG;
                double T = Math.exp(lg0 - lg1);
                assertEquals(pub[a][c], T, 5e-4 * pub[a][c]);
            }
        }
    }

    @Test
    public void rgfmcMatchesConvolutionWithThinkTimes() {
        // Multiclass RGF: iterated residues (Harrison-Coury 2002 Thm 1) with the delay
        // carried by the Bertozzi-McKenna truncation, which neither RGF paper has.
        double[][] L2 = {{0.9, 0.4}, {0.6, 0.7}, {0.3, 1.1}};
        double[][] L3 = {{0.2, 0.4, 1.2}, {0.7, 1.3, 0.3}, {1.4, 0.2, 0.7}};
        double[][] Z2 = {{0.0, 0.0}, {2.0, 0.0}, {2.0, 1.5}};
        double[][] Z3 = {{0.0, 0.0, 0.0}, {1.3, 1.3, 1.3}};
        for (int t = 0; t < Z2.length; t++) {
            assertRgfmc(L2, new double[]{3, 4}, Z2[t]);
        }
        for (int t = 0; t < Z3.length; t++) {
            assertRgfmc(L3, new double[]{2, 3, 2}, Z3[t]);
        }
    }

    private void assertRgfmc(double[][] Ld, double[] N, double[] Z) {
        Matrix L = new Matrix(Ld.length, Ld[0].length);
        for (int i = 0; i < Ld.length; i++) {
            for (int r = 0; r < Ld[0].length; r++) {
                L.set(i, r, Ld[i][r]);
            }
        }
        Matrix Nm = new Matrix(1, N.length);
        Matrix Zm = new Matrix(1, Z.length);
        for (int r = 0; r < N.length; r++) {
            Nm.set(0, r, N[r]);
            Zm.set(0, r, Z[r]);
        }
        double lgCa = Pfqn_ca.pfqn_ca(L, Nm, Zm).lG;
        assertEquals(lgCa, Pfqn_rgfmc.pfqn_rgfmc(L, N, Z).lG, 1e-7);
    }

    @Test
    public void rgfmcIsInsensitiveToTheEliminatedPopulation() {
        // Harrison-Lee sec. 4: a class removed by residues enters only as a pole ORDER,
        // so its population is nearly free -- while it carries no think time, since the
        // truncation reintroduces exactly that dependence.
        double[][] Ld = {{0.9, 0.4}, {0.6, 0.7}, {0.3, 1.1}, {0.5, 0.2}};
        int[] n2 = {20, 200, 2000};
        for (int a = 0; a < n2.length; a++) {
            assertRgfmc(Ld, new double[]{6, n2[a]}, new double[]{1.5, 0.0});
        }
    }

    @Test
    public void rgfmcRefusesRatherThanReturningAWrongConstant() {
        // Near-coincident loads over the eliminated class make the alternating residue
        // sum cancel; the guard must refuse, not answer.
        double[][] Ld = {{0.81733524, 1.24343101, 0.86870538},
                         {1.4732791, 0.38631325, 0.87522251}};
        Matrix L = new Matrix(2, 3);
        for (int i = 0; i < 2; i++) {
            for (int r = 0; r < 3; r++) {
                L.set(i, r, Ld[i][r]);
            }
        }
        try {
            Pfqn_rgfmc.pfqn_rgfmc(L, new double[]{2, 3, 2}, new double[]{0, 0, 0},
                    Pfqn_rgfmc.DEF_TOL, Pfqn_rgfmc.DEF_MAXTERMS, 1.0);
            fail("pfqn_rgfmc must refuse when the residue sum has cancelled away");
        } catch (IllegalStateException expected) {
            assertTrue(expected.getMessage().contains("cancelled"));
        }
    }

    @Test
    public void dncIsExactAtIntegerPopulations() {
        Matrix distinct = col(1 / 20.3, 0.6114 / 10.1, 0.3886 / 1.2);
        Matrix repeated = col(0.5, 0.5, 0.5, 0.2, 0.2, 0.9);
        Matrix[] cases = {distinct, repeated};
        for (int d = 0; d < cases.length; d++) {
            for (int n = 1; n <= 12; n++) {
                double lgCa = Pfqn_ca.pfqn_ca(cases[d], row(n), row(0.0)).lG;
                Pfqn_dnc.Result r = Pfqn_dnc.pfqn_dnc(cases[d], n);
                assertEquals(lgCa, r.lG, 1e-8 * Math.max(1.0, Math.abs(lgCa)));
            }
        }
    }

    @Test
    public void dncAndNintmvaReproduceDowdyGordonFigure5() {
        // Dowdy-Gordon (1984) Sec. 5: DMPavg = 2.5 on the composite one-hour model.
        // Published: DNC 2.999, aMVA 2.991.
        Matrix L = col(1 / 20.3, 0.6114 / 10.1, 0.3886 / 1.2);
        assertEquals(2.999, Pfqn_dnc.pfqn_dnc(L, 2.5).X, 1e-3);
        assertEquals(2.991, Pfqn_nintmva.pfqn_nintmva(L, 2.5, 0.0).X, 1e-3);
        // The interpolant is monotone between the integral points it reproduces.
        double prev = -1;
        for (double n = 2.0; n <= 3.0 + 1e-12; n += 0.25) {
            double x = Pfqn_dnc.pfqn_dnc(L, n).X;
            assertTrue(x > prev, "DNC throughput must increase in the population");
            prev = x;
        }
    }

    @Test
    public void nintmvaMatchesExactMvaAtIntegerPopulations() {
        Matrix L = col(1 / 20.3, 0.6114 / 10.1, 0.3886 / 1.2);
        double[] thinks = {0.0, 5.0};
        for (int z = 0; z < thinks.length; z++) {
            for (int n = 1; n <= 12; n++) {
                // Exact single-class MVA by the same recursion from an empty network.
                double[] q = new double[L.length()];
                double x = 0;
                for (int k = 1; k <= n; k++) {
                    double sumR = 0;
                    double[] r = new double[L.length()];
                    for (int i = 0; i < L.length(); i++) {
                        r[i] = L.get(i) * (1 + q[i]);
                        sumR += r[i];
                    }
                    x = k / (thinks[z] + sumR);
                    for (int i = 0; i < L.length(); i++) {
                        q[i] = x * r[i];
                    }
                }
                Pfqn_nintmva.Result res = Pfqn_nintmva.pfqn_nintmva(L, n, thinks[z]);
                assertEquals(x, res.X, 1e-12);
                for (int i = 0; i < L.length(); i++) {
                    assertEquals(q[i], res.Q.get(i), 1e-12);
                }
            }
        }
    }

    @Test
    public void tayReproducesSurveyExample4() {
        // Tay's Example 4 (CN82 Fig. 2 network) at N = (8,1): the survey's Tay row is
        // X1 = 0.654, X2 = 0.336, L11 = 5.234, L12 = 2.766, L22 = 1.000.
        Matrix L = new Matrix(1, 2);
        L.set(0, 0, 1.0);
        L.set(0, 1, 1.0);
        Pfqn_tay.Result r = Pfqn_tay.pfqn_tay(L, row(8, 1), row(8, 0));
        assertEquals(0.654, r.X.get(0), 1e-3);
        assertEquals(0.336, r.X.get(1), 1e-3);
        assertEquals(5.234, r.X.get(0) * 8, 1e-3);
        assertEquals(2.766, r.Q.get(0, 0), 1e-3);
        assertEquals(1.000, r.Q.get(0, 1), 1e-3);
        // Table B: the arrival-instant queue lengths, not re-solves at N - e_r.
        assertEquals(2.23, r.Qarr[0].get(0, 0), 5e-3);
        assertEquals(1.98, r.Qarr[1].get(0, 0), 5e-3);
        assertEquals(1.00, r.Qarr[0].get(0, 1), 1e-9);
    }

    @Test
    public void tayHandlesEmptyClasses() {
        Matrix L = new Matrix(2, 2);
        L.set(0, 0, 1.0);
        L.set(0, 1, 1.0);
        L.set(1, 0, 2.0);
        L.set(1, 1, 1.0);
        Pfqn_tay.Result r = Pfqn_tay.pfqn_tay(L, row(3, 0), row(0, 0));
        assertEquals(0.0, r.X.get(1), 1e-12);
        assertTrue(r.X.get(0) > 0);
        for (int i = 0; i < 2; i++) {
            assertEquals(0.0, r.Q.get(i, 1), 1e-12);
        }
    }

    @Test
    public void hstReproducesSuriFigure1() {
        // Suri (1983) Figure 1: M = 4, N = 7, workloads (1,1,1,2), station 4 analysed.
        Matrix L = col(1, 1, 1, 2);
        Pfqn_hst.Result r = Pfqn_hst.pfqn_hst(L, 7, 0.0, 3);
        assertEquals(0.48, r.X, 5e-3);
        assertEquals(0.96, r.U, 5e-3);
        double[] pubP = {0.037, 0.058, 0.087, 0.124, 0.166, 0.198, 0.198, 0.132};
        for (int i = 0; i < pubP.length; i++) {
            assertEquals(pubP[i], r.p[i], 1e-3);
        }
        double[] pubC = {0.023, 0.055, 0.097, 0.145, 0.186, 0.193, 0.132};
        for (int i = 0; i < pubC.length; i++) {
            assertEquals(pubC[i], Math.abs(r.c[i]), 1e-3);
        }
        assertEquals(0.831, r.total, 1e-3);
        assertEquals(0.102, r.worst, 1.5e-3);
        // Lemma 3.1: the total equals Q_i(N) - Q_i(N-1).
        double q7 = Pfqn_hst.pfqn_hst(L, 7, 0.0, 3).Q;
        double q6 = Pfqn_hst.pfqn_hst(L, 6, 0.0, 3).Q;
        assertEquals(q7 - q6, r.total, 1e-10);
        // Station 1 of the same system: published total 0.057.
        assertEquals(0.057, Pfqn_hst.pfqn_hst(L, 7, 0.0, 0).total, 1.5e-3);
    }
}
