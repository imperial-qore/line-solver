package jline.api;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import org.junit.jupiter.api.Test;

import jline.api.pfqn.nc.Pfqn_ca;
import jline.api.pfqn.nc.Pfqn_gerasimov;
import jline.io.Ret;
import jline.util.matrix.Matrix;

/**
 * Gerasimov's residue closed form for the normalizing constant, generalized to R classes.
 *
 * <p>A. I. Gerasimov, "On Normalizing Constants in Multiclass Queueing Networks",
 * Operations Research 43(4):704-711, 1995.
 *
 * <p>The paper's own worked example (Figure 1, Table I, and the two closed forms in the
 * Appendix) pins the R = 2 algorithm rather than this port. One superscript is lost in the
 * scan of the second Appendix form: its leading term prints as x21/(x21-x11)^3 where the
 * series it sums requires x21^(N1+4)/(x21-x11)^3, which is what is asserted here. With that
 * reading it agrees with pfqn_ca to 1e-14; with the printed reading it is off by a factor of
 * order one. Everything else is checked against pfqn_ca, which is exact and independent.
 * Values agree across MATLAB, JAR, Python native and C++.
 */
public class PfqnGerasimovTest {

    private static final double Y = 0.05202;    // Table I: x12 = x32
    private static final double X11 = 0.06627;  // Table I: x11

    private static Matrix mat(double[][] a) {
        Matrix m = new Matrix(a.length, a[0].length);
        for (int i = 0; i < a.length; i++) {
            for (int j = 0; j < a[0].length; j++) {
                m.set(i, j, a[i][j]);
            }
        }
        return m;
    }

    private static Matrix vec(double[] a) {
        Matrix m = new Matrix(1, a.length);
        for (int i = 0; i < a.length; i++) {
            m.set(0, i, a[i]);
        }
        return m;
    }

    @Test
    public void testAppendixFormIEqualLoads() {
        // Appendix (i), x11 = x21 = c: G = (y^2/6) c^N1 (N1^3+9N1^2+26N1+18).
        double c = 0.08;
        for (int n1 = 1; n1 <= 6; n1++) {
            Matrix L = mat(new double[][] { { c, Y }, { c, 0 }, { 0, Y } });
            Ret.pfqnNc g = Pfqn_gerasimov.pfqn_gerasimov(L, vec(new double[] { n1, 2 }));
            double gf = (Y * Y / 6) * Math.pow(c, n1) * (n1 * n1 * n1 + 9.0 * n1 * n1 + 26.0 * n1 + 18);
            assertEquals(gf, g.G, Math.abs(gf) * 1e-12);
        }
    }

    @Test
    public void testAppendixFormIIDistinctLoads() {
        // Appendix (ii), x11 != x21, with the leading superscript restored.
        double x21 = 5.0;
        for (int n1 = 1; n1 <= 6; n1++) {
            Matrix L = mat(new double[][] { { X11, Y }, { x21, 0 }, { 0, Y } });
            Ret.pfqnNc g = Pfqn_gerasimov.pfqn_gerasimov(L, vec(new double[] { n1, 2 }));
            double gf = (Y * Y / X11) * (Math.pow(x21, n1 + 4) / Math.pow(x21 - X11, 3)
                    + (n1 * n1 + 7.0 * n1 + 12) * Math.pow(X11, n1 + 2) / (2 * (X11 - x21))
                    - (n1 + 4) * Math.pow(X11, n1 + 3) / Math.pow(X11 - x21, 2)
                    + Math.pow(X11, n1 + 4) / Math.pow(x21 - X11, 3)
                    - Math.pow(x21, n1 + 1));
            assertEquals(gf, g.G, Math.abs(gf) * 1e-7);
        }
    }

    private void check(String label, double[][] L, double[] N, double[] Z) {
        Ret.pfqnNc ge = Pfqn_gerasimov.pfqn_gerasimov(mat(L), vec(N), vec(Z));
        Ret.pfqnNc ca = Pfqn_ca.pfqn_ca(mat(L), vec(N), vec(Z));
        assertEquals(ca.G, ge.G, Math.abs(ca.G) * 1e-12, label);
    }

    @Test
    public void testMatchesConvolution() {
        // The degeneracies the paper's hypotheses exclude are ordinary cases here:
        // coincident poles, tied class-2 demands, stations unvisited by a class and
        // identical station rows.
        check("fig1", new double[][] { { X11, Y }, { 5.0, 0 }, { 0, Y } },
                new double[] { 3, 2 }, new double[] { 0, 0 });
        check("coincident poles x11=x21", new double[][] { { 0.08, Y }, { 0.08, 0 }, { 0, Y } },
                new double[] { 4, 2 }, new double[] { 0, 0 });
        check("tied class-2 demands", new double[][] { { 1, 2 }, { 3, 2 }, { 5, 7 } },
                new double[] { 3, 2 }, new double[] { 0, 0 });
        check("identical station rows", new double[][] { { 1, 2 }, { 1, 2 }, { 5, 7 } },
                new double[] { 3, 3 }, new double[] { 0, 0 });
        check("class-2 unvisited stations", new double[][] { { 1, 0 }, { 3, 0 }, { 5, 7 } },
                new double[] { 4, 2 }, new double[] { 0, 0 });
        check("think time", new double[][] { { 1, 2 }, { 3, 4 } },
                new double[] { 3, 2 }, new double[] { 0.5, 1.5 });
        check("think time one class", new double[][] { { 1, 2 }, { 3, 4 }, { 2, 1 } },
                new double[] { 2, 3 }, new double[] { 1, 0 });
        check("single class", new double[][] { { 1 }, { 2 }, { 3 } },
                new double[] { 6 }, new double[] { 0 });
        check("single class with delay", new double[][] { { 1 }, { 2 }, { 3 } },
                new double[] { 6 }, new double[] { 2 });
        check("empty class", new double[][] { { 1, 2 }, { 3, 4 } },
                new double[] { 4, 0 }, new double[] { 0, 0 });
        check("three classes", new double[][] { { 1, 2, 3 }, { 3, 4, 1 }, { 2, 1, 2 } },
                new double[] { 2, 2, 2 }, new double[] { 0, 0, 0 });
        check("three classes with delay", new double[][] { { 1, 2, 3 }, { 3, 4, 1 }, { 2, 1, 2 } },
                new double[] { 3, 2, 1 }, new double[] { 0.7, 0, 0.3 });
        check("four classes", new double[][] { { 1, 2, 3, 1 }, { 3, 4, 1, 2 }, { 2, 1, 2, 3 } },
                new double[] { 2, 1, 2, 1 }, new double[] { 0, 0, 0, 0 });
        check("four classes with delay", new double[][] { { 1, 2, 3, 1 }, { 3, 4, 1, 2 }, { 2, 1, 2, 3 } },
                new double[] { 2, 1, 2, 1 }, new double[] { 0.3, 0.2, 0, 0.1 });
    }

    @Test
    public void testEliminatedPopulationIsFree() {
        // A population removed by residues enters only as a pole ORDER, so it costs
        // nothing: the answer must stay exact as it grows by three orders of magnitude.
        double[][] L = { { 1, 2 }, { 3, 1 }, { 2, 4 }, { 0.5, 0.7 } };
        // Capped at 2000 because the pfqn_ca REFERENCE is what costs here: it needs
        // 79 s at N2 = 20000 against 0.1 ms for the residue form, which is the point
        // of the test but not something to put in the suite's critical path.
        int[] pops = { 10, 100, 2000 };
        for (int k = 0; k < pops.length; k++) {
            Ret.pfqnNc ge = Pfqn_gerasimov.pfqn_gerasimov(mat(L), vec(new double[] { 6, pops[k] }));
            Ret.pfqnNc ca = Pfqn_ca.pfqn_ca(mat(L), vec(new double[] { 6, pops[k] }),
                    vec(new double[] { 0, 0 }));
            assertEquals(ca.lG, ge.lG, Math.abs(ca.lG) * 1e-12);
        }
    }

    @Test
    public void testMaxtermsRefusesRatherThanTruncates() {
        // A truncated residue sum is not a bound on G, it is a wrong number.
        double[][] L = { { 1.0, 2.0, 3.0 }, { 3.0, 1.0, 2.0 }, { 2.0, 3.0, 1.0 },
                         { 1.5, 2.5, 0.5 }, { 0.5, 1.5, 2.5 }, { 2.5, 0.5, 1.5 } };
        assertThrows(IllegalStateException.class, () -> Pfqn_gerasimov.pfqn_gerasimov(
                mat(L), vec(new double[] { 4, 8, 8 }), vec(new double[] { 0, 0, 0 }), 1e-12, 50));
    }

    @Test
    public void testRejectsNonintegerPopulation() {
        assertThrows(IllegalArgumentException.class, () -> Pfqn_gerasimov.pfqn_gerasimov(
                mat(new double[][] { { 1, 2 } }), vec(new double[] { 1.5, 2 })));
    }

    @Test
    public void testEmptyModelIsUnity() {
        Ret.pfqnNc g = Pfqn_gerasimov.pfqn_gerasimov(mat(new double[][] { { 1, 2 } }),
                vec(new double[] { 0, 0 }));
        assertEquals(1.0, g.G, 1e-14);
        assertTrue(Math.abs(g.lG) < 1e-14);
    }
}
