package jline.solvers.ba;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertTrue;

import org.junit.jupiter.api.Test;

import java.util.List;

import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.RoutingMatrix;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.processes.Cox2;
import jline.lang.processes.Exp;
import jline.solvers.NetworkAvgTable;
import jline.solvers.NetworkBoundsTable;
import jline.solvers.ba.SolverBA;
import jline.solvers.ctmc.SolverCTMC;
import jline.solvers.mva.SolverMVA;

/**
 * Tests for {@link SolverBA}: asymptotic/hierarchical throughput bounds bracket
 * the exact MVA solution and hierarchical families converge with level. Also
 * covers the Marie aggregation-decomposition method ('marie', a SolverMVA
 * method), which is exact for exponential service and close to CTMC for Coxian
 * FCFS service in both single- and multi-class settings.
 *
 * <p>Note: the geometric families (gb/sb) assume zero think time, so the
 * bracket test uses a no-delay cyclic model; Marie tests use think time.
 */
public class SolverBATest {

    private static final double TOL = 1e-9;

    // ---- model builders ----

    /** No-delay single-class closed model: Q0 -> Q1 -> Q2 (cyclic), FCFS exp. */
    private static Network cyclicExp(int N, double s0, double s1, double s2) {
        Network model = new Network("baCyc");
        Queue q0 = new Queue(model, "Q0", SchedStrategy.FCFS);
        Queue q1 = new Queue(model, "Q1", SchedStrategy.FCFS);
        Queue q2 = new Queue(model, "Q2", SchedStrategy.FCFS);
        ClosedClass c = new ClosedClass(model, "C", N, q0);
        q0.setService(c, new Exp(1.0 / s0));
        q1.setService(c, new Exp(1.0 / s1));
        q2.setService(c, new Exp(1.0 / s2));
        model.link(model.serialRouting(q0, q1, q2));
        return model;
    }

    /** No-delay two-class PS closed model (product-form, class-dependent rates). */
    private static Network cyclicPsTwoClass(int n1, int n2) {
        Network model = new Network("baPs2c");
        Queue q0 = new Queue(model, "Q0", SchedStrategy.PS);
        Queue q1 = new Queue(model, "Q1", SchedStrategy.PS);
        ClosedClass ca = new ClosedClass(model, "A", n1, q0);
        ClosedClass cb = new ClosedClass(model, "B", n2, q0);
        q0.setService(ca, new Exp(1.0 / 1.0));
        q0.setService(cb, new Exp(1.0 / 1.5));
        q1.setService(ca, new Exp(1.0 / 1.2));
        q1.setService(cb, new Exp(1.0 / 0.8));
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(ca, ca, model.serialRouting(q0, q1));
        P.set(cb, cb, model.serialRouting(q0, q1));
        model.link(P);
        return model;
    }

    /** Single-class closed model with Delay(Z) and Coxian FCFS service. */
    private static Network singleClassCox(int N, double Z, double s1, double scv1,
                                          double s2, double scv2) {
        Network model = new Network("baCox");
        Delay delay = new Delay(model, "Think");
        Queue q1 = new Queue(model, "Q1", SchedStrategy.FCFS);
        Queue q2 = new Queue(model, "Q2", SchedStrategy.FCFS);
        ClosedClass c = new ClosedClass(model, "C", N, delay);
        delay.setService(c, new Exp(1.0 / Z));
        q1.setService(c, Cox2.fitMeanAndSCV(s1, scv1));
        q2.setService(c, Cox2.fitMeanAndSCV(s2, scv2));
        model.link(model.serialRouting(delay, q1, q2));
        return model;
    }

    /** Single-class closed model with Delay(Z) and exponential FCFS service. */
    private static Network singleClassExp(int N, double Z, double s1, double s2) {
        Network model = new Network("baExp");
        Delay delay = new Delay(model, "Think");
        Queue q1 = new Queue(model, "Q1", SchedStrategy.FCFS);
        Queue q2 = new Queue(model, "Q2", SchedStrategy.FCFS);
        ClosedClass c = new ClosedClass(model, "C", N, delay);
        delay.setService(c, new Exp(1.0 / Z));
        q1.setService(c, new Exp(1.0 / s1));
        q2.setService(c, new Exp(1.0 / s2));
        model.link(model.serialRouting(delay, q1, q2));
        return model;
    }

    /** Two-class closed model, Delay -> Q1(FCFS) -> Q2(FCFS), Coxian per class. */
    private static Network twoClassCox(int n1, int n2) {
        Network model = new Network("ba2cCox");
        Delay delay = new Delay(model, "Think");
        Queue q1 = new Queue(model, "Q1", SchedStrategy.FCFS);
        Queue q2 = new Queue(model, "Q2", SchedStrategy.FCFS);
        ClosedClass ca = new ClosedClass(model, "A", n1, delay);
        ClosedClass cb = new ClosedClass(model, "B", n2, delay);
        delay.setService(ca, new Exp(1.0));
        delay.setService(cb, new Exp(1.0));
        q1.setService(ca, Cox2.fitMeanAndSCV(1.0, 0.5));
        q1.setService(cb, Cox2.fitMeanAndSCV(1.4, 0.5));
        q2.setService(ca, Cox2.fitMeanAndSCV(1.2, 0.5));
        q2.setService(cb, Cox2.fitMeanAndSCV(0.8, 0.5));
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(ca, ca, model.serialRouting(delay, q1, q2));
        P.set(cb, cb, model.serialRouting(delay, q1, q2));
        model.link(P);
        return model;
    }

    /** Two-class closed model with class-independent FCFS rates (product-form). */
    private static Network twoClassExpPF(int n1, int n2, double s1, double s2) {
        Network model = new Network("ba2cPF");
        Delay delay = new Delay(model, "Think");
        Queue q1 = new Queue(model, "Q1", SchedStrategy.FCFS);
        Queue q2 = new Queue(model, "Q2", SchedStrategy.FCFS);
        ClosedClass ca = new ClosedClass(model, "A", n1, delay);
        ClosedClass cb = new ClosedClass(model, "B", n2, delay);
        delay.setService(ca, new Exp(1.0));
        delay.setService(cb, new Exp(1.0));
        q1.setService(ca, new Exp(1.0 / s1));
        q1.setService(cb, new Exp(1.0 / s1));
        q2.setService(ca, new Exp(1.0 / s2));
        q2.setService(cb, new Exp(1.0 / s2));
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(ca, ca, model.serialRouting(delay, q1, q2));
        P.set(cb, cb, model.serialRouting(delay, q1, q2));
        model.link(P);
        return model;
    }

    private static double refTput(NetworkAvgTable t) {
        // Reference station is row 0 in all builders here.
        return t.getTput().get(0);
    }

    // ---- bound bracketing ----

    @Test
    public void testAllBoundFamiliesBracketExact() {
        Network model = cyclicExp(6, 1.0, 1.5, 0.7);
        double exact = refTput(new SolverMVA(model, "exact").getAvgTable());

        String[] fams = {"aba", "bjb", "pb", "gb", "sb", "mwba"};
        for (String fam : fams) {
            SolverBA lo = new SolverBA(model);
            lo.getOptions().method(fam + ".lower");
            double xl = refTput(lo.getAvgTable());
            SolverBA up = new SolverBA(model);
            up.getOptions().method(fam + ".upper");
            double xu = refTput(up.getAvgTable());
            assertTrue(xl <= exact + 1e-6,
                    fam + " lower (" + xl + ") should not exceed exact (" + exact + ")");
            assertTrue(exact <= xu + 1e-6,
                    fam + " upper (" + xu + ") should not fall below exact (" + exact + ")");
            assertTrue(xl <= xu + 1e-9, fam + " lower must not exceed upper");
        }
    }

    @Test
    public void testHierarchicalConvergeAtLevelN() {
        int N = 5;
        Network model = cyclicExp(N, 1.0, 2.0, 0.5);
        double exact = refTput(new SolverMVA(model, "exact").getAvgTable());

        for (String fam : new String[]{"pbh", "cbh"}) {
            SolverBA lo = new SolverBA(model);
            lo.getOptions().method(fam + ".lower");
            lo.getOptions().level = N;
            double xl = refTput(lo.getAvgTable());
            SolverBA up = new SolverBA(model);
            up.getOptions().method(fam + ".upper");
            up.getOptions().level = N;
            double xu = refTput(up.getAvgTable());
            assertEquals(exact, xl, 1e-6, fam + ".lower at level N should equal exact");
            assertEquals(exact, xu, 1e-6, fam + ".upper at level N should equal exact");
        }
    }

    @Test
    public void testCubMbjbBracketMulticlass() {
        Network model = cyclicPsTwoClass(3, 3);
        double x1 = new SolverMVA(model, "exact").getAvgTable().getTput().get(0);
        SolverBA up = new SolverBA(model);
        up.getOptions().method("cub.upper");
        double xu = up.getAvgTable().getTput().get(0);
        SolverBA lo = new SolverBA(model);
        lo.getOptions().method("mbjb.lower");
        double xl = lo.getAvgTable().getTput().get(0);
        assertTrue(xl <= x1 + 1e-6, "mbjb lower (" + xl + ") should bracket exact (" + x1 + ") from below");
        assertTrue(x1 <= xu + 1e-6, "cub upper (" + xu + ") should bracket exact (" + x1 + ") from above");
    }

    // ---- Marie (SolverMVA method): single-class ----

    @Test
    public void testMarieExponentialIsExact() {
        Network model = singleClassExp(6, 2.0, 1.0, 1.5);
        NetworkAvgTable exact = new SolverMVA(model, "exact").getAvgTable();
        NetworkAvgTable m = new SolverMVA(model, "marie").getAvgTable();
        List<Double> qe = exact.getQLen();
        List<Double> qm = m.getQLen();
        for (int i = 0; i < qe.size(); i++) {
            assertEquals(qe.get(i), qm.get(i), TOL,
                    "Marie exponential queue length must equal exact MVA at station " + i);
        }
    }

    @Test
    public void testMarieCoxianCloseToCtmc() {
        Network model = singleClassCox(4, 1.0, 1.0, 0.5, 1.2, 0.5);
        NetworkAvgTable ctmc = new SolverCTMC(model, "keep", false).getAvgTable();
        NetworkAvgTable m = new SolverMVA(model, "marie").getAvgTable();
        List<Double> qc = ctmc.getQLen();
        List<Double> qm = m.getQLen();
        for (int i = 0; i < qc.size(); i++) {
            assertEquals(qc.get(i), qm.get(i), 0.15 + 0.1 * qc.get(i),
                    "Marie Coxian queue length near CTMC at station " + i);
        }
    }

    // ---- Marie (SolverMVA method): multiclass ----

    @Test
    public void testMarieMulticlassExponentialIsExact() {
        // scv==1 and class-independent -> exact MVA dispatch inside Marie.
        Network model = twoClassExpPF(2, 3, 1.0, 1.5);
        NetworkAvgTable exact = new SolverMVA(model, "exact").getAvgTable();
        NetworkAvgTable m = new SolverMVA(model, "marie").getAvgTable();
        List<Double> qe = exact.getQLen();
        List<Double> qm = m.getQLen();
        for (int i = 0; i < qe.size(); i++) {
            assertEquals(qe.get(i), qm.get(i), 1e-6,
                    "Multiclass Marie exponential must equal exact MVA at row " + i);
        }
    }

    @Test
    public void testMarieMulticlassCoxianCloseToCtmc() {
        Network model = twoClassCox(2, 2);
        NetworkAvgTable ctmc = new SolverCTMC(model, "keep", false).getAvgTable();
        NetworkAvgTable m = new SolverMVA(model, "marie").getAvgTable();
        assertNotNull(m);
        List<Double> qc = ctmc.getQLen();
        List<Double> qm = m.getQLen();
        for (int i = 0; i < qc.size(); i++) {
            assertEquals(qc.get(i), qm.get(i), 0.25 + 0.15 * qc.get(i),
                    "Multiclass Marie Coxian near CTMC at row " + i);
        }
    }

    // ---- lr: LP linear-reduction bound ----

    /**
     * The lr LP bound must bracket the exact solution per station on
     * asymmetric models. Asymmetry matters: a symmetric model happens to mask
     * a missing non-negativity restriction on the LP variables, which is the
     * defect that previously made the JAR diverge from MATLAB/Python.
     */
    @Test
    public void testLrLpBracketsExact() {
        Network[] models = {
                cyclicExp(3, 1.0, 0.5, 2.0),          // rates [1, 2, 0.5]
                cyclicExp(3, 1.0, 1.0, 1.0),          // symmetric
                cyclicExp(4, 0.5, 1.0, 1.0 / 3.0),    // rates [2, 1, 3]
                cyclicExp(5, 1.0, 1.0 / 3.0, 1.0 / 0.7) // rates [1, 3, 0.7]
        };
        for (Network model : models) {
            List<Double> exact = new SolverMVA(model, "exact").getAvgTable().getUtil();
            List<Double> lower = new SolverBA(model, "lr.lower").getAvgTable().getUtil();
            List<Double> upper = new SolverBA(model, "lr.upper").getAvgTable().getUtil();
            for (int i = 0; i < exact.size(); i++) {
                assertTrue(lower.get(i) <= exact.get(i) + 1e-6,
                        "lr.lower (" + lower.get(i) + ") must not exceed exact ("
                                + exact.get(i) + ") at station " + i);
                assertTrue(exact.get(i) <= upper.get(i) + 1e-6,
                        "lr.upper (" + upper.get(i) + ") must not fall below exact ("
                                + exact.get(i) + ") at station " + i);
            }
        }
    }

    /**
     * Cross-codebase parity: the asym3 model (3 FCFS queues, rates [1, 2, 0.5],
     * cyclic, N=3) must reproduce the utilizations produced by MATLAB and native
     * Python, which agree to 1e-6.
     */
    @Test
    public void testLrLpMatchesMatlabPythonReference() {
        Network model = cyclicExp(3, 1.0, 0.5, 2.0);
        double[] refLower = {0.436063, 0.218032, 0.872125};
        double[] refUpper = {0.495412, 0.247706, 0.990826};

        List<Double> lower = new SolverBA(model, "lr.lower").getAvgTable().getUtil();
        List<Double> upper = new SolverBA(model, "lr.upper").getAvgTable().getUtil();
        for (int i = 0; i < refLower.length; i++) {
            assertEquals(refLower[i], lower.get(i), 1e-6,
                    "lr.lower parity with MATLAB/Python at station " + i);
            assertEquals(refUpper[i], upper.get(i), 1e-6,
                    "lr.upper parity with MATLAB/Python at station " + i);
        }
    }

    // ---- getBoundsTable ----

    /**
     * The table must carry exactly the numbers getBounds returns, for a
     * two-sided family, and its rows must be ordered lower <= upper.
     */
    @Test
    public void testBoundsTableAgreesWithGetBounds() {
        Network model = cyclicExp(6, 1.0, 1.5, 0.7);
        SolverBA s = new SolverBA(model, "gb.upper");
        SolverBA.Bounds b = s.getBounds();
        NetworkBoundsTable t = new SolverBA(model, "gb.upper").getBoundsTable();

        assertNotNull(t);
        assertEquals(3, t.getQlower().size(), "one row per station-class pair");
        for (int i = 0; i < t.getQlower().size(); i++) {
            assertEquals(b.Qlower.get(i, 0), t.getQlower().get(i), 1e-12,
                    "Qlower row " + i + " must match getBounds");
            assertEquals(b.Qupper.get(i, 0), t.getQupper().get(i), 1e-12,
                    "Qupper row " + i + " must match getBounds");
            assertEquals(b.Tlower.get(i, 0), t.getTlower().get(i), 1e-12,
                    "Tlower row " + i + " must match getBounds");
            assertEquals(b.Tupper.get(i, 0), t.getTupper().get(i), 1e-12,
                    "Tupper row " + i + " must match getBounds");
            assertTrue(t.getQlower().get(i) <= t.getQupper().get(i) + 1e-9,
                    "Qlower must not exceed Qupper at row " + i);
            assertTrue(t.getTlower().get(i) <= t.getTupper().get(i) + 1e-9,
                    "Tlower must not exceed Tupper at row " + i);
        }
    }

    /**
     * A one-sided family must yield NaN on the missing side without erroring,
     * and must not have the row dropped or the value coerced to zero. cub is
     * upper-only, so the lower side is absent.
     */
    @Test
    public void testBoundsTableOneSidedFamilyKeepsNaN() {
        Network model = cyclicPsTwoClass(3, 3);
        NetworkBoundsTable t = new SolverBA(model, "cub.upper").getBoundsTable();

        assertNotNull(t);
        assertTrue(t.getQupper().size() > 0, "one-sided family must still yield rows");
        for (int i = 0; i < t.getQupper().size(); i++) {
            assertTrue(Double.isNaN(t.getQlower().get(i)),
                    "cub has no lower side: Qlower must be NaN at row " + i);
            assertTrue(Double.isNaN(t.getTlower().get(i)),
                    "cub has no lower side: Tlower must be NaN at row " + i);
            assertTrue(!Double.isNaN(t.getQupper().get(i)),
                    "cub upper side must be present at row " + i);
        }
    }

    /**
     * Regression for the option-propagation defect: getBounds re-runs the solver
     * for each side, and if it does not inherit the caller's options the level
     * silently reverts to the default of 2. A hierarchical family would then
     * never tighten however high the user set level. At level = N the bracket
     * must collapse onto the exact solution.
     */
    @Test
    public void testGetBoundsPropagatesLevel() {
        int N = 5;
        Network model = cyclicExp(N, 1.0, 2.0, 0.5);
        double exact = refTput(new SolverMVA(model, "exact").getAvgTable());

        SolverBA coarse = new SolverBA(model, "pbh.upper");
        coarse.getOptions().level = 2;
        SolverBA.Bounds b2 = coarse.getBounds();
        double width2 = b2.Tupper.get(0, 0) - b2.Tlower.get(0, 0);

        SolverBA fine = new SolverBA(model, "pbh.upper");
        fine.getOptions().level = N;
        SolverBA.Bounds bN = fine.getBounds();
        double widthN = bN.Tupper.get(0, 0) - bN.Tlower.get(0, 0);

        assertTrue(width2 > 1e-6, "level 2 bracket should be strictly loose");
        assertTrue(widthN < width2 - 1e-9,
                "level N bracket (" + widthN + ") must be tighter than level 2 (" + width2 + ")");
        assertEquals(exact, bN.Tlower.get(0, 0), 1e-6,
                "pbh.lower at level N reached via getBounds should equal exact");
        assertEquals(exact, bN.Tupper.get(0, 0), 1e-6,
                "pbh.upper at level N reached via getBounds should equal exact");
    }

    /** getBoundsTable inherits the level fix, since it routes through getBounds. */
    @Test
    public void testBoundsTablePropagatesLevel() {
        int N = 5;
        Network model = cyclicExp(N, 1.0, 2.0, 0.5);
        double exact = refTput(new SolverMVA(model, "exact").getAvgTable());

        SolverBA fine = new SolverBA(model, "pbh.upper");
        fine.getOptions().level = N;
        NetworkBoundsTable t = fine.getBoundsTable();

        assertEquals(exact, t.getTlower().get(0), 1e-6,
                "table Tlower at level N should equal exact");
        assertEquals(exact, t.getTupper().get(0), 1e-6,
                "table Tupper at level N should equal exact");
    }
}
