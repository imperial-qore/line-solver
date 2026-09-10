package jline.solvers.ba;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import org.junit.jupiter.api.Test;

import java.util.List;

import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.RoutingMatrix;
import jline.lang.constant.DropStrategy;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.processes.Cox2;
import jline.lang.processes.Exp;
import jline.lang.processes.MMPP2;
import jline.lang.processes.Erlang;
import jline.solvers.NetworkAvgTable;
import jline.solvers.QrfParams;
import jline.solvers.SolverOptions;
import jline.api.pfqn.mva.Pfqn_scb;
import jline.util.matrix.Matrix;
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

    /**
     * The Harel-Namn-Sturm sharp bounds. They are NOT the 'sb' family of the
     * same paper: the upper side extrapolates from the EXACT normalizing
     * constant at populations n &lt;= 7, so it must be at least as tight as sb
     * and exact when the extrapolation point is N itself.
     */
    @Test
    public void testHarelBracketsExactAndIsTighterThanSb() {
        for (int N : new int[]{3, 5, 10, 20}) {
            Network model = cyclicExp(N, 0.5, 1.0 / 3.0, 0.2);
            double exact = refTput(new SolverMVA(model, "exact").getAvgTable());
            SolverBA lo = new SolverBA(model);
            lo.getOptions().method("harel.lower");
            double xl = refTput(lo.getAvgTable());
            SolverBA up = new SolverBA(model);
            up.getOptions().method("harel.upper");
            double xu = refTput(up.getAvgTable());
            assertTrue(xl <= exact + 1e-9, "harel.lower (" + xl + ") exceeds exact at N=" + N);
            assertTrue(exact <= xu + 1e-9, "harel.upper (" + xu + ") below exact at N=" + N);
            SolverBA sbUp = new SolverBA(model);
            sbUp.getOptions().method("sb.upper");
            assertTrue(xu <= refTput(sbUp.getAvgTable()) + 1e-9,
                    "harel.upper must not be looser than sb.upper at N=" + N);
            if (N <= 7) {
                assertEquals(exact, xu, 1e-9,
                        "at N <= 7 the extrapolation point is N itself, so the bound is exact");
            }
        }
    }

    /** Harel forbids think times rather than dropping them. */
    @Test
    public void testHarelRejectsDelayStation() {
        Network model = singleClassExp(5, 2.0, 1.0, 1.5);
        SolverBA up = new SolverBA(model);
        up.getOptions().method("harel.upper");
        assertThrows(RuntimeException.class, up::getAvgTable);
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

    /**
     * The QRF LP tokens DERIVE their blocking parameters, and agree with the
     * hand-built ones.
     *
     * <p>They used to refuse a model that carried no {@code qrfParams}, because
     * the only alternative then was to INVENT the tables and assuming no
     * blocking answers a different model -- roughly 31x farther from exact.
     * {@link jline.api.sn.SnToQrfBlocking} derives them from the model instead,
     * and on this UNBLOCKED cyclic network the derivation is exactly the
     * trivial table {@link #noBlockingParams} spells out: one configuration, no
     * blocked station, capacity N. So the two paths must agree entry for entry,
     * which is what pins the derivation rather than merely exercising it. The
     * bound itself is checked against the exact CTMC in the direction the
     * relaxation guarantees, utilization from above.
     */
    @Test
    public void testQrfLpTokensDeriveBlockingParameters() {
        int N = 3;
        Network model = cyclicExp(N, 1.0, 1.5, 2.0);
        for (String token : new String[]{"qrf.bas", "qrf.rsrd"}) {
            SolverOptions bare = SolverBA.defaultOptions();
            bare.verbose = jline.VerboseLevel.SILENT;
            bare.method = token;
            Matrix derived = new SolverBA(model, bare).getAvgUtil();

            SolverOptions opt = SolverBA.defaultOptions();
            opt.verbose = jline.VerboseLevel.SILENT;
            opt.method = token;
            opt.qrfParams = noBlockingParams(3, N);
            Matrix U = new SolverBA(model, opt).getAvgUtil();
            for (int i = 0; i < 3; i++) {
                assertEquals(U.get(i, 0), derived.get(i, 0), 1e-9,
                        token + " derived tables must match the hand-built no-blocking ones"
                                + " at station " + i);
            }
            Matrix exU = new SolverCTMC(model).getAvgUtil();
            for (int i = 0; i < 3; i++) {
                assertTrue(U.get(i, 0) >= -1e-9 && U.get(i, 0) <= 1.0 + 1e-9,
                        token + " utilization out of [0,1] at station " + i);
                assertTrue(U.get(i, 0) >= exU.get(i, 0) - 1e-6,
                        token + " is not an upper bound at station " + i
                                + ": " + U.get(i, 0) + " < " + exU.get(i, 0));
            }
        }
    }

    /**
     * The ALPHA-FREE arms reject what their single-server formulation cannot model.
     *
     * The refusal NAMES THE TWO ARMS THAT DO SERVE THE MODEL rather than the
     * construct that rules the others out: since the load-dependent arms landed,
     * alpha(i,n) is the rate law of a delay (alpha = n), of a c-server station
     * (alpha = min(n,c)) and of limited load dependence alike, so one message
     * covers all three and points at 'qrf.mmi.ld' / 'qrf.mmi.linear'. Matching on
     * that pointer is what the python twins assert
     * (`test_qrf_alpha_free_arms_reject_delay_stations`,
     * `test_infinite_server_is_rejected_by_the_alpha_free_arms`), and it is the
     * stable half of the text: the words "infinite-server" and "multi-server" are
     * gone from the message precisely because it no longer refuses per construct.
     */
    @Test
    public void testQrfRejectsDelayAndMultiserver() {
        Network delayModel = singleClassExp(3, 1.0, 1.0, 2.0);
        SolverOptions o1 = SolverBA.defaultOptions();
        o1.verbose = jline.VerboseLevel.SILENT;
        o1.method = "qrf.mmi";
        final SolverBA s1 = new SolverBA(delayModel, o1);
        RuntimeException e1 = assertThrows(RuntimeException.class,
                new org.junit.jupiter.api.function.Executable() {
                    public void execute() {
                        s1.getAvgTable();
                    }
                }, "qrf.mmi must reject an infinite-server station");
        assertTrue(e1.getMessage().contains("Use 'qrf.mmi.ld' or 'qrf.mmi.linear'"),
                "expected an infinite-server rejection naming the ld arms, got: "
                        + e1.getMessage());

        Network msModel = cyclicExp(3, 1.0, 1.5, 2.0);
        msModel.getStations().get(1).setNumberOfServers(3);
        SolverOptions o2 = SolverBA.defaultOptions();
        o2.verbose = jline.VerboseLevel.SILENT;
        o2.method = "qrf.mmi";
        final SolverBA s2 = new SolverBA(msModel, o2);
        RuntimeException e2 = assertThrows(RuntimeException.class,
                new org.junit.jupiter.api.function.Executable() {
                    public void execute() {
                        s2.getAvgTable();
                    }
                }, "qrf.mmi must reject a multi-server station");
        assertTrue(e2.getMessage().contains("Use 'qrf.mmi.ld' or 'qrf.mmi.linear'"),
                "expected a multiserver rejection naming the ld arms, got: "
                        + e2.getMessage());
    }

    /**
     * The NLP tokens optimize: they must move off the phase-1 feasible start,
     * and the two objectives must not agree by construction.
     *
     * <p>This is the regression that motivated the reduced-space rewrite. A
     * backend that silently returns the start point reports the SAME numbers
     * for every objective, which is how the defect stayed invisible.
     */
    @Test
    public void testQrfNlpObjectivesDiffer() {
        Network model = cyclicExp(3, 1.0, 1.5, 2.0);
        Matrix mmi = qrfUtil(model, "qrf.mmi");
        Matrix mem = qrfUtil(model, "qrf.mem");
        Matrix ld = qrfUtil(model, "qrf.mmi.ld");
        for (int i = 0; i < 3; i++) {
            assertTrue(mmi.get(i, 0) > 0.0 && mmi.get(i, 0) <= 1.0 + 1e-9,
                    "qrf.mmi utilization out of (0,1] at station " + i);
            assertTrue(mem.get(i, 0) > 0.0 && mem.get(i, 0) <= 1.0 + 1e-9,
                    "qrf.mem utilization out of (0,1] at station " + i);
            // alpha = 1 makes the load-dependent variant the same problem
            assertEquals(mmi.get(i, 0), ld.get(i, 0), 1e-6,
                    "qrf.mmi.ld must equal qrf.mmi at unit load dependence");
        }
        double gap = 0.0;
        for (int i = 0; i < 3; i++) {
            gap = Math.max(gap, Math.abs(mmi.get(i, 0) - mem.get(i, 0)));
        }
        assertTrue(gap > 1e-6,
                "MMI and MEM returned the same point (" + gap + "), which is what "
                        + "a backend that never leaves the feasible start does");
    }

    // ---- scb: Dowdy et al. (1992) single-class bounds ----

    /**
     * scb is the only SolverBA family that does NOT bracket the given model.
     * Its lower side IS the model's exact single-class throughput; what the
     * pair brackets is the multiclass system the single-class model
     * aggregates. Asserting the usual "brackets the exact answer" property
     * would therefore be asserting the wrong thing.
     */
    @Test
    public void testScbBracketsTheMulticlassSystemNotThisModel() {
        Network model = cyclicExp(3, 0.5, 1.0 / 3.0, 0.2);
        double exact = refTput(new SolverMVA(model, "exact").getAvgTable());
        SolverBA lo = new SolverBA(model);
        lo.getOptions().method("scb.lower");
        double xl = refTput(lo.getAvgTable());
        SolverBA up = new SolverBA(model);
        up.getOptions().method("scb.upper");
        double xu = refTput(up.getAvgTable());
        assertEquals(exact, xl, 1e-9, "scb.lower must be the EXACT single-class throughput");
        assertTrue(xu >= xl, "scb.upper below scb.lower");
        // Theorem 3 at N = 3 customers over K = 3 devices: the gap is 2/5.
        double gap = Pfqn_scb.pfqn_scbgap(3, 3);
        assertTrue(xu / xl <= 1.0 + gap / (1.0 - gap) + 1e-9,
                "scb.upper exceeds the Theorem-3 aggregation gap");
        // Corollary 1: one utilization ratio for every device, all capped at 1.
        Matrix Ulo = lo.getAvgUtil();
        Matrix Uhi = up.getAvgUtil();
        for (int i = 0; i < Ulo.getNumRows(); i++) {
            assertEquals(xu / xl, Uhi.get(i, 0) / Ulo.get(i, 0), 1e-9,
                    "utilization ratio is not uniform at station " + i);
            assertTrue(Uhi.get(i, 0) <= 1.0 + 1e-9, "U > 1 at station " + i);
        }
    }

    /** A delay station and a multiclass model are refused by name, not served. */
    @Test
    public void testScbRejectsDelayAndMulticlass() {
        SolverBA withDelay = new SolverBA(singleClassExp(5, 2.0, 1.0, 1.5));
        withDelay.getOptions().method("scb.lower");
        assertThrows(RuntimeException.class, withDelay::getAvgTable);
        SolverBA multiclass = new SolverBA(cyclicPsTwoClass(2, 2));
        multiclass.getOptions().method("scb.upper");
        assertThrows(RuntimeException.class, multiclass::getAvgTable);
    }

    /**
     * auto.lower must stay a bound on THIS model, so scb cannot be one of its
     * candidates: scb.lower equals the exact throughput, and admitting it would
     * collapse the composite lower bound onto the exact answer.
     */
    @Test
    public void testScbIsNotAnAutoCandidate() {
        Network model = cyclicExp(3, 0.5, 1.0 / 3.0, 0.2);
        double exact = refTput(new SolverMVA(model, "exact").getAvgTable());
        SolverBA auto = new SolverBA(model);
        auto.getOptions().method("auto.lower");
        assertTrue(refTput(auto.getAvgTable()) < exact,
                "auto.lower reached the exact answer, so scb leaked into the composite");
        List<String> valid = java.util.Arrays.asList(new SolverBA(model).listValidMethods());
        assertTrue(valid.contains("scb.upper") && valid.contains("scb.lower"));
    }

    /** The published tables of Dowdy et al. (1992), which the API must match. */
    @Test
    public void testScbCombinatorialBoundsMatchThePublishedTables() {
        // Table I, p.200: rows N = 1..5, columns K = 1..5, in percent.
        int[][] tableI = {{0, 0, 0, 0, 0}, {0, 33, 33, 33, 33}, {0, 25, 40, 40, 40},
                          {0, 20, 33, 43, 43}, {0, 17, 29, 38, 44}};
        for (int N = 1; N <= 5; N++) {
            for (int K = 1; K <= 5; K++) {
                assertEquals(tableI[N - 1][K - 1],
                        (int) Math.round(100 * Pfqn_scb.pfqn_scbgap(N, K)),
                        "Table I entry N=" + N + " K=" + K);
            }
        }
        assertTrue(Pfqn_scb.pfqn_scbgap(1000, 1000) < 0.5, "Theorem 3 caps the error at 50%");
        // The undominated form is defined for r <= K and refused past it.
        for (int r = 2; r < 5; r++) {
            assertTrue(Pfqn_scb.pfqn_scbgap(8, 5, r, true) < Pfqn_scb.pfqn_scbgap(8, 5, r, false),
                    "the undominated form must be tighter at r=" + r);
        }
        assertThrows(RuntimeException.class, () -> Pfqn_scb.pfqn_scbgap(8, 5, 6, true));
        // Section 4.7: K = 2 devices, N = 3 customers, one class admits 1.5.
        assertEquals(1.5, Pfqn_scb.pfqn_usumbound(1, 2, 3), 1e-12);
        assertEquals(2.0, Pfqn_scb.pfqn_usumbound(2, 2, 3), 1e-12);
        assertEquals(2, Pfqn_scb.pfqn_minclasses(1.6, 2, 3));
        assertEquals(1, Pfqn_scb.pfqn_minclasses(1.4, 2, 3));
        assertEquals(-1, Pfqn_scb.pfqn_minclasses(2.5, 2, 3), "unattainable by any class structure");
        // Section 2 example: aggregate at 8.152, multiclass counterpart at 8.761.
        Matrix L = new Matrix(3, 1);
        L.set(0, 0, 0.114);
        L.set(1, 0, 0.040);
        L.set(2, 0, 0.062);
        double[] b = Pfqn_scb.pfqn_scb(L, 4);
        assertEquals(8.151894490678735, b[0], 1e-9);
        assertTrue(b[0] <= 8.7615 && b[1] >= 8.7615, "the paper's multiclass X is outside the bracket");
        assertEquals(1.0, b[2 + L.length()], 1e-12, "the single-server cap must bind on the busiest device");
    }

    /**
     * The `cqn_bas_blocking` example: Queue1 declares the BAS drop rule and
     * Queue2 has a buffer of 1, below the population of 2, so the buffer BINDS.
     */
    private static Network basBlocking() {
        Network model = new Network("cqn_bas_blocking");
        Queue q1 = new Queue(model, "Queue1", SchedStrategy.FCFS);
        Queue q2 = new Queue(model, "Queue2", SchedStrategy.FCFS);
        ClosedClass c = new ClosedClass(model, "Class1", 2, q1, 0);
        q1.setService(c, new Exp(1.0));
        q2.setService(c, new Exp(0.8));
        q2.setCap(1);
        q1.setDropRule(c, DropStrategy.BlockingAfterService);
        model.link(model.serialRouting(q1, q2));
        return model;
    }

    /**
     * A blocking-blind family must be REFUSED on a model with a binding finite
     * buffer, not silently answered with the unblocked model's bound. Before
     * the gate, `gb.upper` on cqn_bas_blocking returned QLen 1.28 at Queue2 --
     * a station that can never hold more than one job -- because the cap and
     * the BAS drop rule were never read.
     */
    @Test
    public void testBlockingIsRefusedNotBounded() {
        Network model = basBlocking();
        SolverOptions opt = SolverBA.defaultOptions();
        opt.verbose = jline.VerboseLevel.SILENT;

        // A blind family asked for BY NAME still refuses.
        SolverOptions og = SolverBA.defaultOptions();
        og.verbose = jline.VerboseLevel.SILENT;
        og.method = "gb.upper";
        assertThrows(RuntimeException.class, () -> new SolverBA(model, og).getAvgTable());
        assertThrows(RuntimeException.class, () -> new SolverBA(model, og).getBounds(),
                "getBounds must refuse rather than return an all-null bracket");
        // 'auto.lower' has no blocking counterpart to route to: the analyzer
        // solves qrf.bas in the 'max' direction alone.
        SolverOptions oal = SolverBA.defaultOptions();
        oal.verbose = jline.VerboseLevel.SILENT;
        oal.method = "auto.lower";
        assertThrows(RuntimeException.class, () -> new SolverBA(model, oal).getAvgTable());

        // But the UPPER-side aliases MEAN the QRF BAS bound here. They used to
        // be refused with the blind families, which left this model with no
        // usable SolverBA method at all: the advice was 'qrf.bas', and that in
        // turn demanded hand-built blocking tables. SnToQrfBlocking derives
        // them, so the alias resolves -- the routing SolverMVA performs for
        // 'sqd'.
        SolverOptions oq = SolverBA.defaultOptions();
        oq.verbose = jline.VerboseLevel.SILENT;
        oq.method = "qrf.bas";
        Matrix named = new SolverBA(model, oq).getAvgUtil();
        // Bare "auto" included on purpose: checkDeclaredMethod validates the RAW
        // name, so it only reaches resolveMethod because listAllMethods declares
        // it beside "qr" and "lr". It used to be missing from that list, which
        // made this very assertion fail with "The 'auto' method is unsupported
        // by this solver" -- so this case is the regression test for that entry.
        for (String alias : new String[]{"default", "auto", "auto.upper"}) {
            SolverOptions oa = SolverBA.defaultOptions();
            oa.verbose = jline.VerboseLevel.SILENT;
            oa.method = alias;
            Matrix routed = new SolverBA(model, oa).getAvgUtil();
            for (int i = 0; i < 2; i++) {
                assertEquals(named.get(i, 0), routed.get(i, 0), 1e-9,
                        "'" + alias + "' must resolve to qrf.bas on a blocked model");
            }
        }
        // getBounds refuses the routed default as ONE-SIDED, not as blind:
        // saying "does not support blocking" would contradict the run above.
        assertThrows(RuntimeException.class, () -> new SolverBA(model, opt).getBounds());

        // The narrowed list keeps only the families that model the blocking,
        // plus 'default' -- the one entry whose RESOLVED form ('gb.upper') is
        // blind, which is exactly why runAnalyzer rewrites it.
        List<String> valid = java.util.Arrays.asList(new SolverBA(model, opt).listValidMethods());
        assertTrue(!valid.isEmpty(), "the QRF blocking bounds stay available on this model");
        assertTrue(valid.contains("default"), "the model's own default must be offered back");
        for (int i = 0; i < valid.size(); i++) {
            if (valid.get(i).equals("default")) {
                continue;
            }
            assertTrue(valid.get(i).startsWith("qrf.bas") || valid.get(i).startsWith("qrf.rsrd"),
                    "blocking-blind method still offered: " + valid.get(i));
        }
        // An unblocked model of the same shape is untouched by the gate.
        Network open = cyclicExp(2, 1.0, 1.25, 1.0);
        assertNotNull(new SolverBA(open, opt).getAvgTable());
    }

    private static Matrix qrfUtil(Network model, String method) {
        SolverOptions opt = SolverBA.defaultOptions();
        opt.verbose = jline.VerboseLevel.SILENT;
        opt.method = method;
        return new SolverBA(model, opt).getAvgUtil();
    }

    /** No blocking: one configuration, no blocked station, capacity N. */
    private static QrfParams noBlockingParams(int M, int N) {
        QrfParams qp = new QrfParams();
        qp.f = 1;
        qp.MR = 1;
        qp.ZM = 0;
        qp.MM = new Matrix(1, 2);
        qp.MM1 = new Matrix(1, M);
        qp.BB = new Matrix(1, M);
        qp.ZZ = new int[]{0};
        qp.F = new int[M];
        for (int i = 0; i < M; i++) qp.F[i] = N;
        return qp;
    }

    @Test
    public void testFeatureSetAcceptsTheModelsTheBoundsAreDerivedFor() {
        // SolverBA declared NO feature set, so it inherited Solver.supports,
        // which returns true: every model in the language was accepted,
        // including the ones whose constructs the analyzer has no representation
        // of. MATLAB had the mirror-image defect -- a set naming no service
        // distribution, which refused everything.
        assertTrue(new SolverBA(cyclicExp(4, 1.0, 0.8, 0.6)).supports(
                cyclicExp(4, 1.0, 0.8, 0.6)));

        // A renewal law is admissible whatever its higher moments: the analyzer
        // reads sn.rates and sn.visits and nothing else.
        Network erl = new Network("baErl");
        Queue q0 = new Queue(erl, "Q0", SchedStrategy.FCFS);
        Queue q1 = new Queue(erl, "Q1", SchedStrategy.FCFS);
        ClosedClass c = new ClosedClass(erl, "C", 3, q0);
        q0.setService(c, Erlang.fitMeanAndOrder(1.0, 2L));
        q1.setService(c, new Exp(0.8));
        erl.link(erl.serialRouting(q0, q1));
        assertTrue(new SolverBA(erl).supports(erl));
    }

    /** Single-class closed, one server each, MODULATED service at the first station. */
    private static Network baModulated() {
        Network m = new Network("baMap");
        Queue q0 = new Queue(m, "Q0", SchedStrategy.FCFS);
        Queue q1 = new Queue(m, "Q1", SchedStrategy.FCFS);
        ClosedClass c = new ClosedClass(m, "C", 3, q0);
        q0.setService(c, new MMPP2(0.2, 0.6, 0.1, 0.2));
        q1.setService(c, new Exp(0.8));
        m.link(m.serialRouting(q0, q1));
        return m;
    }

    @Test
    public void testAModulatedServiceProcessReachesMapamvaOnly() {
        // THIS ASSERTED THE OPPOSITE UNTIL 'mapamva' LANDED (2026-09-04), and the
        // premise rather than the rule is what changed. A modulated process was
        // out of the SOLVER envelope because every family here was a function of
        // the demands D = V./rates: the mean rate exists, so the utilization law
        // holds and the formula returns a number, but that number brackets a
        // DIFFERENT system, one whose successive services are independent.
        //
        // MAP-AMVA (Casale-Smirni, DSN 2009) is derived FOR the correlated model,
        // its LP variables being the per-phase QN(i,k) and UN(i,k), so the same
        // reasoning that refuses the others admits it. 'MAP'/'MMPP2' therefore
        // moved INTO the base envelope, and getMethodFeatureSet takes them back
        // from every other family. The direction is forced: a feature set refuses
        // a model for HAVING a construct and never for lacking one, so a law can
        // only be granted to one family by removing it from the rest.
        //
        // Keeping the old assertion would have made the new family unreachable
        // and hidden SolverBA from model.help() on exactly the models it was
        // written for, so the guarantee is re-pinned here at the level where it
        // now lives: the SOLVER accepts the model, and every family but mapamva
        // still refuses it BY NAME.
        Network m = baModulated();
        SolverBA solver = new SolverBA(m);
        assertTrue(solver.supports(m));

        List<String> supported = new java.util.ArrayList<String>();
        for (String name : solver.listValidMethods()) {
            if (solver.supportsModelMethod(name).isEmpty()) {
                supported.add(name);
            }
        }
        java.util.Collections.sort(supported);
        assertEquals(java.util.Arrays.asList("mapamva.lower", "mapamva.upper"), supported);
    }

    @Test
    public void testARenewalFamilyStillRefusesModulatedServiceByName() {
        // The half of the old assertion that must NOT weaken: a demand-parameterized
        // family asked for by name on a modulated model is refused, and the reason
        // names the offending law rather than silently bounding a different system.
        Network m = baModulated();
        SolverBA solver = new SolverBA(m);
        for (String name : new String[]{"gb.upper", "gb.lower", "aba.upper", "lr.upper"}) {
            assertFalse(solver.supportsModelMethod(name).isEmpty(),
                    name + " no longer refuses a modulated service process");
        }
    }
}
