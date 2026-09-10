package jline.solvers.auto;

import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.processes.Erlang;
import jline.lang.processes.Exp;
import jline.lang.processes.MMPP2;
import jline.lang.constant.DropStrategy;
import jline.solvers.NetworkSolver;
import jline.solvers.ctmc.SolverCTMC;
import jline.solvers.mva.SolverMVA;
import jline.solvers.nc.SolverNC;
import jline.lang.RoutingMatrix;
import org.junit.jupiter.api.Test;

import java.lang.reflect.Method;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * SolverAUTO ranking: LDES for simulation, MVA over NC over MAM for
 * analytical solvers, Fluid where a smooth or high-load answer is wanted, and
 * per-metric capability routing. Mirrors
 * line-test.git/test/testsMisc/test_auto_solver_selection.m.
 */
public class SolverAUTOSelectionTest {

    private static Network closedPS(int n) {
        Network model = new Network("closedPS");
        Delay d = new Delay(model, "d");
        Queue q = new Queue(model, "q", SchedStrategy.PS);
        ClosedClass c = new ClosedClass(model, "c", n, d);
        d.setService(c, new Exp(1.0));
        q.setService(c, new Exp(0.5));
        model.link(model.serialRouting(d, q));
        return model;
    }

    private static Network openMAP() {
        Network model = new Network("openMAP");
        Source s = new Source(model, "src");
        Queue q = new Queue(model, "q", SchedStrategy.FCFS);
        Sink k = new Sink(model, "sink");
        OpenClass c = new OpenClass(model, "c");
        s.setArrival(c, new MMPP2(0.2, 0.6, 0.1, 0.2));
        q.setService(c, new Exp(1.0));
        model.link(model.serialRouting(s, q, k));
        return model;
    }

    private static Network openExp() {
        Network model = new Network("openExp");
        Source s = new Source(model, "src");
        Queue q = new Queue(model, "q", SchedStrategy.FCFS);
        Sink k = new Sink(model, "sink");
        OpenClass c = new OpenClass(model, "c");
        s.setArrival(c, new Exp(3.0));
        q.setService(c, new Exp(4.0));
        model.link(model.serialRouting(s, q, k));
        return model;
    }

    private static String selected(Network model, String mode) {
        SolverAUTO solver = new SolverAUTO(model, mode);
        solver.getAvgTable();
        return solver.getSelectedSolverName();
    }

    @Test
    public void testClosedProductFormPrefersMVAOverNC() {
        assertEquals("SolverMVA", selected(closedPS(5), "default"));
    }

    @Test
    public void testLargePopulationPrefersFluid() {
        assertEquals("SolverFluid", selected(closedPS(60), "default"));
    }

    @Test
    public void testOpenProductFormPrefersMVA() {
        // Open classes carry njobs=Inf, which must not count as population: a
        // finite cast would send every open model to the Fluid branch.
        assertEquals("SolverMVA", selected(openExp(), "default"));
        assertEquals(0, new ModelAnalyzer(openExp()).getTotalJobs());
    }

    @Test
    public void testAutocorrelatedArrivalsPreferMAM() {
        assertEquals("SolverMAM", selected(openMAP(), "default"));
    }

    @Test
    public void testExactModeReturnsExactSolver() {
        assertEquals("SolverNC", selected(closedPS(5), "exact"));
    }

    @Test
    public void testSimModePrefersLDESOverSSA() {
        assertEquals("SolverLDES", selected(closedPS(5), "sim"));
    }

    @Test
    public void testJMTIsNeverAnAutoCandidate() {
        // LDES subsumes SolverJMT's feature set, so automatic selection must
        // never reach the external simulator in any mode.
        String[] modes = new String[]{"default", "sim", "fast", "accurate", "exact"};
        for (String mode : modes) {
            assertFalse(selected(closedPS(5), mode).equals("SolverJMT"),
                    "mode " + mode + " selected SolverJMT");
        }
        assertFalse(SolverAUTO.getFeatureSet().inspectFeature("Normal"),
                "AUTO must not advertise the JMT-only Normal distribution");
    }

    @Test
    public void testFastModePrefersMVA() {
        assertEquals("SolverMVA", selected(closedPS(5), "fast"));
    }

    @Test
    public void testAccurateModePrefersFluid() {
        assertEquals("SolverFluid", selected(closedPS(5), "accurate"));
    }

    @Test
    public void testDelayIsNotCountedAsMultiServer() {
        // A Delay has infinite servers; counting it as multiserver made every
        // model with a think time miss the exact product-form branch.
        assertEquals("SolverNC", selected(closedPS(5), "exact"));
    }

    /**
     * The per-metric ranking lives in the private chooseSolver, which the
     * getAvg* path reaches through delegate; reflection asserts the mapping
     * itself rather than the accessor that happens to exercise it.
     *
     * @param model the model under analysis
     * @param mode the AUTO selection token
     * @param methodName the accessor being routed
     * @return the name of the solver AUTO would delegate that accessor to
     */
    private static String routedTo(Network model, String mode, String methodName) throws Exception {
        SolverAUTO solver = new SolverAUTO(model, mode);
        Method chooser = SolverAUTO.class.getDeclaredMethod("chooseSolver", String.class);
        chooser.setAccessible(true);
        NetworkSolver chosen = (NetworkSolver) chooser.invoke(solver, methodName);
        return chosen.getName();
    }

    @Test
    public void testDistributionsAndTransientsRouteToFluid() throws Exception {
        assertEquals("SolverFluid", routedTo(closedPS(5), "default", "getCdfRespT"));
        assertEquals("SolverFluid", routedTo(closedPS(5), "default", "getTranAvg"));
        assertEquals("SolverFluid", routedTo(closedPS(5), "default", "getSensitivityTable"));
    }

    @Test
    public void testProbabilitiesRouteToNCAndTransientProbabilitiesToCTMC() throws Exception {
        assertEquals("SolverNC", routedTo(closedPS(5), "default", "getProb"));
        assertEquals("SolverCTMC", routedTo(closedPS(5), "default", "getTranProb"));
    }

    @Test
    public void testSamplingRoutesToSSAAndAggregateSamplingToLDES() throws Exception {
        assertEquals("SolverSSA", routedTo(closedPS(5), "default", "sample"));
        assertEquals("SolverLDES", routedTo(closedPS(5), "default", "sampleAggr"));
    }

    @Test
    public void testLossTablesRouteToLDESAndMomentsToMVA() throws Exception {
        assertEquals("SolverLDES", routedTo(closedPS(5), "default", "getAvgLossTable"));
        assertEquals("SolverMVA", routedTo(closedPS(5), "default", "getMomentTable"));
    }

    @Test
    public void testSimModePrefersLDESForEveryMetricButSampling() throws Exception {
        assertEquals("SolverLDES", routedTo(closedPS(5), "sim", "getCdfRespT"));
        assertEquals("SolverSSA", routedTo(closedPS(5), "sim", "sample"));
    }

    @Test
    public void testFastAndAccurateModesRouteConsistently() throws Exception {
        assertEquals("SolverMVA", routedTo(closedPS(5), "fast", "getAvgTable"));
        assertEquals("SolverFluid", routedTo(closedPS(5), "accurate", "getAvgTable"));
    }

    /**
     * Eight PS stations, N=400, Erlang-5 service: logNstates is about 200, so
     * the chain cannot be built even though every feature is supported.
     *
     * @return a model whose CTMC state space is intractable
     */
    private static Network intractableCTMC() {
        Network model = new Network("intractableCTMC");
        Queue[] st = new Queue[8];
        for (int i = 0; i < 8; i++) {
            st[i] = new Queue(model, "q" + i, SchedStrategy.PS);
        }
        ClosedClass c = new ClosedClass(model, "c", 400, st[0]);
        for (int i = 0; i < 8; i++) {
            st[i].setService(c, Erlang.fitMeanAndOrder(1.0, 5));
        }
        model.link(model.serialRouting(st[0], st[1], st[2], st[3], st[4], st[5], st[6], st[7]));
        return model;
    }

    /**
     * Multiclass FCFS with unequal per-class rates: no product-form solution.
     *
     * @return a model on which the exact methods cannot run
     */
    private static Network nonProductForm() {
        return nonProductForm(2);
    }

    /**
     * Same non-product-form topology at a chosen per-class population.
     *
     * @param n jobs per class
     * @return the model
     */
    private static Network nonProductForm(int n) {
        Network model = new Network("nonPF");
        Delay d = new Delay(model, "d");
        Queue q = new Queue(model, "q", SchedStrategy.FCFS);
        ClosedClass c1 = new ClosedClass(model, "c1", n, d);
        ClosedClass c2 = new ClosedClass(model, "c2", n, d);
        d.setService(c1, new Exp(1.0));
        d.setService(c2, new Exp(1.0));
        q.setService(c1, new Exp(2.0));
        q.setService(c2, new Exp(0.5));
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(c1, c1, d, q, 1.0);
        P.set(c1, c1, q, d, 1.0);
        P.set(c2, c2, d, q, 1.0);
        P.set(c2, c2, q, d, 1.0);
        model.link(P);
        return model;
    }

    @Test
    public void testSmallPopulationsTakeAnExactSolver() throws Exception {
        // Total population 2n: at or below EXACT_POPULATION_MAX an exact solver
        // is preferred, and on this non-product-form model only CTMC qualifies.
        assertEquals("SolverCTMC", routedTo(nonProductForm(1), "default", "getAvgTable"));
        assertEquals("SolverCTMC", routedTo(nonProductForm(2), "default", "getAvgTable"));
        // Above it the ranking falls back to the approximation.
        assertEquals("SolverMVA", routedTo(nonProductForm(3), "default", "getAvgTable"));
    }

    @Test
    public void testSmallProductFormPopulationStillPrefersMVA() throws Exception {
        // MVA leads NC leads CTMC when all three are exact for the model.
        assertEquals("SolverMVA", routedTo(closedPS(3), "default", "getAvgTable"));
    }

    @Test
    public void testMethodGateRejectsExactWithoutProductForm() {
        Network model = nonProductForm();
        assertFalse(model.hasProductFormSolution());
        assertFalse(new SolverMVA(model).supportsModelMethod("exact").isEmpty());
        assertFalse(new SolverNC(model).supportsModelMethod("exact").isEmpty());
        // The other NC methods fall back to comom, so they stay admissible.
        assertTrue(new SolverNC(model).supportsModelMethod("default").isEmpty());
    }

    @Test
    public void testExactModeRoutesAroundANonProductFormModel() throws Exception {
        assertEquals("SolverCTMC", routedTo(nonProductForm(), "exact", "getAvgTable"));
        assertEquals("SolverNC", routedTo(closedPS(4), "exact", "getAvgTable"));
    }

    @Test
    public void testCTMCAccessorsAlwaysResolveToCTMC() {
        // State space and generator are CTMC-only concepts: AUTO exposes them
        // under the same names and never routes them through the ranking.
        Network model = closedPS(2);
        SolverAUTO auto = new SolverAUTO(model, "default");
        SolverCTMC.StateSpace ss = auto.getStateSpace();
        SolverCTMC.StateSpace ref = new SolverCTMC(model).getStateSpace();
        assertEquals(ref.stateSpace.getNumRows(), ss.stateSpace.getNumRows());
        assertEquals(ref.stateSpace.getNumCols(), ss.stateSpace.getNumCols());

        SolverCTMC.generatorResult gen = auto.getGenerator();
        assertEquals(ss.stateSpace.getNumRows(), gen.infGen.getNumRows());
        assertEquals(ss.stateSpace.getNumRows(), gen.infGen.getNumCols());

        assertTrue(auto.getSymbolicGenerator().eventFilt.size() > 0);
    }

    @Test
    public void testStateSpaceGateAcceptsSmallAndRefusesLarge() {
        assertTrue(SolverCTMC.isStateSpaceTractable(closedPS(3), SolverCTMC.defaultOptions()).ok);
        assertFalse(SolverCTMC.isStateSpaceTractable(intractableCTMC(), SolverCTMC.defaultOptions()).ok);
    }

    @Test
    public void testAutoSkipsCTMCWhenTheChainCannotBeBuilt() throws Exception {
        // getTranProb ranks CTMC alone; an untractable chain must not be proposed.
        assertEquals("SolverCTMC", routedTo(closedPS(3), "default", "getTranProb"));
        assertEquals("SolverFluid", routedTo(intractableCTMC(), "default", "getTranProb"));
    }

    @Test
    public void testAutocorrelatedArrivalsReachFluidUnderFastAndAccurate() throws Exception {
        // MAP and MMPP2 are in the Fluid feature set in all three codebases:
        // the ODE takes the stationary rate and conserves flow.
        assertEquals("SolverFluid", routedTo(openMAP(), "fast", "getAvgTable"));
        assertEquals("SolverFluid", routedTo(openMAP(), "accurate", "getAvgTable"));
    }

    /**
     * cqn_bas_blocking: two FCFS queues, one closed class of 2, a finite buffer
     * on the second and the BAS drop rule on the first.
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

    private static Set<String> familiesOf(String[] tokens) {
        Set<String> fams = new HashSet<String>();
        for (String t : tokens) {
            int dot = t.indexOf('.');
            fams.add(dot < 0 ? t : t.substring(0, dot));
        }
        return fams;
    }

    @Test
    public void testListValidMethodsIsGatedOnTheModel() {
        // listValidMethods USED TO BE the method name universe: it built every family,
        // asked each for its own list and kept everything, so a BAS-blocking
        // model was offered every NC method although SolverNC refuses it method
        // by method. It now applies the gate chooseSolverRanked already applies
        // before delegating.
        String[] valid = new SolverAUTO(basBlocking()).listValidMethods();
        Set<String> fams = familiesOf(valid);
        List<String> tokens = Arrays.asList(valid);
        for (String t : valid) {
            assertFalse(t.startsWith("nc."), "SolverNC refuses a binding buffer: " + t);
            assertFalse(t.startsWith("fluid.") && !t.endsWith(".dae") && !t.endsWith(".default"),
                    "the fluid drift ignores sn.cap outside 'dae': " + t);
        }
        // 'fluid.default' SURVIVES, and it names the same run as 'fluid.dae':
        // on a blocked model SolverFluid resolves 'default' to 'dae' (see
        // blockedResolvesToDae), so gating the literal name would drop one of
        // two method names for one method.
        assertTrue(tokens.contains("fluid.default"));
        assertTrue(tokens.contains("fluid.dae"));
        // A family whose every method is refused loses its bare method name too.
        assertFalse(fams.contains("nc"));
        // MVA stays: solver_mva_analyzer routes a BAS model to solver_sqd.
        assertTrue(tokens.contains("mva.sqd"));
        // And the state-space solvers, which represent the buffer exactly.
        assertTrue(fams.contains("ctmc"));
        assertTrue(fams.contains("ssa"));
        assertTrue(fams.contains("ldes"));
    }

    @Test
    public void testListAllMethodsIsTheUnnarrowedTokenUniverse() {
        // The name check must gate on THIS list, so that asking for a method a
        // candidate refuses gets the candidate's own reason rather than a flat
        // "unsupported by this solver".
        SolverAUTO auto = new SolverAUTO(basBlocking());
        List<String> every = Arrays.asList(auto.listAllMethods());
        List<String> valid = Arrays.asList(auto.listValidMethods());
        assertTrue(every.containsAll(valid));
        assertTrue(every.size() > valid.size());
        assertTrue(every.contains("nc.exact"));
        assertFalse(valid.contains("nc.exact"));
    }

    @Test
    public void testAnUnconstrainedModelKeepsTheWideList() {
        // The narrowing must be model-sensitive, not a blanket trim: drop the cap
        // and the product-form and fluid families come back.
        String[] valid = new SolverAUTO(closedPS(2)).listValidMethods();
        Set<String> fams = familiesOf(valid);
        assertTrue(fams.contains("nc"));
        int fluid = 0;
        for (String t : valid) {
            if (t.startsWith("fluid.")) {
                fluid++;
            }
        }
        assertTrue(fluid > 2, "fluid methods on an unconstrained model: " + fluid);
    }
}
