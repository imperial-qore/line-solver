package jline.solvers.fluid.analyzers;

import java.util.ArrayList;
import java.util.List;

import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.Region;
import jline.lang.constant.DropStrategy;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Node;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.processes.Exp;
import jline.solvers.SolverOptions;
import jline.solvers.fluid.FluidResult;
import jline.solvers.fluid.SolverFluid;
import jline.util.matrix.Matrix;
import org.apache.commons.math3.util.FastMath;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * The fluid {@code dae} method: the min-normal closure solved as one system.
 *
 * <p>WHAT IS ACTUALLY BEING TESTED. {@code dae} is not a new closure -- it is the
 * SAME closure {@code minnormal} computes, taken from the same drift, the same
 * rate factors and the same Lyapunov equation. A test that only checked "does it
 * give a plausible number" would pass on an implementation that quietly
 * reproduced {@code minnormal} and threw the DAE away. The properties below are
 * the ones that distinguish the two:</p>
 *
 * <ul>
 *   <li>population conservation is an EQUATION here, not a consequence of the
 *       drift. {@code minnormal} conserves to integrator tolerance; this
 *       conserves to solver tolerance.</li>
 *   <li>a finite capacity region is solvable. Every other fluid method refuses
 *       one outright, because an ODE has nowhere to put a linear inequality.</li>
 *   <li>the transient runs on {@link jline.util.ode.Rodas} with a singular mass
 *       matrix, so conservation holds along the trajectory too.</li>
 * </ul>
 *
 * <p>and the refusals, which have to name the model feature rather than surface
 * from inside the Newton solve.</p>
 */
public class DaeAnalyzerTest {

    private static Network cqn(int n) {
        Network model = new Network("cqn");
        Delay d = new Delay(model, "Think");
        Queue q = new Queue(model, "Q1", SchedStrategy.PS);
        ClosedClass c = new ClosedClass(model, "C", n, d, 0);
        d.setService(c, new Exp(1.0));
        q.setService(c, new Exp(2.0));
        model.link(model.serialRouting(d, q));
        return model;
    }

    private static Network cqnRegion(int n, int cap, DropStrategy rule) {
        Network model = new Network("fcr");
        Delay d = new Delay(model, "Think");
        Queue q = new Queue(model, "Q1", SchedStrategy.PS);
        ClosedClass c = new ClosedClass(model, "C", n, d, 0);
        d.setService(c, new Exp(1.0));
        q.setService(c, new Exp(2.0));
        model.link(model.serialRouting(d, q));
        List<Node> nodes = new ArrayList<Node>();
        nodes.add(q);
        Region reg = model.addRegion(nodes);
        reg.setGlobalMaxJobs(cap);
        reg.setDropRule(c, rule);
        return model;
    }

    private static SolverOptions opts() {
        SolverOptions o = new SolverOptions();
        o.method = "dae";
        return o;
    }

    /**
     * Every case goes through SolverFluid rather than calling the analyzer
     * directly: the seed integration reads state that {@code runAnalyzer} sets
     * up, so a bare {@code analyze()} hands LSODA neq = 0.
     */
    private static FluidResult solve(Network model) {
        SolverFluid s = new SolverFluid(model, opts());
        s.getAvgQLen();
        return (FluidResult) s.result;
    }

    private static double sum(Matrix m) {
        double s = 0;
        for (int i = 0; i < m.getNumRows(); i++) {
            for (int j = 0; j < m.getNumCols(); j++) {
                s += m.get(i, j);
            }
        }
        return s;
    }

    // ------------------------------------------------------------- closure

    @Test
    public void testDaeReproducesTheMinNormalClosure() {
        // Same closure, different discharge. They must agree to the accuracy the
        // substitution route stops at: its outer loop stops at CoarseTol while
        // the Newton runs to options.tol, so the DAE answer is the tighter of
        // the two and the gap is the substitution error, not a disagreement.
        SolverOptions omn = new SolverOptions();
        omn.method = "minnormal";
        Matrix qmn = new SolverFluid(cqn(5), omn).getAvgQLen();
        Matrix qdae = new SolverFluid(cqn(5), opts()).getAvgQLen();
        for (int i = 0; i < qmn.getNumRows(); i++) {
            assertEquals(qmn.get(i, 0), qdae.get(i, 0), 1e-3);
        }
    }

    @Test
    public void testConservationIsAnEquationNotAConsequence() {
        // The whole point of the algebraic row. A drift whose rows sum to zero
        // conserves the population only as well as the integrator does; writing
        // it as an equation conserves it to the solve's own tolerance.
        int[] pops = {2, 5, 20};
        for (int p = 0; p < pops.length; p++) {
            Matrix qn = new SolverFluid(cqn(pops[p]), opts()).getAvgQLen();
            assertEquals(pops[p], sum(qn), 1e-8,
                    "population is not conserved at N=" + pops[p]);
        }
    }

    @Test
    public void testSaturationIsReachedWithoutAKinkProbe() {
        // At N=20 the queue is saturated and the first-order fixed point sits
        // exactly on the sigma2=0 kink of min(n,c) -- the configuration
        // MinNormalAnalyzer needs its two-sided Jacobian probe for. The DAE
        // seeds the variance POSITIVE and never adopts that iterate, so it must
        // simply solve. Think time 1 at throughput 2 holds 2 jobs.
        Matrix qn = new SolverFluid(cqn(20), opts()).getAvgQLen();
        assertEquals(2.0, qn.get(0, 0), 1e-6);
        assertEquals(18.0, qn.get(1, 0), 1e-6);
    }

    @Test
    public void testTheAnalyzerReportsAConvergedFixedPoint() {
        FluidResult fr = solve(cqn(5));
        assertTrue(fr.daeConverged, "the simultaneous solve did not reach options.tol");
        // the claim is that it reached the tolerance it was GIVEN; the JAR's
        // default options.tol is 1e-4, and pinning a constant here would only
        // record that default rather than the convergence
        assertTrue(fr.daeResidual < opts().tol, "residual " + fr.daeResidual);
        assertTrue(fr.daeConservation < 1e-8, "conservation " + fr.daeConservation);
    }

    // ------------------------------------------------- finite capacity region

    @Test
    public void testACapacityRegionBindsAndConserves() {
        // SolverFLD refuses a finite capacity region for every other method.
        // Here the cap is an algebraic equation, so the region sits exactly ON
        // it and the mass that does not fit is held in the waiting queue rather
        // than vanishing.
        int[] caps = {3, 2};
        for (int c = 0; c < caps.length; c++) {
            FluidResult fr = solve(cqnRegion(8, caps[c], DropStrategy.WaitingQueue));
            assertEquals(1, fr.daeCapacityActive.length, "the cap should bind");
            // the region holds exactly its cap
            assertEquals(caps[c], fr.QN.get(1, 0), 1e-6);
            // and the population is still all there: at a station, or blocked
            assertEquals(8.0, sum(fr.QN) + fr.daeBlocked, 1e-6);
            // a waiting queue drains at a finite positive rate
            assertTrue(fr.daeDrain.get(0, 0) > 0, "the throttle must be positive");
        }
    }

    @Test
    public void testACapThePopulationCannotReachNeverBinds() {
        // A region capped above the whole closed population is not a constraint,
        // and must be pruned rather than left in the active-set loop being
        // tested on every pass.
        FluidResult fr = solve(cqnRegion(8, 50, DropStrategy.WaitingQueue));
        assertEquals(0, fr.daeCapacityActive.length);
        assertEquals(0.0, fr.daeBlocked, 1e-12);
        assertEquals(8.0, sum(fr.QN), 1e-8);
    }

    @Test
    public void testADropRuleThatIsNotAWaitingQueueIsRefusedByName() {
        // Only a waiting queue conserves the population and throttles the flow.
        // DROP destroys the job, which is a different event set, not a
        // constraint on this drift.
        final Network model = cqnRegion(8, 3, DropStrategy.Drop);
        RuntimeException e = assertThrows(RuntimeException.class,
                new org.junit.jupiter.api.function.Executable() {
                    public void execute() {
                        solve(model);
                    }
                });
        assertTrue(e.getMessage().contains("waiting queue"),
                "the refusal must name the rule: " + e.getMessage());
    }

    @Test
    public void testTheRegionTransientHoldsTheCapAndLocatesTheCrossing() {
        // The transient under a cap is a HYBRID DAE: the system switches every time
        // the region fills or drains. Integrating with the binding set frozen would
        // report the unconstrained path through a cap the model declares, so what is
        // tested is that the path NEVER exceeds the cap, that the crossing is
        // located rather than stepped over, and that the trajectory ends where the
        // steady-state solve says it should.
        final Network model = cqnRegion(8, 5, DropStrategy.WaitingQueue);
        final SolverOptions o = opts();
        o.timespan = new double[]{0.0, 40.0};
        SolverFluid s = new SolverFluid(model, o);
        s.getTranAvg();
        FluidResult fr = (FluidResult) s.result;
        Matrix region = fr.QNt[1][0];
        double peak = 0;
        for (int i = 0; i < region.getNumRows(); i++) {
            peak = FastMath.max(peak, region.get(i, 0));
        }
        assertTrue(peak <= 5.0 + 1e-6, "the trajectory left the feasible set at " + peak);
        boolean activated = false;
        for (int i = 0; i < fr.daeCapacitySwitches.size(); i++) {
            activated |= fr.daeCapacitySwitches.get(i)[2] == 1;
        }
        assertTrue(activated, "the cap was reached and must have been activated");
        // and the two formulations describe one model: the drain RATE the steady
        // state solves for and the admitted FLOW the transient carries agree at the
        // fixed point
        FluidResult ss = solve(cqnRegion(8, 5, DropStrategy.WaitingQueue));
        assertEquals(ss.QN.get(1, 0), region.get(region.getNumRows() - 1, 0), 1e-3);
    }

    // ------------------------------------------------------------ station buffers

    private static Network capped(int n, int cap, double mu) {
        Network model = new Network("capped");
        Delay d = new Delay(model, "Think");
        Queue q = new Queue(model, "Q1", SchedStrategy.FCFS);
        q.setCapacity(cap);
        ClosedClass c = new ClosedClass(model, "C", n, d, 0);
        d.setService(c, new Exp(1.0));
        q.setService(c, new Exp(mu));
        model.link(model.serialRouting(d, q));
        return model;
    }

    @Test
    public void testAStationBufferBindsAndHoldsTheBlockedJobUpstream() {
        // A station buffer is the one-station case of the same row -- but NOT of the
        // same model. LINE refuses to lose a closed job (State.arrivalIsLost) and
        // disables the upstream departure instead, so the blocked mass is still AT
        // the upstream station and still counted there: the station queues sum to
        // the whole population and nothing is staged. That is the opposite of a
        // region, whose blocked jobs are reported separately.
        FluidResult fr = solve(capped(8, 2, 2.0));
        assertEquals(1, fr.daeCapacityActive.length, "the buffer should bind");
        assertEquals(2.0, fr.QN.get(1, 0), 1e-6);
        assertEquals(8.0, sum(fr.QN), 1e-6, "a held job is at the upstream station");
        assertEquals(0.0, fr.daeBlocked, 1e-12, "a station buffer has no waiting room");
        // the multiplier is the FRACTION of upstream completions the cap admits
        double alpha = fr.daeDrain.get(0, 0);
        assertTrue(alpha > 0 && alpha < 1, "the admitted fraction must be in (0,1): " + alpha);
    }

    @Test
    public void testAStationBufferThePopulationCannotReachIsInert() {
        // refreshCapacity derives a finite classcap for every closed model, so a cap
        // that cannot bind is the common case, not a corner one.
        FluidResult fr = solve(capped(8, 50, 2.0));
        assertEquals(0, fr.daeCapacityActive.length);
        assertEquals(8.0, sum(fr.QN), 1e-8);
    }

    @Test
    public void testAnArrivalIsLostAtAFullBufferOfAnOpenClass() {
        // The same predicate that HOLDS a closed job LOSES an open one: the external
        // stream is memoryless, so a job that finds the buffer full never enters.
        // In overload the fluid loss rate is exact -- what gets through is the
        // server -- and the multiplier is the admitted fraction.
        Network model = new Network("loss");
        Source src = new Source(model, "S");
        Queue q = new Queue(model, "Q1", SchedStrategy.FCFS);
        Sink k = new Sink(model, "K");
        q.setCapacity(5);
        OpenClass c = new OpenClass(model, "C");
        src.setArrival(c, new Exp(4.0));
        q.setService(c, new Exp(2.0));
        model.link(model.serialRouting(src, q, k));
        FluidResult fr = solve(model);
        assertEquals(1, fr.daeCapacityActive.length);
        assertEquals(5.0, fr.QN.get(1, 0), 1e-6);
        assertEquals(2.0, fr.TN.get(1, 0), 1e-4, "the carried flow is the server");
        assertEquals(0.5, fr.daeDrain.get(0, 0), 1e-3, "half the arrivals are admitted");
    }

    @Test
    public void testDefaultResolvesToDaeOnABlockedModel() {
        // The refusal below is correct and its advice was useless: it told the
        // caller to type the one method the resolution could have picked itself.
        // "default" now stands for "dae" wherever a buffer or a region BINDS and
        // the dae route accepts the model, and for nothing else.
        Network blocked = capped(4, 2, 2.0);
        SolverOptions def = new SolverOptions();
        def.method = "default";
        assertTrue(SolverFluid.blockedResolvesToDae(blocked, blocked.getStruct(false), def));
        assertEquals(2.0, new SolverFluid(blocked, def).getAvgQLen().get(1, 0), 1e-6,
                "the buffer holds exactly its cap");

        // A cap the population cannot reach is not a constraint, so the
        // resolution must be untouched there.
        Network free = capped(4, 50, 2.0);
        SolverOptions def2 = new SolverOptions();
        def2.method = "default";
        assertFalse(SolverFluid.blockedResolvesToDae(free, free.getStruct(false), def2));
    }

    @Test
    public void testEveryOtherFluidMethodRefusesABindingBuffer() {
        // Nothing in the fluid tree reads sn.cap, so a capped station was integrated
        // as an unbounded one and the table reported more jobs in the buffer than
        // the buffer holds. Only the route that can enforce the cap may accept it.
        final Network model = capped(4, 2, 2.0);
        final SolverOptions omn = new SolverOptions();
        omn.method = "minnormal";
        RuntimeException e = assertThrows(RuntimeException.class,
                new org.junit.jupiter.api.function.Executable() {
                    public void execute() {
                        new SolverFluid(model, omn).getAvgQLen();
                    }
                });
        assertTrue(e.getMessage().contains("dae"), "the refusal must name the route: "
                + e.getMessage());
        assertEquals(2.0, new SolverFluid(capped(4, 2, 2.0), opts()).getAvgQLen().get(1, 0), 1e-6);
    }

    @Test
    public void testTwoClassBuffersBindAtOnce() {
        // TWO CAPS BINDING TOGETHER used to be refused outright: one region carried a
        // single throttle, so a second equality had no control to satisfy it. Each
        // active row now carries its own multiplier.
        Network model = new Network("mc");
        Delay d = new Delay(model, "Think");
        Queue q = new Queue(model, "Q1", SchedStrategy.PS);
        ClosedClass a = new ClosedClass(model, "A", 6, d, 0);
        ClosedClass b = new ClosedClass(model, "B", 6, d, 0);
        d.setService(a, new Exp(1.0));
        d.setService(b, new Exp(1.0));
        q.setService(a, new Exp(2.0));
        q.setService(b, new Exp(2.0));
        q.setClassCap(a, 2);
        q.setClassCap(b, 3);
        model.link(model.serialRouting(d, q));
        FluidResult fr = solve(model);
        assertEquals(2, fr.daeCapacityActive.length, "both class buffers must bind");
        for (int c = 0; c < fr.daeCapacityB.getNumRows(); c++) {
            assertEquals(fr.daeCapacityB.get(c, 0), fr.daeCapacityValue.get(c, 0), 1e-6);
        }
        assertTrue(fr.daeConverged, "the two-multiplier system must converge");
    }

    @Test
    public void testTheStationTotalImpliedByItsClassBuffersIsPruned() {
        // LINE derives the station total from the per-class buffers, so a two-class
        // station capped 3 and 3 also declares a total of 6 -- exactly the sum of the
        // two class rows. All three would bind together with rank 2, and the
        // multipliers would be one arbitrary point of a line of solutions.
        Network model = new Network("mc2");
        Delay d = new Delay(model, "Think");
        Queue q = new Queue(model, "Q1", SchedStrategy.PS);
        ClosedClass a = new ClosedClass(model, "A", 6, d, 0);
        ClosedClass b = new ClosedClass(model, "B", 6, d, 0);
        d.setService(a, new Exp(1.0));
        d.setService(b, new Exp(1.0));
        q.setService(a, new Exp(2.0));
        q.setService(b, new Exp(2.0));
        q.setClassCap(a, 3);
        q.setClassCap(b, 3);
        model.link(model.serialRouting(d, q));
        FluidResult fr = solve(model);
        assertEquals(2, fr.daeCapacityLabel.length, "the implied station total must be pruned");
        for (int c = 0; c < fr.daeCapacityLabel.length; c++) {
            assertTrue(fr.daeCapacityLabel[c].contains("class"), fr.daeCapacityLabel[c]);
        }
    }

    // ------------------------------------------------------------- refusals

    @Test
    public void testTheStateCapIsRefusedByName() {
        final SolverOptions o = opts();
        o.config.dae_maxstate = 1;
        final Network model = cqn(5);
        RuntimeException e = assertThrows(RuntimeException.class,
                new org.junit.jupiter.api.function.Executable() {
                    public void execute() {
                        new SolverFluid(model, o).getAvgQLen();
                    }
                });
        assertTrue(e.getMessage().contains("dae_maxstate"), e.getMessage());
    }

    @Test
    public void testDaeIsReachableThroughTheSolverFrontDoor() {
        // The dispatch, not the analyzer: options.method='dae' must resolve, and
        // 'dae' must appear in listValidMethods or the gate rejects it.
        String[] methods = new SolverFluid(cqn(5)).listValidMethods();
        boolean found = false;
        for (int i = 0; i < methods.length; i++) {
            found |= "dae".equals(methods[i]);
        }
        assertTrue(found, "dae is not in listValidMethods");
        Matrix qn = new SolverFluid(cqn(5), opts()).getAvgQLen();
        assertEquals(5.0, sum(qn), 1e-8);
    }

    @Test
    public void testTheFeatureSetDeclaresRegionOnlyForDae() {
        // The DAE method WIDENS the static envelope, which is the one direction
        // the coarse gate cannot express: every other method must go on
        // rejecting a finite capacity region.
        SolverFluid s = new SolverFluid(cqn(5));
        assertTrue(s.getMethodFeatureSet("dae").inspectFeature("Region"),
                "dae must declare Region or the model is rejected before the method is consulted");
        assertTrue(!s.getMethodFeatureSet("minnormal").inspectFeature("Region"),
                "minnormal must keep rejecting a finite capacity region");
    }

    // ----------------------------------------------------------------------
    // The non-hyperbolic fallback ladder: minnormal -> dae -> first order
    // ----------------------------------------------------------------------
    //
    // FluidLyapunov throws when the drift Jacobian at the converged mean is not
    // exponentially stable on range(D). The dominant case is NEUTRAL rather than
    // unstable, and it is an artifact of the alternation: MinNormalAnalyzer must
    // start at sigma2 = 0, where min(n,c) has no derivative, so a saturated or
    // balanced model's first-order fixed point lands on the kink and sits on a
    // continuum of equilibria. DaeAnalyzer seeds the variance positive and never
    // adopts sigma2 = 0, so the same closure has an isolated, hyperbolic fixed
    // point there.

    /**
     * Two identical stations in a closed cycle: the fluid drift is DEGENERATE.
     * Every state with both populations at or above the server count is a fluid
     * equilibrium, so a first-order method has no reason to prefer one point of
     * that continuum over another and lands wherever the integrator stopped.
     */
    private static Network balancedCycle(int n, int nservers, SchedStrategy sched) {
        Network model = new Network("balanced");
        Queue q1 = new Queue(model, "Q1", sched);
        Queue q2 = new Queue(model, "Q2", sched);
        q1.setNumberOfServers(nservers);
        q2.setNumberOfServers(nservers);
        ClosedClass c = new ClosedClass(model, "C", n, q1, 0);
        q1.setService(c, new Exp(1.0));
        q2.setService(c, new Exp(1.0));
        List<Node> route = new ArrayList<Node>();
        route.add(q1);
        route.add(q2);
        model.link(model.serialRouting(route));
        return model;
    }

    @Test
    public void testADeclinedMinNormalLandsOnDaeAndOnTheExactAnswer() {
        // By symmetry the exact answer splits the population evenly, which is
        // what the closure gives once the degeneracy is broken. The first-order
        // fallback this replaces returned [9 1] at N=10 -- a point of the
        // continuum, not the mean.
        int[] pop = {4, 6, 6, 10, 10, 6, 6};
        int[] servers = {1, 1, 2, 1, 2, 1, 2};
        SchedStrategy[] sched = {SchedStrategy.PS, SchedStrategy.PS, SchedStrategy.PS,
                SchedStrategy.PS, SchedStrategy.PS, SchedStrategy.FCFS, SchedStrategy.FCFS};
        for (int k = 0; k < pop.length; k++) {
            SolverOptions o = new SolverOptions();
            o.method = "minnormal";
            // the dae rung converges to options.tol, whose default here is 1e-4;
            // tighten it so the assertion below measures the closure and not the
            // Newton stopping rule
            o.tol = 1e-8;
            Matrix qn = new SolverFluid(balancedCycle(pop[k], servers[k], sched[k]), o).getAvgQLen();
            String at = "N=" + pop[k] + " c=" + servers[k] + " " + sched[k];
            assertEquals(pop[k] / 2.0, qn.get(0, 0), 1e-6, at);
            assertEquals(pop[k] / 2.0, qn.get(1, 0), 1e-6, at);
        }
    }

    @Test
    public void testAGenuinelyUnstableFixedPointWalksPastDaeToFirstOrder() {
        // An overloaded open station has no stationary distribution at all, so no
        // closure has a stationary covariance there and dae declines on the same
        // exception minnormal did. The mean is still reported, by the first-order
        // method, which is the whole reason the last rung exists.
        Network model = new Network("overloaded");
        Source src = new Source(model, "Source");
        Queue q = new Queue(model, "Q1", SchedStrategy.FCFS);
        Sink snk = new Sink(model, "Sink");
        OpenClass c = new OpenClass(model, "C", 0);
        src.setArrival(c, new Exp(1.1));
        q.setService(c, new Exp(1.0));
        List<Node> route = new ArrayList<Node>();
        route.add(src);
        route.add(q);
        route.add(snk);
        model.link(model.serialRouting(route));

        SolverOptions o = new SolverOptions();
        o.method = "minnormal";
        // reaching a finite answer at all IS the assertion: both closure rungs
        // decline, and without the last rung this throws
        Matrix qn = new SolverFluid(model, o).getAvgQLen();
        assertTrue(Double.isFinite(sum(qn)), "the first-order rung must still report a mean");
    }

    @Test
    public void testTheDaeRungIsDeclinedInAdvanceForTheFeaturesItCannotTake() {
        // FluidDaeApplicable is the STATIC difference set between the closures. A
        // rung entered only to be refused a moment later would spend a whole seed
        // integration to learn what the model already declares.
        Network dps = new Network("dps");
        Delay d = new Delay(dps, "Think");
        Queue q = new Queue(dps, "Q1", SchedStrategy.DPS);
        ClosedClass c = new ClosedClass(dps, "C", 4, d, 0);
        d.setService(c, new Exp(1.0));
        q.setService(c, new Exp(2.0));
        List<Node> route = new ArrayList<Node>();
        route.add(d);
        route.add(q);
        dps.link(dps.serialRouting(route));
        SolverOptions o = new SolverOptions();
        String why = jline.solvers.fluid.moments.FluidDaeApplicable
                .reasonToDecline(dps.getStruct(), o);
        assertTrue(why != null && why.contains("covariance between"), String.valueOf(why));

        // the simultaneous solve is quartic where one Lyapunov solve is cubic, so
        // it carries its own cap, lower than moment_maxstate
        SolverOptions capped = new SolverOptions();
        capped.config.dae_maxstate = 1;
        why = jline.solvers.fluid.moments.FluidDaeApplicable
                .reasonToDecline(cqn(5).getStruct(), capped);
        assertTrue(why != null && why.contains("dae_maxstate"), String.valueOf(why));

        // and the model the ladder does take
        why = jline.solvers.fluid.moments.FluidDaeApplicable
                .reasonToDecline(balancedCycle(6, 1, SchedStrategy.PS).getStruct(),
                        new SolverOptions());
        assertEquals(null, why);
    }
}
