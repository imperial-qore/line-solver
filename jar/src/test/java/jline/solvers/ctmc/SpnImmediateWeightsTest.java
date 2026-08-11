package jline.solvers.ctmc;

import jline.VerboseLevel;
import jline.lang.ClosedClass;
import jline.lang.Mode;
import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.RoutingMatrix;
import jline.lang.constant.TimingStrategy;
import jline.lang.nodes.Place;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.nodes.Transition;
import jline.lang.processes.Det;
import jline.lang.processes.Exp;
import jline.lang.processes.Immediate;
import jline.solvers.NetworkAvgTable;
import org.junit.jupiter.api.Test;

import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * GSPN immediate-transition semantics and phase-type expansion of non-Markovian
 * firing times. Mirror of the native-Python
 * {@code python/tests/test_spn_immediate_weights.py}.
 *
 * Neither defect covered here was visible to the existing SPN tests, because
 * those exercise only a single immediate mode (never a conflict) and only
 * exponential firing times.
 *
 * 1. FREE-CHOICE CONFLICT. When two immediate modes are enabled by the same
 *    marking, GSPN semantics require the branching to follow the ratio of their
 *    firing weights. LINE latches a transition's servers with an ENABLE event
 *    emitted at GlobalConstants.Immediate, the same scale as an immediate
 *    firing, so a mode that happened to be latched first could fire before its
 *    competitor was latched at all. Gating the firing on the latched server
 *    count made the branching follow the latching order: with weights 1 and 3
 *    the split came out at 0.34375 instead of 0.25. State.java now gates an
 *    immediate firing on the marking (markDegreeM), as MATLAB does at
 *    afterGlobalEvent.m:315-330.
 *
 * 2. NON-MARKOVIAN FIRING TIMES. SnNonmarkovToPh had no Transition branch at
 *    all, so an SPN firing time kept the single nominal phase installed by the
 *    node layer and the model was silently solved as if the firing time were
 *    exponential with the same mean. MATLAB (sn_nonmarkov_toph.m:152-249) and
 *    native Python (transforms.py:1150-1211) both converted it.
 */
public class SpnImmediateWeightsTest {

    private static final double LAMBDA = 0.5;
    private static final double MU = 1.0;
    private static final double RHO = LAMBDA / MU;
    private static final int CAP = 12;

    // ---------------------------------------------------------------
    // 1. Free-choice conflict between two immediate modes
    // ---------------------------------------------------------------

    /**
     * One token cycling through a free-choice conflict:
     * P1 -> (T1 | T2) -> P2 | P3 -> (T3 | T4) -> P1.
     *
     * In steady state every cycle passes through exactly one of T1, T2, so the
     * throughputs of P2 and P3 are the branch probabilities scaled by the common
     * cycle rate. Their ratio is the firing-weight ratio, independently of which
     * immediate server happens to be latched first.
     */
    private static Network buildFreeChoice(double weight1, double weight2) {
        Network model = new Network("spn_free_choice");
        Place p1 = new Place(model, "P1");
        Place p2 = new Place(model, "P2");
        Place p3 = new Place(model, "P3");
        Transition t1 = new Transition(model, "T1");
        Transition t2 = new Transition(model, "T2");
        Transition t3 = new Transition(model, "T3");
        Transition t4 = new Transition(model, "T4");
        ClosedClass jc = new ClosedClass(model, "C", 1, p1, 0);

        Mode m1 = t1.addMode("m1");
        t1.setDistribution(m1, new Immediate());
        t1.setTimingStrategy(m1, TimingStrategy.IMMEDIATE);
        t1.setFiringWeights(m1, weight1);
        t1.setEnablingConditions(m1, jc, p1, 1);
        t1.setFiringOutcome(m1, jc, p2, 1);

        Mode m2 = t2.addMode("m2");
        t2.setDistribution(m2, new Immediate());
        t2.setTimingStrategy(m2, TimingStrategy.IMMEDIATE);
        t2.setFiringWeights(m2, weight2);
        t2.setEnablingConditions(m2, jc, p1, 1);
        t2.setFiringOutcome(m2, jc, p3, 1);

        Mode m3 = t3.addMode("m3");
        t3.setDistribution(m3, new Exp(1.0));
        t3.setEnablingConditions(m3, jc, p2, 1);
        t3.setFiringOutcome(m3, jc, p1, 1);

        Mode m4 = t4.addMode("m4");
        t4.setDistribution(m4, new Exp(1.0));
        t4.setEnablingConditions(m4, jc, p3, 1);
        t4.setFiringOutcome(m4, jc, p1, 1);

        RoutingMatrix R = model.initRoutingMatrix();
        R.set(jc, jc, p1, t1, 1.0);
        R.set(jc, jc, p1, t2, 1.0);
        R.set(jc, jc, t1, p2, 1.0);
        R.set(jc, jc, t2, p3, 1.0);
        R.set(jc, jc, p2, t3, 1.0);
        R.set(jc, jc, t3, p1, 1.0);
        R.set(jc, jc, p3, t4, 1.0);
        R.set(jc, jc, t4, p1, 1.0);
        model.link(R);
        p1.setState(1);
        p2.setState(0);
        p3.setState(0);
        return model;
    }

    /** Branch probability taken by T1, i.e. Tput(P2) / (Tput(P2) + Tput(P3)). */
    private static double branchProbability(double weight1, double weight2) {
        Network model = buildFreeChoice(weight1, weight2);
        SolverCTMC solver = new SolverCTMC(model, "verbose", VerboseLevel.SILENT);
        NetworkAvgTable table = solver.getAvgTable();
        List<String> names = table.getStationNames();
        double x2 = table.getTput().get(names.indexOf("P2"));
        double x3 = table.getTput().get(names.indexOf("P3"));
        assertTrue(x2 + x3 > 0, "net is deadlocked");
        return x2 / (x2 + x3);
    }

    private void assertWeightRatio(double w1, double w2) {
        assertEquals(w1 / (w1 + w2), branchProbability(w1, w2), 1e-6,
                "branch probability for weights " + w1 + ":" + w2);
    }

    @Test
    public void testBranchingFollowsWeights1to3() {
        assertWeightRatio(1.0, 3.0);
    }

    @Test
    public void testBranchingFollowsWeightsEqual() {
        assertWeightRatio(1.0, 1.0);
    }

    @Test
    public void testBranchingFollowsWeights2to1() {
        assertWeightRatio(2.0, 1.0);
    }

    @Test
    public void testBranchingFollowsWeights1to9() {
        assertWeightRatio(1.0, 9.0);
    }

    /**
     * Guard the specific regression value: gating on the latched server count
     * instead of on the marking produced 0.34375 for weights 1:3. Asserting we
     * are not back at that value stops a reintroduction from hiding behind a
     * loose tolerance.
     */
    @Test
    public void testBranchingIsNotTheLatchingArtifact() {
        double observed = branchProbability(1.0, 3.0);
        assertTrue(Math.abs(observed - 0.34375) > 1e-3,
                "branching reverted to the latching-order artifact (" + observed + ")");
    }

    // ---------------------------------------------------------------
    // 2. Non-Markovian firing time expanded to phase-type
    // ---------------------------------------------------------------

    /**
     * Open SPN forming an M/G/1/K queue: Source -> P1 -> T1 -> Sink. The firing
     * distribution of T1 is the service time, so the mean marking of P1 is the
     * M/G/1 queue length and depends on the service SCV, not only on its mean.
     */
    private static Network buildMG1K(jline.lang.processes.Distribution firingDist) {
        Network model = new Network("spn_mg1k");
        Source source = new Source(model, "source");
        Sink sink = new Sink(model, "sink");
        Place p1 = new Place(model, "P1");
        Transition t1 = new Transition(model, "T1");
        OpenClass jc = new OpenClass(model, "jobs");

        source.setArrival(jc, Exp.fitMean(1.0 / LAMBDA));
        p1.setClassCapacity(jc, CAP);
        Mode fire = t1.addMode("fire");
        t1.setDistribution(fire, firingDist);
        t1.setEnablingConditions(fire, jc, p1, 1);

        RoutingMatrix R = model.initRoutingMatrix();
        R.set(jc, jc, source, p1, 1.0);
        R.set(jc, jc, p1, t1, 1.0);
        R.set(jc, jc, t1, sink, 1.0);
        model.link(R);
        return model;
    }

    private static double solveP1(Network model) {
        SolverCTMC solver = new SolverCTMC(model, "cutoff", CAP + 4, "verbose", VerboseLevel.SILENT);
        NetworkAvgTable table = solver.getAvgTable();
        return table.getQLen().get(table.getStationNames().indexOf("P1"));
    }

    /** Mean number in an M/G/1 system: rho + rho^2 (1 + SCV) / (2 (1 - rho)). */
    private static double pollaczekKhinchine(double scv) {
        return RHO + RHO * RHO * (1.0 + scv) / (2.0 * (1.0 - RHO));
    }

    /**
     * Det(1) and Exp(1) share a mean but differ in SCV (0 against 1), so the
     * M/G/1 queue length differs: 0.75 for M/D/1 against 1.0 for M/M/1. The
     * conversion approximates Det by an Erlang of nonmkvorder phases (20 by
     * default), whose SCV is 1/20, so the reachable target is the
     * Pollaczek-Khinchine value at SCV = 0.05, i.e. 0.7625.
     */
    @Test
    public void testDeterministicFiringMatchesMD1NotMM1() {
        double qDet = solveP1(buildMG1K(new Det(1.0 / MU)));
        double qExp = solveP1(buildMG1K(Exp.fitMean(1.0 / MU)));

        assertTrue(Double.isFinite(qDet), "Det firing produced a non-finite marking");
        // Exponential firing reproduces M/M/1 (truncation at K=12, rho=0.5, is ~1e-4).
        assertEquals(pollaczekKhinchine(1.0), qExp, 2e-3 * pollaczekKhinchine(1.0),
                "Exp firing should reproduce M/M/1");
        // Det firing must land on the Erlang-20 approximation of M/D/1.
        assertEquals(pollaczekKhinchine(1.0 / 20.0), qDet, 5e-3 * pollaczekKhinchine(1.0 / 20.0),
                "Det firing should reproduce M/D/1 via Erlang-20, not M/M/1 ("
                        + pollaczekKhinchine(1.0) + ")");
        assertTrue(qDet < qExp, "lower service variability must not raise the queue");
    }
}
