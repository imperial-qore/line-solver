package jline.solvers.ssa;

import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.constant.SchedStrategy;
import jline.lang.OpenClass;
import jline.lang.OpenSignal;
import jline.lang.RoutingMatrix;
import jline.lang.constant.RemovalPolicy;
import jline.lang.constant.SignalType;
import jline.lang.processes.Immediate;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.nodes.Queue;
import jline.lang.processes.Exp;
import jline.solvers.ctmc.SolverCTMC;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Validates the NRM propensities of the processor-sharing family against the
 * exact CTMC solution.
 * <p>
 * The DPS/GPS rate laws read the per-class weights from sn.schedparam, and the
 * PSPRIO/DPSPRIO/GPSPRIO ones additionally gate on the priority group. A
 * previous implementation collapsed DPS and GPS onto plain PS and dropped the
 * weights entirely, which these tests would not detect unless the weights are
 * strongly asymmetric and the populations exceed one job per class -- hence the
 * fixtures below.
 * </p>
 */
public class SolverSSANrmPsFamilyTest {

    private static final int SAMPLES = 400000;
    private static final int SEED = 23000;
    /** Relative tolerance on the simulated means. */
    private static final double RTOL = 0.03;
    /**
     * Absolute floor, in jobs. A priority-starved class sits at its reference
     * station essentially never (mean queue length ~0.012 in these fixtures), so
     * a purely relative criterion would measure simulation noise on a near-zero
     * quantity rather than the correctness of the rate law.
     */
    private static final double ATOL = 0.01;

    /**
     * Closed 2-class model whose queue runs the given policy, with asymmetric
     * service weights (3:1) so that DPS and GPS cannot coincide with PS.
     */
    private static Network twoClassModel(SchedStrategy sched, int prio1, int prio2) {
        Network model = new Network("nrm_" + sched);
        Delay delay = new Delay(model, "Think");
        Queue queue = new Queue(model, "Q1", sched);
        ClosedClass c1 = new ClosedClass(model, "C1", 2, delay, prio1);
        ClosedClass c2 = new ClosedClass(model, "C2", 2, delay, prio2);
        delay.setService(c1, new Exp(1.0));
        delay.setService(c2, new Exp(2.0));
        queue.setService(c1, new Exp(3.0), 3.0);
        queue.setService(c2, new Exp(5.0), 1.0);
        model.link(model.serialRouting(delay, queue));
        return model;
    }

    /**
     * Three-class model with two classes sharing the urgent priority group, so
     * the weights are live inside the group while the priority gate is live
     * across groups. With one job per class DPS and GPS coincide (n_r equals
     * the activity indicator), so each urgent class carries two jobs.
     */
    private static Network prioWeightModel(SchedStrategy sched) {
        Network model = new Network("prio_" + sched);
        Delay delay = new Delay(model, "Think");
        Queue queue = new Queue(model, "Q1", sched);
        ClosedClass c1 = new ClosedClass(model, "C1", 2, delay, 1);
        ClosedClass c2 = new ClosedClass(model, "C2", 2, delay, 1);
        ClosedClass c3 = new ClosedClass(model, "C3", 1, delay, 2);
        delay.setService(c1, new Exp(1.0));
        delay.setService(c2, new Exp(1.0));
        delay.setService(c3, new Exp(1.0));
        queue.setService(c1, new Exp(1.0), 9.0);
        queue.setService(c2, new Exp(1.0), 1.0);
        queue.setService(c3, new Exp(1.0), 1.0);
        model.link(model.serialRouting(delay, queue));
        return model;
    }

    /**
     * Worst deviation under a mixed criterion |e-s| / (ATOL + |e|), so that
     * entries of negligible magnitude are judged on an absolute scale and
     * larger ones on a relative scale.
     */
    private static double maxRelErr(Matrix exact, Matrix sim) {
        double worst = 0.0;
        for (int i = 0; i < exact.getNumRows(); i++) {
            for (int j = 0; j < exact.getNumCols(); j++) {
                double e = exact.get(i, j);
                double s = sim.get(i, j);
                worst = Math.max(worst, Math.abs(e - s) / (ATOL + Math.abs(e)));
            }
        }
        return worst;
    }

    private static Matrix nrmQLen(Network model) {
        SolverSSA solver = new SolverSSA(model);
        solver.options.method = "nrm";
        solver.options.samples = SAMPLES;
        solver.options.seed = SEED;
        return solver.getAvgQLen();
    }

    private static void assertMatchesCtmc(Network sim, Network exact, String label) {
        Matrix qExact = new SolverCTMC(exact).getAvgQLen();
        Matrix qSim = nrmQLen(sim);
        double err = maxRelErr(qExact, qSim);
        assertTrue(err < RTOL, label + ": NRM deviates from CTMC by " + err
                + " (exact=" + qExact + ", nrm=" + qSim + ")");
    }

    @Test
    public void dpsMatchesCtmc() {
        assertMatchesCtmc(twoClassModel(SchedStrategy.DPS, 0, 0),
                twoClassModel(SchedStrategy.DPS, 0, 0), "DPS");
    }

    @Test
    public void gpsMatchesCtmc() {
        assertMatchesCtmc(twoClassModel(SchedStrategy.GPS, 0, 0),
                twoClassModel(SchedStrategy.GPS, 0, 0), "GPS");
    }

    @Test
    public void lpsMatchesCtmc() {
        assertMatchesCtmc(twoClassModel(SchedStrategy.LPS, 0, 0),
                twoClassModel(SchedStrategy.LPS, 0, 0), "LPS");
    }

    @Test
    public void psprioMatchesCtmc() {
        assertMatchesCtmc(prioWeightModel(SchedStrategy.PSPRIO),
                prioWeightModel(SchedStrategy.PSPRIO), "PSPRIO");
    }

    @Test
    public void dpsprioMatchesCtmc() {
        assertMatchesCtmc(prioWeightModel(SchedStrategy.DPSPRIO),
                prioWeightModel(SchedStrategy.DPSPRIO), "DPSPRIO");
    }

    @Test
    public void gpsprioMatchesCtmc() {
        assertMatchesCtmc(prioWeightModel(SchedStrategy.GPSPRIO),
                prioWeightModel(SchedStrategy.GPSPRIO), "GPSPRIO");
    }

    /**
     * The weights must actually reach the rate law: DPS with a 3:1 weight split
     * must not produce the plain-PS queue lengths.
     */
    @Test
    public void dpsWeightsAreNotIgnored() {
        Matrix qDps = new SolverCTMC(twoClassModel(SchedStrategy.DPS, 0, 0)).getAvgQLen();
        Matrix qPs = new SolverCTMC(twoClassModel(SchedStrategy.PS, 0, 0)).getAvgQLen();
        assertTrue(maxRelErr(qPs, qDps) > 0.05,
                "fixture is degenerate: weighted DPS must differ from PS");
        Matrix qNrm = nrmQLen(twoClassModel(SchedStrategy.DPS, 0, 0));
        assertTrue(maxRelErr(qDps, qNrm) < RTOL, "NRM DPS does not track weighted DPS");
        assertTrue(maxRelErr(qPs, qNrm) > 0.05, "NRM DPS collapsed onto unweighted PS");
    }

    /**
     * GPS splits the weights across active classes rather than across jobs, so
     * it must separate from DPS once a class can hold several jobs.
     */
    @Test
    public void gpsSeparatesFromDps() {
        Matrix qGps = new SolverCTMC(twoClassModel(SchedStrategy.GPS, 0, 0)).getAvgQLen();
        Matrix qDps = new SolverCTMC(twoClassModel(SchedStrategy.DPS, 0, 0)).getAvgQLen();
        assertTrue(maxRelErr(qDps, qGps) > 0.005,
                "fixture is degenerate: GPS must differ from DPS");
        Matrix qNrm = nrmQLen(twoClassModel(SchedStrategy.GPS, 0, 0));
        assertTrue(maxRelErr(qGps, qNrm) < RTOL, "NRM GPS does not track GPS");
    }

    /**
     * Buffered non-preemptive model: distinct means (so SEPT/LEPT have a strict
     * order) and distinct priorities (so HOL has one). The buffered policies all
     * share the FCFS rate law and differ only in which waiting job is promoted.
     */
    private static Network bufferedModel(SchedStrategy sched, int nclasses) {
        double[] svc = {1.0, 4.0, 2.0};   // means 1.000, 0.250, 0.500
        Network model = new Network("buf_" + sched);
        Delay delay = new Delay(model, "Think");
        Queue queue = new Queue(model, "Q1", sched);
        for (int k = 0; k < nclasses; k++) {
            ClosedClass c = new ClosedClass(model, "C" + (k + 1), 2, delay, k + 1);
            delay.setService(c, new Exp(1.0));
            queue.setService(c, new Exp(svc[k]));
        }
        model.link(model.serialRouting(delay, queue));
        return model;
    }

    private static void assertBufferedMatchesCtmc(SchedStrategy sched, int nclasses) {
        String label = sched + "/" + nclasses + "class";
        assertMatchesCtmc(bufferedModel(sched, nclasses), bufferedModel(sched, nclasses), label);
    }

    @Test
    public void siroMatchesCtmc() {
        assertBufferedMatchesCtmc(SchedStrategy.SIRO, 2);
        assertBufferedMatchesCtmc(SchedStrategy.SIRO, 3);
    }

    @Test
    public void holMatchesCtmc() {
        assertBufferedMatchesCtmc(SchedStrategy.HOL, 2);
        assertBufferedMatchesCtmc(SchedStrategy.HOL, 3);
    }

    @Test
    public void septMatchesCtmc() {
        assertBufferedMatchesCtmc(SchedStrategy.SEPT, 2);
        assertBufferedMatchesCtmc(SchedStrategy.SEPT, 3);
    }

    @Test
    public void leptMatchesCtmc() {
        assertBufferedMatchesCtmc(SchedStrategy.LEPT, 2);
        assertBufferedMatchesCtmc(SchedStrategy.LEPT, 3);
    }

    @Test
    public void fcfsAndLcfsStillMatchCtmc() {
        assertBufferedMatchesCtmc(SchedStrategy.FCFS, 3);
        assertBufferedMatchesCtmc(SchedStrategy.LCFS, 3);
    }

    /**
     * SEPT must favour the shortest-mean class, so that class must hold the
     * smallest queue. Guards against a promotion order that merely happens to
     * self-consistently name whichever class it visited first.
     */
    @Test
    public void septFavoursShortestMeanClass() {
        Matrix q = new SolverCTMC(bufferedModel(SchedStrategy.SEPT, 3)).getAvgQLen();
        // class 2 has the shortest mean (0.25); class 1 the longest (1.0)
        assertTrue(q.get(1, 1) < q.get(1, 0), "SEPT: shortest-mean class must not have the largest queue");
        assertTrue(q.get(1, 1) < q.get(1, 2), "SEPT: shortest-mean class must have the smallest queue");
        Matrix qn = nrmQLen(bufferedModel(SchedStrategy.SEPT, 3));
        assertTrue(qn.get(1, 1) < qn.get(1, 0) && qn.get(1, 1) < qn.get(1, 2),
                "SEPT (NRM): shortest-mean class must have the smallest queue");
    }

    /**
     * Plain LCFS is not priority-aware: SchedStrategy.LCFSPRIO is. Setting
     * distinct priorities on an LCFS station must not change any metric.
     */
    @Test
    public void lcfsIgnoresPriorities() {
        Network flat = bufferedModel(SchedStrategy.LCFS, 3);
        Network prio = bufferedModel(SchedStrategy.LCFS, 3);
        // bufferedModel already assigns distinct priorities; build a flat twin
        Network flatTwin = new Network("lcfs_flat");
        Delay d = new Delay(flatTwin, "Think");
        Queue q = new Queue(flatTwin, "Q1", SchedStrategy.LCFS);
        double[] svc = {1.0, 4.0, 2.0};
        for (int k = 0; k < 3; k++) {
            ClosedClass c = new ClosedClass(flatTwin, "C" + (k + 1), 2, d, 0);
            d.setService(c, new Exp(1.0));
            q.setService(c, new Exp(svc[k]));
        }
        flatTwin.link(flatTwin.serialRouting(d, q));

        Matrix qFlat = new SolverCTMC(flatTwin).getAvgQLen();
        Matrix qPrio = new SolverCTMC(prio).getAvgQLen();
        assertTrue(maxRelErr(qFlat, qPrio) < 1e-9,
                "LCFS must ignore class priorities (LCFSPRIO is the priority-aware variant): flat="
                        + qFlat + " prio=" + qPrio);
    }

    /**
     * Two PS queues inside a DROP-rule finite capacity region with a global job
     * cap. The region constrains an aggregate of the member stations' per-class
     * populations, which is linear in the NRM state vector, so admission is a
     * gate on the routing draw; the DROP rule censors the refused transition.
     */
    private static Network fcrModel(int cap) {
        Network model = new Network("fcr");
        Delay delay = new Delay(model, "Think");
        Queue q1 = new Queue(model, "Q1", SchedStrategy.PS);
        Queue q2 = new Queue(model, "Q2", SchedStrategy.PS);
        ClosedClass c1 = new ClosedClass(model, "C1", 3, delay, 0);
        ClosedClass c2 = new ClosedClass(model, "C2", 2, delay, 0);
        delay.setService(c1, new Exp(1.0));
        delay.setService(c2, new Exp(1.0));
        q1.setService(c1, new Exp(2.0));
        q1.setService(c2, new Exp(3.0));
        q2.setService(c1, new Exp(2.5));
        q2.setService(c2, new Exp(1.5));
        model.link(model.serialRouting(delay, q1, q2));
        if (cap > 0) {
            java.util.List<jline.lang.nodes.Node> nodes = new java.util.ArrayList<jline.lang.nodes.Node>();
            nodes.add(q1);
            nodes.add(q2);
            model.addRegion(nodes).setGlobalMaxJobs(cap);
        }
        return model;
    }

    private static double regionOccupancy(Matrix q) {
        double occ = 0.0;
        for (int j = 0; j < q.getNumCols(); j++) {
            occ += q.get(1, j) + q.get(2, j);   // the two queues inside the region
        }
        return occ;
    }

    @Test
    public void dropRuleFcrMatchesCtmc() {
        assertMatchesCtmc(fcrModel(2), fcrModel(2), "FCR cap=2");
        assertMatchesCtmc(fcrModel(3), fcrModel(3), "FCR cap=3");
    }

    /**
     * The cap must actually bind, and it must be respected: otherwise the
     * fixture would pass even with the gate removed entirely.
     */
    @Test
    public void dropRuleFcrCapBinds() {
        double occFree = regionOccupancy(new SolverCTMC(fcrModel(0)).getAvgQLen());
        double occCap2 = regionOccupancy(new SolverCTMC(fcrModel(2)).getAvgQLen());
        assertTrue(occFree > occCap2 + 0.5,
                "fixture is degenerate: the region cap must bind (free=" + occFree + ", cap2=" + occCap2 + ")");
        assertTrue(occCap2 <= 2.0, "region occupancy must respect the cap: " + occCap2);

        double occNrm = regionOccupancy(nrmQLen(fcrModel(2)));
        assertTrue(occNrm <= 2.0, "NRM must never exceed the region cap: " + occNrm);
        assertTrue(Math.abs(occNrm - occCap2) < 0.1,
                "NRM region occupancy must track the exact one (nrm=" + occNrm + ", ctmc=" + occCap2 + ")");
    }

    // ---- balking and reneging -------------------------------------------

    /**
     * M/M/1 with QUEUE_LENGTH balking. A balked job is lost, so the departure
     * rate is untouched and only the arrival outcome differs. Reference values
     * come from the MATLAB CTMC at cutoff=30 on the identical model; a cutoff
     * that low would truncate the unbalked queue, hence the exact M/M/1 value
     * is used for the zero-balking case instead.
     */
    private static Network balkModel(double balkProb) {
        Network model = new Network("balk");
        Source source = new Source(model, "Src");
        Queue queue = new Queue(model, "Q1", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Snk");
        OpenClass oc = new OpenClass(model, "C1", 0);
        source.setArrival(oc, new Exp(1.0));
        queue.setService(oc, new Exp(1.5));
        queue.setNumberOfServers(1);
        model.link(model.serialRouting(source, queue, sink));
        if (balkProb > 0) {
            java.util.List<jline.lang.constant.BalkingThreshold> th =
                    new java.util.ArrayList<jline.lang.constant.BalkingThreshold>();
            th.add(new jline.lang.constant.BalkingThreshold(2, 20, balkProb));
            queue.setBalking(oc, jline.lang.constant.BalkingStrategy.QUEUE_LENGTH, th);
        }
        return model;
    }

    private static Network renegeModel(double theta) {
        Network model = new Network("renege");
        Source source = new Source(model, "Src");
        Queue queue = new Queue(model, "Q1", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Snk");
        OpenClass oc = new OpenClass(model, "C1", 0);
        source.setArrival(oc, new Exp(1.0));
        queue.setService(oc, new Exp(1.2));
        queue.setNumberOfServers(1);
        model.link(model.serialRouting(source, queue, sink));
        if (theta > 0) {
            queue.setPatience(oc, new Exp(theta));   // memoryless patience
        }
        return model;
    }

    /**
     * Balking must both track the exact solution and actually bind: raising the
     * balk probability must shorten the queue monotonically. Without binding,
     * a no-op implementation would pass.
     */
    @Test
    public void balkingMatchesExactAndBinds() {
        // MATLAB CTMC (cutoff=30) on the identical model
        double q0 = nrmQLen(balkModel(0.0)).get(1, 0);
        double q5 = nrmQLen(balkModel(0.5)).get(1, 0);
        double q9 = nrmQLen(balkModel(0.9)).get(1, 0);
        assertTrue(Math.abs(q0 - 2.0) < 0.15, "no balking must recover M/M/1 rho/(1-rho)=2: " + q0);
        assertTrue(Math.abs(q5 - 1.0000) < 0.08, "balkProb=0.5 must match the exact 1.0000: " + q5);
        assertTrue(Math.abs(q9 - 0.7714) < 0.06, "balkProb=0.9 must match the exact 0.7714: " + q9);
        assertTrue(q0 > q5 + 0.3 && q5 > q9 + 0.1,
                "balking must bind monotonically: " + q0 + " > " + q5 + " > " + q9);
    }

    /**
     * Reneging must track the exact solution and bind: a faster abandonment
     * rate must shorten the queue.
     */
    @Test
    public void renegingMatchesExactAndBinds() {
        // MATLAB CTMC (cutoff=30) on the identical model
        double qA = nrmQLen(renegeModel(0.5)).get(1, 0);
        double qB = nrmQLen(renegeModel(2.0)).get(1, 0);
        assertTrue(Math.abs(qA - 1.1256) < 0.08, "theta=0.5 must match the exact 1.1256: " + qA);
        assertTrue(Math.abs(qB - 0.7141) < 0.06, "theta=2.0 must match the exact 0.7141: " + qB);
        assertTrue(qA > qB + 0.2, "faster abandonment must shorten the queue: " + qA + " > " + qB);
    }

    // ---- G-network signals ----------------------------------------------

    /**
     * M/M/1 with a targeted NEGATIVE signal. A signal never joins the station:
     * it removes a job there and is annihilated.
     * <p>
     * The reference values are the MATLAB CTMC (cutoff=12) on identical models,
     * NOT the JAR's own CTMC: the JAR's AfterEventStation signal path is still
     * the older single-removal semantics (no batch distribution, no removal
     * policy), so it is not an equivalent of the MATLAB reference this NRM path
     * was ported from. Cross-solver alignment of the JAR CTMC is open work.
     * </p>
     */
    private static Network signalModel(SchedStrategy sched, double sigRate, SignalType type) {
        Network model = new Network("sig");
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", sched);
        Sink sink = new Sink(model, "Sink");
        OpenClass c1 = new OpenClass(model, "Class1", 0);
        source.setArrival(c1, new Exp(2));
        queue.setService(c1, new Exp(5));
        RoutingMatrix P = model.initRoutingMatrix();
        if (sigRate > 0) {
            OpenSignal s1 = new OpenSignal(model, "Signal1", type);
            s1.forJobClass(c1);
            s1.setRemovalPolicy(RemovalPolicy.RANDOM);
            source.setArrival(s1, new Exp(sigRate));
            queue.setService(s1, new Immediate());
            P.set(c1, c1, source, queue, 1.0);
            P.set(c1, c1, queue, sink, 1.0);
            P.set(s1, s1, source, queue, 1.0);
            P.set(s1, s1, queue, sink, 1.0);
        } else {
            P.set(c1, c1, source, queue, 1.0);
            P.set(c1, c1, queue, sink, 1.0);
        }
        model.link(P);
        return model;
    }

    /**
     * Signals must track the exact solution and bind: a faster signal must
     * shorten the queue. Without the binding check a no-op signal would pass.
     */
    @Test
    public void negativeSignalMatchesExactAndBinds() {
        for (SchedStrategy sched : new SchedStrategy[]{SchedStrategy.PS, SchedStrategy.FCFS}) {
            double q0 = nrmQLen(signalModel(sched, 0.0, SignalType.NEGATIVE)).get(1, 0);
            double q5 = nrmQLen(signalModel(sched, 0.5, SignalType.NEGATIVE)).get(1, 0);
            double q3 = nrmQLen(signalModel(sched, 3.0, SignalType.NEGATIVE)).get(1, 0);
            assertTrue(Math.abs(q0 - 0.6666) < 0.05, sched + " no signal must match 0.6666: " + q0);
            assertTrue(Math.abs(q5 - 0.5714) < 0.05, sched + " rate=0.5 must match 0.5714: " + q5);
            assertTrue(Math.abs(q3 - 0.3333) < 0.04, sched + " rate=3.0 must match 0.3333: " + q3);
            assertTrue(q0 > q5 + 0.05 && q5 > q3 + 0.15,
                    sched + ": a faster signal must shorten the queue: " + q0 + " > " + q5 + " > " + q3);
        }
    }

    /** A catastrophe empties the station, so it must bite harder than a
     *  single-job negative signal at the same rate. */
    @Test
    public void catastropheSignalMatchesExact() {
        double q = nrmQLen(signalModel(SchedStrategy.PS, 1.0, SignalType.CATASTROPHE)).get(1, 0);
        assertTrue(Math.abs(q - 0.4495) < 0.05, "catastrophe must match the exact 0.4495: " + q);
        double qNeg = nrmQLen(signalModel(SchedStrategy.PS, 1.0, SignalType.NEGATIVE)).get(1, 0);
        assertTrue(q < qNeg, "a catastrophe must empty more than a negative signal: " + q + " < " + qNeg);
    }

    /**
     * Priority-awareness is a property of the DECLARED policy, never of the
     * data. Plain LCFSPR must therefore give identical results with flat and
     * with distinct class priorities; LCFSPRPRIO is the priority-aware variant.
     * <p>
     * The second half is the guard that matters: LCFSPRPRIO with those same
     * priorities must still differ from plain LCFSPR. Without it, the first
     * assertion could be satisfied by a priority variant that had silently
     * stopped honouring priorities altogether.
     * </p>
     */
    @Test
    public void lcfsprIgnoresPrioritiesButLcfsprprioDoesNot() {
        Matrix flat = new SolverCTMC(preemptModel(SchedStrategy.LCFSPR, 0, 0)).getAvgQLen();
        Matrix prio = new SolverCTMC(preemptModel(SchedStrategy.LCFSPR, 1, 2)).getAvgQLen();
        assertTrue(maxRelErr(flat, prio) < 1e-9,
                "plain LCFSPR must ignore class priorities (LCFSPRPRIO is the priority-aware"
                        + " variant): flat=" + flat + " prio=" + prio);

        Matrix prioVar = new SolverCTMC(preemptModel(SchedStrategy.LCFSPRPRIO, 1, 2)).getAvgQLen();
        assertTrue(maxRelErr(flat, prioVar) > 0.05,
                "LCFSPRPRIO must honour the priorities plain LCFSPR ignores: " + prioVar);
    }

    private static Network preemptModel(SchedStrategy sched, int prio1, int prio2) {
        Network model = new Network("prio_" + sched);
        Delay delay = new Delay(model, "Think");
        Queue queue = new Queue(model, "Q1", sched);
        ClosedClass c1 = new ClosedClass(model, "C1", 2, delay, prio1);
        ClosedClass c2 = new ClosedClass(model, "C2", 2, delay, prio2);
        delay.setService(c1, new Exp(1.0));
        delay.setService(c2, new Exp(2.0));
        queue.setService(c1, new Exp(3.0));
        queue.setService(c2, new Exp(0.6));
        model.link(model.serialRouting(delay, queue));
        return model;
    }

    // ---- tier 2: LCFSPR, round-robin, retrial ---------------------------
    // Reference values are the MATLAB CTMC (exact) on identical models. Each is
    // paired with a binding control: every one of these features has a fixture
    // in which a no-op implementation still passes.

    private static Network preemptModel(SchedStrategy sched) {
        Network model = new Network("pr");
        Delay delay = new Delay(model, "Think");
        Queue queue = new Queue(model, "Q1", sched);
        ClosedClass c1 = new ClosedClass(model, "C1", 2, delay, 0);
        ClosedClass c2 = new ClosedClass(model, "C2", 2, delay, 0);
        delay.setService(c1, new Exp(1.0));
        delay.setService(c2, new Exp(2.0));
        queue.setService(c1, new Exp(3.0));
        queue.setService(c2, new Exp(0.6));
        model.link(model.serialRouting(delay, queue));
        return model;
    }

    /** LCFSPR: the newest job is always the one in service, so with
     *  class-dependent rates it must separate from non-preemptive LCFS. */
    @Test
    public void lcfsprMatchesExactAndPreemptionBinds() {
        Matrix q = nrmQLen(preemptModel(SchedStrategy.LCFSPR));
        assertTrue(Math.abs(q.get(1, 0) - 1.0373) < 0.05, "LCFSPR c1 must match 1.0373: " + q.get(1, 0));
        assertTrue(Math.abs(q.get(1, 1) - 1.7988) < 0.06, "LCFSPR c2 must match 1.7988: " + q.get(1, 1));
        Matrix ql = nrmQLen(preemptModel(SchedStrategy.LCFS));
        assertTrue(ql.get(1, 0) - q.get(1, 0) > 0.3,
                "preemption must bind: LCFS " + ql.get(1, 0) + " vs LCFSPR " + q.get(1, 0));
    }

    private static Network rrModel(jline.lang.constant.RoutingStrategy strategy, int w) {
        Network model = new Network("rr");
        Delay delay = new Delay(model, "Think");
        jline.lang.nodes.Router router = new jline.lang.nodes.Router(model, "Rtr");
        Queue q1 = new Queue(model, "Q1", SchedStrategy.PS);
        Queue q2 = new Queue(model, "Q2", SchedStrategy.PS);
        ClosedClass c1 = new ClosedClass(model, "C1", 3, delay, 0);
        delay.setService(c1, new Exp(1.0));
        q1.setService(c1, new Exp(2.0));
        q2.setService(c1, new Exp(2.0));
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(c1, c1, delay, router, 1.0);
        P.set(c1, c1, router, q1, 1.0);
        P.set(c1, c1, router, q2, 1.0);
        P.set(c1, c1, q1, delay, 1.0);
        P.set(c1, c1, q2, delay, 1.0);
        model.link(P);
        if (strategy == jline.lang.constant.RoutingStrategy.WRROBIN) {
            router.setRouting(c1, strategy, q1, w);
            router.setRouting(c1, strategy, q2, 1);
        } else {
            router.setRouting(c1, strategy);
        }
        return model;
    }

    /** Deterministic rotation balances better than RAND, so RROBIN must
     *  separate from it; WRROBIN at 1:1 must collapse back onto RROBIN. */
    @Test
    public void roundRobinMatchesExactAndBinds() {
        Matrix qr = nrmQLen(rrModel(jline.lang.constant.RoutingStrategy.RROBIN, 0));
        assertTrue(Math.abs(qr.get(1, 0) - 0.5583) < 0.04, "RROBIN q1 must match 0.5583: " + qr.get(1, 0));
        Matrix qa = nrmQLen(rrModel(jline.lang.constant.RoutingStrategy.RAND, 0));
        assertTrue(qa.get(1, 0) - qr.get(1, 0) > 0.02,
                "RROBIN must balance better than RAND: " + qa.get(1, 0) + " vs " + qr.get(1, 0));

        Matrix qw = nrmQLen(rrModel(jline.lang.constant.RoutingStrategy.WRROBIN, 3));
        assertTrue(Math.abs(qw.get(1, 0) - 1.0107) < 0.06, "WRROBIN w=3:1 q1 must match 1.0107: " + qw.get(1, 0));
        assertTrue(Math.abs(qw.get(2, 0) - 0.2236) < 0.04, "WRROBIN w=3:1 q2 must match 0.2236: " + qw.get(2, 0));
        Matrix qw11 = nrmQLen(rrModel(jline.lang.constant.RoutingStrategy.WRROBIN, 1));
        assertTrue(Math.abs(qw11.get(1, 0) - qr.get(1, 0)) < 0.04,
                "WRROBIN at 1:1 must collapse onto RROBIN: " + qw11.get(1, 0) + " vs " + qr.get(1, 0));
    }

    private static Network retrialModel(double mu) {
        Network model = new Network("rt");
        Source source = new Source(model, "Src");
        Queue queue = new Queue(model, "Q1", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Snk");
        OpenClass oc = new OpenClass(model, "C1", 0);
        source.setArrival(oc, new Exp(1.0));
        queue.setService(oc, new Exp(2.0));
        queue.setNumberOfServers(1);
        model.link(model.serialRouting(source, queue, sink));
        queue.setRetrial(oc, new Exp(mu));
        return model;
    }

    /** A slower retrial rate keeps jobs in orbit longer, so the station
     *  population must rise. A retry that never fires would explode instead. */
    @Test
    public void retrialMatchesExactAndBinds() {
        double fast = nrmQLen(retrialModel(5.0)).get(1, 0);
        double slow = nrmQLen(retrialModel(0.5)).get(1, 0);
        assertTrue(Math.abs(fast - 1.1974) < 0.08, "retrial mu=5.0 must match 1.1974: " + fast);
        assertTrue(Math.abs(slow - 2.9546) < 0.20, "retrial mu=0.5 must match 2.9546: " + slow);
        assertTrue(slow > fast + 1.0, "a slower retrial must lengthen the orbit: " + fast + " -> " + slow);
    }
}
