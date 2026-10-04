/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.ldes;

import jline.VerboseLevel;
import jline.api.qsys.GeoGeo1Convention;
import jline.api.qsys.GeoGeo1Result;
import jline.api.qsys.Qsys_geogeo1;
import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.RoutingMatrix;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.processes.Det;
import jline.lang.processes.Geometric;
import jline.solvers.NetworkAvgTable;
import jline.solvers.ldes.handlers.Solver_ssj;
import org.junit.jupiter.api.Test;
import umontreal.ssj.rng.MRG32k3a;

import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Discrete-time validation of the LDES engine against the exact Geo/Geo/1
 * formulas in {@link Qsys_geogeo1}.
 *
 * A Source with a Geometric(a) interarrival time feeding an FCFS Queue with a
 * Geometric(s) service time is a Geo/Geo/1 queue in the late-arrival
 * delayed-access convention: both distributions are supported on {1,2,...}, so
 * every event lands on the slot lattice and no job can depart in its own arrival
 * slot.
 */
public class SolverLDESGeoGeo1Test {

    private static final int SEED = 23000;

    /** Builds the Geo/Geo/1 model with per-slot probabilities a and s. */
    private static Network geoGeo1Model(double a, double s) {
        Network model = new Network("GeoGeo1");
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");

        OpenClass jobClass = new OpenClass(model, "Class1");
        source.setArrival(jobClass, new Geometric(a));
        queue.setService(jobClass, new Geometric(s));

        RoutingMatrix P = model.initRoutingMatrix();
        P.set(jobClass, jobClass, source, queue, 1.0);
        P.set(jobClass, jobClass, queue, sink, 1.0);
        model.link(P);
        return model;
    }

    private static LDESOptions slottedOptions(int samples) {
        LDESOptions options = new LDESOptions();
        options.verbose = VerboseLevel.SILENT;
        options.seed = SEED;
        options.samples = samples;
        options.setSlotted(true);
        return options;
    }

    // ------------------------------------------------------------------
    // The simulator reproduces the closed form
    // ------------------------------------------------------------------

    /**
     * The observed agreement at this seed and sample count is within 0.4% on all
     * four metrics; the 1.5% band leaves headroom for the sampling noise without
     * losing the ability to detect a genuine regression.
     *
     * <p>This also settles the estimator question the mode raises. LDES reports a
     * time-integrated mean, whereas the discrete-time chain is defined at slot
     * boundaries. The two coincide because the state is constant across each open
     * slot and simultaneous events contribute zero area, and the agreement below
     * is the evidence: a mean queue length biased by the intra-slot event order
     * could not track the closed form this closely.
     */
    @Test
    public void slottedSimulationMatchesGeoGeo1ClosedForm() {
        double[][] grid = {{0.2, 0.5}, {0.1, 0.4}, {0.3, 0.6}};
        for (double[] pt : grid) {
            double a = pt[0];
            double s = pt[1];
            GeoGeo1Result exact = Qsys_geogeo1.qsys_geogeo1(a, s, GeoGeo1Convention.LAS_DA);

            SolverLDES solver = new SolverLDES(geoGeo1Model(a, s), slottedOptions(2000000));
            NetworkAvgTable table = solver.getAvgTable();

            // Row 0 is the Source, row 1 is the Queue.
            double qlen = table.getQLen().get(1);
            double respT = table.getRespT().get(1);
            double util = table.getUtil().get(1);
            double tput = table.getTput().get(1);

            String at = " (a=" + a + ", s=" + s + ")";
            assertRelative(exact.getMeanQueueLength(), qlen, 0.015, "mean queue length" + at);
            assertRelative(exact.getMeanSojournTime(), respT, 0.015, "mean sojourn time" + at);
            assertRelative(exact.getUtilization(), util, 0.015, "utilization" + at);
            assertRelative(exact.getThroughput(), tput, 0.015, "throughput" + at);
        }
    }

    @Test
    public void slottedSimulationSatisfiesLittlesLaw() {
        SolverLDES solver = new SolverLDES(geoGeo1Model(0.2, 0.5), slottedOptions(1000000));
        NetworkAvgTable table = solver.getAvgTable();
        double qlen = table.getQLen().get(1);
        double respT = table.getRespT().get(1);
        double tput = table.getTput().get(1);
        assertRelative(qlen, tput * respT, 0.02, "simulated E[N] = X E[T]");
    }

    @Test
    public void meanSojournNeverFallsBelowOneSlot() {
        // Service is supported on {1,2,...}, so no job can traverse the queue in
        // less than a full slot however light the load.
        SolverLDES solver = new SolverLDES(geoGeo1Model(0.02, 0.9), slottedOptions(400000));
        NetworkAvgTable table = solver.getAvgTable();
        assertTrue(table.getRespT().get(1) >= 1.0 - 1e-9,
                "sojourn must be at least one slot, got " + table.getRespT().get(1));
    }

    // ------------------------------------------------------------------
    // The mode is a semantics guarantee, not a numerical change
    // ------------------------------------------------------------------

    @Test
    public void continuousModeAgreesOnAnAlreadyLatticeModel() {
        // Geometric samples are integral either way, so switching the mode on must
        // not move the estimates: what it adds is the lattice check and the fixed
        // intra-slot event ordering.
        LDESOptions continuous = new LDESOptions();
        continuous.verbose = VerboseLevel.SILENT;
        continuous.seed = SEED;
        continuous.samples = 1000000;

        NetworkAvgTable slottedTable =
                new SolverLDES(geoGeo1Model(0.2, 0.5), slottedOptions(1000000)).getAvgTable();
        NetworkAvgTable continuousTable =
                new SolverLDES(geoGeo1Model(0.2, 0.5), continuous).getAvgTable();

        assertRelative(slottedTable.getQLen().get(1), continuousTable.getQLen().get(1),
                0.02, "queue length is mode independent");
        assertRelative(slottedTable.getRespT().get(1), continuousTable.getRespT().get(1),
                0.02, "sojourn time is mode independent");
    }

    // ------------------------------------------------------------------
    // Lattice enforcement
    // ------------------------------------------------------------------

    @Test
    public void nonLatticeServiceTimeIsRejected() {
        Network model = new Network("GeoDet");
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");
        OpenClass jobClass = new OpenClass(model, "Class1");
        source.setArrival(jobClass, new Geometric(0.2));
        queue.setService(jobClass, new Det(2.5));
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(jobClass, jobClass, source, queue, 1.0);
        P.set(jobClass, jobClass, queue, sink, 1.0);
        model.link(P);

        RuntimeException thrown = assertThrows(RuntimeException.class,
                () -> new SolverLDES(model, slottedOptions(10000)).getAvgTable());
        assertTrue(rootMessage(thrown).contains("slotted mode"),
                "a non-lattice service time must be reported, not rounded; got: "
                        + rootMessage(thrown));
    }

    @Test
    public void nonLatticeServiceTimeIsAcceptedInContinuousMode() {
        Network model = new Network("GeoDetContinuous");
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");
        OpenClass jobClass = new OpenClass(model, "Class1");
        source.setArrival(jobClass, new Geometric(0.2));
        queue.setService(jobClass, new Det(2.5));
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(jobClass, jobClass, source, queue, 1.0);
        P.set(jobClass, jobClass, queue, sink, 1.0);
        model.link(P);

        LDESOptions continuous = new LDESOptions();
        continuous.verbose = VerboseLevel.SILENT;
        continuous.seed = SEED;
        continuous.samples = 50000;
        NetworkAvgTable table = new SolverLDES(model, continuous).getAvgTable();
        assertTrue(table.getRespT().get(1) >= 2.5,
                "continuous mode must still simulate the deterministic 2.5 service");
    }

    // ------------------------------------------------------------------
    // Generator moments
    // ------------------------------------------------------------------

    @Test
    public void shiftedGeometricGenReproducesGeometricMoments() {
        for (double p : new double[]{0.1, 0.25, 0.5, 0.8}) {
            MRG32k3a stream = new MRG32k3a();
            Solver_ssj.ShiftedGeometricGen gen = new Solver_ssj.ShiftedGeometricGen(stream, p);
            int n = 2000000;
            double sum = 0.0;
            double sumSq = 0.0;
            double min = Double.POSITIVE_INFINITY;
            for (int i = 0; i < n; i++) {
                double x = gen.nextDouble();
                assertEquals(x, Math.rint(x), 0.0, "samples must be integral");
                sum += x;
                sumSq += x * x;
                if (x < min) min = x;
            }
            double mean = sum / n;
            double var = sumSq / n - mean * mean;
            double scv = var / (mean * mean);
            assertRelative(1.0 / p, mean, 0.01, "mean must be 1/p for p=" + p);
            assertRelative(1.0 - p, scv, 0.05, "SCV must be 1-p for p=" + p);
            assertEquals(1.0, min, 0.0, "support must start at one for p=" + p);
        }
    }

    @Test
    public void shiftedGeometricGenRejectsInvalidProbability() {
        MRG32k3a stream = new MRG32k3a();
        assertThrows(IllegalArgumentException.class,
                () -> new Solver_ssj.ShiftedGeometricGen(stream, 0.0));
        assertThrows(IllegalArgumentException.class,
                () -> new Solver_ssj.ShiftedGeometricGen(stream, 1.5));
    }

    @Test
    public void degenerateGeometricAlwaysCompletesInOneSlot() {
        MRG32k3a stream = new MRG32k3a();
        Solver_ssj.ShiftedGeometricGen gen = new Solver_ssj.ShiftedGeometricGen(stream, 1.0);
        for (int i = 0; i < 1000; i++) {
            assertEquals(1.0, gen.nextDouble(), 0.0, "p=1 completes in the first slot");
        }
    }

    // ------------------------------------------------------------------
    // Helpers
    // ------------------------------------------------------------------

    private static void assertRelative(double expected, double actual, double relTol, String what) {
        double tol = relTol * Math.max(Math.abs(expected), 1e-12);
        assertEquals(expected, actual, tol,
                what + ": expected " + expected + " within " + (relTol * 100) + "%, got " + actual);
    }

    private static String rootMessage(Throwable t) {
        StringBuilder sb = new StringBuilder();
        Throwable cur = t;
        while (cur != null) {
            if (cur.getMessage() != null) sb.append(cur.getMessage()).append(" | ");
            cur = cur.getCause();
        }
        return sb.toString();
    }
}
