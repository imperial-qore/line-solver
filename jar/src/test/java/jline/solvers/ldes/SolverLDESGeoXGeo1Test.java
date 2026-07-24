/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.ldes;

import jline.VerboseLevel;
import jline.api.qsys.GeoGeo1Convention;
import jline.api.qsys.GeoXGeo1Result;
import jline.api.qsys.Qsys_geoxgeo1;
import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.RoutingMatrix;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.processes.Geometric;
import jline.solvers.NetworkAvgTable;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Discrete-time validation of batch arrivals in the LDES engine against the
 * exact Geo^X/Geo/1 formulas in {@link Qsys_geoxgeo1}.
 *
 * A Source whose interarrival time is Geometric(a) and whose arrival batch is
 * Geometric(beta), feeding an FCFS Queue with Geometric(s) service, is a
 * Geo^X/Geo/1 queue in the late-arrival delayed-access convention.
 */
public class SolverLDESGeoXGeo1Test {

    private static final int SEED = 23000;

    /** Geo^X/Geo/1: batch epochs Geometric(a), batch size Geometric(beta). */
    private static Network geoXGeo1Model(double a, double beta, double s) {
        Network model = new Network("GeoXGeo1");
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");

        OpenClass jobClass = new OpenClass(model, "Class1");
        source.setArrival(jobClass, new Geometric(a));
        source.setArrivalBatch(jobClass, new Geometric(beta));
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

    @Test
    public void slottedSimulationMatchesGeoXGeo1ClosedForm() {
        // (a, beta, s); lambda = a/beta must stay below s.
        double[][] grid = {{0.1, 0.5, 0.9}, {0.2, 0.8, 0.6}, {0.15, 0.4, 0.8}};
        for (double[] pt : grid) {
            double a = pt[0];
            double beta = pt[1];
            double s = pt[2];
            GeoXGeo1Result exact =
                    Qsys_geoxgeo1.qsys_geoxgeo1(a, beta, s, GeoGeo1Convention.LAS_DA);

            NetworkAvgTable table =
                    new SolverLDES(geoXGeo1Model(a, beta, s), slottedOptions(2000000)).getAvgTable();

            double qlen = table.getQLen().get(1);
            double respT = table.getRespT().get(1);
            double util = table.getUtil().get(1);
            double tput = table.getTput().get(1);

            String at = " (a=" + a + ", beta=" + beta + ", s=" + s + ")";
            assertRelative(exact.getMeanQueueLength(), qlen, 0.03, "mean queue length" + at);
            assertRelative(exact.getMeanSojournTime(), respT, 0.03, "mean sojourn time" + at);
            assertRelative(exact.getUtilization(), util, 0.02, "utilization" + at);
            assertRelative(exact.getThroughput(), tput, 0.02, "throughput" + at);
        }
    }

    @Test
    public void throughputCountsJobsNotBatches() {
        // The distinguishing property of a batch stream: the job arrival rate is
        // a*E[X], not a. A solver that ignored the batch law would report a.
        double a = 0.1;
        double beta = 0.25;   // E[X] = 4, so lambda = 0.4
        double s = 0.8;
        NetworkAvgTable table =
                new SolverLDES(geoXGeo1Model(a, beta, s), slottedOptions(1000000)).getAvgTable();
        assertRelative(a / beta, table.getTput().get(1), 0.02, "throughput must be a*E[X]");
        assertTrue(table.getTput().get(1) > 2.0 * a,
                "throughput must exceed the batch epoch rate");
    }

    @Test
    public void degenerateBatchReducesToGeoGeo1() {
        // Geometric(1) always yields a batch of exactly one job, so the run must
        // agree with the single-arrival model.
        double a = 0.2;
        double s = 0.5;
        NetworkAvgTable batched =
                new SolverLDES(geoXGeo1Model(a, 1.0, s), slottedOptions(1000000)).getAvgTable();

        Network single = new Network("GeoGeo1");
        Source source = new Source(single, "Source");
        Queue queue = new Queue(single, "Queue", SchedStrategy.FCFS);
        Sink sink = new Sink(single, "Sink");
        OpenClass jobClass = new OpenClass(single, "Class1");
        source.setArrival(jobClass, new Geometric(a));
        queue.setService(jobClass, new Geometric(s));
        RoutingMatrix P = single.initRoutingMatrix();
        P.set(jobClass, jobClass, source, queue, 1.0);
        P.set(jobClass, jobClass, queue, sink, 1.0);
        single.link(P);
        NetworkAvgTable plain =
                new SolverLDES(single, slottedOptions(1000000)).getAvgTable();

        assertRelative(plain.getQLen().get(1), batched.getQLen().get(1), 0.02,
                "unit batch must match the single-arrival model");
        assertRelative(plain.getRespT().get(1), batched.getRespT().get(1), 0.02,
                "unit batch must match the single-arrival model");
    }

    @Test
    public void batchingRaisesCongestionAtEqualLoad() {
        // Hold lambda = a/beta fixed and enlarge the batches: the simulated
        // sojourn must rise, as the closed form predicts.
        double s = 0.8;
        double lambda = 0.4;
        double previous = 0.0;
        for (double beta : new double[]{1.0, 0.5, 0.25}) {
            double a = lambda * beta;
            NetworkAvgTable table =
                    new SolverLDES(geoXGeo1Model(a, beta, s), slottedOptions(1000000)).getAvgTable();
            double respT = table.getRespT().get(1);
            assertTrue(respT > previous,
                    "larger batches must not reduce the sojourn, beta=" + beta + " got " + respT);
            previous = respT;
        }
    }

    // ------------------------------------------------------------------
    // Model validation
    // ------------------------------------------------------------------

    @Test
    public void batchLawSupportedBelowOneIsRejected() {
        Network model = new Network("BadBatch");
        Source source = new Source(model, "Source");
        new Queue(model, "Queue", SchedStrategy.FCFS);
        new Sink(model, "Sink");
        OpenClass jobClass = new OpenClass(model, "Class1");
        source.setArrival(jobClass, new Geometric(0.2));
        // Bernoulli(0.5) has mean 0.5, so some "batches" would carry no job.
        assertThrows(IllegalArgumentException.class,
                () -> source.setArrivalBatch(jobClass, new jline.lang.processes.Bernoulli(0.5)),
                "a batch must carry at least one job");
    }

    @Test
    public void nullBatchRestoresSingleArrivals() {
        Network model = new Network("ResetBatch");
        Source source = new Source(model, "Source");
        OpenClass jobClass = new OpenClass(model, "Class1");
        source.setArrivalBatch(jobClass, new Geometric(0.5));
        assertEquals(0.5, ((Geometric) source.getArrivalBatch(jobClass)).getRate(), 1e-12);
        source.setArrivalBatch(jobClass, null);
        assertEquals(null, source.getArrivalBatch(jobClass),
                "a null batch law must restore single arrivals");
    }

    // ------------------------------------------------------------------
    // JSON round-trip (the channel MATLAB and Python use)
    // ------------------------------------------------------------------

    @Test
    public void arrivalBatchSurvivesTheModelJson() throws Exception {
        // MATLAB and Python reach the engine only through model.json, so a batch
        // law that does not serialize is a batch law those codebases cannot use.
        Network model = geoXGeo1Model(0.1, 0.5, 0.9);
        java.io.File f = java.io.File.createTempFile("geoxbatch", ".json");
        f.deleteOnExit();
        jline.io.LineModelIO.save(model, f.getAbsolutePath());

        String json = new String(java.nio.file.Files.readAllBytes(f.toPath()), "UTF-8");
        assertTrue(json.contains("arrivalBatch"), "arrivalBatch must appear in the saved JSON");

        Network reloaded = (Network) jline.io.LineModelIO.load(f.getAbsolutePath());
        Source reloadedSource = (Source) reloaded.getNodeByName("Source");
        jline.lang.JobClass reloadedClass = reloaded.getClassByName("Class1");
        jline.lang.processes.DiscreteDistribution batch =
                reloadedSource.getArrivalBatch(reloadedClass);
        assertTrue(batch != null, "batch law must survive the round trip");
        assertEquals(2.0, batch.getMean(), 1e-12, "batch mean E[X] = 1/beta");
        assertEquals(10.0, reloadedSource.getArrivalDistribution(reloadedClass).getMean(),
                1e-12, "interarrival mean must be untouched by the batch round trip");

        // sn.arrivalbatch is what the engine reads, so check it too.
        assertTrue(reloaded.getStruct(true).arrivalbatch.get(0) != null,
                "sn.arrivalbatch must be populated after a reload");
    }

    @Test
    public void reloadedBatchModelReproducesTheClosedForm() throws Exception {
        Network model = geoXGeo1Model(0.1, 0.5, 0.9);
        java.io.File f = java.io.File.createTempFile("geoxsolve", ".json");
        f.deleteOnExit();
        jline.io.LineModelIO.save(model, f.getAbsolutePath());
        Network reloaded = (Network) jline.io.LineModelIO.load(f.getAbsolutePath());

        NetworkAvgTable table =
                new SolverLDES(reloaded, slottedOptions(2000000)).getAvgTable();
        GeoXGeo1Result exact =
                Qsys_geoxgeo1.qsys_geoxgeo1(0.1, 0.5, 0.9, GeoGeo1Convention.LAS_DA);
        assertRelative(exact.getMeanQueueLength(), table.getQLen().get(1), 0.03,
                "reloaded model mean queue length");
        assertRelative(exact.getThroughput(), table.getTput().get(1), 0.02,
                "reloaded model throughput must still be a*E[X]");
    }

    // ------------------------------------------------------------------
    // Helpers
    // ------------------------------------------------------------------

    private static void assertRelative(double expected, double actual, double relTol, String what) {
        double tol = relTol * Math.max(Math.abs(expected), 1e-12);
        assertEquals(expected, actual, tol,
                what + ": expected " + expected + " within " + (relTol * 100) + "%, got " + actual);
    }
}
