/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.mam;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.List;

import org.junit.jupiter.api.Test;

import jline.VerboseLevel;
import jline.api.dqsys.Dqsys_geogeo1;
import jline.api.dqsys.GeoGeo1Convention;
import jline.api.dqsys.GeoGeo1Result;
import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.RoutingMatrix;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.processes.DMAP;
import jline.lang.processes.Det;
import jline.lang.processes.Geometric;
import jline.solvers.NetworkAvgTable;
import jline.util.matrix.Matrix;

/**
 * SolverMAM on a discrete (slotted) time scale.
 *
 * <p>The model is recognized as discrete-time from its distributions alone and
 * solved by the Q-MAM discrete-time queues under the late arrival system with
 * delayed access. Targets are the closed form of Dqsys_geogeo1 and the MATLAB
 * numbers recorded in _kb/06-solver-catalog.md, which LDES slotted independently
 * corroborated.
 */
public class SolverMAMDiscreteTimeTest {

    private static Network geoGeo1(double a, double s) {
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

    private static NetworkAvgTable solve(Network model) {
        MAMOptions options = new MAMOptions();
        options.verbose = VerboseLevel.SILENT;
        SolverMAM solver = new SolverMAM(model, options);
        return solver.getAvgTable();
    }

    @Test
    public void geoGeo1IsSolvedExactly() {
        double a = 0.2, s = 0.5;
        GeoGeo1Result exact = Dqsys_geogeo1.dqsys_geogeo1(a, s, GeoGeo1Convention.LAS_DA);
        NetworkAvgTable t = solve(geoGeo1(a, s));

        List<Double> qlen = t.getQLen();
        List<Double> util = t.getUtil();
        List<Double> respt = t.getRespT();
        List<Double> tput = t.getTput();

        // Row 0 is the Source, row 1 is the Queue
        assertEquals(exact.getMeanQueueLength(), qlen.get(1), 1e-6, "mean queue length");
        assertEquals(exact.getUtilization(), util.get(1), 1e-6, "utilization");
        assertEquals(exact.getMeanSojournTime(), respt.get(1), 1e-6, "mean sojourn time");
        assertEquals(exact.getThroughput(), tput.get(1), 1e-6, "throughput");
    }

    @Test
    public void detServiceStaysOnTheLattice() {
        // MATLAB SolverMAM reports 0.466667 here, LDES slotted 0.466615
        Network model = new Network("GeoDet1");
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");
        OpenClass jobClass = new OpenClass(model, "Class1");
        source.setArrival(jobClass, new Geometric(0.2));
        queue.setService(jobClass, new Det(2.0));
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(jobClass, jobClass, source, queue, 1.0);
        P.set(jobClass, jobClass, queue, sink, 1.0);
        model.link(P);

        NetworkAvgTable t = solve(model);
        assertEquals(0.466667, t.getQLen().get(1), 1e-5, "Geo/Det(2)/1 mean queue length");
        assertEquals(0.4, t.getUtil().get(1), 1e-6, "Geo/Det(2)/1 utilization");
    }

    @Test
    public void dmapArrivalsReachTheStructAndAreSolved() {
        // Before the discrete-time path existed a DMAP was consumed by the
        // continuous machinery, which computes inv(-D0) where the law needs
        // inv(I-D0). MATLAB reports 0.700000 here, LDES slotted 0.698332.
        Matrix D0 = new Matrix(2, 2);
        D0.set(0, 0, 0.5); D0.set(0, 1, 0.2); D0.set(1, 0, 0.1); D0.set(1, 1, 0.6);
        Matrix D1 = new Matrix(2, 2);
        D1.set(0, 0, 0.25); D1.set(0, 1, 0.05); D1.set(1, 0, 0.1); D1.set(1, 1, 0.2);

        Network model = new Network("DmapGeo1");
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");
        OpenClass jobClass = new OpenClass(model, "Class1");
        source.setArrival(jobClass, new DMAP(D0, D1));
        queue.setService(jobClass, new Geometric(0.6));
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(jobClass, jobClass, source, queue, 1.0);
        P.set(jobClass, jobClass, queue, sink, 1.0);
        model.link(P);

        NetworkAvgTable t = solve(model);
        assertEquals(0.700000, t.getQLen().get(1), 1e-5, "DMAP/Geo/1 mean queue length");
        assertEquals(0.5, t.getUtil().get(1), 1e-6, "DMAP/Geo/1 utilization");
    }

    @Test
    public void tandemReproducesTheDiscreteBurkeResult() {
        // The stationary departure stream of a Geo/Geo/1 queue is Bernoulli, so
        // the second queue is EXACT: 0.8 = A(1-A)/(S-A) with A = 0.2, S = 0.4.
        // MATLAB reports 0.5333/0.8000, LDES slotted 0.5363/0.8038.
        Network model = new Network("DTtandem");
        Source source = new Source(model, "Source");
        Queue q1 = new Queue(model, "Queue1", SchedStrategy.FCFS);
        Queue q2 = new Queue(model, "Queue2", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");
        OpenClass jobClass = new OpenClass(model, "Class1");
        source.setArrival(jobClass, new Geometric(0.2));
        q1.setService(jobClass, new Geometric(0.5));
        q2.setService(jobClass, new Geometric(0.4));
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(jobClass, jobClass, source, q1, 1.0);
        P.set(jobClass, jobClass, q1, q2, 1.0);
        P.set(jobClass, jobClass, q2, sink, 1.0);
        model.link(P);

        NetworkAvgTable t = solve(model);
        assertEquals(0.533333, t.getQLen().get(1), 1e-4, "first queue, exact Geo/Geo/1");
        assertEquals(0.800000, t.getQLen().get(2), 1e-3, "second queue, exact by discrete-time Burke");
        assertEquals(0.4, t.getUtil().get(1), 1e-6, "first queue utilization");
        assertEquals(0.5, t.getUtil().get(2), 1e-6, "second queue utilization");
    }

    @Test
    public void continuousModelsStillTakeTheContinuousPath() {
        // A guard against the detection swallowing ordinary models: an
        // exponential Source and service must not be read as slotted.
        Network model = new Network("MM1");
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");
        OpenClass jobClass = new OpenClass(model, "Class1");
        source.setArrival(jobClass, new jline.lang.processes.Exp(0.2));
        queue.setService(jobClass, new jline.lang.processes.Exp(0.5));
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(jobClass, jobClass, source, queue, 1.0);
        P.set(jobClass, jobClass, queue, sink, 1.0);
        model.link(P);

        NetworkAvgTable t = solve(model);
        // M/M/1 with rho = 0.4: E[N] = rho/(1-rho)
        assertEquals(0.4 / 0.6, t.getQLen().get(1), 1e-4, "M/M/1 mean queue length");
        assertTrue(t.getUtil().get(1) > 0.39 && t.getUtil().get(1) < 0.41, "M/M/1 utilization");
    }
}
