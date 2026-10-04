/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.ldes;

import jline.GlobalConstants;
import jline.VerboseLevel;
import jline.lang.JobClass;
import jline.lang.Network;
import jline.lang.NetworkStruct;
import jline.lang.OpenClass;
import jline.lang.RoutingMatrix;
import jline.lang.constant.ProcessType;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.processes.Bernoulli;
import jline.lang.processes.Binomial;
import jline.lang.processes.Distribution;
import jline.lang.processes.Exp;
import jline.lang.processes.Poisson;
import jline.solvers.NetworkAvgTable;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Coverage for the counting distributions (Bernoulli, Binomial, Poisson) used as
 * service times.
 *
 * These three were reachable in the engine's generator factory long before they
 * were declared anywhere else, so a model using one failed in whichever registry
 * happened to be checked first. The tests below exercise each registry in turn:
 * the struct representation, the feature gate, the analyzer's process-type list,
 * and the generator itself.
 *
 * They are supported on {0,1,...}, so the zero atom is the substantive question:
 * in continuous time it denotes a zero-length interval and is resolved to the
 * Immediate constant, and in slotted mode it is rejected because the lattice has
 * no such point.
 */
public class SolverLDESCountingDistTest {

    private static final int SEED = 23000;

    private static Network model(Distribution service) {
        Network model = new Network("counting");
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");
        OpenClass jobClass = new OpenClass(model, "Class1");
        source.setArrival(jobClass, new Exp(0.05));
        queue.setService(jobClass, service);
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(jobClass, jobClass, source, queue, 1.0);
        P.set(jobClass, jobClass, queue, sink, 1.0);
        model.link(P);
        return model;
    }

    private static LDESOptions options(int samples) {
        LDESOptions o = new LDESOptions();
        o.verbose = VerboseLevel.SILENT;
        o.seed = SEED;
        o.samples = samples;
        return o;
    }

    // ------------------------------------------------------------------
    // Struct representation (Station.getServiceRates)
    // ------------------------------------------------------------------

    @Test
    public void structCarriesTheMomentPair() {
        Distribution[] dists = {new Bernoulli(0.4), new Binomial(6, 0.5), new Poisson(3.0)};
        ProcessType[] types = {ProcessType.BERNOULLI, ProcessType.BINOMIAL, ProcessType.POISSON};
        for (int i = 0; i < dists.length; i++) {
            Distribution d = dists[i];
            Network m = model(d);
            NetworkStruct sn = m.getStruct(true);
            JobClass jc = m.getClassByName("Class1");

            // Without the Station arm these fell to the NaN default.
            jline.util.matrix.MatrixCell proc = sn.proc.get(m.getStations().get(1)).get(jc);
            assertEquals(d.getMean(), proc.get(0).get(0, 0), 1e-12,
                    d.getName() + " proc[0] must be the mean");
            assertEquals(d.getSCV(), proc.get(1).get(0, 0), 1e-12,
                    d.getName() + " proc[1] must be the SCV");
            assertEquals(types[i], sn.procid.get(m.getStations().get(1)).get(jc),
                    d.getName() + " procid");
            // sn.rates is populated from getRate(), whose contract is 1/mean.
            assertEquals(1.0 / d.getMean(), sn.rates.get(1), 1e-12,
                    d.getName() + " rate must be 1/mean");
        }
    }

    @Test
    public void poissonRateIsTheReciprocalOfTheMean() {
        // Poisson.getRate used to return lambda, i.e. the MEAN, which put the
        // mean into sn.rates and disagreed with MATLAB's Poisson.getRate.
        Poisson p = new Poisson(4.0);
        assertEquals(4.0, p.getMean(), 1e-12);
        assertEquals(0.25, p.getRate(), 1e-12, "getRate must be 1/mean");
        assertEquals(0.25, p.getSCV(), 1e-12, "Poisson SCV is 1/lambda");
        assertEquals(4.0, p.getVar(), 1e-12, "Poisson variance is lambda");
        assertEquals(0.5, p.getSkewness(), 1e-12, "Poisson skewness is 1/sqrt(lambda)");
    }

    // ------------------------------------------------------------------
    // The model is accepted and simulated
    // ------------------------------------------------------------------

    @Test
    public void countingServiceModelsSolve() {
        // lambda = 0.05 against means 0.4 / 3.0 / 3.0 keeps every case stable.
        Distribution[] dists = {new Bernoulli(0.4), new Binomial(6, 0.5), new Poisson(3.0)};
        for (Distribution d : dists) {
            NetworkAvgTable table = new SolverLDES(model(d), options(200000)).getAvgTable();
            double tput = table.getTput().get(1);
            double util = table.getUtil().get(1);
            assertEquals(0.05, tput, 0.05 * 0.05,
                    d.getName() + " throughput must match the arrival rate");
            // U = X * E[S] for a single server.
            assertEquals(0.05 * d.getMean(), util, Math.max(0.05 * 0.05 * d.getMean(), 1e-3),
                    d.getName() + " utilization must be X*E[S]");
            assertTrue(table.getRespT().get(1) > 0.0,
                    d.getName() + " response time must be measured, not defaulted");
        }
    }

    @Test
    public void zeroAtomBecomesAnImmediateInterval() {
        // Bernoulli(0.05) is 0 in 95% of draws, so a run that survives at all is
        // evidence the zero atom is being resolved rather than read as "never".
        // Mean service is 0.05, so utilization must stay tiny but positive.
        NetworkAvgTable table =
                new SolverLDES(model(new Bernoulli(0.05)), options(200000)).getAvgTable();
        double util = table.getUtil().get(1);
        assertTrue(util >= 0.0 && util < 0.02,
                "near-zero service must give near-zero utilization, got " + util);
        assertEquals(0.05, table.getTput().get(1), 0.05 * 0.1,
                "throughput must still equal the arrival rate");
        assertTrue(1.0 / GlobalConstants.Immediate > 0.0,
                "the immediate constant must be strictly positive to be schedulable");
    }

    // ------------------------------------------------------------------
    // Slotted mode rejects them
    // ------------------------------------------------------------------

    @Test
    public void slottedModeRejectsTheZeroAtom() {
        LDESOptions o = options(50000);
        o.setSlotted(true);
        RuntimeException thrown = assertThrows(RuntimeException.class,
                () -> new SolverLDES(model(new Bernoulli(0.4)), o).getAvgTable());
        String msg = rootMessage(thrown);
        assertTrue(msg.contains("slotted mode"),
                "the slotted rejection must name the mode; got: " + msg);
    }

    // ------------------------------------------------------------------
    // Feature registries
    // ------------------------------------------------------------------

    @Test
    public void allThreeAreDeclaredAndAccepted() {
        for (Distribution d : new Distribution[]{
                new Bernoulli(0.4), new Binomial(6, 0.5), new Poisson(3.0)}) {
            Network m = model(d);
            // Exercises the FeatureSet static name registry (an unregistered
            // name makes setTrue throw) and the LDES featset in one call.
            assertTrue(new SolverLDES(m, options(1000)).supports(m),
                    d.getName() + " must be accepted by the LDES feature gate");
        }
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
