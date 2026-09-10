/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.nc;

import jline.VerboseLevel;
import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.RoutingMatrix;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.processes.Exp;
import jline.io.LineCitations;
import jline.solvers.MethodType;
import jline.solvers.ctmc.SolverCTMC;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;

import java.util.Arrays;
import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * SolverNC's mixed load-dependent route must read the WHOLE rate lattice.
 *
 * <p>An open or mixed limited-load-dependent model goes to the
 * Bruell-Balbo-Afshari effective-capacity MVA (Pfqn_mvaldmx, method
 * <code>ncldmx</code>), which is exact. Two defects made it answer something
 * else:
 *
 * <ul>
 * <li>the rate row was cut at the CLOSED population. Pfqn_ldmx_ec infers the
 * limited-load-dependence level b_i as the first column equal to the last one
 * and treats every rate past it as saturated, so a c-server station was read as
 * saturated at min(n,c) with n&lt;c whenever c exceeded that population -- and
 * with no closed class at all the row collapsed to mu(1), one server.</li>
 * <li>Pfqn_mvaldmx could not run a purely open model at all: pprod on an empty
 * population bound never returns the -1 sentinel, so the lattice recursion
 * never terminated.</li>
 * </ul>
 *
 * <p>The oracles are the exact M/M/c law and SolverCTMC.
 */
public class NcMixedLdMultiserverTest {

    private static final int SERVERS = 3;
    private static final int LLD_WIDTH = 40;

    /** Exact mean number in an M/M/c queue. */
    private static double mmcQLen(double lambda, double mu, int c) {
        double a = lambda / mu;
        double rho = a / c;
        double s = 0.0;
        double fact = 1.0;
        for (int k = 0; k < c; k++) {
            if (k > 0) {
                fact *= k;
            }
            s += Math.pow(a, k) / fact;
        }
        double factc = fact * c;
        double p0 = 1.0 / (s + Math.pow(a, c) / (factc * (1 - rho)));
        return p0 * Math.pow(a, c) * rho / (factc * (1 - rho) * (1 - rho)) + a;
    }

    /**
     * Source -&gt; load-dependent queue -&gt; Sink, with an optional closed chain
     * through a delay. The multiserver is expressed as mu(n) = min(n, c), the
     * only form the load-dependent solver admits.
     */
    private static Network model(double lambda, int closedJobs) {
        Network model = new Network("mixed_ld");
        Source source = new Source(model, "Src");
        Queue queue = new Queue(model, "Q1", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Snk");
        OpenClass oclass = new OpenClass(model, "oc");
        source.setArrival(oclass, new Exp(lambda));
        queue.setService(oclass, new Exp(1.0));
        RoutingMatrix routing;
        if (closedJobs > 0) {
            Delay delay = new Delay(model, "D");
            ClosedClass cclass = new ClosedClass(model, "cc", closedJobs, delay);
            delay.setService(cclass, new Exp(1.0));
            queue.setService(cclass, new Exp(1.0));
            routing = model.initRoutingMatrix();
            routing.set(cclass, cclass, delay, queue, 1.0);
            routing.set(cclass, cclass, queue, delay, 1.0);
        } else {
            routing = model.initRoutingMatrix();
        }
        routing.set(oclass, oclass, source, queue, 1.0);
        routing.set(oclass, oclass, queue, sink, 1.0);
        model.link(routing);
        Matrix alpha = new Matrix(1, LLD_WIDTH);
        for (int n = 1; n <= LLD_WIDTH; n++) {
            alpha.set(0, n - 1, Math.min(n, SERVERS));
        }
        queue.setLoadDependence(alpha);
        return model;
    }

    @Test
    public void testPurelyOpenMultiserverIsTheMmcQueue() {
        // no closed class: the rate row must still describe c servers, not one
        double lambda = 1.5;
        SolverNC solver = new SolverNC(model(lambda, 0), "method", "exact",
                "verbose", VerboseLevel.SILENT);
        assertEquals(mmcQLen(lambda, 1.0, SERVERS), solver.getAvgQLen().get(1, 0), 1e-9);
        assertEquals(lambda, solver.getAvgTput().get(1, 0), 1e-9);
    }

    @Test
    public void testMixedModelWithServersAboveTheClosedPopulation() {
        // c = 3 with one closed job: the row used to be cut to mu(1)
        double lambda = 0.4;
        SolverNC nc = new SolverNC(model(lambda, 1), "method", "exact",
                "verbose", VerboseLevel.SILENT);
        SolverCTMC ctmc = new SolverCTMC(model(lambda, 1), "cutoff", 16,
                "verbose", VerboseLevel.SILENT);
        Matrix ncQ = nc.getAvgQLen();
        Matrix ctmcQ = ctmc.getAvgQLen();
        for (int r = 0; r < 2; r++) {
            assertEquals(ctmcQ.get(1, r), ncQ.get(1, r), 1e-3 * (1 + Math.abs(ctmcQ.get(1, r))),
                    "QLen of class " + r + " at the load-dependent queue");
        }
    }

    @Test
    public void testNcldmxIsClassifiedAsAnExactMethod() {
        // the route evaluates the product form itself, so the banner must say so
        assertEquals(MethodType.EXACT_DET, MethodType.of("SolverNC", "ncldmx"));
        assertEquals(MethodType.EXACT_DET, MethodType.of("SolverNC", "default/ncldmx"));
    }

    @Test
    public void testNcldmxCarriesItsCitation() {
        List<LineCitations.Citation> refs =
                LineCitations.citationsFor(Arrays.asList("ncldmx"));
        assertEquals(1, refs.size());
        assertEquals("AfBaBr84", refs.get(0).key);
        assertTrue(refs.get(0).ref.contains("Bruell"), refs.get(0).ref);
    }
}
