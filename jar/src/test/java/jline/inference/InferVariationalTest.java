/*
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 */
package jline.inference;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import org.junit.jupiter.api.Test;

import java.util.Arrays;

import jline.inference.api.Infer_variational;
import jline.inference.api.VariationalOptions;
import jline.inference.api.VariationalResult;
import jline.inference.api.VariationalSpec;
import jline.inference.lang.ParamEstimator;
import jline.inference.lang.SampledMetric;
import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.constant.MetricType;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Station;
import jline.lang.processes.Exp;
import jline.util.matrix.Matrix;

/**
 * Variational inference for Markovian queueing networks (Perez-Casale, AAP
 * 53(3), 2021).
 *
 * <p>The estimator carries no random-number stream, so the fixture below must
 * reproduce the MATLAB and Python implementations digit for digit; the golden
 * values were taken from the MATLAB reference.</p>
 */
public class InferVariationalTest {

    /** Closed two-station loop of the paper's Section 6.1, on fixed readings. */
    private static VariationalSpec closedLoopSpec() {
        int N = 10;
        double lam = 0.5;
        double[] obsQ = {1, 2, 3, 2, 4, 3, 5, 4, 3, 4};
        VariationalSpec spec = new VariationalSpec();
        spec.arcs = new int[][]{{1, 2, 1}, {2, 1, 1}};
        spec.x0 = new double[][]{{N}, {0}};
        spec.sched = new int[]{0, 1};
        spec.nservers = new double[]{1, 1};
        spec.routeprob = new double[]{1, 1};
        spec.arcparam = new int[]{0, 1};
        spec.arcrate = new double[]{lam, Double.NaN};
        spec.alpha0 = new double[]{2.0};
        spec.beta0 = new double[]{1.0};
        spec.obsTimes = new double[10];
        spec.obsData = new double[10][2];
        for (int k = 0; k < 10; k++) {
            spec.obsTimes[k] = k + 1;
            spec.obsData[k][0] = N - obsQ[k];
            spec.obsData[k][1] = obsQ[k];
        }
        spec.obsRange = new double[]{N, N};
        spec.epsilon = 0.1;
        spec.capacity = new double[]{N, N};
        return spec;
    }

    @Test
    public void testClosedLoopMatchesMatlab() {
        VariationalOptions opt = new VariationalOptions();
        opt.ngrid = 51;
        opt.nsamples = 32;
        opt.ymax = 60;
        opt.iterMax = 5;
        opt.tol = 0.0;
        opt.delta = 1e-3;
        VariationalResult out = Infer_variational.infer_variational(closedLoopSpec(), opt);

        assertEquals(20.8738918926, out.alpha[0], 1e-8);
        assertEquals(10.4687500000, out.beta[0], 1e-8);
        assertEquals(1.9939240017, out.rates[0], 1e-8);
        assertEquals(0.5015230000, out.meanServiceTime[0], 1e-6);
        assertEquals(5, out.iter);

        double[] expectedBound = {-89.392514468, -164.0463074709, -163.7108480833,
                -174.6500856726, -163.2336393829};
        for (int i = 0; i < expectedBound.length; i++) {
            assertEquals(expectedBound[i], out.bound[i], 1e-6);
        }

        int ny = out.Y[0][0].length;
        double ey1 = 0.0;
        double ey2 = 0.0;
        for (int y = 0; y < ny; y++) {
            ey1 += out.Y[0][out.Y[0].length - 1][y] * y;
            ey2 += out.Y[1][out.Y[1].length - 1][y] * y;
        }
        assertEquals(22.7916783314, ey1, 1e-8);
        assertEquals(18.8738918926, ey2, 1e-8);

        double[] last = out.qlen[out.qlen.length - 1];
        assertEquals(6.0822135613, last[0], 1e-8);
        assertEquals(3.9177864387, last[1], 1e-8);
        // the two stations hold the whole closed population at every epoch
        assertEquals(10.0, last[0] + last[1], 1e-9);
    }

    /**
     * With a single transition the mean field is exact, so the marginal must
     * reproduce the transient of the underlying birth process. An infinite
     * server emptying at rate lam gives E[Y(t)] = N(1-exp(-lam t)).
     */
    @Test
    public void testSingleTransitionIsExact() {
        int N = 50;
        double lam = 0.1;
        double tmax = 20.0;
        VariationalSpec spec = new VariationalSpec();
        spec.arcs = new int[][]{{1, 0, 1}};
        spec.x0 = new double[][]{{N}};
        spec.sched = new int[]{0};
        spec.nservers = new double[]{1};
        spec.routeprob = new double[]{1};
        spec.arcparam = new int[]{0};
        spec.arcrate = new double[]{lam};
        spec.alpha0 = new double[0];
        spec.beta0 = new double[0];
        spec.obsTimes = new double[]{tmax};
        spec.obsData = new double[][]{{Double.NaN}};
        spec.obsRange = new double[]{N};
        spec.epsilon = 0.2;
        spec.capacity = new double[]{N};

        VariationalOptions opt = new VariationalOptions();
        opt.ngrid = 201;
        opt.nsamples = 50;
        opt.iterMax = 4;
        opt.tmax = tmax;
        VariationalResult out = Infer_variational.infer_variational(spec, opt);

        int ny = out.Y[0][0].length;
        double ey = 0.0;
        for (int y = 0; y < ny; y++) {
            ey += out.Y[0][out.Y[0].length - 1][y] * y;
        }
        assertEquals(N * (1 - Math.exp(-lam * tmax)), ey, 2e-2);
        assertTrue(out.tailmass < 1e-6, "truncation must not carry mass");
    }

    /**
     * The same fixture driven through ParamEstimator, so that the translation
     * from the network to the transition set is pinned to the api golden: same
     * transitions in the same order, same priors, same observations.
     */
    @Test
    public void testEstimatorReachesTheApiFixture() {
        int N = 10;
        double[] obsQ = {1, 2, 3, 2, 4, 3, 5, 4, 3, 4};
        Network model = new Network("vi");
        Delay think = new Delay(model, "Think");
        Queue queue = new Queue(model, "Q", SchedStrategy.FCFS);
        ClosedClass cl = new ClosedClass(model, "C", N, think, 0);
        think.setService(cl, new Exp(0.5));
        queue.setService(cl, new Exp(2.0));
        model.link(model.serialRouting(think, queue));

        double[] t = new double[10];
        double[] qd = new double[10];
        double[] dd = new double[10];
        for (int k = 0; k < 10; k++) {
            t[k] = k + 1;
            qd[k] = obsQ[k];
            dd[k] = N - obsQ[k];
        }
        ParamEstimator pe = new ParamEstimator(model);
        pe.addSamples(new SampledMetric(MetricType.QLen, t, dd, think, cl));
        pe.addSamples(new SampledMetric(MetricType.QLen, t, qd, queue, cl));
        pe.options.method = "vi";
        pe.options.epsilon = 0.1;
        pe.options.priorShape = 2.0;  // with the model rate 2.0 this gives Gamma(2,1)
        pe.options.variational.ngrid = 51;
        pe.options.variational.nsamples = 32;
        pe.options.variational.ymax = 60;
        pe.options.variational.iterMax = 5;
        pe.options.variational.tol = 0.0;
        pe.options.variational.delta = 1e-3;

        Matrix est = pe.estimateAt(Arrays.asList((Station) queue));
        assertEquals(20.8738918926, pe.options.posteriorAlpha[0], 1e-8);
        assertEquals(10.4687500000, pe.options.posteriorBeta[0], 1e-8);
        assertEquals(10.4687500000 / 20.8738918926, est.get(0, 0), 1e-9);
        // the estimate has replaced the service time the model started with
        assertEquals(est.get(0, 0), queue.getService(cl).getMean(), 1e-9);
    }
}
