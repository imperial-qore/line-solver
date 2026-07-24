/*
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 */
package jline.inference;

import jline.lang.RoutingMatrix;
import jline.lang.Network;
import jline.lang.ClosedClass;
import jline.lang.Event;
import jline.lang.constant.EventType;
import jline.lang.constant.MetricType;
import jline.lang.constant.SchedStrategy;
import jline.lang.constant.SolverType;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Station;
import jline.lang.processes.Exp;
import jline.solvers.NetworkAvgTable;
import jline.solvers.SolverOptions;
import jline.solvers.mva.SolverMVA;
import jline.solvers.ssa.SolverSSA;
import jline.solvers.ssa.SampleNodeState;
import jline.inference.lang.ConditionEvent;
import jline.inference.lang.EstimatorOptions;
import jline.inference.lang.ParamEstimator;
import jline.inference.lang.SampledMetric;
import jline.util.Maths;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

import static org.junit.jupiter.api.Assertions.*;

/** One test per online-parameter-estimation method (mcmc/ubr/qmle/ubo/ekf/mlps/fmlps/erps/mle/gibbs). */
public class InferenceEstimatorsTest {

    @BeforeAll
    public static void setup() {
        Maths.setRandomNumbersMatlab(true);
    }

    @AfterAll
    public static void teardown() {
        Maths.setRandomNumbersMatlab(false);
    }

    @Test
    public void testMcmcClosed() throws Exception {
        Random rng = new Random(1);
        double trueDemand = 0.1;

        Network model = new Network("model");
        Delay delay = new Delay(model, "Delay");
        Queue queue = new Queue(model, "Queue1", SchedStrategy.PS);
        ClosedClass class1 = new ClosedClass(model, "Class1", 1, delay, 0);

        delay.setService(class1, Exp.fitMean(1.0));
        queue.setService(class1, Exp.fitMean(trueDemand));

        RoutingMatrix P = model.initRoutingMatrix();
        P.set(class1, class1, delay, queue, 1.0);
        P.set(class1, class1, queue, delay, 1.0);
        model.link(P);

        // Get true steady-state queue length from MVA
        SolverMVA solverMva = new SolverMVA(model);
        Matrix trueQLen = solverMva.getAvgQLen();
        double trueQLenQueue = trueQLen.get(1, 0); // Queue station index = 1 (0-based)

        // Generate model-consistent queue length samples
        int n = 5000;
        double[] ts = new double[n];
        for (int i = 0; i < n; i++) {
            ts[i] = (double) (i + 1);
        }
        double[] qlenSamples = new double[n];
        for (int i = 0; i < n; i++) {
            qlenSamples[i] = trueQLenQueue + rng.nextDouble() * 0.005 - 0.0025;
        }

        // Reset service for estimation
        queue.setService(class1, new Exp(Double.NaN));

        EstimatorOptions options = ParamEstimator.defaultOptions();
        options.method = "mcmc";
        ParamEstimator se = new ParamEstimator(model, options);

        se.addSamples(new SampledMetric(MetricType.QLen, ts, qlenSamples, queue));
        se.interpolate();
        List<Station> nodes = new ArrayList<Station>();
        nodes.add(queue);
        Matrix estVal = se.estimateAt(nodes);

        double est = estVal.get(0, 0);
        double relErr = Math.abs(est - trueDemand) / trueDemand;
        assertTrue(relErr < 0.10,
            "MCMC: estimated " + est + ", relative error " + (relErr * 100) + "% exceeds 10%");

        SolverMVA solver = new SolverMVA(model);
        NetworkAvgTable avgTable = solver.getAvgTable();
        assertNotNull(avgTable);
    }

    @Test
    public void testUbrClosed() throws Exception {
        Random rng = new Random(1);

        Network model = new Network("model");
        Delay delay = new Delay(model, "Delay");
        Queue queue = new Queue(model, "Queue1", SchedStrategy.PS);
        ClosedClass class1 = new ClosedClass(model, "Class1", 1, delay, 0);
        ClosedClass class2 = new ClosedClass(model, "Class2", 2, delay, 0);

        delay.setService(class1, Exp.fitMean(1.0));
        delay.setService(class2, Exp.fitMean(1.0));
        queue.setService(class1, new Exp(Double.NaN));
        queue.setService(class2, new Exp(Double.NaN));

        RoutingMatrix P = model.initRoutingMatrix();
        P.set(class1, class1, delay, queue, 1.0);
        P.set(class1, class1, queue, delay, 1.0);
        P.set(class2, class2, delay, queue, 1.0);
        P.set(class2, class2, queue, delay, 1.0);
        model.link(P);

        int n = 1000;
        double[] ts = new double[n];
        double[] arvr1 = new double[n];
        double[] arvr2 = new double[n];
        double[] util = new double[n];
        for (int i = 0; i < n; i++) {
            ts[i] = (double) (i + 1);
            arvr1[i] = 1.5 + rng.nextDouble() * 0.1;
            arvr2[i] = 2.0 + rng.nextDouble() * 0.1;
            util[i] = 0.2 * arvr1[i] + 0.4 * arvr2[i];
        }

        EstimatorOptions options = ParamEstimator.defaultOptions();
        options.method = "ubr";
        ParamEstimator se = new ParamEstimator(model, options);

        se.addSamples(new SampledMetric(MetricType.ArvR, ts, arvr1, queue, class1));
        se.addSamples(new SampledMetric(MetricType.ArvR, ts, arvr2, queue, class2));
        se.addSamples(new SampledMetric(MetricType.Util, ts, util, queue));
        se.interpolate();
        List<Station> nodes = Arrays.<Station>asList(queue);
        Matrix estVal = se.estimateAt(nodes);

        double[] trueDemands = new double[]{0.2, 0.4};
        for (int r = 0; r < trueDemands.length; r++) {
            double relErr = Math.abs(estVal.get(0, r) - trueDemands[r]) / trueDemands[r];
            assertTrue(relErr < 0.10,
                "UBR: class " + (r + 1) + " estimated " + estVal.get(0, r) + ", " +
                "relative error " + (relErr * 100) + "% exceeds 10%");
        }

        // Verify model solves after estimation
        SolverMVA solver = new SolverMVA(model);
        Object avgTable = solver.getAvgTable();
        assertNotNull(avgTable);
    }

    @Test
    public void testQmleClosed() throws Exception {
        Random rng = new Random(1);
        double[] trueDemands = new double[]{0.2, 0.4};

        Network model = new Network("model");
        Delay delay = new Delay(model, "Delay");
        Queue queue = new Queue(model, "Queue1", SchedStrategy.PS);
        ClosedClass class1 = new ClosedClass(model, "Class1", 2, delay, 0);
        ClosedClass class2 = new ClosedClass(model, "Class2", 3, delay, 0);

        delay.setService(class1, Exp.fitMean(1.0));
        delay.setService(class2, Exp.fitMean(1.0));
        queue.setService(class1, Exp.fitMean(trueDemands[0]));
        queue.setService(class2, Exp.fitMean(trueDemands[1]));

        RoutingMatrix P = model.initRoutingMatrix();
        P.set(class1, class1, delay, queue, 1.0);
        P.set(class1, class1, queue, delay, 1.0);
        P.set(class2, class2, delay, queue, 1.0);
        P.set(class2, class2, queue, delay, 1.0);
        model.link(P);

        // Get true steady-state queue lengths from MVA
        SolverMVA solverMva = new SolverMVA(model);
        Matrix trueQLen = solverMva.getAvgQLen();
        double trueQLen1 = trueQLen.get(1, 0); // Queue station, Class 1
        double trueQLen2 = trueQLen.get(1, 1); // Queue station, Class 2

        // Generate model-consistent queue length samples
        int n = 5000;
        double[] ts = new double[n];
        double[] qlen1 = new double[n];
        double[] qlen2 = new double[n];
        for (int i = 0; i < n; i++) {
            ts[i] = (double) (i + 1);
            qlen1[i] = trueQLen1 + rng.nextDouble() * 0.02 - 0.01;
            qlen2[i] = trueQLen2 + rng.nextDouble() * 0.02 - 0.01;
        }

        // Reset service for estimation
        queue.setService(class1, new Exp(Double.NaN));
        queue.setService(class2, new Exp(Double.NaN));

        EstimatorOptions options = ParamEstimator.defaultOptions();
        options.method = "qmle";
        ParamEstimator se = new ParamEstimator(model, options);

        se.addSamples(new SampledMetric(MetricType.QLen, ts, qlen1, queue, class1));
        se.addSamples(new SampledMetric(MetricType.QLen, ts, qlen2, queue, class2));
        se.interpolate();
        List<Station> nodes = Arrays.<Station>asList(queue);
        Matrix estVal = se.estimateAt(nodes);

        for (int r = 0; r < trueDemands.length; r++) {
            double relErr = Math.abs(estVal.get(0, r) - trueDemands[r]) / trueDemands[r];
            assertTrue(relErr < 0.10,
                "QMLE: class " + (r + 1) + " estimated " + estVal.get(0, r) + ", " +
                "relative error " + (relErr * 100) + "% exceeds 10%");
        }

        SolverMVA solver = new SolverMVA(model);
        Object avgTable = solver.getAvgTable();
        assertNotNull(avgTable);
    }

    @Test
    public void testUboClosed() throws Exception {
        Random rng = new Random(1);

        Network model = new Network("model");
        Delay delay = new Delay(model, "Delay");
        Queue queue = new Queue(model, "Queue1", SchedStrategy.PS);
        ClosedClass class1 = new ClosedClass(model, "Class1", 1, delay, 0);
        ClosedClass class2 = new ClosedClass(model, "Class2", 3, delay, 0);

        delay.setService(class1, Exp.fitMean(1.0));
        delay.setService(class2, Exp.fitMean(1.0));
        queue.setService(class1, new Exp(Double.NaN));
        queue.setService(class2, new Exp(Double.NaN));

        RoutingMatrix P = model.initRoutingMatrix();
        P.set(class1, class1, delay, queue, 1.0);
        P.set(class1, class1, queue, delay, 1.0);
        P.set(class2, class2, delay, queue, 1.0);
        P.set(class2, class2, queue, delay, 1.0);
        model.link(P);

        int n = 1000;
        double[] ts = new double[n];
        double[] arvr1 = new double[n];
        double[] arvr2 = new double[n];
        double[] util = new double[n];
        double[] respt1 = new double[n];
        double[] respt2 = new double[n];
        for (int i = 0; i < n; i++) {
            ts[i] = (double) (i + 1);
            arvr1[i] = 2.0 - rng.nextDouble() * 0.15;
            arvr2[i] = 1.0 - rng.nextDouble() * 0.15;
            util[i] = 0.1 * arvr1[i] + 0.3 * arvr2[i];
            respt1[i] = 0.1 / (1.0 - util[i]);
            respt2[i] = 0.3 / (1.0 - util[i]);
        }

        EstimatorOptions options = ParamEstimator.defaultOptions();
        options.method = "ubo";
        ParamEstimator se = new ParamEstimator(model, options);

        se.addSamples(new SampledMetric(MetricType.ArvR, ts, arvr1, queue, class1));
        se.addSamples(new SampledMetric(MetricType.ArvR, ts, arvr2, queue, class2));
        se.addSamples(new SampledMetric(MetricType.RespT, ts, respt1, queue, class1));
        se.addSamples(new SampledMetric(MetricType.RespT, ts, respt2, queue, class2));
        se.addSamples(new SampledMetric(MetricType.Util, ts, util, queue));
        se.interpolate();
        List<Station> nodes = Arrays.<Station>asList(queue);
        Matrix estVal = se.estimateAt(nodes);

        double[] trueDemands = new double[]{0.1, 0.3};
        // estVal may have multiple rows; use last row
        int lastRow = estVal.getNumRows() - 1;
        for (int r = 0; r < trueDemands.length; r++) {
            double relErr = Math.abs(estVal.get(lastRow, r) - trueDemands[r]) / trueDemands[r];
            assertTrue(relErr < 0.10,
                "UBO: class " + (r + 1) + " estimated " + estVal.get(lastRow, r) + ", " +
                "relative error " + (relErr * 100) + "% exceeds 10%");
        }

        SolverMVA solver = new SolverMVA(model);
        Object avgTable = solver.getAvgTable();
        assertNotNull(avgTable);
    }

    @Test
    public void testEkfClosed() throws Exception {
        Random rng = new Random(1);
        double trueDemand = 0.3;

        Network model = new Network("model");
        Delay delay = new Delay(model, "Delay");
        Queue queue = new Queue(model, "Queue1", SchedStrategy.PS);
        ClosedClass class1 = new ClosedClass(model, "Class1", 5, delay, 0);

        delay.setService(class1, Exp.fitMean(1.0));
        queue.setService(class1, Exp.fitMean(trueDemand));

        RoutingMatrix P = model.initRoutingMatrix();
        P.set(class1, class1, delay, queue, 1.0);
        P.set(class1, class1, queue, delay, 1.0);
        model.link(P);

        // Get true steady-state metrics from MVA
        SolverMVA solverMva = new SolverMVA(model);
        Matrix trueRespT = solverMva.getAvgRespT();
        Matrix trueUtil = solverMva.getAvgUtil();
        Matrix trueTput = solverMva.getAvgTput();

        int stIdx = queue.getStationIdx();
        double trueR = trueRespT.get(stIdx, 0);
        double trueU = trueUtil.get(stIdx, 0);
        double trueX = trueTput.get(stIdx, 0);

        // Generate noisy dataset
        int n = 1000;
        double[] ts = new double[n];
        for (int i = 0; i < n; i++) {
            ts[i] = (double) (i + 1);
        }
        double noiseScale = 0.05;
        double[] arvr = new double[n];
        double[] respt = new double[n];
        double[] util = new double[n];
        for (int i = 0; i < n; i++) {
            arvr[i] = trueX + (rng.nextDouble() - 0.5) * noiseScale * trueX;
        }
        for (int i = 0; i < n; i++) {
            respt[i] = trueR + (rng.nextDouble() - 0.5) * noiseScale * trueR;
        }
        for (int i = 0; i < n; i++) {
            util[i] = trueU + (rng.nextDouble() - 0.5) * noiseScale * trueU;
        }

        // Reset service for estimation
        queue.setService(class1, new Exp(Double.NaN));

        EstimatorOptions options = ParamEstimator.defaultOptions();
        options.method = "ekf";
        ParamEstimator se = new ParamEstimator(model, options);

        se.addSamples(new SampledMetric(MetricType.ArvR, ts, arvr, queue, class1));
        se.addSamples(new SampledMetric(MetricType.RespT, ts, respt, queue, class1));
        se.addSamples(new SampledMetric(MetricType.Util, ts, util, queue));
        List<Station> nodes = new ArrayList<Station>();
        nodes.add(queue);
        Matrix estVal = se.estimateAt(nodes);

        double est = estVal.get(0, 0);
        double relErr = Math.abs(est - trueDemand) / trueDemand;
        assertTrue(relErr < 0.10,
            "EKF: estimated " + est + ", relative error " + (relErr * 100) + "% exceeds 10%");

        SolverMVA solver = new SolverMVA(model);
        NetworkAvgTable avgTable = solver.getAvgTable();
        assertNotNull(avgTable);
    }

    @Test
    public void testMlpsClosed() throws Exception {
        double trueDemand = 0.5;

        Network model = new Network("model");
        Delay delay = new Delay(model, "Delay");
        Queue queue = new Queue(model, "Queue1", SchedStrategy.PS);
        ClosedClass class1 = new ClosedClass(model, "Class1", 1, delay, 0);

        delay.setService(class1, Exp.fitMean(1.0));
        queue.setService(class1, Exp.fitMean(trueDemand));

        RoutingMatrix P = model.initRoutingMatrix();
        P.set(class1, class1, delay, queue, 1.0);
        P.set(class1, class1, queue, delay, 1.0);
        model.link(P);

        // Generate trace data from SSA simulation
        SolverOptions ssaOptions = new SolverOptions(SolverType.SSA);
        ssaOptions.seed = 1;
        ssaOptions.samples = 10000;
        SolverSSA solverSsa = new SolverSSA(model, ssaOptions);
        SampleNodeState samplePath = solverSsa.sample(queue, 10000);

        // Extract arrival and departure times at the queue for class 1
        ArrayList<Double> arvTimes = new ArrayList<Double>();
        ArrayList<Double> depTimes = new ArrayList<Double>();
        for (Event ev : samplePath.event) {
            if (ev.getNode() == queue.getNodeIndex() && ev.getJobClass() == class1.getIndex() - 1) {
                if (ev.getEvent() == EventType.ARV) {
                    arvTimes.add(ev.getT());
                } else if (ev.getEvent() == EventType.DEP) {
                    depTimes.add(ev.getT());
                }
            }
        }
        int n = Math.min(arvTimes.size(), depTimes.size());
        double[] arrivalTimes = new double[n];
        double[] responseTimes = new double[n];
        for (int i = 0; i < n; i++) {
            arrivalTimes[i] = arvTimes.get(i);
            responseTimes[i] = depTimes.get(i) - arvTimes.get(i);
        }

        // Reset service for estimation
        queue.setService(class1, new Exp(Double.NaN));

        // Create trace-format SampledMetric objects
        SampledMetric arvData = new SampledMetric(MetricType.ArvR, arrivalTimes, arrivalTimes, queue, class1);
        arvData.setTrace();
        SampledMetric rtData = new SampledMetric(MetricType.RespT, arrivalTimes, responseTimes, queue, class1);
        rtData.setTrace();

        EstimatorOptions options = ParamEstimator.defaultOptions();
        options.method = "mlps";
        ParamEstimator se = new ParamEstimator(model, options);
        se.addSamples(arvData);
        se.addSamples(rtData);
        se.interpolate();
        List<Station> nodes = Arrays.<Station>asList(queue);
        Matrix estVal = se.estimateAt(nodes);

        double est = estVal.get(0, 0);
        double relErr = Math.abs(est - trueDemand) / trueDemand;
        assertTrue(relErr < 0.10,
            "MLPS: estimated " + est + ", relative error " + (relErr * 100) + "% exceeds 10%");

        SolverMVA solver = new SolverMVA(model);
        Object avgTable = solver.getAvgTable();
        assertNotNull(avgTable);
    }

    @Test
    public void testFmlpsClosed() throws Exception {
        double trueDemand = 0.5;

        Network model = new Network("model");
        Delay delay = new Delay(model, "Delay");
        Queue queue = new Queue(model, "Queue1", SchedStrategy.PS);
        ClosedClass class1 = new ClosedClass(model, "Class1", 1, delay, 0);

        delay.setService(class1, Exp.fitMean(1.0));
        queue.setService(class1, Exp.fitMean(trueDemand));

        RoutingMatrix P = model.initRoutingMatrix();
        P.set(class1, class1, delay, queue, 1.0);
        P.set(class1, class1, queue, delay, 1.0);
        model.link(P);

        // Generate trace data from SSA simulation
        SolverOptions ssaOptions = new SolverOptions(SolverType.SSA);
        ssaOptions.seed = 1;
        // Runtime is dominated by infer_fmlps, whose cost is linear in the trace
        // length (one fluid ODE likelihood solve per trace point per fmincon
        // evaluation). 1000 samples yield n=500 trace points, ~23 s, 4.7% error.
        ssaOptions.samples = 1000;
        SolverSSA solverSsa = new SolverSSA(model, ssaOptions);
        SampleNodeState samplePath = solverSsa.sample(queue, 1000);

        // Extract arrival and departure times at the queue for class 1
        List<Double> arvTimes = new ArrayList<Double>();
        List<Double> depTimes = new ArrayList<Double>();
        for (Event ev : samplePath.event) {
            if (ev.getNode() == queue.getNodeIndex() && ev.getJobClass() == class1.getIndex() - 1) {
                if (ev.getEvent() == EventType.ARV) {
                    arvTimes.add(ev.getT());
                } else if (ev.getEvent() == EventType.DEP) {
                    depTimes.add(ev.getT());
                }
            }
        }
        int n = Math.min(arvTimes.size(), depTimes.size());
        double[] arrivalTimes = new double[n];
        double[] responseTimes = new double[n];
        for (int i = 0; i < n; i++) {
            arrivalTimes[i] = arvTimes.get(i);
        }
        for (int i = 0; i < n; i++) {
            responseTimes[i] = depTimes.get(i) - arvTimes.get(i);
        }

        // Reset service for estimation
        queue.setService(class1, new Exp(Double.NaN));

        // Create trace-format SampledMetric objects
        SampledMetric arvData = new SampledMetric(MetricType.ArvR, arrivalTimes, arrivalTimes, queue, class1);
        arvData.setTrace();
        SampledMetric rtData = new SampledMetric(MetricType.RespT, arrivalTimes, responseTimes, queue, class1);
        rtData.setTrace();

        EstimatorOptions options = ParamEstimator.defaultOptions();
        options.method = "fmlps";
        ParamEstimator se = new ParamEstimator(model, options);
        se.addSamples(arvData);
        se.addSamples(rtData);
        se.interpolate();
        List<Station> nodes = new ArrayList<Station>();
        nodes.add(queue);
        Matrix estVal = se.estimateAt(nodes);

        double est = estVal.get(0, 0);
        double relErr = Math.abs(est - trueDemand) / trueDemand;
        assertTrue(relErr < 0.10,
            "FMLPS: estimated " + est + ", relative error " + (relErr * 100) + "% exceeds 10%");

        SolverMVA solver = new SolverMVA(model);
        NetworkAvgTable avgTable = solver.getAvgTable();
        assertNotNull(avgTable);
    }

    @Test
    public void testErpsClosed() throws Exception {
        Random rng = new Random(1);

        Network model = new Network("model");
        Delay delay = new Delay(model, "Delay");
        Queue queue = new Queue(model, "Queue1", SchedStrategy.PS);
        ClosedClass class1 = new ClosedClass(model, "Class1", 1, delay, 0);
        ClosedClass class2 = new ClosedClass(model, "Class2", 3, delay, 0);

        delay.setService(class1, Exp.fitMean(1.0));
        delay.setService(class2, Exp.fitMean(1.0));
        queue.setService(class1, new Exp(Double.NaN));
        queue.setService(class2, new Exp(Double.NaN));

        RoutingMatrix P = model.initRoutingMatrix();
        P.set(class1, class1, delay, queue, 1.0);
        P.set(class1, class1, queue, delay, 1.0);
        P.set(class2, class2, delay, queue, 1.0);
        P.set(class2, class2, queue, delay, 1.0);
        model.link(P);

        int n = 1000;
        double[] ts = new double[n];
        for (int i = 0; i < n; i++) {
            ts[i] = (double) (i + 1);
        }
        double[] arvr1 = new double[n];
        double[] arvr2 = new double[n];
        for (int i = 0; i < n; i++) {
            arvr1[i] = 1.0 - rng.nextDouble() * 0.15;
        }
        for (int i = 0; i < n; i++) {
            arvr2[i] = 2.0 - rng.nextDouble() * 0.15;
        }
        double[] util = new double[n];
        double[] respt1 = new double[n];
        double[] respt2 = new double[n];
        double[] aqlen1 = new double[n];
        double[] aqlen2 = new double[n];
        for (int i = 0; i < n; i++) {
            util[i] = 0.1 * arvr1[i] + 0.3 * arvr2[i];
        }
        for (int i = 0; i < n; i++) {
            respt1[i] = 0.1 / (1.0 - util[i]);
        }
        for (int i = 0; i < n; i++) {
            respt2[i] = 0.3 / (1.0 - util[i]);
        }
        for (int i = 0; i < n; i++) {
            aqlen1[i] = 1.0 + util[i] / (1.0 - util[i]);
        }
        for (int i = 0; i < n; i++) {
            aqlen2[i] = 1.0 + util[i] / (1.0 - util[i]);
        }

        EstimatorOptions options = ParamEstimator.defaultOptions();
        options.method = "erps";
        ParamEstimator se = new ParamEstimator(model, options);

        SampledMetric aql1 = new SampledMetric(MetricType.QLen, ts, aqlen1, queue);
        aql1.setConditional(new ConditionEvent(queue, class1, EventType.ARV));
        SampledMetric aql2 = new SampledMetric(MetricType.QLen, ts, aqlen2, queue);
        aql2.setConditional(new ConditionEvent(queue, class2, EventType.ARV));

        se.addSamples(aql1);
        se.addSamples(aql2);
        se.addSamples(new SampledMetric(MetricType.RespT, ts, respt1, queue, class1));
        se.addSamples(new SampledMetric(MetricType.RespT, ts, respt2, queue, class2));
        se.interpolate();
        List<Station> nodes = new ArrayList<Station>();
        nodes.add(queue);
        Matrix estVal = se.estimateAt(nodes);

        double[] trueDemands = new double[]{0.1, 0.3};
        for (int r = 0; r < trueDemands.length; r++) {
            double relErr = Math.abs(estVal.get(0, r) - trueDemands[r]) / trueDemands[r];
            assertTrue(relErr < 0.10,
                "ERPS: class " + (r + 1) + " estimated " + estVal.get(0, r) + ", " +
                "relative error " + (relErr * 100) + "% exceeds 10%");
        }

        SolverMVA solver = new SolverMVA(model);
        NetworkAvgTable avgTable = solver.getAvgTable();
        assertNotNull(avgTable);
    }

    @Test
    public void testMleClosed() throws Exception {
        Random rng = new Random(1);
        double[] trueDemands = new double[]{0.1, 0.3};

        Network model = new Network("model");
        Delay delay = new Delay(model, "Delay");
        Queue queue = new Queue(model, "Queue1", SchedStrategy.PS);
        ClosedClass class1 = new ClosedClass(model, "Class1", 1, delay, 0);
        ClosedClass class2 = new ClosedClass(model, "Class2", 3, delay, 0);

        delay.setService(class1, Exp.fitMean(1.0));
        delay.setService(class2, Exp.fitMean(1.0));
        queue.setService(class1, Exp.fitMean(trueDemands[0]));
        queue.setService(class2, Exp.fitMean(trueDemands[1]));

        RoutingMatrix P = model.initRoutingMatrix();
        P.set(class1, class1, delay, queue, 1.0);
        P.set(class1, class1, queue, delay, 1.0);
        P.set(class2, class2, delay, queue, 1.0);
        P.set(class2, class2, queue, delay, 1.0);
        model.link(P);

        // Get true steady-state metrics from MVA
        SolverMVA solverMva = new SolverMVA(model);
        Matrix trueRespT = solverMva.getAvgRespT();
        Matrix trueUtil = solverMva.getAvgUtil();
        Matrix trueTput = solverMva.getAvgTput();

        int stIdx = queue.getStationIdx();
        double trueR1 = trueRespT.get(stIdx, 0);
        double trueR2 = trueRespT.get(stIdx, 1);
        double trueU = trueUtil.get(stIdx, 0) + trueUtil.get(stIdx, 1);
        double trueX1 = trueTput.get(stIdx, 0);
        double trueX2 = trueTput.get(stIdx, 1);

        // Generate noisy dataset
        int n = 1000;
        double[] ts = new double[n];
        for (int i = 0; i < n; i++) {
            ts[i] = (double) (i + 1);
        }
        double noiseScale = 0.05;
        double[] arvr1 = new double[n];
        double[] arvr2 = new double[n];
        double[] respt1 = new double[n];
        double[] respt2 = new double[n];
        double[] utilSamples = new double[n];
        for (int i = 0; i < n; i++) {
            arvr1[i] = trueX1 + (rng.nextDouble() - 0.5) * noiseScale * trueX1;
        }
        for (int i = 0; i < n; i++) {
            arvr2[i] = trueX2 + (rng.nextDouble() - 0.5) * noiseScale * trueX2;
        }
        for (int i = 0; i < n; i++) {
            respt1[i] = trueR1 + (rng.nextDouble() - 0.5) * noiseScale * trueR1;
        }
        for (int i = 0; i < n; i++) {
            respt2[i] = trueR2 + (rng.nextDouble() - 0.5) * noiseScale * trueR2;
        }
        for (int i = 0; i < n; i++) {
            utilSamples[i] = trueU + (rng.nextDouble() - 0.5) * noiseScale * trueU;
        }

        // Reset service for estimation
        queue.setService(class1, new Exp(Double.NaN));
        queue.setService(class2, new Exp(Double.NaN));

        EstimatorOptions options = ParamEstimator.defaultOptions();
        options.method = "mle";
        ParamEstimator se = new ParamEstimator(model, options);

        se.addSamples(new SampledMetric(MetricType.ArvR, ts, arvr1, queue, class1));
        se.addSamples(new SampledMetric(MetricType.ArvR, ts, arvr2, queue, class2));
        se.addSamples(new SampledMetric(MetricType.RespT, ts, respt1, queue, class1));
        se.addSamples(new SampledMetric(MetricType.RespT, ts, respt2, queue, class2));
        se.addSamples(new SampledMetric(MetricType.Util, ts, utilSamples, queue));
        se.interpolate();
        List<Station> nodes = new ArrayList<Station>();
        nodes.add(queue);
        Matrix estVal = se.estimateAt(nodes);

        for (int r = 0; r < trueDemands.length; r++) {
            double relErr = Math.abs(estVal.get(0, r) - trueDemands[r]) / trueDemands[r];
            assertTrue(relErr < 0.10,
                "MLE: class " + (r + 1) + " estimated " + estVal.get(0, r) + ", " +
                "relative error " + (relErr * 100) + "% exceeds 10%");
        }

        SolverMVA solver = new SolverMVA(model);
        NetworkAvgTable avgTable = solver.getAvgTable();
        assertNotNull(avgTable);
    }

    @Test
    public void testGibbsClosed() throws Exception {
        double trueDemand = 0.3;

        Network model = new Network("model");
        Delay delay = new Delay(model, "Delay");
        Queue queue = new Queue(model, "Queue1", SchedStrategy.PS);
        ClosedClass class1 = new ClosedClass(model, "Class1", 1, delay, 0);

        delay.setService(class1, Exp.fitMean(1.0));
        queue.setService(class1, Exp.fitMean(trueDemand));

        RoutingMatrix P = model.initRoutingMatrix();
        P.set(class1, class1, delay, queue, 1.0);
        P.set(class1, class1, queue, delay, 1.0);
        model.link(P);

        // Generate trace data from SSA simulation
        SolverOptions ssaOptions = new SolverOptions(SolverType.SSA);
        ssaOptions.seed = 1;
        ssaOptions.samples = 10000;
        SolverSSA solverSsa = new SolverSSA(model, ssaOptions);
        SampleNodeState samplePath = solverSsa.sample(queue, 10000);

        // Extract arrival and departure times at the queue for class 1
        List<Double> arvTimes = new ArrayList<Double>();
        List<Double> depTimes = new ArrayList<Double>();
        for (Event ev : samplePath.event) {
            if (ev.getNode() == queue.getNodeIndex() && ev.getJobClass() == class1.getIndex() - 1) {
                if (ev.getEvent() == EventType.ARV) {
                    arvTimes.add(ev.getT());
                } else if (ev.getEvent() == EventType.DEP) {
                    depTimes.add(ev.getT());
                }
            }
        }
        int n = Math.min(arvTimes.size(), depTimes.size());
        double[] arrivalTimes = new double[n];
        double[] responseTimes = new double[n];
        for (int i = 0; i < n; i++) {
            arrivalTimes[i] = arvTimes.get(i);
        }
        for (int i = 0; i < n; i++) {
            responseTimes[i] = depTimes.get(i) - arvTimes.get(i);
        }

        // Compute throughput from inter-departure times
        double[] depTimesArr = new double[n];
        for (int i = 0; i < n; i++) {
            depTimesArr[i] = depTimes.get(i);
        }
        double sumInterDep = 0.0;
        for (int i = 1; i < n; i++) {
            sumInterDep += depTimesArr[i] - depTimesArr[i - 1];
        }
        double meanInterDep = sumInterDep / (n - 1);
        double tputVal = 1.0 / meanInterDep;
        double[] tputSamples = new double[n];
        for (int i = 0; i < n; i++) {
            tputSamples[i] = tputVal;
        }

        // Reset service for estimation
        queue.setService(class1, new Exp(Double.NaN));

        // Create trace-format SampledMetric objects
        SampledMetric arvData = new SampledMetric(MetricType.ArvR, arrivalTimes, arrivalTimes, queue, class1);
        arvData.setTrace();
        SampledMetric rtData = new SampledMetric(MetricType.RespT, arrivalTimes, responseTimes, queue, class1);
        rtData.setTrace();
        SampledMetric tputData = new SampledMetric(MetricType.Tput, arrivalTimes, tputSamples, queue, class1);

        EstimatorOptions options = ParamEstimator.defaultOptions();
        options.method = "gibbs";
        ParamEstimator se = new ParamEstimator(model, options);
        se.addSamples(arvData);
        se.addSamples(rtData);
        se.addSamples(tputData);
        se.interpolate();
        List<Station> nodes = new ArrayList<Station>();
        nodes.add(queue);
        Matrix estVal = se.estimateAt(nodes);

        double est = estVal.get(0, 0);
        double relErr = Math.abs(est - trueDemand) / trueDemand;
        assertTrue(relErr < 0.10,
            "Gibbs: estimated " + est + ", relative error " + (relErr * 100) + "% exceeds 10%");

        SolverMVA solver = new SolverMVA(model);
        NetworkAvgTable avgTable = solver.getAvgTable();
        assertNotNull(avgTable);
    }
}
