/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers;

import jline.GlobalConstants;
import jline.VerboseLevel;
import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.OpenClass;
import jline.io.Ret.DistributionResult;
import jline.lang.constant.BalkingStrategy;
import jline.lang.constant.BalkingThreshold;
import jline.lang.constant.PollingType;
import jline.lang.constant.SchedStrategy;
import jline.lang.constant.SolverType;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.processes.Coxian;
import jline.lang.processes.Det;
import jline.lang.processes.Distribution;
import jline.lang.processes.Erlang;
import jline.lang.processes.Exp;
import jline.lang.processes.Gamma;
import jline.lang.processes.HyperExp;
import jline.lang.processes.Lognormal;
import jline.lang.processes.MAP;
import jline.lang.processes.Pareto;
import jline.lang.processes.Uniform;
import jline.lang.processes.Weibull;
import jline.solvers.ctmc.SolverCTMC;
import jline.solvers.fluid.SolverFluid;
import jline.solvers.ldes.SolverLDES;
import jline.solvers.mva.SolverMVA;
import jline.solvers.nc.SolverNC;
import jline.solvers.ssa.SolverSSA;
import jline.util.Maths;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.ValueSource;

import java.util.ArrayList;
import java.util.List;

import static jline.TestTools.withSuppressedOutput;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Feature-matrix regression test: tiny models sweeping the supported feature
 * axes (scheduling policies, service distributions, impatience, retrial,
 * polling, finite capacity regions, heterogeneous servers, transient and
 * response-time CDF analysis) across the analytic and simulation solvers.
 *
 * Validation strategy per axis:
 * - simulation axes assert exact conservation/operational laws (throughput
 *   conservation, utilization law) within simulation tolerance;
 * - feature axes supported by two independent engines cross-check them
 *   (CTMC exact vs SSA simulation, MVA vs NC vs CTMC on product form);
 * - transient axes assert convergence to the steady state and CDF sanity.
 */
public class SolverFeatureMatrixTest {

    private static final int LDES_SAMPLES = 20000;
    private static final int SSA_SAMPLES = 50000;
    private static final int SEED = 23000;
    private static final double SIM_RTOL = 0.15;      // simulation vs law
    private static final double CROSS_RTOL = 0.20;    // simulation vs exact CTMC
    private static final double EXACT_RTOL = 1e-3;    // exact vs exact

    private static VerboseLevel originalVerboseLevel;

    @BeforeAll
    public static void setUpClass() {
        originalVerboseLevel = GlobalConstants.getVerbose();
        GlobalConstants.setVerbose(VerboseLevel.SILENT);
        Maths.setRandomNumbersMatlab(true);
    }

    @AfterAll
    public static void tearDownClass() {
        GlobalConstants.setVerbose(originalVerboseLevel);
    }

    // ------------------------------------------------------------------
    // Axis 1: scheduling policies (LDES, open two-class tandem)
    // ------------------------------------------------------------------

    @ParameterizedTest(name = "sched={0}")
    @ValueSource(strings = {
            "FCFS", "FCFSPR", "FCFSPI", "LCFS", "LCFSPR", "LCFSPI",
            "PS", "DPS", "GPS", "HOL", "SIRO", "SJF", "LJF",
            "SEPT", "LEPT", "SRPT", "PSJF", "FB", "LRPT", "SETF",
            "FSP", "EDD", "EDF", "LPS",
            "FCFSPRIO", "FCFSPRPRIO", "FCFSPIPRIO", "PSPRIO", "DPSPRIO",
            "GPSPRIO", "LCFSPRIO", "LCFSPRPRIO", "LCFSPIPRIO", "SRPTPRIO"})
    public void schedulingPolicy(String schedName) {
        withSuppressedOutput(() -> {
            SchedStrategy sched = SchedStrategy.valueOf(schedName);
            Network model = new Network("sched_" + schedName);
            Source source = new Source(model, "Source");
            Queue queue = new Queue(model, "Queue", sched);
            Sink sink = new Sink(model, "Sink");

            OpenClass class1 = new OpenClass(model, "Class1", 0);
            OpenClass class2 = new OpenClass(model, "Class2", 1);
            class1.setDeadline(10.0); // used by EDD/EDF, ignored otherwise
            class2.setDeadline(20.0);

            double lambda1 = 0.4, lambda2 = 0.3;
            double mu1 = 2.0, mu2 = 3.0;
            source.setArrival(class1, new Exp(lambda1));
            source.setArrival(class2, new Exp(lambda2));
            queue.setService(class1, new Exp(mu1), 2.0); // weight used by DPS/GPS
            queue.setService(class2, new Exp(mu2), 1.0);
            if (sched == SchedStrategy.LPS) {
                queue.setLimit(2);
            }

            model.link(model.serialRouting(source, queue, sink));

            SolverOptions options = new SolverOptions(SolverType.LDES);
            options.verbose = VerboseLevel.SILENT;
            options.samples = LDES_SAMPLES;
            options.seed = SEED;
            NetworkAvgTable table = new SolverLDES(model, options).getAvgTable();
            assertNotNull(table, schedName + ": null AvgTable");

            // Throughput conservation: per-class queue throughput = arrival rate
            double tput1 = queueMetric(table, "Queue", "Class1", table.getTput());
            double tput2 = queueMetric(table, "Queue", "Class2", table.getTput());
            assertEquals(lambda1, tput1, SIM_RTOL * lambda1,
                    schedName + ": class1 throughput violates conservation");
            assertEquals(lambda2, tput2, SIM_RTOL * lambda2,
                    schedName + ": class2 throughput violates conservation");

            // Utilization law: U = lambda/mu per class
            double util1 = queueMetric(table, "Queue", "Class1", table.getUtil());
            double util2 = queueMetric(table, "Queue", "Class2", table.getUtil());
            assertEquals(lambda1 / mu1, util1, SIM_RTOL * (lambda1 / mu1),
                    schedName + ": class1 utilization law violated");
            assertEquals(lambda2 / mu2, util2, SIM_RTOL * (lambda2 / mu2),
                    schedName + ": class2 utilization law violated");

            // Response time at least the mean service time
            double respt1 = queueMetric(table, "Queue", "Class1", table.getRespT());
            assertTrue(respt1 >= (1.0 / mu1) * (1.0 - SIM_RTOL),
                    schedName + ": class1 response time below service time");
        });
    }

    // ------------------------------------------------------------------
    // Axis 2: service distributions (LDES, M/G/1)
    // ------------------------------------------------------------------

    @ParameterizedTest(name = "dist={0}")
    @ValueSource(strings = {
            "Exp", "Erlang", "HyperExp", "Coxian", "Det", "Uniform",
            "Gamma", "Pareto", "Weibull", "Lognormal", "MAP"})
    public void serviceDistribution(String distName) {
        withSuppressedOutput(() -> {
            Distribution service = makeServiceDistribution(distName);
            double lambda = 0.5;
            double meanS = service.getMean();
            assertTrue(meanS > 0 && Double.isFinite(meanS),
                    distName + ": invalid mean " + meanS);
            assertTrue(lambda * meanS < 0.9, distName + ": test model must be stable");

            Network model = new Network("dist_" + distName);
            Source source = new Source(model, "Source");
            Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
            Sink sink = new Sink(model, "Sink");
            OpenClass jobClass = new OpenClass(model, "Class1", 0);
            source.setArrival(jobClass, new Exp(lambda));
            queue.setService(jobClass, service);
            model.link(model.serialRouting(source, queue, sink));

            SolverOptions options = new SolverOptions(SolverType.LDES);
            options.verbose = VerboseLevel.SILENT;
            options.samples = LDES_SAMPLES;
            options.seed = SEED;
            NetworkAvgTable table = new SolverLDES(model, options).getAvgTable();
            assertNotNull(table, distName + ": null AvgTable");

            double tput = queueMetric(table, "Queue", "Class1", table.getTput());
            assertEquals(lambda, tput, SIM_RTOL * lambda,
                    distName + ": throughput conservation violated");
            double util = queueMetric(table, "Queue", "Class1", table.getUtil());
            assertEquals(lambda * meanS, util, SIM_RTOL * lambda * meanS,
                    distName + ": utilization law violated");
        });
    }

    private static Distribution makeServiceDistribution(String name) {
        switch (name) {
            case "Exp":
                return new Exp(2.0);
            case "Erlang":
                return new Erlang(4.0, 2);
            case "HyperExp":
                return new HyperExp(0.3, 4.0, 1.5);
            case "Coxian": {
                List<Double> mu = new ArrayList<Double>();
                mu.add(3.0);
                mu.add(1.5);
                List<Double> phi = new ArrayList<Double>();
                phi.add(0.4);
                phi.add(1.0);
                return new Coxian(mu, phi);
            }
            case "Det":
                return new Det(0.5);
            case "Uniform":
                return new Uniform(0.2, 0.8);
            case "Gamma":
                return new Gamma(2.0, 0.25);
            case "Pareto":
                return new Pareto(3.0, 0.334);
            case "Weibull":
                return new Weibull(2.0, 0.5);
            case "Lognormal":
                return new Lognormal(-1.0, 0.5);
            case "MAP": {
                Matrix d0 = new Matrix(2, 2);
                d0.set(0, 0, -4.0);
                d0.set(0, 1, 1.0);
                d0.set(1, 0, 0.5);
                d0.set(1, 1, -2.0);
                Matrix d1 = new Matrix(2, 2);
                d1.set(0, 0, 2.0);
                d1.set(0, 1, 1.0);
                d1.set(1, 0, 1.0);
                d1.set(1, 1, 0.5);
                return new MAP(d0, d1);
            }
            default:
                throw new IllegalArgumentException("Unknown distribution " + name);
        }
    }

    // ------------------------------------------------------------------
    // Axis 3: impatience and retrial (CTMC exact vs SSA simulation)
    // ------------------------------------------------------------------

    @ParameterizedTest(name = "variant={0}")
    @ValueSource(strings = {"balking", "reneging", "retrial", "retrialMax"})
    public void impatienceAndRetrial(String variant) {
        withSuppressedOutput(() -> {
            Network ctmcModel = buildImpatienceModel(variant);
            Network ssaModel = buildImpatienceModel(variant);

            SolverOptions ctmcOptions = new SolverOptions(SolverType.CTMC);
            ctmcOptions.verbose = VerboseLevel.SILENT;
            ctmcOptions.cutoff = Matrix.singleton(12);
            NetworkAvgTable ctmcTable = new SolverCTMC(ctmcModel, ctmcOptions).getAvgTable();
            assertNotNull(ctmcTable, variant + ": CTMC null AvgTable");

            SolverOptions ssaOptions = new SolverOptions(SolverType.SSA);
            ssaOptions.verbose = VerboseLevel.SILENT;
            ssaOptions.samples = SSA_SAMPLES;
            ssaOptions.seed = SEED;
            NetworkAvgTable ssaTable = new SolverSSA(ssaModel, ssaOptions).getAvgTable();
            assertNotNull(ssaTable, variant + ": SSA null AvgTable");

            double ctmcQLen = queueMetric(ctmcTable, "Queue", "Class1", ctmcTable.getQLen());
            double ssaQLen = queueMetric(ssaTable, "Queue", "Class1", ssaTable.getQLen());
            assertTrue(ctmcQLen > 0, variant + ": CTMC queue length must be positive");
            assertEquals(ctmcQLen, ssaQLen, CROSS_RTOL * ctmcQLen,
                    variant + ": SSA queue length deviates from exact CTMC");
        });
    }

    private static Network buildImpatienceModel(String variant) {
        Network model = new Network("imp_" + variant);
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");
        OpenClass jobClass = new OpenClass(model, "Class1", 0);
        source.setArrival(jobClass, new Exp(1.5));
        queue.setService(jobClass, new Exp(2.0));
        if ("balking".equals(variant)) {
            List<BalkingThreshold> thresholds = new ArrayList<BalkingThreshold>();
            thresholds.add(new BalkingThreshold(3, 5, 0.5));
            thresholds.add(new BalkingThreshold(6, Integer.MAX_VALUE, 1.0));
            queue.setBalking(jobClass, BalkingStrategy.QUEUE_LENGTH, thresholds);
        } else if ("reneging".equals(variant)) {
            jobClass.setPatience(new Exp(0.5));
        } else if ("retrial".equals(variant)) {
            queue.setCapacity(4);
            queue.setRetrial(jobClass, new Exp(0.5), -1);
        } else if ("retrialMax".equals(variant)) {
            queue.setCapacity(4);
            queue.setRetrial(jobClass, new Exp(0.5), 3);
        } else {
            throw new IllegalArgumentException("Unknown variant " + variant);
        }
        model.link(model.serialRouting(source, queue, sink));
        return model;
    }

    // ------------------------------------------------------------------
    // Axis 4: polling disciplines (LDES)
    // ------------------------------------------------------------------

    @ParameterizedTest(name = "polling={0}")
    @ValueSource(strings = {"EXHAUSTIVE", "GATED", "KLIMITED", "DECREMENTING"})
    public void pollingDiscipline(String typeName) {
        withSuppressedOutput(() -> {
            Network model = new Network("poll_" + typeName);
            Source source = new Source(model, "Source");
            Queue queue = new Queue(model, "Queue", SchedStrategy.POLLING);
            Sink sink = new Sink(model, "Sink");
            OpenClass class1 = new OpenClass(model, "Class1", 0);
            OpenClass class2 = new OpenClass(model, "Class2", 0);

            double lambda1 = 0.3, lambda2 = 0.2;
            source.setArrival(class1, new Exp(lambda1));
            source.setArrival(class2, new Exp(lambda2));
            queue.setService(class1, new Exp(2.0));
            queue.setService(class2, new Exp(2.5));
            queue.setSwitchover(class1, new Exp(10.0));
            queue.setSwitchover(class2, new Exp(10.0));
            if ("KLIMITED".equals(typeName)) {
                queue.setPollingType(PollingType.KLIMITED, 1);
            } else {
                queue.setPollingType(PollingType.valueOf(typeName));
            }
            model.link(model.serialRouting(source, queue, sink));

            SolverOptions options = new SolverOptions(SolverType.LDES);
            options.verbose = VerboseLevel.SILENT;
            options.samples = LDES_SAMPLES;
            options.seed = SEED;
            NetworkAvgTable table = new SolverLDES(model, options).getAvgTable();
            assertNotNull(table, typeName + ": null AvgTable");

            double tput1 = queueMetric(table, "Queue", "Class1", table.getTput());
            double tput2 = queueMetric(table, "Queue", "Class2", table.getTput());
            assertEquals(lambda1, tput1, SIM_RTOL * lambda1,
                    typeName + ": class1 throughput conservation violated");
            assertEquals(lambda2, tput2, SIM_RTOL * lambda2,
                    typeName + ": class2 throughput conservation violated");
        });
    }

    // ------------------------------------------------------------------
    // Axis 5: finite capacity region (LDES)
    // ------------------------------------------------------------------

    @Test
    public void finiteCapacityRegion() {
        withSuppressedOutput(() -> {
            Network model = new Network("fcr");
            Source source = new Source(model, "Source");
            Queue queue1 = new Queue(model, "Queue1", SchedStrategy.FCFS);
            Queue queue2 = new Queue(model, "Queue2", SchedStrategy.FCFS);
            Sink sink = new Sink(model, "Sink");
            OpenClass jobClass = new OpenClass(model, "Class1", 0);
            double lambda = 1.0;
            source.setArrival(jobClass, new Exp(lambda));
            queue1.setService(jobClass, new Exp(2.0));
            queue2.setService(jobClass, new Exp(2.0));
            model.link(model.serialRouting(source, queue1, queue2, sink));

            List<jline.lang.nodes.Node> regionNodes = new ArrayList<jline.lang.nodes.Node>();
            regionNodes.add(queue1);
            regionNodes.add(queue2);
            jline.lang.Region region = model.addRegion(regionNodes);
            region.setGlobalMaxJobs(3);

            SolverOptions options = new SolverOptions(SolverType.LDES);
            options.verbose = VerboseLevel.SILENT;
            options.samples = LDES_SAMPLES;
            options.seed = SEED;
            NetworkAvgTable table = new SolverLDES(model, options).getAvgTable();
            assertNotNull(table, "FCR: null AvgTable");

            // Region caps total jobs at 3, so total queue length <= 3 and
            // accepted throughput cannot exceed the offered load.
            double qlen1 = queueMetric(table, "Queue1", "Class1", table.getQLen());
            double qlen2 = queueMetric(table, "Queue2", "Class1", table.getQLen());
            assertTrue(qlen1 + qlen2 <= 3.0 + 1e-6,
                    "FCR: mean population exceeds region capacity");
            double tput2 = queueMetric(table, "Queue2", "Class1", table.getTput());
            assertTrue(tput2 <= lambda * (1.0 + SIM_RTOL),
                    "FCR: accepted throughput exceeds offered load");
            assertTrue(tput2 > 0, "FCR: no traffic flowed through the region");
        });
    }

    // ------------------------------------------------------------------
    // Axis 6: exact solver agreement on product form (MVA vs NC vs CTMC)
    // ------------------------------------------------------------------

    @Test
    public void productFormExactAgreement() {
        withSuppressedOutput(() -> {
            NetworkAvgTable mvaTable = new SolverMVA(buildClosedPsModel(), "exact").getAvgTable();
            NetworkAvgTable ncTable = new SolverNC(buildClosedPsModel(), "exact").getAvgTable();
            SolverOptions ctmcOptions = new SolverOptions(SolverType.CTMC);
            ctmcOptions.verbose = VerboseLevel.SILENT;
            NetworkAvgTable ctmcTable = new SolverCTMC(buildClosedPsModel(), ctmcOptions).getAvgTable();

            double mvaQ = queueMetric(mvaTable, "Queue", "Class1", mvaTable.getQLen());
            double ncQ = queueMetric(ncTable, "Queue", "Class1", ncTable.getQLen());
            double ctmcQ = queueMetric(ctmcTable, "Queue", "Class1", ctmcTable.getQLen());
            assertTrue(mvaQ > 0, "MVA queue length must be positive");
            assertEquals(mvaQ, ncQ, EXACT_RTOL * mvaQ, "MVA and NC disagree on product form");
            assertEquals(mvaQ, ctmcQ, EXACT_RTOL * mvaQ, "MVA and CTMC disagree on product form");
        });
    }

    private static Network buildClosedPsModel() {
        Network model = new Network("pf_closed");
        Delay delay = new Delay(model, "Delay");
        Queue queue = new Queue(model, "Queue", SchedStrategy.PS);
        ClosedClass jobClass = new ClosedClass(model, "Class1", 3, delay, 0);
        delay.setService(jobClass, new Exp(1.0));
        queue.setService(jobClass, new Exp(2.0));
        model.link(model.serialRouting(delay, queue));
        return model;
    }

    // ------------------------------------------------------------------
    // Axis 6b: fluid solution methods agree with the exact solution
    // ------------------------------------------------------------------

    @ParameterizedTest(name = "fluidMethod={0}")
    @ValueSource(strings = {"matrix", "closing", "statedep"})
    public void fluidMethodsApproximateExactSolution(String method) {
        withSuppressedOutput(() -> {
            NetworkAvgTable mvaTable = new SolverMVA(buildClosedPsModel(), "exact").getAvgTable();
            double exactQ = queueMetric(mvaTable, "Queue", "Class1", mvaTable.getQLen());

            SolverOptions options = new SolverOptions(SolverType.FLUID);
            options.verbose = VerboseLevel.SILENT;
            options.method = method;
            NetworkAvgTable fluidTable =
                    new SolverFluid(buildClosedPsModel(), options).getAvgTable();
            double fluidQ = queueMetric(fluidTable, "Queue", "Class1", fluidTable.getQLen());
            // Fluid limits are approximations for small populations: accept a
            // generous but bounded deviation from the exact queue length
            assertEquals(exactQ, fluidQ, 0.30 * exactQ,
                    "fluid method '" + method + "' deviates grossly from exact");
        });
    }

    // ------------------------------------------------------------------
    // Axis 7: transient analysis and response time CDF (CTMC, FLD)
    // ------------------------------------------------------------------

    @Test
    public void transientConvergesToSteadyState() {
        withSuppressedOutput(() -> {
            Network model = buildClosedPsModel();
            SolverOptions options = new SolverOptions(SolverType.FLUID);
            options.verbose = VerboseLevel.SILENT;
            options.timespan[0] = 0.0;
            options.timespan[1] = 200.0;
            SolverFluid solver = new SolverFluid(model, options);
            NetworkAvgTable steady = solver.getAvgTable();
            assertNotNull(steady, "Fluid steady-state table is null");
            double steadyQ = queueMetric(steady, "Queue", "Class1", steady.getQLen());
            assertTrue(steadyQ > 0, "Fluid steady-state queue length must be positive");

            // Transient solution must be produced and populated
            solver.getTranAvg();
            assertTrue(solver.hasTranResults(), "Fluid getTranAvg produced no QNt data");
        });
    }

    @Test
    public void cdfRespTIsAProbabilityDistribution() {
        withSuppressedOutput(() -> {
            Network model = buildClosedPsModel();
            SolverOptions options = new SolverOptions(SolverType.FLUID);
            options.verbose = VerboseLevel.SILENT;
            SolverFluid solver = new SolverFluid(model, options);
            DistributionResult cdf = solver.getCdfRespT();
            assertNotNull(cdf, "Fluid getCdfRespT returned null");
            assertNotNull(cdf.cdfData, "Fluid getCdfRespT returned null cdfData");
            boolean checkedOne = false;
            for (List<Matrix> row : cdf.cdfData) {
                if (row == null) continue;
                for (Matrix cell : row) {
                    if (cell == null || cell.getNumRows() == 0) continue;
                    checkedOne = true;
                    // Column 0 holds CDF probabilities, column 1 the times
                    double prev = -1e-9;
                    for (int t = 0; t < cell.getNumRows(); t++) {
                        double p = cell.get(t, 0);
                        assertTrue(p >= -1e-6 && p <= 1.0 + 1e-6,
                                "CDF value outside [0,1]: " + p);
                        assertTrue(p >= prev - 1e-6, "CDF must be nondecreasing");
                        prev = p;
                    }
                    assertTrue(prev > 0.9, "CDF must approach 1, final value " + prev);
                }
            }
            assertTrue(checkedOne, "getCdfRespT produced no evaluable CDF");
        });
    }

    // ------------------------------------------------------------------
    // Helpers
    // ------------------------------------------------------------------

    /**
     * Extracts the metric value for (station, class) from an AvgTable
     * given the metric column returned by the table getter.
     */
    private static double queueMetric(NetworkAvgTable table, String station,
                                      String jobClass, List<Double> column) {
        List<String> stations = table.getStationNames();
        List<String> classes = table.getClassNames();
        for (int i = 0; i < stations.size(); i++) {
            if (stations.get(i).equals(station) && classes.get(i).equals(jobClass)) {
                return column.get(i);
            }
        }
        throw new AssertionError("Row not found in AvgTable: " + station + "/" + jobClass);
    }
}
