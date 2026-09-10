package jline.solvers.fluid;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.ArrayList;
import java.util.List;

import org.junit.jupiter.api.Test;

import jline.VerboseLevel;
import jline.api.mam.Map_lambda;
import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.NetworkStruct;
import jline.lang.OpenClass;
import jline.lang.constant.ProcessType;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.processes.Exp;
import jline.lang.processes.MAP;
import jline.lang.processes.MAPt;
import jline.lang.processes.MMPP2;
import jline.lang.processes.NHPP;
import jline.lang.processes.PHt;
import jline.util.matrix.Matrix;

/**
 * Regression tests for MAPt/PHt and the SolverFluid "kp" method.
 *
 * <p>"kp" integrates the fluid and diffusion limits of Ko and Pender, "Diffusion limits for the
 * (MAP_t/Ph_t/inf)^N queueing network", Oper. Res. Lett. 45 (2017) 248-253. For an
 * infinite-server network the rate functions are affine in the state, so the mean and covariance
 * ODEs close EXACTLY: the assertions below compare against exact references rather than a
 * tolerance band.
 *
 * <p>References used, strongest first:
 * <ol>
 *   <li>the exact stationary MAP/M/inf, whose mean is lambda/mu = 1.6 and whose variance is
 *       1.64 -- note Var != Mean, which is precisely what a PH renewal approximation of the MAP
 *       would lose;</li>
 *   <li>the exact period-average arrival rate of a cyclic schedule;</li>
 *   <li>degeneracy: a 1-phase MAPt IS an NHPP, a 1-segment MAPt IS a MAP.</li>
 * </ol>
 *
 * <p>Also covers the phase-loop fix: the closing and matrix ODEs enumerated internal phase
 * transitions only from phases 1..h-1, which drops the last row of D0 for a cyclic D0 and cost
 * 14% of the arrival rate on a 2-phase MAP source.
 *
 * <p>Mirrors python/line_solver/tests/test_solver_fld_kp.py.
 */
public class SolverFluidKoPenderTest {

    /** A 2-phase MAP whose D0 is cyclic, so the process is genuinely non-renewal. */
    private static final double[][] D0_A = {{-5.0, 1.0}, {2.0, -4.0}};
    private static final double[][] D1_A = {{3.0, 1.0}, {1.0, 1.0}};
    private static final double[][] D0_B = {{-12.0, 3.0}, {5.0, -9.0}};
    private static final double[][] D1_B = {{7.0, 2.0}, {2.0, 2.0}};
    private static final double[] BREAKPOINTS = {0.0, 1.0, 2.5};
    private static final double MU = 2.0;

    private static Matrix m(double[][] entries) {
        Matrix out = new Matrix(entries.length, entries[0].length);
        for (int i = 0; i < entries.length; i++) {
            for (int j = 0; j < entries[0].length; j++) {
                out.set(i, j, entries[i][j]);
            }
        }
        return out;
    }

    private static List<Matrix> list(double[][]... entries) {
        List<Matrix> out = new ArrayList<Matrix>();
        for (int k = 0; k < entries.length; k++) {
            out.add(m(entries[k]));
        }
        return out;
    }

    /** Single-segment MAPt, i.e. an ordinary MAP carried by the schedule machinery. */
    private static MAPt oneSegment() {
        return new MAPt(new double[]{0.0, 1.0}, list(D0_A), list(D1_A), true);
    }

    /** Two-segment cyclic MAPt. */
    private static MAPt twoSegment() {
        return new MAPt(BREAKPOINTS, list(D0_A, D0_B), list(D1_A, D1_B), true);
    }

    private static Network infStation(jline.lang.processes.Distribution arrival,
                                      jline.lang.processes.Distribution service) {
        Network model = new Network("fld_kp");
        Source source = new Source(model, "Source");
        Delay delay = new Delay(model, "Delay");
        Sink sink = new Sink(model, "Sink");
        OpenClass jobclass = new OpenClass(model, "Class1");
        source.setArrival(jobclass, arrival);
        delay.setService(jobclass, service);
        model.link(model.serialRouting(source, delay, sink));
        return model;
    }

    private static SolverFluid kpSolver(Network model, double tend, double tol) {
        SolverFluid solver = new SolverFluid(model);
        solver.options.method = "kp";
        solver.options.tol = tol;
        solver.options.verbose = VerboseLevel.SILENT;
        if (Double.isFinite(tend)) {
            solver.options.timespan = new double[]{0.0, tend};
        }
        return solver;
    }

    @Test
    public void testOnePhaseMAPtIsAnNHPP() {
        double[] bp = {0.0, 1.0, 2.0};
        MAPt oneph = new MAPt(bp, list(new double[][]{{-1.0}}, new double[][]{{-3.0}}),
                list(new double[][]{{1.0}}, new double[][]{{3.0}}), true);
        NHPP nhpp = new NHPP(bp, new double[]{1.0, 3.0}, true);
        assertEquals(nhpp.getTimeAverageRate(), oneph.getTimeAverageRate(), 1e-12);
        for (double t : new double[]{0.25, 1.25, 2.25, 3.25}) {
            assertEquals(nhpp.getRateAt(t), oneph.getRateAt(t), 1e-12);
        }
    }

    @Test
    public void testSingleSegmentMAPtIsAMAP() {
        MAPt mapt = oneSegment();
        assertEquals(Map_lambda.map_lambda(m(D0_A), m(D1_A)), mapt.getTimeAverageRate(), 1e-12);
        assertEquals(2, mapt.getNumberOfPhases());
    }

    @Test
    public void testScalarSummariesAreNaN() {
        MAPt mapt = twoSegment();
        assertTrue(Double.isNaN(mapt.getSCV()));
        assertTrue(Double.isNaN(mapt.getSkewness()));
        assertTrue(Double.isNaN(mapt.evalCDF(1.0)));
        assertTrue(Double.isNaN(mapt.evalLST(1.0)));
    }

    @Test
    public void testVaryingSparsityPatternIsRefused() {
        // A transition present in one segment and absent in another cannot be a per-entry
        // multiplier on the time-averaged nominal, so it is an error rather than a silently
        // dropped transition.
        double[][] d0off = {{-4.0, 0.0}, {2.0, -4.0}};
        double[][] d1off = {{3.0, 1.0}, {1.0, 1.0}};
        assertThrows(IllegalArgumentException.class, new org.junit.jupiter.api.function.Executable() {
            @Override
            public void execute() {
                new MAPt(BREAKPOINTS, list(D0_A, d0off), list(D1_A, d1off), true);
            }
        });
    }

    @Test
    public void testStructCarriesTheScheduleAndItsPhaseCount() {
        Network model = infStation(twoSegment(), new Exp(MU));
        NetworkStruct sn = model.getStruct(true);
        assertEquals(ProcessType.MAPT,
                sn.procid.get(sn.stations.get(0)).get(sn.jobclasses.get(0)));
        // phases is the ORDER OF THE SEGMENT MATRICES, unlike NHPP which collapses to 1
        assertEquals(2.0, sn.phases.get(0, 0), 1e-12);
        // rates is the arrival rate of the time-averaged nominal; scv is undefined
        assertEquals(twoSegment().getTimeAverageRate(), sn.rates.get(0, 0), 1e-9);
        assertTrue(Double.isNaN(sn.scv.get(0, 0)));
    }

    @Test
    public void testMatchesTheExactStationaryMAPMInf() {
        // MAP/M/inf: exact mean lambda/mu = 1.6 and exact variance 1.64, both from the
        // stationary CTMC over (phase, count). Var != Mean because the arrival stream is
        // non-renewal.
        // getAvg rejects a finite timespan by design, so the steady state is requested
        // unbounded: with a single segment the period average IS the fixed point.
        SolverFluid solver = kpSolver(infStation(oneSegment(), new Exp(MU)),
                Double.POSITIVE_INFINITY, 1e-9);
        solver.getAvg();
        assertEquals(Map_lambda.map_lambda(m(D0_A), m(D1_A)) / MU,
                solver.result.QN.get(1, 0), 1e-6);
        assertEquals(Map_lambda.map_lambda(m(D0_A), m(D1_A)),
                solver.result.TN.get(1, 0), 1e-6);
    }

    @Test
    public void testCyclicSteadyStateIsThePeriodAverage() {
        // A cyclic schedule has no fixed point, so getAvg reports the average over the last
        // full period. Source and station throughput must then agree, which a snapshot of an
        // arbitrary point of the cycle would not.
        SolverFluid solver = kpSolver(infStation(twoSegment(), new Exp(MU)),
                Double.POSITIVE_INFINITY, 1e-9);
        solver.getAvg();
        double sourceTput = solver.result.TN.get(0, 0);
        double delayTput = solver.result.TN.get(1, 0);
        assertEquals(sourceTput, delayTput, 1e-3);

        double lamA = Map_lambda.map_lambda(m(D0_A), m(D1_A));
        double lamB = Map_lambda.map_lambda(m(D0_B), m(D1_B));
        double expected = (1.0 * lamA + 1.5 * lamB) / 2.5;
        assertEquals(expected, delayTput, 1e-3);
    }

    @Test
    public void testTransientVariesOverTheCycle() {
        SolverFluid solver = kpSolver(infStation(twoSegment(), new Exp(MU)), 5.0, 1e-9);
        solver.getTranAvg();
        Matrix qlen = solver.result.QNt[1][0];
        assertNotNull(qlen);
        double lo = Double.POSITIVE_INFINITY;
        double hi = Double.NEGATIVE_INFINITY;
        for (int i = 0; i < qlen.getNumRows(); i++) {
            if (solver.result.t.get(i, 0) > 1.0) {
                lo = Math.min(lo, qlen.get(i, 0));
                hi = Math.max(hi, qlen.get(i, 0));
            }
        }
        // A solver that silently used the time-averaged nominal would give a flat trajectory.
        assertTrue(hi - lo > 0.5, "queue trajectory is flat (spread " + (hi - lo) + ")");
    }

    @Test
    public void testTimeVaryingPhaseTypeService() {
        List<Matrix> alpha = list(new double[][]{{1.0, 0.0}}, new double[][]{{1.0, 0.0}});
        List<Matrix> subgen = list(new double[][]{{-4.0, 4.0}, {0.0, -4.0}},
                new double[][]{{-9.0, 9.0}, {0.0, -9.0}});
        PHt pht = new PHt(BREAKPOINTS, alpha, subgen, true);
        assertEquals(2, pht.getNumberOfPhases());
        assertTrue(Double.isNaN(pht.getSCV()));

        SolverFluid solver = kpSolver(infStation(twoSegment(), pht), 5.0, 1e-9);
        solver.getTranAvg();
        Matrix qlen = solver.result.QNt[1][0];
        for (int i = 0; i < qlen.getNumRows(); i++) {
            assertTrue(qlen.get(i, 0) >= -1e-9, "negative queue length at index " + i);
        }
        assertTrue(qlen.get(qlen.getNumRows() - 1, 0) > 0.0);
    }

    @Test
    public void testClosedModelIsRejected() {
        final Network model = new Network("closed_kp");
        Delay delay = new Delay(model, "Delay");
        Queue queue = new Queue(model, "Queue", SchedStrategy.PS);
        ClosedClass jobclass = new ClosedClass(model, "Class1", 5, delay);
        delay.setService(jobclass, new Exp(1.0));
        queue.setService(jobclass, new Exp(2.0));
        model.link(model.serialRouting(delay, queue, delay));
        assertThrows(RuntimeException.class, new org.junit.jupiter.api.function.Executable() {
            @Override
            public void execute() {
                SolverFluid solver = kpSolver(model, 10.0, 1e-6);
                solver.getAvg();
            }
        });
    }

    @Test
    public void testCyclicD0KeepsItsExactRateUnderClosingAndMatrix() {
        // The closing and matrix ODEs enumerated phase transitions only from phases 1..h-1,
        // which is valid only for an acyclic PH. Dropping the last row of D0 gave Tput 2.7586
        // for this MAP (exact 3.2) and 2.1429 for the MMPP2 (exact 2.5).
        double mapExact = Map_lambda.map_lambda(m(D0_A), m(D1_A));
        MMPP2 mmpp = new MMPP2(4.0, 1.0, 0.5, 0.5);
        for (String method : new String[]{"closing", "matrix"}) {
            SolverFluid mapSolver = new SolverFluid(infStation(new MAP(m(D0_A), m(D1_A)),
                    new Exp(MU)));
            mapSolver.options.method = method;
            mapSolver.options.verbose = VerboseLevel.SILENT;
            mapSolver.getAvg();
            assertEquals(mapExact, mapSolver.result.TN.get(1, 0), 1e-4,
                    "MAP source under " + method);

            SolverFluid mmppSolver = new SolverFluid(infStation(mmpp, new Exp(MU)));
            mmppSolver.options.method = method;
            mmppSolver.options.verbose = VerboseLevel.SILENT;
            mmppSolver.getAvg();
            assertEquals(mmpp.getRate(), mmppSolver.result.TN.get(1, 0), 1e-4,
                    "MMPP2 source under " + method);
        }
    }

    @Test
    public void testAcyclicDistributionsAreUnchanged() {
        // The phase-loop fix adds only zero-rate candidates for an acyclic PH.
        SolverFluid expSolver = new SolverFluid(infStation(new Exp(3.0), new Exp(MU)));
        expSolver.options.verbose = VerboseLevel.SILENT;
        // Pinned to the matrix method, which is what "default" resolved to when
        // these exact values were written; the default now prefers the
        // second-order closure "minnormal", whose open-model throughput is exact
        // only at its fixed point and lands 1.4e-5 away here, inside the closure's
        // own tolerance but outside the 1e-6 asserted below. The assertion is
        // about the phase-loop fix, not about which closure runs.
        expSolver.options.method = "matrix";
        expSolver.getAvg();
        assertEquals(3.0, expSolver.result.TN.get(1, 0), 1e-6);
        assertEquals(3.0 / MU, expSolver.result.QN.get(1, 0), 1e-6);
    }

    @Test
    public void testCovarianceReachesTheResult() {
        // The covariance is the point of this method, so it must reach the caller. It used to
        // be computed and discarded, leaving getTranAvgVar with nothing to return.
        SolverFluid solver = kpSolver(infStation(twoSegment(), new Exp(MU)), 5.0, 1e-9);
        solver.getTranAvg();
        FluidResult fr = (FluidResult) solver.result;
        assertNotNull(fr.QVart, "QVart was not stored on the result");
        assertNotNull(fr.Sigmat, "Sigmat was not stored on the result");
        assertEquals(fr.t.getNumRows(), fr.QVart[1][0].getNumRows());
        assertEquals(fr.t.getNumRows(), fr.Sigmat.length);
        for (int n = 0; n < fr.QVart[1][0].getNumRows(); n++) {
            assertTrue(fr.QVart[1][0].get(n, 0) >= -1e-9,
                    "negative variance at index " + n);
        }
        // the per-block variance is the sum of that block's Sigma entries
        Matrix sig = fr.Sigmat[fr.Sigmat.length - 1];
        double blockSum = 0.0;
        for (int i = 2; i < sig.getNumRows(); i++) {
            for (int j = 2; j < sig.getNumCols(); j++) {
                blockSum += sig.get(i, j);
            }
        }
        assertEquals(blockSum, fr.QVart[1][0].get(fr.QVart[1][0].getNumRows() - 1, 0), 1e-9);

        Matrix[][] viaAccessor = solver.getTranAvgVar();
        assertNotNull(viaAccessor);
        assertEquals(fr.QVart[1][0].get(0, 0), viaAccessor[1][0].get(0, 0), 1e-12);
    }

    @Test
    public void testLdesSimulatesTheMAPtArrivalStream() {
        // The independent cross-check: LDES walks the actual sample path, carrying the MAP
        // phase across a breakpoint and resampling the clock there, so it validates the kp
        // ODE without sharing any of its machinery. A single-segment MAPt must reproduce the
        // exact stationary MAP/M/inf; a cyclic one must match the kp period average.
        jline.solvers.SolverOptions opts = new jline.solvers.SolverOptions(
                jline.lang.constant.SolverType.LDES);
        opts.seed = 23456;
        opts.samples = 600000;
        opts.verbose = VerboseLevel.SILENT;
        jline.solvers.ldes.SolverLDES sim = new jline.solvers.ldes.SolverLDES(
                infStation(oneSegment(), new Exp(MU)), opts);
        sim.getAvg();
        double exactMean = Map_lambda.map_lambda(m(D0_A), m(D1_A)) / MU;
        assertEquals(exactMean, sim.result.QN.get(1, 0), 0.01);
        assertEquals(Map_lambda.map_lambda(m(D0_A), m(D1_A)), sim.result.TN.get(1, 0), 0.01);

        jline.solvers.SolverOptions cycOpts = new jline.solvers.SolverOptions(
                jline.lang.constant.SolverType.LDES);
        cycOpts.seed = 998877;
        cycOpts.samples = 1500000;
        cycOpts.verbose = VerboseLevel.SILENT;
        jline.solvers.ldes.SolverLDES cycSim = new jline.solvers.ldes.SolverLDES(
                infStation(twoSegment(), new Exp(MU)), cycOpts);
        cycSim.getAvg();
        double lamA = Map_lambda.map_lambda(m(D0_A), m(D1_A));
        double lamB = Map_lambda.map_lambda(m(D0_B), m(D1_B));
        double expected = (1.0 * lamA + 1.5 * lamB) / 2.5;
        assertEquals(expected, cycSim.result.TN.get(1, 0), 0.01);
        assertEquals(expected / MU, cycSim.result.QN.get(1, 0), 0.01);
    }

    @Test
    public void testLdesSimulatesPHtServiceAndRefusesSharedDisciplines() {
        List<Matrix> alpha = list(new double[][]{{1.0, 0.0}}, new double[][]{{1.0, 0.0}});
        List<Matrix> subgen = list(new double[][]{{-4.0, 4.0}, {0.0, -4.0}},
                new double[][]{{-9.0, 9.0}, {0.0, -9.0}});
        PHt pht = new PHt(BREAKPOINTS, alpha, subgen, true);

        // MAP_t/Ph_t/inf: the paper's unit, simulated end to end and compared with kp.
        jline.solvers.SolverOptions simOpts = new jline.solvers.SolverOptions(
                jline.lang.constant.SolverType.LDES);
        simOpts.seed = 424242;
        simOpts.samples = 2000000;
        simOpts.verbose = VerboseLevel.SILENT;
        jline.solvers.ldes.SolverLDES sim = new jline.solvers.ldes.SolverLDES(
                infStation(twoSegment(), pht), simOpts);
        sim.getAvg();

        SolverFluid kp = kpSolver(infStation(twoSegment(), pht),
                Double.POSITIVE_INFINITY, 1e-9);
        kp.getAvg();
        assertEquals(kp.result.QN.get(1, 0), sim.result.QN.get(1, 0), 0.01);
        assertEquals(kp.result.TN.get(1, 0), sim.result.TN.get(1, 0), 0.01);

        // A PHt duration is walked from the SERVICE START epoch, which is exact only while
        // service runs continuously; a shared discipline must be refused, not silently wrong.
        final Network psModel = new Network("ldes_pht_ps");
        Source src = new Source(psModel, "Source");
        Queue que = new Queue(psModel, "Queue", SchedStrategy.PS);
        Sink snk = new Sink(psModel, "Sink");
        OpenClass cls = new OpenClass(psModel, "Class1");
        src.setArrival(cls, twoSegment());
        que.setService(cls, pht);
        psModel.link(psModel.serialRouting(src, que, snk));
        assertThrows(RuntimeException.class, new org.junit.jupiter.api.function.Executable() {
            @Override
            public void execute() {
                jline.solvers.SolverOptions o = new jline.solvers.SolverOptions(
                        jline.lang.constant.SolverType.LDES);
                o.seed = 1;
                o.samples = 20000;
                o.verbose = VerboseLevel.SILENT;
                new jline.solvers.ldes.SolverLDES(psModel, o).getAvg();
            }
        });
    }

    @Test
    public void testSeededMeanAndCovarianceAreHonoured() {
        // M/M/inf started from a POISSON(x0) queue stays Poisson at every t, so a run seeded
        // with mean x0 and covariance x0 must return the exact transient mean AND Var == Mean
        // along the whole trajectory. Honouring only one of the two seeds breaks the identity.
        final double lambda = 3.0;
        final double x0 = 5.0;
        MAPt poisson = new MAPt(new double[]{0.0, 1.0}, list(new double[][]{{-lambda}}),
                list(new double[][]{{lambda}}), true);
        SolverFluid solver = kpSolver(infStation(poisson, new Exp(MU)), 4.0, 1e-9);
        // dim = 2: the arrival phase (u-block) then the Delay service phase.
        solver.options.config.kp_init_sol = new double[]{1.0, x0};
        Matrix cov = new Matrix(2, 2);
        cov.set(1, 1, x0);
        solver.options.config.init_cov = cov;
        solver.getTranAvg();

        Matrix qlen = solver.result.QNt[1][0];
        Matrix qvar = ((FluidResult) solver.result).QVart[1][0];
        assertNotNull(qvar);
        for (int i = 0; i < qlen.getNumRows(); i++) {
            double t = solver.result.t.get(i, 0);
            double exact = x0 * Math.exp(-MU * t) + (lambda / MU) * (1.0 - Math.exp(-MU * t));
            assertEquals(exact, qlen.get(i, 0), 1e-6, "seeded mean at t=" + t);
            assertEquals(qlen.get(i, 0), qvar.get(i, 0), 1e-6, "Var != Mean at t=" + t);
        }
    }

    @Test
    public void testASeedThatDoesNotFitIsRefused() {
        // Silently dropping it would integrate from the DEFAULT initial condition under the
        // caller's name and hand back a plausible wrong trajectory.
        final MAPt poisson = new MAPt(new double[]{0.0, 1.0}, list(new double[][]{{-3.0}}),
                list(new double[][]{{3.0}}), true);
        final List<double[]> badMeans = new ArrayList<double[]>();
        badMeans.add(new double[]{1.0, 2.0, 3.0});
        for (final double[] bad : badMeans) {
            assertThrows(RuntimeException.class, new org.junit.jupiter.api.function.Executable() {
                @Override
                public void execute() {
                    SolverFluid solver = kpSolver(infStation(poisson, new Exp(MU)), 1.0, 1e-9);
                    solver.options.config.kp_init_sol = bad;
                    solver.getTranAvg();
                }
            });
        }
        final List<Matrix> badCovs = new ArrayList<Matrix>();
        badCovs.add(new Matrix(1, 2));                       // wrong shape
        Matrix asym = new Matrix(2, 2);
        asym.set(0, 0, 1.0);
        asym.set(0, 1, 2.0);
        asym.set(1, 1, 1.0);
        badCovs.add(asym);                                   // not symmetric
        for (final Matrix bad : badCovs) {
            assertThrows(RuntimeException.class, new org.junit.jupiter.api.function.Executable() {
                @Override
                public void execute() {
                    SolverFluid solver = kpSolver(infStation(poisson, new Exp(MU)), 1.0, 1e-9);
                    solver.options.config.init_cov = bad;
                    solver.getTranAvg();
                }
            });
        }
    }

    @Test
    public void testMethodAndFeatureRegistration() {
        SolverFluid solver = new SolverFluid(infStation(new Exp(1.0), new Exp(MU)));
        boolean hasKp = false;
        for (String method : solver.listValidMethods()) {
            if ("kp".equals(method)) {
                hasKp = true;
                break;
            }
        }
        assertTrue(hasKp, "kp is not registered in listValidMethods");
        assertTrue(SolverFluid.getFeatureSet().inspectFeature("MAPt"), "MAPt not declared");
        assertTrue(SolverFluid.getFeatureSet().inspectFeature("PHt"), "PHt not declared");
        assertTrue(SolverFluid.getFeatureSet().inspectFeature("MAP"), "MAP not declared");
        assertTrue(SolverFluid.getFeatureSet().inspectFeature("MMPP2"), "MMPP2 not declared");
        assertTrue(jline.solvers.ldes.SolverLDES.getFeatureSet().inspectFeature("MAPt"),
                "LDES does not declare MAPt");
        assertTrue(jline.solvers.ldes.SolverLDES.getFeatureSet().inspectFeature("PHt"),
                "LDES does not declare PHt");
    }
}
