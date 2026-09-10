/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.env;

import jline.io.MAPQN2RENV;
import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.processes.Exp;
import jline.lang.processes.MAP;
import jline.lang.processes.MMPP2;
import jline.solvers.SolverResult;
import jline.solvers.ctmc.SolverCTMC;
import jline.solvers.mva.SolverMVA;
import jline.solvers.nc.SolverNC;
import jline.util.matrix.Matrix;
import jline.VerboseLevel;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Verifies the MAP/MMPP random-environment fallback of NetworkSolver: a solver
 * that cannot consume a non-renewal process solves the model through its
 * environment image instead of rejecting it.
 *
 * <p>Targets are the MATLAB reference values of the same models
 * (matlab @NetworkSolver/mapEnvApprox.m).
 */
public class MapEnvFallbackTest {

    /** Closed network, MMPP2 service: Delay(Exp 1) -> Queue(MMPP2 1,10,.2,.3), N=5. */
    private static Network mmppClosed() {
        Network model = new Network("mmppClosed");
        Delay delay = new Delay(model, "Think");
        Queue queue = new Queue(model, "Q1", SchedStrategy.FCFS);
        ClosedClass jobclass = new ClosedClass(model, "C1", 5, delay);
        delay.setService(jobclass, new Exp(1.0));
        queue.setService(jobclass, new MMPP2(1.0, 10.0, 0.2, 0.3));
        model.link(model.serialRouting(delay, queue));
        return model;
    }

    /** Closed network, general MAP service (non-diagonal D1), N=4. */
    private static Network genMapClosed() {
        Network model = new Network("genMap");
        Delay delay = new Delay(model, "Think");
        Queue queue = new Queue(model, "Q1", SchedStrategy.FCFS);
        ClosedClass jobclass = new ClosedClass(model, "C1", 4, delay);
        delay.setService(jobclass, new Exp(1.0));
        Matrix D0 = new Matrix(2, 2);
        D0.set(0, 0, -3.0); D0.set(0, 1, 0.5); D0.set(1, 0, 0.2); D0.set(1, 1, -2.0);
        Matrix D1 = new Matrix(2, 2);
        D1.set(0, 0, 2.0); D1.set(0, 1, 0.5); D1.set(1, 0, 1.0); D1.set(1, 1, 0.8);
        queue.setService(jobclass, new MAP(D0, D1));
        model.link(model.serialRouting(delay, queue));
        return model;
    }

    @Test
    public void mmppServiceIsInterceptedByMVA() {
        Network model = mmppClosed();
        SolverMVA solver = new SolverMVA(model);
        solver.options.verbose = VerboseLevel.SILENT;
        assertTrue(solver.needsMapEnv(solver.options),
                "MVA has no MAP support on the resolved method, so the fallback must fire");
        SolverResult res = solver.getAvg();
        // MATLAB reference: the 'auto' timescale test picks the rate-averaged
        // limit here, X = 3.4427 (exact CTMC 3.0952).
        assertTrue(solver.result.method.contains("env.avg"), "reported method was " + solver.result.method);
        assertEquals(3.4427, res.TN.get(0, 0), 1e-3);
        assertEquals(1.5573, res.QN.get(1, 0), 1e-3);
    }

    @Test
    public void mmppServiceIsInterceptedByNC() {
        Network model = mmppClosed();
        SolverNC solver = new SolverNC(model);
        solver.options.verbose = VerboseLevel.SILENT;
        assertTrue(solver.needsMapEnv(solver.options));
        SolverResult res = solver.getAvg();
        assertEquals(3.4427, res.TN.get(0, 0), 1e-3);
    }

    @Test
    public void ctmcSolvesTheSameModelNatively() {
        Network model = mmppClosed();
        SolverCTMC solver = new SolverCTMC(model);
        solver.options.verbose = VerboseLevel.SILENT;
        assertFalse(solver.needsMapEnv(solver.options), "CTMC declares MAP/MMPP2, so nothing may be intercepted");
        SolverResult res = solver.getAvg();
        assertEquals(3.0952, res.TN.get(0, 0), 1e-3);
    }

    @Test
    public void forcedLimitsMatchMatlab() {
        Network model = mmppClosed();
        SolverMVA dec = new SolverMVA(model);
        dec.options.verbose = VerboseLevel.SILENT;
        dec.options.config.map_env_method = "dec";
        assertEquals(2.3424, dec.getAvg().TN.get(0, 0), 1e-3);

        SolverMVA avg = new SolverMVA(mmppClosed());
        avg.options.verbose = VerboseLevel.SILENT;
        avg.options.config.map_env_method = "avg";
        assertEquals(3.4427, avg.getAvg().TN.get(0, 0), 1e-3);
    }

    @Test
    public void mapEnvOffRestoresTheRejection() {
        Network model = mmppClosed();
        SolverMVA solver = new SolverMVA(model);
        solver.options.verbose = VerboseLevel.SILENT;
        solver.options.config.map_env = "off";
        assertFalse(solver.needsMapEnv(solver.options));
        assertThrows(RuntimeException.class, solver::getAvg);
    }

    @Test
    public void generalMapImageIsIntensityMatched() {
        Network model = genMapClosed();
        MAPQN2RENV.RenvImage image = MAPQN2RENV.map2renvImage(model, null);
        assertEquals(2, image.nstages);
        assertFalse(image.isMMPP, "a non-diagonal D1 is not an MMPP");

        SolverMVA solver = new SolverMVA(model);
        solver.options.verbose = VerboseLevel.SILENT;
        assertEquals(1.9317, solver.getAvg().TN.get(0, 0), 1e-3);
    }

    @Test
    public void twoModulatedProcessesGiveAProductStageSpace() {
        Network model = new Network("twoMaps");
        Delay delay = new Delay(model, "Think");
        Queue q1 = new Queue(model, "Q1", SchedStrategy.FCFS);
        Queue q2 = new Queue(model, "Q2", SchedStrategy.FCFS);
        ClosedClass jobclass = new ClosedClass(model, "C1", 3, delay);
        delay.setService(jobclass, new Exp(1.0));
        q1.setService(jobclass, new MMPP2(1.0, 4.0, 0.5, 0.5));
        q2.setService(jobclass, new MMPP2(2.0, 6.0, 0.4, 0.6));
        model.link(model.serialRouting(delay, q1, q2));

        MAPQN2RENV.RenvImage image = MAPQN2RENV.map2renvImage(model, null);
        assertEquals(4, image.nstages);
        assertEquals(2, image.orders[0]);
        assertEquals(2, image.orders[1]);

        SolverMVA solver = new SolverMVA(model);
        solver.options.verbose = VerboseLevel.SILENT;
        assertEquals(1.5041, solver.getAvg().TN.get(0, 0), 1e-3);
    }

    @Test
    public void transientCapabilityGatesTheMeanFieldCoupling() {
        Network model = mmppClosed();
        SolverMVA mva = new SolverMVA(model);
        assertFalse(mva.supportsTransientAnalysis());
        assertTrue(new SolverCTMC(mmppClosed()).supportsTransientAnalysis());

        mva.options.verbose = VerboseLevel.SILENT;
        mva.options.config.map_env_method = "meanfield";
        assertThrows(RuntimeException.class, mva::getAvg);
    }
}
