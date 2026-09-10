/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.env;

import jline.VerboseLevel;
import jline.lang.ClosedClass;
import jline.lang.Env;
import jline.lang.Network;
import jline.lang.constant.SchedStrategy;
import jline.lang.constant.SolverType;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.processes.Exp;
import jline.solvers.NetworkSolver;
import jline.solvers.SolverResult;
import jline.solvers.SolverOptions;
import jline.solvers.ctmc.SolverCTMC;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

/**
 * Regression: the ENV state-vector analyzer must report the same utilization as the plain
 * CTMC for a load-dependent station.
 *
 * <p>{@code Ctmc_avg_from_pi}, which backs the ENV state-vector path, had drifted from
 * {@code Solver_ctmc_analyzer}: in its load/class-dependent branch it accumulated the bare
 * per-class capacity share {@code nir(k)*schedparam/(nir*schedparam')} without weighting by
 * the current scaling {@code lldnow} or normalizing by the effective capacity
 * {@code ceff = max(nservers, max lldscaling)}. That is a P(busy)-style value, not a
 * busy-server fraction, and it overstated the utilization of a load-dependent station.</p>
 *
 * <p>The oracle used here needs no external reference: an environment whose stages are
 * IDENTICAL cannot change anything, so it must reproduce the plain CTMC solution of that
 * one model exactly. Queue lengths already matched before the fix, which isolates the
 * defect to the utilization branch.</p>
 */
public class SolverEnvLldUtilTest {

    private static final double TOL = 1e-9;
    private static final int N = 3;
    private static final int C = 2;

    /** Closed model with a delay and a load-dependent PS queue, alpha(n) = min(n,c). */
    private static Network build() {
        Network m = new Network("base");
        Delay d = new Delay(m, "D");
        Queue q = new Queue(m, "Q", SchedStrategy.PS);
        ClosedClass jc = new ClosedClass(m, "C", N, d);
        d.setService(jc, new Exp(1.0));
        q.setService(jc, new Exp(2.0));
        Matrix lld = new Matrix(1, N);
        for (int n = 1; n <= N; n++) {
            lld.set(0, n - 1, Math.min(n, C));
        }
        q.setLoadDependence(lld);
        m.link(m.serialRouting(d, q));
        return m;
    }

    @Test
    public void envStatevecMatchesCtmcOnIdenticalLoadDependentStages() {
        Matrix want = new SolverCTMC(build()).getAvgUtil();

        int E = 2;
        Env envModel = new Env("IdenticalStages", E);
        envModel.addStage(0, "S1", "operational", build());
        envModel.addStage(1, "S2", "operational", build());
        envModel.addTransition(0, 1, new Exp(0.5));
        envModel.addTransition(1, 0, new Exp(0.5));

        SolverOptions options = new SolverOptions(SolverType.ENV);
        options.method = "statevec";
        options.verbose = VerboseLevel.SILENT;

        NetworkSolver[] solvers = new NetworkSolver[E];
        for (int e = 0; e < E; e++) {
            SolverOptions ctmcOptions = new SolverOptions(SolverType.CTMC);
            ctmcOptions.verbose = VerboseLevel.SILENT;
            ctmcOptions.timespan[0] = 0;
            ctmcOptions.timespan[1] = 1e6;
            solvers[e] = new SolverCTMC(envModel.getModel(e));
            solvers[e].options = ctmcOptions;
        }

        SolverENV envSolver = new SolverENV(envModel, solvers, options);
        envSolver.getAvg();
        Matrix got = envSolver.result.UN;

        // The environment is a no-op, so every station's utilization must agree exactly.
        for (int i = 0; i < want.getNumRows(); i++) {
            assertEquals(want.get(i, 0), got.get(i, 0), TOL,
                    "ENV statevec Util must match the plain CTMC at station " + i);
        }
    }
}
