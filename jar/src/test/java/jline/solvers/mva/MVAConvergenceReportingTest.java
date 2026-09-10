/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.mva;

import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.RoutingMatrix;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.processes.Exp;
import jline.solvers.SolverOptions;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * The AMVA iteration count and convergence flag must reach an API caller.
 *
 * <p>Two defects made them unreadable. {@code MVAResult} redeclared {@code iter},
 * shadowing {@link jline.solvers.SolverResult#iter}: {@code setAvgResults} writes
 * the base field through a {@code SolverResult} reference while the analyzer wrote
 * the subclass one, so a caller that cast the result to {@code MVAResult} -- which
 * is what reading an MVA-specific field forces -- saw 0 on every method. And
 * {@code runAnalyzer} never carried the analyzer's {@code converged} flag onto the
 * solver's own result container, which is a different object, so the flag was
 * always null downstream.
 *
 * <p>The count alone cannot decide convergence on the load-dependent (multiserver)
 * route: {@code Solver_amvald} increments {@code totiter} inside the nested
 * per-chain and inner loops, so it aggregates inner sweeps and saturates the
 * budget by construction on a solve whose outer residual is exactly zero. That is
 * why the flag exists and why it has to be readable. See G31(b) in line-gaps.md.
 */
public class MVAConvergenceReportingTest {

    /** Multiclass, multiserver closed model: dispatches through the amvald route. */
    private static Network multiserverClosedModel() {
        Network model = new Network("amvald_convergence");
        Delay think = new Delay(model, "Think");
        int nstations = 4;
        int[] nservers = {1, 2, 3, 2};
        List<Queue> queues = new ArrayList<Queue>();
        for (int i = 0; i < nstations; i++) {
            Queue q = new Queue(model, "Q" + (i + 1), SchedStrategy.PS);
            q.setNumberOfServers(nservers[i]);
            queues.add(q);
        }
        int nclasses = 3;
        int[] njobs = {6, 8, 5};
        List<ClosedClass> classes = new ArrayList<ClosedClass>();
        for (int r = 0; r < nclasses; r++) {
            classes.add(new ClosedClass(model, "C" + (r + 1), njobs[r], think, 0));
        }
        for (int r = 0; r < nclasses; r++) {
            think.setService(classes.get(r), new Exp(1.0));
            for (int i = 0; i < nstations; i++) {
                queues.get(i).setService(classes.get(r), new Exp(1.0 / (0.3 + 0.2 * i + 0.1 * r)));
            }
        }
        RoutingMatrix P = model.initRoutingMatrix();
        for (int r = 0; r < nclasses; r++) {
            ClosedClass c = classes.get(r);
            P.set(c, c, think, queues.get(0), 1.0);
            for (int i = 0; i < nstations - 1; i++) {
                P.set(c, c, queues.get(i), queues.get(i + 1), 1.0);
            }
            P.set(c, c, queues.get(nstations - 1), think, 1.0);
        }
        model.link(P);
        return model;
    }

    private static SolverMVA solved(String method) {
        return solved(method, null, 0);
    }

    /**
     * @param multiserver config rule, or null for the default. "softmin"
     *        dispatches to Solver_amvald unconditionally, so it pins the flag to
     *        that route independently of how the default rule happens to dispatch.
     * @param iterMax     iteration budget, or 0 to leave the default.
     */
    private static SolverMVA solved(String method, String multiserver, int iterMax) {
        SolverOptions options = new SolverOptions(jline.lang.constant.SolverType.MVA);
        options.method = method;
        options.verbose = jline.VerboseLevel.SILENT;
        if (multiserver != null) {
            options.config.multiserver = multiserver;
        }
        if (iterMax > 0) {
            options.iter_max = iterMax;
        }
        SolverMVA solver = new SolverMVA(multiserverClosedModel(), options);
        solver.getAvgTable();
        return solver;
    }

    @Test
    @DisplayName("the amvald route reports a nonzero iteration count on the result")
    public void iterationCountReachesTheResult() {
        for (String method : new String[]{"lin", "qd", "bs"}) {
            SolverMVA solver = solved(method);
            assertTrue(solver.result.iter > 0,
                    "method '" + method + "' reported iter=" + solver.result.iter
                            + "; the count never reached the result container");
        }
    }

    @Test
    @DisplayName("MVAResult no longer shadows SolverResult.iter")
    public void castingToMvaResultReadsTheSameCount() {
        SolverMVA solver = solved("lin");
        MVAResult mvaResult = (MVAResult) solver.result;
        assertEquals(solver.result.iter, mvaResult.iter,
                "the cast result reports a different iteration count, i.e. the field is shadowed again");
    }

    @Test
    @DisplayName("the amvald convergence flag is carried onto the solver result")
    public void convergenceFlagReachesTheResult() {
        SolverMVA solver = solved("lin", "softmin", 0);
        MVAResult mvaResult = (MVAResult) solver.result;
        assertNotNull(mvaResult.converged,
                "the amvald route reports a convergence flag; it did not reach the result");
        assertTrue(mvaResult.converged.booleanValue(),
                "this model converges: an outer residual within tolerance must not report as non-converged");
    }

    @Test
    @DisplayName("the default multiserver rule carries the flag too, not only softmin")
    public void defaultMultiserverRuleAlsoCarriesTheFlag() {
        // This model reaches Solver_amvald under the default rule as well, so the
        // flag must survive there. A null flag is still legitimate in general --
        // the single-loop handlers (the linearizer family, bs, aql) report none,
        // and for them the count is a sound signal because it is not an aggregate
        // -- but it must not be null merely because a layer forgot to copy it,
        // which is what made every AMVA solve look flagless.
        SolverMVA solver = solved("lin");
        MVAResult mvaResult = (MVAResult) solver.result;
        assertNotNull(mvaResult.converged, "the convergence flag was lost on the default route");
        assertTrue(mvaResult.converged.booleanValue(),
                "this model converges; reporting non-convergence would raise a spurious warning");
        assertTrue(solver.result.iter > 0, "no iteration count reported");
    }
}
