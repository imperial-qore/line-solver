/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.env;

import jline.VerboseLevel;
import jline.lang.Env;
import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.constant.SchedStrategy;
import jline.lang.constant.SolverType;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.processes.Exp;
import jline.solvers.NetworkSolver;
import jline.solvers.SolverOptions;
import jline.solvers.ctmc.SolverCTMC;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

/**
 * Regression: the ENV state-vector analyzer must be CUTOFF-INVARIANT when one environment
 * stage is individually unstable.
 *
 * <p>{@code svPre} seeded every stage's entry distribution with THAT STAGE'S OWN stationary
 * law. A stage whose arrival rate is at or above its own service rate has no stationary law,
 * so {@code ctmc_solve_reducible} returned the stationary law of the TRUNCATED generator --
 * mass piled against the truncation wall, mean growing linearly with the cutoff (uniform,
 * mean N/2, for a critical M/M/1 truncated at N).</p>
 *
 * <p>The fixed point itself was never wrong: chaining the resolvent is exactly the stationary
 * equation of the joint (queue, stage) chain, so it contracts to the right answer. What grew
 * with the cutoff was the number of sweeps needed to drain the seed, which made the reported
 * number drift AWAY from the truth as the cutoff was RAISED -- the natural response to a
 * suspect answer made it worse, and a small cutoff looking right was a coincidence.</p>
 *
 * <p>Here the Slow stage is critical (lambda = mu = 0.8) while the system is stable on average
 * (mu_bar = 1.4), the regime that exposed it. Both cutoffs sit far past the real mass, so
 * anything that moves between them is the seed and not the model.</p>
 *
 * <p>THE SWEEP BUDGET IS PART OF THE TEST. With the fixed seed the fixed point is reached in
 * 162 sweeps at EVERY cutoff, so iter_max = 200 is a comfortable 19% margin; with the old
 * per-stage seed it needs 293 sweeps at cutoff 40 and 456 at cutoff 100, so neither finishes.
 * Raising iter_max far enough lets the old seed converge too and SILENTLY DEFANGS this test --
 * it would still pass while testing nothing. Measured on the reference iteration, the old seed
 * under this budget gives 1.478346 at cutoff 40 and 1.600286 at cutoff 100 (throughput 0.808905
 * against an arrival rate of 0.8), versus an exact 1.478318.</p>
 */
public class SolverEnvStatevecCutoffTest {

    /** Exact joint (queue, stage) CTMC answer for this environment. */
    private static final double EXACT_Q = 1.478318;
    private static final double LAMBDA = 0.8;

    /** Open M/M/1: Source -> FCFS Queue -> Sink, one class. */
    private static Network mm1(double lambda, double mu) {
        Network m = new Network("mm1");
        Source src = new Source(m, "Source");
        Queue q = new Queue(m, "Q", SchedStrategy.FCFS);
        Sink snk = new Sink(m, "Sink");
        OpenClass c = new OpenClass(m, "C");
        src.setArrival(c, new Exp(lambda));
        q.setService(c, new Exp(mu));
        m.link(m.serialRouting(src, q, snk));
        return m;
    }

    /** Total queue length and the queue's throughput at the given per-class cutoff. */
    private static double[] solve(int cutoff) {
        int E = 2;
        Env envModel = new Env("Modes", E);
        envModel.addStage(0, "Slow", "degraded", mm1(LAMBDA, 0.8));
        envModel.addStage(1, "Fast", "operational", mm1(LAMBDA, 2.0));
        envModel.addTransition(0, 1, new Exp(1.0));
        envModel.addTransition(1, 0, new Exp(1.0));

        SolverOptions options = new SolverOptions(SolverType.ENV);
        options.method = "statevec";
        options.verbose = VerboseLevel.SILENT;
        // See the class comment: this budget is chosen so the OLD seed cannot finish.
        options.iter_max = 200;
        options.iter_tol = 1e-8;

        NetworkSolver[] solvers = new NetworkSolver[E];
        for (int e = 0; e < E; e++) {
            SolverOptions ctmcOptions = new SolverOptions(SolverType.CTMC);
            ctmcOptions.verbose = VerboseLevel.SILENT;
            ctmcOptions.cutoff = Matrix.singleton(cutoff);
            ctmcOptions.timespan[0] = 0;
            ctmcOptions.timespan[1] = 1e6;
            solvers[e] = new SolverCTMC(envModel.getModel(e));
            solvers[e].options = ctmcOptions;
        }

        SolverENV envSolver = new SolverENV(envModel, solvers, options);
        envSolver.getAvg();
        Matrix QN = envSolver.result.QN;
        Matrix TN = envSolver.result.TN;
        // Station 1 is the Queue; station 0 is the Source, whose queue length is
        // infinite by convention and must NOT be summed into the total.
        double q = 0.0;
        double t = 0.0;
        for (int r = 0; r < QN.getNumCols(); r++) {
            q += QN.get(1, r);
            t += TN.get(1, r);
        }
        return new double[]{q, t};
    }

    @Test
    public void criticalStageIsCutoffInvariant() {
        double[] small = solve(40);
        double[] large = solve(100);

        // Flow balance at the queue is the cheap tell the defect violated by up to 33%:
        // in steady state the queue must clear exactly what arrives.
        assertEquals(LAMBDA, small[1], 1e-5, "queue throughput at cutoff 40 must equal the arrival rate");
        assertEquals(LAMBDA, large[1], 1e-5, "queue throughput at cutoff 100 must equal the arrival rate");

        assertEquals(EXACT_Q, small[0], 1e-3, "ENV statevec QLen at cutoff 40");
        assertEquals(EXACT_Q, large[0], 1e-3, "ENV statevec QLen at cutoff 100");

        // The invariance itself, independent of how close either is to the truth.
        assertEquals(small[0], large[0], 1e-5, "raising the cutoff must not move the answer");
    }
}
