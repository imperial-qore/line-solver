package jline.solvers.env;

import jline.lang.Environment;
import jline.lang.constant.SchedStrategy;
import jline.lang.constant.SolverType;
import jline.lang.layered.Activity;
import jline.lang.layered.Entry;
import jline.lang.layered.LayeredNetwork;
import jline.lang.layered.Processor;
import jline.lang.layered.Task;
import jline.lang.processes.Exp;
import jline.solvers.Solver;
import jline.solvers.SolverOptions;
import jline.solvers.SolverResult;
import jline.solvers.ctmc.SolverCTMC;
import jline.solvers.ln.LNTranAvgResult;
import jline.solvers.ln.SolverFactory;
import jline.solvers.ln.SolverLN;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.Timeout;

import static org.junit.jupiter.api.Assertions.*;

/**
 * SolverENV running over LayeredNetwork (LQN) stages. The environment
 * alternates between fast/slow-database LQN stages; the environment-averaged
 * metrics are validated against the single-stage LQN solutions.
 */
public class SolverEnvLqnTest {

    private static final SolverFactory CTMC = m -> new SolverCTMC(m);

    private static LayeredNetwork buildLqn(String name, double dbMean) {
        LayeredNetwork model = new LayeredNetwork(name);
        Processor P1 = new Processor(model, "CP", 1, SchedStrategy.PS);
        Processor P2 = new Processor(model, "DP", 1, SchedStrategy.PS);
        Task T1 = new Task(model, "CT", 3, SchedStrategy.REF).on(P1);
        T1.setThinkTime(new Exp(1.0 / 5.0));
        Task T2 = new Task(model, "DT", 1, SchedStrategy.INF).on(P2);
        Entry E1 = new Entry(model, "CE").on(T1);
        Entry E2 = new Entry(model, "DE").on(T2);
        new Activity(model, "CA", new Exp(1.0)).on(T1).boundTo(E1).synchCall(E2, 2.0);
        new Activity(model, "DA", new Exp(1.0 / dbMean)).on(T2).boundTo(E2).repliesTo(E2);
        return model;
    }

    private static SolverOptions lnOptions(double T) {
        SolverOptions o = new SolverOptions(SolverType.LN);
        o.timespan = new double[]{0.0, T};
        o.iter_max = 12;
        o.verbose = jline.VerboseLevel.SILENT;
        return o;
    }

    private static SolverOptions envOptions() {
        SolverOptions o = new SolverOptions(SolverType.ENV);
        o.iter_max = 4;
        o.iter_tol = 0.05;
        o.verbose = jline.VerboseLevel.SILENT;
        return o;
    }

    /** Sum of the finite tail (last time point) values over a dense block. */
    private static double tailSum(Matrix[][] cells) {
        double s = 0.0;
        if (cells == null) return 0.0;
        for (Matrix[] row : cells) {
            if (row == null) continue;
            for (Matrix c : row) {
                if (c != null && c.getNumRows() > 0) {
                    double v = c.get(c.getNumRows() - 1, 0);
                    if (!Double.isNaN(v) && !Double.isInfinite(v)) s += v;
                }
            }
        }
        return s;
    }

    private static double sumFinite(Matrix m) {
        double s = 0.0;
        for (int i = 0; i < m.getNumRows(); i++) {
            for (int j = 0; j < m.getNumCols(); j++) {
                double v = m.get(i, j);
                if (!Double.isNaN(v) && !Double.isInfinite(v)) s += v;
            }
        }
        return s;
    }

    /** {aggregate Q, aggregate T} steady tail for one LQN stage. */
    private static double[] stageAgg(LayeredNetwork model, SolverFactory f, double T) {
        SolverLN s = new SolverLN(model, f, lnOptions(T));
        LNTranAvgResult r = s.getTranAvg();
        return new double[]{tailSum(r.QNt), tailSum(r.TNt)};
    }

    private static double envTputSum(SolverENV env) {
        env.getAvg();
        return sumFinite(env.result.TN);
    }

    // -----------------------------------------------------------------------

    @Test
    @Timeout(600)
    public void lqnInEnvironmentIsBracketed() {
        double T = 8.0;
        LayeredNetwork up = buildLqn("LQN_UP", 0.8);
        LayeredNetwork down = buildLqn("LQN_DOWN", 3.0);

        Environment renv = new Environment("DBReliability", 2);
        renv.addStage(0, "UP", "operational", up);
        renv.addStage(1, "DOWN", "degraded", down);
        renv.addTransition(0, 1, new Exp(0.2));
        renv.addTransition(1, 0, new Exp(1.0));

        Solver[] solvers = new Solver[]{
            new SolverLN(up, CTMC, lnOptions(T)),
            new SolverLN(down, CTMC, lnOptions(T))
        };
        SolverENV envSolver = new SolverENV(renv, solvers, envOptions());
        envSolver.getAvg();
        SolverResult res = envSolver.result;

        assertNotNull(res.QN);
        assertNotNull(res.TN);
        assertFalse(res.QN.hasNaN(), "ENV aggregate QN has NaN");

        double qEnv = res.QN.elementSum();
        double tEnv = sumFinite(res.TN);
        double[] upA = stageAgg(up, CTMC, T);
        double[] downA = stageAgg(down, CTMC, T);

        assertTrue(Math.abs(qEnv - upA[0]) <= 0.05 * upA[0],
                String.format("population not conserved: env=%.4f up=%.4f", qEnv, upA[0]));
        double lo = Math.min(upA[1], downA[1]), hi = Math.max(upA[1], downA[1]);
        double tol = 0.02 * Math.max(1.0, hi);
        assertTrue(tEnv >= lo - tol && tEnv <= hi + tol,
                String.format("throughput not bracketed: env=%.4f in [%.4f, %.4f]", tEnv, lo, hi));
    }

    /** Two-stage UP/DOWN environment throughput at switch rates a (UP->DOWN)
     *  and b (DOWN->UP); probEnv(UP)=b/(a+b). */
    private static double env2StageTput(double a, double b, double T) {
        LayeredNetwork up = buildLqn("UP", 0.8);
        LayeredNetwork down = buildLqn("DOWN", 3.0);
        Environment renv = new Environment("R", 2);
        renv.addStage(0, "UP", "operational", up);
        renv.addStage(1, "DOWN", "degraded", down);
        renv.addTransition(0, 1, new Exp(a));
        renv.addTransition(1, 0, new Exp(b));
        Solver[] solvers = new Solver[]{
            new SolverLN(up, CTMC, lnOptions(T)),
            new SolverLN(down, CTMC, lnOptions(T))
        };
        return envTputSum(new SolverENV(renv, solvers, envOptions()));
    }

    /**
     * Quantitative coupling check: the environment-averaged throughput increases
     * monotonically with the stationary probability of the fast (UP) stage.
     * probEnv(UP)=b/(a+b), so raising b/a spends more time in UP and must raise
     * the aggregate throughput, while both settings stay bracketed by the
     * single-stage solutions. This validates the probEnv-weighted coupling.
     */
    @Test
    @Timeout(600)
    public void envThroughputMonotoneInStageProbability() {
        double T = 8.0;
        double tLowPUp = env2StageTput(1.0, 0.2, T);   // probEnv(UP)=0.167
        double tHighPUp = env2StageTput(0.2, 1.0, T);  // probEnv(UP)=0.833

        double tUp = stageAgg(buildLqn("UP", 0.8), CTMC, T)[1];
        double tDown = stageAgg(buildLqn("DOWN", 3.0), CTMC, T)[1];
        double lo = Math.min(tUp, tDown), hi = Math.max(tUp, tDown);
        double tol = 0.03 * Math.max(1.0, hi);

        assertTrue(tHighPUp > tLowPUp + 1e-3,
                String.format("env throughput not monotone in probEnv(UP): lowP=%.4f highP=%.4f",
                        tLowPUp, tHighPUp));
        assertTrue(tLowPUp >= lo - tol && tHighPUp <= hi + tol,
                String.format("env throughputs out of [%.4f, %.4f]: low=%.4f high=%.4f",
                        lo, hi, tLowPUp, tHighPUp));
    }

    /** Three-stage environment (UP / MID / DOWN): env-average is bracketed by
     *  the fastest and slowest single-stage solutions. Exercises the E>2
     *  stage-coupling matrix. */
    @Test
    @Timeout(600)
    public void threeStageEnvironmentIsBracketed() {
        double T = 8.0;
        LayeredNetwork up = buildLqn("UP", 0.8);
        LayeredNetwork mid = buildLqn("MID", 1.6);
        LayeredNetwork down = buildLqn("DOWN", 3.0);

        Environment renv = new Environment("R3", 3);
        renv.addStage(0, "UP", "operational", up);
        renv.addStage(1, "MID", "degraded", mid);
        renv.addStage(2, "DOWN", "failed", down);
        renv.addTransition(0, 1, new Exp(0.3));
        renv.addTransition(1, 2, new Exp(0.3));
        renv.addTransition(2, 1, new Exp(0.6));
        renv.addTransition(1, 0, new Exp(0.6));

        Solver[] solvers = new Solver[]{
            new SolverLN(up, CTMC, lnOptions(T)),
            new SolverLN(mid, CTMC, lnOptions(T)),
            new SolverLN(down, CTMC, lnOptions(T))
        };
        SolverENV envSolver = new SolverENV(renv, solvers, envOptions());
        envSolver.getAvg();

        assertFalse(envSolver.result.QN.hasNaN(), "3-stage ENV QN has NaN");
        double tEnv = sumFinite(envSolver.result.TN);
        double tUp = stageAgg(up, CTMC, T)[1];
        double tDown = stageAgg(down, CTMC, T)[1];
        double lo = Math.min(tUp, tDown), hi = Math.max(tUp, tDown);
        double tol = 0.03 * Math.max(1.0, hi);
        assertTrue(tEnv >= lo - tol && tEnv <= hi + tol,
                String.format("3-stage throughput not bracketed: env=%.4f in [%.4f, %.4f]", tEnv, lo, hi));
    }

}
