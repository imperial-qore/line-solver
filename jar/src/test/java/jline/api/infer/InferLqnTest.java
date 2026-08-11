package jline.api.infer;

import jline.VerboseLevel;
import jline.lang.constant.SchedStrategy;
import jline.lang.layered.Activity;
import jline.lang.layered.Entry;
import jline.lang.layered.LayeredNetwork;
import jline.lang.layered.Processor;
import jline.lang.layered.Task;
import jline.lang.processes.Exp;
import jline.solvers.LayeredNetworkAvgTable;
import jline.solvers.SolverOptions;
import jline.solvers.ln.LNOptions;
import jline.solvers.ln.SolverLN;
import jline.util.matrix.Matrix;

import org.junit.jupiter.api.Test;

import java.util.Arrays;
import java.util.List;
import java.util.function.Function;

import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Tests for JAR LQN parameter identification via the Extended Kalman Filter
 * ({@link InferLqn}). Method: Zheng, Yang, Woodside, Litoiu, Iszlai, CASCON
 * 2005. Two hidden parameters are estimated (reference-task think time and the
 * P2 host demand of activity AS3) from [R(E1), U(P1), U(P2)].
 */
class InferLqnTest {

    private static LayeredNetwork buildLqn() {
        LayeredNetwork m = new LayeredNetwork("paramident_LQN");
        Processor P1 = new Processor(m, "P1", 1, SchedStrategy.PS);
        Processor P2 = new Processor(m, "P2", 1, SchedStrategy.PS);
        Task T1 = new Task(m, "T1", 5, SchedStrategy.REF).on(P1);
        T1.setThinkTime(Exp.fitMean(0.5));
        Task T2 = new Task(m, "T2", 5, SchedStrategy.FCFS).on(P1);
        T2.setThinkTime(Exp.fitMean(1.0 / 3.0));
        Task T3 = new Task(m, "T3", 3, SchedStrategy.FCFS).on(P2);
        T3.setThinkTime(Exp.fitMean(0.25));
        Entry E1 = new Entry(m, "E1").on(T1);
        Entry E2 = new Entry(m, "E2").on(T2);
        Entry E3 = new Entry(m, "E3").on(T3);
        new Activity(m, "AS1", Exp.fitMean(0.1)).on(T1).boundTo(E1).synchCall(E2, 1);
        new Activity(m, "AS2", Exp.fitMean(0.05)).on(T2).boundTo(E2).synchCall(E3, 5).repliesTo(E2);
        new Activity(m, "AS3", Exp.fitMean(0.02)).on(T3).boundTo(E3).repliesTo(E3);
        return m;
    }

    private static Matrix solveObs(LayeredNetwork model, List<ObsSpec> obsSpec) {
        SolverOptions so = new LNOptions();
        so.verbose = VerboseLevel.SILENT;
        LayeredNetworkAvgTable table = (LayeredNetworkAvgTable) new SolverLN(model, so).getEnsembleAvg();
        return InferLqn.getObs(table, obsSpec);
    }

    @Test
    public void testNoiseFreeRecovery() {
        List<ParamSpec> paramSpec = Arrays.asList(
                ParamSpec.think("T1"), ParamSpec.hostDemand("AS3"));
        List<ObsSpec> obsSpec = Arrays.asList(
                ObsSpec.respT("E1"), ObsSpec.util("P1"), ObsSpec.util("P2"));

        Matrix aTrue = new Matrix(2, 1);
        aTrue.set(0, 0, 0.5);
        aTrue.set(1, 0, 1.0 / 50.0);

        // noise-free measurement (constant), single column tiled over the sequence.
        // Each EKF step costs 1 + paramSpec.size() SolverLN solves (base + one
        // forward-difference column per parameter); on a constant noise-free Z the
        // filter is converged by step 4, so 8 steps keep margin on every assertion
        // below. Matches the MATLAB/Python tests.
        LayeredNetwork truth = buildLqn();
        InferLqn.setParams(truth, paramSpec, aTrue);
        Matrix z = solveObs(truth, obsSpec);
        int nsteps = 8;
        Matrix Z = new Matrix(z.getNumRows(), nsteps);
        for (int c = 0; c < nsteps; c++) {
            for (int r = 0; r < z.getNumRows(); r++) {
                Z.set(r, c, z.get(r, 0));
            }
        }

        InferLqnOptions opt = new InferLqnOptions();
        opt.a0 = new Matrix(2, 1);
        opt.a0.set(0, 0, 0.8);
        opt.a0.set(1, 0, 1.0 / 35.0);
        opt.QFac = 1e-3;
        opt.gammaT = 1.0;
        opt.aTrue = aTrue;

        LayeredNetwork model = buildLqn();
        InferLqnResult res = InferLqn.inferLqn(model, paramSpec, obsSpec, Z, opt);

        double zHat = res.ahat.get(0, nsteps - 1);
        double sdHat = res.ahat.get(1, nsteps - 1);
        assertTrue(Math.abs(zHat - 0.5) / 0.5 < 1e-2, "think-time relerr too large: " + zHat);
        assertTrue(Math.abs(sdHat - 0.02) / 0.02 < 5e-3, "host-demand relerr too large: " + sdHat);
        assertTrue(res.Er < 1e-1, "prediction RMS too large: " + res.Er);
    }

    @Test
    public void testUtilMonotoneInDemand() {
        List<ParamSpec> paramSpec = Arrays.asList(
                ParamSpec.think("T1"), ParamSpec.hostDemand("AS3"));
        List<ObsSpec> obsSpec = Arrays.asList(ObsSpec.util("P2"));

        LayeredNetwork model = buildLqn();
        Matrix aLo = new Matrix(2, 1);
        aLo.set(0, 0, 0.5);
        aLo.set(1, 0, 1.0 / 50.0);
        InferLqn.setParams(model, paramSpec, aLo);
        double zLo = solveObs(model, obsSpec).get(0, 0);

        Matrix aHi = new Matrix(2, 1);
        aHi.set(0, 0, 0.5);
        aHi.set(1, 0, 1.0 / 25.0);
        InferLqn.setParams(model, paramSpec, aHi);
        double zHi = solveObs(model, obsSpec).get(0, 0);

        assertTrue(zHi > zLo, "larger demand must give larger utilization");
    }

    @Test
    public void testJacobianAnalytic() {
        Function<Matrix, Matrix> hfun = new Function<Matrix, Matrix>() {
            @Override
            public Matrix apply(Matrix a) {
                Matrix h = new Matrix(3, 1);
                h.set(0, 0, a.get(0, 0) * a.get(0, 0));
                h.set(1, 0, 2.0 * a.get(1, 0));
                h.set(2, 0, a.get(0, 0) * a.get(1, 0));
                return h;
            }
        };
        Matrix a = new Matrix(2, 1);
        a.set(0, 0, 3.0);
        a.set(1, 0, 5.0);
        InferLqn.JacobianResult jr = InferLqn.jacobian(hfun, a, 1e-6, 1e-9);
        double[][] jexact = {{2 * 3.0, 0.0}, {0.0, 2.0}, {5.0, 3.0}};
        for (int r = 0; r < 3; r++) {
            for (int c = 0; c < 2; c++) {
                assertTrue(Math.abs(jr.H.get(r, c) - jexact[r][c]) < 1e-3,
                        "jacobian mismatch at (" + r + "," + c + ")");
            }
        }
    }
}
