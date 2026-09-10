package jline.solvers.fluid;

import jline.api.sym.SymEngine;
import jline.api.sym.SymEngines;
import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.processes.Exp;
import jline.solvers.SolverOptions;
import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Tag;
import org.junit.jupiter.api.Test;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import static org.junit.jupiter.api.Assertions.*;
import static org.junit.jupiter.api.Assumptions.assumeTrue;

/**
 * Symbolic drift and Jacobian of the mean-field ODE system, through the
 * line-sage-rest service.
 *
 * <p>The drift is checked against the field the solver actually integrates,
 * recomputed here in double precision, so a divergence in the exported
 * expression shows up as a number rather than as a differently spelled
 * formula. Skipped when no symbolic backend can be resolved.</p>
 */
@Tag("remote")
public class SolverFluidSymbolicSageTest {

    private static SymEngine engine;

    @BeforeAll
    public static void setUp() {
        engine = SymEngines.resolve("auto");
        if (engine == null) {
            System.err.println("[sage] No symbolic backend available; skipping. Start one with "
                    + "'docker run -d -p 8080:8080 " + SymEngines.DOCKER_IMAGE + "'.");
        }
    }

    @AfterAll
    public static void tearDown() {
        SymEngines.stopContainer();
    }

    /** Closed model: Delay Exp(1) -> PS Queue Exp(2) with 2 servers, N = 4. */
    private static Network model() {
        Network model = new Network("FluidSymbolic");
        Delay delay = new Delay(model, "Think");
        Queue queue = new Queue(model, "Queue1", SchedStrategy.PS);
        queue.setNumberOfServers(2);
        ClosedClass jobs = new ClosedClass(model, "Class1", 4, delay);
        delay.setService(jobs, new Exp(1.0));
        queue.setService(jobs, new Exp(2.0));
        model.link(model.serialRouting(delay, queue));
        return model;
    }

    private static SolverOptions options(boolean smooth) {
        SolverOptions options = new SolverOptions(jline.lang.constant.SolverType.FLUID);
        options.method = "matrix";
        if (smooth) {
            options.config.pstar = new ArrayList<Double>();
            options.config.pstar.add(8.0);
        }
        return options;
    }

    @Test
    public void minScaledDriftIsRefusedNotOneSidedlyDifferentiated() {
        SolverFluid solver = new SolverFluid(model(), options(false));
        RuntimeException e = assertThrows(RuntimeException.class, solver::getSymbolicDrift);
        assertTrue(e.getMessage().contains("min(n_i, S_i)"),
                "the refusal must name the non-differentiable factor, got: " + e.getMessage());
    }

    @Test
    public void smoothDriftMatchesTheIntegratedField() throws IOException {
        assumeTrue(engine != null);
        SolverFluid solver = new SolverFluid(model(), options(true));
        FluidODEsExporter.SymODEs sys =
                FluidODEsExporter.build(solver.getStruct(), options(true));
        List<String> rhs = FluidODEsExporter.symbolicDrift(sys);
        List<String> vars = FluidODEsExporter.stateVariables(sys);
        int n = sys.nstates;
        assertEquals(n, rhs.size());

        double[] x = new double[n];
        for (int s = 0; s < n; s++) {
            x[s] = 0.5 + 0.37 * (s + 1);
        }
        Map<String, Double> assignment = new HashMap<String, Double>();
        for (int s = 0; s < n; s++) {
            assignment.put(vars.get(s), x[s]);
        }
        double[] symbolic = engine.eval(rhs, assignment);
        double[] reference = pnormDrift(sys, x);
        for (int s = 0; s < n; s++) {
            assertEquals(reference[s], symbolic[s], 1e-9,
                    "state " + s + " of the exported drift disagrees with the integrated field");
        }
    }

    @Test
    public void jacobianMatchesACentralDifferenceOfTheField() {
        assumeTrue(engine != null);
        SolverFluid solver = new SolverFluid(model(), options(true));
        FluidODEsExporter.SymODEs sys =
                FluidODEsExporter.build(solver.getStruct(), options(true));
        int n = sys.nstates;
        String[][] J = solver.getJacobian();
        assertEquals(n, J.length);

        double[] x = new double[n];
        for (int s = 0; s < n; s++) {
            x[s] = 0.5 + 0.37 * (s + 1);
        }
        Map<String, Double> assignment = new HashMap<String, Double>();
        List<String> vars = FluidODEsExporter.stateVariables(sys);
        for (int s = 0; s < n; s++) {
            assignment.put(vars.get(s), x[s]);
        }
        double h = 1e-6;
        for (int i = 0; i < n; i++) {
            double[] row;
            try {
                row = engine.eval(java.util.Arrays.asList(J[i]), assignment);
            } catch (IOException e) {
                throw new RuntimeException(e);
            }
            for (int j = 0; j < n; j++) {
                double[] xp = x.clone();
                double[] xm = x.clone();
                xp[j] += h;
                xm[j] -= h;
                double numeric = (pnormDrift(sys, xp)[i] - pnormDrift(sys, xm)[i]) / (2 * h);
                assertEquals(numeric, row[j], 1e-6,
                        "J[" + i + "][" + j + "] disagrees with the differenced field");
            }
        }
    }

    /**
     * The field the p-norm method integrates, in double precision:
     * dx/dt = W' * (x .* ghat) + Alambda with
     * ghat_s = 1/(1 + (n_i/S_i)^p_i)^(1/p_i), and theta = 0 at a Source.
     * This is pnorm_ode of the MATLAB matrix analyzer, transcribed.
     */
    private static double[] pnormDrift(FluidODEsExporter.SymODEs sys, double[] x) {
        int n = sys.nstates;
        double[] theta = new double[n];
        for (int s = 0; s < n; s++) {
            if (sys.isSource[s]) {
                continue;
            }
            double ni = jline.GlobalConstants.FineTol;
            for (int k = 0; k < n; k++) {
                if (sys.stateStation[k] == sys.stateStation[s]) {
                    ni += x[k];
                }
            }
            double S = sys.S[sys.stateStation[s]];
            double p = sys.pstar[sys.stateStation[s]];
            double ghat = (S > 0 && p > 0)
                    ? 1.0 / Math.pow(1 + Math.pow(ni / S, p), 1.0 / p) : 1.0;
            theta[s] = x[s] * ghat;
        }
        double[] f = new double[n];
        for (int s = 0; s < n; s++) {
            double acc = sys.Alambda[s];
            for (int t = 0; t < n; t++) {
                acc += sys.W.get(t, s) * theta[t];
            }
            f[s] = acc;
        }
        return f;
    }
}
