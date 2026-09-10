package jline.solvers.ctmc;

import jline.api.sym.SageRestEngine;
import jline.api.sym.SymEngine;
import jline.api.sym.SymEngines;
import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.processes.Exp;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Tag;
import org.junit.jupiter.api.Test;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import static jline.api.mc.Ctmc_solve.ctmc_solve;
import static org.junit.jupiter.api.Assertions.*;
import static org.junit.jupiter.api.Assumptions.assumeTrue;

/**
 * Symbolic CTMC analysis through the line-sage-rest service.
 *
 * <p>Backend resolution goes through {@link SymEngines#resolve}, the same path
 * production code takes, so the test also covers the service discovery and the
 * container lifecycle. Every test is skipped when no backend can be resolved,
 * which is the case on a machine with no Docker and no service running.</p>
 *
 * <p>To run against a hand started service instead:
 * <pre>docker run --rm -p 8080:8080 imperialqore/line-sage-rest:latest</pre></p>
 */
@Tag("remote")
public class SolverCTMCSymbolicSageTest {

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

    /**
     * Two station closed model: Delay Exp(1), FCFS Queue Exp(2), N = 3.
     *
     * <p>This is the oracle the symbolic CTMC tests use in all three codebases.
     * Its stationary distribution is proportional to the balance product, and
     * at x1 = 1, x2 = 2 it is [0.15789, 0.31579, 0.31579, 0.21053] over the
     * MATLAB state ordering.</p>
     *
     * @return the model
     */
    private static Network oracleModel() {
        Network model = new Network("SymbolicOracle");
        Delay delay = new Delay(model, "Delay");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        ClosedClass jobs = new ClosedClass(model, "Class1", 3, delay);
        delay.setService(jobs, new Exp(1.0));
        queue.setService(jobs, new Exp(2.0));
        model.link(model.serialRouting(delay, queue));
        return model;
    }

    @Test
    public void serviceIdentifiesItself() throws IOException {
        assumeTrue(engine != null);
        assertTrue(engine.isAvailable());
        assertEquals("sage", engine.name());
        assertTrue(((SageRestEngine) engine).info().has("sage_version"),
                "a line-sage-rest service must report its Sage version");
    }

    @Test
    public void symbolicSolutionMatchesTheNumericOne() {
        assumeTrue(engine != null);
        SolverCTMC solver = new SolverCTMC(oracleModel());
        SolverCTMC.symbolicGeneratorResult sym = solver.getSymbolicGenerator();
        SymEngine.CTMCSolution sol = solver.getSymbolicSolution(engine);

        int n = sym.stateSpace.getNumRows();
        assertEquals(n, sol.pi.size());
        assertEquals(1, sol.nConnComp);

        // Each event filtration was normalized by its own minimum positive
        // rate, so substituting that rate for the event's symbol recovers the
        // model's own generator. The rate has to come from the raw filtration:
        // the one carried in the symbolic result is the normalized one, whose
        // minimum positive entry is 1 by construction.
        Matrix[] rawFilt = new Matrix[sym.symbols.size()];
        for (int e = 0; e < sym.symbols.size(); e++) {
            rawFilt[e] = solver.getGenerator().eventFilt.get(e);
        }
        Map<String, Double> assignment = new HashMap<String, Double>();
        List<String> symbols = sym.activeSymbols();
        for (int k = 0; k < symbols.size(); k++) {
            int e = indexOfSymbol(sym, symbols.get(k));
            assignment.put(symbols.get(k), minPositiveRate(rawFilt[e]));
        }

        double[] symbolicPi;
        try {
            symbolicPi = engine.eval(sol.pi, assignment);
        } catch (IOException e) {
            throw new RuntimeException(e);
        }

        Matrix Q = sym.evalInfGen(assignment);
        Matrix numericPi = ctmc_solve(Q);

        double sum = 0;
        for (int i = 0; i < n; i++) {
            assertEquals(numericPi.get(0, i), symbolicPi[i], 1e-12,
                    "state " + i + " disagrees between the symbolic and the numeric solve");
            sum += symbolicPi[i];
        }
        assertEquals(1.0, sum, 1e-12);
    }

    @Test
    public void symbolicSolutionIsTheKnownRationalFunction() {
        assumeTrue(engine != null);
        SolverCTMC solver = new SolverCTMC(oracleModel());
        SymEngine.CTMCSolution sol = solver.getSymbolicSolution(engine);

        // Every entry is printed over one common denominator: without the
        // primitive normal form the same probability comes back scaled
        // differently from state to state, because a rational constant is a
        // unit in the field and so never enters a gcd.
        String den = null;
        for (String p : sol.pi) {
            int slash = p.indexOf('/');
            assertTrue(slash > 0, "expected a rational function, got " + p);
            String d = p.substring(slash + 1);
            if (den == null) {
                den = d;
            } else {
                assertEquals(den, d, "entries do not share a denominator");
            }
        }
        // The balance product of the oracle, up to the state ordering.
        assertEquals("(6*x1^3 + 6*x1^2*x2 + 3*x1*x2^2 + x2^3)", den);
        assertEquals(sol.pi.size(), sol.num.size());
        assertEquals("6*x1^3 + 6*x1^2*x2 + 3*x1*x2^2 + x2^3", sol.den);
    }

    @Test
    public void evaluationSubstitutesExactly() throws IOException {
        assumeTrue(engine != null);
        Map<String, Double> assignment = new HashMap<String, Double>();
        assignment.put("x1", 1.0);
        assignment.put("x2", 2.0);
        double[] values = engine.eval(
                Arrays.asList("x1/(x1 + x2)", "x2/(x1 + x2)", "x1^2 + x2"), assignment);
        assertEquals(1.0 / 3.0, values[0], 1e-15);
        assertEquals(2.0 / 3.0, values[1], 1e-15);
        assertEquals(3.0, values[2], 1e-15);
    }

    @Test
    public void simplifyAndDifferentiate() throws IOException {
        assumeTrue(engine != null);
        List<String> cancelled = engine.simplify(
                Arrays.asList("(x1^2 - 1)/(x1 - 1)"), "cancel");
        assertEquals("x1 + 1", cancelled.get(0));
        List<String> second = engine.diff(Arrays.asList("x1^3"), "x1", 2);
        assertEquals("6*x1", second.get(0));
    }

    @Test
    public void sensitivityIsExact() throws IOException {
        assumeTrue(engine != null);
        // Two state chain: pi = [x2, x1]/(x1 + x2), reward 1 on the second
        // state, so d(E[r])/dx1 = x2/(x1 + x2)^2, which is 2/9 at x1=1, x2=2.
        String[][] Q = {{"-x1", "x1"}, {"x2", "-x2"}};
        SymEngine.Sensitivity s = engine.ctmcSensitivity(Q,
                Arrays.asList("x1", "x2"), "x1", Arrays.asList("0", "1"));
        assertEquals(2, s.dpi.size());
        Map<String, Double> assignment = new HashMap<String, Double>();
        assignment.put("x1", 1.0);
        assignment.put("x2", 2.0);
        double[] values = engine.eval(Arrays.asList(s.S), assignment);
        assertEquals(2.0 / 9.0, values[0], 1e-12);
    }

    @Test
    public void fluidJacobianIsExact() throws IOException {
        assumeTrue(engine != null);
        SymEngine.FluidODEs odes = engine.fluidODEs(
                Arrays.asList("-a*x + b*y", "a*x - b*y"),
                Arrays.asList("x", "y"),
                Arrays.asList("jacobian", "latex"));
        assertEquals("-a", odes.jacobian[0][0]);
        assertEquals("b", odes.jacobian[0][1]);
        assertEquals("a", odes.jacobian[1][0]);
        assertEquals("-b", odes.jacobian[1][1]);
        assertEquals(2, odes.latex.size());
    }

    @Test
    public void malformedRequestsAreRefusedNotGuessed() {
        assumeTrue(engine != null);
        // A symbol that was never declared must be an error: silently treating
        // it as a new variable would return an expression for a different chain.
        String[][] Q = {{"-x1", "x9"}, {"x1", "-x1"}};
        assertThrows(IOException.class,
                () -> engine.solveCTMC(Q, new ArrayList<String>(Arrays.asList("x1"))));
        // An absorbing chain has no unique stationary distribution.
        String[][] absorbing = {{"-x1", "x1"}, {"0", "0"}};
        assertThrows(IOException.class,
                () -> engine.solveCTMC(absorbing, new ArrayList<String>(Arrays.asList("x1"))));
    }

    private static int indexOfSymbol(SolverCTMC.symbolicGeneratorResult sym, String name) {
        for (int e = 0; e < sym.symbols.size(); e++) {
            if (name.equals(sym.symbols.get(e))) {
                return e;
            }
        }
        throw new IllegalArgumentException("unknown symbol " + name);
    }

    private static double minPositiveRate(Matrix filtration) {
        double min = Double.POSITIVE_INFINITY;
        for (int i = 0; i < filtration.getNumRows(); i++) {
            for (int j = 0; j < filtration.getNumCols(); j++) {
                double v = filtration.get(i, j);
                if (v > 0 && v < min) {
                    min = v;
                }
            }
        }
        return min;
    }
}
