package jline.solvers.fluid;

import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.RoutingMatrix;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.processes.Erlang;
import jline.lang.processes.Exp;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Tests for SolverFluid.exportODEs, the LaTeX export of the mean-field ODE
 * system. Expected strings are cross-validated against the MATLAB
 * SolverFLD.exportODEs output, which is numerically verified against the
 * ODE right-hand sides integrated by the solver.
 */
public class SolverFluidExportODEsTest {

    private static Network buildClosedExpModel() {
        Network model = new Network("D");
        Delay delay = new Delay(model, "Delay1");
        Queue queue = new Queue(model, "Queue1", SchedStrategy.PS);
        ClosedClass c1 = new ClosedClass(model, "C1", 4, delay);
        delay.setService(c1, new Exp(1.0));
        queue.setService(c1, new Exp(2.0));
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(c1, Network.serialRouting(delay, queue));
        model.link(P);
        return model;
    }

    private static Network buildDpsModel() {
        Network model = new Network("B");
        Delay delay = new Delay(model, "Delay1");
        Queue qdps = new Queue(model, "QueueDPS", SchedStrategy.DPS);
        ClosedClass b1 = new ClosedClass(model, "B1", 2, delay);
        ClosedClass b2 = new ClosedClass(model, "B2", 3, delay);
        delay.setService(b1, new Exp(1.0));
        delay.setService(b2, new Exp(0.5));
        qdps.setService(b1, Erlang.fitMeanAndOrder(1.0, 2), 2.0);
        qdps.setService(b2, new Exp(1.0), 1.0);
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(b1, Network.serialRouting(delay, qdps));
        P.set(b2, Network.serialRouting(delay, qdps));
        model.link(P);
        return model;
    }

    private static Network buildOpenModel() {
        Network model = new Network("C");
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue1", SchedStrategy.PS);
        Sink sink = new Sink(model, "Sink");
        OpenClass oc = new OpenClass(model, "OC");
        source.setArrival(oc, new Exp(0.5));
        queue.setService(oc, new Exp(1.0));
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(oc, Network.serialRouting(source, queue, sink));
        model.link(P);
        return model;
    }

    @Test
    public void testClosingScalarExportClosedModel() {
        SolverFluid solver = new SolverFluid(buildClosedExpModel());
        solver.options.method = "closing";
        String tex = solver.exportODEs("", "scalar");
        assertTrue(tex.contains("% method: closing"));
        assertTrue(tex.contains("% form: dx/dt = J*r(x)"));
        assertTrue(tex.contains("% nstates: 2"));
        assertTrue(tex.contains("% nevents: 2"));
        assertTrue(tex.contains("% STATE 1 station=Delay1 class=C1 phase=1"));
        assertTrue(tex.contains("% STATE 2 station=Queue1 class=C1 phase=1"));
        assertTrue(tex.contains("n_{2}(\\mathbf{x}) &= x_{2}"));
        assertTrue(tex.contains("g_{2}(\\mathbf{x}) &= \\frac{\\min(n_{2}(\\mathbf{x}),\\, 1)}{n_{2}(\\mathbf{x})}"));
        assertTrue(tex.contains("\\frac{\\mathrm{d}x_{1}}{\\mathrm{d}t} &= -x_{1} + 2\\,x_{2}\\,g_{2}(\\mathbf{x})\\\\"));
        assertTrue(tex.contains("\\frac{\\mathrm{d}x_{2}}{\\mathrm{d}t} &= x_{1} - 2\\,x_{2}\\,g_{2}(\\mathbf{x})"));
        assertTrue(tex.contains("\\mathbf{x}(0) = \\begin{pmatrix} 4 & 0 \\end{pmatrix}^{\\top}"));
    }

    @Test
    public void testClosingScalarExportDpsModel() {
        SolverFluid solver = new SolverFluid(buildDpsModel());
        solver.options.method = "closing";
        String tex = solver.exportODEs("", "scalar");
        // DPS: weights normalized to (2/3, 1/3) and folded into the coefficients;
        // the shares divide the capacity min(n_2, S_2), with no additive seed
        assertTrue(tex.contains("\\tilde{n}_{2}(\\mathbf{x}) &= 0.66666667\\,(x_{3} + x_{4}) + 0.33333333\\,(x_{5})"));
        assertTrue(tex.contains("g_{2}(\\mathbf{x}) &= \\frac{\\min(n_{2}(\\mathbf{x}),\\, 1)}{\\tilde{n}_{2}(\\mathbf{x})}"));
        assertTrue(tex.contains("\\frac{\\mathrm{d}x_{1}}{\\mathrm{d}t} &= -x_{1} + 1.3333333\\,x_{4}\\,g_{2}(\\mathbf{x})\\\\"));
        assertTrue(tex.contains("\\frac{\\mathrm{d}x_{5}}{\\mathrm{d}t} &= 0.5\\,x_{2} - 0.33333333\\,x_{5}\\,g_{2}(\\mathbf{x})"));
        assertTrue(tex.contains("\\mathbf{x}(0) = \\begin{pmatrix} 2 & 3 & 0 & 0 & 0 \\end{pmatrix}^{\\top}"));
    }

    @Test
    public void testStatedepAndSoftminExports() {
        Network model = buildDpsModel();
        SolverFluid solver = new SolverFluid(model);
        solver.options.method = "statedep";
        String tex = solver.exportODEs("", "scalar");
        // piecewise DPS factor with per-class definition
        assertTrue(tex.contains("g_{2,1}(\\mathbf{x}) &= \\begin{cases} 1 & n_{2}(\\mathbf{x}) \\le 1\\\\"));
        assertTrue(tex.contains("x_{3}\\,g_{2,1}(\\mathbf{x})"));

        SolverFluid solver2 = new SolverFluid(buildDpsModel());
        solver2.options.method = "softmin";
        String tex2 = solver2.exportODEs("", "scalar");
        assertTrue(tex2.contains("% method: softmin"));
    }

    @Test
    public void testMatrixNotationOpenModel() {
        SolverFluid solver = new SolverFluid(buildOpenModel());
        solver.options.method = "matrix";
        String tex = solver.exportODEs("", "matrix");
        assertTrue(tex.contains("% form: dx/dt = W^T*theta(x) + lambda"));
        assertTrue(tex.contains("\\frac{\\mathrm{d}\\mathbf{x}}{\\mathrm{d}t} = W^{\\top}\\,\\theta(\\mathbf{x}) + \\boldsymbol{\\lambda}"));
        // Source state theta = 0, arrivals via lambda
        assertTrue(tex.contains("\\theta(\\mathbf{x}) = \\begin{bmatrix} 0 \\\\ x_{2}\\,g_{2}(\\mathbf{x}) \\end{bmatrix}"));
        assertTrue(tex.contains("\\boldsymbol{\\lambda} = \\begin{pmatrix} 0 & 0.5 \\end{pmatrix}^{\\top}"));
    }

    @Test
    public void testUnsupportedMethodAndNotationRejected() {
        SolverFluid solver = new SolverFluid(buildClosedExpModel());
        solver.options.method = "mfq";
        assertThrows(RuntimeException.class, new org.junit.jupiter.api.function.Executable() {
            @Override
            public void execute() {
                solver.exportODEs("", "scalar");
            }
        });
        SolverFluid solver2 = new SolverFluid(buildClosedExpModel());
        solver2.options.method = "closing";
        assertThrows(RuntimeException.class, new org.junit.jupiter.api.function.Executable() {
            @Override
            public void execute() {
                solver2.exportODEs("", "vector");
            }
        });
        // statedep does not support open models
        SolverFluid solver3 = new SolverFluid(buildOpenModel());
        solver3.options.method = "statedep";
        assertThrows(RuntimeException.class, new org.junit.jupiter.api.function.Executable() {
            @Override
            public void execute() {
                solver3.exportODEs("", "scalar");
            }
        });
    }

    @Test
    public void testScalarAndMatrixShareStructure() {
        SolverFluid solver = new SolverFluid(buildClosedExpModel());
        solver.options.method = "closing";
        String scalar = solver.exportODEs("", "scalar");
        String matrix = solver.exportODEs("", "matrix");
        assertEquals(headerOf(scalar, "% notation: scalar"), headerOf(matrix, "% notation: matrix"));
        assertTrue(matrix.contains("\\frac{\\mathrm{d}\\mathbf{x}}{\\mathrm{d}t} = J\\,r(\\mathbf{x})"));
        assertTrue(matrix.contains("r_{1}(\\mathbf{x}) &= x_{1}\\\\"));
        assertTrue(matrix.contains("r_{2}(\\mathbf{x}) &= 2\\,x_{2}\\,g_{2}(\\mathbf{x})"));
        assertFalse(scalar.contains("stoichiometry"));
    }

    private static String headerOf(String tex, String notationLine) {
        // headers must agree except for the notation line
        StringBuilder sb = new StringBuilder();
        String[] lines = tex.split("\n");
        for (int i = 0; i < lines.length; i++) {
            if (!lines[i].startsWith("%")) {
                break;
            }
            if (lines[i].startsWith("% notation:")) {
                assertEquals(notationLine, lines[i]);
                continue;
            }
            sb.append(lines[i]).append("\n");
        }
        return sb.toString();
    }
}
