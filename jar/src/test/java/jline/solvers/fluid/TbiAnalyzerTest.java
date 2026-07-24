package jline.solvers.fluid;

import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.constant.SchedStrategy;
import jline.lang.constant.SolverType;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.processes.Erlang;
import jline.lang.processes.Exp;
import jline.solvers.SolverOptions;
import jline.solvers.SolverResult;
import jline.util.matrix.Matrix;
import jline.VerboseLevel;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

public class TbiAnalyzerTest {

    private static Network closedCyclicModel() {
        // Closed cyclic network alternating 3 Delay and 3 Queue (PS) stations,
        // 30 jobs, one class. One delay carries Erlang service, the rest Exp.
        Network model = new Network("tbi_closed_cyclic");

        Delay delay1 = new Delay(model, "Delay1");
        Queue queue1 = new Queue(model, "Queue1", SchedStrategy.PS);
        Delay delay2 = new Delay(model, "Delay2");
        Queue queue2 = new Queue(model, "Queue2", SchedStrategy.PS);
        Delay delay3 = new Delay(model, "Delay3");
        Queue queue3 = new Queue(model, "Queue3", SchedStrategy.PS);

        ClosedClass class1 = new ClosedClass(model, "Class1", 30, delay1, 0);

        delay1.setService(class1, Erlang.fitMeanAndOrder(0.5, 2));
        queue1.setService(class1, new Exp(2.0));
        delay2.setService(class1, new Exp(1.0));
        queue2.setService(class1, new Exp(1.5));
        delay3.setService(class1, new Exp(1.0));
        queue3.setService(class1, new Exp(2.5));

        model.link(Network.serialRouting(delay1, queue1, delay2, queue2, delay3, queue3));
        return model;
    }

    private static Matrix solveQLen(String method) {
        Network model = closedCyclicModel();
        SolverOptions options = new SolverOptions(SolverType.FLUID);
        options.verbose = VerboseLevel.SILENT;
        options.method = method;
        options.iter_max = 200;
        SolverFluid solver = new SolverFluid(model, options);
        solver.runAnalyzer();
        SolverResult result = solver.result;
        // Per-station queue length aggregated over classes.
        Matrix qn = result.QN;
        Matrix perStation = new Matrix(qn.getNumRows(), 1);
        for (int i = 0; i < qn.getNumRows(); i++) {
            double s = 0.0;
            for (int k = 0; k < qn.getNumCols(); k++) {
                s += qn.get(i, k);
            }
            perStation.set(i, 0, s);
        }
        return perStation;
    }

    @Test
    public void testTbiMatchesClosingOnClosedCyclic() {
        Matrix qClosing = solveQLen("closing");
        Matrix qTbi = solveQLen("tbi");

        assertTrue(qClosing.getNumRows() == qTbi.getNumRows(),
                "TBI and closing must report the same number of stations");

        double totalClosing = 0.0;
        double totalTbi = 0.0;
        for (int i = 0; i < qClosing.getNumRows(); i++) {
            double gap = Math.abs(qClosing.get(i, 0) - qTbi.get(i, 0));
            assertTrue(gap < 0.05,
                    "Per-station QLen gap at station " + i + " is " + gap
                            + " (closing=" + qClosing.get(i, 0) + ", tbi=" + qTbi.get(i, 0) + ")");
            totalClosing += qClosing.get(i, 0);
            totalTbi += qTbi.get(i, 0);
        }
        // Closed-population conservation: total fluid mass must equal 30.
        assertTrue(Math.abs(totalClosing - 30.0) < 0.5, "closing total mass " + totalClosing);
        assertTrue(Math.abs(totalTbi - 30.0) < 0.5, "tbi total mass " + totalTbi);
    }

    @Test
    public void testTbiRejectsOpenModel() {
        Network model = new Network("tbi_open");
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue1", SchedStrategy.PS);
        Sink sink = new Sink(model, "Sink");

        OpenClass oclass = new OpenClass(model, "OpenClass1");
        source.setArrival(oclass, new Exp(1.0));
        queue.setService(oclass, new Exp(2.0));

        model.link(Network.serialRouting(source, queue, sink));

        SolverOptions options = new SolverOptions(SolverType.FLUID);
        options.verbose = VerboseLevel.SILENT;
        options.method = "tbi";
        SolverFluid solver = new SolverFluid(model, options);
        assertThrows(RuntimeException.class, solver::runAnalyzer);
    }
}
