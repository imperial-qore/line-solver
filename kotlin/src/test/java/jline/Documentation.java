package jline;

import jline.solvers.jmt.SolverJMT;
import org.junit.jupiter.api.Test;

import jline.lang.*;
import jline.lang.constant.*;
import jline.lang.nodes.*;
import jline.lang.distributions.*;
import jline.lang.nodes.Queue;

import jline.solvers.SolverOptions;
import jline.solvers.NetworkAvgTable;

import static org.junit.jupiter.api.Assertions.assertEquals;

public class Documentation {
    /*
    @Test
    public void ex_chap2_1() {

        Network model = new Network("M/M/1");
        // Block 1: nodes
        Source mySource = new Source(model, "mySource");
        Queue myQueue = new Queue(model, "myQueue", SchedStrategy.FCFS);
        Sink mySink = new Sink(model, "mySink");
        // Block 2: classes
        OpenClass myClass = new OpenClass(model, "myClass", 0);
        mySource.setArrival(myClass, new Exp(1.0)); // (mySource,myClass)
        myQueue.setService(myClass, new Exp(2.0)); // (myQueue,myClass)
        // Block 3: topology
        model.link(model.serialRouting(mySource, myQueue, mySink));
        // Block 4: solution
        SolverOptions options = new SolverOptions(SolverType.Fluid);
        options.seed = 23000;
        SolverJMT solver = new SolverJMT(model, options);
        NetworkAvgTable avgTable = solver.getAvgTable();
        avgTable.print();
    }

     */
}
