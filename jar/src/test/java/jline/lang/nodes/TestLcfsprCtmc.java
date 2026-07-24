package jline.lang.nodes;

import jline.lang.*;
import jline.lang.constant.*;
import jline.lang.nodes.*;
import jline.lang.processes.*;
import jline.solvers.ctmc.*;
import jline.solvers.*;
import jline.util.matrix.*;
import jline.VerboseLevel;

public class TestLcfsprCtmc {
    public static void main(String[] args) {
        System.out.println("=== Testing cqn_bcmp_theorem_lcfspr model with CTMC ===");
        
        // Create the model
        Network model = new Network("myModel");
        
        // Block 1: nodes
        Delay node1 = new Delay(model, "Delay");
        Queue node2 = new Queue(model, "Queue1", SchedStrategy.LCFSPR);
        
        // Block 2: classes
        ClosedClass jobclass1 = new ClosedClass(model, "Class1", 2, node1, 0);
        ClosedClass jobclass2 = new ClosedClass(model, "Class2", 2, node1, 0);
        
        node1.setService(jobclass1, new Erlang(3, 2));
        node1.setService(jobclass2, new HyperExp(0.5, 3.0, 10.0));
        node2.setService(jobclass1, Exp.fitMean(1.00));
        node2.setService(jobclass2, Exp.fitMean(1.00));
        
        // Block 3: topology
        RoutingMatrix routingMatrix = model.initRoutingMatrix();
        routingMatrix.set(jobclass1, jobclass1, node1, node2, 1.00);
        routingMatrix.set(jobclass1, jobclass1, node2, node1, 1.00);
        routingMatrix.set(jobclass2, jobclass2, node1, node2, 1.00);
        routingMatrix.set(jobclass2, jobclass2, node2, node1, 1.00);
        model.link(routingMatrix);
        
        System.out.println("Model created successfully");
        System.out.println("Nodes: " + model.getNumberOfNodes());
        System.out.println("Classes: " + model.getNumberOfClasses());
        
        // Solve with CTMC
        SolverCTMC solver = new SolverCTMC(model);
        solver.options.verbose = VerboseLevel.DEBUG;
        
        System.out.println("\nRunning CTMC solver...");
        NetworkAvgTable avgTable = solver.getAvgTable();
        
        System.out.println("\n=== CTMC Results ===");
        avgTable.print();
        
        System.out.println("\n=== Expected Values (from MATLAB) ===");
        System.out.println("Delay, Class1: QLen=0.30861, Tput=0.46291");
        System.out.println("Delay, Class2: QLen=0.11625, Tput=0.53653");
        System.out.println("Queue1, Class1: QLen=1.69139, Tput=0.46291");
        System.out.println("Queue1, Class2: QLen=1.88375, Tput=0.53653");
    }
}
