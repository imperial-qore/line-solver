package jline.solvers.ssa;

import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.RoutingMatrix;
import jline.lang.constant.SchedStrategy;
import jline.lang.constant.SolverType;
import jline.lang.distributions.APH;
import jline.lang.distributions.Exp;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Router;
import jline.solvers.SolverOptions;
import jline.solvers.mva.SolverMVA;
import jline.solvers.NetworkAvgTable;
import org.junit.jupiter.api.Test;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;

public class SolverSSAClosedExamplesTest {


	@Test
	public void test_example_closedModel_1() throws IllegalAccessException {
		Network model = new Network("example_closedModel_1");

		// Block 1: nodes			
		Delay node1 = new Delay(model, "Delay");
		Queue node2 = new Queue(model, "Queue1", SchedStrategy.FCFS);

		// Block 2: classes
		ClosedClass jobclass1 = new ClosedClass(model, "Class1", 10, node1, 0);
		
		node1.setService(jobclass1, Exp.fitMean(1.000000)); // (Delay,Class1)
		node2.setService(jobclass1, Exp.fitMean(1.500000)); // (Queue1,Class1)

		// Block 3: topology	
		RoutingMatrix routingMatrix = new RoutingMatrix(model,
				Collections.singletonList(jobclass1),
			 Arrays.asList(node1, node2));
	
		routingMatrix.addConnection(jobclass1, jobclass1, node1, node1, 0.700000); // (Delay,Class1) -> (Delay,Class1)
		routingMatrix.addConnection(jobclass1, jobclass1, node1, node2, 0.300000); // (Delay,Class1) -> (Queue1,Class1)
		routingMatrix.addConnection(jobclass1, jobclass1, node2, node1, 1.000000); // (Queue1,Class1) -> (Delay,Class1)

		model.link(routingMatrix);

		SolverOptions options = new SolverOptions(SolverType.SSA);
		options.seed = 23000;
		SolverSSA solver = new SolverSSA(model, options);

		NetworkAvgTable avgTable = solver.getAvgTable();
		//avgTable.print(options);

		List<Double> QLen = avgTable.get(0);
		assertEquals(2.3229476172278685, QLen.get(0), 1e-13);
		assertEquals(7.6765052100706920, QLen.get(1), 1e-13);

		List<Double> Util = avgTable.get(1);
		assertEquals(2.3229476172278685, Util.get(0), 1e-13);
		assertEquals(0.9998816876197125, Util.get(1), 1e-13);

		List<Double> RespT = avgTable.get(2);
		assertEquals(1.0118384964373368, RespT.get(0), 1e-13);
		assertEquals(11.5385682729017560, RespT.get(1), 1e-13);

		List<Double> ResidT = avgTable.get(3);
		assertEquals(1.0118384964373368, ResidT.get(0), 1e-13);
		assertEquals(3.4615704818705266, ResidT.get(1), 1e-13);

		List<Double> Tput = avgTable.get(4);
		assertEquals(2.2900449730000120, Tput.get(0), 1e-13);
		assertEquals(0.6649185948607406, Tput.get(1), 1e-13);
	}
}
