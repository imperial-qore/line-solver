/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.ldes;

import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.RoutingMatrix;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.processes.Exp;
import jline.lang.processes.HyperExp;
import jline.solvers.NetworkAvgTable;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;

import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;

/**
 * Java/Kotlin parity tests for the LDES solver.
 *
 * <p>The LDES simulator is deterministic for a fixed (seed, samples) pair, and
 * the Java ({@code jline.jar}) and Kotlin ({@code jline-kotlin.jar}) builds are
 * line-for-line mirrors that share the same SSJ-based event engine and seeding,
 * so for identical models they produce <em>byte-identical</em> metrics.
 *
 * <p>Each test runs LDES on a small, self-contained model at {@code seed=23000}
 * and {@code samples=100000} and asserts the QLen/Util/Tput/RespT vectors equal
 * the golden values captured from the canonical Java run. The <strong>same</strong>
 * golden values are asserted by the mirrored {@code SolverLDESParityTest.kt} in
 * the Kotlin module: if either implementation drifts, its parity test fails.
 *
 * <p>The tolerance only guards against floating-point reassociation across JVMs;
 * any genuine Java/Kotlin divergence (different event order, RNG, or algorithm)
 * compounds far beyond it.
 *
 * @see SolverLDES
 * @see SolverLDESTestFixtures
 */
public class SolverLDESParityTest extends SolverLDESTestFixtures {

	/** Fixed sample budget for the parity golden values. */
	private static final int PARITY_SAMPLES = 100000;

	/** Absolute tolerance: deterministic match expected, this only absorbs FP reassociation. */
	private static final double PARITY_TOL = 1e-6;

	// ---------------------------------------------------------------- models

	/** Closed repairmen: infinite-server Delay + FCFS Queue, one class, 10 jobs. */
	private static Network closedRepairmen() {
		Network model = new Network("closedRepairmen");
		Delay node1 = new Delay(model, "Delay");
		Queue node2 = new Queue(model, "Queue1", SchedStrategy.FCFS);
		ClosedClass jobclass1 = new ClosedClass(model, "Class1", 10, node1, 0);
		node1.setService(jobclass1, Exp.fitMean(1.00));
		node2.setService(jobclass1, Exp.fitMean(1.50));
		RoutingMatrix P = model.initRoutingMatrix();
		P.set(jobclass1, jobclass1, node1, node1, 0.70);
		P.set(jobclass1, jobclass1, node1, node2, 0.30);
		P.set(jobclass1, jobclass1, node2, node1, 1.00);
		model.link(P);
		return model;
	}

	/** Open tandem: Source + Delay (HyperExp) + FCFS Queue + Sink, one open class. */
	private static Network openBasic() {
		Network model = new Network("openBasic");
		Delay node1 = new Delay(model, "Delay");
		Queue node2 = new Queue(model, "Queue1", SchedStrategy.FCFS);
		Source node3 = new Source(model, "Source");
		Sink node4 = new Sink(model, "Sink");
		OpenClass jobclass1 = new OpenClass(model, "Class1", 0);
		node1.setService(jobclass1, new HyperExp(0.5, 3.0, 10.0));
		node2.setService(jobclass1, new Exp(1.0));
		node3.setArrival(jobclass1, new Exp(0.1));
		RoutingMatrix P = model.initRoutingMatrix();
		P.set(jobclass1, jobclass1, new Matrix("[0,1,0,0; 0,0,0,1; 1,0,0,0; 0,0,0,0]"));
		model.link(P);
		return model;
	}

	/** Closed multiclass: Delay + PS Queue, two classes (2 and 1 jobs). */
	private static Network closedMulticlassPS() {
		Network model = new Network("closedMulticlassPS");
		Delay d = new Delay(model, "Delay");
		Queue q = new Queue(model, "Queue1", SchedStrategy.PS);
		ClosedClass c1 = new ClosedClass(model, "Class1", 2, d, 0);
		ClosedClass c2 = new ClosedClass(model, "Class2", 1, d, 0);
		d.setService(c1, Exp.fitMean(1.0));
		d.setService(c2, Exp.fitMean(2.0));
		q.setService(c1, Exp.fitMean(1.5));
		q.setService(c2, Exp.fitMean(0.8));
		RoutingMatrix P = model.initRoutingMatrix();
		P.set(c1, c1, d, q, 1.0);
		P.set(c1, c1, q, d, 1.0);
		P.set(c2, c2, d, q, 1.0);
		P.set(c2, c2, q, d, 1.0);
		model.link(P);
		return model;
	}

	/** Closed multiserver: Delay + 2-server FCFS Queue, one class, 5 jobs. */
	private static Network closedMultiserver() {
		Network model = new Network("closedMultiserver");
		Delay d = new Delay(model, "Delay");
		Queue q = new Queue(model, "Queue1", SchedStrategy.FCFS);
		q.setNumberOfServers(2);
		ClosedClass c1 = new ClosedClass(model, "Class1", 5, d, 0);
		d.setService(c1, Exp.fitMean(1.0));
		q.setService(c1, Exp.fitMean(1.5));
		RoutingMatrix P = model.initRoutingMatrix();
		P.set(c1, c1, d, q, 1.0);
		P.set(c1, c1, q, d, 1.0);
		model.link(P);
		return model;
	}

	// ------------------------------------------------------------- assertion

	private void assertParity(String name, Network model,
			double[] qlen, double[] util, double[] tput, double[] respt) {
		LDESOptions opts = createTestOptions(PARITY_SAMPLES);
		SolverLDES solver = new SolverLDES(model, opts);
		NetworkAvgTable t = solver.getAvgTable();
		assertVector(name, "QLen", qlen, t.getQLen());
		assertVector(name, "Util", util, t.getUtil());
		assertVector(name, "Tput", tput, t.getTput());
		assertVector(name, "RespT", respt, t.getRespT());
	}

	private void assertVector(String name, String metric, double[] golden, List<Double> actual) {
		assertEquals(golden.length, actual.size(),
				name + ": " + metric + " row count mismatch (Java/Kotlin parity)");
		for (int i = 0; i < golden.length; i++) {
			assertEquals(golden[i], actual.get(i), PARITY_TOL,
					name + ": " + metric + "[" + i + "] Java/Kotlin LDES parity mismatch");
		}
	}

	// ------------------------------------------------------------------ tests

	@Test
	@DisplayName("LDES parity: closed repairmen (FCFS, single class)")
	public void parityClosedRepairmen() {
		assertParity("closedRepairmen", closedRepairmen(),
				new double[]{2.21852361727186, 7.781476382728091},
				new double[]{2.2141866895358495, 0.9999455981753816},
				new double[]{2.2141866895358495, 0.6651423999204868},
				new double[]{1.0025783708598637, 11.679842689056642});
	}

	@Test
	@DisplayName("LDES parity: open tandem (HyperExp delay, FCFS queue)")
	public void parityOpenBasic() {
		assertParity("openBasic", openBasic(),
				new double[]{0.021545948652226407, 0.11049675717207369, 0.0},
				new double[]{0.021709189403355814, 0.09968098339610598, 0.0},
				new double[]{0.10019625878471913, 0.10019625878471913, 0.1},
				new double[]{0.21500903101137978, 1.1028696263184417, 0.0});
	}

	@Test
	@DisplayName("LDES parity: closed multiclass PS (two classes)")
	public void parityClosedMulticlassPS() {
		assertParity("closedMulticlassPS", closedMulticlassPS(),
				new double[]{0.499304632415569, 0.5092223728947529, 1.5006953675844212, 0.4907776271052469},
				new double[]{0.49347584257729576, 0.5100643991875738, 0.7339629098835129, 0.20402037703624515},
				new double[]{0.49347584257729576, 0.2550321995937869, 0.49346085743130635, 0.25504718473977633},
				new double[]{1.0121617987071256, 1.9962920287197428, 3.0405861890553547, 1.9241679650669947});
	}

	@Test
	@DisplayName("LDES parity: closed multiserver (2 servers, FCFS)")
	public void parityClosedMultiserver() {
		assertParity("closedMultiserver", closedMultiserver(),
				new double[]{1.3045512077856398, 3.695448792214329},
				new double[]{1.296950764907132, 0.9779835211003793},
				new double[]{1.296950764907132, 1.296898835986193},
				new double[]{1.0062645062177846, 2.848432870789654});
	}

}
