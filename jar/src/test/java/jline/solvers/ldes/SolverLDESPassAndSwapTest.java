/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.ldes;

import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.processes.Exp;
import jline.solvers.NetworkAvgTable;
import jline.solvers.SolverOptions;
import jline.solvers.ctmc.SolverCTMC;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;

import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;

/**
 * LDES tests for pass-and-swap (PAS) / order-independent (OI) queues.
 *
 * <p>A PAS station completes the whole station at the total state-dependent rate
 * mu(c) of its ordered job list (Dorsman and Gardner 2024); the completing
 * position is drawn proportionally to the marginal rate increments and the
 * departing job is selected by the pass-and-swap scan over the swapping graph.
 * The exact stationary distribution is the OI product form, computed by the CTMC
 * solver. JMT does not support PAS, so the LDES output is validated against CTMC
 * (rather than against a JMT t-test like the other LDES tests).
 *
 * <p>For an M/M/K OI queue every in-service job occupies exactly one server, so
 * the per-class utilization conventions of CTMC and LDES coincide and all four
 * metrics agree within simulation tolerance.
 *
 * @see SolverLDES
 * @see jline.examples.java.advanced.PassAndSwapExample
 */
public class SolverLDESPassAndSwapTest {

	private static final int SAMPLES = 1_000_000;
	private static final int SEED = 23000;
	/** Absolute tolerance absorbing Monte-Carlo error at the chosen sample budget. */
	private static final double SIM_TOL = 2e-2;

	/** M/M/K order-independent queue: mu(c) = min(n, K), empty swap graph. */
	private static Network mmkModel() {
		final int K = 2;
		double[] lam = {0.7, 0.5};
		int C = 4;
		Network model = new Network("PASmmk");
		Source src = new Source(model, "Source");
		Queue q = new Queue(model, "PASQueue", SchedStrategy.PAS);
		Sink snk = new Sink(model, "Sink");
		OpenClass c1 = new OpenClass(model, "Class1");
		OpenClass c2 = new OpenClass(model, "Class2");
		src.setArrival(c1, new Exp(lam[0]));
		src.setArrival(c2, new Exp(lam[1]));
		q.setService((Matrix c) -> (double) Math.min(c.getNumCols(), K));
		q.setSwapGraph(new Matrix(2, 2));
		q.setNumberOfServers(K);
		q.setCap(C);
		model.link(Network.serialRouting(src, q, snk));
		return model;
	}

	private static NetworkAvgTable ctmc(Network model, int cutoff) {
		SolverOptions opt = SolverCTMC.defaultOptions();
		opt.cutoff = Matrix.singleton(cutoff);
		return new SolverCTMC(model, opt).getAvgTable();
	}

	private static NetworkAvgTable ldes(Network model) {
		LDESOptions opts = new LDESOptions();
		opts.seed = SEED;
		opts.samples = SAMPLES;
		return new SolverLDES(model, opts).getAvgTable();
	}

	@Test
	@DisplayName("LDES PAS M/M/K matches CTMC (PASQueue rows)")
	public void pasMMKMatchesCTMC() {
		NetworkAvgTable ref = ctmc(mmkModel(), 4);
		NetworkAvgTable sim = ldes(mmkModel());
		// Rows: Source (2 classes) then PASQueue (2 classes). Compare the
		// PASQueue rows only; the Source-row throughput uses the raw arrival
		// rate in LDES vs the post-loss rate in CTMC (a Source-node convention
		// difference unrelated to PAS).
		int nclasses = 2;
		int from = ref.getQLen().size() - nclasses;
		assertVector("QLen", tail(ref.getQLen(), from), tail(sim.getQLen(), from));
		assertVector("Util", tail(ref.getUtil(), from), tail(sim.getUtil(), from));
		assertVector("Tput", tail(ref.getTput(), from), tail(sim.getTput(), from));
		assertVector("RespT", tail(ref.getRespT(), from), tail(sim.getRespT(), from));
	}

	private static List<Double> tail(List<Double> v, int from) {
		return v.subList(from, v.size());
	}

	@Test
	@DisplayName("LDES PAS is deterministic for fixed (seed, samples)")
	public void pasDeterministic() {
		List<Double> a = ldes(mmkModel()).getQLen();
		List<Double> b = ldes(mmkModel()).getQLen();
		assertEquals(a.size(), b.size());
		for (int i = 0; i < a.size(); i++) {
			assertEquals(a.get(i), b.get(i), 0.0,
					"LDES PAS must be deterministic at fixed seed/samples (QLen[" + i + "])");
		}
	}

	private void assertVector(String metric, List<Double> ref, List<Double> sim) {
		assertEquals(ref.size(), sim.size(), metric + ": row count mismatch");
		for (int i = 0; i < ref.size(); i++) {
			assertEquals(ref.get(i), sim.get(i), SIM_TOL,
					"PAS LDES vs CTMC " + metric + "[" + i + "] mismatch");
		}
	}
}
