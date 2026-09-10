/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.ssa;

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
 * SSA tests for pass-and-swap (PAS) / order-independent (OI) queues.
 *
 * <p>The stochastic-simulation engine reuses the shared {@code AfterEventStation}
 * PAS handler. A PAS station departs at the total state-dependent rate mu(c) of
 * its ordered job list (Dorsman and Gardner 2024); in simulation mode the handler
 * collapses the per-position completions of a departing class into a single
 * sampled transition whose rate is the class-r departure rate, so the sum over
 * classes reproduces mu(c). The exact stationary distribution is the OI product
 * form, computed by the CTMC solver, against which the SSA output is validated.
 *
 * @see SolverSSA
 * @see jline.examples.java.advanced.PassAndSwapExample
 */
public class SolverSSAPassAndSwapTest {

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

	private static NetworkAvgTable ssa(Network model) {
		SolverOptions opt = SolverSSA.defaultOptions();
		opt.seed = SEED;
		opt.samples = SAMPLES;
		return new SolverSSA(model, opt).getAvgTable();
	}

	@Test
	@DisplayName("SSA PAS M/M/K matches CTMC (PASQueue rows)")
	public void pasMMKMatchesCTMC() {
		NetworkAvgTable ref = ctmc(mmkModel(), 4);
		NetworkAvgTable sim = ssa(mmkModel());
		// Rows: Source (2 classes) then PASQueue (2 classes). Compare the
		// PASQueue rows only; the Source-row throughput uses the raw arrival
		// rate convention that differs between solvers and is unrelated to PAS.
		int nclasses = 2;
		int from = ref.getQLen().size() - nclasses;
		assertVector("QLen", tail(ref.getQLen(), from), tail(sim.getQLen(), from));
		assertVector("Util", tail(ref.getUtil(), from), tail(sim.getUtil(), from));
		assertVector("Tput", tail(ref.getTput(), from), tail(sim.getTput(), from));
		assertVector("RespT", tail(ref.getRespT(), from), tail(sim.getRespT(), from));
	}

	@Test
	@DisplayName("SSA PAS is deterministic for fixed (seed, samples)")
	public void pasDeterministic() {
		List<Double> a = ssa(mmkModel()).getQLen();
		List<Double> b = ssa(mmkModel()).getQLen();
		assertEquals(a.size(), b.size());
		for (int i = 0; i < a.size(); i++) {
			assertEquals(a.get(i), b.get(i), 0.0,
					"SSA PAS must be deterministic at fixed seed/samples (QLen[" + i + "])");
		}
	}

	private static List<Double> tail(List<Double> v, int from) {
		return v.subList(from, v.size());
	}

	private void assertVector(String metric, List<Double> ref, List<Double> sim) {
		assertEquals(ref.size(), sim.size(), metric + ": row count mismatch");
		for (int i = 0; i < ref.size(); i++) {
			assertEquals(ref.get(i), sim.get(i), SIM_TOL,
					"PAS SSA vs CTMC " + metric + "[" + i + "] mismatch");
		}
	}
}
