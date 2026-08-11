/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers.ldes;

import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.constant.SchedStrategy;
import jline.lang.constant.ServerType;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.processes.Exp;
import jline.lang.processes.NHPP;
import jline.solvers.NetworkAvgTable;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Correctness tests for the cyclic NHPP rate schedule in the LDES engine.
 *
 * <p>These assert analytic invariants rather than captured golden values, so
 * they stay meaningful without regeneration. The invariants are the ones that
 * exposed the original defects:
 *
 * <ul>
 * <li><b>Flow balance.</b> In an open, infinite-capacity, stable queue no jobs
 * are lost, so the departure rate equals the arrival rate. The engine used to
 * sample an interarrival from the rate current at the previous arrival and keep
 * that sample across a phase boundary, which starved every fast phase; the
 * error grew as boundaries became more frequent (21% low at short durations)
 * and vanished as they became rare, so a single-schedule test could not have
 * caught it. The schedule sweep below is therefore load-bearing.</li>
 *
 * <li><b>Utilization law.</b> U = lambda * E[S] holds for any stable open
 * single-server queue irrespective of how bursty the arrival stream is.</li>
 *
 * <li><b>Switching limits.</b> A time-varying service rate mu(t) must converge
 * to the time-average-rate exponential as the phase durations vanish, and to
 * the pointwise-stationary (PSA) mixture of per-phase exponentials as they grow
 * without bound. These bracket the truth and need no simulator to compute.</li>
 * </ul>
 */
public class SolverLDESNHPPTest extends SolverLDESTestFixtures {

	private static final int SAMPLES = 400000;
	/** Simulation tolerance: the metrics below are sample means. */
	private static final double SIM_TOL = 0.01;

	// ------------------------------------------------------------- fixtures

	/**
	 * A cyclic NHPP holding {@code rates[i]} for {@code durations[i]} time units,
	 * i.e. {@code NHPP(cumsum([0, durations]), rates, true)}.
	 */
	private static NHPP cyclic(double[] rates, double[] durations) {
		double[] breakpoints = new double[durations.length + 1];
		breakpoints[0] = 0.0;
		for (int i = 0; i < durations.length; i++) {
			breakpoints[i + 1] = breakpoints[i] + durations[i];
		}
		return new NHPP(breakpoints, rates, true);
	}

	/** Open M/M/1 with a cyclic NHPP arrival schedule and Exp(mu) service. */
	private static Network arrivalSchedule(double[] rates, double[] durations, double mu) {
		Network model = new Network("cyclic-arrival");
		Source source = new Source(model, "Source");
		Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
		Sink sink = new Sink(model, "Sink");
		OpenClass oclass = new OpenClass(model, "OpenClass", 0);
		source.setArrival(oclass, cyclic(rates, durations));
		queue.setService(oclass, new Exp(mu));
		model.link(model.serialRouting(source, queue, sink));
		return model;
	}

	/** Open M/Mt/1 with Poisson arrivals and a cyclic NHPP service schedule. */
	private static Network serviceSchedule(double lambda, double[] rates,
			double[] durations, SchedStrategy sched) {
		Network model = new Network("cyclic-service");
		Source source = new Source(model, "Source");
		Queue queue = new Queue(model, "Queue", sched);
		Sink sink = new Sink(model, "Sink");
		OpenClass oclass = new OpenClass(model, "OpenClass", 0);
		source.setArrival(oclass, new Exp(lambda));
		queue.setService(oclass, cyclic(rates, durations));
		model.link(model.serialRouting(source, queue, sink));
		return model;
	}

	private static NetworkAvgTable solve(Network model) {
		return new SolverLDES(model, createTestOptions(SAMPLES)).getAvgTable();
	}

	private static double timeAverageRate(double[] rates, double[] durations) {
		double mass = 0.0;
		double span = 0.0;
		for (int i = 0; i < rates.length; i++) {
			mass += rates[i] * durations[i];
			span += durations[i];
		}
		return mass / span;
	}

	// ------------------------------------------------------- arrival lambda(t)

	@Test
	@DisplayName("NHPP arrivals: flow balance and utilization law hold at every phase length")
	public void testArrivalScheduleInvariants() {
		double[] rates = {2.0, 8.0, 4.0};
		double mu = 100.0;
		double lambda = timeAverageRate(rates, new double[]{3.0, 1.0, 2.0});
		// Phase lengths spanning four orders of magnitude relative to the mean
		// interarrival time. The pre-fix engine was 21% low at the short end and
		// within 0.3% at the long end, so only the sweep discriminates.
		double[] scales = {0.1, 1.0, 10.0, 100.0};
		for (int s = 0; s < scales.length; s++) {
			double c = scales[s];
			double[] durations = {3.0 * c, 1.0 * c, 2.0 * c};
			NetworkAvgTable t = solve(arrivalSchedule(rates, durations, mu));
			String at = " (phase scale " + c + ")";
			// Row 1 is the Queue; row 0 is the Source.
			double tput = t.getTput().get(1);
			double util = t.getUtil().get(1);
			assertEquals(lambda, tput, SIM_TOL * lambda,
					"departure rate must equal the time-average arrival rate" + at);
			assertEquals(lambda / mu, util, SIM_TOL * (lambda / mu),
					"utilization law U = lambda*S must hold" + at);
		}
	}

	@Test
	@DisplayName("NHPP arrivals: a single-phase schedule reproduces Exp arrivals")
	public void testArrivalDegeneratesToExp() {
		double lambda = 3.5;
		double mu = 10.0;
		NetworkAvgTable cyclic =
				solve(arrivalSchedule(new double[]{lambda}, new double[]{1.0}, mu));
		NetworkAvgTable exp = solve(arrivalSchedule(
				new double[]{lambda, lambda}, new double[]{1.0, 2.0}, mu));
		// A constant schedule is a homogeneous Poisson process however it is cut
		// into phases, so the phase structure must not affect the result.
		assertEquals(cyclic.getQLen().get(1), exp.getQLen().get(1),
				SIM_TOL * Math.max(1e-9, cyclic.getQLen().get(1)),
				"a constant rate must not depend on how the cycle is partitioned");
	}

	// ------------------------------------------------------- service mu(t)

	@Test
	@DisplayName("NHPP service: a constant schedule reproduces Exp service")
	public void testServiceDegeneratesToExp() {
		double lambda = 1.0;
		double mu = 11.0;
		NetworkAvgTable cyclic =
				solve(serviceSchedule(lambda, new double[]{mu, mu},
						new double[]{1.0, 1.0}, SchedStrategy.FCFS));
		double util = cyclic.getUtil().get(1);
		// Regression: an NHPP service used to reach sn.proc as {NaN, NaN}
		// and simulate as zero service time, reporting QLen = Util = 0.
		assertTrue(util > 0.0, "service schedule must not simulate as zero service time");
		assertEquals(lambda / mu, util, SIM_TOL * (lambda / mu),
				"utilization law must hold for a constant service schedule");
	}

	@Test
	@DisplayName("NHPP service: converges to the fast- and slow-switching limits")
	public void testServiceSwitchingLimits() {
		double lambda = 1.0;
		double[] rates = {2.0, 20.0};

		// Fast switching: mu(t) self-averages within one service, so the station
		// behaves as Exp(timeAverageRate) and U -> lambda/rateBar.
		double rateBar = timeAverageRate(rates, new double[]{1.0, 1.0});
		NetworkAvgTable fast = solve(serviceSchedule(lambda, rates,
				new double[]{1e-4, 1e-4}, SchedStrategy.FCFS));
		assertEquals(lambda / rateBar, fast.getUtil().get(1), SIM_TOL * (lambda / rateBar),
				"fast-switching limit must match the time-average-rate exponential");

		// Slow switching (PSA): a service sees a single frozen rate, so
		// E[S] -> mean of the per-phase means and U -> lambda*E[S].
		double psaMeanService = 0.5 * (1.0 / rates[0]) + 0.5 * (1.0 / rates[1]);
		NetworkAvgTable slow = solve(serviceSchedule(lambda, rates,
				new double[]{1e3, 1e3}, SchedStrategy.FCFS));
		assertEquals(lambda * psaMeanService, slow.getUtil().get(1),
				SIM_TOL * (lambda * psaMeanService),
				"slow-switching limit must match the pointwise-stationary mixture");

		// The two limits must be genuinely different, else the test is vacuous.
		assertTrue(psaMeanService > 1.5 / rateBar,
				"fixture must separate the two limits");
	}

	/**
	 * Disciplines under which a time-varying service rate is now simulated
	 * holistically via the operational-time transform. FCFS is the reference;
	 * PS shares the time-varying rate, and FCFSPR/LCFSPR/SRPT interrupt and
	 * resume service, so each exercises a different residual path.
	 */
	private static final SchedStrategy[] SCHEDULES = {
		SchedStrategy.FCFS, SchedStrategy.PS, SchedStrategy.FCFSPR,
		SchedStrategy.LCFSPR, SchedStrategy.SRPT
	};

	@Test
	@DisplayName("NHPP service: flow balance holds under every scheduling discipline")
	public void testServiceScheduleFlowBalanceAcrossDisciplines() {
		double lambda = 1.0;
		double[] rates = {2.0, 20.0};
		// A stable, infinite-capacity, work-conserving station loses no jobs, so
		// its departure rate equals lambda whatever the discipline. An interrupted
		// service resumed with the wrong work would change the busy period and so
		// the achieved throughput; flow balance is the discipline-free check.
		for (int s = 0; s < SCHEDULES.length; s++) {
			SchedStrategy sched = SCHEDULES[s];
			NetworkAvgTable t = solve(serviceSchedule(lambda, rates,
					new double[]{1.0, 1.0}, sched));
			assertEquals(lambda, t.getTput().get(1), SIM_TOL * lambda,
					"departure rate must equal lambda under " + sched);
		}
	}

	@Test
	@DisplayName("NHPP service: fast-switching limit is reproduced under every discipline")
	public void testServiceScheduleFastLimitAcrossDisciplines() {
		double lambda = 1.0;
		double[] rates = {2.0, 20.0};
		double rateBar = timeAverageRate(rates, new double[]{1.0, 1.0});
		// mu(t) self-averages within one service, so the station behaves as
		// Exp(rateBar) and U -> lambda/rateBar irrespective of the discipline. A
		// residual carried in wall time rather than cumulative intensity would
		// miss this limit under PS and the preemptive policies.
		for (int s = 0; s < SCHEDULES.length; s++) {
			SchedStrategy sched = SCHEDULES[s];
			NetworkAvgTable t = solve(serviceSchedule(lambda, rates,
					new double[]{1e-4, 1e-4}, sched));
			assertEquals(lambda / rateBar, t.getUtil().get(1),
					SIM_TOL * (lambda / rateBar),
					"fast-switching utilization must match Exp(rateBar) under " + sched);
		}
	}

	@Test
	@DisplayName("NHPP service: slow-switching (PSA) limit is reproduced under every discipline")
	public void testServiceScheduleSlowLimitAcrossDisciplines() {
		double lambda = 1.0;
		double[] rates = {2.0, 20.0};
		// Slow switching: a whole service sees a single frozen rate, so E[S] ->
		// mean of the per-phase means and U -> lambda*E[S]. Utilization of a
		// work-conserving single-server station is discipline-invariant, so every
		// discipline must land on the same pointwise-stationary value.
		double psaMeanService = 0.5 * (1.0 / rates[0]) + 0.5 * (1.0 / rates[1]);
		for (int s = 0; s < SCHEDULES.length; s++) {
			SchedStrategy sched = SCHEDULES[s];
			NetworkAvgTable t = solve(serviceSchedule(lambda, rates,
					new double[]{1e3, 1e3}, sched));
			assertEquals(lambda * psaMeanService, t.getUtil().get(1),
					SIM_TOL * (lambda * psaMeanService),
					"slow-switching utilization must match the PSA mixture under " + sched);
		}
	}

	@Test
	@DisplayName("NHPP service: utilization is discipline-invariant at intermediate phase scales")
	public void testServiceScheduleUtilizationDisciplineInvariant() {
		double lambda = 1.0;
		double[] rates = {2.0, 20.0};
		// For a work-conserving single-server queue the busy fraction depends only
		// on the offered-work process, not the discipline: unfinished operational
		// time drains at rate mu(t) whenever the server is busy, however jobs are
		// ordered or shared. PS and the preemptive policies must reproduce the
		// FCFS utilization at a phase scale where neither switching limit applies.
		double[] durations = {0.3, 0.3};
		double reference = solve(serviceSchedule(lambda, rates, durations,
				SchedStrategy.FCFS)).getUtil().get(1);
		for (int s = 0; s < SCHEDULES.length; s++) {
			SchedStrategy sched = SCHEDULES[s];
			NetworkAvgTable t = solve(serviceSchedule(lambda, rates, durations, sched));
			assertEquals(reference, t.getUtil().get(1), 2.0 * SIM_TOL * reference,
					"utilization must be discipline-invariant under " + sched);
		}
	}

	// -------------------------------------------------- unsupported combinations
	//
	// Neither restriction below is expressible in the feature set, which carries a
	// single boolean "NHPP" entry, so both surface only as a run-time exception.
	// These tests pin the messages so a refactor cannot silently turn either into
	// a wrong sample path.

	@Test
	@DisplayName("NHPP service with heterogeneous server types is rejected, not silently mis-sampled")
	public void testHeteroServerNhppServiceRejected() {
		Network model = new Network("hetero-nhpp-service");
		Source source = new Source(model, "Source");
		Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
		Sink sink = new Sink(model, "Sink");
		OpenClass oclass = new OpenClass(model, "OpenClass", 0);
		source.setArrival(oclass, new Exp(1.0));
		ServerType fast = new ServerType("Fast", 1);
		ServerType slow = new ServerType("Slow", 1);
		queue.addServerType(fast);
		queue.addServerType(slow);
		NHPP schedule = cyclic(new double[]{2.0, 20.0}, new double[]{0.3, 0.3});
		queue.setService(oclass, schedule);
		queue.setService(oclass, fast, schedule);
		queue.setService(oclass, slow, new Exp(4.0));
		model.link(model.serialRouting(source, queue, sink));

		// generateHeteroServiceTime has no NHPP branch and would sample a zero
		// duration; the engine must refuse rather than produce that sample path.
		RuntimeException e = assertThrows(RuntimeException.class,
				new org.junit.jupiter.api.function.Executable() {
					@Override
					public void execute() {
						solve(model);
					}
				});
		assertTrue(e.getMessage().contains("heterogeneous server types"),
				"expected the heterogeneous-server rejection, got: " + e.getMessage());
	}

	@Test
	@DisplayName("Non-cyclic NHPP in steady state is rejected, not averaged over an arbitrary horizon")
	public void testNonCyclicScheduleSteadyStateRejected() {
		Network model = new Network("noncyclic-nhpp");
		Source source = new Source(model, "Source");
		Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
		Sink sink = new Sink(model, "Sink");
		OpenClass oclass = new OpenClass(model, "OpenClass", 0);
		// cyclic = false: the intensity is zero past breakpoints[end], so the
		// steady state is the empty system and a sample-count-driven run would
		// report whatever horizon it happened to stop at.
		source.setArrival(oclass, new NHPP(new double[]{0, 3, 4, 6},
				new double[]{2, 8, 4}, false));
		queue.setService(oclass, new Exp(20.0));
		model.link(model.serialRouting(source, queue, sink));

		RuntimeException e = assertThrows(RuntimeException.class,
				new org.junit.jupiter.api.function.Executable() {
					@Override
					public void execute() {
						solve(model);
					}
				});
		assertTrue(e.getMessage().contains("non-cyclic NHPP"),
				"expected the non-cyclic steady-state rejection, got: " + e.getMessage());
	}
}
