/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.solvers.ag;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import jline.solvers.ag.handlers.RCATModel;
import jline.util.matrix.Matrix;

import static jline.io.InputOutput.line_error;
import static jline.io.InputOutput.line_warning;

/**
 * Execution backend of the reversed-rate fixed point.
 *
 * <p>WHAT MAKES THIS SAFE IS THE DECOMPOSITION, NOT THE SCHEDULING. Agent k's
 * generator is
 *
 * <pre>Q_k(x) = L_k + sum_{c passive at k} x_c Pb_c</pre>
 *
 * so an agent reads the rest of the model only through the scalar reversed rates
 * x, and it writes only its own slot of the sweep's output. The sweep is Jacobi
 * -- every x_a is read off the PREVIOUS sweep's stationary vectors and only then
 * do the agents re-solve -- so the agent order is immaterial. A parallel or
 * distributed sweep therefore produces the SAME iterates as the serial one.
 *
 * <p>How far that survives floating point depends on who runs the agent solve.
 * {@code threads} is BIT-IDENTICAL to {@code serial}: same code, same process,
 * only the evaluation order differs and the order is what does not matter.
 * {@code cluster} is bit-identical only when the worker runs the same
 * implementation as the coordinator -- the wire is exact, since JSON round-trips
 * a double without loss, but the stationary vector comes back from the WORKER's
 * solve, so a Python coordinator driving a Java ag-worker agrees to a few ulp
 * (measured: 1.1e-16 on the M/M/1 tandem) rather than bit for bit. That is the
 * ordinary cross-codebase difference, not a protocol defect, but it means a
 * cluster run must not regenerate a seeded golden a same-language run will read.
 *
 * <p>The three backends differ only in who evaluates an agent:
 * <ul>
 *   <li>{@code serial}: the caller's thread, in agent order (the reference).</li>
 *   <li>{@code threads}: a fixed pool, one task per agent, barrier per sweep.</li>
 *   <li>{@code cluster}: agents partitioned over ag-worker processes; the static
 *       half of each agent is shipped once, and only x crosses the wire per
 *       sweep.</li>
 * </ul>
 *
 * <p>A cluster worker that is missing, slow or broken is NOT fatal: its agents
 * are solved on the coordinator through the same agent path, so the answer is
 * the run's answer either way and only the wall clock changes.</p>
 */
public final class AgExec implements AutoCloseable {

    private final String mode;
    private final ExecutorService pool;
    private final List<AgWorkerClient> workers;
    /** Agent indices owned by each worker; the tail beyond the workers is local. */
    private final List<List<Integer>> owns;
    private boolean assigned;

    private AgExec(String mode, ExecutorService pool, List<AgWorkerClient> workers) {
        this.mode = mode;
        this.pool = pool;
        this.workers = workers;
        this.owns = new ArrayList<List<Integer>>();
    }

    /**
     * Resolve the backend named by the options. Returns null for {@code serial},
     * which is the caller's own loop and needs no object; that keeps the default
     * path free of any scheduling machinery at all.
     */
    public static AgExec create(AGOptions options) {
        if (options == null) return null;
        String mode = (options.exec == null) ? AGOptions.EXEC_SERIAL : options.exec.toLowerCase();
        if (AGOptions.EXEC_SERIAL.equals(mode)) {
            return null;
        }
        if (isParallel(mode)) {
            int n = options.nworkers > 0
                    ? options.nworkers
                    : Runtime.getRuntime().availableProcessors();
            return new AgExec(mode, Executors.newFixedThreadPool(n), null);
        }
        if (AGOptions.EXEC_CLUSTER.equals(mode)) {
            if (options.endpoints == null || options.endpoints.isEmpty()) {
                line_error("AgExec", "The 'cluster' execution backend needs worker endpoints: set "
                        + "AGOptions.endpoints to a list of \"host:port\" strings, each one an "
                        + "ag-worker started with 'java -cp jline.jar jline.cli.AgWorker -p <port>'.");
            }
            List<AgWorkerClient> ws = new ArrayList<AgWorkerClient>();
            for (int i = 0; i < options.endpoints.size(); i++) {
                ws.add(new AgWorkerClient(options.endpoints.get(i), options.workerTimeout));
            }
            return new AgExec(mode, null, ws);
        }
        // "threads" was this backend's name until 2026-08-19. Name the rename
        // rather than reporting a backend that still exists as unknown.
        String hint = "threads".equals(mode)
                ? " ('threads' was renamed to 'parallel', alias 'para')" : "";
        line_error("AgExec", "Unknown AG execution backend '" + mode
                + "'. Use 'serial', 'parallel' (alias 'para') or 'cluster'." + hint);
        return null;
    }

    /** True for either accepted spelling of the local-thread-pool backend. */
    public static boolean isParallel(String mode) {
        return AGOptions.EXEC_PARALLEL.equals(mode) || AGOptions.EXEC_PARA.equals(mode);
    }

    public String mode() {
        return mode;
    }

    /**
     * Evaluate every agent at the reversed rates x, writing agent k's generator
     * into {@code Qs[k]} and its stationary vector into {@code pis[k]}.
     */
    public void sweep(Matrix x, Matrix[] Aa, Matrix[] Pb, Matrix[] L,
                      int[] ACT, int[] PSV, int numProcesses, int numActions,
                      int[] N, RCATModel rcat, Matrix[] Qs, Matrix[] pis) {
        if (isParallel(mode)) {
            sweepThreads(x, Aa, Pb, L, ACT, PSV, numProcesses, numActions, N, rcat, Qs, pis);
        } else {
            sweepCluster(x, Aa, Pb, L, ACT, PSV, numProcesses, numActions, N, rcat, Qs, pis);
        }
    }

    private void sweepThreads(final Matrix x, final Matrix[] Aa, final Matrix[] Pb, final Matrix[] L,
                              final int[] ACT, final int[] PSV, int numProcesses, final int numActions,
                              final int[] N, final RCATModel rcat, final Matrix[] Qs, final Matrix[] pis) {
        List<Future<?>> tasks = new ArrayList<Future<?>>(numProcesses);
        for (int kk = 0; kk < numProcesses; kk++) {
            final int k = kk;
            tasks.add(pool.submit(new Callable<Void>() {
                public Void call() {
                    Qs[k] = AgAgent.generator(k, x, Aa, Pb, L, ACT, PSV, numActions, N);
                    pis[k] = AgAgent.stationary(Qs[k], rcat, k);
                    return null;
                }
            }));
        }
        for (int i = 0; i < tasks.size(); i++) {
            try {
                tasks.get(i).get();
            } catch (Exception e) {
                // An agent solve that throws is a defect in the agent, not in the
                // scheduling, so it must surface as itself rather than as a
                // half-filled sweep that fails later somewhere unrelated.
                throw new RuntimeException("AG agent " + i + " failed on the thread pool: "
                        + e.getMessage(), e);
            }
        }
    }

    private void sweepCluster(Matrix x, Matrix[] Aa, Matrix[] Pb, Matrix[] L,
                              int[] ACT, int[] PSV, int numProcesses, int numActions,
                              int[] N, RCATModel rcat, Matrix[] Qs, Matrix[] pis) {
        if (!assigned) {
            partition(numProcesses);
            for (int w = 0; w < workers.size(); w++) {
                AgWorkerClient c = workers.get(w);
                if (!c.isLive()) continue;
                try {
                    c.assign(owns.get(w), Aa, Pb, L, ACT, PSV, numActions, N, rcat);
                } catch (Exception e) {
                    line_warning("AgExec", "AG worker " + c.endpoint() + " refused the assignment ("
                            + e.getMessage() + "); its " + owns.get(w).size()
                            + " agent(s) run on the coordinator instead.");
                    c.kill();
                }
            }
            assigned = true;
        }

        // The generator is rebuilt here in any case: the metrics stage reads it,
        // and rebuilding is cheaper than transporting it.
        for (int k = 0; k < numProcesses; k++) {
            Qs[k] = AgAgent.generator(k, x, Aa, Pb, L, ACT, PSV, numActions, N);
        }

        boolean[] pending = new boolean[numProcesses];
        for (int k = 0; k < numProcesses; k++) pending[k] = true;

        for (int w = 0; w < workers.size(); w++) {
            AgWorkerClient c = workers.get(w);
            List<Integer> mine = owns.get(w);
            if (!c.isLive() || mine.isEmpty()) continue;
            try {
                List<Matrix> got = c.sweep(x, mine);
                for (int t = 0; t < mine.size(); t++) {
                    int k = mine.get(t).intValue();
                    pis[k] = got.get(t);
                    pending[k] = false;
                }
            } catch (Exception e) {
                line_warning("AgExec", "AG worker " + c.endpoint() + " failed mid-sweep ("
                        + e.getMessage() + "); its " + mine.size()
                        + " agent(s) are solved on the coordinator for the rest of the run.");
                c.kill();
            }
        }

        for (int k = 0; k < numProcesses; k++) {
            if (pending[k]) {
                pis[k] = AgAgent.stationary(Qs[k], rcat, k);
            }
        }
    }

    /**
     * Round-robin over the worker list in agent index order, computed before any
     * connection is attempted so that a dead worker does not shift the others'
     * agents. The partition is a pure function of (agent count, worker count), so
     * a rerun assigns the same agents to the same workers.
     */
    private void partition(int numProcesses) {
        owns.clear();
        for (int w = 0; w < workers.size(); w++) owns.add(new ArrayList<Integer>());
        for (int k = 0; k < numProcesses; k++) {
            owns.get(k % workers.size()).add(Integer.valueOf(k));
        }
    }

    public void close() {
        if (pool != null) pool.shutdownNow();
        if (workers != null) {
            for (int w = 0; w < workers.size(); w++) workers.get(w).close();
        }
    }
}
