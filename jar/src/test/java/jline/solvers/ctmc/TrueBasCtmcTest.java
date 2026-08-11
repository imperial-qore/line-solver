package jline.solvers.ctmc;

import jline.VerboseLevel;
import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.RoutingMatrix;
import jline.lang.constant.DropStrategy;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Queue;
import jline.lang.processes.Exp;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

/**
 * Exact validation of true Blocking-After-Service (BAS) semantics in SolverCTMC.
 *
 * Under true BAS a job that completes service at a finite-capacity station whose
 * downstream is full holds its server (it is not re-served) and moves the instant the
 * destination frees. This differs from repetitive service (RS), where the held job is
 * re-served and, service being exponential and hence memoryless, only leaves at rate
 * mu once the destination frees. The two therefore have distinct throughputs, which is
 * what pins the semantics here.
 *
 * Oracle: closed 2-station cyclic network, mu=[1,1], K=[4,4]. The exact BAS throughputs
 * are 0.8333 (N=5), 0.8000 (N=6) and 0.7500 (N=7); RS gives 0.750 and 0.667 for N=5,6.
 */
public class TrueBasCtmcTest {

    private static final double TOL = 1e-3;

    private static Network buildCyclic(String name, double[] mu, int[] K, int N, DropStrategy rule) {
        Network model = new Network(name);
        Queue[] q = new Queue[mu.length];
        for (int i = 0; i < mu.length; i++) {
            q[i] = new Queue(model, "Q" + (i + 1), SchedStrategy.FCFS);
        }
        ClosedClass jobs = new ClosedClass(model, "Jobs", N, q[0]);
        for (int i = 0; i < mu.length; i++) {
            q[i].setService(jobs, new Exp(mu[i]));
            q[i].setCapacity(K[i]);
            q[i].setDropRule(jobs, rule);
        }
        RoutingMatrix P = model.initRoutingMatrix();
        for (int i = 0; i < mu.length; i++) {
            P.set(jobs, jobs, q[i], q[(i + 1) % mu.length], 1.0);
        }
        model.link(P);
        return model;
    }

    private static double tput(Network model) {
        return new SolverCTMC(model, "verbose", VerboseLevel.SILENT).getAvgTable().getTput().get(0);
    }

    @Test
    public void trueBasMatchesExactThroughput() {
        double[] mu = {1.0, 1.0};
        int[] K = {4, 4};
        int[] Ns = {5, 6, 7};
        double[] exact = {0.8333, 0.8000, 0.7500};
        for (int i = 0; i < Ns.length; i++) {
            assertEquals(exact[i],
                    tput(buildCyclic("bas" + Ns[i], mu, K, Ns[i], DropStrategy.BlockingAfterService)),
                    TOL, "true BAS throughput, N=" + Ns[i]);
        }
    }

    @Test
    public void reServiceOnRejectionKeepsRepetitiveServiceSemantics() {
        double[] mu = {1.0, 1.0};
        int[] K = {4, 4};
        int[] Ns = {5, 6};
        double[] rsExact = {0.750, 0.667};
        for (int i = 0; i < Ns.length; i++) {
            assertEquals(rsExact[i],
                    tput(buildCyclic("rs" + Ns[i], mu, K, Ns[i], DropStrategy.ReServiceOnRejection)),
                    TOL, "RS throughput must not follow the BAS branch, N=" + Ns[i]);
        }
    }

    @Test
    public void nonBlockingClosedNetworkIsUnaffected() {
        // Symmetric 2-station cyclic network with ample capacity: the chain is uniform
        // over its N+1 states, so X = mu * (1 - 1/(N+1)) = 0.75 for N=3. Guards against
        // the blocked marker leaking into models that do not use BAS.
        assertEquals(0.75,
                tput(buildCyclic("pf", new double[]{1.0, 1.0}, new int[]{100, 100}, 3,
                        DropStrategy.WaitingQueue)),
                TOL, "non-BAS closed network must be unchanged");
    }
}
