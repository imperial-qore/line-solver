package jline.solvers.mva;

import jline.GlobalConstants;
import jline.VerboseLevel;
import jline.api.npfqn.Npfqn_sqd;
import jline.io.Ret;
import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.RoutingMatrix;
import jline.lang.constant.DropStrategy;
import jline.lang.constant.SchedStrategy;
import jline.lang.constant.SolverType;
import jline.lang.nodes.Queue;
import jline.lang.processes.Exp;
import jline.solvers.NetworkAvgTable;
import jline.solvers.SolverOptions;
import jline.solvers.wrappers.jmt.SolverJMT;
import jline.solvers.ldes.LDESOptions;
import jline.solvers.ldes.SolverLDES;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;

import java.util.List;

/**
 * Oracle harness: compares the SQD (Npfqn_sqd) BAS approximation against JMT and LDES
 * simulation on the finite-buffer closed networks that regress in FiniteBufferBASTest.
 * Not a gate — prints Q1 throughput sweeps for accuracy assessment.
 */
public class MgsOracleTest {

    private static final int    SEED    = 23000;
    private static final int    SAMPLES = 200000;

    static Network buildClosed(String name, double[] mu, int[] K, double[][] r, int N) {
        Network model = new Network(name);
        Queue[] q = new Queue[mu.length];
        for (int i = 0; i < mu.length; i++) {
            q[i] = new Queue(model, "Q" + (i + 1), SchedStrategy.FCFS);
        }
        ClosedClass jobs = new ClosedClass(model, "Jobs", N, q[0]);
        for (int i = 0; i < mu.length; i++) {
            q[i].setService(jobs, new Exp(mu[i]));
            q[i].setCapacity(K[i]);
            q[i].setDropRule(jobs, DropStrategy.BlockingAfterService);
            q[i].setNumberOfServers(1);
        }
        RoutingMatrix P = new RoutingMatrix(model, java.util.Arrays.asList(jobs),
                java.util.Arrays.asList(q));
        for (int i = 0; i < mu.length; i++)
            for (int j = 0; j < mu.length; j++)
                if (r[i][j] > 0) P.set(jobs, jobs, q[i], q[j], r[i][j]);
        model.link(P);
        return model;
    }

    static double[][] cyc(int m) {
        double[][] r = new double[m][m];
        for (int i = 0; i < m; i++) r[i][(i + 1) % m] = 1.0;
        return r;
    }

    private static double sqdTput(Network model, int N) {
        Ret.pfqnMVA res = Npfqn_sqd.npfqn_sqd(model.getStruct(false), N,
                0, false, Npfqn_sqd.NeighborMode.DOWNSTREAM, Npfqn_sqd.V1Policy.COMPOUND, null);
        return res.X.get(0, 0);
    }

    /** Feasible BAS initial marginal: greedily fill stations up to capacity K. */
    private static void initFeasible(Network model, int[] K, int N) {
        int m = K.length;
        Matrix n = new Matrix(m, 1);
        int left = N;
        for (int i = 0; i < m && left > 0; i++) {
            int put = Math.min(K[i], left);
            n.set(i, 0, put);
            left -= put;
        }
        model.initFromMarginal(n);
    }

    private static double jmtTput(Network model, int[] K, int N) {
        initFeasible(model, K, N);
        SolverOptions o = new SolverOptions(SolverType.JMT);
        o.verbose = VerboseLevel.SILENT;
        o.samples = SAMPLES;
        o.seed = SEED;
        NetworkAvgTable t = new SolverJMT(model, o).getAvgTable();
        return t.getTput().get(0);
    }

    private static double ldesTput(Network model, int[] K, int N) {
        initFeasible(model, K, N);
        LDESOptions o = new LDESOptions();
        o.verbose = VerboseLevel.SILENT;
        o.samples = SAMPLES;
        o.seed = SEED;
        NetworkAvgTable t = new SolverLDES(model, o).getAvgTable();
        return t.getTput().get(0);
    }

    private void sweep(String tag, double[] mu, int[] K, double[][] r, int[] Ns, double[] ref) {
        for (int idx = 0; idx < Ns.length; idx++) {
            int N = Ns[idx];
            String nm = tag + "_N" + N;
            double sqd = sqdTput(buildClosed(nm + "a", mu, K, r, N), N);
            try { jmtTput(buildClosed(nm + "b", mu, K, r, N), K, N); } catch (Throwable e) { }
            try { ldesTput(buildClosed(nm + "c", mu, K, r, N), K, N); } catch (Throwable e) { }
        }
    }

    @Test
    public void debugCap() {
        Network m = buildClosed("dbg", new double[]{1,1}, new int[]{4,4}, cyc(2), 5);
        m.getStruct(false);
        int[][] trials = {{4,1},{3,2},{2,3},{3,1},{4,0}};
        StringBuilder sb = new StringBuilder();
        for (int[] t : trials) {
            Network mm = buildClosed("dbg"+t[0]+t[1], new double[]{1,1}, new int[]{4,4}, cyc(2), t[0]+t[1]);
            Matrix n = new Matrix(2,1); n.set(0,0,t[0]); n.set(1,0,t[1]);
            String res;
            try { mm.initFromMarginal(n); res = "OK"; }
            catch (Throwable e) { res = "THREW " + e.getMessage(); }
            sb.append("[").append(t[0]).append(",").append(t[1]).append("] -> ").append(res).append("\n");
        }
        try { java.nio.file.Files.write(java.nio.file.Paths.get("/tmp/mgs_dbg.txt"), sb.toString().getBytes()); }
        catch (Exception ex) { }
    }

    @Test
    public void oracle() {
        GlobalConstants.setVerbose(VerboseLevel.SILENT);
        int[] n17 = {1,2,3,4,5,6,7};
        sweep("Table1", new double[]{1,1}, new int[]{4,4}, cyc(2), n17,
                new double[]{0.500,0.667,0.750,0.800,0.833,0.800,0.750});
        int[] n19 = {1,2,3,4,5,6,7,8,9};
        sweep("Table2", new double[]{2,4}, new int[]{4,6}, cyc(2), n19,
                new double[]{1.333,1.714,1.867,1.936,1.968,1.968,1.968,1.936,1.867});
        int[] n = {10,20,30};
        sweep("Table8", new double[]{2,2,2,2,2,2,2,2}, new int[]{4,4,4,4,4,4,4,4}, cyc(8), n,
                new double[]{1.175,1.380,1.037});
        sweep("Table9", new double[]{2,8,5,2.5,2,4,1.25,5}, new int[]{5,2,3,5,4,3,7,3}, cyc(8), n,
                new double[]{1.194,1.232,1.072});
    }
}
