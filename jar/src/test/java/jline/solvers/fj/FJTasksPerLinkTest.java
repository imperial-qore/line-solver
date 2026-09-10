package jline.solvers.fj;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import jline.api.fj.FJ_ordstat_exp;
import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.RoutingMatrix;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Fork;
import jline.lang.nodes.Join;
import jline.lang.nodes.Queue;
import jline.lang.processes.Exp;
import jline.solvers.mva.SolverMVA;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;

/**
 * Fork.setTasksPerLink(w) under the MMT fork-join transform.
 *
 * <p>A fork with tasksPerLink = w sends w IDENTICAL tasks down each of its B links, so a firing
 * creates w*B siblings. The join synchronises on the order statistic of that many branch times,
 * each branch REPLICATED w times, and NOT on w times the order statistic of B: the latter is
 * w*H_B/mu where the answer is H_(w*B)/mu, already 3.0/mu against 2.083/mu at B = w = 2. A branch
 * station also holds w jobs per circulating parent, so the chain population is not its bound.</p>
 *
 * <p>The reference values are SolverMVA on the same model in MATLAB, native python and the C++
 * port, which agree with the JAR to every printed digit; SolverJMT and SolverLDES (which simulate
 * tasksPerLink directly and agree with each other to 0.3%) put the true throughput about ten
 * percent higher, the same bias an ordinary fork-join carries under this transform.</p>
 */
public class FJTasksPerLinkTest {

    private static final double REL = 1e-5;

    private static Network build(int N, int w, int nb) {
        Network model = new Network("FJ-TPL");
        Delay delay = new Delay(model, "Delay");
        Queue[] qs = new Queue[nb];
        for (int b = 0; b < nb; b++) {
            qs[b] = new Queue(model, "Queue" + (b + 1), SchedStrategy.PS);
        }
        Fork fork = new Fork(model, "Fork");
        Join join = new Join(model, "Join", fork);
        ClosedClass c = new ClosedClass(model, "class1", N, delay);
        delay.setService(c, new Exp(1.0));
        for (int b = 0; b < nb; b++) {
            qs[b].setService(c, new Exp(2.0));
        }
        fork.setTasksPerLink(w);
        RoutingMatrix P = model.initRoutingMatrix();
        P.set(c, c, delay, fork, 1.0);
        for (int b = 0; b < nb; b++) {
            P.set(c, c, fork, qs[b], 1.0);
            P.set(c, c, qs[b], join, 1.0);
        }
        P.set(c, c, join, delay, 1.0);
        model.link(P);
        return model;
    }

    private static double tput(Network model) {
        SolverMVA solver = new SolverMVA(model);
        solver.getAvg();
        return solver.result.TN.get(0, 0);
    }

    @Test
    public void joinSynchronisesOnAllSiblings() {
        double[][] ref = {{1.167124, 0.702141, 0.507307}, {1.058612, 0.662789, 0.487276}};
        for (int bi = 0; bi < 2; bi++) {
            int nb = bi + 2;
            for (int wi = 0; wi < 3; wi++) {
                int w = wi + 1;
                assertEquals(ref[bi][wi], tput(build(3, w, nb)), ref[bi][wi] * REL,
                        "branches=" + nb + " tasksPerLink=" + w);
            }
        }
    }

    @Test
    public void oneTaskPerLinkIsTheOrdinaryForkJoin() {
        // w = 1 replicates a branch to itself and scales the capacity by one, so a plain
        // fork must answer exactly as it did before tasksPerLink was handled.
        for (int nb = 2; nb <= 3; nb++) {
            Network model = build(3, 1, nb);
            double expected = (nb == 2) ? 1.167124 : 1.058612;
            assertEquals(expected, tput(model), expected * REL);
            assertEquals(3.0, model.getStruct().classcap.get(0, 0), 1e-9);
        }
    }

    @Test
    public void moreTasksPerLinkLowersThroughput() {
        for (int nb = 2; nb <= 3; nb++) {
            double x1 = tput(build(3, 1, nb));
            double x2 = tput(build(3, 2, nb));
            double x3 = tput(build(3, 3, nb));
            assertTrue(x1 > x2 && x2 > x3, "throughput must fall as tasksPerLink grows");
        }
    }

    @Test
    public void branchCapacityScalesWithTasksPerLink() {
        for (int nb = 2; nb <= 3; nb++) {
            for (int w = 1; w <= 3; w++) {
                Network model = build(3, w, nb);
                assertEquals(3.0 * w, model.getStruct().classcap.get(0, 0), 1e-9,
                        "branches=" + nb + " tasksPerLink=" + w);
            }
        }
    }

    @Test
    public void replicatingIsNotScalingTheOrderStatistic() {
        // w*E[X_(B)] is strictly larger than E[X_(w*B)] on identical branches, so the two
        // cannot both be the synchronisation instant.
        Matrix ri = new Matrix(1, 2);
        ri.set(0, 0, 0.5);
        ri.set(0, 1, 0.5);
        Matrix rep = new Matrix(1, 4);
        for (int i = 0; i < 4; i++) {
            rep.set(0, i, 0.5);
        }
        double scaled = 2.0 * FJ_ordstat_exp.fj_ordstat_exp(ri, 2);
        double replicated = FJ_ordstat_exp.fj_ordstat_exp(rep, 4);
        assertEquals(2.0 * 0.5 * (1 + 1.0 / 2), scaled, 1e-12);
        assertEquals(0.5 * (1 + 1.0 / 2 + 1.0 / 3 + 1.0 / 4), replicated, 1e-12);
        assertTrue(replicated < scaled);
    }
}
