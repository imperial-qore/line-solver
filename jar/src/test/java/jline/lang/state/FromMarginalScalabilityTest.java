package jline.lang.state;

import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.NetworkStruct;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.processes.Exp;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.Timeout;

import static org.junit.jupiter.api.Assertions.assertEquals;

/**
 * Regression test for the FCFS initial-state construction.
 *
 * Building the default initial state used to enumerate every distinct permutation
 * of the FCFS buffer, which is factorial in the population: the model below has
 * 14!/(2!*3!^4) = 33,634,368 orderings and could not be built at all. Only one
 * ordering is ever kept, so it is now produced directly.
 */
public class FromMarginalScalabilityTest {

    /** Builds a Delay + FCFS Queue closed model with R classes of N jobs, referenced at the queue. */
    private static Network buildModel(int R, int N) {
        Network model = new Network("fcfs-scalability");
        Delay delay = new Delay(model, "Delay");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        queue.setNumberOfServers(1);
        for (int r = 0; r < R; r++) {
            ClosedClass jobClass = new ClosedClass(model, "Class" + (r + 1), N, queue, 0);
            delay.setService(jobClass, new Exp(1.0));
            queue.setService(jobClass, new Exp(0.5));
        }
        model.link(Network.serialRouting(delay, queue));
        return model;
    }

    @Test
    @Timeout(60)
    public void fcfsInitialStateIsBuiltWithoutEnumeratingPermutations() {
        int R = 5, N = 3;
        NetworkStruct sn = buildModel(R, N).getStruct(true);

        // Stateful index 2 is the Queue (index 1 is the Delay); its initial state is
        // [buffer classes, per-class phase counts of the jobs in service].
        Matrix state = sn.state.get(sn.stations.get(1));
        assertEquals(1, state.getNumRows());
        // 14 buffered jobs (one class-1 job is in the single server) + R phase columns
        assertEquals(R * N - 1 + R, state.getNumCols());

        // Buffer classes come out in descending order, the single lexicographic
        // representative of the R*N-1 job orderings.
        double[] expected = {5, 5, 5, 4, 4, 4, 3, 3, 3, 2, 2, 2, 1, 1, 1, 0, 0, 0, 0};
        for (int j = 0; j < expected.length; j++) {
            assertEquals(expected[j], state.get(0, j), 0.0, "state column " + j);
        }
    }
}
