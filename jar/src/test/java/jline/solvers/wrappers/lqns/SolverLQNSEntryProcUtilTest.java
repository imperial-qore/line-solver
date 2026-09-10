package jline.solvers.wrappers.lqns;

import jline.examples.java.basic.LayeredModel;
import jline.lang.layered.LayeredNetwork;
import jline.solvers.LayeredNetworkAvgTable;
import org.junit.jupiter.api.Test;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assumptions.assumeTrue;

/**
 * An entry's Util is the sum over ITS OWN activities, not the raw .lqxo attribute.
 *
 * <p>lqns credits host work to whichever level declares the host demand. In the
 * activity-graph ({@code task-activities}) form -- the only form LINE's writer
 * ever emits -- an entry declares none, so lqns writes a literal
 * {@code proc-utilization="0"} on every {@code result-entry} and puts the work on
 * the {@code result-activity} rows. Read verbatim, every entry's Util came back 0
 * while the model plainly runs.</p>
 *
 * <p>{@code lqn_twotasks} is the case that pins the rule down: T2 hosts TWO
 * entries, so the entry sums (three unit-demand activities against one) are NOT
 * the task's proc-utilization. A fix that copied the task value, or that summed
 * every activity of the task into each of its entries, would pass a one-entry
 * task and fail here.</p>
 */
public class SolverLQNSEntryProcUtilTest {

    /** Util by node name, from one lqns run of the two-task model. */
    private static Map<String, Double> utilByName() throws Exception {
        LayeredNetwork model = LayeredModel.lqn_twotasks();
        SolverLQNS solver = new SolverLQNS(model);
        LayeredNetworkAvgTable table = (LayeredNetworkAvgTable) solver.getAvgTable();
        List<String> names = table.getNodeNames();
        List<Double> util = table.getUtil();
        Map<String, Double> out = new HashMap<String, Double>();
        for (int i = 0; i < names.size(); i++) {
            out.put(names.get(i), util.get(i));
        }
        return out;
    }

    @Test
    public void entryProcUtilSumsItsOwnActivities() throws Exception {
        assumeTrue(SolverLQNS.isAvailable(), "no lqns binary on the PATH");
        Map<String, Double> u = utilByName();

        double chain = u.get("A20") + u.get("A21") + u.get("A22");
        assertEquals(chain, u.get("E2"), 1e-6 * Math.max(1.0, chain));
        assertEquals(u.get("A3"), u.get("E3"), 1e-6);
        assertEquals(u.get("A1"), u.get("E1"), 1e-6);

        // The regression itself: read verbatim, every entry comes back 0.
        assertTrue(u.get("E1") > 0.0, "E1 utilization must not be the raw 0");
        assertTrue(u.get("E2") > 0.0, "E2 utilization must not be the raw 0");
        assertTrue(u.get("E3") > 0.0, "E3 utilization must not be the raw 0");
    }

    @Test
    public void entryProcUtilIsNotTheTaskValue() throws Exception {
        assumeTrue(SolverLQNS.isAvailable(), "no lqns binary on the PATH");
        Map<String, Double> u = utilByName();

        // Three unit-demand activities against one, at the same throughput.
        assertEquals(3.0 * u.get("E3"), u.get("E2"), 1e-3 * u.get("E2"));
        // Together they are the task, so neither entry can be carrying its value.
        assertEquals(u.get("T2"), u.get("E2") + u.get("E3"), 1e-6 * u.get("T2"));
        assertTrue(u.get("E2") > 1.1 * u.get("E3"));
    }
}
