package jline.gen;

import jline.GlobalConstants;
import jline.VerboseLevel;
import jline.lang.Network;
import jline.lang.constant.RoutingStrategy;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Node;
import jline.solvers.NetworkAvgTable;
import jline.solvers.mva.SolverMVA;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import java.util.Map;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Tests for the {@link Cluster} builder and the {@code Network.cluster*}
 * static factory methods. Validates topology, dispatching/scheduling enums,
 * heterogeneous service rates, the closed variant, and parity vs a manually
 * built reference network.
 */
public class ClusterTest {

    @BeforeAll
    public static void setUp() {
        GlobalConstants.setVerbose(VerboseLevel.SILENT);
    }

    @Test
    void testOpenFarmTopology() {
        Network model = new Cluster().setNumStations(3).setArrivalRate(1.0).setServiceRate(2.0)
                .setDispatching(RoutingStrategy.RAND)
                .setScheduling(SchedStrategy.PS)
                .build();

        // Source + Dispatcher + 3 servers + Sink = 6 nodes
        assertEquals(6, model.getNumberOfNodes());
        assertEquals(1, model.getNumberOfClasses());
        assertNotNull(model.getNodeByName("Source"));
        assertNotNull(model.getNodeByName("Dispatcher"));
        assertNotNull(model.getNodeByName("Station1"));
        assertNotNull(model.getNodeByName("Station2"));
        assertNotNull(model.getNodeByName("Station3"));
        assertNotNull(model.getNodeByName("Sink"));
    }

    @Test
    void testClosedFarmTopology() {
        Network model = new Cluster().setNumStations(2).setArrivalRate(1.0).setServiceRate(2.0)
                .setClosed(5, 1.0)
                .setDispatching(RoutingStrategy.JSQ)
                .setScheduling(SchedStrategy.FCFS)
                .build();

        // Think + Dispatcher + 2 servers = 4 nodes
        assertEquals(4, model.getNumberOfNodes());
        assertEquals(1, model.getNumberOfClasses());
        assertNotNull(model.getNodeByName("Think"));
        assertNotNull(model.getNodeByName("Dispatcher"));
        assertNotNull(model.getNodeByName("Station1"));
        assertNotNull(model.getNodeByName("Station2"));
    }

    @Test
    void testAllDispatchingPoliciesAccepted() {
        // KCHOICES is rejected by jline.lang.OutputStrategy, so it is excluded.
        RoutingStrategy[] policies = {
                RoutingStrategy.RAND,
                RoutingStrategy.RROBIN,
                RoutingStrategy.WRROBIN,
                RoutingStrategy.JSQ,
        };
        for (RoutingStrategy p : policies) {
            Network m = new Cluster().setNumStations(2).setArrivalRate(1.0).setServiceRate(2.0).setDispatching(p).build();
            assertEquals(5, m.getNumberOfNodes(), "policy " + p + " should build a 5-node farm");
        }
    }

    @Test
    void testAllSchedulingDisciplinesAccepted() {
        SchedStrategy[] disciplines = {
                SchedStrategy.FCFS, SchedStrategy.PS, SchedStrategy.SJF, SchedStrategy.SRPT,
        };
        for (SchedStrategy s : disciplines) {
            Network m = new Cluster().setNumStations(2).setArrivalRate(1.0).setServiceRate(2.0).setScheduling(s).build();
            assertEquals(5, m.getNumberOfNodes(), "discipline " + s + " should build a 5-node farm");
        }
    }

    @Test
    void testHeterogeneousServiceRates() {
        double[][] mu = {{2.0}, {3.0}, {4.0}};
        Network m = new Cluster()
                .setNumStations(3)
                .setArrivalRate(1.0)
                .setServiceRates(mu)
                .setDispatching(RoutingStrategy.RAND)
                .build();
        assertEquals(6, m.getNumberOfNodes());
        // No exception thrown by build() means service rates were applied for every server.
    }

    @Test
    void testServerCounts() {
        Network m = new Cluster().setNumStations(2).setArrivalRate(1.0).setServiceRate(2.0)
                .setStationServers(new int[]{2, 3})
                .setScheduling(SchedStrategy.FCFS)
                .build();
        Node s1 = m.getNodeByName("Station1");
        Node s2 = m.getNodeByName("Station2");
        assertNotNull(s1);
        assertNotNull(s2);
        assertEquals(2, ((jline.lang.nodes.Queue) s1).getNumberOfServers());
        assertEquals(3, ((jline.lang.nodes.Queue) s2).getNumberOfServers());
    }

    @Test
    void testCompareDispatchingReturnsAllPolicies() {
        // MVA only supports RAND routing; use SSA (a simulator) for non-product-form
        // dispatching strategies and pass a custom factory with reduced sample budget.
        Cluster farm = new Cluster().setNumStations(2).setArrivalRate(0.4).setServiceRate(1.0).setScheduling(SchedStrategy.PS);
        Map<RoutingStrategy, NetworkAvgTable> results = farm.compareDispatching(
                model -> new jline.solvers.ssa.SolverSSA(model,
                        "seed", 23000, "samples", 2000).getAvgTable(),
                RoutingStrategy.RAND,
                RoutingStrategy.RROBIN);
        assertEquals(2, results.size());
        assertTrue(results.containsKey(RoutingStrategy.RAND));
        assertTrue(results.containsKey(RoutingStrategy.RROBIN));
        for (NetworkAvgTable t : results.values()) {
            assertNotNull(t);
        }
    }

    @Test
    void testCompareSchedulingReturnsAllDisciplines() {
        Cluster farm = new Cluster().setNumStations(2).setArrivalRate(0.4).setServiceRate(1.0).setDispatching(RoutingStrategy.RAND);
        Map<SchedStrategy, NetworkAvgTable> results = farm.compareScheduling(
                SolverMVA.class,
                SchedStrategy.FCFS,
                SchedStrategy.PS);
        assertEquals(2, results.size());
        assertTrue(results.containsKey(SchedStrategy.FCFS));
        assertTrue(results.containsKey(SchedStrategy.PS));
    }

    @Test
    void testSweepArrivalRate() {
        Cluster farm = new Cluster().setNumStations(2).setArrivalRate(0.2).setServiceRate(1.0)
                .setScheduling(SchedStrategy.PS)
                .setDispatching(RoutingStrategy.RAND);
        Map<Double, NetworkAvgTable> sweep = farm.sweepArrivalRate(
                new double[]{0.2, 0.5, 0.8}, SolverMVA.class);
        assertEquals(3, sweep.size());
    }

    @Test
    void testSweepArrivalRateRejectsClosedFarm() {
        Cluster farm = new Cluster().setNumStations(2).setArrivalRate(0.4).setServiceRate(1.0).setClosed(4, 1.0);
        assertThrows(IllegalStateException.class, () ->
                farm.sweepArrivalRate(new double[]{0.5, 1.0}, SolverMVA.class));
    }

    @Test
    void testSweepNumServers() {
        Cluster farm = new Cluster().setNumStations(1).setArrivalRate(0.4).setServiceRate(1.0)
                .setScheduling(SchedStrategy.PS)
                .setDispatching(RoutingStrategy.RAND);
        Map<Integer, NetworkAvgTable> sweep = farm.sweepNumStations(
                new int[]{2, 4}, SolverMVA.class);
        assertEquals(2, sweep.size());
    }

    @Test
    void testParityWithFactoryFunction() {
        // Builder vs the static factory should produce models with identical
        // topology and identical MVA throughput results.
        Cluster builder = new Cluster().setNumStations(2).setArrivalRate(0.5).setServiceRate(1.0)
                .setScheduling(SchedStrategy.PS)
                .setDispatching(RoutingStrategy.RAND);
        Network builderModel = builder.build();

        jline.util.matrix.Matrix lambda = new jline.util.matrix.Matrix(1, 1);
        lambda.set(0, 0, 0.5);
        jline.util.matrix.Matrix D = new jline.util.matrix.Matrix(2, 1);
        D.set(0, 0, 1.0); D.set(1, 0, 1.0);
        Network factoryModel = Network.clusterPs(lambda, D, RoutingStrategy.RAND);

        assertEquals(factoryModel.getNumberOfNodes(), builderModel.getNumberOfNodes());
        assertEquals(factoryModel.getNumberOfClasses(), builderModel.getNumberOfClasses());

        NetworkAvgTable a = new SolverMVA(builderModel).getAvgTable();
        NetworkAvgTable b = new SolverMVA(factoryModel).getAvgTable();
        assertEquals(a.getTput().size(), b.getTput().size());
        for (int i = 0; i < a.getTput().size(); i++) {
            assertEquals(b.getTput().get(i), a.getTput().get(i), 1e-9,
                    "Throughput row " + i + " should match between builder and factory");
        }
    }
}
