package jline.lang;

import jline.lang.constant.HeteroSchedPolicy;
import jline.lang.constant.ServerType;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.processes.Distribution;
import jline.lang.processes.Exp;
import org.junit.jupiter.api.Test;

import java.util.List;
import java.util.Map;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Unit tests for heterogeneous server support.
 * Tests ServerType, HeteroSchedPolicy, and Queue heterogeneous server methods.
 */
public class HeteroServerTest {

    // ========== ServerType Tests ==========

    @Test
    void testServerTypeCreation() {
        ServerType st = new ServerType("Fast", 2);
        assertNotNull(st);
        assertEquals("Fast", st.getName());
        assertEquals(2, st.getNumOfServers());
        assertTrue(st.getCompatibleClasses().isEmpty());
    }

    @Test
    void testServerTypeSetNumOfServers() {
        ServerType st = new ServerType("Test", 1);
        assertEquals(1, st.getNumOfServers());
        st.setNumOfServers(5);
        assertEquals(5, st.getNumOfServers());
    }

    @Test
    void testServerTypeCompatibility() {
        Network model = new Network("TestModel");
        Source source = new Source(model, "Source");
        Sink sink = new Sink(model, "Sink");

        OpenClass class1 = new OpenClass(model, "Class1");
        OpenClass class2 = new OpenClass(model, "Class2");

        ServerType st = new ServerType("Fast", 2);

        // Initially no classes are compatible
        assertFalse(st.isCompatible(class1));
        assertFalse(st.isCompatible(class2));

        // Add class1 as compatible
        st.addCompatibleClass(class1);
        assertTrue(st.isCompatible(class1));
        assertFalse(st.isCompatible(class2));
        assertEquals(1, st.getCompatibleClasses().size());

        // Add class2 as compatible
        st.addCompatibleClass(class2);
        assertTrue(st.isCompatible(class1));
        assertTrue(st.isCompatible(class2));
        assertEquals(2, st.getCompatibleClasses().size());
    }

    @Test
    void testServerTypeParentQueue() {
        Network model = new Network("TestModel");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);

        ServerType st = new ServerType("Fast", 2);
        assertNull(st.getParentQueue());

        st.setParentQueue(queue);
        assertEquals(queue, st.getParentQueue());
    }

    // ========== HeteroSchedPolicy Tests ==========

    @Test
    void testHeteroSchedPolicyValues() {
        assertNotNull(HeteroSchedPolicy.ORDER);
        assertNotNull(HeteroSchedPolicy.ALIS);
        assertNotNull(HeteroSchedPolicy.ALFS);
        assertNotNull(HeteroSchedPolicy.FAIRNESS);
        assertNotNull(HeteroSchedPolicy.FSF);
        assertNotNull(HeteroSchedPolicy.RAIS);
    }

    @Test
    void testHeteroSchedPolicyFromText() {
        assertEquals(HeteroSchedPolicy.ORDER, HeteroSchedPolicy.fromText("ORDER"));
        assertEquals(HeteroSchedPolicy.ALIS, HeteroSchedPolicy.fromText("ALIS"));
        assertEquals(HeteroSchedPolicy.ALFS, HeteroSchedPolicy.fromText("ALFS"));
        assertEquals(HeteroSchedPolicy.FAIRNESS, HeteroSchedPolicy.fromText("FAIRNESS"));
        assertEquals(HeteroSchedPolicy.FSF, HeteroSchedPolicy.fromText("FSF"));
        assertEquals(HeteroSchedPolicy.RAIS, HeteroSchedPolicy.fromText("RAIS"));
    }

    @Test
    void testHeteroSchedPolicyFromTextCaseInsensitive() {
        assertEquals(HeteroSchedPolicy.FSF, HeteroSchedPolicy.fromText("fsf"));
        assertEquals(HeteroSchedPolicy.ALIS, HeteroSchedPolicy.fromText("Alis"));
    }

    @Test
    void testHeteroSchedPolicyFromTextInvalid() {
        assertThrows(IllegalArgumentException.class, () -> HeteroSchedPolicy.fromText("INVALID"));
        assertThrows(IllegalArgumentException.class, () -> HeteroSchedPolicy.fromText(""));
        assertThrows(IllegalArgumentException.class, () -> HeteroSchedPolicy.fromText(null));
    }

    @Test
    void testHeteroSchedPolicyToText() {
        assertEquals("ORDER", HeteroSchedPolicy.ORDER.toText());
        assertEquals("ALIS", HeteroSchedPolicy.ALIS.toText());
        assertEquals("ALFS", HeteroSchedPolicy.ALFS.toText());
        assertEquals("FAIRNESS", HeteroSchedPolicy.FAIRNESS.toText());
        assertEquals("FSF", HeteroSchedPolicy.FSF.toText());
        assertEquals("RAIS", HeteroSchedPolicy.RAIS.toText());
    }

    // ========== Queue Heterogeneous Server Tests ==========

    @Test
    void testQueueIsHeterogeneousInitiallyFalse() {
        Network model = new Network("TestModel");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);

        assertFalse(queue.isHeterogeneous());
        assertTrue(queue.getServerTypes().isEmpty());
    }

    @Test
    void testQueueAddServerType() {
        Network model = new Network("TestModel");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);

        ServerType fast = new ServerType("Fast", 2);
        ServerType slow = new ServerType("Slow", 3);

        queue.addServerType(fast);
        assertTrue(queue.isHeterogeneous());
        assertEquals(1, queue.getServerTypes().size());
        assertEquals(queue, fast.getParentQueue());

        queue.addServerType(slow);
        assertEquals(2, queue.getServerTypes().size());
        assertEquals(queue, slow.getParentQueue());
    }

    @Test
    void testQueueAddServerTypeNull() {
        Network model = new Network("TestModel");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);

        assertThrows(IllegalArgumentException.class, () -> queue.addServerType(null));
    }

    @Test
    void testQueueSetHeteroSchedPolicy() {
        Network model = new Network("TestModel");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);

        // Default policy is ORDER
        assertEquals(HeteroSchedPolicy.ORDER, queue.getHeteroSchedPolicy());

        queue.setHeteroSchedPolicy(HeteroSchedPolicy.FSF);
        assertEquals(HeteroSchedPolicy.FSF, queue.getHeteroSchedPolicy());

        queue.setHeteroSchedPolicy(HeteroSchedPolicy.ALIS);
        assertEquals(HeteroSchedPolicy.ALIS, queue.getHeteroSchedPolicy());
    }

    @Test
    void testQueueSetHeteroService() {
        Network model = new Network("TestModel");
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");

        OpenClass jobClass = new OpenClass(model, "Jobs");
        source.setArrival(jobClass, new Exp(1.0));

        ServerType fast = new ServerType("Fast", 2);
        ServerType slow = new ServerType("Slow", 3);

        queue.addServerType(fast);
        queue.addServerType(slow);

        // Set service for fast servers
        Exp fastService = new Exp(2.0);
        queue.setService(jobClass, fast, fastService);

        // Set service for slow servers
        Exp slowService = new Exp(1.0);
        queue.setService(jobClass, slow, slowService);

        // Verify services are stored
        Distribution retrievedFast = queue.getService(jobClass, fast);
        Distribution retrievedSlow = queue.getService(jobClass, slow);

        assertNotNull(retrievedFast);
        assertNotNull(retrievedSlow);
    }

    @Test
    void testQueueSetHeteroServiceAutoCompatibility() {
        Network model = new Network("TestModel");
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");

        OpenClass jobClass = new OpenClass(model, "Jobs");
        source.setArrival(jobClass, new Exp(1.0));

        ServerType st = new ServerType("Fast", 2);
        queue.addServerType(st);

        // Initially not compatible
        assertFalse(st.isCompatible(jobClass));

        // Setting service should auto-add compatibility
        queue.setService(jobClass, st, new Exp(2.0));
        assertTrue(st.isCompatible(jobClass));
    }

    @Test
    void testQueueSetHeteroServiceServerTypeNotAdded() {
        Network model = new Network("TestModel");
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);

        OpenClass jobClass = new OpenClass(model, "Jobs");
        source.setArrival(jobClass, new Exp(1.0));

        ServerType st = new ServerType("Fast", 2);
        // Note: NOT calling queue.addServerType(st)

        assertThrows(IllegalArgumentException.class, () ->
            queue.setService(jobClass, st, new Exp(2.0)));
    }

    @Test
    void testQueueGetHeteroServiceDistributions() {
        Network model = new Network("TestModel");
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");

        OpenClass class1 = new OpenClass(model, "Class1");
        OpenClass class2 = new OpenClass(model, "Class2");
        source.setArrival(class1, new Exp(1.0));
        source.setArrival(class2, new Exp(0.5));

        ServerType fast = new ServerType("Fast", 2);
        ServerType slow = new ServerType("Slow", 3);

        queue.addServerType(fast);
        queue.addServerType(slow);

        queue.setService(class1, fast, new Exp(2.0));
        queue.setService(class1, slow, new Exp(1.0));
        queue.setService(class2, slow, new Exp(1.5));

        Map<ServerType, Map<JobClass, Distribution>> distributions = queue.getHeteroServiceDistributions();

        assertNotNull(distributions);
        assertEquals(2, distributions.size());
        assertTrue(distributions.containsKey(fast));
        assertTrue(distributions.containsKey(slow));

        // Fast server only serves class1
        assertEquals(1, distributions.get(fast).size());
        assertTrue(distributions.get(fast).containsKey(class1));

        // Slow server serves both classes
        assertEquals(2, distributions.get(slow).size());
        assertTrue(distributions.get(slow).containsKey(class1));
        assertTrue(distributions.get(slow).containsKey(class2));
    }

    // ========== Integration Tests ==========

    @Test
    void testHeteroServerModelBuilding() {
        Network model = new Network("HeteroModel");

        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "HeteroQueue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");

        OpenClass class1 = new OpenClass(model, "Class1");
        OpenClass class2 = new OpenClass(model, "Class2");

        source.setArrival(class1, new Exp(1.5));
        source.setArrival(class2, new Exp(1.0));

        // Create server types with different compatibilities
        ServerType fastServers = new ServerType("Fast", 2);
        ServerType slowServers = new ServerType("Slow", 3);

        // Fast only serves Class1, Slow serves both
        fastServers.addCompatibleClass(class1);
        slowServers.addCompatibleClass(class1);
        slowServers.addCompatibleClass(class2);

        queue.addServerType(fastServers);
        queue.addServerType(slowServers);
        queue.setHeteroSchedPolicy(HeteroSchedPolicy.FSF);

        queue.setService(class1, fastServers, new Exp(2.0));
        queue.setService(class1, slowServers, new Exp(1.0));
        queue.setService(class2, slowServers, new Exp(1.5));

        model.link(Network.serialRouting(source, queue, sink));

        // Validate model
        assertTrue(queue.isHeterogeneous());
        assertEquals(2, queue.getServerTypes().size());
        assertEquals(HeteroSchedPolicy.FSF, queue.getHeteroSchedPolicy());

        // Validate compatibility
        assertTrue(fastServers.isCompatible(class1));
        assertFalse(fastServers.isCompatible(class2));
        assertTrue(slowServers.isCompatible(class1));
        assertTrue(slowServers.isCompatible(class2));
    }

    @Test
    void testHeteroServerTotalServerCount() {
        Network model = new Network("TestModel");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);

        ServerType st1 = new ServerType("Type1", 2);
        ServerType st2 = new ServerType("Type2", 3);
        ServerType st3 = new ServerType("Type3", 5);

        queue.addServerType(st1);
        queue.addServerType(st2);
        queue.addServerType(st3);

        // Total servers should be 2 + 3 + 5 = 10
        List<ServerType> types = queue.getServerTypes();
        int totalServers = 0;
        for (ServerType st : types) {
            totalServers += st.getNumOfServers();
        }
        assertEquals(10, totalServers);
    }

    @Test
    void testCompatibilityMatrix() {
        Network model = new Network("TestModel");
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");

        OpenClass classA = new OpenClass(model, "ClassA");
        OpenClass classB = new OpenClass(model, "ClassB");
        OpenClass classC = new OpenClass(model, "ClassC");

        source.setArrival(classA, new Exp(1.0));
        source.setArrival(classB, new Exp(1.0));
        source.setArrival(classC, new Exp(1.0));

        ServerType type1 = new ServerType("Type1", 1);
        ServerType type2 = new ServerType("Type2", 1);

        // Type1 serves A and B
        type1.addCompatibleClass(classA);
        type1.addCompatibleClass(classB);

        // Type2 serves B and C
        type2.addCompatibleClass(classB);
        type2.addCompatibleClass(classC);

        queue.addServerType(type1);
        queue.addServerType(type2);

        // Verify compatibility matrix
        assertTrue(type1.isCompatible(classA));
        assertTrue(type1.isCompatible(classB));
        assertFalse(type1.isCompatible(classC));

        assertFalse(type2.isCompatible(classA));
        assertTrue(type2.isCompatible(classB));
        assertTrue(type2.isCompatible(classC));
    }
}
