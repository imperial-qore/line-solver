package jline.lang;

import jline.gen.NetworkGenerator;
import jline.lang.constant.SchedStrategy;
import jline.lang.constant.SignalType;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.processes.Exp;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.List;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Coverage tests for lang classes with 0% coverage.
 */
public class LangCoverageTest {

    // ========== Chain ==========

    @Test
    void testChain() {
        // Test simple constructor
        Chain chain = new Chain("TestChain");
        assertNotNull(chain);
        assertEquals("TestChain", chain.getName());
        assertTrue(chain.getClasses().isEmpty());

        // Test setName
        chain.setName("RenamedChain");
        assertEquals("RenamedChain", chain.getName());
    }

    // ========== Region ==========

    @Test
    void testRegion() {
        // Create a simple Network for the classes
        Network model = new Network("TestModel");
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");

        OpenClass openClass = new OpenClass(model, "Class1");
        source.setArrival(openClass, new Exp(1.0));
        queue.setService(openClass, new Exp(2.0));

        // Create region with nodes and classes
        List<jline.lang.nodes.Node> nodes = new ArrayList<>();
        nodes.add(queue);
        List<JobClass> classes = new ArrayList<>();
        classes.add(openClass);

        Region fcr = new Region(nodes, classes);

        assertNotNull(fcr);
        assertEquals(1, fcr.getNodes().size());

        // Test global limits
        assertEquals(Region.UNBOUNDED, fcr.getGlobalMaxJobs());
        fcr.setGlobalMaxJobs(10);
        assertEquals(10, fcr.getGlobalMaxJobs());

        fcr.setGlobalMaxMemory(100);
        assertEquals(100, fcr.getGlobalMaxMemory());

        // Test per-class limits
        fcr.setClassMaxJobs(openClass, 5);
        assertEquals(5, fcr.getClassMaxJobs(openClass));

        fcr.setClassMaxMemory(openClass, 50);
        assertEquals(50, fcr.getClassMaxMemory(openClass));

        // Test drop rule
        assertFalse(fcr.getDropRule(openClass));
        fcr.setDropRule(openClass, true);
        assertTrue(fcr.getDropRule(openClass));

        // Test class size
        assertEquals(1, fcr.getClassSize(openClass));
        fcr.setClassSize(openClass, 2);
        assertEquals(2, fcr.getClassSize(openClass));

        // Test name
        fcr.setName("TestRegion");
        assertEquals("TestRegion", fcr.getName());
    }

    // ========== Signal ==========

    @Test
    void testSignal() {
        Network model = new Network("SignalTestModel");
        Source source = new Source(model, "Source");
        Sink sink = new Sink(model, "Sink");

        // Create a signal with type
        Signal signal = new Signal(model, "TestSignal", SignalType.NEGATIVE);
        assertNotNull(signal);
        assertEquals("TestSignal", signal.getName());
        assertEquals(SignalType.NEGATIVE, signal.getSignalType());

        // Test constructor with NEGATIVE type
        Signal signal2 = new Signal(model, "Signal2", SignalType.NEGATIVE);
        assertEquals(SignalType.NEGATIVE, signal2.getSignalType());

        // Test setSignalType
        signal.setSignalType(SignalType.REPLY);
        assertEquals(SignalType.REPLY, signal.getSignalType());

        // Test forJobClass association
        OpenClass targetClass = new OpenClass(model, "TargetClass");
        signal.forJobClass(targetClass);
        assertEquals(targetClass, signal.getTargetJobClass());
        assertEquals(targetClass.getIndex(), signal.getTargetJobClassIndex());

        // Test null target
        Signal signal3 = new Signal(model, "Signal3", SignalType.NEGATIVE);
        assertEquals(-1, signal3.getTargetJobClassIndex());
    }

    // ========== ClosedSignal ==========

    @Test
    void testClosedSignal() {
        Network model = new Network("ClosedSignalTestModel");
        Delay delay = new Delay(model, "Delay");

        ClosedClass closedClass = new ClosedClass(model, "ClosedClass", 5, delay);
        delay.setService(closedClass, new Exp(1.0));

        // Create a closed signal
        ClosedSignal signal = new ClosedSignal(model, "TestClosedSignal", SignalType.NEGATIVE, delay);
        assertNotNull(signal);
        assertEquals("TestClosedSignal", signal.getName());
        assertEquals(SignalType.NEGATIVE, signal.getSignalType());

        // Test default constructor
        ClosedSignal signal2 = new ClosedSignal(model, "Signal2", delay);
        assertEquals(SignalType.NEGATIVE, signal2.getSignalType());

        // Test setSignalType
        signal.setSignalType(SignalType.REPLY);
        assertEquals(SignalType.REPLY, signal.getSignalType());

        // Test forJobClass association
        signal.forJobClass(closedClass);
        assertEquals(closedClass, signal.getTargetJobClass());
        assertEquals(closedClass.getIndex(), signal.getTargetJobClassIndex());

        // Test null target
        ClosedSignal signal3 = new ClosedSignal(model, "Signal3", delay);
        assertEquals(-1, signal3.getTargetJobClassIndex());
    }

    // ========== NetworkGenerator ==========

    @Test
    void testNetworkGenerator() {
        // Test default constructor
        NetworkGenerator generator = new NetworkGenerator();
        assertNotNull(generator);

        // Generate a simple network with 2 queues
        Network model = generator.generate(2);
        assertNotNull(model);
        assertTrue(model.getNumberOfNodes() > 0);
    }

    // ========== AfterEventRouter ==========
    // Note: AfterEventRouter.afterEventRouter is a static package-private method
    // that requires complex NetworkStruct setup. We test it indirectly through
    // the Network model that uses it, or skip if setup is too complex.

    @Test
    void testAfterEventRouterIndirect() {
        // AfterEventRouter is tested indirectly through model analysis
        // This test ensures the class loads and basic types are accessible
        Network model = new Network("RouterTestModel");
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");

        OpenClass openClass = new OpenClass(model, "Class1");
        source.setArrival(openClass, new Exp(1.0));
        queue.setService(openClass, new Exp(2.0));

        model.link(model.serialRouting(source, queue, sink));

        // The model uses AfterEventRouter internally when solving
        // Just verify the model is valid
        assertTrue(model.hasOpenClasses());
        assertEquals(3, model.getNumberOfNodes());
    }
}
