package jline.lang;

import jline.lang.constant.ActivityPrecedenceType;
import jline.lang.constant.EventType;
import jline.lang.constant.NodeType;
import jline.lang.constant.ImpatienceType;
import jline.lang.constant.SchedStrategy;
import jline.lang.constant.SolverType;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.processes.Exp;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Coverage tests for constants and simple classes with 0% coverage.
 */
public class ConstantsCoverageTest {

    // ========== ActivityPrecedenceType ==========

    @Test
    void testActivityPrecedenceType() {
        // Test ID constants
        assertEquals(1, ActivityPrecedenceType.ID_PRE_SEQ);
        assertEquals(2, ActivityPrecedenceType.ID_PRE_AND);
        assertEquals(3, ActivityPrecedenceType.ID_PRE_OR);
        assertEquals(11, ActivityPrecedenceType.ID_POST_SEQ);
        assertEquals(12, ActivityPrecedenceType.ID_POST_AND);
        assertEquals(13, ActivityPrecedenceType.ID_POST_OR);
        assertEquals(14, ActivityPrecedenceType.ID_POST_LOOP);
        assertEquals(15, ActivityPrecedenceType.ID_POST_CACHE);

        // Test String constants
        assertEquals("pre", ActivityPrecedenceType.PRE_SEQ);
        assertEquals("pre-AND", ActivityPrecedenceType.PRE_AND);
        assertEquals("pre-OR", ActivityPrecedenceType.PRE_OR);
        assertEquals("post", ActivityPrecedenceType.POST_SEQ);
        assertEquals("post-AND", ActivityPrecedenceType.POST_AND);
        assertEquals("post-OR", ActivityPrecedenceType.POST_OR);
        assertEquals("post-LOOP", ActivityPrecedenceType.POST_LOOP);
        assertEquals("post-CACHE", ActivityPrecedenceType.POST_CACHE);
    }

    // ========== ImpatienceType ==========

    @Test
    void testImpatienceType() {
        // Test enum values exist
        assertNotNull(ImpatienceType.RENEGING);
        assertNotNull(ImpatienceType.BALKING);
        assertNotNull(ImpatienceType.RETRIAL);

        // Test getID
        assertEquals(1, ImpatienceType.RENEGING.getID());
        assertEquals(2, ImpatienceType.BALKING.getID());
        assertEquals(3, ImpatienceType.RETRIAL.getID());

        // Test fromID
        assertEquals(ImpatienceType.RENEGING, ImpatienceType.fromID(1));
        assertEquals(ImpatienceType.BALKING, ImpatienceType.fromID(2));
        assertEquals(ImpatienceType.RETRIAL, ImpatienceType.fromID(3));

        // Test fromID with invalid ID
        assertThrows(IllegalArgumentException.class, () -> ImpatienceType.fromID(99));

        // Test toText
        assertEquals("reneging", ImpatienceType.toText(ImpatienceType.RENEGING));
        assertEquals("balking", ImpatienceType.toText(ImpatienceType.BALKING));
        assertEquals("retrial", ImpatienceType.toText(ImpatienceType.RETRIAL));

        // Test values()
        ImpatienceType[] values = ImpatienceType.values();
        assertEquals(3, values.length);
    }

    // ========== DisabledClass ==========

    @Test
    void testDisabledClass() throws Exception {
        // Create a simple network
        Network model = new Network("DisabledTestModel");
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");

        // Create an open class first for arrivals
        OpenClass openClass = new OpenClass(model, "OpenClass");
        source.setArrival(openClass, new Exp(1.0));
        queue.setService(openClass, new Exp(2.0));

        model.link(model.serialRouting(source, queue, sink));

        // Create a disabled class - it loops at the queue
        DisabledClass disabledClass = new DisabledClass(model, "DisabledClass", queue);

        assertNotNull(disabledClass);
        assertEquals("DisabledClass", disabledClass.getName());
        assertEquals(queue, disabledClass.getReferenceStation());
    }

    // ========== SelfLoopingClass ==========

    @Test
    void testSelfLoopingClass() {
        // Create a simple closed network
        Network model = new Network("SelfLoopTestModel");
        Delay delay = new Delay(model, "Delay");

        // Create a self-looping class with priority
        SelfLoopingClass selfLooping1 = new SelfLoopingClass(model, "SelfLoop1", 5, delay, 1);
        assertNotNull(selfLooping1);
        assertEquals("SelfLoop1", selfLooping1.getName());
        assertEquals(delay, selfLooping1.getReferenceStation());
        assertEquals(5, selfLooping1.getPopulation());

        // Create a self-looping class without priority (default 0)
        Network model2 = new Network("SelfLoopTestModel2");
        Delay delay2 = new Delay(model2, "Delay2");
        SelfLoopingClass selfLooping2 = new SelfLoopingClass(model2, "SelfLoop2", 10, delay2);
        assertNotNull(selfLooping2);
        assertEquals(10, selfLooping2.getPopulation());
    }

    // ========== Env (deprecated wrapper for Environment) ==========

    @Test
    void testEnv() {
        // Test deprecated Env class (wrapper for Environment)
        Env env = new Env("TestEnv", 2);
        assertNotNull(env);
        assertEquals("TestEnv", env.getName());
        assertNotNull(env.env);  // The environment array
        assertEquals(2, env.env.length);
    }

    // ========== SolverType ==========

    @Test
    void testSolverType() {
        // Test all enum values exist
        SolverType[] values = SolverType.values();
        assertEquals(14, values.length);

        // Test specific values
        assertNotNull(SolverType.AUTO);
        assertNotNull(SolverType.CTMC);
        assertNotNull(SolverType.DES);
        assertNotNull(SolverType.ENV);
        assertNotNull(SolverType.FLUID);
        assertNotNull(SolverType.JMT);
        assertNotNull(SolverType.LN);
        assertNotNull(SolverType.LQNS);
        assertNotNull(SolverType.MAM);
        assertNotNull(SolverType.MVA);
        assertNotNull(SolverType.NC);
        assertNotNull(SolverType.POSTERIOR);
        assertNotNull(SolverType.QNS);
        assertNotNull(SolverType.SSA);

        // Test valueOf
        assertEquals(SolverType.MVA, SolverType.valueOf("MVA"));
        assertEquals(SolverType.CTMC, SolverType.valueOf("CTMC"));

        // Test name
        assertEquals("AUTO", SolverType.AUTO.name());
        assertEquals("DES", SolverType.DES.name());

        // Test ordinal
        assertEquals(0, SolverType.AUTO.ordinal());
    }

    // ========== NodeType ==========

    @Test
    void testNodeType() {
        // Test all enum values exist
        NodeType[] values = NodeType.values();
        assertEquals(13, values.length);

        // Test specific values
        assertNotNull(NodeType.Transition);
        assertNotNull(NodeType.Place);
        assertNotNull(NodeType.Fork);
        assertNotNull(NodeType.Router);
        assertNotNull(NodeType.Cache);
        assertNotNull(NodeType.Logger);
        assertNotNull(NodeType.ClassSwitch);
        assertNotNull(NodeType.Delay);
        assertNotNull(NodeType.Source);
        assertNotNull(NodeType.Sink);
        assertNotNull(NodeType.Join);
        assertNotNull(NodeType.Queue);
        assertNotNull(NodeType.Region);

        // Test valueOf
        assertEquals(NodeType.Queue, NodeType.valueOf("Queue"));
        assertEquals(NodeType.Delay, NodeType.valueOf("Delay"));
        assertEquals(NodeType.Source, NodeType.valueOf("Source"));

        // Test name
        assertEquals("Queue", NodeType.Queue.name());
        assertEquals("Fork", NodeType.Fork.name());
        assertEquals("Join", NodeType.Join.name());

        // Test ordinal
        assertEquals(0, NodeType.Transition.ordinal());
    }

    // ========== EventType ==========

    @Test
    void testEventType() {
        // Test all enum values exist
        EventType[] values = EventType.values();
        assertEquals(11, values.length);

        // Test specific values
        assertNotNull(EventType.INIT);
        assertNotNull(EventType.LOCAL);
        assertNotNull(EventType.ARV);
        assertNotNull(EventType.DEP);
        assertNotNull(EventType.PHASE);
        assertNotNull(EventType.READ);
        assertNotNull(EventType.STAGE);
        assertNotNull(EventType.ENABLE);
        assertNotNull(EventType.FIRE);
        assertNotNull(EventType.PRE);
        assertNotNull(EventType.POST);

        // Test valueOf
        assertEquals(EventType.INIT, EventType.valueOf("INIT"));
        assertEquals(EventType.ARV, EventType.valueOf("ARV"));
        assertEquals(EventType.DEP, EventType.valueOf("DEP"));

        // Test name
        assertEquals("INIT", EventType.INIT.name());
        assertEquals("ARV", EventType.ARV.name());
        assertEquals("DEP", EventType.DEP.name());

        // Test ordinal
        assertEquals(0, EventType.INIT.ordinal());
        assertEquals(2, EventType.ARV.ordinal());
        assertEquals(3, EventType.DEP.ordinal());
    }
}
