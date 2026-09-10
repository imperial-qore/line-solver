package jline.lang;

import jline.examples.java.basic.ClassSwitchingModel;
import jline.lang.constant.RoutingStrategy;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Node;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.processes.Distribution;
import jline.lang.processes.Exp;
import jline.util.NamedParam;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;

import java.util.List;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Covers the MATLAB lang methods ported into the JAR to close the class-level
 * parity gap recorded in _kb/07-cross-language-parity.md: Model.attribute,
 * JobClass.summary, RoutingMatrix.getCell/set-by-index/rtnodes2rtorig,
 * Node.link/isStation/hasClassSwitching/summary/copy and the Distribution
 * parameter management.
 */
public class MatlabPortedApiTest {

    private static Network buildOpenModel() {
        Network model = new Network("ported");
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");
        OpenClass oclass = new OpenClass(model, "Class1", 0);
        source.setArrival(oclass, Exp.fitMean(2.0));
        queue.setService(oclass, Exp.fitMean(1.0));
        model.link(Network.serialRouting(oclass, source, queue, sink));
        return model;
    }

    // ========== Model ==========

    @Test
    void testModelAttributeIsAMetadataContainer() {
        Model model = new Model("plain");
        assertNotNull(model.getAttribute());
        assertFalse(model.getAttribute().has("layer"));
        model.getAttribute().put("layer", 3);
        assertTrue(model.getAttribute().has("layer"));
        assertEquals(3, model.getAttribute().get("layer"));
        assertTrue(model.getAttribute().keys().contains("layer"));
        assertEquals(3, model.getAttribute().remove("layer"));
        assertFalse(model.getAttribute().has("layer"));

        ModelAttribute replacement = new ModelAttribute();
        replacement.put("k", "v");
        model.setAttribute(replacement);
        assertEquals("v", model.getAttribute().get("k"));
    }

    @Test
    void testNetworkAttributeIsTheModelAttribute() {
        Network model = buildOpenModel();
        // Network carries the typed container, and Model must see the same object
        assertSame(model.getAttribute(), ((Model) model).getAttribute());
        model.getAttribute().put("origin", "test");
        assertEquals("test", ((Model) model).getAttribute().get("origin"));
    }

    @Test
    void testModelVersionIsTrimmed() {
        Model model = new Model("plain");
        assertNotNull(model.getVersion());
        assertEquals(model.getVersion().trim(), model.getVersion());
        assertEquals("plain", model.getName());
    }

    @Test
    void testModelIsCopyable() {
        Model model = new Model("plain");
        model.getAttribute().put("k", 7);
        Model clone = model.copy();
        assertNotSame(model, clone);
        assertEquals("plain", clone.getName());
        assertEquals(7, clone.getAttribute().get("k"));
        clone.setName("other");
        assertEquals("plain", model.getName());
    }

    // ========== JobClass ==========

    @Test
    void testJobClassSummary() {
        Network model = buildOpenModel();
        JobClass jobClass = model.getClasses().get(0);
        jobClass.summary(); // one-line form, must not throw
        jobClass.printSummary();
        assertEquals("Class1", jobClass.getName());
    }

    // ========== RoutingMatrix ==========

    @Test
    void testRoutingMatrixSetByIndexAndGetCell() {
        Network model = new Network("rm");
        Delay delay = new Delay(model, "Delay");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        ClosedClass cclass1 = new ClosedClass(model, "Class1", 1, delay);
        ClosedClass cclass2 = new ClosedClass(model, "Class2", 0, delay);
        delay.setService(cclass1, Exp.fitMean(1.0));
        delay.setService(cclass2, Exp.fitMean(1.0));
        queue.setService(cclass1, Exp.fitMean(0.5));
        queue.setService(cclass2, Exp.fitMean(0.5));

        RoutingMatrix routing = model.initRoutingMatrix();

        // set(classIndex, matrix): the whole node-by-node block of class 1.
        // The int-indexed family is 1-based, as get(int,int) and MATLAB are.
        Matrix block = new Matrix(2, 2);
        block.set(0, 1, 1.0);
        block.set(1, 0, 1.0);
        routing.set(1, block);
        assertEquals(1.0, routing.get(1, 1).get(0, 1), 1e-14);
        assertEquals(1.0, routing.get(1, 1).get(1, 0), 1e-14);

        // set(classIndex1, classIndex2, nodeIndex1, nodeIndex2, value): one entry
        routing.set(2, 2, 1, 2, 1.0);
        routing.set(2, 2, 2, 1, 1.0);
        assertEquals(1.0, routing.get(2, 2).get(0, 1), 1e-14);
        assertEquals(1.0, routing.get(2, 2, 1, 2), 1e-14);

        // a class switch declared entry by entry
        routing.set(1, 2, 2, 1, 0.0);
        assertEquals(0.0, routing.get(1, 2).get(1, 0), 1e-14);

        List<List<Matrix>> cell = routing.getCell();
        assertEquals(2, cell.size());
        assertEquals(2, cell.get(0).size());
        assertEquals(1.0, cell.get(0).get(0).get(0, 1), 1e-14);
        // the cell is returned by value, as in MATLAB
        cell.get(0).get(0).set(0, 1, 0.25);
        assertEquals(1.0, routing.get(1, 1).get(0, 1), 1e-14);
    }

    @Test
    void testRoutingMatrixSetRejectsOutOfRangeIndices() {
        Network model = buildOpenModel();
        RoutingMatrix routing = model.initRoutingMatrix();
        assertThrows(RuntimeException.class, () -> routing.set(1, 5, 1, 2, 1.0));
        assertThrows(RuntimeException.class, () -> routing.set(1, 1, 1, 99, 1.0));
        assertThrows(RuntimeException.class, () -> routing.set(0, 1, 1, 2, 1.0));
    }

    /**
     * The reference values are the nonzeros MATLAB
     * {@code RoutingMatrix.rtnodes2rtorig} returns on the same model, one row
     * per nonzero as {class1, class2, node1, node2, value} in MATLAB 1-based
     * indices. Note this is NOT {@code sn.rtorig}: the stochastic complement
     * also fills the rows of stations a class never visits, which is why
     * {@code refreshRoutingMatrix.m} does not use it to populate that field.
     */
    private static final double[][] CS_MULTI_DIAMOND_RTORIG = {
            {1, 1, 1, 2, 1.0},
            {1, 1, 2, 2, 0.2},
            {1, 1, 3, 5, 1.0},
            {1, 1, 4, 5, 1.0},
            {1, 2, 2, 3, 0.3},
            {1, 3, 2, 4, 0.5},
            {2, 2, 1, 2, 1.0},
            {2, 2, 2, 2, 1.0 / 3.0},
            {2, 2, 2, 3, 1.0 / 3.0},
            {2, 2, 2, 4, 1.0 / 3.0},
            {2, 2, 3, 5, 1.0},
            {2, 2, 4, 5, 1.0},
            {3, 3, 1, 2, 1.0},
            {3, 3, 2, 2, 1.0 / 3.0},
            {3, 3, 2, 3, 1.0 / 3.0},
            {3, 3, 2, 4, 1.0 / 3.0},
            {3, 3, 3, 5, 1.0},
            {3, 3, 4, 5, 1.0},
    };

    @Test
    void testRtnodes2rtorigMatchesMatlab() {
        Network model = ClassSwitchingModel.cs_multi_diamond();
        NetworkStruct sn = model.getStruct(true);

        // the comparison is index-by-index, so the node order must be the one
        // MATLAB reported: Source, Queue 0, Queue 1, Queue 2, Sink, CS_, CS_
        assertEquals(7, sn.nnodes);
        assertEquals(3, sn.nclasses);
        assertTrue(sn.nodenames.get(5).startsWith("CS_"), "node 6 must be the first class switch node");

        RoutingMatrix.RtOrigResult result = RoutingMatrix.rtnodes2rtorig(sn);
        assertEquals(3, result.rtorigcell.size());
        int csshift = result.rtorigcell.get(0).get(0).getNumRows();
        assertEquals(5, csshift);

        // every entry MATLAB reports as nonzero, and nothing else
        for (int r = 0; r < 3; r++) {
            for (int s = 0; s < 3; s++) {
                Matrix actual = result.rtorigcell.get(r).get(s);
                for (int i = 0; i < csshift; i++) {
                    for (int j = 0; j < csshift; j++) {
                        double want = 0.0;
                        for (double[] ref : CS_MULTI_DIAMOND_RTORIG) {
                            if ((int) ref[0] == r + 1 && (int) ref[1] == s + 1
                                    && (int) ref[2] == i + 1 && (int) ref[3] == j + 1) {
                                want = ref[4];
                            }
                        }
                        assertEquals(want, actual.get(i, j), 1e-12,
                                String.format("rtorig mismatch at class (%d,%d) node (%d,%d)",
                                        r + 1, s + 1, i + 1, j + 1));
                    }
                }
            }
        }
        assertNotNull(result.rtorig);
        assertEquals(csshift * sn.nclasses, result.rtorig.getNumRows());
    }

    // ========== Node ==========

    @Test
    void testNodePredicatesAndLink() {
        Network model = new Network("nodes");
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");
        OpenClass oclass = new OpenClass(model, "Class1", 0);
        source.setArrival(oclass, Exp.fitMean(2.0));
        queue.setService(oclass, Exp.fitMean(1.0));

        assertTrue(queue.isStation());
        assertTrue(source.isStation());
        assertFalse(sink.isStation());
        assertFalse(queue.hasClassSwitching());

        source.link(queue);
        queue.link(sink);
        assertTrue(model.getConnectionMatrix().get(model.getNodeIndex(source),
                model.getNodeIndex(queue)) > 0);

        queue.summary(); // must not throw
    }

    @Test
    void testClassSwitchNodeReportsClassSwitching() {
        Network model = ClassSwitchingModel.cs_implicit();
        boolean seen = false;
        for (Node node : model.getNodes()) {
            if (node instanceof jline.lang.nodes.ClassSwitch) {
                assertTrue(node.hasClassSwitching());
                assertFalse(node.isStation());
                seen = true;
            }
        }
        assertTrue(seen, "the model must carry a ClassSwitch node");
    }

    @Test
    void testNodeCopySharesTheModelAndClonesTheSections() {
        Network model = buildOpenModel();
        Queue queue = (Queue) model.getNodeByName("Queue");
        Queue clone = queue.copy();

        assertNotSame(queue, clone);
        assertSame(queue.model, clone.model); // model handle is shared, as in MATLAB
        assertEquals(queue.getName(), clone.getName());
        assertNotSame(queue.getOutput(), clone.getOutput());
        assertNotSame(queue.getServer(), clone.getServer());
        // the sections' contents stay shared: same job class objects
        assertSame(queue.getServer().getServiceDistribution(model.getClasses().get(0)),
                clone.getServer().getServiceDistribution(model.getClasses().get(0)));
    }

    @Test
    void testSetRoutingWithSingleParameter() {
        Network model = new Network("sq");
        Source source = new Source(model, "Source");
        Queue queue1 = new Queue(model, "Queue1", SchedStrategy.FCFS);
        Queue queue2 = new Queue(model, "Queue2", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");
        OpenClass oclass = new OpenClass(model, "Class1", 0);
        source.setArrival(oclass, Exp.fitMean(2.0));
        queue1.setService(oclass, Exp.fitMean(1.0));
        queue2.setService(oclass, Exp.fitMean(1.0));

        source.setRouting(oclass, RoutingStrategy.SQ, 2);
        source.setRouting(oclass, RoutingStrategy.RAND, null);

        assertThrows(IllegalArgumentException.class,
                () -> source.setRouting(oclass, RoutingStrategy.SQ, "two"));
        assertThrows(IllegalArgumentException.class,
                () -> source.setRouting(oclass, RoutingStrategy.SQ, 0));
    }

    // ========== Distribution ==========

    @Test
    void testDistributionParameterManagement() {
        Distribution dist = Exp.fitMean(2.0);
        assertTrue(dist.getNumParams() >= 1);
        assertTrue(dist.hasParam(1));
        assertNotNull(dist.getParams());
        assertEquals(dist.getNumParams(), dist.getParams().size());

        NamedParam first = dist.getParam(1);
        assertNotNull(first.getName());
        assertSame(first, dist.getParam(first.getName()));
        assertNull(dist.getParam("no_such_param"));

        // a matrix-valued parameter must survive as the object it is
        Matrix payload = new Matrix(2, 2);
        payload.set(0, 0, 1.5);
        dist.setParam(2, "block", payload);
        assertSame(payload, dist.getParam(2).getValue());
        assertTrue(dist.hasParam(2));
        assertFalse(dist.hasParam(9));

        assertThrows(RuntimeException.class, () -> dist.setParam(0, "bad", 1.0));
        assertThrows(RuntimeException.class, () -> dist.setParam(1, "", 1.0));
    }

    @Test
    void testDistributionIsCopyable() {
        Distribution dist = Exp.fitMean(2.0);
        Distribution clone = dist.copy();
        assertNotSame(dist, clone);
        assertEquals(dist.getMean(), clone.getMean(), 1e-14);
        assertFalse(clone.isImmediate());
    }
}
