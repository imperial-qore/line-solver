/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang;

import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.*;
import jline.lang.processes.Erlang;
import jline.lang.processes.Exp;
import jline.lang.processes.HyperExp;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;

/**
 * Tests that verify model.copy() preserves the printout format.
 * Models are created inline to avoid dependencies on the examples package.
 */
public class ModelCopyTest {

    /**
     * Helper method to capture printout from a model
     */
    private String captureModelPrint(Network model) {
        ByteArrayOutputStream baos = new ByteArrayOutputStream();
        PrintStream ps = new PrintStream(baos);
        PrintStream originalOut = System.out;

        try {
            System.setOut(ps);
            model.getStruct().print();
            return baos.toString();
        } finally {
            System.setOut(originalOut);
        }
    }

    /**
     * Creates a closed network with single customer for response time CDF analysis.
     */
    private Network createClosedModel() {
        Network model = new Network("model");
        Delay delay = new Delay(model, "Delay");
        Queue queue = new Queue(model, "Queue2", SchedStrategy.PS);
        ClosedClass jobClass = new ClosedClass(model, "Class1", 1, delay, 0);
        delay.setService(jobClass, new Exp(1.0/0.1));
        queue.setService(jobClass, Erlang.fitMeanAndSCV(1, 1.0/3.0));
        RoutingMatrix routingMatrix = model.initRoutingMatrix();
        Matrix circularRouting = new Matrix(2, 2);
        circularRouting.set(0, 1, 1.0);
        circularRouting.set(1, 0, 1.0);
        routingMatrix.set(jobClass, circularRouting);
        model.link(routingMatrix);
        return model;
    }

    /**
     * Creates a closed network with three classes and class switching.
     */
    private Network createClosedThreeclassesModel() {
        Network model = new Network("model");
        Delay delay = new Delay(model, "Delay");
        Queue queue = new Queue(model, "Queue2", SchedStrategy.PS);
        ClosedClass class1 = new ClosedClass(model, "Class1", 1, delay, 0);
        ClosedClass class2 = new ClosedClass(model, "Class2", 0, delay, 0);
        ClosedClass class3 = new ClosedClass(model, "Class3", 0, delay, 0);
        class1.setCompletes(false);
        delay.setService(class1, new Exp(1.0/1.0));
        delay.setService(class2, new Exp(1.0/1.0));
        delay.setService(class3, new Exp(1.0/1.0));
        queue.setService(class1, new Exp(1.0/1.0));
        queue.setService(class2, Erlang.fitMeanAndOrder(4.0, 2));
        queue.setService(class3, new Exp(1.0/0.01));
        RoutingMatrix routingMatrix = model.initRoutingMatrix();
        Matrix p11 = new Matrix(2, 2);
        p11.set(0, 1, 1.0);
        p11.set(1, 0, 0.0);
        routingMatrix.set(class1, class1, p11);
        Matrix p12 = new Matrix(2, 2);
        p12.set(0, 0, 0.0);
        p12.set(1, 0, 1.0);
        routingMatrix.set(class1, class2, p12);
        Matrix p21 = new Matrix(2, 2);
        p21.set(0, 0, 0.0);
        p21.set(1, 0, 1.0);
        routingMatrix.set(class2, class1, p21);
        Matrix p22 = new Matrix(2, 2);
        p22.set(0, 1, 1.0);
        p22.set(1, 0, 0.0);
        routingMatrix.set(class2, class2, p22);
        Matrix p33 = new Matrix(2, 2);
        p33.set(0, 1, 1.0);
        p33.set(1, 0, 1.0);
        routingMatrix.set(class3, p33);
        model.link(routingMatrix);
        return model;
    }

    /**
     * Creates an open network with two classes.
     */
    private Network createOpenTwoclassesModel() {
        Network model = new Network("myModel");
        Source source = new Source(model, "Source");
        Queue queue1 = new Queue(model, "Queue1", SchedStrategy.FCFS);
        Queue queue2 = new Queue(model, "Queue2", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");
        OpenClass class1 = new OpenClass(model, "Class1", 0);
        OpenClass class2 = new OpenClass(model, "Class2", 0);
        source.setArrival(class1, Exp.fitMean(4.0));
        source.setArrival(class2, Exp.fitMean(4.0));
        queue1.setService(class1, Exp.fitMean(1.0));
        queue1.setService(class2, Exp.fitMean(1.0));
        queue2.setService(class1, Exp.fitMean(1.0));
        queue2.setService(class2, Exp.fitMean(1.0));
        RoutingMatrix routingMatrix = model.initRoutingMatrix();
        routingMatrix.set(class1, class1, source, queue1, 1.0);
        routingMatrix.set(class1, class2, queue1, queue2, 1.0);
        routingMatrix.set(class2, class1, queue2, sink, 1.0);
        routingMatrix.set(class2, class2, source, queue1, 1.0);
        routingMatrix.set(class2, class1, queue1, queue2, 1.0);
        routingMatrix.set(class1, class2, queue2, sink, 1.0);
        model.link(routingMatrix);
        return model;
    }

    /**
     * Creates a closed network with different service distributions.
     */
    private Network createDistribModel() {
        Network model = new Network("model");
        Delay delay = new Delay(model, "Delay");
        Queue queue = new Queue(model, "Queue1", SchedStrategy.PS);
        ClosedClass class1 = new ClosedClass(model, "Class1", 1, delay, 0);
        ClosedClass class2 = new ClosedClass(model, "Class2", 3, delay, 0);
        delay.setService(class1, Exp.fitMean(1.0));
        queue.setService(class1, Exp.fitMean(2.0));
        delay.setService(class2, Erlang.fitMeanAndOrder(4.0, 2));
        queue.setService(class2, HyperExp.fitMeanAndSCV(5.0, 30.0));
        RoutingMatrix routingMatrix = model.initRoutingMatrix();
        routingMatrix.set(class1, Network.serialRouting(delay, queue));
        routingMatrix.set(class2, Network.serialRouting(delay, queue));
        model.link(routingMatrix);
        return model;
    }

    /**
     * Creates a closed network with population N=4.
     */
    private Network createPopulationsN4Model() {
        Network model = new Network("model");
        Delay delay = new Delay(model, "Delay");
        Queue queue1 = new Queue(model, "Queue1", SchedStrategy.PS);
        Queue queue2 = new Queue(model, "Queue2", SchedStrategy.PS);
        ClosedClass jobClass = new ClosedClass(model, "Class1", 4, delay, 0);
        jobClass.setCompletes(false);
        delay.setService(jobClass, new Exp(1.0/1.0));
        queue1.setService(jobClass, new Exp(1.0/2.0));
        queue2.setService(jobClass, new Exp(1.0/2.0));
        RoutingMatrix routingMatrix = model.initRoutingMatrix();
        Matrix circularRouting = new Matrix(3, 3);
        circularRouting.set(0, 1, 1.0);
        circularRouting.set(1, 2, 1.0);
        circularRouting.set(2, 0, 1.0);
        routingMatrix.set(jobClass, circularRouting);
        model.link(routingMatrix);
        return model;
    }

    @Test
    public void testCdfRespTClosedCopy() {
        Network original = createClosedModel();
        assertNotNull(original, "Original model should not be null");

        Network copy = original.copy();
        assertNotNull(copy, "Copied model should not be null");

        String originalPrint = captureModelPrint(original);
        String copyPrint = captureModelPrint(copy);

        assertEquals(originalPrint, copyPrint, "Original and copied model printouts should be identical");
    }

    @Test
    public void testCdfRespTClosedThreeclassesCopy() {
        Network original = createClosedThreeclassesModel();
        assertNotNull(original, "Original model should not be null");

        Network copy = original.copy();
        assertNotNull(copy, "Copied model should not be null");

        String originalPrint = captureModelPrint(original);
        String copyPrint = captureModelPrint(copy);

        assertEquals(originalPrint, copyPrint, "Original and copied model printouts should be identical");
    }

    @Test
    public void testCdfRespTOpenTwoclassesCopy() {
        Network original = createOpenTwoclassesModel();
        assertNotNull(original, "Original model should not be null");

        Network copy = original.copy();
        assertNotNull(copy, "Copied model should not be null");

        String originalPrint = captureModelPrint(original);
        String copyPrint = captureModelPrint(copy);

        assertEquals(originalPrint, copyPrint, "Original and copied model printouts should be identical");
    }

    @Test
    public void testCdfRespTDistribCopy() {
        Network original = createDistribModel();
        assertNotNull(original, "Original model should not be null");

        Network copy = original.copy();
        assertNotNull(copy, "Copied model should not be null");

        String originalPrint = captureModelPrint(original);
        String copyPrint = captureModelPrint(copy);

        assertEquals(originalPrint, copyPrint, "Original and copied model printouts should be identical");
    }

    @Test
    public void testCdfRespTPopulationsCopy() {
        Network original = createPopulationsN4Model();
        assertNotNull(original, "Original model should not be null");

        Network copy = original.copy();
        assertNotNull(copy, "Copied model should not be null");

        String originalPrint = captureModelPrint(original);
        String copyPrint = captureModelPrint(copy);

        assertEquals(originalPrint, copyPrint, "Original and copied model printouts should be identical");
    }
}