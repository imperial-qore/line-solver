/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api.sn;

import jline.io.Ret;
import jline.lang.ClosedClass;
import jline.lang.Network;
import jline.lang.NetworkStruct;
import jline.lang.OpenClass;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Delay;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.processes.Exp;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Validation of the NetworkStruct utility APIs (jline.api.sn) on models with
 * known structure:
 * - the SnHas* predicates must reflect the model topology exactly;
 * - product-form parameter extraction must return the constructed demands.
 */
public class SnApiTest {

    private static final double TOL = 1e-12;

    /** Closed PS+Delay cycle: 3 jobs, think 1.0, demand 0.5. */
    private static Network closedModel() {
        Network model = new Network("sn_closed");
        Delay delay = new Delay(model, "Delay");
        Queue queue = new Queue(model, "Queue", SchedStrategy.PS);
        ClosedClass jobClass = new ClosedClass(model, "Class1", 3, delay, 0);
        delay.setService(jobClass, new Exp(1.0));
        queue.setService(jobClass, new Exp(2.0));
        model.link(model.serialRouting(delay, queue));
        return model;
    }

    /** Open M/M/2 with two priority classes. */
    private static Network openModel() {
        Network model = new Network("sn_open");
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        queue.setNumberOfServers(2);
        Sink sink = new Sink(model, "Sink");
        OpenClass class1 = new OpenClass(model, "Class1", 0);
        OpenClass class2 = new OpenClass(model, "Class2", 1);
        source.setArrival(class1, new Exp(0.5));
        source.setArrival(class2, new Exp(0.3));
        queue.setService(class1, new Exp(2.0));
        queue.setService(class2, new Exp(3.0));
        model.link(model.serialRouting(source, queue, sink));
        return model;
    }

    @Test
    public void predicatesReflectClosedModelTopology() {
        NetworkStruct sn = closedModel().getStruct(true);
        assertTrue(SnHasClosedClasses.snHasClosedClasses(sn), "model has a closed class");
        assertFalse(SnHasOpenClasses.snHasOpenClasses(sn), "model has no open class");
        assertTrue(SnHasINF.snHasINF(sn), "Delay station is INF scheduling");
        assertFalse(SnHasFCFS.snHasFCFS(sn), "no FCFS station in PS+Delay cycle");
        assertFalse(SnHasDPS.snHasDPS(sn), "no DPS station");
        assertFalse(SnHasForkJoin.snHasForkJoin(sn), "no fork-join");
        assertFalse(SnHasClassSwitching.snHasClassSwitching(sn), "no class switching");
        // Infinite servers are delays, not multiserver queues: the predicate
        // ranges over finite nservers only (as the Python sn_has_multi_server
        // does)
        assertFalse(SnHasMultiServer.snHasMultiServer(sn),
                "Delay is INF, and the PS queue is single-server");
        assertFalse(SnHasPriorities.snHasPriorities(sn), "single class: no priorities");
        assertTrue(SnHasProductForm.snHasProductForm(sn),
                "PS+Delay closed cycle is product form");
    }

    /** cqn_bas_blocking: BAS at Queue1, a one-job buffer at Queue2, two jobs circulating. */
    private static Network basModel() {
        Network model = new Network("sn_bas");
        Queue q1 = new Queue(model, "Queue1", SchedStrategy.FCFS);
        Queue q2 = new Queue(model, "Queue2", SchedStrategy.FCFS);
        ClosedClass jobClass = new ClosedClass(model, "Class1", 2, q1, 0);
        q1.setService(jobClass, new Exp(1.0));
        q2.setService(jobClass, new Exp(0.8));
        q2.setCap(1);
        q1.setDropRule(jobClass, jline.lang.constant.DropStrategy.BlockingAfterService);
        model.link(model.serialRouting(q1, q2));
        return model;
    }

    /** Single-station M/M/1/K with tail drop, the one truncated shape that keeps product form. */
    private static Network mm1kModel() {
        Network model = new Network("sn_mm1k");
        Source source = new Source(model, "Source");
        Queue queue = new Queue(model, "Queue", SchedStrategy.FCFS);
        queue.setNumberOfServers(1);
        queue.setCapacity(3);
        Sink sink = new Sink(model, "Sink");
        OpenClass jobClass = new OpenClass(model, "Class1", 0);
        source.setArrival(jobClass, new Exp(0.8));
        queue.setService(jobClass, new Exp(1.0));
        model.link(model.serialRouting(source, queue, sink));
        return model;
    }

    @Test
    public void bindingBufferExcludesProductForm() {
        // Until 2026-08-16 no conjunct of snHasProductForm read sn.cap/sn.classcap/
        // sn.droprule, so this model reported a product form its truncation does not have
        NetworkStruct sn = basModel().getStruct(true);
        assertTrue(SnHasBlocking.snHasBlocking(sn), "a one-job buffer holding two jobs binds");
        assertFalse(SnHasProductForm.snHasProductForm(sn), "BAS blocking is not product form");
    }

    @Test
    public void nonBindingCapacityKeepsProductForm() {
        // capacity at least the reachable population can never refuse a job: setCap(N)
        // on an N-job closed model is a no-op, and product form stays exact
        Network model = closedModel();
        ((Queue) model.getNodeByName("Queue")).setCap(3);
        NetworkStruct sn = model.getStruct(true);
        assertFalse(SnHasBlocking.snHasBlocking(sn), "capacity 3 with 3 jobs never binds");
        assertTrue(SnHasProductForm.snHasProductForm(sn), "PS+Delay closed cycle stays product form");
    }

    @Test
    public void singleStationLossSystemIsExempt() {
        NetworkStruct sn = mm1kModel().getStruct(true);
        assertTrue(SnIsMm1kLoss.snIsMm1kLoss(sn), "Source-Queue-Sink M/M/1/K with tail drop");
        assertFalse(SnHasBlocking.snHasBlocking(sn), "its truncated geometric IS a product form");
        assertTrue(SnHasProductForm.snHasProductForm(sn));
    }

    @Test
    public void predicatesReflectOpenModelTopology() {
        NetworkStruct sn = openModel().getStruct(true);
        assertTrue(SnHasOpenClasses.snHasOpenClasses(sn), "model has open classes");
        assertFalse(SnHasClosedClasses.snHasClosedClasses(sn), "no closed class");
        assertTrue(SnHasFCFS.snHasFCFS(sn), "queue is FCFS");
        assertTrue(SnHasMultiServer.snHasMultiServer(sn), "queue has 2 servers");
        assertTrue(SnHasPriorities.snHasPriorities(sn), "classes have distinct priorities");
        assertFalse(SnHasForkJoin.snHasForkJoin(sn), "no fork-join");
    }

    @Test
    public void productFormParamsMatchConstruction() {
        NetworkStruct sn = closedModel().getStruct(true);
        Ret.snGetProductFormParams params =
                SnGetProductFormParams.snGetProductFormParams(sn);
        assertNotNull(params, "product form params must not be null");
        // Single chain with 3 jobs
        assertEquals(3.0, params.N.get(0), TOL, "chain population");
        // Think time at the Delay: 1.0; queue demand: 0.5
        assertEquals(1.0, params.Z.elementSum(), TOL, "aggregate think time");
        assertEquals(0.5, params.D.elementSum(), TOL, "aggregate service demand");
    }

    @Test
    public void chainDemandsAggregateServiceDemands() {
        NetworkStruct sn = closedModel().getStruct(true);
        Ret.snGetDemands demands = SnGetDemandsChain.snGetDemandsChain(sn);
        assertNotNull(demands, "chain demands must not be null");
        assertNotNull(demands.Dchain, "chain demand matrix must not be null");
        // One chain visiting Delay (think 1.0) and Queue (demand 0.5): the
        // chain demand across stations sums to 1.5
        assertEquals(1.5, demands.Dchain.elementSum() + 0.0, 1e-9,
                "total chain demand must equal think time plus service demand");
    }
}
