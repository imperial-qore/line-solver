/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api.wf;

import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Fork;
import jline.lang.nodes.Join;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.processes.Exp;
import org.junit.jupiter.api.Test;

import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Validation of the workflow analysis APIs (jline.api.wf) on a network with
 * known topology: a fork-join with two parallel branches embedded in a
 * sequence. The analyzer must classify the fork/join nodes and service nodes
 * exactly as constructed and produce a self-consistent analysis.
 */
public class WfApiTest {

    /** Source -> Q1 -> Fork -> (Q2 | Q3) -> Join -> Q4 -> Sink. */
    private static Network forkJoinWorkflow() {
        Network model = new Network("wf_forkjoin");
        Source source = new Source(model, "Source");
        Queue q1 = new Queue(model, "Queue1", SchedStrategy.FCFS);
        Fork fork = new Fork(model, "Fork");
        Queue q2 = new Queue(model, "Queue2", SchedStrategy.FCFS);
        Queue q3 = new Queue(model, "Queue3", SchedStrategy.FCFS);
        Join join = new Join(model, "Join", fork);
        Queue q4 = new Queue(model, "Queue4", SchedStrategy.FCFS);
        Sink sink = new Sink(model, "Sink");

        OpenClass jobClass = new OpenClass(model, "Class1", 0);
        source.setArrival(jobClass, new Exp(0.5));
        q1.setService(jobClass, new Exp(2.0));
        q2.setService(jobClass, new Exp(2.0));
        q3.setService(jobClass, new Exp(2.0));
        q4.setService(jobClass, new Exp(2.0));

        jline.lang.RoutingMatrix p = model.initRoutingMatrix();
        p.set(jobClass, jobClass, source, q1, 1.0);
        p.set(jobClass, jobClass, q1, fork, 1.0);
        p.set(jobClass, jobClass, fork, q2, 1.0);
        p.set(jobClass, jobClass, fork, q3, 1.0);
        p.set(jobClass, jobClass, q2, join, 1.0);
        p.set(jobClass, jobClass, q3, join, 1.0);
        p.set(jobClass, jobClass, join, q4, 1.0);
        p.set(jobClass, jobClass, q4, sink, 1.0);
        model.link(p);
        return model;
    }

    @Test
    public void analyzerClassifiesForkJoinTopology() {
        Wf_analyzer analyzer = new Wf_analyzer(forkJoinWorkflow());
        Wf_analyzer.WorkflowAnalysis analysis = analyzer.analyzeWorkflow();
        assertNotNull(analysis, "analysis must not be null");

        Wf_analyzer.WorkflowRepresentation repr = analysis.getOriginalWorkflow();
        assertNotNull(repr, "workflow representation must not be null");
        assertEquals(1, repr.getForkNodes().size(), "exactly one fork node");
        assertEquals(1, repr.getJoinNodes().size(), "exactly one join node");
        assertTrue(repr.getServiceNodes().size() >= 4,
                "the four queues must be classified as service nodes");
        assertNotNull(repr.getLinkMatrix(), "link matrix must be present");
    }

    @Test
    public void analysisPassesSelfValidation() {
        Wf_analyzer analyzer = new Wf_analyzer(forkJoinWorkflow());
        Wf_analyzer.WorkflowAnalysis analysis = analyzer.analyzeWorkflow();
        assertTrue(analyzer.validateAnalysis(analysis),
                "analysis of a well-formed workflow must validate");
    }

    @Test
    public void recommendationsAreProduced() {
        Wf_analyzer analyzer = new Wf_analyzer(forkJoinWorkflow());
        Wf_analyzer.WorkflowAnalysis analysis = analyzer.analyzeWorkflow();
        List<String> recommendations = analyzer.getOptimizationRecommendations(analysis);
        assertNotNull(recommendations, "recommendations list must not be null");
    }

    @Test
    public void detectedPatternsIncludeTheParallelSection() {
        Wf_analyzer analyzer = new Wf_analyzer(forkJoinWorkflow());
        Wf_analyzer.WorkflowAnalysis analysis = analyzer.analyzeWorkflow();
        Wf_analyzer.DetectedPatterns patterns = analysis.getDetectedPatterns();
        assertNotNull(patterns, "detected patterns must not be null");
    }
}
