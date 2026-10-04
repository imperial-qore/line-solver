/*
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 */
package jline.inference.lang;

import jline.lang.JobClass;
import jline.lang.nodes.Node;

import java.io.Serializable;
import java.util.Objects;

/**
 * Conditioning event for sampled metrics.
 */
public final class ConditionEvent implements Serializable {
    private static final long serialVersionUID = 1L;

    public final Node node;
    public final JobClass jobClass;
    public final Object event;

    public ConditionEvent(Node node, JobClass jobClass, Object event) {
        this.node = node;
        this.jobClass = jobClass;
        this.event = event;
    }

    public Node getNode() { return node; }
    public JobClass getJobClass() { return jobClass; }
    public Object getEvent() { return event; }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof ConditionEvent)) return false;
        ConditionEvent that = (ConditionEvent) o;
        return Objects.equals(node, that.node) && Objects.equals(jobClass, that.jobClass)
                && Objects.equals(event, that.event);
    }

    @Override
    public int hashCode() {
        return Objects.hash(node, jobClass, event);
    }

    @Override
    public String toString() {
        return "ConditionEvent(node=" + node + ", class=" + jobClass + ", event=" + event + ")";
    }
}
