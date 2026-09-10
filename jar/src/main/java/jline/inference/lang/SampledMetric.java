/*
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 */
package jline.inference.lang;

import jline.lang.Copyable;
import jline.lang.JobClass;
import jline.lang.constant.MetricType;
import jline.lang.nodes.Node;

/**
 * Observed data for a metric.
 */
public class SampledMetric implements Copyable {
    private static final long serialVersionUID = 1L;

    public final MetricType type;
    public double[] t;
    public double[] data;
    public final Node node;
    public final JobClass jobClass;

    public ConditionEvent cond;
    private SampledFormat format = SampledFormat.TIMESERIES;

    public SampledMetric(MetricType type, double[] t, double[] data, Node node, JobClass jobClass) {
        this.type = type;
        this.t = t;
        this.data = data;
        this.node = node;
        this.jobClass = jobClass;
    }

    public SampledMetric(MetricType type, double[] t, double[] data, Node node) {
        this(type, t, data, node, null);
    }

    public MetricType getType() { return type; }
    public double[] getT() { return t; }
    public void setT(double[] t) { this.t = t; }
    public double[] getData() { return data; }
    public void setData(double[] data) { this.data = data; }
    public Node getNode() { return node; }
    public JobClass getJobClass() { return jobClass; }
    public ConditionEvent getCond() { return cond; }
    public SampledFormat getFormat() { return format; }

    /** Set a conditioning event on this metric. */
    public void setConditional(ConditionEvent event) {
        this.cond = event;
    }

    /** Switch the data format to trace (per-request). */
    public void setTrace() {
        this.format = SampledFormat.TRACE;
    }

    /** Check if this metric applies to all classes (aggregate). */
    public boolean isAggregate() {
        return jobClass == null;
    }

    /** Check if this metric has a conditioning event. */
    public boolean isConditional() {
        return cond != null;
    }

    /** Check if this metric is in trace format. */
    public boolean isTrace() {
        return format == SampledFormat.TRACE;
    }

    @Override
    @SuppressWarnings("unchecked")
    public <T extends Copyable> T copy() {
        SampledMetric c = new SampledMetric(type, t.clone(), data.clone(), node, jobClass);
        c.cond = cond;
        c.format = format;
        return (T) c;
    }
}
