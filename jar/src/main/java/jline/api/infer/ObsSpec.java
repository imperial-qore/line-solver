package jline.api.infer;

/**
 * Specification of an observed LQN performance metric.
 *
 * <p>Used by {@link InferLqn} to select which measured quantity forms a
 * component of the observation vector z = h(a): a response time, a utilization,
 * a throughput or a queue length of a named LQN element.</p>
 *
 * Copyright (c) 2012-2026, Imperial College London. All rights reserved.
 */
public class ObsSpec {

    /** Kind of performance metric observed. */
    public enum Metric { RESPT, UTIL, TPUT, QLEN }

    public final Metric metric;
    public final String name;

    public ObsSpec(Metric metric, String name) {
        this.metric = metric;
        this.name = name;
    }

    /** Observe the response time of the named element (entry or reference task). */
    public static ObsSpec respT(String name) {
        return new ObsSpec(Metric.RESPT, name);
    }

    /** Observe the utilization of the named element (typically a processor). */
    public static ObsSpec util(String name) {
        return new ObsSpec(Metric.UTIL, name);
    }

    /** Observe the throughput of the named element. */
    public static ObsSpec tput(String name) {
        return new ObsSpec(Metric.TPUT, name);
    }

    /** Observe the queue length of the named element. */
    public static ObsSpec qLen(String name) {
        return new ObsSpec(Metric.QLEN, name);
    }
}
