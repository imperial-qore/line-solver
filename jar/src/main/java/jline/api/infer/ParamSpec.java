package jline.api.infer;

/**
 * Specification of a hidden LQN parameter to identify.
 *
 * <p>Used by {@link InferLqn} to name which Layered Queueing Network parameters
 * are estimated: an activity host demand ({@code HOSTDEM}) or a task think time
 * ({@code THINK}).</p>
 *
 * Copyright (c) 2012-2026, Imperial College London. All rights reserved.
 */
public class ParamSpec {

    /** Kind of parameter being estimated. */
    public enum Type { HOSTDEM, THINK }

    public final Type type;
    public final String name;

    public ParamSpec(Type type, String name) {
        this.type = type;
        this.name = name;
    }

    /** Estimate the mean host demand of the named activity. */
    public static ParamSpec hostDemand(String activityName) {
        return new ParamSpec(Type.HOSTDEM, activityName);
    }

    /** Estimate the mean think time of the named task. */
    public static ParamSpec think(String taskName) {
        return new ParamSpec(Type.THINK, taskName);
    }
}
