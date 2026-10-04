/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.layered;

import jline.lang.constant.SchedStrategy;

/**
 * FunctionTask is the former name of SetupTask, retained for backward compatibility.
 *
 * <p>Setup and delay-off times are not specific to serverless (function-as-a-service)
 * platforms, so the class carrying them is now named after the modelling primitive
 * rather than after that application domain.
 *
 * @see SetupTask
 * @deprecated Use SetupTask, or a plain Task with setSetupTime() and setDelayOffTime().
 */
@Deprecated
public class FunctionTask extends SetupTask {

    /**
     * Constructor for FunctionTask.
     *
     * @param model The LayeredNetwork model this task belongs to
     * @param name The name of the task
     * @param multiplicity The number of servers that can run concurrently
     * @param scheduling The scheduling strategy for the servers
     */
    public FunctionTask(LayeredNetwork model, String name, int multiplicity, SchedStrategy scheduling) {
        super(model, name, multiplicity, scheduling);
    }

    /**
     * Constructor with default scheduling strategy (FCFS).
     */
    public FunctionTask(LayeredNetwork model, String name, int multiplicity) {
        this(model, name, multiplicity, SchedStrategy.FCFS);
    }

    /**
     * Constructor with default multiplicity and scheduling.
     */
    public FunctionTask(LayeredNetwork model, String name) {
        this(model, name, 1, SchedStrategy.FCFS);
    }

    @Override
    public FunctionTask on(Processor parent) {
        super.on(parent);
        return this;
    }

    @Override
    public FunctionTask setThinkTime(jline.lang.processes.Distribution thinkTime) {
        super.setThinkTime(thinkTime);
        return this;
    }

    @Override
    public FunctionTask setThinkTime(double thinkTime) {
        super.setThinkTime(thinkTime);
        return this;
    }

    @Override
    public FunctionTask setSetupTime(jline.lang.processes.Distribution setupTime) {
        super.setSetupTime(setupTime);
        return this;
    }

    @Override
    public FunctionTask setSetupTime(double setupTime) {
        super.setSetupTime(setupTime);
        return this;
    }

    @Override
    public FunctionTask setDelayOffTime(jline.lang.processes.Distribution delayOffTime) {
        super.setDelayOffTime(delayOffTime);
        return this;
    }

    @Override
    public FunctionTask setDelayOffTime(double delayOffTime) {
        super.setDelayOffTime(delayOffTime);
        return this;
    }
}
