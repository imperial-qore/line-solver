/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.layered;

import jline.lang.constant.SchedStrategy;

/**
 * FunctionTask is an alias for Task, provided for backward compatibility.
 *
 * <p>All setup/delayoff functionality has been moved to the base Task class.
 * Any Task can now have setup time (cold start delay) and delay-off time
 * (teardown delay) configured via setSetupTime() and setDelayOffTime().
 *
 * <p>This class is retained for backward compatibility with existing code
 * that uses FunctionTask to model serverless functions or tasks with
 * initialization overhead.
 *
 * @see Task
 * @deprecated Use Task directly with setSetupTime() and setDelayOffTime() methods.
 */
@Deprecated
public class FunctionTask extends Task {

    /**
     * Constructor for FunctionTask.
     *
     * @param model The LayeredNetwork model this task belongs to
     * @param name The name of the function task
     * @param multiplicity The number of function instances that can run concurrently
     * @param scheduling The scheduling strategy for the function instances
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

    // Override methods to return FunctionTask for method chaining (backward compatibility)

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
