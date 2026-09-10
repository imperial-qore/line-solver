/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.layered;

import jline.lang.constant.SchedStrategy;

/**
 * SetupTask is a Task whose servers are switched off while idle.
 *
 * <p>A server that resumes from the off state pays a setup (activation) time
 * before serving the request that woke it up, and it stays available for a
 * delay-off (idle) period after emptying its queue before switching off again.
 * These are the setup and close-down times of a server with vacations; they
 * model on-demand virtual machines and containers, power-managed servers under
 * a timeout policy, warm-up delays, and serverless cold start / keep-alive.
 *
 * <p>Both times are declared on the base Task class through setSetupTime() and
 * setDelayOffTime(), so this subclass is a naming convenience: a plain Task
 * carrying either time behaves identically.
 *
 * @see Task
 */
public class SetupTask extends Task {

    /**
     * Constructor for SetupTask.
     *
     * @param model The LayeredNetwork model this task belongs to
     * @param name The name of the task
     * @param multiplicity The number of servers that can run concurrently
     * @param scheduling The scheduling strategy for the servers
     */
    public SetupTask(LayeredNetwork model, String name, int multiplicity, SchedStrategy scheduling) {
        super(model, name, multiplicity, scheduling);
    }

    /**
     * Constructor with default scheduling strategy (FCFS).
     */
    public SetupTask(LayeredNetwork model, String name, int multiplicity) {
        this(model, name, multiplicity, SchedStrategy.FCFS);
    }

    /**
     * Constructor with default multiplicity and scheduling.
     */
    public SetupTask(LayeredNetwork model, String name) {
        this(model, name, 1, SchedStrategy.FCFS);
    }

    // Override methods to return SetupTask for method chaining

    @Override
    public SetupTask on(Processor parent) {
        super.on(parent);
        return this;
    }

    @Override
    public SetupTask setThinkTime(jline.lang.processes.Distribution thinkTime) {
        super.setThinkTime(thinkTime);
        return this;
    }

    @Override
    public SetupTask setThinkTime(double thinkTime) {
        super.setThinkTime(thinkTime);
        return this;
    }

    @Override
    public SetupTask setSetupTime(jline.lang.processes.Distribution setupTime) {
        super.setSetupTime(setupTime);
        return this;
    }

    @Override
    public SetupTask setSetupTime(double setupTime) {
        super.setSetupTime(setupTime);
        return this;
    }

    @Override
    public SetupTask setDelayOffTime(jline.lang.processes.Distribution delayOffTime) {
        super.setDelayOffTime(delayOffTime);
        return this;
    }

    @Override
    public SetupTask setDelayOffTime(double delayOffTime) {
        super.setDelayOffTime(delayOffTime);
        return this;
    }
}
