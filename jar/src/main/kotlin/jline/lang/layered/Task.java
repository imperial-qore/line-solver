/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.layered;

import jline.GlobalConstants;
import jline.lang.constant.SchedStrategy;
import jline.lang.processes.Distribution;
import jline.lang.processes.Exp;
import jline.lang.processes.Immediate;

import java.util.ArrayList;
import java.util.List;

import static jline.io.InputOutputKt.line_error;
import static jline.io.InputOutputKt.mfilename;

/**
 * A Task represents a software component or process in a layered queueing network that can host
 * services (Entry objects) and runs on a Processor (Host).
 * 
 * <p>Tasks are fundamental building blocks in layered queueing networks that encapsulate:
 * <ul>
 * <li>Service entries that define the interfaces accessible to other tasks</li>
 * <li>Activities that represent internal processing steps</li>
 * <li>Activity precedences that define the execution order and relationships</li>
 * <li>Resource demands and scheduling policies</li>
 * </ul>
 * 
 * <p>Key characteristics:
 * <ul>
 * <li><b>Multiplicity:</b> Number of concurrent task instances that can execute</li>
 * <li><b>Scheduling:</b> How the task schedules its internal processing (FCFS, PS, etc.)</li>
 * <li><b>Think time:</b> Time between completing one request and starting the next</li>
 * <li><b>Reference tasks:</b> Special tasks that generate workload (external customers)</li>
 * </ul>
 * 
 * <p>Tasks can be reference tasks (SchedStrategy.REF) that represent external workload
 * generators, or regular tasks that provide services to other tasks in the system.
 * 
 * @see Entry
 * @see Activity
 * @see ActivityPrecedence
 * @see Processor
 * @see LayeredNetwork
 */
public class Task extends LayeredNetworkElement {
    protected Processor parent;
    protected int multiplicity;
    protected int replication;
    protected SchedStrategy scheduling;//Enum
    protected int priority = 0;  // Priority level (0 = default/no priority)
    protected String fanInSource = "";  // Source task for fan-in
    protected int fanInValue = 0;  // Fan-in value (load distribution count)
    protected Distribution thinkTime;
    protected double thinkTimeMean;
    protected double thinkTimeSCV;
    protected Distribution setupTime;
    protected double setupTimeMean;
    protected double setupTimeSCV;
    protected Distribution delayOffTime;
    protected double delayOffTimeMean;
    protected double delayOffTimeSCV;
    protected List<Entry> entries;
    protected List<Activity> activities;
    protected List<ActivityPrecedence> precedences;
    private List<Entry> replyEntry;

    public Task(LayeredNetwork model, String name, int multiplicity, SchedStrategy scheduling, Distribution thinkTime) {
        super(name);
        this.parent = null;
        this.setReplication(1);
        this.model = model;
        this.multiplicity = multiplicity;
        this.scheduling = scheduling;
        this.setThinkTime(thinkTime);
        this.setSetupTime(Immediate.getInstance());
        this.setDelayOffTime(Immediate.getInstance());
        this.entries = new ArrayList<>();
        this.activities = new ArrayList<>();
        this.precedences = new ArrayList<>();
        this.replyEntry = new ArrayList<>();
        // link within model
        model.tasks.put(model.tasks.size(), this);
        model.nodes.put(model.tasks.size(), this);
        if (scheduling == SchedStrategy.REF) {
            model.reftasks.put(model.reftasks.size(), this);
        }
    }

    public Task(LayeredNetwork model, String name, int multiplicity, SchedStrategy scheduling) {
        this(model, name, multiplicity, scheduling, Immediate.getInstance());
    }

    public Task(LayeredNetwork model, String name, int multiplicity) {
        this(model, name, multiplicity, SchedStrategy.INF, Immediate.getInstance());
    }

    public Task(LayeredNetwork model, String name) {
        this(model, name, 1, SchedStrategy.INF, Immediate.getInstance());
    }

    public Task addActivity(Activity newActivity) {
        newActivity.setParent(this);
        this.activities.add(newActivity);
        return this;
    }

    public Task addEntry(Entry newEntry) {
        this.entries.add(newEntry);
        return this;
    }

    public Task addPrecedence(ActivityPrecedence newPrec) {
        this.precedences.add(newPrec);
        return this;
    }

    public Task addPrecedence(List<ActivityPrecedence> newPrec) {
        this.precedences.addAll(newPrec);
        return this;
    }

    public Task addPrecedence(ActivityPrecedence[] newPrec) {
        for (ActivityPrecedence ap : newPrec) {
            this.precedences.add(ap);
        }
        return this;
    }

    public double getMeanHostDemand(String entryName) {
        double meanHostDemand = -1;
        for (int j = 0; j < this.entries.size(); j++) {
            if (this.entries.get(j).getName().equals(entryName)) {
                break;
            }
        }
        return meanHostDemand;
    }

    public Task on(Processor parent) {
        if (this.parent != null) {
            line_error(mfilename(new Object() {
            }), "Parent processor already defined.");
        }
        this.parent = parent;
        this.parent.addTask(this);
        return this;
    }

    public Task removeActivity(int index) {
        this.activities.remove(index);
        return this;
    }

    public Task setActivity(Activity newActivity, int index) {
        this.activities.set(index, newActivity);
        return this;
    }

    public Task setAsReferenceTask() {
        this.scheduling = SchedStrategy.REF;
        return this;
    }

    public Task setReplication(int replication) {
        this.replication = replication;
        return this;
    }

    public Task setReplyEntry(List<Entry> replyEntry) {
        this.replyEntry = replyEntry;
        return this;
    }

    public Task setThinkTime(Distribution thinkTime) {
        this.thinkTime = thinkTime;
        this.thinkTimeMean = thinkTime.getMean();
        this.thinkTimeSCV = thinkTime.getSCV();
        return this;
    }

    public Task setThinkTime(double thinkTime) {
        if (thinkTime <= GlobalConstants.Zero) {
            this.thinkTime = Immediate.getInstance();
            this.thinkTimeMean = GlobalConstants.Zero;
            this.thinkTimeSCV = GlobalConstants.Zero;
        } else {
            this.thinkTime = new Exp(1 / thinkTime);
            this.thinkTimeMean = thinkTime;
            this.thinkTimeSCV = 1.0;
        }
        return this;
    }

    /**
     * Set the setup time (cold start time) for the task.
     *
     * @param setupTime The setup time distribution
     * @return This Task instance for method chaining
     */
    public Task setSetupTime(Distribution setupTime) {
        this.setupTime = setupTime;
        this.setupTimeMean = setupTime.getMean();
        this.setupTimeSCV = setupTime.getSCV();
        return this;
    }

    /**
     * Set the setup time using a mean value (creates exponential distribution).
     *
     * @param setupTime The mean setup time
     * @return This Task instance for method chaining
     */
    public Task setSetupTime(double setupTime) {
        if (setupTime <= GlobalConstants.FineTol) {
            this.setupTime = Immediate.getInstance();
            this.setupTimeMean = GlobalConstants.FineTol;
            this.setupTimeSCV = GlobalConstants.FineTol;
        } else {
            this.setupTime = new Exp(1.0 / setupTime);
            this.setupTimeMean = setupTime;
            this.setupTimeSCV = 1.0;
        }
        return this;
    }

    /**
     * Set the delay-off time (teardown time) for the task.
     *
     * @param delayOffTime The delay-off time distribution
     * @return This Task instance for method chaining
     */
    public Task setDelayOffTime(Distribution delayOffTime) {
        this.delayOffTime = delayOffTime;
        this.delayOffTimeMean = delayOffTime.getMean();
        this.delayOffTimeSCV = delayOffTime.getSCV();
        return this;
    }

    /**
     * Set the delay-off time using a mean value (creates exponential distribution).
     *
     * @param delayOffTime The mean delay-off time
     * @return This Task instance for method chaining
     */
    public Task setDelayOffTime(double delayOffTime) {
        if (delayOffTime <= GlobalConstants.FineTol) {
            this.delayOffTime = Immediate.getInstance();
            this.delayOffTimeMean = GlobalConstants.FineTol;
            this.delayOffTimeSCV = GlobalConstants.FineTol;
        } else {
            this.delayOffTime = new Exp(1.0 / delayOffTime);
            this.delayOffTimeMean = delayOffTime;
            this.delayOffTimeSCV = 1.0;
        }
        return this;
    }

    // Getters for setup time
    public Distribution getSetupTime() {
        return setupTime;
    }

    public double getSetupTimeMean() {
        return setupTimeMean;
    }

    public double getSetupTimeSCV() {
        return setupTimeSCV;
    }

    // Getters for delay-off time
    public Distribution getDelayOffTime() {
        return delayOffTime;
    }

    public double getDelayOffTimeMean() {
        return delayOffTimeMean;
    }

    public double getDelayOffTimeSCV() {
        return delayOffTimeSCV;
    }

    public int getMultiplicity() {
        return multiplicity;
    }

    public int getReplication() {
        return replication;
    }

    public SchedStrategy getScheduling() {
        return scheduling;
    }

    public double getThinkTimeMean() {
        return thinkTimeMean;
    }

    public double getThinkTimeSCV() {
        return thinkTimeSCV;
    }

    public Host getParent() {
        return parent;
    }

    public void setMultiplicity(int multiplicity) {
        this.multiplicity = multiplicity;
    }

    public void setScheduling(SchedStrategy scheduling) {
        this.scheduling = scheduling;
    }

    public List<ActivityPrecedence> getPrecedences() {
        return precedences;
    }

    /**
     * Check if this task has setup/delayoff configured (i.e., non-trivial values).
     *
     * @return true if setup time or delay-off time is configured with non-immediate values
     */
    public boolean hasSetupDelayoff() {
        boolean hasSetup = setupTime != null && !(setupTime instanceof Immediate)
                           && setupTimeMean > GlobalConstants.FineTol;
        boolean hasDelayoff = delayOffTime != null && !(delayOffTime instanceof Immediate)
                              && delayOffTimeMean > GlobalConstants.FineTol;
        return hasSetup || hasDelayoff;
    }

    /**
     * Returns the list of entries for this task.
     *
     * @return the list of entries
     */
    public List<Entry> getEntries() {
        return entries;
    }

    /**
     * Returns the list of activities for this task.
     *
     * @return the list of activities
     */
    public List<Activity> getActivities() {
        return activities;
    }

    /**
     * Returns the processor (host) that this task runs on.
     *
     * @return the parent processor
     */
    public Processor getProcessor() {
        return parent;
    }

    public Task setPriority(int priority) {
        this.priority = priority;
        return this;
    }

    public int getPriority() {
        return this.priority;
    }

    public Task setFanIn(String source, int value) {
        this.fanInSource = source;
        this.fanInValue = value;
        return this;
    }

    public String getFanInSource() {
        return this.fanInSource;
    }

    public int getFanInValue() {
        return this.fanInValue;
    }

}
