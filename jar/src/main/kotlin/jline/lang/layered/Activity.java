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
import jline.util.matrix.Matrix;

import java.util.HashMap;
import java.util.Locale;
import java.util.Map;
import java.util.Objects;

import static jline.io.InputOutputKt.line_error;
import static jline.io.InputOutputKt.mfilename;

/**
 * An Activity represents an individual processing step or service operation within a Task
 * in a layered queueing network.
 * 
 * <p>Activities are the fundamental units of work that define what processing happens
 * when a request is served by a task. They encapsulate:
 * <ul>
 * <li>Resource demands (CPU time, I/O operations, etc.)</li>
 * <li>Service calls to other tasks (synchronous and asynchronous)</li>
 * <li>Processing logic and execution flow</li>
 * <li>Binding relationships with entries that trigger them</li>
 * </ul>
 * 
 * <p>Key characteristics:
 * <ul>
 * <li><b>Host demand:</b> The computational resources required on the processor</li>
 * <li><b>Service calls:</b> Requests made to other tasks during execution</li>
 * <li><b>Call ordering:</b> Whether calls are made deterministically or stochastically</li>
 * <li><b>Entry binding:</b> Which entry point triggers this activity</li>
 * <li><b>Reply activity:</b> Whether this activity sends replies back to callers</li>
 * </ul>
 * 
 * <p>Activities can make both synchronous calls (blocking until response) and
 * asynchronous calls (fire-and-forget) to other task entries. They are connected
 * through ActivityPrecedence relationships to form complex execution graphs.
 * 
 * @see Task
 * @see Entry
 * @see ActivityPrecedence
 * @see Distribution
 */
public class Activity extends LayeredNetworkElement {
    protected Distribution hostDemand;
    protected double hostDemandMean;
    protected double hostDemandSCV;
    protected Task parent;//
    protected String boundToEntry;
    protected String callOrder;
    protected Map<Integer, String> syncCallDests = new HashMap<>();
    protected Matrix syncCallMeans = new Matrix(1, 1, 0);
    protected Map<Integer, String> asyncCallDests = new HashMap<>();
    protected Matrix asyncCallMeans = new Matrix(1, 1, 0);
    protected Matrix scheduling = new Matrix(0, 0, 0);
    protected Distribution thinkTime;
    protected double thinkTimeMean;
    protected double thinkTimeSCV;
    protected int phase = 1;  // Phase number (1 or 2), default=1

    public Activity(LayeredNetwork model, String name, Distribution hostDemand, String boundToEntry, String callOrder) {
        super(name);

        this.setHostDemand(hostDemand);
        this.boundToEntry = boundToEntry;
        this.setCallOrder(callOrder);

        // Initialize default think time
        this.thinkTime = Immediate.getInstance();
        this.thinkTimeMean = GlobalConstants.Zero;
        this.thinkTimeSCV = GlobalConstants.Zero;

        model.activities.put(model.activities.size(), this);
        model.nodes.put(model.activities.size(), this);
        this.model = model;
    }

    public Activity(LayeredNetwork model, String name, Distribution hostDemand, String boundToEntry) {
        this(model, name, hostDemand, boundToEntry, "STOCHASTIC");
    }

    public Activity(LayeredNetwork model, String name, Distribution hostDemand) {
        this(model, name, hostDemand, "", "STOCHASTIC");
    }

    public Activity(LayeredNetwork model, String name) {
        this(model, name, Immediate.getInstance(), "", "STOCHASTIC");
    }

    public Activity asynchCall(Entry asynchCallDest, double asynchCallMean) {
        this.asyncCallDests.put(asyncCallDests.size(), asynchCallDest.getName());
        if (this.asyncCallDests.size() == 1) {
            this.asyncCallMeans.set(0, asynchCallMean);
        } else {
            this.asyncCallMeans = this.asyncCallMeans.concatCols(Matrix.singleton(asynchCallMean));
        }
        return this;
    }

    public Activity asynchCall(String asynchCallDest, double asynchCallMean) {
        this.asyncCallDests.put(asyncCallDests.size(), asynchCallDest);
        if (this.asyncCallDests.size() == 1) {
            this.asyncCallMeans.set(0, asynchCallMean);
        } else {
            this.asyncCallMeans = this.asyncCallMeans.concatCols(Matrix.singleton(asynchCallMean));
        }
        return this;
    }

    public Activity asynchCall(Entry asynchCallDest) {
        this.asyncCallDests.put(asyncCallDests.size(), asynchCallDest.getName());
        if (this.asyncCallDests.size() == 1) {
            this.asyncCallMeans.set(0, 1.0);
        } else {
            this.asyncCallMeans = this.asyncCallMeans.concatCols(Matrix.singleton(1.0));
        }
        return this;
    }

    public Activity asynchCall(String asynchCallDest) {
        this.asyncCallDests.put(asyncCallDests.size(), asynchCallDest);
        if (this.asyncCallDests.size() == 1) {
            this.asyncCallMeans.set(0, 1.0);
        } else {
            this.asyncCallMeans = this.asyncCallMeans.concatCols(Matrix.singleton(1.0));
        }
        return this;
    }

    public Activity boundTo(Entry entry) {
        if (this.parent == null) {
            line_error(mfilename(new Object() {
            }), "Activity must have a parent task before binding to an entry. Use activity.on(task) first.");
        }
        if (entry.getParent() != null && entry.getParent() != this.parent) {
            line_error(mfilename(new Object() {
            }), "Activity and entry must belong to the same task. Activity parent: " +
                this.parent.getName() + ", Entry parent: " + entry.getParent().getName());
        }
        this.boundToEntry = entry.getName();
        // Find this activity's index in the model
        Integer activityIndex = findActivityIndex();
        if (activityIndex != null) {
            entry.boundToActivity.put(activityIndex, this.getName());
        }
        return this;
    }

    public Activity boundTo(String entry) {
        this.boundToEntry = entry;
        return this;
    }

    public Activity on(Task parent) {
        if (this.parent != null) {
            line_error(mfilename(new Object() {
            }), "Parent task already defined.");
        }
        parent.addActivity(this);
        this.parent = parent;
        return this;
    }

    public Activity repliesTo(Entry entry) {
        Integer activityIndex = findActivityIndex();
        if (activityIndex != null) {
            if (this.parent != null) {
                if (Objects.requireNonNull(this.parent.scheduling) == SchedStrategy.REF) {
                    line_error(mfilename(new Object() {
                    }), "Activities in reference tasks cannot reply.");
                } else {
                    entry.replyActivity.put(activityIndex, this.getName());
                }
            } else {
                entry.replyActivity.put(activityIndex, this.getName());
            }
        }
        return this;
    }

    public Activity setCallOrder(String callOrder) {
        if (callOrder.equals("STOCHASTIC") || callOrder.equals("DETERMINISTIC")) {
            this.callOrder = callOrder.toUpperCase(Locale.ROOT);
        } else {
            this.callOrder = "STOCHASTIC";
        }

        return this;
    }

    public void setHostDemand(double hostDemand) {
        if (hostDemand <= GlobalConstants.Zero) {
            this.hostDemand = Immediate.getInstance();
            this.hostDemandMean = 1e-8;
            this.hostDemandSCV = 1e-8;
        } else {
            this.hostDemand = new Exp(1 / hostDemand);
            this.hostDemandMean = hostDemand;
            this.hostDemandSCV = 1.0;
        }
    }

    public void setHostDemand(Distribution hostDemand) {
        this.hostDemand = hostDemand;
        this.hostDemandMean = hostDemand.getMean();
        this.hostDemandSCV = hostDemand.getSCV();
    }

    public void setParent(Task parent) {
        this.parent = parent;
    }

    public void setThinkTime(Distribution thinkTime) {
        this.thinkTime = thinkTime;
        this.thinkTimeMean = thinkTime.getMean();
        this.thinkTimeSCV = thinkTime.getSCV();
    }

    public void setThinkTime(double thinkTime) {
        if (thinkTime <= GlobalConstants.Zero) {
            this.thinkTime = Immediate.getInstance();
            this.thinkTimeMean = GlobalConstants.Zero;
            this.thinkTimeSCV = GlobalConstants.Zero;
        } else {
            this.thinkTime = new Exp(1 / thinkTime);
            this.thinkTimeMean = thinkTime;
            this.thinkTimeSCV = 1.0;
        }
    }

    public Distribution getThinkTime() {
        return thinkTime;
    }

    public double getThinkTimeMean() {
        return thinkTimeMean;
    }

    public double getThinkTimeSCV() {
        return thinkTimeSCV;
    }

    public Distribution getHostDemand() {
        return hostDemand;
    }

    public double getHostDemandMean() {
        return hostDemandMean;
    }

    public double getHostDemandSCV() {
        return hostDemandSCV;
    }

    public String getCallOrder() {
        return callOrder;
    }

    public String getBoundToEntry() {
        return boundToEntry;
    }

    /**
     * Get the phase number for this activity.
     * Phase 1: activities before the reply is sent
     * Phase 2: activities after the reply is sent (post-reply processing)
     * @return The phase number (1 or 2)
     */
    public int getPhase() {
        return phase;
    }

    /**
     * Set the phase number for this activity.
     * Phase 1: activities before the reply is sent
     * Phase 2: activities after the reply is sent (post-reply processing)
     * @param phase The phase number (must be 1 or 2)
     * @return This activity for method chaining
     */
    public Activity setPhase(int phase) {
        if (phase < 1 || phase > 2) {
            line_error(mfilename(new Object() {}), "Phase must be 1 or 2.");
        }
        this.phase = phase;
        return this;
    }

    public Activity synchCall(Entry synchCallDest, double synchCallMean) {
        this.syncCallDests.put(syncCallDests.size(), synchCallDest.getName());
        if (this.syncCallDests.size() == 1) {
            this.syncCallMeans.set(0, synchCallMean);
        } else {
            this.syncCallMeans = this.syncCallMeans.concatCols(Matrix.singleton(synchCallMean));
        }
        return this;
    }

    public Activity synchCall(String synchCallDest, double synchCallMean) {
        this.syncCallDests.put(syncCallDests.size(), synchCallDest);
        if (this.syncCallDests.size() == 1) {
            this.syncCallMeans.set(0, synchCallMean);
        } else {
            this.syncCallMeans = this.syncCallMeans.concatCols(Matrix.singleton(synchCallMean));
        }
        return this;
    }

    public Activity synchCall(Entry synchCallDest) {
        this.syncCallDests.put(syncCallDests.size(), synchCallDest.getName());
        if (this.syncCallDests.size() == 1) {
            this.syncCallMeans.set(0, 1.0);
        } else {
            this.syncCallMeans = this.syncCallMeans.concatCols(Matrix.singleton(1.0));
        }
        return this;
    }

    public Activity synchCall(String synchCallDest) {
        this.syncCallDests.put(syncCallDests.size(), synchCallDest);
        if (this.syncCallDests.size() == 1) {
            this.syncCallMeans.set(0, 1.0);
        } else {
            this.syncCallMeans = this.syncCallMeans.concatCols(Matrix.singleton(1.0));
        }
        return this;
    }

    public Task getParent() {
        return parent;
    }

    public Map<Integer, String> getSyncCallDests() {
        return syncCallDests;
    }

    public Matrix getSyncCallMeans() {
        return syncCallMeans;
    }

    public Activity setSyncCallMeans(Matrix syncCallMeans) {
        if (syncCallMeans == null) {
            throw new IllegalArgumentException("syncCallMeans cannot be null");
        }
        this.syncCallMeans = syncCallMeans;
        return this;
    }

    public Activity setSyncCallMeans(double syncCallMean) {
        if (this.syncCallDests.size() == 0) {
            throw new IllegalStateException("Cannot set syncCallMeans: no synchronous call destinations defined");
        }
        if (this.syncCallDests.size() == 1) {
            this.syncCallMeans.set(0, syncCallMean);
        } else {
            this.syncCallMeans = this.syncCallMeans.concatCols(Matrix.singleton(syncCallMean)); // same logic as in synchCall()
        }
        return this;
    }

    public Map<Integer, String> getAsyncCallDests() {
        return asyncCallDests;
    }

    public Matrix getAsyncCallMeans() {
        return asyncCallMeans;
    }

    public Activity setAsyncCallMeans(Matrix asyncCallMeans) {
        if (asyncCallMeans == null) {
            throw new IllegalArgumentException("asyncCallMeans cannot be null");
        }
        this.asyncCallMeans = asyncCallMeans;
        return this;
    }

    public Activity setAsyncCallMeans(double asyncCallMean) {
        if (this.asyncCallDests.size() == 0) {
            throw new IllegalStateException("Cannot set asyncCallMeans: no asynchronous call destinations defined");
        }
        if (this.asyncCallDests.size() == 1) {
            this.asyncCallMeans.set(0, asyncCallMean);
        } else {
            this.asyncCallMeans = this.asyncCallMeans.concatCols(Matrix.singleton(asyncCallMean)); // same logic as in asynchCall()
        }
        return this;
    }

    /**
     * Helper method to find this activity's index in the model
     * @return The index of this activity in the model, or null if not found
     */
    private Integer findActivityIndex() {
        if (this.model != null) {
            for (Map.Entry<Integer, Activity> entry : this.model.activities.entrySet()) {
                if (entry.getValue() == this) {
                    return entry.getKey();
                }
            }
        }
        return null;
    }

}
