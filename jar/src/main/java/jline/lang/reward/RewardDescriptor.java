/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.reward;

import jline.lang.JobClass;
import jline.lang.NetworkStruct;
import jline.lang.nodes.Node;
import jline.util.matrix.Matrix;

/**
 * A {@link RewardFunction} that also records what the reward measures.
 *
 * A RewardDescriptor delegates evaluation to the wrapped reward function, so it can be
 * passed to {@code Network.setReward} exactly like a plain lambda, while additionally
 * exposing the structural description ({@link #getKind()}, {@link #getNode()},
 * {@link #getJobClass()}) that allows the reward to be serialized declaratively.
 *
 * A descriptor of kind {@link Kind#Custom} wraps an arbitrary user function and is
 * deliberately NOT serializable: the writer warns and omits it rather than emitting a
 * reward it cannot reproduce on reload.
 *
 * @see Reward
 */
public class RewardDescriptor implements RewardFunction {

    private static final long serialVersionUID = 1L;

    /** Structural kind of a reward. The names double as the JSON {@code type} values. */
    public enum Kind {
        QLen,
        Util,
        Blocking,
        Custom
    }

    private final Kind kind;
    private final Node node;
    private final JobClass jobClass;
    private final RewardFunction delegate;

    /**
     * Creates a reward descriptor.
     *
     * @param kind     structural kind of the reward
     * @param node     node the reward refers to, or null when not node-scoped
     * @param jobClass job class the reward refers to, or null when not class-scoped
     * @param delegate the reward function actually evaluated
     */
    public RewardDescriptor(Kind kind, Node node, JobClass jobClass, RewardFunction delegate) {
        if (kind == null) {
            throw new IllegalArgumentException("Reward kind must not be null");
        }
        if (delegate == null) {
            throw new IllegalArgumentException("Reward descriptor must wrap a reward function");
        }
        this.kind = kind;
        this.node = node;
        this.jobClass = jobClass;
        this.delegate = delegate;
    }

    @Override
    public double compute(Matrix state, NetworkStruct sn) {
        return this.delegate.compute(state, sn);
    }

    /**
     * @return the structural kind of this reward
     */
    public Kind getKind() {
        return this.kind;
    }

    /**
     * @return the node this reward refers to, or null when not node-scoped
     */
    public Node getNode() {
        return this.node;
    }

    /**
     * @return the job class this reward refers to, or null when not class-scoped
     */
    public JobClass getJobClass() {
        return this.jobClass;
    }

    /**
     * @return the wrapped reward function
     */
    public RewardFunction getDelegate() {
        return this.delegate;
    }

    @Override
    public String toString() {
        String nodeName = (this.node == null) ? "null" : this.node.getName();
        String className = (this.jobClass == null) ? "null" : this.jobClass.getName();
        return "RewardDescriptor(kind=" + this.kind + ", node=" + nodeName + ", jobClass=" + className + ")";
    }
}
