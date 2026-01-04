/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.layered;

import jline.lang.processes.Distribution;
import jline.util.matrix.Matrix;

import java.util.HashMap;
import java.util.Map;

/**
 * An Entry represents a service interface exposed by a Task in a layered queueing network.
 *
 * <p>Entries define the points where external requests can enter a task and receive service.
 * They act as the public interface of a task, similar to methods or endpoints in software systems.
 * Each entry can receive requests from other tasks or from external arrival sources.
 *
 * <p>Key characteristics:
 * <ul>
 * <li><b>Service interface:</b> Defines what service the task provides to callers</li>
 * <li><b>Arrival distribution:</b> For open systems, the arrival process of external requests</li>
 * <li><b>Activity binding:</b> Can be bound to specific activities that implement the service</li>
 * <li><b>Reply activity:</b> Designates which activity sends the response back to the caller</li>
 * </ul>
 *
 * <p>Entries are essential for modeling service-oriented architectures, web services,
 * database systems, and any system where tasks provide distinct services to clients.
 *
 * @see Task
 * @see Activity
 * @see LayeredNetwork
 */
public class Entry extends LayeredNetworkElement {
    private Distribution arrival;
    protected Task parent;
    protected String type = "PH1PH2";  // Entry type (default: PH1PH2, can be NONE, etc.)
    protected Map<Integer, String> boundToActivity = new HashMap<Integer, String>();
    protected Map<Integer, String> replyActivity = new HashMap<Integer, String>();
    protected Matrix scheduling = new Matrix(0, 0, 0);
    protected Map<Integer, String> forwardingDests = new HashMap<Integer, String>();
    protected Matrix forwardingProbs = new Matrix(0, 0, 0);

    public Entry(LayeredNetwork model, String name) {
        super(name);
        this.arrival = null;
        model.entries.put(model.entries.size(), this);
        model.nodes.put(model.entries.size(), this);
        this.model = model;
    }

    public Entry on(Task newParent) {
        newParent.addEntry(this);
        this.parent = newParent;
        return this;
    }

    /**
     * Sets the open arrival process for this entry.
     *
     * @param arrival The arrival distribution (e.g., Exp, Erlang, HyperExp)
     * @return this Entry for method chaining
     */
    public Entry setArrival(Distribution arrival) {
        this.arrival = arrival;
        return this;
    }

    /**
     * Gets the open arrival distribution for this entry.
     *
     * @return The arrival distribution, or null if not set
     */
    public Distribution getArrival() {
        return this.arrival;
    }

    /**
     * Add a forwarding call to another entry with a specified probability.
     * Forwarding allows this entry to redirect the reply to another entry
     * instead of replying directly to the original caller.
     *
     * @param dest The destination entry to forward to
     * @param prob The probability of forwarding (0.0 to 1.0)
     * @return This entry for method chaining
     */
    public Entry forward(Entry dest, double prob) {
        validateForwardingProbability(prob);
        int idx = this.forwardingDests.size();
        this.forwardingDests.put(idx, dest.getName());
        if (this.forwardingProbs.isEmpty()) {
            this.forwardingProbs = new Matrix(1, 1, 1);
            this.forwardingProbs.set(0, 0, prob);
        } else {
            this.forwardingProbs = this.forwardingProbs.concatCols(Matrix.singleton(prob));
        }
        return this;
    }

    /**
     * Add a forwarding call to another entry by name with a specified probability.
     *
     * @param destName The name of the destination entry
     * @param prob The probability of forwarding (0.0 to 1.0)
     * @return This entry for method chaining
     */
    public Entry forward(String destName, double prob) {
        validateForwardingProbability(prob);
        int idx = this.forwardingDests.size();
        this.forwardingDests.put(idx, destName);
        if (this.forwardingProbs.isEmpty()) {
            this.forwardingProbs = new Matrix(1, 1, 1);
            this.forwardingProbs.set(0, 0, prob);
        } else {
            this.forwardingProbs = this.forwardingProbs.concatCols(Matrix.singleton(prob));
        }
        return this;
    }

    /**
     * Add a forwarding call to another entry with probability 1.0.
     *
     * @param dest The destination entry to forward to
     * @return This entry for method chaining
     */
    public Entry forward(Entry dest) {
        return forward(dest, 1.0);
    }

    /**
     * Add a forwarding call to another entry by name with probability 1.0.
     *
     * @param destName The name of the destination entry
     * @return This entry for method chaining
     */
    public Entry forward(String destName) {
        return forward(destName, 1.0);
    }

    /**
     * Get the map of forwarding destinations.
     *
     * @return Map of forwarding destination entry names
     */
    public Map<Integer, String> getForwardingDests() {
        return this.forwardingDests;
    }

    /**
     * Get the matrix of forwarding probabilities.
     *
     * @return Matrix of forwarding probabilities
     */
    public Matrix getForwardingProbs() {
        return this.forwardingProbs;
    }

    public Map<Integer, String> getBoundToActivity() {
        return this.boundToActivity;
    }

    public Map<Integer, String> getReplyActivity() {
        return this.replyActivity;
    }

    public Task getParent() {
        return this.parent;
    }

    /**
     * Validate that a forwarding probability is valid and that the sum of all
     * forwarding probabilities does not exceed 1.0.
     *
     * @param prob The probability to validate
     * @throws IllegalArgumentException if probability is invalid or sum exceeds 1.0
     */
    protected void validateForwardingProbability(double prob) {
        if (prob < 0.0 || prob > 1.0) {
            throw new IllegalArgumentException("Forwarding probability must be between 0.0 and 1.0, got: " + prob);
        }
        double currentSum = this.forwardingProbs.isEmpty() ? 0.0 : this.forwardingProbs.sumRows(0);
        if (currentSum + prob > 1.0 + 1e-6) { // Small tolerance for floating point errors
            throw new IllegalArgumentException(
                String.format("Sum of forwarding probabilities would exceed 1.0 (current: %.6f, adding: %.6f)",
                             currentSum, prob));
        }
    }

    /**
     * Sets the type of this entry.
     *
     * @param type The entry type (e.g., "PH1PH2", "NONE")
     * @return this Entry for method chaining
     */
    public Entry setType(String type) {
        this.type = type;
        return this;
    }

    /**
     * Gets the type of this entry.
     *
     * @return The entry type
     */
    public String getType() {
        return this.type;
    }
}

