/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.constant;

import java.io.Serializable;

/**
 * Represents a balking threshold for queue-length based balking.
 *
 * A threshold specifies that when the queue length is within a certain range,
 * the customer will balk (refuse to join) with a given probability.
 *
 * Example: BalkingThreshold(5, 10, 0.3) means:
 * - When queue length is between 5 and 10 (inclusive)
 * - Customer balks with 30% probability
 *
 * Thresholds are evaluated in order; the first matching range determines the probability.
 */
public class BalkingThreshold implements Serializable {

    private final int minJobs;
    private final int maxJobs;
    private final double probability;

    /**
     * Creates a balking threshold for a queue length range.
     *
     * @param minJobs minimum queue length (inclusive) for this threshold
     * @param maxJobs maximum queue length (inclusive), use Integer.MAX_VALUE for unbounded
     * @param probability probability of balking when queue length is in range [0.0, 1.0]
     * @throws IllegalArgumentException if probability is not in [0, 1] or minJobs > maxJobs
     */
    public BalkingThreshold(int minJobs, int maxJobs, double probability) {
        if (probability < 0.0 || probability > 1.0) {
            throw new IllegalArgumentException("Balking probability must be in [0.0, 1.0], got: " + probability);
        }
        if (minJobs > maxJobs) {
            throw new IllegalArgumentException("minJobs (" + minJobs + ") cannot exceed maxJobs (" + maxJobs + ")");
        }
        if (minJobs < 0) {
            throw new IllegalArgumentException("minJobs cannot be negative: " + minJobs);
        }
        this.minJobs = minJobs;
        this.maxJobs = maxJobs;
        this.probability = probability;
    }

    /**
     * Creates a balking threshold for all queue lengths >= minJobs.
     *
     * @param minJobs minimum queue length (inclusive)
     * @param probability probability of balking
     */
    public BalkingThreshold(int minJobs, double probability) {
        this(minJobs, Integer.MAX_VALUE, probability);
    }

    /**
     * Get the minimum queue length for this threshold.
     * @return minimum queue length (inclusive)
     */
    public int getMinJobs() {
        return minJobs;
    }

    /**
     * Get the maximum queue length for this threshold.
     * @return maximum queue length (inclusive)
     */
    public int getMaxJobs() {
        return maxJobs;
    }

    /**
     * Get the balking probability for this threshold.
     * @return probability in [0.0, 1.0]
     */
    public double getProbability() {
        return probability;
    }

    /**
     * Check if a given queue length falls within this threshold's range.
     *
     * @param queueLength the current queue length
     * @return true if queueLength is in [minJobs, maxJobs]
     */
    public boolean matches(int queueLength) {
        return queueLength >= minJobs && queueLength <= maxJobs;
    }

    @Override
    public String toString() {
        if (maxJobs == Integer.MAX_VALUE) {
            return "BalkingThreshold[" + minJobs + "+, p=" + probability + "]";
        }
        return "BalkingThreshold[" + minJobs + "-" + maxJobs + ", p=" + probability + "]";
    }
}
