/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang;

import jline.lang.constant.EventType;
import jline.lang.nodes.Node;
import jline.util.Pair;
import jline.util.SerializableFunction;
import jline.util.matrix.Matrix;

import java.io.Serializable;
import java.util.Map;

/**
 * Class abstracting an event within a Network model
 */
public class Event implements Serializable {

    protected int node;
    protected EventType event;
    protected int jobclass;
    protected double prob;
    protected SerializableFunction<Pair<Map<Node, Matrix>, Map<Node, Matrix>>, Double> probFun;
    protected Matrix state;
    protected double t;
    protected double job;

    /*
     * prob = NaN if not set
     * t = NaN if not set
     * job = NaN if not set
     * state = new Matrix(0,0) if not set
     */
    /**
     * Creates a new Event with basic parameters. Other properties are initialized with default values:
     * probability = NaN, state = empty matrix, time = NaN, job = NaN, probFun = null.
     * 
     * @param event the type of event that occurred
     * @param node the node index where the event occurred
     * @param jobclass the job class index associated with the event
     */
    public Event(EventType event, int node, int jobclass) {
        this.event = event;
        this.node = node;
        this.jobclass = jobclass;
        this.prob = Double.NaN;
        this.state = new Matrix(0,0);
        this.t = Double.NaN;
        this.job = Double.NaN;
        this.probFun = null;
    }


    /*
     * prob = NaN if not set
     * t = NaN if not set
     * job = NaN if not set
     * state = new JLineMatrix(0,0) if not set
     */
    /**
     * Creates a new Event with all basic parameters specified.
     * 
     * @param event the type of event that occurred
     * @param node the node index where the event occurred  
     * @param jobclass the job class index associated with the event
     * @param prob the probability associated with this event
     * @param state the system state matrix when this event occurred
     * @param t the time when this event occurred
     * @param job the job identifier associated with this event
     */
    public Event(EventType event, int node, int jobclass, double prob, Matrix state, double t, double job) {
        this.event = event;
        this.node = node;
        this.jobclass = jobclass;
        this.prob = prob;
        this.state = state;
        this.t = t;
        this.job = job;
        this.probFun = null;
    }

    /*
     * The input probability might be a function
     */
    /**
     * Creates a new Event where probability is determined by a function of system state.
     * The probability parameter is set to NaN since it will be computed dynamically.
     * 
     * @param event the type of event that occurred
     * @param node the node index where the event occurred
     * @param jobclass the job class index associated with the event
     * @param probFun function that computes probability based on current system state
     * @param state the system state matrix when this event occurred
     * @param t the time when this event occurred
     * @param job the job identifier associated with this event
     */
    public Event(EventType event, int node, int jobclass, SerializableFunction<Pair<Map<Node, Matrix>, Map<Node, Matrix>>, Double> probFun, Matrix state, double t, double job) {
        this.event = event;
        this.node = node;
        this.jobclass = jobclass;
        this.prob = Double.NaN;
        this.state = state;
        this.t = t;
        this.job = job;
        this.probFun = probFun;
    }

    /**
     * Gets the type of this event.
     * 
     * @return the event type
     */
    public EventType getEvent() {
        return event;
    }

    /**
     * Sets the type of this event.
     * 
     * @param event the new event type
     */
    public void setEvent(EventType event) {
        this.event = event;
    }

    /**
     * Gets the job identifier associated with this event.
     * 
     * @return the job identifier, or NaN if not set
     */
    public double getJob() {
        return job;
    }

    /**
     * Sets the job identifier for this event.
     * 
     * @param job the job identifier
     */
    public void setJob(double job) {
        this.job = job;
    }

    /**
     * Gets the job class index associated with this event.
     * 
     * @return the job class index
     */
    public int getJobClass() {
        return jobclass;
    }

    /**
     * Sets the job class index for this event.
     * 
     * @param jobclass the job class index
     */
    public void setJobClass(int jobclass) {
        this.jobclass = jobclass;
    }

    /**
     * Gets the node index where this event occurred.
     * 
     * @return the node index
     */
    public int getNode() {
        return node;
    }

    /**
     * Sets the node index where this event occurred.
     * 
     * @param node the node index
     */
    public void setNode(int node) {
        this.node = node;
    }

    /**
     * Computes the probability of this event using the probability function and current system state.
     * This method requires that a probability function was set during construction.
     * 
     * @param state the current system state as a pair of node-to-matrix mappings
     * @return the computed probability
     * @throws NullPointerException if no probability function was set
     */
    public double getProb(Pair<Map<Node, Matrix>, Map<Node, Matrix>> state) {
        return this.probFun.apply(state);
    }

    /**
     * Gets the fixed probability value for this event.
     * 
     * @return the probability, or NaN if not set or if using a probability function
     */
    public double getProb() {
        return prob;
    }

    /**
     * Sets the fixed probability value for this event.
     * 
     * @param prob the probability value
     */
    public void setProb(double prob) {
        this.prob = prob;
    }

    /**
     * Gets the probability function used to dynamically compute event probability.
     * 
     * @return the probability function, or null if using fixed probability
     */
    public SerializableFunction<Pair<Map<Node, Matrix>, Map<Node, Matrix>>, Double> getProbFun() {
        return probFun;
    }

    /**
     * Sets the probability function to dynamically compute event probability based on system state.
     * 
     * @param probFun function that takes system state and returns probability
     */
    public void setProbFun(SerializableFunction<Pair<Map<Node, Matrix>, Map<Node, Matrix>>, Double> probFun) {
        this.probFun = probFun;
    }

    /**
     * Gets the system state matrix when this event occurred.
     * 
     * @return the state matrix, or empty matrix if not set
     */
    public Matrix getState() {
        return state;
    }

    /**
     * Sets the system state matrix for this event.
     * 
     * @param state the system state matrix
     */
    public void setState(Matrix state) {
        this.state = state;
    }

    /**
     * Gets the time when this event occurred.
     * 
     * @return the event time, or NaN if not set
     */
    public double getT() {
        return t;
    }

    /**
     * Sets the time when this event occurred.
     * 
     * @param t the event time
     */
    public void setT(double t) {
        this.t = t;
    }
}
