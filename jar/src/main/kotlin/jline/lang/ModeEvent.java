/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang;

import jline.lang.constant.EventType;
import jline.util.matrix.Matrix;

import java.io.Serializable;

/**
 * A mode event occurring in a Network.
 */
public class ModeEvent implements Serializable {
    
    private int node;
    private EventType event;
    private int mode;
    private double weight;
    private double prob;
    private Matrix state; // state information when the event occurs (optional)
    private double t; // timestamp when the event occurs (optional)
    private double job; // job id (optional)
    
    public ModeEvent(EventType event, int node, int mode) {
        this(event, node, mode, 1.0);
    }
    
    public ModeEvent(EventType event, int node, int mode, double weight) {
        this(event, node, mode, weight, Double.NaN);
    }
    
    public ModeEvent(EventType event, int node, int mode, double weight, double prob) {
        this(event, node, mode, weight, prob, null);
    }
    
    public ModeEvent(EventType event, int node, int mode, double weight, double prob, Matrix state) {
        this(event, node, mode, weight, prob, state, Double.NaN);
    }
    
    public ModeEvent(EventType event, int node, int mode, double weight, double prob, Matrix state, double t) {
        this(event, node, mode, weight, prob, state, t, Double.NaN);
    }
    
    public ModeEvent(EventType event, int node, int mode, double weight, double prob, Matrix state, double t, double job) {
        this.node = node;
        this.event = event;
        this.mode = mode;
        this.weight = weight;
        this.prob = prob;
        this.state = state;
        this.t = t;
        this.job = job;
    }
    
    // Getters and setters
    public int getNode() {
        return node;
    }
    
    public void setNode(int node) {
        this.node = node;
    }
    
    public EventType getEvent() {
        return event;
    }
    
    public void setEvent(EventType event) {
        this.event = event;
    }
    
    public int getMode() {
        return mode;
    }
    
    public void setMode(int mode) {
        this.mode = mode;
    }
    
    public double getWeight() {
        return weight;
    }
    
    public void setWeight(double weight) {
        this.weight = weight;
    }
    
    public double getProb() {
        return prob;
    }
    
    public void setProb(double prob) {
        this.prob = prob;
    }
    
    public Matrix getState() {
        return state;
    }
    
    public void setState(Matrix state) {
        this.state = state;
    }
    
    public double getT() {
        return t;
    }
    
    public void setT(double t) {
        this.t = t;
    }
    
    public double getJob() {
        return job;
    }
    
    public void setJob(double job) {
        this.job = job;
    }
}