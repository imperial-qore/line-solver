package jline.lang;

import java.io.Serializable;
import java.util.Map;

import jline.lang.constant.EventType;
import jline.lang.nodes.Node;
import jline.util.SerializableFunction;
import jline.util.Matrix;
import jline.util.Pair;

/**
 * Class abstracting an event within a Network model
 */
public class NetworkEvent implements Serializable {
	
	protected int nodeIdx;
	protected EventType event;
	protected int jobclassIdx;
	protected double prob;
	protected SerializableFunction<Pair<Map<Node, Matrix>, Map<Node, Matrix>>, Double> probFun;
	protected Matrix state;
	protected double t;
	protected double job;
	
	/*
	 * prob = NaN if not set
	 * t = NaN if not set
	 * job = NaN if not set
	 * state = new JLineMatrix(0,0) if not set
	 */
	public NetworkEvent(EventType event, int nodeIdx, int jobclassIdx, double prob, Matrix state, double t, double job) {
		this.event = event;
		this.nodeIdx = nodeIdx;
		this.jobclassIdx = jobclassIdx;
		this.prob = prob;
		this.state = state;
		this.t = t;
		this.job = job;
		this.probFun = null;
	}
	
	/*
	 * The input probability might be a function
	 */
	public NetworkEvent(EventType event, int nodeIdx, int jobclassIdx, SerializableFunction<Pair<Map<Node, Matrix>, Map<Node, Matrix>>, Double> probFun, Matrix state, double t, double job) {
		this.event = event;
		this.nodeIdx = nodeIdx;
		this.jobclassIdx = jobclassIdx;
		this.prob = Double.NaN;
		this.state = state;
		this.t = t;
		this.job = job;
		this.probFun = probFun;
	}
	
	public int getNodeIdx() {
		return nodeIdx;
	}
	
	public void setNodeIdx(int nodeIdx) {
		this.nodeIdx = nodeIdx;
	}

	public EventType getEvent() {
		return event;
	}

	public void setEvent(EventType event) {
		this.event = event;
	}

	public int getJobclassIdx() {
		return jobclassIdx;
	}

	public void setJobclassIdx(int jobclassIdx) {
		this.jobclassIdx = jobclassIdx;
	}

	public double getProb() {
		return prob;
	}
	
	public double getProb(Pair<Map<Node, Matrix>, Map<Node, Matrix>> state) {
		return this.probFun.apply(state);
	}

	public void setProb(double prob) {
		this.prob = prob;
	}

	public SerializableFunction<Pair<Map<Node, Matrix>, Map<Node, Matrix>>, Double> getProbFun() {
		return probFun;
	}

	public void setProbFun(SerializableFunction<Pair<Map<Node, Matrix>, Map<Node, Matrix>>, Double> probFun) {
		this.probFun = probFun;
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
