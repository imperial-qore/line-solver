package jline.lang;

import java.io.Serializable;
import java.util.List;
import java.util.Map;

import jline.lang.constant.ProcessType;

import jline.lang.constant.TimingStrategy;
import jline.lang.constant.JoinStrategy;
import jline.lang.constant.ReplacementStrategy;
import jline.util.Matrix;

/**
 * Class for the nodeparam field within NetworkStruct
 */
public class NodeParam implements Serializable {
	//TODO: this data structure is a blob and may need reimplementation

	// Cache
	public Matrix[][] accost;
	public Matrix hitclass;
	public Matrix itemcap;
	public Matrix missclass;
	public int nitems;
		public Map<Integer, List<Double>> pread;
	public ReplacementStrategy rpolicy;

	//Fork
	public double fanOut = Double.NaN;
	
	//Join
	public Map<JobClass, JoinStrategy> joinStrategy = null;
	public Map<JobClass, Double> fanIn = null;
	public Map<JobClass, Double> joinRequired;

	//RoutingStrategy WRROBIN:
	public Map<JobClass, Matrix> weights = null;
	
	//RoutingStrategy RROBIN, WRROBIN:
	public Map<JobClass, Matrix> outlinks = null;

	//RoutingStrategy KCHOICES:
	public Map<JobClass, Matrix> withMemory = null;
	public Map<JobClass, Integer> k;

	public List<Matrix> enabling = null;

	public List<Matrix> inhibiting = null;

	public List<String> modenames;

	public Matrix nmodeservers;

	public Map<JobClass, TimingStrategy> firingid;

	public List<Matrix> firing = null;

	public Map<JobClass, ProcessType> firingprocid;

	public Map<JobClass, Map<Integer, Matrix>> firingproc;

	public List<Integer> firingphases;

	public Matrix firingprio;

	public Matrix fireweight;

	public String fileName = null;

	public String filePath = null;

	public String startTime = null;
	public String loggerName = null;

	public String timestamp = null;

	public String jobID = null;

	public String jobClass = null;

	public String timeSameClass = null;

	public String timeAnyClass = null;

	public int nmodes;

	public boolean isEmpty() {
		return Double.isNaN(fanOut) && joinStrategy == null && fanIn == null && weights == null && outlinks == null;
	}
}
