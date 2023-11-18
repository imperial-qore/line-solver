package jline.solvers.mva;

import jline.lang.FeatureSet;
import jline.lang.Network;
import jline.lang.NetworkStruct;
import jline.lang.constant.SolverType;
import jline.solvers.NetworkSolver;
import jline.solvers.SolverOptions;

public class SolverMVA extends NetworkSolver {
	
	public SolverMVA(Network model, SolverOptions options) {
		super(model, "SolverMVA", options);
		this.sn = model.getStruct(false);
		this.result = new SolverMVAResult();
	}
	
	public SolverMVA(Network model) {
		super(model, "SolverMVA", SolverMVA.defaultOptions());
		this.sn = model.getStruct(false);
		this.result = new SolverMVAResult();
	}
	
	public NetworkStruct getStruct() {
		if (this.sn == null) 
			this.sn = this.model.getStruct(false);
		return this.sn;
	}
	
	public void setStruct(NetworkStruct sn) {
		this.sn = sn;
	}

	@Override
	public void runAnalyzer() throws IllegalAccessException {
		if (this.model == null)
			throw new RuntimeException("Model is not provided");
		if (this.sn == null) 
			this.sn = this.model.getStruct(false);
		if (this.options == null)
			this.options = new SolverOptions(SolverType.MVA);
		
		MVARunner runner = new MVARunner(this);
		this.result = runner.run();
	}

	/**
	 * Returns the feature set supported by the MVA solver
	 * @return - the feature set supported by the MVA solver
	 */
	public static FeatureSet getFeatureSet(){
		FeatureSet s = new FeatureSet();
		// TODO: update with the features supported by JLINE. These are the features supported by LINE.
		String[] features = {"Sink","Source",
		"ClassSwitch","Delay","Queue",
		"APH","Coxian","Erlang","Exp","HyperExp",
		"Pareto","Weibull","Lognormal","Uniform","Det",
		"StatelessClassSwitcher","InfiniteServer","SharedServer","Buffer","Dispatcher",
		"CacheClassSwitcher","Cache",
		"Server","RandomSource","ServiceTunnel",
		"SchedStrategy_INF","SchedStrategy_PS",
		"SchedStrategy_DPS","SchedStrategy_FCFS","SchedStrategy_SIRO","SchedStrategy_HOL",
		"SchedStrategy_LCFSPR",
		"Fork","Forker","Join","Joiner",
		"RoutingStrategy_PROB","RoutingStrategy_RAND",
		"ClosedClass","OpenClass","Replayer"};
		s.setTrue(features);
		return s;
	}

	/**
	 * Checks whether the given model is supported by the MVA solver
	 * @param model - the network model
	 * @return - true if the model is supported, false otherwise
	 */
	public static boolean supports(Network model){
		FeatureSet featUsed = model.getUsedLangFeatures();
		FeatureSet featSupported = SolverMVA.getFeatureSet();
		return FeatureSet.supports(featSupported, featUsed);
	}

	public static SolverOptions defaultOptions() {
		return new SolverOptions(SolverType.MVA);
	}
}
