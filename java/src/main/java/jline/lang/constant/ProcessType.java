package jline.lang.constant;

import java.io.Serializable;

/**
 *  Constants for specifying a stochastic process
 */
public enum ProcessType implements Serializable {
	EXP,
	ERLANG,
	DISABLED,
	IMMEDIATE,
	HYPEREXP,
	APH,
	COXIAN,
	PH,
	MAP,
	UNIFORM,
	DET,
	GAMMA,
	PARETO,
	WEIBULL,
	LOGNORMAL,
	MMPP2,
	REPLAYER,
	TRACE,
	COX2,
	BINOMIAL,
	POISSON;

	public static String toText(ProcessType type){
		switch (type){
			case EXP:
				return "Exp";
			case ERLANG:
				return "Erlang";
			case HYPEREXP:
				return "HyperExp";
			case PH:
				return "PH";
			case APH:
				return "APH";
			case MAP:
				return "MAP";
			case UNIFORM:
				return "Uniform";
			case DET:
				return "Det";
			case COXIAN:
				return "Coxian";
			case GAMMA:
				return "Gamma";
			case PARETO:
				return "Pareto";
			case MMPP2:
				return "MMPP2";
			case REPLAYER:
			case TRACE:
				return "Replayer";
			case IMMEDIATE:
				return "Immediate";
			case DISABLED:
				return "Disabled";
			case COX2:
				return "Cox2";
			case WEIBULL:
				return "Weibull";
			case LOGNORMAL:
				return "Lognormal";
			case POISSON:
				return "Poisson";
			case BINOMIAL:
				return "Binomial";
			default:
				return type.name();
		}
	}
}