package jline.solvers.ag.handlers;

import jline.lang.processes.DiscreteDistribution;

/** Information about an action in an RCAT model. */
public final class ActionInfo {
    public final int fromStation;
    public final int fromClass;
    public final int toStation;
    public final int toClass;
    public final double prob;
    public final boolean isNegative;
    public final boolean isCatastrophe;
    public final DiscreteDistribution removalDistribution;

    public ActionInfo(int fromStation, int fromClass, int toStation, int toClass, double prob) {
        this(fromStation, fromClass, toStation, toClass, prob, false, false, null);
    }

    public ActionInfo(int fromStation, int fromClass, int toStation, int toClass, double prob, boolean isNegative) {
        this(fromStation, fromClass, toStation, toClass, prob, isNegative, false, null);
    }

    public ActionInfo(int fromStation, int fromClass, int toStation, int toClass, double prob,
                      boolean isNegative, boolean isCatastrophe, DiscreteDistribution removalDistribution) {
        this.fromStation = fromStation;
        this.fromClass = fromClass;
        this.toStation = toStation;
        this.toClass = toClass;
        this.prob = prob;
        this.isNegative = isNegative;
        this.isCatastrophe = isCatastrophe;
        this.removalDistribution = removalDistribution;
    }
}
