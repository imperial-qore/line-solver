package jline.solvers.ssa.strategies;

import jline.solvers.ssa.strategies.TauLeapingOrderStrategy;
import jline.solvers.ssa.strategies.TauLeapingStateStrategy;
import jline.solvers.ssa.strategies.TauLeapingVarType;

public class TauLeapingType {
    protected TauLeapingVarType varType;
    protected TauLeapingOrderStrategy orderStrategy;
    protected TauLeapingStateStrategy stateStrategy;
    protected double tau;
    protected double initialTau;

    public TauLeapingType(TauLeapingVarType varType,
                          TauLeapingOrderStrategy orderStrategy,
                          TauLeapingStateStrategy stateStrategy,
                          double tau) {
        this.varType = varType;
        this.orderStrategy = orderStrategy;
        this.stateStrategy = stateStrategy;
        this.tau = tau;
        this.initialTau = tau;
    }

    public TauLeapingVarType getVarType() {
        return this.varType;
    }

    public TauLeapingStateStrategy getStateStrategy() {
        return this.stateStrategy;
    }

    public TauLeapingOrderStrategy getOrderStrategy() {
        return this.orderStrategy;
    }

    public double getTau() {
        return this.tau;
    }

    public void setTau(double tau) {
        this.tau = tau;
    }

    public void resetTau() {
        this.tau = initialTau;
    }

    public String stringRep() {
        return String.format("%s_%s_%s_%f", varType.name(), orderStrategy.name(), stateStrategy.name(), tau);
    }
}
