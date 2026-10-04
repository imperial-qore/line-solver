package jline.opt.variables;

import jline.lang.Network;
import jline.lang.layered.LayeredNetwork;

/**
 * Abstract base for LayeredNetwork (LQN) decision variables. Flat decision
 * variables mutate a {@link Network}; LQN variables mutate a
 * {@link LayeredNetwork}, so the flat {@link #apply(Network, Object)} is invalid
 * here and throws, while {@link #applyLqn(LayeredNetwork, Object)} carries the
 * LQN mutation. The evaluator branches on {@code problem.isLayered()} and calls
 * {@code applyLqn} for LQN variables. Mirrors the LQN variable classes in
 * native-Python {@code line_solver.opt.variables}.
 */
public abstract class LqnDecisionVariable extends DecisionVariable {

    protected LqnDecisionVariable(String name) {
        super(name);
    }

    /**
     * Invalid for LQN variables: they apply to a LayeredNetwork, not a flat
     * Network. The evaluator's layered path calls {@link #applyLqn} instead.
     */
    @Override
    public final void apply(Network model, Object value) {
        throw new UnsupportedOperationException(
                "LQN decision variable '" + name + "' cannot be applied to a flat Network; "
                        + "use applyLqn(LayeredNetwork, value).");
    }

    /** Apply a decoded value to a (copied) LayeredNetwork model. */
    public abstract void applyLqn(LayeredNetwork model, Object value);
}
