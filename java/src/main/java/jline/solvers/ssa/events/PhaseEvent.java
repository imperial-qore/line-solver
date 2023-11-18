package jline.solvers.ssa.events;

import jline.solvers.SolverOptions;

public class PhaseEvent extends Event{
    public PhaseEvent() {
        super();
    }

    public long getNPhases() {
        return 1;
    }
}
