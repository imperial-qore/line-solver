package jline.solvers.ssa.strategies;

public enum TauLeapingStateStrategy {
    Cutoff,                // cutoff number between 0 and [capacity]
    TimeWarp,              // try again if it reaches an illegal state
    TauTimeWarp,           // Time warp with adjustments to tau value at every step (up if legal, down if illegal)
    TwoTimes,              // cycle through each event twice if it doesn't reach a legal state
}
