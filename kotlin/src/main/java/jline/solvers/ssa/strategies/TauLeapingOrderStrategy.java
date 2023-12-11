package jline.solvers.ssa.strategies;

public enum TauLeapingOrderStrategy {
    RandomEvent,
    RandomEventFixed,
    InOrder,        // update in the order that was given
    DirectedGraph,  // attempt to create a directed graph of the network, fail if not possible
                    // consider feedbacks in this, cause chaos
    DirectedCycle   // create a directed graph and cycle through start node
}
