package jline.solvers.ag.handlers;

import java.util.List;

import jline.util.matrix.Matrix;

/**
 * RCAT (Reversed Compound Agent Theorem) model representation.
 *
 * Every component is a QBD: NLEV levels of MPH phases each, the level being the
 * queue length and the phase the pair (arrival phase, service phase). LEVEL
 * gives the level index of every state, SVCRATE the service completion rate out
 * of every state (zero on level 0, where no server runs) and SVCDOWN the same
 * rate per phase at a busy level. With exponential processes MPH is 1 and the
 * QBD is the scalar birth-death chain this analyzer built before.
 */
public final class RCATModel {
    public final Matrix[][] R;
    public final Matrix AP;
    public final Matrix processMap;
    public final List<ActionInfo> actionMap;
    public final int[] N;
    public final int[] nlev;
    public final int[] mph;
    public final int[][] level;
    public final double[][] svcrate;
    public final double[][] svcdown;

    public RCATModel(Matrix[][] R, Matrix AP, Matrix processMap, List<ActionInfo> actionMap, int[] N,
                     int[] nlev, int[] mph, int[][] level, double[][] svcrate, double[][] svcdown) {
        this.R = R;
        this.AP = AP;
        this.processMap = processMap;
        this.actionMap = actionMap;
        this.N = N;
        this.nlev = nlev;
        this.mph = mph;
        this.level = level;
        this.svcrate = svcrate;
        this.svcdown = svcdown;
    }
}
