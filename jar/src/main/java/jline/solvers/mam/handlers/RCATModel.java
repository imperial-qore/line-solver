package jline.solvers.mam.handlers;

import java.util.List;

import jline.util.matrix.Matrix;

/** RCAT (Reversed Compound Agent Theorem) model representation. */
public final class RCATModel {
    public final Matrix[][] R;
    public final Matrix AP;
    public final Matrix processMap;
    public final List<ActionInfo> actionMap;
    public final int[] N;

    public RCATModel(Matrix[][] R, Matrix AP, Matrix processMap, List<ActionInfo> actionMap, int[] N) {
        this.R = R;
        this.AP = AP;
        this.processMap = processMap;
        this.actionMap = actionMap;
        this.N = N;
    }
}
