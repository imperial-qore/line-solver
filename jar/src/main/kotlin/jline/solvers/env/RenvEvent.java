package jline.solvers.env;

import jline.lang.Event;
import jline.lang.constant.EventType;
import jline.util.Pair;
import jline.util.matrix.Matrix;

public class RenvEvent extends Event {
    public Event event;
    private final Pair<Integer, Integer> envTransition;

    public RenvEvent(int nodeIdx, int jobclassIdx, double prob, Matrix state, double t, double job, Pair<Integer, Integer> envTransition) {
        super(EventType.STAGE, nodeIdx, jobclassIdx, prob, state, t, job);
        this.envTransition = envTransition;

    }
}
