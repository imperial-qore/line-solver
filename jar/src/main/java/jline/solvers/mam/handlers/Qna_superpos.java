package jline.solvers.mam.handlers;

import jline.api.da.Da_traffic_superpos;
import jline.util.matrix.Matrix;

/** Legacy entry point; the algorithm lives in {@link Da_traffic_superpos}. */
public final class Qna_superpos {
    private Qna_superpos() {}

    public static double qna_superpos(Matrix lambda, Matrix a2) {
        return Da_traffic_superpos.da_traffic_superpos(lambda, a2);
    }
}
