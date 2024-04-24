package jline.util;

import java.util.List;
import java.util.Map;

public class UniqueRowResult {

    public final Matrix sortedMatrix;
    public final Matrix vi;
    public final Map<Integer, List<Integer>> vj;


    public UniqueRowResult(Matrix sortedMatrix, Matrix vi, Map<Integer, List<Integer>> vj) {
        this.sortedMatrix = sortedMatrix;
        this.vi = vi;
        this.vj = vj;
    }
}
