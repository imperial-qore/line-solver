package jline.api.trace;

import jline.util.matrix.Matrix;

public final class Mtrace_joint {
    private Mtrace_joint() {}

    /**
     * Given a multi-class trace, computes the empirical class-dependent joint
     * moments that estimate E[ ( X^(a)_j )^i(1) (X^(a)_(j+l) )^i(2) ]
     * for all classes a.
     *
     * @param T the inter-event times
     * @param A the class of each event
     * @param i a vector specifying the power of each variable in the joint moment
     * @return the joint moment of each class
     */
    public static Matrix mtrace_joint(double[] T, int[] A, int[] i) {
        // number of classes
        int C = 0;
        for (int v : A) {
            if (v > C) C = v;
        }

        // number of events
        int N = A.length;

        // result
        double[] JM = new double[C];

        for (int a = 1; a <= C; a++) {
            // Count events of class a, excluding the first and last event
            int Na = 0;
            for (int j = 1; j < N - 1; j++) {
                if (A[j] == a) Na++;
            }

            double tmp = 0.0;
            for (int j = 0; j < N - 2; j++) {
                if (A[j + 1] == a) {
                    tmp += Math.pow(T[j], i[0]) * Math.pow(T[j + 1], i[1]);
                }
            }

            JM[a - 1] = (Na > 0) ? tmp / Na : 0.0;
        }

        return new Matrix(JM);
    }
}
