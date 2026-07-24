/**
 * @file Absorbing Phase-type distribution assembly from parameters
 *
 * Constructs APH(2) transition matrices from specified rates and transition probabilities.
 * Used for building phase-type distributions from fitted or prescribed parameters.
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Aph2_assemble {
    private Aph2_assemble() {}

    /**
     * Assembles an acyclic phase-type (APH) distribution with two phases (APH(2)) using the given parameters.
     *
     * <p>This method constructs the transition matrices {@code D0} and {@code D1} for an APH(2) distribution
     * based on the parameters: {@code l1} and {@code l2} (the rates of the exponential phases), and {@code p1}
     * (the probability of transitioning between the phases). The resulting distribution has two phases, each
     * represented by an exponential distribution, with transitions between the phases as specified by the
     * given parameters.
     *
     * @param l1 the rate of the first exponential phase
     * @param l2 the rate of the second exponential phase
     * @param p1 the probability of transitioning from the first phase to the second phase
     * @return a MatrixCell containing the transition matrices {@code D0} and {@code D1} of the APH(2) distribution
     */
    public static MatrixCell aph2_assemble(double l1, double l2, double p1) {
        MatrixCell APH = new MatrixCell();
        Matrix D0 = new Matrix(2, 2, 4);
        D0.set(0, 0, -1 / l1);
        D0.set(0, 1, 1 / l1 * p1);
        D0.set(1, 0, 0);
        D0.set(1, 1, -1 / l2);

        Matrix D1 = new Matrix(2, 2, 4);
        D1.set(0, 0, 1 / l1 * (1 - p1));
        D1.set(0, 1, 0);
        D1.set(1, 0, 1 / l2);
        D1.set(1, 1, 0);

        APH.set(0, D0);
        APH.set(1, D1);

        return APH;
    }
}
