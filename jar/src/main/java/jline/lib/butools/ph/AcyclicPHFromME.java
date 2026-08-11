/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 * Ported from BUTools V2.0
 */
package jline.lib.butools.ph;

import jline.lib.butools.reptrans.ExtendToMarkovian;
import jline.lib.butools.reptrans.MarkovianRepresentation;
import jline.lib.butools.reptrans.SimilarityMatrix;
import jline.lib.butools.reptrans.TransformToAcyclic;
import jline.util.matrix.Matrix;

public final class AcyclicPHFromME {
    private AcyclicPHFromME() {}

    /**
     * Transforms an arbitrary matrix-exponential representation
     * to an acyclic phase-type representation.
     *
     * @param alpha Initial vector of the distribution (shape 1 x N)
     * @param A Matrix parameter of the distribution (shape N x N)
     * @param maxSize The maximum number of phases for the result. The default value is 100.
     * @param prec Vector and matrix entries smaller than the precision are considered to be zeros.
     * @return PHRepresentation containing the initial probability vector and
     *         transient generator matrix of the Markovian acyclic representation.
     * @throws IllegalArgumentException if no Markovian acyclic representation has been found.
     */
    public static PHRepresentation acyclicPHFromME(Matrix alpha, Matrix A, int maxSize, double prec) {
        Matrix G = TransformToAcyclic.transformToAcyclic(A, maxSize, prec);

        Matrix T = SimilarityMatrix.similarityMatrix(A, G);
        Matrix gamma = alpha.mult(T);

        boolean hasNeg = false;
        for (int i = 0; i < gamma.length(); i++) {
            if (gamma.get(0, i) < -prec) {
                hasNeg = true;
                break;
            }
        }

        if (hasNeg) {
            MarkovianRepresentation result = ExtendToMarkovian.extendToMarkovian(gamma, G, maxSize, prec);
            if (!CheckPHRepresentation.checkPHRepresentation(result.getBeta(), result.getB(), prec)) {
                throw new IllegalArgumentException("AcyclicPHFromME: No acyclic representation found up to the given size and precision!");
            }
            return new PHRepresentation(result.getBeta(), result.getB());
        } else {
            if (!CheckPHRepresentation.checkPHRepresentation(gamma, G, prec)) {
                throw new IllegalArgumentException("AcyclicPHFromME: Result is not a valid PH representation!");
            }
            return new PHRepresentation(gamma, G);
        }
    }

    public static PHRepresentation acyclicPHFromME(Matrix alpha, Matrix A, int maxSize) {
        return acyclicPHFromME(alpha, A, maxSize, 1e-14);
    }

    public static PHRepresentation acyclicPHFromME(Matrix alpha, Matrix A) {
        return acyclicPHFromME(alpha, A, 100, 1e-14);
    }
}
