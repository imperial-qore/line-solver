/**
 * @file Markovian Arrival MAP with Marked arrivals coefficient computation
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import jline.util.matrix.Matrix;

public final class Mamap2m_coefficients {
    private Mamap2m_coefficients() {}

    /**
     * Result tuple holding G, U, Y matrices.
     */
    public static final class Mamap2mCoefficients {
        public final Matrix first;
        public final Matrix second;
        public final Matrix third;

        public Mamap2mCoefficients(Matrix first, Matrix second, Matrix third) {
            this.first = first;
            this.second = second;
            this.third = third;
        }
    }

    /**
     * Returns the coefficients used in the direct and inverse formulas for
     * fitting a MAMAP(2,m) in first canonical form (gamma > 0).
     */
    public static Mamap2mCoefficients mamap2m_can1_coefficients(double h1, double h2, double r1, double r2) {
        Matrix G = Matrix.zeros(15, 1);
        Matrix U = Matrix.zeros(12, 1);
        Matrix Y = Matrix.zeros(3, 1);

        G.set(0, 0, 1.0 - r1 / (r2 * (r1 - 1.0) + 1.0));
        G.set(1, 0, -(r1 * (r2 - 1.0)) / (r1 * r2 - r2 + 1.0));
        G.set(2, 0, (r1 * r2) / (r1 * r2 - r2 + 1.0));
        G.set(3, 0, (r1 * (r1 - 1.0)) / (r2 * (r1 - 1.0) + 1.0) - r1 + 1.0);
        G.set(4, 0, -(r1 * (r1 - 1.0) * (r2 - 1.0) * (r2 - 2.0)) / (r1 * r2 - r2 + 1.0));
        G.set(5, 0, (r1 * r2 * (r1 - 1.0) * (r2 - 1.0)) / (r1 * r2 - r2 + 1.0));
        G.set(6, 0, (r1 * r1 * (r2 - 1.0) * (r2 - 1.0)) / (r2 * (r1 - 1.0) + 1.0));
        G.set(7, 0, -(r1 * r2 * (r1 + 1.0) * (r2 - 1.0)) / (r1 * r2 - r2 + 1.0));
        G.set(8, 0, (r1 * r2 * r2) / (r1 * r2 - r2 + 1.0));
        G.set(9, 0, h1 - (h1 * r1) / (r2 * (r1 - 1.0) + 1.0));
        G.set(10, 0, -(r1 * (r2 - 1.0) * (h1 + h2 - h1 * r2)) / (r1 * r2 - r2 + 1.0));
        G.set(11, 0, (r1 * r2 * (h1 + h2 - h1 * r2)) / (r1 * r2 - r2 + 1.0));
        G.set(12, 0, ((h1 + h2 * r1) * (r1 - 1.0) * (r2 - 1.0)) / (r1 * r2 - r2 + 1.0));
        G.set(13, 0, -(r1 * (h1 + h2 * r1) * (r2 - 1.0)) / (r1 * r2 - r2 + 1.0));
        G.set(14, 0, (h2 * r1 * r2) / (r1 * r2 - r2 + 1.0));

        double a = (r1 * r2 - r2 + 1.0);
        U.set(0, 0, a * a);
        U.set(1, 0, -(r1 * r2 - r2 + 1.0) * (2.0 * h1 - h1 * r1 - 2.0 * h1 * r2 + 3.0 * h2 * r1 - h2 * r1 * r1 + h2 * r1 * r1 * r2 + h1 * r1 * r2 - h2 * r1 * r2));
        double t1 = (h1 - h2 + h2 * r1);
        U.set(2, 0, r1 * (r2 - 1.0) * t1 * t1);
        U.set(3, 0, (r1 * r2 - r2 + 1.0) * (h2 * h2 * r1 - h1 * h1 * r2 + h1 * h1 + h1 * h2 * r1 - h1 * h2 * r1 * r2));
        U.set(4, 0, -r1 * (r2 - 1.0) * (r1 * r2 - r2 + 1.0) * (h1 - h2 + h2 * r1));
        U.set(5, 0, r1 * (r2 - 1.0) * (h1 - h1 * r2 + h2 * r1) * (h1 - h2 + h2 * r1));
        U.set(6, 0, a * a);
        U.set(7, 0, -(r1 * r2 - r2 + 1.0) * (2.0 * h1 - 2.0 * h1 * r2 + h2 * r1 - h1 * r1 * r2 * r2 + h1 * r1 * r2 + h2 * r1 * r2));
        double t2 = (h2 - h1 * r2);
        U.set(8, 0, r1 * t2 * t2 * (r2 - 1.0));
        U.set(9, 0, (r1 * r2 - r2 + 1.0) * (h2 * h2 * r1 - h1 * h1 * r2 + h1 * h1 + h1 * h2 * r1 - h1 * h2 * r1 * r2));
        U.set(10, 0, -r1 * (h2 - h1 * r2) * (r2 - 1.0) * (r1 * r2 - r2 + 1.0));
        U.set(11, 0, r1 * (h2 - h1 * r2) * (r2 - 1.0) * (h1 - h1 * r2 + h2 * r1));

        Y.set(0, 0, G.get(0, 0) * G.get(10, 0) * G.get(14, 0) - G.get(0, 0) * G.get(11, 0) * G.get(13, 0)
                - G.get(1, 0) * G.get(9, 0) * G.get(14, 0) + G.get(1, 0) * G.get(11, 0) * G.get(12, 0)
                + G.get(2, 0) * G.get(9, 0) * G.get(13, 0) - G.get(2, 0) * G.get(10, 0) * G.get(12, 0));
        Y.set(1, 0, G.get(2, 0) * G.get(12, 0) - G.get(0, 0) * G.get(14, 0));
        Y.set(2, 0, G.get(9, 0) * G.get(2, 0) - G.get(11, 0) * G.get(0, 0));

        return new Mamap2mCoefficients(G, U, Y);
    }

    /**
     * Returns the coefficients used in the direct and inverse formulas for
     * fitting a MAMAP(2,m) in second canonical form (gamma < 0).
     */
    public static Mamap2mCoefficients mamap2m_can2_coefficients(double h1, double h2, double r1, double r2) {
        Matrix E = Matrix.zeros(14, 1);
        Matrix V = Matrix.zeros(12, 1);
        Matrix Z = Matrix.zeros(3, 1);

        E.set(0, 0, 1.0 - 1.0 / (r2 * (r1 - 1.0) - r1 + 2.0));
        E.set(1, 0, -(r2 - 1.0) / (r1 * (r2 - 1.0) - r2 + 2.0));
        E.set(2, 0, r2 / (r1 * (r2 - 1.0) - r2 + 2.0));
        E.set(3, 0, (r2 - 2.0) / (r1 * (r2 - 1.0) - r2 + 2.0) - r2 + 2.0);
        E.set(4, 0, r2 - r2 / (r1 * (r2 - 1.0) - r2 + 2.0));
        E.set(5, 0, -(r1 * (r2 - 1.0) * (r2 - 1.0)) / (r1 + r2 - r1 * r2 - 2.0));
        E.set(6, 0, -r2 - (r2 * (2.0 * r2 - 3.0)) / (r1 * (r2 - 1.0) - r2 + 2.0));
        E.set(7, 0, r2 * r2 / (r2 * (r1 - 1.0) - r1 + 2.0));
        E.set(8, 0, h1 - h1 / (r2 * (r1 - 1.0) - r1 + 2.0));
        E.set(9, 0, h1 * (r2 - 1.0) - ((r2 - 1.0) * (2.0 * h1 + h2 - h1 * r2)) / (r1 * (r2 - 1.0) - r2 + 2.0));
        E.set(10, 0, (r2 * (2.0 * h1 + h2 - h1 * r2)) / (r1 * (r2 - 1.0) - r2 + 2.0) - h1 * r2);
        E.set(11, 0, h2 - h2 / (r2 * (r1 - 1.0) - r1 + 2.0));
        E.set(12, 0, ((h1 + h2 * r1) * (r2 - 1.0)) / (r1 + r2 - r1 * r2 - 2.0));
        E.set(13, 0, (h2 * r2) / (r1 * (r2 - 1.0) - r2 + 2.0));

        double aa = (r1 + r2 - r1 * r2 - 2.0);
        V.set(0, 0, -aa * aa);
        V.set(1, 0, -(r1 + r2 - r1 * r2 - 2.0) * (2.0 * h1 + 2.0 * h2 - h1 * r2 - h2 * r2 + h2 * r1 * r2));
        V.set(2, 0, h2 * (2.0 * h1 - h1 * r2 + h2 * r1) * (r1 + r2 - r1 * r2 - 2.0));
        double t1 = (h1 - h2 + h2 * r1);
        V.set(3, 0, (r2 - 1.0) * t1 * t1);
        V.set(4, 0, (h1 - h2 + h2 * r1) * (2.0 * r2 - r1 * r2 + r1 * r2 * r2 - r2 * r2));
        V.set(5, 0, -(h1 * r2 + h2 * r2 - h1 * r2 * r2) * (h1 - h2 + h2 * r1));
        V.set(6, 0, -aa * aa);
        V.set(7, 0, -(r1 + r2 - r1 * r2 - 2.0) * (2.0 * h1 + 2.0 * h2 - h1 * r2 - h2 * r2 + h1 * r1 * r2 * r2 - h1 * r1 * r2));
        V.set(8, 0, h1 * (r1 + r2 - r1 * r2 - 2.0) * (2.0 * h2 + h1 * r1 - h2 * r2 + h1 * r1 * r2 * r2 - 2.0 * h1 * r1 * r2));
        double t2 = (h1 - h2 - h1 * r1 + h1 * r1 * r2);
        V.set(9, 0, (r2 - 1.0) * t2 * t2);
        V.set(10, 0, -r2 * (h1 - h2 - h1 * r1 + h1 * r1 * r2) * (r1 + r2 - r1 * r2 - 2.0));
        V.set(11, 0, -r2 * (h1 + h2 - h1 * r2) * (h1 - h2 - h1 * r1 + h1 * r1 * r2));

        Z.set(0, 0, E.get(9, 0) * E.get(11, 0) * E.get(2, 0) - E.get(9, 0) * E.get(13, 0) * E.get(0, 0)
                - E.get(10, 0) * E.get(11, 0) * E.get(1, 0) + E.get(10, 0) * E.get(12, 0) * E.get(0, 0)
                - E.get(12, 0) * E.get(2, 0) * E.get(8, 0) + E.get(13, 0) * E.get(1, 0) * E.get(8, 0));
        Z.set(1, 0, E.get(11, 0) * E.get(1, 0) - E.get(12, 0) * E.get(0, 0));
        Z.set(2, 0, E.get(9, 0) * E.get(0, 0) - E.get(1, 0) * E.get(8, 0));

        return new Mamap2mCoefficients(E, V, Z);
    }
}
