/**
 * @file Marked Markovian Arrival Process random generation
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

import java.util.Random;

public final class Mmap_rand {
    private Mmap_rand() {}

    /**
     * Generates a random MMAP with a given order and number of classes.
     *
     * @param order   the number of phases (order) in the MAP
     * @param classes the number of different classes (types) of arrivals in the MMAP
     * @return a MatrixCell representing the MMAP
     */
    public static MatrixCell mmap_rand(int order, int classes) {
        MatrixCell MMAP = new MatrixCell();
        for (int c = 0; c < 2 + classes; c++) {
            MMAP.set(c, new Matrix(order, order));
        }
        MatrixCell MAP = Map_rand.map_rand(order);

        MMAP.set(0, MAP.get(0));
        MMAP.set(1, MAP.get(1));

        Random rand = new Random();
        for (int i = 0; i < order; i++) {
            double[] p = new double[classes];
            double sum = 0.0;

            for (int c = 0; c < classes; c++) {
                p[c] = rand.nextDouble();
                sum += p[c];
            }

            for (int c = 0; c < classes; c++) {
                p[c] /= sum;
            }

            for (int j = 0; j < order; j++) {
                for (int c = 0; c < classes; c++) {
                    MMAP.set(2 + c, MMAP.get(1).copy());
                    MMAP.get(2 + c).scaleEq(p[c]);
                }
            }
        }

        return MMAP;
    }
}
