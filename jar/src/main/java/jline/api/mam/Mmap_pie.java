/**
 * @file Marked Markovian Arrival Process class-wise steady-state analysis
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

/**
 * Computes the steady-state probability vector for each class in an MMAP.
 */
public final class Mmap_pie {
    private Mmap_pie() {}

    public static Matrix mmap_pie(MatrixCell mmap) {
        if (mmap.size() < 2) {
            throw new IllegalArgumentException("MMAP must have at least D0 and D1");
        }

        int m = mmap.size() - 2;

        if (m <= 0) {
            Matrix result = new Matrix(1, 1);
            result.set(0, 0, 1.0);
            return result;
        }

        MatrixCell mapCell = new MatrixCell(2);
        mapCell.set(0, mmap.get(0));
        mapCell.set(1, mmap.get(1));
        Matrix piMap = Map_pie.map_pie(mapCell);

        Matrix classProbabilities = new Matrix(1, m);
        double totalRate = 0.0;
        double[] classRates = new double[m];

        for (int i = 0; i < m; i++) {
            Matrix Di = mmap.get(2 + i);
            double classRate = 0.0;
            for (int j = 0; j < Di.getNumRows(); j++) {
                for (int k = 0; k < Di.getNumCols(); k++) {
                    classRate += piMap.get(j) * Di.get(j, k);
                }
            }
            classRates[i] = classRate;
            totalRate += classRate;
        }

        if (totalRate > 0) {
            for (int i = 0; i < m; i++) {
                classProbabilities.set(0, i, classRates[i] / totalRate);
            }
        } else {
            for (int i = 0; i < m; i++) {
                classProbabilities.set(0, i, 1.0 / m);
            }
        }
        return classProbabilities;
    }

    public static Matrix mmap_pie(Matrix[] mmap) {
        return mmap_pie(new MatrixCell(mmap));
    }
}
