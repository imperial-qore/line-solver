package jline.api.mam;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Mmap_sigma2 {
    private Mmap_sigma2() {}

    /**
     * Computes two-step class transition probabilities for an MMAP.
     *
     * @return a 3D array of class-transition probabilities
     */
    public static Double[][][] mmap_sigma2(MatrixCell mmap) {
        int C = mmap.size() - 2;

        Double[][][] sigma = new Double[C][C][C];
        for (int i = 0; i < C; i++) {
            for (int j = 0; j < C; j++) {
                for (int k = 0; k < C; k++) {
                    sigma[i][j][k] = 0.0;
                }
            }
        }

        Matrix alpha = Map_pie.map_pie(mmap);

        Matrix negD0inv = mmap.get(0).scale(-1.0).inv();

        for (int i = 0; i < C; i++) {
            Matrix starti = alpha.mult(negD0inv).mult(mmap.get(i + 2));

            for (int j = 0; j < C; j++) {
                Matrix startj = starti.mult(negD0inv).mult(mmap.get(j + 2));

                for (int h = 0; h < C; h++) {
                    Matrix result = startj.mult(negD0inv).mult(mmap.get(h + 2));
                    double sum = 0.0;
                    for (int row = 0; row < result.getNumRows(); row++) {
                        for (int col = 0; col < result.getNumCols(); col++) {
                            sum += result.get(row, col);
                        }
                    }
                    sigma[i][j][h] = sum;
                }
            }
        }

        return sigma;
    }

    /**
     * MatrixCell overload.
     */
    public static MatrixCell mmap_sigma2_cell(MatrixCell mmap) {
        int C = mmap.size() - 2;
        Double[][][] sigma = mmap_sigma2(mmap);

        MatrixCell result = new MatrixCell(C);
        for (int i = 0; i < C; i++) {
            Matrix matrix = new Matrix(C, C);
            for (int j = 0; j < C; j++) {
                for (int h = 0; h < C; h++) {
                    matrix.set(j, h, sigma[i][j][h]);
                }
            }
            result.set(i, matrix);
        }

        return result;
    }
}
