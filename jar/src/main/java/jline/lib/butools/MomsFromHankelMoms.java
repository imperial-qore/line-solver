package jline.lib.butools;

import java.util.ArrayList;
import java.util.List;

import jline.util.matrix.Matrix;

public final class MomsFromHankelMoms {
    private MomsFromHankelMoms() {}

    /**
     * Constructs a Hankel matrix from the given first column and last row.
     */
    private static Matrix hankelMH(double[] col, double[] row) {
        int n = col.length;
        int m = row.length;
        Matrix H = Matrix.zeros(n, m);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                if (i + j < n) {
                    H.set(i, j, col[i + j]);
                } else {
                    H.set(i, j, row[i + j - n + 1]);
                }
            }
        }
        return H;
    }

    /**
     * Returns the raw moments given the Hankel moments.
     *
     * Reference: http://en.wikipedia.org/wiki/Stieltjes_moment_problem
     */
    public static Matrix momsFromHankelMoms(Matrix hm) {
        int n = hm.length();
        List<Double> mList = new ArrayList<Double>();
        mList.add(hm.get(0));

        for (int i = 1; i < n; i++) {
            double[] mArr = new double[mList.size()];
            for (int k = 0; k < mList.size(); k++) mArr[k] = mList.get(k);
            int N;
            Matrix H;

            if (i % 2 == 0) {
                N = i / 2 + 1;
                double[] col = new double[N];
                for (int k = 0; k < N; k++) col[k] = mArr[k];
                double[] row = new double[N];
                for (int k = 0; k < N - 1; k++) {
                    row[k] = mArr[N - 1 + k];
                }
                row[N - 1] = 0.0;
                H = hankelMH(col, row);
            } else {
                N = (i + 1) / 2 + 1;
                double[] col = new double[N];
                col[0] = 1.0;
                for (int k = 1; k < N; k++) {
                    col[k] = mArr[k - 1];
                }
                double[] row = new double[N];
                for (int k = 0; k < N - 1; k++) {
                    row[k] = mArr[N - 2 + k];
                }
                row[N - 1] = 0.0;
                H = hankelMH(col, row);
            }

            // Solve for the unknown moment using cofactor expansion along the last row
            double h = hm.get(i);
            Matrix rH = Matrix.zeros(N - 1, N); // first N-1 rows of H
            for (int r = 0; r < N - 1; r++) {
                for (int c = 0; c < N; c++) {
                    rH.set(r, c, H.get(r, c));
                }
            }

            double lastCofactor = 0.0;
            for (int j = 0; j < N; j++) {
                Matrix sub = Matrix.zeros(N - 1, N - 1);
                for (int r = 0; r < N - 1; r++) {
                    int cc = 0;
                    for (int c = 0; c < N; c++) {
                        if (c != j) {
                            sub.set(r, cc, rH.get(r, c));
                            cc++;
                        }
                    }
                }
                double sign = ((N + j - 1) % 2 == 0) ? 1.0 : -1.0;
                double cofactor = sign * sub.det();

                if (j < N - 1) {
                    h -= cofactor * H.get(N - 1, j);
                } else {
                    lastCofactor = cofactor;
                }
            }

            mList.add(h / lastCofactor);
        }

        Matrix m = new Matrix(hm.getNumRows(), hm.getNumCols(), n);
        for (int i = 0; i < n; i++) {
            m.set(i, mList.get(i));
        }
        return m;
    }
}
