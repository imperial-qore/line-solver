package jline.lib.butools;

import jline.util.matrix.Matrix;

public final class CheckMoments {
    private CheckMoments() {}

    /**
     * Checks if the given moment sequence is valid in the sense
     * that it belongs to a distribution with support (0,inf).
     *
     * This procedure checks the determinant of Delta_n
     * and Delta_n^(1) according to the Stieltjes moment problem.
     *
     * @param m The (raw) moments to check (starts with the first moment).
     *          Its length must be odd.
     * @param prec Entries with absolute value less than prec are considered to be zeros.
     * @return The result of the check.
     */
    public static boolean checkMoments(Matrix m, double prec) {
        if (m.length() % 2 == 0) {
            throw new IllegalArgumentException("CheckMoments: the number of moments must be odd!");
        }

        // Prepend 1.0 to moments (zeroth moment)
        Matrix moments = new Matrix(1, m.length() + 1, m.length() + 1);
        moments.set(0, 1.0);
        for (int i = 0; i < m.length(); i++) {
            moments.set(i + 1, m.get(i));
        }
        int N = (moments.length() / 2) - 1;

        for (int n = 0; n <= N; n++) {
            // Create Hankel matrices
            Matrix H = hankelMatrix(moments, n + 1, 0, n, 2 * n);
            Matrix H0 = hankelMatrix(moments, n + 1, 1, n + 1, 2 * n + 1);

            if (H.det() < -prec || H0.det() < -prec) {
                return false;
            }
        }

        return true;
    }

    public static boolean checkMoments(Matrix m) {
        return checkMoments(m, 1e-14);
    }

    /**
     * Creates a Hankel matrix from the given moment vector.
     */
    private static Matrix hankelMatrix(Matrix moments, int size, int firstRowStart, int lastColStart, int lastColEnd) {
        Matrix result = Matrix.zeros(size, size);

        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                int idx = i + j;
                if (idx < size) {
                    result.set(i, j, moments.get(firstRowStart + idx));
                } else {
                    result.set(i, j, moments.get(lastColStart + idx - size + 1));
                }
            }
        }

        return result;
    }
}
