/**
 * @file Markovian Arrival Process autocorrelation function analysis
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Map_acf {
    private Map_acf() {}

    /**
     * Computes the autocorrelation function (ACF) for a given MAP at multiple lags.
     */
    public static Matrix map_acf(Matrix D0, Matrix D1, Matrix lags) {
        Matrix P = Map_embedded.map_embedded(D0, D1);
        Matrix x = Map_piq.map_piq(D0, D1);
        x.scaleEq(Map_lambda.map_lambda(D0, D1));

        Matrix neg_D0 = D0.copy();
        neg_D0.scaleEq(-1.0);
        Matrix y = neg_D0.inv().sumRows();
        Matrix acfCoeffs = new Matrix(1, lags.length(), lags.length());
        for (int i = 0; i < lags.length(); i++) {
            acfCoeffs.set(i, x.mult(Matrix.pow(P, (int) lags.get(i))).mult(y).get(0));
        }
        for (int i = 0; i < acfCoeffs.length(); i++) {
            acfCoeffs.set(i, acfCoeffs.get(i) - 1);
        }

        acfCoeffs.scaleEq(1 / Map_scv.map_scv(D0, D1));
        return acfCoeffs;
    }

    public static Matrix map_acf(MatrixCell MAP, Matrix lags) {
        return map_acf(MAP.get(0), MAP.get(1), lags);
    }

    public static Matrix map_acf(Matrix D0, Matrix D1) {
        Matrix LAGS = new Matrix(1, 1, 1);
        LAGS.set(0, 0, 1);
        return map_acf(D0, D1, LAGS);
    }

    public static Matrix map_acf(MatrixCell MAP) {
        return map_acf(MAP.get(0), MAP.get(1));
    }

    public static Matrix map_acf(Matrix D0, Matrix D1, int lag) {
        return map_acf(D0, D1, Matrix.singleton((double) lag));
    }
}
