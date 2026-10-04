package jline.api.mam;

import jline.util.matrix.Matrix;

/**
 * True if a process (D0, D1) is a renewal process.
 *
 * <p>A renewal process embedded as a MAP has D1 = t * alpha (rank one): the
 * phase entered after an event does not depend on the phase the process was in
 * at that event, so successive interevent times are independent and the process
 * is fully described by its marginal.</p>
 *
 * <p>The test is on the algebraic structure alone, so it holds equally for a
 * MAP, a matrix-exponential (ME) and a rational arrival process (RAP): an ME
 * renewal stream satisfies it, a correlated MAP or RAP does not. That is what
 * makes it the right guard in front of any closed form that reads only the
 * marginal (a PH pair (pie, D0), a mean and an SCV), because such a form is
 * exact for a renewal input and silently discards the correlation otherwise.</p>
 *
 * <p>Mirrors MATLAB {@code mam_is_renewal_map.m} and python
 * {@code _is_renewal_map} in the MAM handler.</p>
 */
public final class Mam_is_renewal_map {

    private Mam_is_renewal_map() {
    }

    public static boolean mam_is_renewal_map(Matrix D0, Matrix D1) {
        int ns = D0.getNumRows();
        if (ns == 1) {
            return true;   // a one-phase process is memoryless, hence renewal
        }
        Matrix tExit = new Matrix(ns, 1);
        for (int i = 0; i < ns; i++) {
            double s = 0.0;
            for (int j = 0; j < ns; j++) {
                s += D1.get(i, j);
            }
            tExit.set(i, 0, s);
        }
        Matrix alpha = Map_pie.map_pie(D0, D1);
        double diff = 0.0;
        double normD1 = 0.0;
        for (int i = 0; i < ns; i++) {
            for (int j = 0; j < ns; j++) {
                double rankOne = tExit.get(i, 0) * alpha.get(j);
                double d = D1.get(i, j) - rankOne;
                diff += d * d;
                normD1 += D1.get(i, j) * D1.get(i, j);
            }
        }
        return Math.sqrt(diff) < 1e-9 * Math.max(1.0, Math.sqrt(normD1));
    }
}
