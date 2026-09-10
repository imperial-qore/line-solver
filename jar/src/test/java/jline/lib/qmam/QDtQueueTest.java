package jline.lib.qmam;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;

import org.junit.jupiter.api.Test;

import jline.util.matrix.Matrix;

/**
 * Discrete-time Q-MAM queues, pinned to the closed forms and to the MATLAB
 * numbers recorded in _kb/06-solver-catalog.md.
 *
 * <p>The Geo/Geo/1 targets come from dqsys_geogeo1 under the LAS_DA
 * convention: with A = 0.2 and S = 0.5, E[N] = A(1-A)/(S-A) = 0.5333... and
 * P[empty] = 1 - A/S = 0.6. The DMAP row reproduces the MATLAB
 * SolverMAM discrete-time path, which LDES slotted independently measured at
 * 0.698332 over 400k slots.
 */
public class QDtQueueTest {

    private static final double TOL = 1e-8;

    private static Matrix row(double... v) {
        Matrix m = new Matrix(1, v.length);
        for (int i = 0; i < v.length; i++) {
            m.set(0, i, v[i]);
        }
        return m;
    }

    private static Matrix sq(int n, double... v) {
        Matrix m = new Matrix(n, n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                m.set(i, j, v[i * n + j]);
            }
        }
        return m;
    }

    @Test
    public void geoGeo1MatchesClosedForm() {
        double a = 0.2, s = 0.5;
        DTQueueResult r = Q_DT_MAP_MAP_1.qDtMapMap1(
                sq(1, 1 - a), sq(1, a), sq(1, 1 - s), sq(1, s), null);

        assertEquals(a * (1 - a) / (s - a), r.getMeanQueueLength(), 1e-6,
                "mean number in system under LAS-DA");
        assertEquals(1 - a / s, r.getQueueLength().get(0, 0), 1e-6, "probability of an empty system");
        assertEquals(a / s, r.getUtilization(), 1e-6, "utilization equals the load");
    }

    @Test
    public void phPh1AgreesWithMapMap1OnAGeometricPair() {
        double a = 0.2, s = 0.5;
        DTQueueResult viaMap = Q_DT_MAP_MAP_1.qDtMapMap1(
                sq(1, 1 - a), sq(1, a), sq(1, 1 - s), sq(1, s), null);
        DTQueueResult viaPh = Q_DT_PH_PH_1.qDtPhPh1(
                row(1.0), sq(1, 1 - a), row(1.0), sq(1, 1 - s), null);

        assertEquals(viaMap.getMeanQueueLength(), viaPh.getMeanQueueLength(), TOL,
                "the renewal DPH route must reproduce the DMAP route exactly");
    }

    @Test
    public void detServiceMatchesTheMatlabPath() {
        // Geo(0.2) arrivals, Det(2) service as a two-phase DPH: alpha = e1,
        // T(1,2) = 1. MATLAB SolverMAM reports 0.466667, LDES slotted 0.466615.
        double a = 0.2;
        Matrix T = sq(1, 1 - a);
        Matrix S = new Matrix(2, 2);
        S.set(0, 1, 1.0);
        DTQueueResult r = Q_DT_PH_PH_1.qDtPhPh1(row(1.0), T, row(1.0, 0.0), S, null);

        assertEquals(0.466667, r.getMeanQueueLength(), 1e-5, "Geo/Det(2)/1 mean queue length");
        assertEquals(0.4, r.getUtilization(), 1e-6, "Geo/Det(2)/1 utilization");
    }

    @Test
    public void dmapArrivalsMatchTheMatlabPath() {
        // The MATLAB SolverMAM discrete-time path reports 0.700000 with a
        // utilization of 0.5 on this pair; LDES slotted measured 0.698332.
        Matrix D0 = sq(2, 0.5, 0.2, 0.1, 0.6);
        Matrix D1 = sq(2, 0.25, 0.05, 0.1, 0.2);
        double s = 0.6;
        DTQueueResult r = Q_DT_MAP_MAP_1.qDtMapMap1(D0, D1, sq(1, 1 - s), sq(1, s), null);

        assertEquals(0.700000, r.getMeanQueueLength(), 1e-5, "DMAP/Geo/1 mean queue length");
        assertEquals(0.5, r.getUtilization(), 1e-6, "DMAP/Geo/1 utilization");
    }

    @Test
    public void overloadIsRefusedRatherThanReturned() {
        assertThrows(RuntimeException.class, () -> Q_DT_MAP_MAP_1.qDtMapMap1(
                sq(1, 0.4), sq(1, 0.6), sq(1, 0.7), sq(1, 0.3), null),
                "a load above one must be refused, not reported");
    }
}
