package jline.api.qsys;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import org.junit.jupiter.api.Test;

import jline.util.matrix.Matrix;

/**
 * The exact MAP/PH/c queue.
 *
 * THE ORACLES ARE INDEPENDENT OF THE IMPLEMENTATION. With Poisson arrivals and
 * exponential service the model collapses to M/M/c, whose queue length, delay
 * probability, mean wait, second moment and waiting-time CCDF all have
 * closed forms (Erlang-C); those are checked to 1e-9. With a correlated MMPP2
 * arrival and exponential service it must reproduce qsys_mapmc, which solves a
 * DIFFERENT chain (Q-MAM's level-dependent QBD, phases = arrival phases only)
 * and is itself validated against the MATLAB reference.
 */
public class QsysMapPhcTest {

    private static Matrix rowVec(double... v) {
        Matrix m = new Matrix(1, v.length);
        for (int i = 0; i < v.length; i++) m.set(0, i, v[i]);
        return m;
    }

    private static Matrix square(int n, double... v) {
        Matrix m = new Matrix(n, n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) m.set(i, j, v[i * n + j]);
        }
        return m;
    }

    @Test
    public void mmcCollapseMatchesErlangC() {
        final double lam = 1.2, mu = 1.0;
        final int c = 3;
        final Matrix wPoints = rowVec(0.1, 0.5, 1.0);
        QsysMapPhcResult r = Qsys_mapphc.qsys_mapphc(
                square(1, -lam), square(1, lam), rowVec(1.0), square(1, -mu), c, 500, 3, wPoints);

        final double a = lam / mu, rho = a / c;
        double sum = 0.0, fact = 1.0;
        for (int n = 0; n < c; n++) {
            if (n > 0) fact *= n;
            sum += Math.pow(a, n) / fact;
        }
        final double factC = fact * c;
        final double p0 = 1.0 / (sum + Math.pow(a, c) / (factC * (1 - rho)));
        final double cErl = Math.pow(a, c) / (factC * (1 - rho)) * p0;
        final double lq = cErl * rho / (1 - rho);

        assertEquals(lq + a, r.getMeanQueueLength(), 1e-9);
        assertEquals(lq / lam, r.getMeanWaitingTime(), 1e-9);
        assertEquals(cErl, r.getProbWait(), 1e-9);
        assertEquals(2 * cErl / Math.pow(c * mu - lam, 2), r.getWaitingTimeMoments().get(0, 1), 1e-9);
        for (int i = 0; i < 3; i++) {
            final double t = wPoints.get(0, i);
            assertEquals(cErl * Math.exp(-(c * mu - lam) * t), r.getWaitingTimeCCDF().get(0, i), 1e-9);
        }
        assertEquals(1, r.getPhaseCount());
    }

    @Test
    public void correlatedArrivalMatchesTheIndependentMapMcAlgorithm() {
        // MMPP2(1.8, 0.4, 0.15, 0.25)
        final Matrix D0 = square(2, -(1.8 + 0.15), 0.15, 0.25, -(0.4 + 0.25));
        final Matrix D1 = square(2, 1.8, 0.0, 0.0, 0.4);
        final int c = 3;
        QsysMapPhcResult r = Qsys_mapphc.qsys_mapphc(D0, D1, rowVec(1.0), square(1, -1.0), c);
        QsysMapPhResult ref = Qsys_mapmc.qsys_mapmc(D0, D1, 1.0, c);
        assertEquals(ref.getMeanQueueLength(), r.getMeanQueueLength(), 1e-6);
        assertEquals(ref.getMeanWaitingTime(), r.getMeanWaitingTime(), 1e-6);
        // the MATLAB reference reports this number through its own Q-MAM path
        assertEquals(1.521875, r.getMeanQueueLength(), 1e-5);
    }

    @Test
    public void phaseCountIsTheMultisetCountNotThePowerCount() {
        // Erlang-2 service, so ms = 2; with c = 3 the ordered space is 2^3 = 8
        // and the multiset space is binomial(2+3-1,3) = 4.
        final Matrix alpha = rowVec(1.0, 0.0);
        final Matrix S = square(2, -2.0, 2.0, 0.0, -2.0);
        QsysMapPhcResult r = Qsys_mapphc.qsys_mapphc(square(1, -0.9), square(1, 0.9), alpha, S, 3);
        assertEquals(4, r.getPhaseCount());
        // Little's law on the servers ties the two means together
        final double lambda = 0.9;
        assertEquals(r.getMeanQueueLength() / lambda, r.getMeanSojournTime(), 1e-8);
        assertTrue(r.getUtilization() > 0 && r.getUtilization() < 1);
    }
}
