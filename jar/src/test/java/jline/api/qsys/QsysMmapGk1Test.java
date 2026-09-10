package jline.api.qsys;

import static org.junit.jupiter.api.Assertions.assertEquals;

import java.util.ArrayList;
import java.util.List;

import org.junit.jupiter.api.Test;

import jline.lang.processes.Det;
import jline.lang.processes.Distribution;
import jline.lang.processes.Erlang;
import jline.lang.processes.Exp;
import jline.lang.processes.HyperExp;
import jline.lang.processes.Uniform;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

/**
 * Per-type waiting times of the MMAP[K]/G[K]/1 queue.
 *
 * THE ORACLES ARE INDEPENDENT OF THE IMPLEMENTATION. With one type, Poisson
 * arrivals and exponential service the queue is M/M/1, whose waiting time is an
 * atom plus an exponential tail in closed form; with deterministic service it is
 * M/D/1 and Pollaczek-Khinchine applies. The MULTICLASS case is checked against
 * the MATLAB reference, which reaches it through the same construction but is
 * itself validated against BuTools MMAPPH1FCFS (agreement 4e-16 on phase-type
 * service) and against JMT at 8e6 samples on general service (5e-4 per type).
 */
public class QsysMmapGk1Test {

    private static MatrixCell cell(Matrix... ms) {
        MatrixCell c = new MatrixCell(ms.length);
        for (int i = 0; i < ms.length; i++) c.set(i, ms[i]);
        return c;
    }

    private static Matrix scalar(double v) {
        Matrix m = new Matrix(1, 1);
        m.set(0, 0, v);
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
    public void mm1CollapseMatchesTheClosedForm() {
        final double lam = 0.6, mu = 1.0;
        List<Distribution> svc = new ArrayList<Distribution>();
        svc.add(new Exp(mu));
        Matrix pts = new Matrix(1, 2);
        pts.set(0, 0, 0.5);
        pts.set(0, 1, 2.0);
        QsysMmapGk1Result r = Qsys_mmapgk1.qsys_mmapgk1(
                cell(scalar(-lam), scalar(lam), scalar(lam)), svc, pts, 3, 1e-12, 10000);
        final double rho = lam / mu;
        assertEquals(rho / (mu - lam), r.getMeanWaitingTime().get(0, 0), 1e-9);
        assertEquals(2 * rho / Math.pow(mu - lam, 2), r.getWaitMoments().get(0, 1), 1e-9);
        assertEquals(1 - rho, r.getIdleVector().elementSum(), 1e-12);
        assertEquals(1 - rho * Math.exp(-(mu - lam) * 0.5), r.getWaitCDF().get(0, 0), 1e-6);
        assertEquals(1 - rho * Math.exp(-(mu - lam) * 2.0), r.getWaitCDF().get(0, 1), 1e-6);
    }

    @Test
    public void md1MatchesPollaczekKhinchine() {
        final double d = 0.8, lam = 0.9;
        List<Distribution> svc = new ArrayList<Distribution>();
        svc.add(new Det(d));
        QsysMmapGk1Result r = Qsys_mmapgk1.qsys_mmapgk1(
                cell(scalar(-lam), scalar(lam), scalar(lam)), svc);
        assertEquals(lam * d * d / (2 * (1 - lam * d)), r.getMeanWaitingTime().get(0, 0), 1e-9);
    }

    @Test
    public void twoTypesWithPhaseTypeServiceMatchTheReference() {
        final double l1 = 0.3, l2 = 0.25;
        List<Distribution> svc = new ArrayList<Distribution>();
        svc.add(Erlang.fitMeanAndSCV(1.0, 0.5));
        svc.add(HyperExp.fitMeanAndSCV(0.8, 3.0));
        QsysMmapGk1Result r = Qsys_mmapgk1.qsys_mmapgk1(
                cell(scalar(-(l1 + l2)), scalar(l1 + l2), scalar(l1), scalar(l2)), svc);
        assertEquals(2.09, r.getMeanSojournTime().get(0, 0), 1e-8);
        assertEquals(1.89, r.getMeanSojournTime().get(0, 1), 1e-8);
    }

    @Test
    public void twoTypesWithGeneralServiceAndCorrelatedArrivalsMatchTheReference() {
        // MMPP2(1.2, 0.3, 0.2, 0.4) marked by phase: type 1 in phase 1, type 2 in phase 2
        Matrix D0 = square(2, -(1.2 + 0.2), 0.2, 0.4, -(0.3 + 0.4));
        Matrix D1 = square(2, 1.2, 0.0, 0.0, 0.0);
        Matrix D2 = square(2, 0.0, 0.0, 0.0, 0.3);
        List<Distribution> svc = new ArrayList<Distribution>();
        svc.add(new Det(0.5));
        svc.add(new Uniform(0.1, 0.9));
        QsysMmapGk1Result r = Qsys_mmapgk1.qsys_mmapgk1(
                cell(D0, D1.add(D2), D1, D2), svc);
        assertEquals(0.319282, r.getMeanWaitingTime().get(0, 0), 1e-5);
        assertEquals(0.105970, r.getMeanWaitingTime().get(0, 1), 1e-5);
    }
}
