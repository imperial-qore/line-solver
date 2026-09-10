/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.qsys;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.junit.jupiter.api.Test;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * MAP/G/1/K family against the M/M/1/K closed form and against MATLAB
 * (matlab/src/api/qsys/qsys_mapg1k.m and companions), which is the reference
 * implementation.
 */
public class QsysMapG1kTest {

    private static Matrix mat(double[][] a) {
        Matrix m = new Matrix(a.length, a[0].length);
        for (int i = 0; i < a.length; i++) {
            for (int j = 0; j < a[0].length; j++) {
                m.set(i, j, a[i][j]);
            }
        }
        return m;
    }

    /** The arrival MAP of MATLAB cases B, C and D. */
    private static Matrix d0() {
        return mat(new double[][]{{-1.5, 0.5}, {0.2, -0.9}});
    }

    private static Matrix d1() {
        return mat(new double[][]{{1.0, 0.0}, {0.0, 0.7}});
    }

    @Test
    public void mm1kMatchesTheClosedForm() {
        // M/M/1/K, lambda = 2, mu = 3: p_K = (1-rho) rho^K/(1-rho^(K+1))
        double lambda = 2.0;
        double mu = 3.0;
        double rho = lambda / mu;
        for (int K = 1; K <= 8; K++) {
            QsysMapG1kResult r = Qsys_mapg1k.qsys_mapg1k(mat(new double[][]{{-lambda}}),
                    mat(new double[][]{{lambda}}), QsysServiceLaw.exponential(mu), K);
            double exact = (1.0 - rho) * Math.pow(rho, K) / (1.0 - Math.pow(rho, K + 1));
            assertEquals(exact, r.pK, 1e-11, "pK at K=" + K);
            // Every level of the truncated geometric, not just the boundary
            double norm = (1.0 - rho) / (1.0 - Math.pow(rho, K + 1));
            for (int l = 0; l <= K; l++) {
                assertEquals(norm * Math.pow(rho, l), r.plevel.get(0, l), 1e-11, "plevel at K=" + K);
            }
            assertEquals(1.0 - r.throughput * r.meanServiceTime, r.p0, 1e-12, "p0 = 1 - T S");
            assertEquals(1.0 - r.throughput / r.lambda, r.lossProbability, 1e-14);
        }
    }

    @Test
    public void mm1k5AgreesWithMatlab() {
        QsysMapG1kResult r = Qsys_mapg1k.qsys_mapg1k(mat(new double[][]{{-2.0}}),
                mat(new double[][]{{2.0}}), QsysServiceLaw.gamma(1.0, 1.0 / 3.0), 5);
        assertEquals(0.0481203007518798, r.pK, 1e-12);
        assertEquals(0.365413533834586, r.p0, 1e-12);
        assertEquals(1.90375939849624, r.throughput, 1e-12);
        assertEquals(1.42255639097744, r.meanQueueLength, 1e-12);
        assertEquals(0.0481203007518797, r.lossProbability, 1e-12);
    }

    @Test
    public void mapDeterministic1K4AgreesWithMatlab() {
        QsysMapG1kResult r = Qsys_mapg1k.qsys_mapg1k(d0(), d1(), QsysServiceLaw.deterministic(0.5), 4);
        assertEquals(0.785714285714286, r.lambda, 1e-13);
        assertEquals(0.00298542064568351, r.pK, 1e-13);
        assertEquals(0.608427817222949, r.p0, 1e-13);
        assertEquals(0.783144365554101, r.throughput, 1e-13);
        assertEquals(0.518092851090509, r.meanQueueLength, 1e-13);
        double[] plevel = {0.60842781722294936, 0.28950631452107134, 0.08059648884418355,
                0.018483958766112186, 0.00298542064568351};
        for (int l = 0; l < plevel.length; l++) {
            assertEquals(plevel[l], r.plevel.get(0, l), 1e-13, "plevel[" + l + "]");
        }
    }

    @Test
    public void mapPhaseType1K6AgreesWithMatlab() {
        Matrix alpha = mat(new double[][]{{0.6, 0.4}});
        Matrix T = mat(new double[][]{{-3.0, 1.0}, {0.0, -2.0}});
        QsysMapG1kResult r = Qsys_mapg1k.qsys_mapg1k(d0(), d1(), QsysServiceLaw.phaseType(alpha, T), 6);
        assertEquals(0.5, r.meanServiceTime, 1e-14);
        assertEquals(0.00257223043669522, r.pK, 1e-13);
        assertEquals(0.60821501825072, r.p0, 1e-13);
        assertEquals(0.78356996349856, r.throughput, 1e-13);
        assertEquals(0.644309452450642, r.meanQueueLength, 1e-13);
    }

    @Test
    public void densityPathReproducesTheGammaPath() {
        // Same Gamma(2, 0.25) law given as a closed family and as a raw density:
        // two independent code paths, the second by quadrature.
        final double al = 2.0;
        final double th = 0.25;
        QsysServiceLaw closed = QsysServiceLaw.gamma(al, th);
        QsysServiceLaw byDensity = QsysServiceLaw.density(x ->
                Math.exp((al - 1.0) * Math.log(x) - x / th - al * Math.log(th)
                        - org.apache.commons.math3.special.Gamma.logGamma(al)));
        QsysMapG1kResult a = Qsys_mapg1k.qsys_mapg1k(d0(), d1(), closed, 5);
        QsysMapG1kResult b = Qsys_mapg1k.qsys_mapg1k(d0(), d1(), byDensity, 5);
        assertEquals(a.meanServiceTime, b.meanServiceTime, 1e-9);
        assertEquals(a.pK, b.pK, 1e-9);
        assertEquals(a.p0, b.p0, 1e-9);
        assertEquals(a.meanQueueLength, b.meanQueueLength, 1e-8);
    }

    @Test
    public void mmapPerClassLossAgreesWithMatlab() {
        List<Matrix> D1c = new ArrayList<Matrix>(Arrays.asList(
                mat(new double[][]{{0.8, 0.0}, {0.0, 0.2}}),
                mat(new double[][]{{0.2, 0.0}, {0.0, 0.5}})));
        QsysMmapG1kResult r = Qsys_mmapg1k.qsys_mmapg1k(d0(), D1c,
                QsysServiceLaw.gamma(2.0, 0.25), 5);
        assertEquals(0.3714285714285715, r.lambda.get(0, 0), 1e-13);
        assertEquals(0.41428571428571431, r.lambda.get(0, 1), 1e-13);
        assertEquals(0.0036665139707426175, r.lossRatio.get(0, 0), 1e-13);
        assertEquals(0.0023764014868457519, r.lossRatio.get(0, 1), 1e-13);
        assertEquals(0.37006672338229568, r.throughput.get(0, 0), 1e-13);
        assertEquals(0.41330120509830681, r.throughput.get(0, 1), 1e-13);
        assertEquals(0.00298627284286945, r.lossAggregate, 1e-13);
        assertEquals(0.581221776230446, r.meanQueueLength, 1e-12);
        // The two classes share a rate ordering but not a loss ratio: this is the
        // effect an aggregate-only analysis cannot express.
        assertTrue(r.lossRatio.get(0, 0) > r.lossRatio.get(0, 1));
    }

    @Test
    public void perflowAgreesWithMatlab() {
        List<MatrixCell> maps = new ArrayList<MatrixCell>();
        maps.add(new MatrixCell(mat(new double[][]{{-1.0, 0.2}, {0.3, -0.6}}),
                mat(new double[][]{{0.7, 0.1}, {0.1, 0.2}})));
        maps.add(new MatrixCell(mat(new double[][]{{-0.8}}), mat(new double[][]{{0.8}})));
        maps.add(new MatrixCell(mat(new double[][]{{-2.0, 0.4}, {0.1, -1.2}}),
                mat(new double[][]{{1.6, 0.0}, {0.2, 0.9}})));
        QsysMapG1kPerflowResult r = Qsys_mapg1k_perflow.qsys_mapg1k_perflow(maps,
                QsysServiceLaw.gamma(0.25, 0.4), 8);
        double[] lam = {0.58571428571428574, 0.80000000000000004, 1.3142857142857145};
        double[] tput = {0.5845756713286363, 0.79861017778923693, 1.3116030804338918};
        double[] loss = {0.001943975780377083, 0.0017372777634538927, 0.0020411344524737807};
        for (int n = 0; n < 3; n++) {
            assertEquals(lam[n], r.lambda.get(0, n), 1e-13, "lambda[" + n + "]");
            assertEquals(tput[n], r.throughput.get(0, n), 1e-11, "throughput[" + n + "]");
            assertEquals(loss[n], r.lossRatio.get(0, n), 1e-11, "lossRatio[" + n + "]");
        }
        assertEquals(0.10000000000000001, r.meanServiceTime, 1e-13);
        assertEquals(0.27000000000000002, r.rho, 1e-13);
        assertEquals(0.0019300260919389747, r.lossAggregate, 1e-11);
    }

    @Test
    public void retrialFixedPointAgreesWithMatlab() {
        Qsys_mmcc_retrial_fp.Result a = Qsys_mmcc_retrial_fp.qsys_mmcc_retrial_fp(8.0, 1.0, 10);
        assertEquals(0.23617844093433224, a.blocProb, 1e-14);
        assertEquals(2.4736504291936243, a.r, 1e-13);
        assertEquals(66, a.niter);
        Qsys_mmcc_retrial_fp.Result b = Qsys_mmcc_retrial_fp.qsys_mmcc_retrial_fp(3.0, 2.0, 4);
        assertEquals(0.055418042680785297, b.blocProb, 1e-14);
        assertEquals(0.17600815550208007, b.r, 1e-14);
        assertEquals(14, b.niter);
        // Erlang-B with no retrials is the c-server loss system
        assertEquals(0.0 + erlangBReference(1.0, 1), Qsys_mmcc_retrial_fp.erlangB(1.0, 1), 1e-15);
    }

    /** B(a,1) = a/(1+a), the M/M/1/1 loss probability. */
    private static double erlangBReference(double a, int c) {
        return a / (1.0 + a);
    }
}
