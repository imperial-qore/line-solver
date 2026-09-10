package jline.api.mam;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.ArrayList;
import java.util.List;

import org.junit.jupiter.api.Test;

import jline.lib.smc.MG1_ETAQA;
import jline.solvers.mam.handlers.BMAPMAP1Result;
import jline.solvers.mam.handlers.MAPBMAP1Result;
import jline.solvers.mam.handlers.Solver_mam_bmap_map_1;
import jline.solvers.mam.handlers.Solver_mam_map_bmap_1;
import jline.util.matrix.Matrix;

/**
 * Validates the ETAQA mean measures of the BMAP/MAP/1 (M/G/1-type) and
 * MAP/BMAP/1 (GI/M/1-type) queues against closed forms and against MATLAB.
 *
 * Both entry points used to answer with SUBSTITUTES rather than MAMSolver's
 * ETAQA: {@code mg1_g_etaqa} ran a plain functional iteration where the
 * reference runs Bini-Meini cyclic reduction, {@code mg1_qlen_etaqa} reported
 * sum(pi1) + 2 sum(pi*) -- which assumes every level at or above 2 IS level 2
 * -- and the GI/M/1 side never reached GIM1_R_ETAQA at all. Two live defects
 * came out with them: {@code MG1_pi_ETAQA} read the column vector returned by
 * {@code sumRows()} as a row and threw "Outside of matrix bounds" for every
 * boundary with more than one row, and its drift test used the COLUMN sums
 * where the reference forms the ROW sums.
 *
 * The oracles are of two kinds. A unit-batch BMAP with Poisson arrivals and
 * exponential service is an M/M/1, whose level is geometric, so the mean and
 * the higher moments are exact numbers; a genuine batch of size 1 or 2 is an
 * M[X]/M/1, which the unit-batch case cannot distinguish. The four-phase cases
 * are pinned to MATLAB, INCLUDING the negative GI/M/1 mean: MATLAB's
 * GIM1_qlen_ETAQA initializes its accumulator with the SCALAR A(3) where the
 * third BLOCK is meant, which corrupts the last column of the moment system
 * for more than one phase. That is reproduced deliberately, and pi and R --
 * which the defect does not touch -- are pinned alongside it.
 */
public class BmapEtaqaTest {

    static Matrix m2(double a, double b, double c, double d) {
        Matrix m = new Matrix(2, 2);
        m.set(0, 0, a);
        m.set(0, 1, b);
        m.set(1, 0, c);
        m.set(1, 1, d);
        return m;
    }

    static Matrix m1(double a) {
        Matrix m = new Matrix(1, 1);
        m.set(0, 0, a);
        return m;
    }

    static Matrix arrD0() { return m2(-1.4, 0.2, 0.3, -0.8); }

    static Matrix arrD1() { return m2(1.2, 0.0, 0.0, 0.5); }

    static Matrix svcD0() { return m2(-3.0, 0.5, 0.1, -2.0); }

    static Matrix svcD1() { return m2(2.0, 0.5, 1.4, 0.5); }

    static List<Matrix> list(Matrix... ms) {
        List<Matrix> out = new ArrayList<Matrix>();
        for (int i = 0; i < ms.length; i++) {
            out.add(ms[i]);
        }
        return out;
    }

    @Test
    public void unitBatchBmapIsAnMM1() {
        BMAPMAP1Result r = Solver_mam_bmap_map_1.solver_mam_bmap_map_1(
                list(m1(-0.6), m1(0.6)), m1(-1.0), m1(1.0));
        double rho = 0.6;
        assertEquals(rho / (1 - rho), r.meanQueueLength, 1e-10);
        assertEquals(rho, r.utilization, 1e-12);
        assertEquals(0.6, r.throughput, 1e-12);
        assertEquals(2.5, r.meanResponseTime, 1e-10);
        // A rank-one A0 sends MG1_EG home with G = 1 before cyclic reduction runs.
        assertEquals(1.0, r.G.get(0, 0), 1e-12);
        // The aggregates are (1-rho, rho(1-rho), rho^2).
        assertEquals(0.4, r.pi.get(0, 0), 1e-10);
        assertEquals(0.24, r.pi.get(0, 1), 1e-10);
        assertEquals(0.36, r.pi.get(0, 2), 1e-10);
    }

    @Test
    public void genuineBatchMatchesTheMxM1MeanQueueLength() {
        // Batch sizes 1 and 2 with equal probability, batch rate 0.3, lambda 0.45.
        BMAPMAP1Result r = Solver_mam_bmap_map_1.solver_mam_bmap_map_1(
                list(m1(-0.3), m1(0.15), m1(0.15)), m1(-1.0), m1(1.0));
        double rho = 0.45, ex = 1.5, exx1 = 1.0;
        double exact = rho / (1 - rho) + rho * (exx1 / ex) / (2 * (1 - rho));
        assertEquals(exact, r.meanQueueLength, 1e-10);
        assertEquals(0.45, r.throughput, 1e-12);
    }

    @Test
    public void bmapMap1MatchesMatlabOnAMapInput() {
        BMAPMAP1Result r = Solver_mam_bmap_map_1.solver_mam_bmap_map_1(
                list(arrD0(), arrD1().scale(0.5), arrD1().scale(0.5)), svcD0(), svcD1());
        assertEquals(2.496123595502, r.meanQueueLength, 1e-9);
        assertEquals(0.610619469027, r.utilization, 1e-10);
        assertEquals(1.808785214132, r.meanResponseTime, 1e-9);
        assertEquals(1.38, r.throughput, 1e-12);

        double mass = 0.0;
        for (int j = 0; j < r.pi.getNumCols(); j++) {
            mass += r.pi.get(0, j);
        }
        assertEquals(1.0, mass, 1e-12);
        assertEquals(0.118941068513, r.pi.get(0, 0), 1e-9);
        assertEquals(0.052543965796, r.pi.get(0, 4), 1e-9);
        assertEquals(0.188514965691, r.pi.get(0, 8), 1e-9);

        // G is stochastic: a positive recurrent chain leaves a level downwards
        // with probability one.
        for (int i = 0; i < r.G.getNumRows(); i++) {
            double s = 0.0;
            for (int j = 0; j < r.G.getNumCols(); j++) {
                s += r.G.get(i, j);
            }
            assertEquals(1.0, s, 1e-9);
        }
        assertEquals(0.686039381261, r.G.get(0, 0), 1e-9);
        assertEquals(0.217998974669, r.G.get(3, 3), 1e-9);
    }

    @Test
    public void higherMomentsMatchMatlab() {
        int ma = 2, ms = 2, m = ma * ms, K = 2;
        Matrix d0 = arrD0(), d1 = arrD1().scale(0.5);
        Matrix s0 = svcD0(), s1 = svcD1();
        Matrix Ima = Matrix.eye(ma), Ims = Matrix.eye(ms);

        Matrix A = new Matrix(m, m * (K + 2));
        put(A, Ima.kron(s1), 0);
        put(A, d0.kron(Ims).add(Ima.kron(s0)), m);
        put(A, d1.kron(Ims), 2 * m);
        put(A, d1.kron(Ims), 3 * m);

        Matrix B = new Matrix(m, m * (K + 1));
        put(B, d0.kron(Ims).add(Ima.kron(s0.add(s1))), 0);
        put(B, d1.kron(Ims), m);
        put(B, d1.kron(Ims), 2 * m);

        Matrix G = MG1_ETAQA.mg1_g_etaqa(A);
        Matrix pi = MG1_ETAQA.mg1_pi_etaqa(B, A, G, null);
        assertEquals(2.496123595502, MG1_ETAQA.mg1_qlen_etaqa(B, A, pi, 1), 1e-9);
        assertEquals(17.835504350850, MG1_ETAQA.mg1_qlen_etaqa(B, A, pi, 2), 1e-8);
        assertEquals(191.796380027147, MG1_ETAQA.mg1_qlen_etaqa(B, A, pi, 3), 1e-7);
    }

    static void put(Matrix dest, Matrix block, int colOffset) {
        for (int i = 0; i < block.getNumRows(); i++) {
            for (int j = 0; j < block.getNumCols(); j++) {
                dest.set(i, colOffset + j, block.get(i, j));
            }
        }
    }

    @Test
    public void unitBatchServiceIsAnMM1() {
        MAPBMAP1Result r = Solver_mam_map_bmap_1.solver_mam_map_bmap_1(
                m1(-0.6), m1(0.6), list(m1(-1.0), m1(1.0)));
        assertEquals(1.5, r.meanQueueLength, 1e-10);
        assertEquals(2.5, r.meanResponseTime, 1e-10);
        assertEquals(0.6, r.R.get(0, 0), 1e-10);
        assertEquals(0.4, r.pi.get(0, 0), 1e-10);
        assertEquals(0.36, r.pi.get(0, 2), 1e-10);
    }

    @Test
    public void singlePhaseBatchService() {
        MAPBMAP1Result r = Solver_mam_map_bmap_1.solver_mam_map_bmap_1(
                m1(-0.9), m1(0.9), list(m1(-1.0), m1(0.5), m1(0.5)));
        assertEquals(2.0611000442, r.meanQueueLength, 1e-9);
        assertEquals(0.6, r.utilization, 1e-12);
        assertEquals(0.6733200531, r.R.get(0, 0), 1e-9);
    }

    @Test
    public void mapBmap1MatchesMatlabIncludingTheNegativeMean() {
        MAPBMAP1Result r = Solver_mam_map_bmap_1.solver_mam_map_bmap_1(
                arrD0(), arrD1(), list(svcD0(), svcD1().scale(0.5), svcD1().scale(0.5)));
        // Negative, and negative in MATLAB too, to twelve digits: see the class
        // comment and GIM1_ETAQA.
        assertEquals(-3.434197005219, r.meanQueueLength, 1e-9);
        assertEquals(0.271386430678, r.utilization, 1e-10);
        assertEquals(0.92, r.throughput, 1e-12);

        double mass = 0.0;
        for (int j = 0; j < r.pi.getNumCols(); j++) {
            mass += r.pi.get(0, j);
        }
        assertEquals(1.0, mass, 1e-12);
        assertEquals(0.227094812985, r.pi.get(0, 0), 1e-9);
        assertEquals(0.082609049829, r.pi.get(0, 4), 1e-9);
        assertEquals(0.050296137186, r.pi.get(0, 8), 1e-9);

        assertEquals(0.319301430437, r.R.get(0, 0), 1e-9);
        assertEquals(0.187111846865, r.R.get(3, 3), 1e-9);
        assertTrue(r.R.get(0, 0) > 0.0);
    }
}
