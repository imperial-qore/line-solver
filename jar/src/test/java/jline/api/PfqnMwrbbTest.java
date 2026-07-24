/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api;

import jline.api.pfqn.Pfqn_mwrbb;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Validates the Majumdar-Woodside robust box bounds
 * ({@link Pfqn_mwrbb#pfqn_mwrbb}) against the numerical results reported by
 * S. Majumdar and C.M. Woodside, "Robust bounds and throughput guarantees for
 * closed multiclass queueing networks", Performance Evaluation 32 (1998)
 * 101-136.
 *
 * <p>The reference model QNM1 (Section 3.1, Fig. 1-2) has two FIFO devices and
 * two closed classes with N1=2, N2=3, zero think time, and
 * V=[[1,1],[1,10]], S=[[0.9,0.9],[0.1,0.1]] (rows = devices, cols = classes).
 * The paper's robust box bounds are X1-=0.4, X1+=0.74075, X2-=0.37037,
 * X2+=0.71111 (the +0.74075 value is the paper's 5-digit rounding of
 * 0.740741).
 */
public class PfqnMwrbbTest {

    private static final double TOL = 1e-4;

    @Test
    public void testQNM1RobustBoxBounds() {
        Matrix V = new Matrix(2, 2);
        V.set(0, 0, 1.0); V.set(0, 1, 1.0);
        V.set(1, 0, 1.0); V.set(1, 1, 10.0);
        Matrix S = new Matrix(2, 2);
        S.set(0, 0, 0.9); S.set(0, 1, 0.9);
        S.set(1, 0, 0.1); S.set(1, 1, 0.1);
        Matrix N = new Matrix(1, 2);
        N.set(0, 0, 2.0); N.set(0, 1, 3.0);
        Matrix Z = new Matrix(1, 2);

        Pfqn_mwrbb.Result r = Pfqn_mwrbb.pfqn_mwrbb(V, S, N, Z);

        assertEquals(0.40000, r.Xlo.get(0, 0), TOL, "X1-");
        assertEquals(0.74074, r.Xup.get(0, 0), TOL, "X1+");
        assertEquals(0.37037, r.Xlo.get(0, 1), TOL, "X2-");
        assertEquals(0.71111, r.Xup.get(0, 1), TOL, "X2+");

        // the lower bound must not exceed the upper bound for any class
        for (int c = 0; c < 2; c++) {
            assertTrue(r.Xlo.get(0, c) <= r.Xup.get(0, c) + TOL,
                    "lower bound must not exceed upper bound for class " + c);
        }
    }

    /**
     * Validates the mixed-discipline lower bounds (Lemmas 2-5, Theorem 2)
     * against MQNM2, Table 2 (four classes, four devices, N=1 each, device 2 =
     * FIFO, device 3 = non-preemptive priority, device 4 = PS, device 1 varying
     * over FIFO/NPP/PP/PS). Priorities (LINE convention, lower = higher):
     * class 1 = 0, classes 2,3 = 1, class 4 = 2.
     */
    @Test
    public void testMQNM2MixedDisciplines() {
        double[][] Va = {
                {1, 1, 1, 1.0 / 3}, {0, 1, 1, 1.0 / 3},
                {0, 1, 1, 1.0 / 3}, {0, 4, 4, 2}};
        double[][] Sa = {
                {1, 2, 3, 10}, {0, 2, 3, 10}, {0, 2, 3, 10}, {0, 2, 3, 10}};
        Matrix V = new Matrix(4, 4);
        Matrix S = new Matrix(4, 4);
        for (int k = 0; k < 4; k++) {
            for (int c = 0; c < 4; c++) {
                V.set(k, c, Va[k][c]);
                S.set(k, c, Sa[k][c]);
            }
        }
        Matrix N = new Matrix(1, 4);
        for (int c = 0; c < 4; c++) N.set(0, c, 1.0);
        Matrix Z = new Matrix(1, 4);
        Z.set(0, 0, 2); Z.set(0, 1, 5); Z.set(0, 2, 5); Z.set(0, 3, 20);
        Matrix prio = new Matrix(1, 4);
        prio.set(0, 0, 0); prio.set(0, 1, 1); prio.set(0, 2, 1); prio.set(0, 3, 2);

        // device1 discipline code, then expected class lower bounds (Table 2)
        int[] dev1 = {Pfqn_mwrbb.FIFO, Pfqn_mwrbb.NPPRIO, Pfqn_mwrbb.PPPRIO, Pfqn_mwrbb.PS};
        double[][] expLo = {
                {0.2376, 0.01778, 0.01404, 0.00832},
                {0.2376, 0.01061, 0.008333, 0.002459},
                {0.3333, 0.01212, 0.009524, 0.002459},
                {0.2376, 0.01844, 0.01376, 0.007666}};
        for (int di = 0; di < 4; di++) {
            Matrix sched = new Matrix(4, 1);
            sched.set(0, 0, dev1[di]);
            sched.set(1, 0, Pfqn_mwrbb.FIFO);
            sched.set(2, 0, Pfqn_mwrbb.NPPRIO);
            sched.set(3, 0, Pfqn_mwrbb.PS);
            Pfqn_mwrbb.Result r = Pfqn_mwrbb.pfqn_mwrbb(V, S, N, Z, sched, prio);
            for (int c = 0; c < 4; c++) {
                assertEquals(expLo[di][c], r.Xlo.get(0, c), TOL,
                        "MQNM2 dev1=" + di + " class " + c);
            }
            // upper bounds are discipline-independent
            assertEquals(0.3333, r.Xup.get(0, 0), TOL);
            assertEquals(0.05263, r.Xup.get(0, 1), TOL);
            assertEquals(0.03846, r.Xup.get(0, 2), TOL);
            assertEquals(0.02, r.Xup.get(0, 3), TOL);
        }
    }

    /**
     * Validates the population dependence of the priority lower bounds against
     * MQNM1, Table 1 (device 1 = preemptive priority, device 2 = FIFO, device
     * 3 = non-preemptive priority, device 4 = PS; N1=N3=N4=1, N2 swept).
     */
    @Test
    public void testMQNM1PopulationSweep() {
        double[][] Va = {
                {1, 1, 1, 1.0 / 3}, {0, 1, 1, 1.0 / 3},
                {0, 1, 1, 1.0 / 3}, {0, 4, 4, 2}};
        double[][] Sa = {
                {1, 2, 3, 10}, {0, 2, 3, 10}, {0, 2, 3, 10}, {0, 2, 3, 10}};
        Matrix V = new Matrix(4, 4);
        Matrix S = new Matrix(4, 4);
        for (int k = 0; k < 4; k++) {
            for (int c = 0; c < 4; c++) {
                V.set(k, c, Va[k][c]);
                S.set(k, c, Sa[k][c]);
            }
        }
        Matrix Z = new Matrix(1, 4);
        Z.set(0, 0, 2); Z.set(0, 1, 5); Z.set(0, 2, 5); Z.set(0, 3, 20);
        Matrix prio = new Matrix(1, 4);
        prio.set(0, 0, 0); prio.set(0, 1, 1); prio.set(0, 2, 1); prio.set(0, 3, 2);
        Matrix sched = new Matrix(4, 1);
        sched.set(0, 0, Pfqn_mwrbb.PPPRIO);
        sched.set(1, 0, Pfqn_mwrbb.FIFO);
        sched.set(2, 0, Pfqn_mwrbb.NPPRIO);
        sched.set(3, 0, Pfqn_mwrbb.PS);

        int[] n2s = {1, 2, 3};
        double[][] expLo = {
                {0.3333, 0.01212, 0.009524, 0.002459},
                {0.3333, 0.01839, 0.007207, 0.0001321},
                {0.3333, 0.02222, 0.005952, 0.0}};
        double[][] expUp = {
                {0.3333, 0.05263, 0.03846, 0.02},
                {0.3333, 0.1053, 0.03846, 0.02},
                {0.3333, 0.1161, 0.03846, 0.02}};
        for (int t = 0; t < n2s.length; t++) {
            Matrix N = new Matrix(1, 4);
            N.set(0, 0, 1); N.set(0, 1, n2s[t]); N.set(0, 2, 1); N.set(0, 3, 1);
            Pfqn_mwrbb.Result r = Pfqn_mwrbb.pfqn_mwrbb(V, S, N, Z, sched, prio);
            for (int c = 0; c < 4; c++) {
                assertEquals(expLo[t][c], r.Xlo.get(0, c), TOL, "N2=" + n2s[t] + " Xlo class " + c);
                assertEquals(expUp[t][c], r.Xup.get(0, c), TOL, "N2=" + n2s[t] + " Xup class " + c);
            }
        }
    }

    /**
     * QNM1 with all stations set to the ABA full-contention code must give a
     * lower (looser) throughput guarantee than FIFO, since ABA forces the
     * arrival-contention probability P_cm = 1 (no rate-ratio discount).
     */
    @Test
    public void testAbaDisciplineLooserThanFifo() {
        Matrix V = new Matrix(2, 2);
        V.set(0, 0, 1.0); V.set(0, 1, 1.0); V.set(1, 0, 1.0); V.set(1, 1, 10.0);
        Matrix S = new Matrix(2, 2);
        S.set(0, 0, 0.9); S.set(0, 1, 0.9); S.set(1, 0, 0.1); S.set(1, 1, 0.1);
        Matrix N = new Matrix(1, 2); N.set(0, 0, 2.0); N.set(0, 1, 3.0);
        Matrix Z = new Matrix(1, 2);
        Matrix schedFifo = new Matrix(2, 1);
        schedFifo.set(0, 0, Pfqn_mwrbb.FIFO); schedFifo.set(1, 0, Pfqn_mwrbb.FIFO);
        Matrix schedAba = new Matrix(2, 1);
        schedAba.set(0, 0, Pfqn_mwrbb.ABA); schedAba.set(1, 0, Pfqn_mwrbb.ABA);

        Pfqn_mwrbb.Result rf = Pfqn_mwrbb.pfqn_mwrbb(V, S, N, Z, schedFifo, null);
        Pfqn_mwrbb.Result ra = Pfqn_mwrbb.pfqn_mwrbb(V, S, N, Z, schedAba, null);
        assertEquals(0.31579, ra.Xlo.get(0, 1), TOL, "class 2 ABA lower");
        for (int c = 0; c < 2; c++) {
            assertTrue(ra.Xlo.get(0, c) <= rf.Xlo.get(0, c) + 1e-12,
                    "ABA must not exceed FIFO lower for class " + c);
        }
    }

    /**
     * A single-class closed network must reduce to the Muntz-Wong asymptotic
     * bounds (paper Section 2.2), consistent with the well-known single-class
     * result X+ = min(1/Dmax, N/(Z+sum D)), X- = N/(Z+N sum D).
     */
    @Test
    public void testSingleClassReducesToAsymptotic() {
        double[] D = {0.9, 0.1};
        int Np = 4;
        double Zc = 1.0;
        Matrix V = new Matrix(2, 1);
        V.set(0, 0, 1.0); V.set(1, 0, 1.0);
        Matrix S = new Matrix(2, 1);
        S.set(0, 0, D[0]); S.set(1, 0, D[1]);
        Matrix N = new Matrix(1, 1); N.set(0, 0, Np);
        Matrix Z = new Matrix(1, 1); Z.set(0, 0, Zc);

        Pfqn_mwrbb.Result r = Pfqn_mwrbb.pfqn_mwrbb(V, S, N, Z);

        double Dsum = D[0] + D[1];
        double Dmax = Math.max(D[0], D[1]);
        double xUpRef = Math.min(1.0 / Dmax, Np / (Zc + Dsum));
        double xLoRef = Np / (Zc + Np * Dsum);
        assertEquals(xUpRef, r.Xup.get(0, 0), TOL, "single-class upper");
        assertEquals(xLoRef, r.Xlo.get(0, 0), TOL, "single-class lower");
    }
}
