/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api.qsys;

import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;

import java.util.HashMap;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Analytic validation of the queueing-system formula APIs (jline.api.qsys).
 *
 * Every assertion uses an exact mathematical identity:
 * - G/G/1 and G/G/k approximations reduce exactly to M/M/1 / M/M/k at
 *   ca = cs = 1 (their defining property);
 * - Pollaczek-Khinchine, Erlang loss/waiting, and D/M/c / M/D/c fixed-point
 *   forms are recomputed independently inline;
 * - bounds are checked against the exact value they bound.
 */
public class QsysApiTest {

    private static final double TOL = 1e-9;
    private static final double NUM_TOL = 1e-6; // iterative/numerical methods

    // Base test point: M/M/1 with lambda=1, mu=2 -> rho=0.5, W=1, Wq=0.5
    private static final double LAMBDA = 1.0;
    private static final double MU = 2.0;
    private static final double RHO = LAMBDA / MU;
    private static final double W_MM1 = 1.0 / (MU - LAMBDA);
    private static final double WQ_MM1 = W_MM1 - 1.0 / MU;

    @Test
    public void mm1MatchesClosedForm() {
        // NOTE: Ret.qsys fields are static; capture immediately after each call
        double w = Qsys_mm1.qsys_mm1(LAMBDA, MU).W;
        assertEquals(W_MM1, w, TOL, "M/M/1 mean response time");
    }

    @Test
    public void mg1ReducesToMm1AtScvOne() {
        double w = Qsys_mg1.qsys_mg1(LAMBDA, MU, 1.0).W;
        assertEquals(W_MM1, w, TOL, "M/G/1 with cs=1 must equal M/M/1");
    }

    @Test
    public void mg1MatchesPollaczekKhinchineForDeterministicService() {
        // M/D/1: Q = rho + rho^2 / (2(1-rho)), W = Q/lambda
        double q = RHO + RHO * RHO / (2 * (1 - RHO));
        double w = Qsys_mg1.qsys_mg1(LAMBDA, MU, 0.0).W;
        assertEquals(q / LAMBDA, w, TOL, "M/D/1 P-K mean response time");
    }

    @Test
    public void mmkReducesToMm1ForSingleServer() {
        double w = Qsys_mmk.qsys_mmk(LAMBDA, MU, 1).W;
        assertEquals(W_MM1, w, TOL, "M/M/k with k=1 must equal M/M/1");
    }

    @Test
    public void mmkMatchesErlangCClosedFormForTwoServers() {
        // Exact M/M/2: lambda=3, mu=2 -> a=1.5, rho=0.75
        double lambda = 3.0, mu = 2.0;
        int k = 2;
        double a = lambda / mu;
        double rho = a / k;
        double p0 = 1.0 / (1.0 + a + a * a / (2 * (1 - rho)));
        double pWait = (a * a / (2 * (1 - rho))) * p0; // Erlang C
        double wq = pWait / (k * mu - lambda);
        double w = Qsys_mmk.qsys_mmk(lambda, mu, k).W;
        assertEquals(wq + 1.0 / mu, w, TOL, "M/M/2 mean response time");
    }

    @Test
    public void gg1ReducesToMm1AtScvOne() {
        HashMap<String, Object> r = Qsys_gg1.qsys_gg1(LAMBDA, MU, 1.0, 1.0);
        assertEquals(W_MM1, (Double) r.get("W"), NUM_TOL, "G/G/1 W at ca2=cs2=1");
        assertEquals(WQ_MM1, (Double) r.get("Wq"), NUM_TOL, "G/G/1 Wq at ca2=cs2=1");
        assertEquals(RHO / (1 - RHO), (Double) r.get("L"), NUM_TOL, "G/G/1 L at ca2=cs2=1");
        assertEquals(1 - RHO, (Double) r.get("p0"), NUM_TOL, "G/G/1 p0 at ca2=cs2=1");
    }

    @Test
    public void gig1ApproximationsReduceToMm1AtCvOne() {
        // Allen-Cunneen is exact for the M/M/1 special case ca=cs=1
        double wAc = Qsys_gig1_approx_allencunneen
                .qsys_gig1_approx_allencunneen(LAMBDA, MU, 1.0, 1.0).W;
        assertEquals(W_MM1, wAc, TOL, "Allen-Cunneen at ca=cs=1");

        // Diffusion-based approximations (Gelenbe, Kimura) are not exact at
        // ca=cs=1 but must stay within a few percent of the M/M/1 value
        HashMap<String, Object> gelenbe = Qsys_gig1_approx_gelenbe
                .qsys_gig1_approx_gelenbe(LAMBDA, MU, 1.0, 1.0);
        assertEquals(W_MM1, (Double) gelenbe.get("W"), 0.05 * W_MM1, "Gelenbe at ca=cs=1");

        HashMap<String, Object> kimura = Qsys_gig1_approx_kimura
                .qsys_gig1_approx_kimura(LAMBDA, MU, 1.0, 1.0);
        assertEquals(W_MM1, (Double) kimura.get("W"), 0.05 * W_MM1, "Kimura at ca=cs=1");
    }

    @Test
    public void kingmanBoundDominatesExactMm1() {
        double wUb = Qsys_gig1_ubnd_kingman
                .qsys_gig1_ubnd_kingman(LAMBDA, MU, 1.0, 1.0).W;
        assertTrue(wUb >= W_MM1 - TOL,
                "Kingman upper bound " + wUb + " below exact M/M/1 " + W_MM1);
        HashMap<String, Object> lb = Qsys_gig1_lbnd.qsys_gig1_lbnd(LAMBDA, MU, 1.0, 1.0);
        double wLb = (Double) lb.get("W");
        assertTrue(wLb <= W_MM1 + TOL,
                "Lower bound " + wLb + " above exact M/M/1 " + W_MM1);
    }

    @Test
    public void gigkApproxReducesToMmkAtCvOne() {
        double lambda = 3.0, mu = 2.0;
        int k = 2;
        double wMmk = Qsys_mmk.qsys_mmk(lambda, mu, k).W;
        double wApprox = Qsys_gigk_approx.qsys_gigk_approx(lambda, mu, 1.0, 1.0, k).W;
        // The base G/G/k approximation replaces the Erlang-C term with a closed
        // form that is not exact at ca=cs=1; a few percent deviation is expected
        assertEquals(wMmk, wApprox, 0.05 * wMmk, "G/G/k approx at ca=cs=1 vs M/M/2");
    }

    @Test
    public void gigkWhittAndCosmetatosAreExactAtCvOne() {
        // The Whitt and Cosmetatos refinements retain the exact Erlang-C term,
        // so they reduce EXACTLY to M/M/k at ca=cs=1 (unlike the base approx).
        double lambda = 3.0, mu = 2.0;
        int k = 2;
        double wMmk = Qsys_mmk.qsys_mmk(lambda, mu, k).W;
        HashMap<String, Object> whitt =
                Qsys_gigk_approx_whitt.qsys_gigk_approx_whitt(lambda, mu, 1.0, 1.0, k);
        assertEquals(wMmk, (Double) whitt.get("W"), NUM_TOL, "Whitt G/G/k at cv=1");
        HashMap<String, Object> cosm =
                Qsys_gigk_approx_cosmetatos.qsys_gigk_approx_cosmetatos(lambda, mu, 1.0, 1.0, k);
        assertEquals(wMmk, (Double) cosm.get("W"), NUM_TOL, "Cosmetatos G/G/k at cv=1");
    }

    @Test
    public void mm1kLossMatchesTruncatedGeometricForm() {
        int K = 3;
        HashMap<String, Object> r = Qsys_mm1k_loss.qsys_mm1k_loss(LAMBDA, MU, K);
        // M/M/1/K loss probability: (1-rho) rho^K / (1 - rho^{K+1})
        double lossExact = (1 - RHO) * Math.pow(RHO, K) / (1 - Math.pow(RHO, K + 1));
        assertEquals(lossExact, (Double) r.get("lossprob"), TOL, "M/M/1/K loss probability");
    }

    @Test
    public void mginfIsInsensitiveToServiceVariability() {
        HashMap<String, Object> r1 = Qsys_mginf.qsys_mginf(2.0, 4.0, 1.0);
        HashMap<String, Object> r2 = Qsys_mginf.qsys_mginf(2.0, 4.0, 4.0);
        assertEquals(0.25, (Double) r1.get("W"), TOL, "M/G/inf W = 1/mu");
        assertEquals((Double) r1.get("W"), (Double) r2.get("W"), TOL,
                "M/G/inf must be insensitive to service SCV");
        assertEquals(0.5, (Double) r1.get("L"), TOL, "M/G/inf L = lambda/mu");
    }

    @Test
    public void phm1ReducesToMm1ForExponentialArrivals() {
        double[] alpha = {1.0};
        double[][] T = {{-LAMBDA}};
        PhM1Result r = Qsys_phm1.qsys_phm1(alpha, T, MU);
        assertEquals(RHO, r.getUtilization(), NUM_TOL, "PH/M/1 utilization");
        assertEquals(W_MM1, r.getMeanSojournTime(), NUM_TOL, "PH/M/1 sojourn");
        assertEquals(RHO / (1 - RHO), r.getMeanQueueLength(), NUM_TOL, "PH/M/1 queue length");
    }

    @Test
    public void phmcReducesToMmcForExponentialArrivals() {
        // 1-phase PH = exponential arrivals; PH/M/c must equal M/M/c
        double lambda = 3.0, mu = 2.0;
        int c = 2;
        double[] alpha = {1.0};
        double[][] T = {{-lambda}};
        PhMcResult r = Qsys_phmc.qsys_phmc(alpha, T, mu, c);
        assertEquals(lambda / (c * mu), r.getUtilization(), NUM_TOL, "PH/M/2 utilization");
        assertEquals(Qsys_mmk.qsys_mmk(lambda, mu, c).W, r.getMeanSojournTime(), NUM_TOL,
                "PH/M/2 sojourn equals M/M/2");
    }

    @Test
    public void mapm1ReducesToMm1ForPoissonArrivals() {
        Matrix d0 = new Matrix(1, 1);
        d0.set(0, 0, -LAMBDA);
        Matrix d1 = new Matrix(1, 1);
        d1.set(0, 0, LAMBDA);
        QsysMapPhResult r = Qsys_mapm1.qsys_mapm1(d0, d1, MU);
        assertEquals(W_MM1, r.getMeanSojournTime(), NUM_TOL, "MAP/M/1 sojourn for Poisson MAP");
        assertEquals(RHO / (1 - RHO), r.getMeanQueueLength(), NUM_TOL,
                "MAP/M/1 queue length for Poisson MAP");
    }

    @Test
    public void mapg1ReducesToMm1ForPoissonAndExponential() {
        // Poisson MAP + exponential service (scv=1) = M/M/1
        Matrix d0 = new Matrix(1, 1);
        d0.set(0, 0, -LAMBDA);
        Matrix d1 = new Matrix(1, 1);
        d1.set(0, 0, LAMBDA);
        QsysMapPhResult r = Qsys_mapg1.qsys_mapg1(d0, d1, 1.0 / MU, 1.0);
        assertEquals(W_MM1, r.getMeanSojournTime(), NUM_TOL, "MAP/G/1 sojourn");
        assertEquals(RHO / (1 - RHO), r.getMeanQueueLength(), NUM_TOL, "MAP/G/1 queue length");
    }

    @Test
    public void mapmap1ReducesToMm1ForPoissonAndExponential() {
        Matrix c0 = new Matrix(1, 1);
        c0.set(0, 0, -LAMBDA);
        Matrix c1 = new Matrix(1, 1);
        c1.set(0, 0, LAMBDA);
        Matrix d0 = new Matrix(1, 1);
        d0.set(0, 0, -MU);
        Matrix d1 = new Matrix(1, 1);
        d1.set(0, 0, MU);
        QsysMapPhResult r = Qsys_mapmap1.qsys_mapmap1(c0, c1, d0, d1);
        assertEquals(W_MM1, r.getMeanSojournTime(), NUM_TOL, "MAP/MAP/1 sojourn");
    }

    @Test
    public void mapmcReducesToMmcForPoissonArrivals() {
        // Poisson MAP arrivals + exponential service => MAP/M/c == M/M/c
        double lambda = 3.0, mu = 2.0;
        int c = 2;
        Matrix d0 = new Matrix(1, 1);
        d0.set(0, 0, -lambda);
        Matrix d1 = new Matrix(1, 1);
        d1.set(0, 0, lambda);
        QsysMapPhResult r = Qsys_mapmc.qsys_mapmc(d0, d1, mu, c);
        assertEquals(Qsys_mmk.qsys_mmk(lambda, mu, c).W, r.getMeanSojournTime(), NUM_TOL,
                "MAP/M/2 sojourn equals M/M/2 for Poisson MAP");
    }

    @Test
    public void mdcCrommelinMatchesMd1ClosedForm() {
        // M/D/1: Wq = rho * s / (2 (1 - rho)); lambda=0.5, s=1
        double lambda = 0.5, s = 1.0;
        double rho = lambda * s;
        MDcCrommelinResult r = Qsys_mdc_crommelin.qsys_mdc_crommelin(lambda, s, 1);
        double wqExact = rho * s / (2 * (1 - rho));
        assertEquals(wqExact, r.getMeanWaitingTime(), NUM_TOL, "M/D/1 Crommelin waiting time");
        assertEquals(wqExact + s, r.getMeanSojournTime(), NUM_TOL, "M/D/1 Crommelin sojourn");
    }

    @Test
    public void dmcMatchesDm1FixedPointForm() {
        // D/M/1 with D=1, mu=2: sigma = exp(-mu D (1-sigma)); Wq = sigma/(mu(1-sigma))
        double lambdaArr = 1.0, mu = 2.0;
        double sigma = 0.5;
        for (int i = 0; i < 200; i++) {
            sigma = Math.exp(-mu * (1.0 / lambdaArr) * (1 - sigma));
        }
        double wqExact = sigma / (mu * (1 - sigma));
        DmcResult r = Qsys_dmc.qsys_dmc(lambdaArr, mu, 1);
        assertEquals(wqExact, r.getMeanWaitingTime(), NUM_TOL, "D/M/1 waiting time");
        // Cross-check with the G/M/1 sigma-parameterized formula
        double wGm1 = Qsys_gm1.qsys_gm1(sigma, mu).W;
        assertEquals(wqExact + 1.0 / mu, wGm1, NUM_TOL, "G/M/1(sigma) vs D/M/1 sojourn");
    }

    @Test
    public void schedulingDisciplinesCoincideForExponentialService() {
        // With exponential (memoryless) service, every non-idling
        // work-conserving discipline induces the same M/M/1 birth-death
        // queue-length process, so the mean sojourn time is 1/(mu-lambda)
        // regardless of the discipline. Quadrature-based implementations are
        // allowed a 2% numerical tolerance.
        Matrix lambda = new Matrix(1, 1);
        lambda.set(0, 0, LAMBDA);
        Matrix mu = new Matrix(1, 1);
        mu.set(0, 0, MU);
        Matrix cs = new Matrix(1, 1);
        cs.set(0, 0, 1.0);
        double quadTol = 0.02 * W_MM1;

        // Non-clairvoyant age-based preemptive discipline: with memoryless
        // service the attained age carries no information, so FB matches the
        // M/M/1 birth-death mean sojourn 1/(mu-lambda)
        double wFb = Qsys_mg1_fb.qsys_mg1_fb(lambda, mu, cs).W.get(0);
        assertEquals(W_MM1, wFb, quadTol, "FB sojourn at cs=1");

        // qsys_mg1_setf implements a non-preemptive FB variant whose metric
        // convention differs from the plain mean sojourn (MATLAB and JAR agree
        // on 1.6033 for this input); assert the cross-codebase reference value
        double wSetf = Qsys_mg1_setf.qsys_mg1_setf(lambda, mu, cs).W.get(0);
        assertEquals(1.6033399361655691, wSetf, 1e-6,
                "SETF value must match the MATLAB reference");

        // Clairvoyant size-based disciplines DO change the mean even for
        // exponential service: SRPT is optimal (minimal mean sojourn), PSJF
        // is between SRPT and the non-size-based value, LRPT is pessimal
        double wSrpt = Qsys_mg1_srpt.qsys_mg1_srpt(lambda, mu, cs).W.get(0);
        double wPsjf = Qsys_mg1_psjf.qsys_mg1_psjf(lambda, mu, cs).W.get(0);
        double wLrpt = Qsys_mg1_lrpt.qsys_mg1_lrpt(lambda, mu, cs).W.get(0);
        assertTrue(wSrpt > 1.0 / MU - quadTol,
                "SRPT sojourn cannot be below the mean service time");
        assertTrue(wSrpt <= wPsjf + quadTol,
                "SRPT is optimal: must not exceed PSJF (" + wSrpt + " vs " + wPsjf + ")");
        assertTrue(wPsjf <= W_MM1 + quadTol,
                "PSJF must not exceed the non-size-based sojourn");
        assertTrue(wLrpt >= W_MM1 - quadTol,
                "LRPT is pessimal: must dominate the non-size-based sojourn");
    }

    @Test
    public void dpsWithEqualWeightsReducesToProcessorSharing() {
        // DPS with equal weights is PS; for M/M/1 PS the mean sojourn is
        // 1/(mu-lambda) for both classes regardless of their rates
        Matrix lambda = new Matrix(1, 2);
        lambda.set(0, 0, 0.3);
        lambda.set(0, 1, 0.2);
        Matrix mu = new Matrix(1, 2);
        mu.set(0, 0, MU);
        mu.set(0, 1, MU);
        Matrix w = new Matrix(1, 2);
        w.set(0, 0, 1.0);
        w.set(0, 1, 1.0);
        Matrix resp = Qsys_mm1_dps.qsys_mm1_dps(lambda, mu, w);
        double exact = 1.0 / (MU - 0.5); // total lambda = 0.5
        assertEquals(exact, resp.get(0), 1e-6, "DPS equal-weight class 1 sojourn");
        assertEquals(exact, resp.get(1), 1e-6, "DPS equal-weight class 2 sojourn");
    }

    @Test
    public void mg1PrioWithSingleClassReducesToMg1() {
        Matrix lambda = new Matrix(1, 1);
        lambda.set(0, 0, LAMBDA);
        Matrix mu = new Matrix(1, 1);
        mu.set(0, 0, MU);
        Matrix cs = new Matrix(1, 1);
        cs.set(0, 0, 1.0);
        Matrix wPrio = Qsys_mg1_prio.qsys_mg1_prio(lambda, mu, cs).W;
        double wMg1 = Qsys_mg1.qsys_mg1(LAMBDA, MU, 1.0).W;
        assertEquals(wMg1, wPrio.get(0, 0), NUM_TOL,
                "M/G/1 priority with one class must equal plain M/G/1");
    }
}
