package jline.solvers.mam;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import org.junit.jupiter.api.Test;

import jline.lang.Network;
import jline.lang.OpenClass;
import jline.lang.constant.SchedStrategy;
import jline.lang.nodes.Queue;
import jline.lang.nodes.Sink;
import jline.lang.nodes.Source;
import jline.lang.processes.Distribution;
import jline.lang.processes.Exp;
import jline.lang.processes.MAP;
import jline.lang.processes.RAP;
import jline.util.matrix.Matrix;

/**
 * Solver-level validation of Rational Arrival Process (RAP) components: RAP
 * arrival with exponential service, exponential arrival with RAP service, and
 * RAP arrival with RAP service, solved through {@link SolverMAM}.
 *
 * Two oracles are used and neither is a tolerance chosen to fit.
 *
 * 1. The MAP control is an exact identity. A RAP whose (H0, H1) happen to be
 *    nonnegative IS a Markovian arrival process, so SolverMAM applied to the
 *    RAP object must return the same numbers as SolverMAM applied to the
 *    equivalent MAP object. Both reach the exact MAP/MAP/1 queue, so what this
 *    pins down is the RAP side of the sn marshalling - that a RAP object
 *    delivers its (H0, H1) to the solver unaltered and is recognised as
 *    Markovian when it is - rather than an agreement between two algorithms.
 *    It is an identity, so the tolerance is numerical only, and the LDES
 *    interval in the MATLAB companion test is what ties the value to reality.
 *
 * 2. For a genuinely non-Markovian RAP (H0 carries a negative off-diagonal
 *    entry, so the representation is not a MAP) there is no closed form. The
 *    reference values below were produced by MATLAB, which is the ground truth
 *    of this codebase, and independently confirmed to sit inside a SolverLDES
 *    95% confidence interval built from 20 replications of 1e6 samples. The
 *    LDES RAP sampler is the conditional-vector recursion of
 *    jline.api.mam.Rap_sample.RapSampler, which shares no code with the matrix
 *    analytic path, so that confirmation is a real independent check.
 *
 * The non-Markovian RAP used here is
 *   H0 = [-9.9 -0.2; 0.1 -1.0] * 30/119,  H1 = [9.7 0.4; 0 0.9] * 30/119
 * which has mean 1, SCV 4.5448028674 and lag-1 autocorrelation 0.3432018450.
 * Its dominant left eigenvector is strictly positive ([1, 88.98]), which is
 * what keeps the conditional density of one sign over the reachable set and so
 * makes the representation a genuine point process rather than merely an
 * algebraically admissible (H0, H1) pair.
 *
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 */
public class RapQueueingModelTest {

    /** Identity between a RAP object and the equivalent MAP object: numerical only. */
    private static final double IDENTITY_TOL = 1e-9;

    /**
     * Agreement with the MATLAB reference. The two codebases run the same
     * algorithm on the same matrices, so they agree far beyond this; the margin
     * absorbs the level-truncation cutoff of the QBD, which is the only place
     * the ports can legitimately differ.
     */
    private static final double MATLAB_TOL = 1e-6;

    private static final double NM_SCALE = 30.0 / 119.0;

    private static Matrix mat(double[][] a) {
        Matrix m = new Matrix(a.length, a[0].length);
        for (int i = 0; i < a.length; i++) {
            for (int j = 0; j < a[i].length; j++) {
                m.set(i, j, a[i][j]);
            }
        }
        return m;
    }

    /** MMPP2 scaled to unit rate. Nonnegative, so this pair is a MAP. */
    private static Matrix mapH0() {
        return mat(new double[][] {{-10.1 / 5.5, 0.1 / 5.5}, {0.1 / 5.5, -1.1 / 5.5}});
    }

    private static Matrix mapH1() {
        return mat(new double[][] {{10.0 / 5.5, 0.0}, {0.0, 1.0 / 5.5}});
    }

    /** Genuinely non-Markovian RAP at unit rate: H0 has a negative off-diagonal. */
    private static Matrix nmH0(double scale) {
        return mat(new double[][] {{-9.9 * NM_SCALE * scale, -0.2 * NM_SCALE * scale},
                                   {0.1 * NM_SCALE * scale, -1.0 * NM_SCALE * scale}});
    }

    private static Matrix nmH1(double scale) {
        return mat(new double[][] {{9.7 * NM_SCALE * scale, 0.4 * NM_SCALE * scale},
                                   {0.0, 0.9 * NM_SCALE * scale}});
    }

    /**
     * Partner of (nmH0, nmH1) with the SAME H0 and H1 = (-H0*e)*pie, where pie
     * is the arrival-embedded equilibrium vector of the correlated process. The
     * embedded transition matrix (-H0)^-1 H1 is then rank one, so the process is
     * a renewal stream with lag-k autocorrelation zero at every lag, while the
     * marginal - and therefore the mean and the SCV - is unchanged.
     */
    private static Matrix nmH1Uncorrelated(double scale) {
        return mat(new double[][] {{8.232773109243697 * NM_SCALE * scale, 1.867226890756303 * NM_SCALE * scale},
                                   {0.7336134453781513 * NM_SCALE * scale, 0.16638655462184874 * NM_SCALE * scale}});
    }

    private static Matrix scaled(Matrix m, double c) {
        Matrix out = new Matrix(m.getNumRows(), m.getNumCols());
        for (int i = 0; i < m.getNumRows(); i++) {
            for (int j = 0; j < m.getNumCols(); j++) {
                out.set(i, j, m.get(i, j) * c);
            }
        }
        return out;
    }

    private static Network build(Distribution arrival, Distribution service) {
        Network m = new Network("rapqn");
        Source src = new Source(m, "Src");
        Queue q = new Queue(m, "Q1", SchedStrategy.FCFS);
        Sink sink = new Sink(m, "Snk");
        OpenClass cls = new OpenClass(m, "C1", 0);
        src.setArrival(cls, arrival);
        q.setService(cls, service);
        m.link(Network.serialRouting(src, q, sink));
        return m;
    }

    private static double qlen(Network m) throws Exception {
        SolverMAM solver = new SolverMAM(m);
        return solver.getAvgQLen().get(1, 0);
    }

    // ---------------------------------------------------------------------
    // 1. MAP control: the RAP object and the MAP object describe one process
    // ---------------------------------------------------------------------

    @Test
    public void rapArrivalWithNonnegativeMatricesMatchesMapArrival() throws Exception {
        double qMap = qlen(build(new MAP(mapH0(), mapH1()), new Exp(2.0)));
        double qRap = qlen(build(new RAP(mapH0(), mapH1()), new Exp(2.0)));
        assertEquals(qMap, qRap, IDENTITY_TOL,
                "a RAP with nonnegative (H0,H1) is a MAP, so the two must agree exactly");
        // Guard against the identity holding because both branches returned a
        // degenerate value: the queue must be a genuine correlated-arrival one,
        // well above the M/M/1 value of 1.0 at the same utilization.
        assertTrue(qMap > 3.0, "correlated MAP arrival must inflate the queue well past M/M/1");
    }

    @Test
    public void rapServiceWithNonnegativeMatricesMatchesMapService() throws Exception {
        double qMap = qlen(build(new Exp(1.0), new MAP(scaled(mapH0(), 2.0), scaled(mapH1(), 2.0))));
        double qRap = qlen(build(new Exp(1.0), new RAP(scaled(mapH0(), 2.0), scaled(mapH1(), 2.0))));
        assertEquals(qMap, qRap, IDENTITY_TOL,
                "RAP service and MAP service describe the same process here");
        assertTrue(qMap > 8.0, "correlated service must inflate the queue past M/M/1");
    }

    @Test
    public void rapArrivalAndRapServiceMatchTheMapPair() throws Exception {
        double qMap = qlen(build(new MAP(mapH0(), mapH1()),
                                 new MAP(scaled(mapH0(), 2.0), scaled(mapH1(), 2.0))));
        double qRap = qlen(build(new RAP(mapH0(), mapH1()),
                                 new RAP(scaled(mapH0(), 2.0), scaled(mapH1(), 2.0))));
        assertEquals(qMap, qRap, IDENTITY_TOL,
                "the RAP/RAP pair and the MAP/MAP pair describe the same queue");
        assertTrue(qMap > 16.0, "both streams correlated: the queue must be large");
    }

    // ---------------------------------------------------------------------
    // 2. Genuinely non-Markovian RAP against the MATLAB reference
    // ---------------------------------------------------------------------

    @Test
    public void nonMarkovianRapServiceMatchesMatlab() throws Exception {
        // MATLAB: SolverMAM(build(Exp(1), RAP(H0*2, H1*2))).getAvgTable().QLen(2)
        // = 6.36312174319464, which is qbd_raprap1({-1,1},{H0*2,H1*2}) to 12
        // digits, i.e. the model does reach the RAP/RAP/1 QBD and is not being
        // answered by a phase-type surrogate.
        double q = qlen(build(new Exp(1.0), new RAP(nmH0(2.0), nmH1(2.0))));
        assertEquals(6.36312174319464, q, MATLAB_TOL,
                "non-Markovian RAP service must reproduce the MATLAB reference");
    }

    @Test
    public void nonMarkovianRapServiceBeatsTheExponentialBaseline() throws Exception {
        // Sanity floor rather than a golden: a correlated non-Markovian service
        // at rho = 0.5 must sit far above the M/M/1 value of 1.0. An
        // implementation that collapsed the RAP to its rate would still satisfy
        // the equality above if the reference were ever regenerated from it, so
        // the floor is what makes that equality load-bearing.
        assertTrue(qlen(build(new Exp(1.0), new RAP(nmH0(2.0), nmH1(2.0)))) > 5.0,
                "a correlated non-Markovian service must inflate the queue past M/M/1");
    }

    // The RAP-arrival-with-RAP-service case is deliberately NOT pinned to a
    // value here. SolverMAM currently returns 10.80992052415988 for it in this
    // codebase and 10.80992052415996 in MATLAB, but SolverLDES puts it at
    // 16.13581322 +/- 0.10309226 over 20 replications of 2e6 samples, so the
    // matrix analytic answer is about a third low and pinning it would enshrine
    // a defect. The autocorrelation test below is what covers this model until
    // that is resolved.

    // ---------------------------------------------------------------------
    // 3. Representation validity agrees with MATLAB and Python
    // ---------------------------------------------------------------------

    @Test
    public void tiedDominantEigenvalueIsAccepted() {
        // The arrival block Ca of the Bean and Nielsen (2010) example has
        // eigenvalues -1 and -1 +/- i: the dominant real eigenvalue ties in real
        // part with a complex pair. BuTools, and therefore MATLAB
        // CheckRAPRepresentation.m and the Python port, accept this and only
        // note the tie. An earlier version of the JAR port picked the single
        // eigenvalue of smallest modulus of the real part and demanded that it
        // be real, which rejected this pair - and did so nondeterministically,
        // since the tie was broken by whatever order eig() returned. The
        // constructor throws on rejection, so reaching the assertion is the test.
        Matrix ca = mat(new double[][] {{-1.0, 0.0, 0.0},
                                        {-2.0 / 3, -1.0, 1.0},
                                        {2.0 / 3, -1.0, -1.0}});
        Matrix da = mat(new double[][] {{14.0 / 5, -9.0 / 10, -9.0 / 10},
                                        {26.0 / 15, -8.0 / 15, -8.0 / 15},
                                        {58.0 / 15, -19.0 / 15, -19.0 / 15}});
        RAP rap = new RAP(ca, da);
        assertTrue(rap.getNumberOfPhases() == 3, "the representation must be accepted with its 3 phases");
    }

    // ---------------------------------------------------------------------
    // 4. The autocorrelation must survive the sn marshalling
    // ---------------------------------------------------------------------

    @Test
    public void arrivalAutocorrelationChangesTheQueueLength() throws Exception {
        RAP correlated = new RAP(nmH0(1.0), nmH1(1.0));
        RAP uncorrelated = new RAP(nmH0(1.0), nmH1Uncorrelated(1.0));

        // The two processes are constructed to have the same marginal, so the
        // mean and the SCV must agree to machine precision. If they did not, a
        // difference in the solved queue length would prove nothing.
        assertEquals(correlated.getMean(), uncorrelated.getMean(), 1e-12,
                "the two arrival processes must share their mean");
        assertEquals(correlated.getSCV(), uncorrelated.getSCV(), 1e-12,
                "the two arrival processes must share their SCV");
        assertEquals(0.3432018450, correlated.getACF(new double[] {1.0}).get(0), 1e-9,
                "correlated arrival lag-1 autocorrelation");
        assertEquals(0.0, uncorrelated.getACF(new double[] {1.0}).get(0), 1e-12,
                "the renewal partner must have zero lag-1 autocorrelation");

        // Same service in both models, so any difference in the queue length is
        // attributable to the arrival autocorrelation alone. If the correlation
        // were dropped anywhere between the RAP object and the solver the two
        // would come out equal.
        double qCorr = qlen(build(correlated, new RAP(nmH0(2.0), nmH1(2.0))));
        double qUncorr = qlen(build(uncorrelated, new RAP(nmH0(2.0), nmH1(2.0))));
        // SolverMAM gives 10.8099 against 7.7742 here, a ratio of 1.3905;
        // SolverLDES puts the same ratio at 2.0013 (16.13581322 +/- 0.10309226
        // against 8.06276774 +/- 0.05791355), so the solver understates the
        // effect but does carry it. The threshold is 1.25: far enough above 1
        // that it cannot be met by numerical drift, which is at the 1e-6 level
        // here, and far enough below the observed 1.39 to be stable, while
        // still failing outright if the correlation were dropped in the
        // marshalling - in which case the two models would return equal values.
        assertTrue(qCorr > 1.25 * qUncorr,
                "arrival autocorrelation must inflate the queue, got "
                        + qCorr + " against " + qUncorr);
    }
}
