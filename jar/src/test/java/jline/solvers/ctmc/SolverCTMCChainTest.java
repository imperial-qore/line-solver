package jline.solvers.ctmc;

import jline.lang.processes.MarkovChain;
import jline.lang.processes.MarkovProcess;
import jline.solvers.SolverOptions;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Regression tests for chain mode, i.e. SolverCTMC built from a user-supplied
 * MarkovProcess (CTMC) or MarkovChain (DTMC) instead of a Network. The solver
 * then skips state-space generation and solves the given generator directly.
 *
 * The reference is the closed-form stationary vector of the two-state chains
 * used here. The MATLAB twin is line-test.git/test/testsCTMC/test_ctmc_chain.m and the Python
 * twin is python/tests/test_ctmc_chain.py.
 */
public class SolverCTMCChainTest {

    private static final double TOL = 1e-9;

    private static Matrix generator() {
        // Two-state CTMC with pi = [2/7, 5/7]
        Matrix Q = new Matrix(2, 2);
        Q.set(0, 0, -0.5);
        Q.set(0, 1, 0.5);
        Q.set(1, 0, 0.2);
        Q.set(1, 1, -0.2);
        return Q;
    }

    private static Matrix transitionMatrix() {
        // Two-state DTMC with pi = [0.8, 0.2]
        Matrix P = new Matrix(2, 2);
        P.set(0, 0, 0.9);
        P.set(0, 1, 0.1);
        P.set(1, 0, 0.4);
        P.set(1, 1, 0.6);
        return P;
    }

    @Test
    public void ctmcChainStationaryMatchesClosedForm() {
        SolverCTMC solver = new SolverCTMC(new MarkovProcess(generator()));
        Matrix pi = solver.getProbSys().probability;
        assertTrue(solver.isChainSolver());
        assertFalse(solver.isDiscreteChain());
        assertEquals(2.0 / 7.0, pi.get(0), TOL);
        assertEquals(5.0 / 7.0, pi.get(1), TOL);
    }

    @Test
    public void ctmcChainReturnsTheGivenGenerator() {
        Matrix Q = generator();
        SolverCTMC solver = new SolverCTMC(new MarkovProcess(Q));
        Matrix infGen = solver.getGenerator().infGen;
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                assertEquals(Q.get(i, j), infGen.get(i, j), TOL);
            }
        }
        // No state space attached, so the state space is the state indices
        Matrix space = solver.getStateSpace().stateSpace;
        assertEquals(2, space.getNumRows());
        assertEquals(1.0, space.get(0, 0), TOL);
        assertEquals(2.0, space.get(1, 0), TOL);
    }

    @Test
    public void ctmcChainProbOfSingleState() {
        Matrix space = new Matrix(2, 1);
        space.set(0, 0, 0);
        space.set(1, 0, 1);
        MarkovProcess ctmc = new MarkovProcess(generator(), true, space);
        SolverCTMC solver = new SolverCTMC(ctmc);
        Matrix state = new Matrix(1, 1);
        state.set(0, 0, 1);
        assertEquals(5.0 / 7.0, solver.getProb(state), TOL);
        state.set(0, 0, 7);
        assertThrows(RuntimeException.class, () -> solver.getProb(state));
    }

    // Cross-codebase reference: pi(1) of the CTMC and pi after 5 steps of the
    // DTMC, both from the uniform start. MATLAB, JAR and Python assert these
    // same numbers.
    private static final double[] PI_CTMC_T1 = {0.392125422241, 0.607874577759};
    private static final double[] PI_DTMC_K5 = {0.790625, 0.209375};

    @Test
    public void chainTransientMatchesTheCrossCodebaseReference() {
        SolverOptions ctmcOptions = SolverCTMC.defaultOptions();
        ctmcOptions.timespan[0] = 0;
        ctmcOptions.timespan[1] = 1;
        Matrix pit = new SolverCTMC(new MarkovProcess(generator()), ctmcOptions).getTranProbSys().probability;
        int last = pit.getNumRows() - 1;
        assertEquals(PI_CTMC_T1[0], pit.get(last, 0), 1e-5);
        assertEquals(PI_CTMC_T1[1], pit.get(last, 1), 1e-5);

        SolverOptions dtmcOptions = SolverCTMC.defaultOptions();
        dtmcOptions.timespan[0] = 0;
        dtmcOptions.timespan[1] = 5;
        Matrix pitd = new SolverCTMC(new MarkovChain(transitionMatrix()), dtmcOptions).getTranProbSys().probability;
        assertEquals(PI_DTMC_K5[0], pitd.get(5, 0), 1e-9);
        assertEquals(PI_DTMC_K5[1], pitd.get(5, 1), 1e-9);

        // the container methods answer the same distribution as the solver
        Matrix piC = new MarkovProcess(generator()).transientProb(null, 1.0);
        assertEquals(PI_CTMC_T1[0], piC.get(0), 1e-5);
        Matrix piD = new MarkovChain(transitionMatrix()).transientProb(null, 5);
        assertEquals(PI_DTMC_K5[0], piD.get(5, 0), 1e-9);
    }

    @Test
    public void ctmcChainTransientConvergesToStationary() {
        SolverOptions options = SolverCTMC.defaultOptions();
        options.timespan[0] = 0;
        options.timespan[1] = 100;
        SolverCTMC solver = new SolverCTMC(new MarkovProcess(generator()), options);
        Matrix pit = solver.getTranProbSys().probability;
        int last = pit.getNumRows() - 1;
        assertEquals(2.0 / 7.0, pit.get(last, 0), 1e-3);
        assertEquals(5.0 / 7.0, pit.get(last, 1), 1e-3);
    }

    @Test
    public void dtmcChainStationaryMatchesClosedForm() {
        MarkovChain dtmc = new MarkovChain(transitionMatrix());
        SolverCTMC solver = new SolverCTMC(dtmc);
        Matrix pi = solver.getProbSys().probability;
        assertTrue(solver.isDiscreteChain());
        assertEquals(0.8, pi.get(0), TOL);
        assertEquals(0.2, pi.get(1), TOL);
        // MarkovChain.solve is the twin of MarkovProcess.solve
        Matrix piDirect = dtmc.solve();
        assertEquals(0.8, piDirect.get(0), TOL);
        assertEquals(0.2, piDirect.get(1), TOL);
    }

    @Test
    public void dtmcChainGeneratorIsTheUniformizedOne() {
        Matrix P = transitionMatrix();
        SolverCTMC solver = new SolverCTMC(new MarkovChain(P));
        Matrix infGen = solver.getGenerator().infGen;
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                assertEquals(P.get(i, j) - (i == j ? 1 : 0), infGen.get(i, j), TOL);
            }
        }
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                assertEquals(P.get(i, j), solver.getTransMat().get(i, j), TOL);
            }
        }
    }

    @Test
    public void dtmcChainTransientAdvancesOneStepPerUnitTime() {
        SolverOptions options = SolverCTMC.defaultOptions();
        options.timespan[0] = 0;
        options.timespan[1] = 5;
        SolverCTMC solver = new SolverCTMC(new MarkovChain(transitionMatrix()), options);
        Matrix pit = solver.getTranProbSys().probability;
        assertEquals(6, pit.getNumRows());
        // Uniform start advanced five steps of P
        Matrix pik = new Matrix(1, 2);
        pik.set(0, 0, 0.5);
        pik.set(0, 1, 0.5);
        for (int k = 0; k < 5; k++) {
            pik = pik.mult(transitionMatrix());
        }
        assertEquals(pik.get(0, 0), pit.get(5, 0), TOL);
        assertEquals(pik.get(0, 1), pit.get(5, 1), TOL);
    }

    @Test
    public void chainSampleStaysInTheStateSpace() {
        SolverCTMC solver = new SolverCTMC(new MarkovProcess(generator()));
        jline.io.Ret.SampleResult sample = solver.sampleSys(50);
        Matrix states = (Matrix) sample.state;
        assertEquals(50, states.getNumRows());
        for (int k = 0; k < 50; k++) {
            double s = states.get(k, 0);
            assertTrue(s == 1.0 || s == 2.0);
        }
        for (int k = 1; k < 50; k++) {
            assertTrue(sample.t.get(k, 0) > sample.t.get(k - 1, 0));
        }
    }

    @Test
    public void markovProcessApiMethods() {
        MarkovProcess ctmc = new MarkovProcess(generator());
        // transient endpoint and time average both converge to stationary
        Matrix piT = ctmc.transientProb(null, 200.0);
        assertEquals(2.0 / 7.0, piT.get(0), 1e-6);
        assertEquals(5.0 / 7.0, piT.get(1), 1e-6);
        Matrix avg = ctmc.timeAverage(null, 200.0);
        assertEquals(2.0 / 7.0, avg.get(0), 1e-2);
        // the constructor repairs the diagonal, so the carried generator is valid
        assertTrue(ctmc.isFeasible());
        // the embedded jump chain drops the holding times, so its stationary
        // vector differs from the CTMC one
        MarkovChain embedded = ctmc.toEmbedded();
        assertEquals(0.0, embedded.getTransMat().get(0, 0), TOL);
        assertEquals(1.0, embedded.getTransMat().get(0, 1), TOL);
        assertEquals(0.5, embedded.solve().get(0), 1e-9);
        // sensitivity of the stationary vector conserves total probability
        Matrix dQ = new Matrix(2, 2);
        dQ.set(0, 0, -1.0);
        dQ.set(0, 1, 1.0);
        Matrix dpi = ctmc.sens(dQ);
        assertEquals(0.0, dpi.get(0) + dpi.get(1), 1e-8);
        // stochastic complement of a single state is a 1x1 generator, and the
        // full form carries the blocks it was built from
        java.util.List<Double> keep = new java.util.ArrayList<Double>();
        keep.add(0.0);
        assertEquals(1, ctmc.stochComp(keep).getNumRows());
        jline.solvers.ctmc.SolverCTMC.StochCompResult blocks = ctmc.stochCompFull(keep);
        assertEquals(-0.5, blocks.Q11.get(0, 0), TOL);
        assertEquals(0.5, blocks.Q12.get(0, 0), TOL);
        assertEquals(0.2, blocks.Q21.get(0, 0), TOL);
        assertEquals(-0.2, blocks.Q22.get(0, 0), TOL);
        // T = Q12*inv(-Q22)*Q21 is the return path, and the complement is Q11 + T
        assertEquals(blocks.Q12.get(0, 0) / 0.2 * blocks.Q21.get(0, 0), blocks.T.get(0, 0), TOL);
        assertEquals(blocks.Q11.get(0, 0) + blocks.T.get(0, 0), blocks.S.get(0, 0), TOL);
    }

    @Test
    public void markovChainApiMethods() {
        MarkovChain dtmc = new MarkovChain(transitionMatrix());
        Matrix pi0 = new Matrix(1, 2);
        pi0.set(0, 0, 0.5);
        pi0.set(0, 1, 0.5);
        Matrix pit = dtmc.transientProb(pi0, 5);
        assertEquals(6, pit.getNumRows());
        Matrix pik = new Matrix(pi0);
        for (int k = 0; k < 5; k++) {
            pik = pik.mult(transitionMatrix());
        }
        assertEquals(pik.get(0, 0), pit.get(5, 0), TOL);
        // from state 0 the target state 1 is hit after 1/0.1 steps on average
        Matrix h = dtmc.hittingTime(new int[]{1});
        assertEquals(10.0, h.get(0, 0), 1e-9);
        assertEquals(0.0, h.get(1, 0), TOL);
        assertTrue(dtmc.isFeasible());
        java.util.List<Integer> keep = new java.util.ArrayList<Integer>();
        keep.add(0);
        assertEquals(1, dtmc.stochComp(keep).getNumRows());
    }

    @Test
    public void markovProcessApiExtras() {
        MarkovProcess ctmc = new MarkovProcess(generator());
        // the two transient engines agree
        Matrix unif = ctmc.transientProb(null, 1.0, "unif");
        Matrix fox = ctmc.transientProb(null, 1.0, "foxglynn");
        assertEquals(unif.get(0), fox.get(0), 1e-9);
        // the relative solution is the stationary vector scaled by p(refstate)
        Matrix rel = ctmc.solveRelative(0);
        Matrix pi = ctmc.solve();
        double total = rel.get(0) + rel.get(1);
        assertEquals(1.0, rel.get(0), TOL);
        assertEquals(pi.get(0), rel.get(0) / total, 1e-9);
        // aggregation of a nearly-decomposable chain recovers the exact vector
        double[][] qb = {{-1.001, 1.0, 0.001, 0.0}, {1.0, -1.001, 0.0, 0.001},
                {0.001, 0.0, -1.001, 1.0}, {0.0, 0.001, 1.0, -1.001}};
        Matrix Qb = new Matrix(4, 4);
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                Qb.set(i, j, qb[i][j]);
            }
        }
        MarkovProcess block = new MarkovProcess(Qb);
        Matrix exact = block.solve();
        java.util.List<java.util.List<Integer>> MS = new java.util.ArrayList<java.util.List<Integer>>();
        MS.add(java.util.Arrays.asList(0, 1));
        MS.add(java.util.Arrays.asList(2, 3));
        String[] methods = {"courtois", "kms", "takahashi"};
        for (int m = 0; m < methods.length; m++) {
            jline.util.Triple<Matrix, Double, Double> agg = block.aggregate(MS, methods[m], null);
            for (int i = 0; i < 4; i++) {
                assertEquals(exact.get(i), agg.getFirst().get(i), 1e-3);
            }
            assertTrue(agg.getSecond() <= agg.getThird());
        }
        java.util.List<java.util.List<Integer>> MSS = new java.util.ArrayList<java.util.List<Integer>>();
        MSS.add(java.util.Arrays.asList(0, 1));
        jline.util.Triple<Matrix, Double, Double> multi = block.aggregate(MS, "multi", MSS);
        assertEquals(exact.get(0), multi.getFirst().get(0), 1e-3);
    }

    @Test
    public void markovChainApiExtras() {
        MarkovChain dtmc = new MarkovChain(transitionMatrix());
        jline.api.mc.Dtmc_stochcomp.DtmcStochCompResult blocks =
                dtmc.stochCompFull(java.util.Arrays.asList(0));
        assertEquals(0.9, blocks.P11.get(0, 0), TOL);
        assertEquals(0.1, blocks.P12.get(0, 0), TOL);
        assertEquals(0.4, blocks.P21.get(0, 0), TOL);
        assertEquals(0.6, blocks.P22.get(0, 0), TOL);
        // S = P11 + P12*inv(I-P22)*P21, a 1x1 stochastic matrix here
        assertEquals(blocks.P11.get(0, 0) + blocks.P12.get(0, 0) / (1 - blocks.P22.get(0, 0))
                * blocks.P21.get(0, 0), blocks.S.get(0, 0), TOL);
        // read as the randomized image of a CTMC, the chain relaxes to its
        // stationary vector as t grows
        Matrix piT = dtmc.transientUnif(null, 500.0);
        assertEquals(0.8, piT.get(0), 1e-6);
        assertEquals(0.2, piT.get(1), 1e-6);
    }

    @Test
    public void averageMetricsAreRefusedInChainMode() {
        SolverCTMC solver = new SolverCTMC(new MarkovProcess(generator()));
        assertThrows(RuntimeException.class, solver::getAvgTable);
        assertThrows(RuntimeException.class, solver::getAvg);
    }
}
