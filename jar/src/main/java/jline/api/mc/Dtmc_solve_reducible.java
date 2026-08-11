package jline.api.mc;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import jline.io.Ret;
import jline.util.Pair;
import jline.util.graph.DirectedGraph;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Dtmc_solve_reducible {
    private Dtmc_solve_reducible() {}

    /**
     * Result class for DTMC solve reducible.
     */
    public static final class DtmcSolveReducibleResult {
        public final Matrix pi;
        public final List<Matrix> pis;
        public final Matrix pi0;
        public final List<List<Integer>> scc;
        public final List<Boolean> isrec;
        public final Matrix Pl;
        public final Matrix pil;

        public DtmcSolveReducibleResult(Matrix pi, List<Matrix> pis, Matrix pi0,
                                         List<List<Integer>> scc, List<Boolean> isrec,
                                         Matrix Pl, Matrix pil) {
            this.pi = pi;
            this.pis = pis;
            this.pi0 = pi0;
            this.scc = scc;
            this.isrec = isrec;
            this.Pl = Pl;
            this.pil = pil;
        }
    }

    private static Map<String, Object> defaultOptions() {
        Map<String, Object> opts = new HashMap<String, Object>();
        opts.put("tol", 1e-12);
        return opts;
    }

    /**
     * Estimate limiting distribution for a DTMC that may have reducible components.
     *
     * @param P DTMC transition matrix
     * @param pin Initial vector (may be null)
     * @param options Solution options including tolerance
     * @return Pair of (limiting distribution, SCC information)
     */
    public static Pair<Matrix, List<List<Integer>>> dtmc_solve_reducible(
            Matrix P, Matrix pin, Map<String, Object> options) {
        Pair<int[], boolean[]> sccInfo = stronglyconncomp(P);
        int[] scc = sccInfo.getLeft();
        boolean[] isrec = sccInfo.getRight();

        int maxScc = -1;
        for (int v : scc) {
            if (v > maxScc) maxScc = v;
        }
        int numSCC = maxScc + 1;

        if (numSCC == 1) {
            Matrix piResult = Dtmc_solve.dtmc_solve(P);
            List<Integer> all = new ArrayList<Integer>();
            for (int i = 0; i < P.getNumRows(); i++) all.add(i);
            List<List<Integer>> sccList = new ArrayList<List<Integer>>();
            sccList.add(all);
            return new Pair<Matrix, List<List<Integer>>>(piResult, sccList);
        }

        // Create lumped transition matrix
        Matrix Pl = Matrix.zeros(numSCC, numSCC);
        List<List<Integer>> sccIdx = new ArrayList<List<Integer>>();
        for (int i = 0; i < numSCC; i++) {
            List<Integer> idx = new ArrayList<Integer>();
            for (int j = 0; j < P.getNumRows(); j++) {
                if (scc[j] == i) idx.add(j);
            }
            sccIdx.add(idx);
        }

        // Compute transition probabilities between SCCs
        for (int i = 0; i < numSCC; i++) {
            for (int j = 0; j < numSCC; j++) {
                if (i != j) {
                    double sum = 0.0;
                    for (Integer si : sccIdx.get(i)) {
                        for (Integer sj : sccIdx.get(j)) {
                            sum += P.get(si, sj);
                        }
                    }
                    Pl.set(i, j, sum);
                }
            }
        }

        // Make lumped matrix stochastic
        Pl = Dtmc_makestochastic.dtmc_makestochastic(Pl);

        // Ensure recurrent SCCs have self-loops
        for (int i = 0; i < numSCC; i++) {
            if (isrec[i]) {
                Pl.set(i, i, 1.0);
            }
        }

        // Compute initial distribution for lumped chain
        Matrix pinl;
        if (pin == null) {
            double[] temp = new double[numSCC];
            for (int i = 0; i < numSCC; i++) temp[i] = 1.0;
            // Check for zero columns (absorbing states)
            for (int j = 0; j < P.getNumCols(); j++) {
                double colSum = 0.0;
                for (int i = 0; i < P.getNumRows(); i++) {
                    colSum += P.get(i, j);
                }
                if (colSum < 1e-12) {
                    int sccIndex = scc[j];
                    temp[sccIndex] = 0.0;
                }
            }
            double sum = 0.0;
            for (double v : temp) sum += v;
            pinl = new Matrix(1, numSCC);
            for (int i = 0; i < numSCC; i++) {
                pinl.set(0, i, temp[i] / sum);
            }
        } else {
            double[] temp = new double[numSCC];
            for (int i = 0; i < numSCC; i++) {
                double sum = 0.0;
                for (Integer idx : sccIdx.get(i)) {
                    sum += pin.get(0, idx);
                }
                temp[i] = sum;
            }
            pinl = new Matrix(1, numSCC);
            for (int i = 0; i < numSCC; i++) {
                pinl.set(0, i, temp[i]);
            }
        }

        // Compute limiting distribution
        Matrix piResult = Matrix.zeros(1, P.getNumRows());
        List<Matrix> pis = new ArrayList<Matrix>();

        // Compute spectral decomposition of lumped matrix
        Matrix PI = computeLimitingMatrix(Pl);

        for (int i = 0; i < numSCC; i++) {
            if (pinl.get(0, i) > 0) {
                Matrix pi0i = Matrix.zeros(1, numSCC);
                pi0i.set(0, i, 1.0);
                Matrix pili = pi0i.mult(PI);

                Matrix pisi = Matrix.zeros(1, P.getNumRows());
                for (int j = 0; j < numSCC; j++) {
                    if (pili.get(0, j) > 0 && !sccIdx.get(j).isEmpty()) {
                        Matrix indices = new Matrix(sccIdx.get(j).size(), 1);
                        for (int k = 0; k < sccIdx.get(j).size(); k++) {
                            indices.set(k, 0, sccIdx.get(j).get(k).doubleValue());
                        }
                        Matrix subP = P.getSubMatrix(indices, indices);
                        Matrix subPi = Dtmc_solve.dtmc_solve(subP);
                        for (int k = 0; k < sccIdx.get(j).size(); k++) {
                            int stateIdx = sccIdx.get(j).get(k);
                            pisi.set(0, stateIdx, pili.get(0, j) * subPi.get(0, k));
                        }
                    }
                }
                pis.add(pisi);

                for (int k = 0; k < P.getNumRows(); k++) {
                    piResult.set(0, k, piResult.get(0, k) + pisi.get(0, k) * pinl.get(0, i));
                }
            }
        }

        // Handle single transient SCC case
        List<Integer> transientStates = new ArrayList<Integer>();
        for (int i = 0; i < isrec.length; i++) {
            if (!isrec[i]) transientStates.add(i);
        }
        if (transientStates.size() == 1 && pin == null) {
            return new Pair<Matrix, List<List<Integer>>>(pis.get(transientStates.get(0)), sccIdx);
        }

        return new Pair<Matrix, List<List<Integer>>>(piResult, sccIdx);
    }

    public static Pair<Matrix, List<List<Integer>>> dtmc_solve_reducible(Matrix P) {
        return dtmc_solve_reducible(P, null, defaultOptions());
    }

    public static Pair<Matrix, List<List<Integer>>> dtmc_solve_reducible(Matrix P, Matrix pin) {
        return dtmc_solve_reducible(P, pin, defaultOptions());
    }

    /**
     * Full version that returns all computed values.
     */
    public static DtmcSolveReducibleResult dtmc_solve_reducible_full(
            Matrix P, Matrix pin, Map<String, Object> options) {
        Pair<Matrix, List<List<Integer>>> r = dtmc_solve_reducible(P, pin, options);
        Matrix pi = r.getLeft();
        List<List<Integer>> scc = r.getRight();

        List<Matrix> pis = new ArrayList<Matrix>();
        pis.add(pi);

        List<Boolean> isrec = new ArrayList<Boolean>();
        for (int i = 0; i < scc.size(); i++) isrec.add(true);

        return new DtmcSolveReducibleResult(pi, pis, Matrix.zeros(1, 1), scc, isrec, P, pi);
    }

    public static DtmcSolveReducibleResult dtmc_solve_reducible_full(Matrix P) {
        return dtmc_solve_reducible_full(P, null, defaultOptions());
    }

    public static DtmcSolveReducibleResult dtmc_solve_reducible_full(Matrix P, Matrix pin) {
        return dtmc_solve_reducible_full(P, pin, defaultOptions());
    }

    private static Pair<int[], boolean[]> stronglyconncomp(Matrix P) {
        DirectedGraph graph = new DirectedGraph(P);
        DirectedGraph.SCCResult result = graph.stronglyconncomp();

        // Convert 1-based indexing to 0-based indexing for SCC assignments
        int[] scc = new int[result.I.length];
        for (int i = 0; i < result.I.length; i++) scc[i] = result.I[i] - 1;

        return new Pair<int[], boolean[]>(scc, result.recurrent);
    }

    private static Matrix computeLimitingMatrix(Matrix P) {
        // Pre-check: if eigenvector matrix is singular, use power method directly
        try {
            Ret.Eigs eigs = P.eigvec();
            Matrix V = eigs.vectors;
            try {
                V.inv();
            } catch (Exception e) {
                return computeLimitingMatrixPowerMethod(P);
            }
        } catch (Exception e) {
            return computeLimitingMatrixPowerMethod(P);
        }

        try {
            Ret.SpectralDecomposition spectral = Matrix.spectd(P);
            MatrixCell projectors = spectral.projectors;

            Matrix PI = Matrix.zeros(P.getNumRows(), P.getNumCols());
            Matrix spectrum = spectral.spectrum;

            for (int e = 0; e < spectrum.getNumRows(); e++) {
                if (Math.abs(spectrum.get(e, e) - 1.0) < 1e-12) {
                    Matrix proj = (Matrix) projectors.get(e);
                    PI = PI.add(proj);
                }
            }

            if (PI.hasNaN()) {
                return computeLimitingMatrixPowerMethod(P);
            }

            return PI;
        } catch (Exception e) {
            return computeLimitingMatrixPowerMethod(P);
        }
    }

    private static Matrix computeLimitingMatrixPowerMethod(Matrix P) {
        int n = P.getNumRows();
        Matrix Pk = new Matrix(P);
        int maxIter = 1000;
        double tol = 1e-10;

        for (int iter = 0; iter < maxIter; iter++) {
            Matrix Pk1 = Pk.mult(P);

            double maxDiff = 0.0;
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    double diff = Math.abs(Pk1.get(i, j) - Pk.get(i, j));
                    if (diff > maxDiff) maxDiff = diff;
                }
            }

            Pk = Pk1;

            if (maxDiff < tol) {
                break;
            }
        }

        return Pk;
    }
}
