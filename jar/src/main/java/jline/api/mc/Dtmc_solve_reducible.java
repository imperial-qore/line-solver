package jline.api.mc;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import jline.io.Ret;
import jline.util.Pair;
import jline.util.graph.DirectedGraph;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;
import jline.util.matrix.MatrixEntry;

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

        // Compute transition probabilities between SCCs, over the nonzeros of P: the pairwise
        // sum over SCC members read every (state, state) entry, O(n^2) on a sparse P
        double[][] plSum = new double[numSCC][];
        for (Iterator<MatrixEntry> it = P.nonZeroIterator(); it.hasNext(); ) {
            MatrixEntry e = it.next();
            int i = scc[e.row], j = scc[e.col];
            if (i != j) {
                if (plSum[i] == null) plSum[i] = new double[numSCC];
                plSum[i][j] += e.value;
            }
        }
        for (int i = 0; i < numSCC; i++) {
            if (plSum[i] == null) continue;
            for (int j = 0; j < numSCC; j++) {
                if (plSum[i][j] != 0) Pl.set(i, j, plSum[i][j]);
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
            double[] colSums = new double[P.getNumCols()];
            for (Iterator<MatrixEntry> it = P.nonZeroIterator(); it.hasNext(); ) {
                MatrixEntry e = it.next();
                colSums[e.col] += e.value;
            }
            for (int j = 0; j < P.getNumCols(); j++) {
                double colSum = colSums[j];
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

        // One row per SCC, so that the position of a row in `pis` IS its SCC index.
        // The rows used to be appended only when pinl(i) > 0, which compacted the list
        // and silently desynchronized it from the SCC numbering; the "single transient
        // SCC" branch below then indexed `pis` with an SCC id and threw
        // IndexOutOfBounds as soon as a skipped component preceded the transient one.
        // Each SCC's internal stationary vector depends only on that SCC, not on where
        // the chain started, so solve it ONCE per SCC and reuse it across the rows
        // below. Solving it inside the (i, j) loop costs numSCC^2 solves for numSCC
        // distinct answers.
        Matrix[] sccPi = new Matrix[numSCC];
        for (int j = 0; j < numSCC; j++) {
            if (!sccIdx.get(j).isEmpty()) {
                Matrix indices = new Matrix(sccIdx.get(j).size(), 1);
                for (int k = 0; k < sccIdx.get(j).size(); k++) {
                    indices.set(k, 0, sccIdx.get(j).get(k).doubleValue());
                }
                sccPi[j] = Dtmc_solve.dtmc_solve(P.getSubMatrix(indices, indices));
            }
        }

        for (int i = 0; i < numSCC; i++) {
            // e_i * PI, read as the row it is rather than as numSCC row-by-matrix products
            Matrix pili = PI.getRow(i);

            Matrix pisi = Matrix.zeros(1, P.getNumRows());
            for (int j = 0; j < numSCC; j++) {
                if (pili.get(0, j) > 0 && sccPi[j] != null) {
                    for (int k = 0; k < sccIdx.get(j).size(); k++) {
                        int stateIdx = sccIdx.get(j).get(k);
                        pisi.set(0, stateIdx, pili.get(0, j) * sccPi[j].get(0, k));
                    }
                }
            }
            pis.add(pisi);

            if (pinl.get(0, i) > 0) {
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
        Matrix PI = computeLimitingMatrixAcyclic(P);
        if (PI != null) {
            return PI;
        }
        return computeLimitingMatrixSpectral(P);
    }

    /**
     * lim P^n of a lumped SCC chain, exactly, by absorption probabilities.
     *
     * <p>The lumped chain is a condensation: its recurrent states are absorbing
     * and the others form a DAG, so P is triangular up to a permutation with
     * eigenvalues 0 and 1 only. Whenever a transient lump leads to another the
     * eigenvector matrix is defective, which is where the spectral route is
     * both slow (a dense eigendecomposition, hours at n = 6131) and abandoned
     * by the reference for the power method. The limit is instead the
     * absorption matrix: PI(a,a) = 1 for absorbing a, and PI(i,:) =
     * sum_{j != i} P(i,j) PI(j,:) / (1 - P(i,i)) for a transient i, solved in
     * reverse topological order. The division keeps it exact under the
     * round-off self-loop dtmc_makestochastic may leave on a transient row.
     *
     * @return the limiting matrix, or null if P is not of that form
     */
    private static Matrix computeLimitingMatrixAcyclic(Matrix P) {
        int n = P.getNumRows();
        List<List<Integer>> succ = new ArrayList<List<Integer>>();
        List<List<Double>> prob = new ArrayList<List<Double>>();
        double[] diag = new double[n];
        for (int i = 0; i < n; i++) {
            succ.add(new ArrayList<Integer>());
            prob.add(new ArrayList<Double>());
        }
        for (Iterator<MatrixEntry> it = P.nonZeroIterator(); it.hasNext(); ) {
            MatrixEntry e = it.next();
            if (e.value == 0) continue;
            if (e.row == e.col) {
                diag[e.row] = e.value;
            } else {
                succ.get(e.row).add(e.col);
                prob.get(e.row).add(e.value);
            }
        }
        boolean[] absorbing = new boolean[n];
        for (int i = 0; i < n; i++) {
            if (succ.get(i).isEmpty()) {
                if (Math.abs(diag[i] - 1.0) > 1e-12) return null;
                absorbing[i] = true;
            } else if (diag[i] >= 1.0 - 1e-12) {
                return null;   // a self-loop of 1 with an exit is not a stochastic row
            }
        }
        // Iterative DFS postorder: a row is solved after all of its successors
        int[] mark = new int[n];   // 0 unvisited, 1 on the stack, 2 done
        List<Map<Integer, Double>> row = new ArrayList<Map<Integer, Double>>();
        for (int i = 0; i < n; i++) row.add(null);
        int[] stack = new int[n];
        int[] cursor = new int[n];
        for (int s = 0; s < n; s++) {
            if (mark[s] != 0) continue;
            int top = 0;
            stack[top] = s;
            cursor[s] = 0;
            mark[s] = 1;
            while (top >= 0) {
                int u = stack[top];
                if (cursor[u] < succ.get(u).size()) {
                    int v = succ.get(u).get(cursor[u]++);
                    if (mark[v] == 1) return null;   // a cycle: not a condensation
                    if (mark[v] == 0) {
                        mark[v] = 1;
                        cursor[v] = 0;
                        stack[++top] = v;
                    }
                    continue;
                }
                Map<Integer, Double> r = new HashMap<Integer, Double>();
                if (absorbing[u]) {
                    r.put(u, 1.0);
                } else {
                    double scale = 1.0 / (1.0 - diag[u]);
                    for (int k = 0; k < succ.get(u).size(); k++) {
                        double p = prob.get(u).get(k) * scale;
                        for (Map.Entry<Integer, Double> e : row.get(succ.get(u).get(k)).entrySet()) {
                            Double old = r.get(e.getKey());
                            r.put(e.getKey(), (old == null ? 0.0 : old) + p * e.getValue());
                        }
                    }
                }
                row.set(u, r);
                mark[u] = 2;
                top--;
            }
        }
        int nnz = 0;
        for (int i = 0; i < n; i++) nnz += row.get(i).size();
        Matrix PI = new Matrix(n, n, nnz);
        for (int i = 0; i < n; i++) {
            for (Map.Entry<Integer, Double> e : row.get(i).entrySet()) {
                if (e.getValue() != 0) PI.set(i, e.getKey(), e.getValue());
            }
        }
        return PI;
    }

    private static Matrix computeLimitingMatrixSpectral(Matrix P) {
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
