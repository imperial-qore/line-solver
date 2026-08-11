package jline.api.mc;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import jline.util.Pair;
import jline.util.graph.DirectedGraph;
import jline.util.matrix.Matrix;

/**
 * Solve reducible CTMCs via direct block decomposition on the generator matrix.
 *
 * Algorithm (based on SMART's computeInfinityDistribution):
 *   1. Decompose states into transient and recurrent classes via SCC detection
 *   2. For transient states: solve n * Q_tt = -p0_t for expected sojourn
 *   3. Compute hitting probabilities: h = n * Q_ta + p0_r
 *   4. For each recurrent class: solve pi_c * Q_cc = 0, scale by hitting prob
 */
public final class Ctmc_solve_reducible_blkdecomp {
    private Ctmc_solve_reducible_blkdecomp() {}

    private static Map<String, Object> defaultOptions() {
        Map<String, Object> opts = new HashMap<String, Object>();
        opts.put("tol", 1e-12);
        return opts;
    }

    public static Pair<Matrix, List<List<Integer>>> ctmc_solve_reducible_blkdecomp(Matrix Q) {
        return ctmc_solve_reducible_blkdecomp(Q, null, defaultOptions());
    }

    public static Pair<Matrix, List<List<Integer>>> ctmc_solve_reducible_blkdecomp(Matrix Q, Matrix pi0) {
        return ctmc_solve_reducible_blkdecomp(Q, pi0, defaultOptions());
    }

    public static Pair<Matrix, List<List<Integer>>> ctmc_solve_reducible_blkdecomp(Matrix Q, Matrix pi0,
                                                                                   Map<String, Object> options) {
        Matrix Qmat = Ctmc_makeinfgen.ctmc_makeinfgen(Q.copy());
        int N = Qmat.getNumRows();

        // Build adjacency from off-diagonal positive entries
        Matrix Adj = new Matrix(N, N);
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                if (i != j && Qmat.get(i, j) > 0) Adj.set(i, j, 1.0);
            }
        }

        DirectedGraph graph = new DirectedGraph(Adj);
        DirectedGraph.SCCResult sccResult = graph.stronglyconncomp();
        int[] scc = new int[N];
        for (int i = 0; i < N; i++) {
            scc[i] = sccResult.I[i] - 1;
        }
        boolean[] isrec = sccResult.recurrent;
        int numSCC = isrec.length;

        // Irreducible case
        if (numSCC == 1) {
            Matrix pi = Ctmc_solve.ctmc_solve(Qmat);
            List<Integer> all = new ArrayList<Integer>();
            for (int i = 0; i < N; i++) all.add(i);
            List<List<Integer>> result = new ArrayList<List<Integer>>();
            result.add(all);
            return new Pair<Matrix, List<List<Integer>>>(pi, result);
        }

        // Build SCC index sets
        @SuppressWarnings("unchecked")
        List<Integer>[] sccIdx = (List<Integer>[]) new List[numSCC];
        for (int i = 0; i < numSCC; i++) {
            List<Integer> list = new ArrayList<Integer>();
            for (int s = 0; s < N; s++) {
                if (scc[s] == i) list.add(s);
            }
            Collections.sort(list);
            sccIdx[i] = list;
        }

        // Classify SCCs
        List<Integer> transSccIds = new ArrayList<Integer>();
        List<Integer> recSccIds = new ArrayList<Integer>();
        for (int i = 0; i < numSCC; i++) {
            if (!isrec[i]) transSccIds.add(i);
            else recSccIds.add(i);
        }

        // Gather ordered state indices
        List<Integer> transStates = new ArrayList<Integer>();
        for (int i : transSccIds) transStates.addAll(sccIdx[i]);
        Collections.sort(transStates);
        List<Integer> recStates = new ArrayList<Integer>();
        for (int i : recSccIds) recStates.addAll(sccIdx[i]);
        Collections.sort(recStates);
        int nt = transStates.size();
        int nr = recStates.size();

        // Extract Q sub-blocks
        Matrix Q_tt = null;
        Matrix Q_ta = null;
        if (nt > 0 && nr > 0) {
            Matrix transIdx = new Matrix(nt, 1);
            for (int i = 0; i < nt; i++) transIdx.set(i, 0, (double) transStates.get(i));
            Matrix recIdx = new Matrix(nr, 1);
            for (int i = 0; i < nr; i++) recIdx.set(i, 0, (double) recStates.get(i));
            Q_tt = Qmat.getSubMatrix(transIdx, transIdx);
            Q_ta = Qmat.getSubMatrix(transIdx, recIdx);
        }

        // Compute per-SCC limiting distributions
        Matrix[] pis = new Matrix[numSCC];
        for (int i = 0; i < numSCC; i++) pis[i] = Matrix.zeros(1, N);

        for (int s = 0; s < numSCC; s++) {
            List<Integer> classStates = sccIdx[s];
            int classSize = classStates.size();

            double[] hit = new double[nr];

            if (nt > 0 && Q_tt != null && Q_ta != null) {
                double[] p0_t = new double[nt];
                for (int idx = 0; idx < transStates.size(); idx++) {
                    int state = transStates.get(idx);
                    if (classStates.contains(state)) {
                        p0_t[idx] = 1.0 / classSize;
                    }
                }

                boolean anyPositive = false;
                for (double v : p0_t) {
                    if (v > 0.0) { anyPositive = true; break; }
                }

                if (anyPositive) {
                    Matrix negP0Col = new Matrix(nt, 1);
                    for (int i = 0; i < nt; i++) negP0Col.set(i, 0, -p0_t[i]);

                    Matrix sojournCol = new Matrix(nt, 1);
                    // Above the dispatch threshold the transient block is what the direct
                    // factorization cannot hold; the direct solve stays the fallback.
                    boolean solved = false;
                    if (nt > Ctmc_solve.GMRES_MIN_STATES) {
                        Ctmc_gmres.GmresResult g =
                                Ctmc_gmres.ctmc_gmres(Q_tt.transpose(), negP0Col, 0.0, 0, 0, null);
                        if (g.flag == 0) {
                            sojournCol = g.x;
                            solved = true;
                        }
                    }
                    if (!solved) {
                        Matrix.solveDirect(Q_tt.transpose(), negP0Col, sojournCol);
                    }

                    Matrix sojournRow = sojournCol.transpose();
                    Matrix hitMatrix = sojournRow.mult(Q_ta);
                    for (int i = 0; i < nr; i++) hit[i] = hitMatrix.get(0, i);
                }
            }

            for (int idx = 0; idx < recStates.size(); idx++) {
                int state = recStates.get(idx);
                if (classStates.contains(state)) {
                    hit[idx] += 1.0 / classSize;
                }
            }

            for (int c : recSccIds) {
                List<Integer> idxC = sccIdx[c];

                double reachprob = 0.0;
                for (int state : idxC) {
                    int loc = recStates.indexOf(state);
                    if (loc >= 0) reachprob += hit[loc];
                }

                if (reachprob < 1e-15) continue;

                if (idxC.size() == 1) {
                    pis[s].set(0, idxC.get(0), reachprob);
                } else {
                    Matrix indices = new Matrix(idxC.size(), 1);
                    for (int i = 0; i < idxC.size(); i++) indices.set(i, 0, (double) idxC.get(i));
                    Matrix Qcc = Qmat.getSubMatrix(indices, indices);
                    Matrix piC = Ctmc_solve.ctmc_solve(Qcc);
                    for (int k = 0; k < idxC.size(); k++) {
                        int stateIdx = idxC.get(k);
                        pis[s].set(0, stateIdx, piC.get(0, k) * reachprob);
                    }
                }
            }
        }

        // Compute initial SCC probabilities for weighted average
        double[] pinl = new double[numSCC];
        if (pi0 == null) {
            for (int i = 0; i < numSCC; i++) pinl[i] = 1.0;
            for (int j = 0; j < N; j++) {
                double colSum = 0.0;
                for (int i = 0; i < N; i++) colSum += Math.abs(Qmat.get(i, j));
                if (colSum < 1e-12) pinl[scc[j]] = 0.0;
            }
            double totalPinl = 0.0;
            for (double v : pinl) totalPinl += v;
            if (totalPinl > 0) {
                for (int i = 0; i < numSCC; i++) pinl[i] /= totalPinl;
            } else {
                for (int i = 0; i < numSCC; i++) pinl[i] = 1.0 / numSCC;
            }
        } else {
            for (int i = 0; i < numSCC; i++) {
                double sum = 0.0;
                for (int idx : sccIdx[i]) sum += pi0.get(0, idx);
                pinl[i] = sum;
            }
        }

        // Weighted average over starting SCCs
        Matrix pi = Matrix.zeros(1, N);
        for (int i = 0; i < numSCC; i++) {
            if (pinl[i] > 0) {
                for (int k = 0; k < N; k++) {
                    pi.set(0, k, pi.get(0, k) + pis[i].get(0, k) * pinl[i]);
                }
            }
        }

        // Special case: single transient SCC without explicit initial distribution
        if (transSccIds.size() == 1 && pi0 == null) {
            for (int k = 0; k < N; k++) {
                pi.set(0, k, pis[transSccIds.get(0)].get(0, k));
            }
        }

        // Normalize
        double total = 0.0;
        for (int k = 0; k < N; k++) total += pi.get(0, k);
        if (total > 0) {
            for (int k = 0; k < N; k++) pi.set(0, k, pi.get(0, k) / total);
        }

        List<List<Integer>> sccLists = new ArrayList<List<Integer>>();
        for (List<Integer> l : sccIdx) {
            sccLists.add(new ArrayList<Integer>(l));
        }
        return new Pair<Matrix, List<List<Integer>>>(pi, sccLists);
    }

    public static Ctmc_solve_reducible_blkdecompResult ctmc_solve_reducible_blkdecomp_full(Matrix Q) {
        return ctmc_solve_reducible_blkdecomp_full(Q, null, defaultOptions());
    }

    public static Ctmc_solve_reducible_blkdecompResult ctmc_solve_reducible_blkdecomp_full(Matrix Q, Matrix pi0) {
        return ctmc_solve_reducible_blkdecomp_full(Q, pi0, defaultOptions());
    }

    public static Ctmc_solve_reducible_blkdecompResult ctmc_solve_reducible_blkdecomp_full(Matrix Q, Matrix pi0,
                                                                                          Map<String, Object> options) {
        Matrix Qmat = Ctmc_makeinfgen.ctmc_makeinfgen(Q.copy());
        int N = Qmat.getNumRows();

        Matrix Adj = new Matrix(N, N);
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                if (i != j && Qmat.get(i, j) > 0) Adj.set(i, j, 1.0);
            }
        }

        DirectedGraph graph = new DirectedGraph(Adj);
        DirectedGraph.SCCResult sccResult = graph.stronglyconncomp();
        int[] scc = new int[N];
        for (int i = 0; i < N; i++) scc[i] = sccResult.I[i] - 1;
        boolean[] isrec = sccResult.recurrent;
        int numSCC = isrec.length;

        Pair<Matrix, List<List<Integer>>> baseResult = ctmc_solve_reducible_blkdecomp(Q, pi0, options);
        Matrix pi = baseResult.getLeft();
        List<List<Integer>> sccLists = baseResult.getRight();

        @SuppressWarnings("unchecked")
        List<Integer>[] sccIdx = (List<Integer>[]) new List[numSCC];
        for (int i = 0; i < numSCC; i++) {
            List<Integer> list = new ArrayList<Integer>();
            for (int s = 0; s < N; s++) {
                if (scc[s] == i) list.add(s);
            }
            Collections.sort(list);
            sccIdx[i] = list;
        }
        List<Matrix> pisList = new ArrayList<Matrix>();
        for (int s = 0; s < numSCC; s++) {
            Matrix p0s = Matrix.zeros(1, N);
            List<Integer> classStates = sccIdx[s];
            for (int state : classStates) {
                p0s.set(0, state, 1.0 / classStates.size());
            }
            Pair<Matrix, List<List<Integer>>> stepResult = ctmc_solve_reducible_blkdecomp(Q, p0s, options);
            pisList.add(stepResult.getLeft());
        }

        List<Boolean> isrecList = new ArrayList<Boolean>();
        for (boolean b : isrec) isrecList.add(b);

        return new Ctmc_solve_reducible_blkdecompResult(
                pi,
                pisList,
                pi0 != null ? pi0 : Matrix.zeros(1, N),
                sccLists,
                isrecList
        );
    }
}
