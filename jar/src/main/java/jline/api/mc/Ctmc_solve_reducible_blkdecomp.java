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
 * Algorithm:
 *   1. Decompose states into transient and recurrent classes via SCC detection
 *   2. For transient states: solve sojourn * Q_tt = -p0_t for expected sojourn
 *   3. Compute hitting probabilities: hit = sojourn * Q_ta + p0_r
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

    public static Pair<Matrix, List<List<Integer>>> ctmc_solve_reducible_blkdecomp(Matrix Q, Matrix pin) {
        return ctmc_solve_reducible_blkdecomp(Q, pin, defaultOptions());
    }

    public static Pair<Matrix, List<List<Integer>>> ctmc_solve_reducible_blkdecomp(Matrix Q, Matrix pin,
                                                                                   Map<String, Object> options) {
        Matrix Qmat = Ctmc_makeinfgen.ctmc_makeinfgen(Q.copy());
        int N = Qmat.getNumRows();

        // Build adjacency from off-diagonal positive entries
        Matrix Adj = new Matrix(N, N);
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                // arc by MAGNITUDE, never by sign: an ME generator embeds genuinely negative off-diagonals -- see _kb/11-conventions-and-gotchas.md
                if (i != j && Math.abs(Qmat.get(i, j)) > jline.GlobalConstants.ArcTol) Adj.set(i, j, 1.0);
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

        // Compute per-SCC limiting distributions, each from a uniform start in its SCC
        Matrix[] pis = new Matrix[numSCC];
        for (int s = 0; s < numSCC; s++) {
            double[] p0 = new double[N];
            List<Integer> classStates = sccIdx[s];
            for (int k = 0; k < classStates.size(); k++) {
                p0[classStates.get(k)] = 1.0 / classStates.size();
            }
            pis[s] = absorbLimiting(p0, Qmat, N, nt, nr, transStates, recStates, sccIdx, recSccIds, Q_tt, Q_ta);
        }

        Matrix pi;
        if (pin == null) {
            // No initial vector: mix the uniform-start rows, weighting the SCCs equally.
            // An SCC holding a state whose column of Q is entirely zero is EXCLUDED: such a
            // state has no transition in and none out, so it is isolated and carries no
            // dynamics to start from. (Not the same as absorbing, which has incoming
            // transitions and a zero off-diagonal ROW.)
            double[] pinl = new double[numSCC];
            for (int i = 0; i < numSCC; i++) pinl[i] = 1.0;
            for (int j = 0; j < N; j++) {
                double colSum = 0.0;
                for (int i = 0; i < N; i++) colSum += Math.abs(Qmat.get(i, j));
                if (colSum < 1e-12) pinl[scc[j]] = 0.0;
            }
            double totalPinl = 0.0;
            for (int i = 0; i < numSCC; i++) totalPinl += pinl[i];
            if (totalPinl > 0) {
                for (int i = 0; i < numSCC; i++) pinl[i] /= totalPinl;
            } else {
                for (int i = 0; i < numSCC; i++) pinl[i] = 1.0 / numSCC;
            }

            pi = Matrix.zeros(1, N);
            for (int i = 0; i < numSCC; i++) {
                if (pinl[i] > 0) {
                    for (int k = 0; k < N; k++) {
                        pi.set(0, k, pi.get(0, k) + pis[i].get(0, k) * pinl[i]);
                    }
                }
            }

            // Special case: single transient SCC without explicit initial distribution
            if (transSccIds.size() == 1) {
                pi = pis[transSccIds.get(0)].copy();
            }
        } else {
            // An initial vector is available, so the exact absorption probabilities can be
            // computed from it directly. Lumping pin onto its SCCs and mixing the pis rows
            // would instead assume a uniform start within each SCC, which differs from the
            // truth whenever pin puts mass on a transient SCC holding more than one state
            // (states of the same transient SCC reach the recurrent classes with different
            // probabilities).
            double[] p0 = new double[N];
            double total0 = 0.0;
            for (int k = 0; k < N; k++) {
                p0[k] = pin.get(0, k);
                total0 += p0[k];
            }
            if (total0 > 0) {
                for (int k = 0; k < N; k++) p0[k] /= total0;
            }
            pi = absorbLimiting(p0, Qmat, N, nt, nr, transStates, recStates, sccIdx, recSccIds, Q_tt, Q_ta);
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


    /**
     * Limiting distribution reached from the initial vector p0: the mass absorbed in each
     * recurrent class (BSCC) redistributed over that class according to its own stationary
     * vector. Transient states receive zero.
     */
    private static Matrix absorbLimiting(double[] p0, Matrix Qmat, int N, int nt, int nr,
                                         List<Integer> transStates, List<Integer> recStates,
                                         List<Integer>[] sccIdx, List<Integer> recSccIds,
                                         Matrix Q_tt, Matrix Q_ta) {
        Matrix piv = Matrix.zeros(1, N);

        // Absorption probabilities into the recurrent states
        double[] hit = new double[nr];
        if (nt > 0 && Q_tt != null && Q_ta != null) {
            double[] p0_t = new double[nt];
            boolean anyPositive = false;
            for (int idx = 0; idx < nt; idx++) {
                p0_t[idx] = p0[transStates.get(idx)];
                if (p0_t[idx] > 0.0) anyPositive = true;
            }
            if (anyPositive) {
                // Solve sojourn * Q_tt = -p0_t for expected sojourn in transient states.
                // Q_tt is non-singular (Hurwitz) for transient states.
                Matrix negP0Col = new Matrix(nt, 1);
                for (int i = 0; i < nt; i++) negP0Col.set(i, 0, -p0_t[i]);
                Matrix sojourn = new Matrix(nt, 1);
                // above the dispatch threshold the transient block is what the direct factorization cannot hold; the direct solve stays the fallback
                boolean solved = false;
                if (nt > Ctmc_solve.GMRES_MIN_STATES) {
                    Ctmc_gmres.GmresResult g =
                            Ctmc_gmres.ctmc_gmres(Q_tt.transpose(), negP0Col, 0.0, 0, 0, null);
                    if (g.flag == 0) {
                        sojourn = g.x;
                        solved = true;
                    } else {
                        // Short-recurrence retry before the cubic factorization, as in Ctmc_solve.
                        Ctmc_bicgstab.BicgstabResult bs =
                                Ctmc_bicgstab.ctmc_bicgstab(Q_tt.transpose(), negP0Col, 0.0, 0, null);
                        if (bs.flag == 0) {
                            sojourn = bs.x;
                            solved = true;
                        }
                    }
                }
                if (!solved) {
                    Matrix.solveDirect(Q_tt.transpose(), negP0Col, sojourn);
                }
                Matrix hitMatrix = sojourn.transpose().mult(Q_ta);
                for (int i = 0; i < nr; i++) hit[i] = hitMatrix.get(0, i);
            }
        }

        // Add initial mass already in recurrent states
        for (int idx = 0; idx < nr; idx++) {
            hit[idx] += p0[recStates.get(idx)];
        }

        // Solve steady state per recurrent class, scaled by hitting probability
        for (int ci = 0; ci < recSccIds.size(); ci++) {
            List<Integer> idx_c = sccIdx[recSccIds.get(ci)];
            double reachprob = 0.0;
            for (int k = 0; k < idx_c.size(); k++) {
                int loc = recStates.indexOf(idx_c.get(k));
                if (loc >= 0) reachprob += hit[loc];
            }
            if (reachprob < 1e-15) continue;

            if (idx_c.size() == 1) {
                // Absorbing state: hitting probability IS the final probability
                piv.set(0, idx_c.get(0), reachprob);
            } else {
                Matrix indices = new Matrix(idx_c.size(), 1);
                for (int i = 0; i < idx_c.size(); i++) indices.set(i, 0, (double) idx_c.get(i));
                Matrix Q_cc = Qmat.getSubMatrix(indices, indices);
                Matrix pi_c = Ctmc_solve.ctmc_solve(Q_cc);
                for (int k = 0; k < idx_c.size(); k++) {
                    piv.set(0, idx_c.get(k), pi_c.get(0, k) * reachprob);
                }
            }
        }

        return piv;
    }

    public static Ctmc_solve_reducible_blkdecompResult ctmc_solve_reducible_blkdecomp_full(Matrix Q) {
        return ctmc_solve_reducible_blkdecomp_full(Q, null, defaultOptions());
    }

    public static Ctmc_solve_reducible_blkdecompResult ctmc_solve_reducible_blkdecomp_full(Matrix Q, Matrix pin) {
        return ctmc_solve_reducible_blkdecomp_full(Q, pin, defaultOptions());
    }

    public static Ctmc_solve_reducible_blkdecompResult ctmc_solve_reducible_blkdecomp_full(Matrix Q, Matrix pin,
                                                                                          Map<String, Object> options) {
        Matrix Qmat = Ctmc_makeinfgen.ctmc_makeinfgen(Q.copy());
        int N = Qmat.getNumRows();

        Matrix Adj = new Matrix(N, N);
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                // arc by MAGNITUDE, never by sign: an ME generator embeds genuinely negative off-diagonals -- see _kb/11-conventions-and-gotchas.md
                if (i != j && Math.abs(Qmat.get(i, j)) > jline.GlobalConstants.ArcTol) Adj.set(i, j, 1.0);
            }
        }

        DirectedGraph graph = new DirectedGraph(Adj);
        DirectedGraph.SCCResult sccResult = graph.stronglyconncomp();
        int[] scc = new int[N];
        for (int i = 0; i < N; i++) scc[i] = sccResult.I[i] - 1;
        boolean[] isrec = sccResult.recurrent;
        int numSCC = isrec.length;

        Pair<Matrix, List<List<Integer>>> baseResult = ctmc_solve_reducible_blkdecomp(Q, pin, options);
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
                pin != null ? pin : Matrix.zeros(1, N),
                sccLists,
                isrecList
        );
    }
}
