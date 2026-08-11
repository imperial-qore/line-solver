package jline.api.mc;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import jline.util.Pair;
import jline.util.matrix.Matrix;

public final class Ctmc_solve_reducible {
    private Ctmc_solve_reducible() {}

    /**
     * Result of full reducible CTMC solve.
     */
    public static final class CtmcSolveReducibleResult {
        public final Matrix pi;
        public final List<Matrix> pis;
        public final Matrix pi0;
        public final List<List<Integer>> scc;
        public final List<Boolean> isrec;

        public CtmcSolveReducibleResult(Matrix pi, List<Matrix> pis, Matrix pi0,
                                        List<List<Integer>> scc, List<Boolean> isrec) {
            this.pi = pi;
            this.pis = pis;
            this.pi0 = pi0;
            this.scc = scc;
            this.isrec = isrec;
        }
    }

    private static Map<String, Object> defaultOptions() {
        Map<String, Object> opts = new HashMap<String, Object>();
        opts.put("tol", 1e-12);
        return opts;
    }

    /**
     * Solve reducible CTMCs by converting to DTMC via randomization.
     *
     * @param Q Infinitesimal generator matrix (possibly reducible)
     * @param pi0 Initial state distribution (may be null)
     * @param options Solution options
     * @return Pair of (pi: steady-state distribution, scc: strongly connected components info)
     */
    public static Pair<Matrix, List<List<Integer>>> ctmc_solve_reducible(
            Matrix Q, Matrix pi0, Map<String, Object> options) {
        // Convert CTMC to DTMC via randomization
        Pair<Matrix, Double> rand = Ctmc_randomization.ctmc_randomization(Q);
        Matrix P = rand.getLeft();

        // Solve the reducible DTMC
        Pair<Matrix, List<List<Integer>>> dtmcResult =
                Dtmc_solve_reducible.dtmc_solve_reducible(P, pi0, options);
        Matrix pi = dtmcResult.getLeft();
        List<List<Integer>> scc = dtmcResult.getRight();

        return new Pair<Matrix, List<List<Integer>>>(pi, scc);
    }

    public static Pair<Matrix, List<List<Integer>>> ctmc_solve_reducible(Matrix Q) {
        return ctmc_solve_reducible(Q, null, defaultOptions());
    }

    public static Pair<Matrix, List<List<Integer>>> ctmc_solve_reducible(Matrix Q, Matrix pi0) {
        return ctmc_solve_reducible(Q, pi0, defaultOptions());
    }

    /**
     * Alternative signature that returns additional information.
     */
    public static CtmcSolveReducibleResult ctmc_solve_reducible_full(
            Matrix Q, Matrix pi0, Map<String, Object> options) {
        // Convert CTMC to DTMC via randomization
        Pair<Matrix, Double> rand = Ctmc_randomization.ctmc_randomization(Q);
        Matrix P = rand.getLeft();

        // Solve the reducible DTMC (full implementation)
        Dtmc_solve_reducible.DtmcSolveReducibleResult result =
                Dtmc_solve_reducible.dtmc_solve_reducible_full(P, pi0, options);

        return new CtmcSolveReducibleResult(result.pi, result.pis, result.pi0, result.scc, result.isrec);
    }

    public static CtmcSolveReducibleResult ctmc_solve_reducible_full(Matrix Q) {
        return ctmc_solve_reducible_full(Q, null, defaultOptions());
    }

    public static CtmcSolveReducibleResult ctmc_solve_reducible_full(Matrix Q, Matrix pi0) {
        return ctmc_solve_reducible_full(Q, pi0, defaultOptions());
    }
}
