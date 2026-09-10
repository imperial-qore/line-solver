/**
 * @file Continuous-time Markov chain steady-state solver
 *
 * Computes the steady-state probability distribution for CTMCs by solving the linear
 * system pi * Q = 0 with normalization pi * 1 = 1. Handles reducible chains by decomposing
 * into strongly connected components and solving each component separately.
 *
 * @since LINE 3.0
 */
package jline.api.mc;

import java.util.Iterator;
import java.util.Set;

import jline.GlobalConstants;
import jline.solvers.SolverOptions;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixEntry;

public final class Ctmc_solve {
    private Ctmc_solve() {}

    /**
     * Order above which the direct sparse factorization is abandoned in favour of GMRES.
     * Below it the direct solve is both faster and exact, so the iterative path buys
     * nothing; above it the factorization fill-in is what limits the tractable model.
     */
    public static final int GMRES_MIN_STATES = 6000;

    /**
     * Column mass below which a state counts as ISOLATED and is eliminated. The mass is
     * a computed sum, so testing it against an exact zero makes the answer depend on the
     * order the terms were accumulated -- i.e. on the BLAS kernel, i.e. on the CPU. Same
     * value native python uses in {@code api/mc/ctmc.py}.
     */
    private static final double ISOLATED_TOL = 1e-12;

    /**
     * Return the steady-state probability of a CTMC, choosing the solution method by
     * size alone.
     *
     * @param Q Infinitesimal generator of the CTMC
     * @return Steady-state probability vector
     */
    public static Matrix ctmc_solve(Matrix Q) {
        return ctmc_solve(Q, null);
    }

    /**
     * Return the steady-state probability of a CTMC.
     *
     * @param Q       Infinitesimal generator of the CTMC
     * @param options Solver options; {@code method} may force "gmres" or "direct", and
     *                {@code config.gmres_restart} sets the Krylov subspace dimension.
     *                May be null, in which case the defaults apply.
     * @return Steady-state probability vector
     */
    public static Matrix ctmc_solve(Matrix Q, SolverOptions options) {
        Matrix Qmat = Q.copy();
        if (Qmat.length() == 1) {
            Matrix p = new Matrix(1, 1, 1);
            p.set(0, 0, 1);
            return p;
        }

        Qmat = Ctmc_makeinfgen.ctmc_makeinfgen(Qmat);
        int n = Qmat.length();

        // B = abs(Q+Q')>0
        Matrix B = Qmat.add(1.0, Qmat.transpose());
        B.absEq();
        int[] colIndexesB = B.getColIndexes();
        double[] nzvB = B.getNonZeroValues();
        for (int colIdx = 0; colIdx < B.getNumCols(); colIdx++) {
            int col1 = colIndexesB[colIdx];
            int col2 = colIndexesB[colIdx + 1];

            for (int i = col1; i < col2; i++) {
                if (nzvB[i] > GlobalConstants.Zero) nzvB[i] = 1.0;
            }
        }

        // [nConnComp, connComp] = weaklyconncomp(B);
        Set<Set<Integer>> sets = Matrix.weaklyConnect(B, null);
        if (sets.size() > 1) {
            Matrix p = new Matrix(1, n);

            for (Set<Integer> set_c : sets) {
                // Qc = Q(connComp==c,connComp==c);
                Matrix Qc = new Matrix(set_c.size(), set_c.size());
                int Qc_row = 0;
                int Qc_col = 0;
                for (int q_row : set_c) {
                    for (int q_col : set_c) {
                        Qc.set(Qc_row, Qc_col++, Qmat.get(q_row, q_col));
                    }
                    Qc_row++;
                    Qc_col = 0;
                }
                // Qc = ctmc_makeinfgen(Qc);
                Qc = Ctmc_makeinfgen.ctmc_makeinfgen(Qc);
                // p(connComp==c) = ctmc_solve(Qc);
                Matrix ctmc_solve_Qc = ctmc_solve(Qc, options);
                int idx = 0;
                for (int i : set_c) {
                    p.set(0, i, ctmc_solve_Qc.get(0, idx++));
                }
            }
            p.divide(p.sumRows(0), p, true);
            return p;
        }

        if (Qmat.getNonZeroLength() == 0) {
            // No transitions at all: every distribution satisfies p*Q=0, so the
            // stationary distribution is not unique and uniform is as good as any.
            Matrix p = new Matrix(1, n);
            p.fill(1.0 / n);
            return p;
        }

        Matrix p = new Matrix(1, n);
        Matrix b = new Matrix(n, 1);
        Matrix nnzel = new Matrix(1, n);
        for (int i = 0; i < n; i++) nnzel.set(0, i, i);
        Matrix Qnnz = Qmat.copy();
        Matrix bnnz = b.copy();
        Matrix Qnnz_1 = Qmat.copy();
        Matrix bnnz_1 = bnnz.copy();

        boolean isReducible = false;
        boolean goon = true;
        while (goon) {
            // nnzel = find(sum(abs(Qnnz),1) > ISOLATED_TOL);
            //
            // AN ISOLATED STATE IS DROPPED; AN ABSORBING ONE IS NOT. The column mass
            // counts a state's inflow plus its own outflow through the diagonal, so it
            // vanishes exactly for a state with neither -- isolated, carrying no
            // stationary mass. This test used to ALSO require a nonzero ROW mass, which
            // vanishes for an ABSORBING state: the one state the mass ends up in.
            // Dropping it left its feeders with nothing to flow into, makeinfgen
            // re-zeroed their diagonals, and the elimination cascaded until nothing was
            // left and this method refused a chain whose stationary distribution is
            // unique (Q = [0 0; 1 -1] has pi = [1 0]).
            //
            // The row mass is also a COMPUTED SUM tested against an exact zero, so the
            // same generator assembled through a different BLAS kernel took opposite
            // branches on two CPUs. In MATLAB that cost SolverMAM's dec.source.mmap a
            // host-dependent answer on the self-looping sanity models; see
            // _kb/06-solver-catalog.md. Native python has always tested the column mass
            // alone (api/mc/ctmc.py).
            Matrix Qnnz_abs = Qnnz.copy();
            Qnnz_abs.absEq();
            Matrix Qnnz_abs_sum_col = Qnnz_abs.sumCols();
            Matrix find_res = new Matrix(1, Qnnz_abs_sum_col.getNumCols());
            for (int i = 0; i < Qnnz_abs_sum_col.getNumCols(); i++) {
                if (Qnnz_abs_sum_col.get(i) > ISOLATED_TOL) find_res.set(0, i, 1);
            }
            nnzel = find_res.find().transpose();

            if (nnzel.length() < n && !isReducible) {
                isReducible = true;
                // if (nargin > 1 && options.verbose == 2) % debug
                // fprintf(1,'ctmc_solve: the infinitesimal generator is reducible.\n');
                // end
            }

            // Qnnz = Qnnz(nnzel, nnzel);
            //
            // Walk Qnnz's NONZEROS through an inverse index instead of probing all
            // nnzel^2 cells. The old form issued a sparse get -- a binary search in a
            // CSC column -- for every pair, so extracting a submatrix with nnz entries
            // cost Theta(nnzel^2 log). Same entries, same values, same != 0.0 filter
            // (a stored explicit zero must still be dropped), so the result is
            // bit-identical; the insertion order becomes column-major, which is also
            // the order CSC wants.
            int nsel = nnzel.getNumCols();
            int[] inv = new int[Qnnz.getNumRows()];
            java.util.Arrays.fill(inv, -1);
            for (int i = 0; i < nsel; i++) {
                inv[(int) nnzel.get(0, i)] = i;
            }
            Matrix new_Qnnz = new Matrix(nsel, nsel, Qnnz.getNonZeros());
            Iterator<MatrixEntry> qit = Qnnz.nonZeroIterator();
            while (qit.hasNext()) {
                MatrixEntry e = qit.next();
                int i = inv[e.row];
                if (i < 0) continue;
                int j = inv[e.col];
                // Copy all non-zero values including negative diagonal elements
                if (j >= 0 && e.value != 0.0) new_Qnnz.set(i, j, e.value);
            }
            Qnnz = new_Qnnz;

            // bnnz = bnnz(nnzel);
            Matrix new_bnnz = new Matrix(nsel, 1, nsel);
            for (int i = 0; i < nsel; i++) {
                new_bnnz.set(i, 0, bnnz.get((int) nnzel.get(0, i), 0));
            }
            bnnz = new_bnnz;

            // Qnnz = ctmc_makeinfgen(Qnnz);
            Qnnz = Ctmc_makeinfgen.ctmc_makeinfgen(Qnnz);

            if ((Qnnz.getNumCols() * Qnnz.getNumRows() == Qnnz_1.getNumCols() * Qnnz_1.getNumRows())
                    && (bnnz.getNumCols() * bnnz.getNumRows() == bnnz_1.getNumCols() * bnnz_1.getNumRows())) {
                goon = false;
            } else {
                Qnnz_1 = Qnnz.copy();
                bnnz_1 = bnnz.copy();
                nnzel = new Matrix(1, Qnnz.length());
                for (int i = 0; i < Qnnz.length(); i++) nnzel.set(0, i, i);
            }
        }

        if (Qnnz == null || Qnnz.isEmpty()) {
            // Every state was ISOLATED -- no inflow and no outflow anywhere -- so the
            // elimination emptied the generator. Filling p with a uniform vector here
            // does NOT satisfy p*Q=0 (it is not a stationary distribution, just a shape
            // of the right size), and a caller cannot tell it apart from a real answer:
            // a generator missing all its arrivals reads back as a plausible mean of
            // cutoff/2. Fail instead. A chain with SEVERAL recurrent classes has no
            // unique stationary distribution without an initial vector either and
            // belongs in ctmc_solve_reducible; a chain with ONE absorbing state is no
            // longer refused, its distribution being the point mass the elimination
            // used to throw away.
            throw new RuntimeException("The infinitesimal generator has no connected state: every state was "
                    + "eliminated as isolated. This generator admits no unique stationary distribution. It "
                    + "usually means the generator is malformed -- e.g. a state with no outgoing transitions "
                    + "that absorbs the whole chain, as happens when a class of transitions was dropped while "
                    + "building it. Use ctmc_solve_reducible for a genuinely absorbing chain.");
        }

        // Qnnz(:,end) = 1;
        for (int i = 0; i < Qnnz.getNumRows(); i++) Qnnz.set(i, Qnnz.getNumCols() - 1, 1.0);

        // bnnz(end) = 1;
        bnnz.set(bnnz.getNumRows() - 1, 0, 1.0);

        // see _kb/03-api-layer.md for rationale
        Matrix Qt = Qnnz.transpose();
        boolean gmresRequested = options != null && "gmres".equalsIgnoreCase(options.method);
        boolean bicgstabRequested = options != null && "bicgstab".equalsIgnoreCase(options.method);
        boolean directRequested = options != null && "direct".equalsIgnoreCase(options.method);
        if (gmresRequested || bicgstabRequested
                || (!directRequested && Qnnz.getNumRows() > GMRES_MIN_STATES)) {
            int restart = 0;
            int maxit = 0;
            if (options != null) {
                if (options.config != null && options.config.gmres_restart > 0) {
                    restart = options.config.gmres_restart;
                }
                if (options.iter_max > 0) {
                    int r = restart > 0 ? restart : Math.min(Qnnz.getNumRows(), Ctmc_gmres.GMRES_DEFAULT_RESTART);
                    maxit = Math.min((int) Math.ceil((double) Qnnz.getNumRows() / r), options.iter_max);
                }
            }
            // GMRES(m) first: its residual is monotone and it is the more robust of the
            // two. The way it fails on a generator is stagnation, the useful subspace
            // being wider than the restart window, and a short-recurrence method has no
            // restart to stagnate on, so BiCGSTAB is tried before the direct solve rather
            // than instead of it. The direct solve is cubic at this size, so the second
            // iterative attempt is cheap against what it may avoid.
            if (!bicgstabRequested) {
                Ctmc_gmres.GmresResult g = Ctmc_gmres.ctmc_gmres(Qt, bnnz, 0.0, restart, maxit, null);
                if (g.flag == 0) {
                    for (int i = 0; i < nnzel.getNumCols(); i++) {
                        p.set(0, (int) nnzel.get(0, i), g.x.get(i, 0));
                    }
                    return p;
                }
            }
            int bmaxit = options != null && options.iter_max > 0 ? options.iter_max : 0;
            Ctmc_bicgstab.BicgstabResult bs = Ctmc_bicgstab.ctmc_bicgstab(Qt, bnnz, 0.0, bmaxit, null);
            if (bs.flag == 0) {
                for (int i = 0; i < nnzel.getNumCols(); i++) {
                    p.set(0, (int) nnzel.get(0, i), bs.x.get(i, 0));
                }
                return p;
            }
        }

        // p(nnzel) = Qnnz' \ bnnz;
        Matrix solve_res = new Matrix(bnnz.getNumRows(), 1);
        // Use solveSafe to match MATLAB behavior - returns NaN for singular matrices instead of throwing
        Matrix.solveSafe(Qt, bnnz, solve_res);
        for (int i = 0; i < nnzel.getNumCols(); i++) p.set(0, (int) nnzel.get(0, i), solve_res.get(i, 0));

        if (p.hasNaN()) {
            // B = abs(Qnnz+Qnnz')>0;
            B = Qnnz.add(1.0, Qnnz.transpose());
            B.absEq();
            int[] colIndexesB2 = B.getColIndexes();
            double[] nzvB2 = B.getNonZeroValues();
            for (int colIdx = 0; colIdx < B.getNumCols(); colIdx++) {
                int col1 = colIndexesB2[colIdx];
                int col2 = colIndexesB2[colIdx + 1];

                for (int i = col1; i < col2; i++) {
                    if (nzvB2[i] > 0) nzvB2[i] = 1.0;
                }
            }

            // [nConnComp, connComp] = weaklyconncomp(B);
            sets = Matrix.weaklyConnect(B, null);
            if (sets.size() > 1) {
                // p(nnzel) = zeros(1,n);
                for (int i = 0; i < nnzel.getNumCols(); i++) p.remove(0, (int) nnzel.get(0, i));

                for (Set<Integer> set_c : sets) {
                    // Qc = Q(connComp==c,connComp==c);
                    Matrix Qc = new Matrix(set_c.size(), set_c.size());
                    int Qc_row = 0;
                    int Qc_col = 0;
                    for (int q_row : set_c) {
                        for (int q_col : set_c) {
                            Qc.set(Qc_row, Qc_col++, Qmat.get(q_row, q_col));
                        }
                        Qc_row++;
                        Qc_col = 0;
                    }
                    Qc = Ctmc_makeinfgen.ctmc_makeinfgen(Qc);
                    // p(intersect(find(connComp==c),nnzel)) = ctmc_solve(Qc);
                    Matrix ctmc_solve_Qc = ctmc_solve(Qc, options);
                    int idx = 0;
                    for (int i = 0; i < nnzel.getNumCols(); i++) {
                        int nnzelValue = (int) nnzel.get(0, i);
                        if (set_c.contains(nnzelValue)) p.set(0, nnzelValue, ctmc_solve_Qc.get(0, idx++));
                    }
                }
                p.divide(p.sumRows(0), p, true);
                return p;
            }
        }
        return p;
    }
}
