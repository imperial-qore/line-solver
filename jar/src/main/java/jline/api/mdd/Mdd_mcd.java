package jline.api.mdd;

import java.util.ArrayList;
import java.util.List;

/**
 * Miner-Ciardo-Donatelli approximate stationary analysis.
 *
 * <p>Solve a structured CTMC whose EXACT reachable state space is stored in a
 * decision diagram, by building and iterating K level-CTMCs (a
 * decision-diagram-guided aggregation), after A.S. Miner, G. Ciardo,
 * S. Donatelli, "Using the exact state space of a Markov model to compute
 * approximate stationary measures", SIGMETRICS 2000.</p>
 *
 * <p>The method never forms the |S|-state generator or probability vector. It
 * keeps one CTMC per decision-diagram level k, over states M_k = {(p,i_k)} with
 * p a level-k node and i_k a local state on a non-null arc, and iterates the
 * coupled system to a fixed point. The single approximation (Eq. 5) is
 * Pr{i_k | alpha} = Pr{i_k | p}: the local-state law at level k depends only on
 * the node p, not the full path above it, which the exact reachability the
 * diagram encodes justifies. For product-form models the method is EXACT (paper
 * Sec. 5), so on a single-class closed QN it reproduces SolverCTMC.</p>
 *
 * <p><b>Orientation.</b> The paper indexes levels K (top/root) down to 1
 * (bottom/terminal); {@link MDD} uses level 0 as the root. This class works in
 * the paper's orientation with 0-based indices, so paper level k (0 = bottom)
 * maps to MDD level K-1-k and to station K-1-k. Getting this backwards silently
 * mislabels every per-station metric.</p>
 */
public class Mdd_mcd {

    private Mdd_mcd() {}

    /** Run the aggregation with the default level-iteration knobs. */
    public static MddMcdResult mdd_mcd(MddStruct mdds, MddDescriptor desc) {
        return mdd_mcd(mdds, desc, new MddMcdOptions());
    }

    /**
     * Approximate stationary measures by decision-diagram-guided aggregation.
     *
     * @param mdds the reachable set, in MDD orientation
     * @param desc the Kronecker rate descriptor
     * @param options level-iteration knobs
     */
    public static MddMcdResult mdd_mcd(MddStruct mdds, MddDescriptor desc, MddMcdOptions options) {
        if (options == null) {
            options = new MddMcdOptions();
        }
        final int K = mdds.K;

        // ---- paper orientation: paper level k <-> MDD level K-1-k = station K-1-k
        int[][][] Pnode = new int[K][][];
        int[] nn = new int[K];
        int[] dom = new int[K];
        for (int k = 0; k < K; k++) {
            int oL = K - 1 - k;
            Pnode[k] = mdds.node[oL];
            nn[k] = mdds.nnodes[oL];
            dom[k] = mdds.domain[oL];
        }

        // ---- per (event, paper level) local matrices W_k^e; an untouched level
        // carries the identity, which is supplied here rather than stored
        final int E = desc.events.size();
        MddLocalMatrix[][] W = new MddLocalMatrix[E][K];
        for (int e = 0; e < E; e++) {
            MddEvent ev = desc.events.get(e);
            for (int k = 0; k < K; k++) {
                W[e][k] = null;
            }
            for (int t = 0; t < ev.lev.length; t++) {
                W[e][K - 1 - ev.lev[t]] = ev.W[t];   // station level -> paper level
            }
            for (int k = 0; k < K; k++) {
                if (W[e][k] == null) {
                    W[e][k] = MddLocalMatrix.identity(dom[k]);
                }
            }
        }

        // ---- level-k CTMC state sets M_k = {(p, v) : arc p[v] non-null}
        int[][][] Mrows = new int[K][][];
        int[][] Midx = new int[K][];
        int[] levelSizes = new int[K];
        for (int k = 0; k < K; k++) {
            List<int[]> rows = new ArrayList<int[]>();
            // column-major order, matching the MATLAB find() the reference uses
            for (int v = 0; v < dom[k]; v++) {
                for (int p = 0; p < nn[k]; p++) {
                    int ch = Pnode[k][p][v];
                    boolean live = (k == 0) ? (ch == MDD.TERM_TRUE) : (ch > 0);
                    if (live) {
                        rows.add(new int[]{p + 1, v});
                    }
                }
            }
            Mrows[k] = rows.toArray(new int[rows.size()][]);
            levelSizes[k] = rows.size();
            Midx[k] = new int[nn[k] * dom[k]];
            for (int r = 0; r < Mrows[k].length; r++) {
                Midx[k][(Mrows[k][r][0] - 1) * dom[k] + Mrows[k][r][1]] = r + 1;
            }
        }

        // ---- initialise level stationary vectors and node marginals
        double[][][] counts = pathCounts(mdds, K);
        double[][] above = counts[0];
        double[][] pik;
        if (options.initpik != null) {
            pik = new double[K][];
            for (int k = 0; k < K; k++) {
                pik[k] = options.initpik[k].clone();
            }
        } else {
            pik = uniformInit(mdds, Mrows, K, above, counts[1]);
        }
        double[][] Prp = new double[K][];
        for (int k = 0; k < K; k++) {
            Prp[k] = nodeMarginal(Mrows[k], pik[k], nn[k]);
        }

        // ---- fixed-point iteration (Fig. 3, procedure Solve)
        int iters = 0;
        boolean converged = false;
        double delta = Double.POSITIVE_INFINITY;
        for (int it = 1; it <= options.maxiter; it++) {
            iters = it;
            double[][] piold = new double[K][];
            for (int k = 0; k < K; k++) {
                piold[k] = pik[k].clone();
            }

            // ComputeBs, bottom-up:
            // b_k^e[p] = sum_v Pr{v|p} * b_{k-1}^e[p[v]] * lambda
            double[][][] bcell = new double[K][][];
            for (int k = 0; k < K; k++) {
                double[][] bk = new double[nn[k]][E];
                int[][] rows = Mrows[k];
                for (int r = 0; r < rows.length; r++) {
                    int p = rows[r][0];
                    int v = rows[r][1];
                    if (Prp[k][p - 1] <= 0) {
                        continue;
                    }
                    double adjust = pik[k][r] / Prp[k][p - 1];       // Pr{v|p}
                    for (int e = 0; e < E; e++) {
                        double le = W[e][k].rowSum[v];
                        if (le == 0) {
                            continue;                                // not locally enabled
                        }
                        double down = 1.0;                           // terminal ONE
                        if (k > 0) {
                            down = bcell[k - 1][Pnode[k][p - 1][v] - 1][e];
                        }
                        bk[p - 1][e] += adjust * down * le;
                    }
                }
                bcell[k] = bk;
            }

            // top-down: ComputeAs(k) then SolveLevel(k)
            double[][][][] Acell = new double[K][][][];
            Acell[K - 1] = new double[E][][];
            for (int e = 0; e < E; e++) {
                Acell[K - 1][e] = identity(nn[K - 1]);
            }
            for (int k = K - 1; k >= 0; k--) {
                if (k < K - 1) {
                    Acell[k] = computeAs(k, Acell[k + 1], Pnode, pik, W, Mrows, nn, E);
                }
                // SolveLevel(k): assemble R_k (Eq. 6), solve pi_k Q_k = 0
                double[][] Rk = computeMC(k, Acell[k], bcell, Pnode, W, Mrows, Midx,
                        levelSizes, dom, E);
                double[][] Qk = generator(Rk);
                pik[k] = solveStat(Qk);
                Prp[k] = nodeMarginal(Mrows[k], pik[k], nn[k]);
            }

            // A diverged iterate must not be read as converged, which is what a
            // NaN-skipping maximum would do.
            delta = 0;
            for (int k = 0; k < K; k++) {
                double dk = 0;
                for (int r = 0; r < pik[k].length; r++) {
                    dk = Math.max(dk, Math.abs(pik[k][r] - piold[k][r]));
                }
                if (Double.isNaN(dk) || Double.isInfinite(dk)) {
                    throw new RuntimeException("mdd_mcd: level " + (k + 1)
                            + " iterate is not finite at iteration " + it + "; the level-"
                            + (k + 1) + " CTMC did not yield a proper stationary vector.");
                }
                delta = Math.max(delta, dk);
            }
            if (delta < options.tol) {
                converged = true;
                break;
            }
        }
        if (!converged) {
            throw new RuntimeException(String.format(
                    "mdd_mcd: the coupled level iteration did not converge in %d sweeps "
                    + "(last change %.3e against tol %.3e); the level marginals returned would "
                    + "not be a fixed point. Raise maxiter or relax tol.",
                    options.maxiter, delta, options.tol));
        }

        // ---- performance measures from the per-level marginals. QLen is the
        // mean local value and is defined for any descriptor (jobs at a station,
        // tokens in a place); X and U need the queueing parameters.
        boolean isQN = desc.mu != null && desc.mu.length > 0 && desc.servers != null;
        double[] QLen = new double[K];
        double[] X = isQN ? new double[K] : null;
        double[] U = isQN ? new double[K] : null;
        for (int s = 0; s < K; s++) {
            int k = K - 1 - s;                       // paper level of station s
            int[][] rows = Mrows[k];
            double[] pk = pik[k];
            double q = 0;
            double busySum = 0;
            for (int r = 0; r < rows.length; r++) {
                double v = (desc.valuemap != null)
                        ? desc.valuemap[s][rows[r][1]] : rows[r][1];
                q += v * pk[r];
                if (isQN) {
                    busySum += Math.min(v, desc.servers[s]) * pk[r];
                }
            }
            QLen[s] = q;
            if (isQN) {
                X[s] = desc.mu[s] * busySum;
                U[s] = Double.isInfinite(desc.servers[s]) ? q : busySum / desc.servers[s];
            }
        }

        // The level chains are coupled only through rates, so nothing in the
        // iteration forces the marginals to describe the same population; a fixed
        // point that does not is a wrong answer, not an approximation, and must
        // not be returned. The test is a conservation law of the model: the
        // closed population for a QN, a place invariant w'*m = const for a net.
        double[] winv = desc.invariantWeights;
        double vinv = desc.invariantValue;
        if (winv == null && desc.N > 0) {
            winv = new double[K];
            for (int k = 0; k < K; k++) {
                winv[k] = 1.0;                       // closed QN: total population
            }
            vinv = desc.N;
        }
        if (winv != null) {
            double got = 0;
            for (int s = 0; s < K; s++) {
                got += winv[s] * QLen[s];
            }
            if (Math.abs(got - vinv) > 1e-6 * Math.max(1.0, Math.abs(vinv))) {
                throw new RuntimeException(String.format(
                        "mdd_mcd: the level marginals converged to an invariant value of %.6g "
                        + "against the model value %g, so the fixed point reached is degenerate "
                        + "(the level chains are mutually inconsistent). Supply initpik with a "
                        + "consistent starting law.", got, vinv));
            }
        }

        MddMcdResult out = new MddMcdResult();
        out.QLen = QLen;
        out.X = X;
        out.U = U;
        out.pik = pik;
        out.Mrows = Mrows;
        out.levelSizes = levelSizes;
        out.iters = iters;
        out.pathsPerLevel = new double[K];
        boolean noAgg = true;
        for (int k = 0; k < K; k++) {
            int oL = K - 1 - k;
            double mx = 1.0;
            for (int p = 0; p < above[oL].length; p++) {
                mx = Math.max(mx, above[oL][p]);
            }
            out.pathsPerLevel[k] = mx;
            if (mx > 1.0 + 1e-12) {
                noAgg = false;
            }
        }
        out.noAggregation = noAgg;

        if (options.verbose) {
            StringBuilder sb = new StringBuilder();
            int total = 0;
            int max = 0;
            for (int k = 0; k < K; k++) {
                sb.append(levelSizes[k]);
                if (k < K - 1) {
                    sb.append(' ');
                }
                total += levelSizes[k];
                max = Math.max(max, levelSizes[k]);
            }
            System.out.format("%nMDD-MCD approximate aggregation: %d levels, %d fixed-point iters%n",
                    K, iters);
            System.out.format("  level-CTMC sizes |M_k| = [%s] (max %d, total %d)%n",
                    sb.toString(), max, total);
        }
        return out;
    }

    // -----------------------------------------------------------------------
    private static double[] nodeMarginal(int[][] rows, double[] pk, int nnodes) {
        double[] pr = new double[nnodes];
        for (int r = 0; r < rows.length; r++) {
            pr[rows[r][0] - 1] += pk[r];
        }
        return pr;
    }

    private static double[][] identity(int n) {
        double[][] I = new double[n][n];
        for (int i = 0; i < n; i++) {
            I[i][i] = 1.0;
        }
        return I;
    }

    private static double[][] generator(double[][] R) {
        int n = R.length;
        double[][] Q = new double[n][n];
        for (int i = 0; i < n; i++) {
            double s = 0;
            for (int j = 0; j < n; j++) {
                Q[i][j] = R[i][j];
                s += R[i][j];
            }
            Q[i][i] -= s;
        }
        return Q;
    }

    /**
     * Path counts per MDD level: above[oL][p] is the number of distinct
     * root-to-p paths, |A(p)| in the paper's notation, and below[oL][p] the
     * number of accepted states under p.
     *
     * <p>Both are O(#nodes) and serve two purposes: the uniform initialisation,
     * and the exactness certificate. The single approximation is
     * Pr{i_k | alpha} = Pr{i_k | p}, so when a node is reached by exactly one
     * path, conditioning on the node IS conditioning on the path and the
     * identity is exact. If that holds at every node the fixed point is the
     * exact stationary law, with no reference solve needed to know it. That test
     * is SUFFICIENT, not necessary: a product-form model is exact too however
     * much its diagram shares. Note also that max |A(p)| = 1 means no node is
     * shared, i.e. the diagram compresses nothing, so exactness by this route
     * and a useful saving are mutually exclusive.</p>
     *
     * @return {above, below}
     */
    private static double[][][] pathCounts(MddStruct mdds, int K) {
        double[][] below = new double[K][];
        double[][] above = new double[K][];
        for (int oL = K - 1; oL >= 0; oL--) {
            double[] nb = new double[mdds.nnodes[oL]];
            for (int p = 0; p < mdds.nnodes[oL]; p++) {
                double s = 0;
                for (int v = 0; v < mdds.domain[oL]; v++) {
                    int ch = mdds.node[oL][p][v];
                    if (oL == K - 1) {
                        if (ch == MDD.TERM_TRUE) {
                            s += 1.0;
                        }
                    } else if (ch > 0) {
                        s += below[oL + 1][ch - 1];
                    }
                }
                nb[p] = s;
            }
            below[oL] = nb;
        }
        for (int oL = 0; oL < K; oL++) {
            above[oL] = new double[mdds.nnodes[oL]];
        }
        above[0][mdds.root - 1] = 1.0;
        for (int oL = 0; oL < K - 1; oL++) {
            for (int p = 0; p < mdds.nnodes[oL]; p++) {
                double w = above[oL][p];
                if (w == 0) {
                    continue;
                }
                for (int v = 0; v < mdds.domain[oL]; v++) {
                    int ch = mdds.node[oL][p][v];
                    if (ch > 0) {
                        above[oL + 1][ch - 1] += w;
                    }
                }
            }
        }
        return new double[][][]{above, below};
    }

    /**
     * Uniform law over the EXACT reachable set, projected onto each level.
     *
     * <p>Pr{(p,v)} = (paths root-&gt;p) * (states below arc p[v]) / |S|. A flat
     * law over M_k instead treats level states as equiprobable irrespective of
     * how many global states they stand for, which breaks the population
     * invariant the diagram encodes; from about K=8 the coupled iteration then
     * descends into the basin of the DEGENERATE empty-population fixed point
     * (all mass on local state 0 at every level, a true fixed point since no
     * station can then emit) and converges to it with zero residual. The
     * projection below is consistent across levels by construction, so the
     * iteration starts inside the physical simplex.</p>
     */
    private static double[][] uniformInit(MddStruct mdds, int[][][] Mrows, int K,
                                          double[][] above, double[][] below) {
        double[][] pik = new double[K][];
        for (int k = 0; k < K; k++) {
            int oL = K - 1 - k;                      // paper level k is MDD level K-1-k
            int[][] rows = Mrows[k];
            double[] w = new double[rows.length];
            double sum = 0;
            for (int r = 0; r < rows.length; r++) {
                int p = rows[r][0];
                int v = rows[r][1];
                double val;
                if (oL == K - 1) {
                    val = above[oL][p - 1];          // a TRUE arc stands for one state
                } else {
                    int ch = mdds.node[oL][p - 1][v];
                    val = above[oL][p - 1] * below[oL + 1][ch - 1];
                }
                w[r] = val;
                sum += val;
            }
            for (int r = 0; r < rows.length; r++) {
                w[r] /= sum;
            }
            pik[k] = w;
        }
        return pik;
    }

    /**
     * ComputeAs(k): A_k^e from A_{k+1}^e, the "from above" contribution (Fig. 3).
     *
     * <p>The adjust denominator Pr{p[v]} is the FROM-ABOVE marginal of the child
     * node, Pr{p} = sum over parents of pi_{k+1}. Using it (rather than the
     * level-k CTMC marginal, which only equals it at convergence) makes adjust a
     * proper conditional Pr{(parent,arc)|child} and pins the inter-level node
     * marginals, removing the spurious fixed points.</p>
     */
    private static double[][][] computeAs(int k, double[][][] Aup, int[][][] Pnode,
                                          double[][] pik, MddLocalMatrix[][] W,
                                          int[][][] Mrows, int[] nn, int E) {
        int[][] rows1 = Mrows[k + 1];
        double[] PrAbove = new double[nn[k]];
        for (int r = 0; r < rows1.length; r++) {
            int child = Pnode[k + 1][rows1[r][0] - 1][rows1[r][1]];
            if (child > 0) {
                PrAbove[child - 1] += pik[k + 1][r];
            }
        }
        double[][][] Ak = new double[E][][];
        for (int e = 0; e < E; e++) {
            Ak[e] = new double[nn[k]][nn[k]];
        }
        for (int r = 0; r < rows1.length; r++) {
            int p = rows1[r][0];
            int v = rows1[r][1];
            int childp = Pnode[k + 1][p - 1][v];     // p[v]: node at level k
            if (childp <= 0 || PrAbove[childp - 1] <= 0) {
                continue;
            }
            double adjust = pik[k + 1][r] / PrAbove[childp - 1];
            for (int e = 0; e < E; e++) {
                int[] wcols = W[e][k + 1].cols[v];
                if (wcols.length == 0) {
                    continue;
                }
                double[] wvals = W[e][k + 1].vals[v];
                double[] arow = Aup[e][p - 1];
                for (int wi = 0; wi < wcols.length; wi++) {
                    int w = wcols[wi];
                    double wv = wvals[wi];
                    for (int q = 0; q < arow.length; q++) {
                        if (arow[q] == 0) {
                            continue;
                        }
                        int childq = Pnode[k + 1][q][w];
                        if (childq <= 0) {
                            continue;                // q[w] null
                        }
                        Ak[e][childp - 1][childq - 1] += arow[q] * wv * adjust;
                    }
                }
            }
        }
        return Ak;
    }

    /**
     * ComputeMC(k): level-k rate matrix (Eq. 6),
     * R_k^e[(p,i),(q,j)] = A_k^e[p,q] * W_k^e[i,j] * b_{k-1}^e[p[i]].
     */
    private static double[][] computeMC(int k, double[][][] Ak, double[][][] bcell,
                                        int[][][] Pnode, MddLocalMatrix[][] W,
                                        int[][][] Mrows, int[][] Midx, int[] levelSizes,
                                        int[] dom, int E) {
        int nm = levelSizes[k];
        int[][] rows = Mrows[k];
        double[][] Rk = new double[nm][nm];
        for (int r = 0; r < nm; r++) {
            int p = rows[r][0];
            int v = rows[r][1];
            for (int e = 0; e < E; e++) {
                int[] wcols = W[e][k].cols[v];
                if (wcols.length == 0) {
                    continue;
                }
                double bfac = 1.0;                   // terminal ONE
                if (k > 0) {
                    bfac = bcell[k - 1][Pnode[k][p - 1][v] - 1][e];
                }
                if (bfac == 0) {
                    continue;
                }
                double[] wvals = W[e][k].vals[v];
                double[] arow = Ak[e][p - 1];
                for (int wi = 0; wi < wcols.length; wi++) {
                    int w = wcols[wi];
                    double wv = wvals[wi];
                    for (int q = 0; q < arow.length; q++) {
                        if (arow[q] == 0) {
                            continue;
                        }
                        int di = Midx[k][q * dom[k] + w];
                        if (di == 0) {
                            continue;                // (q,w) not in M_k
                        }
                        Rk[r][di - 1] += arow[q] * wv * bfac;
                    }
                }
            }
        }
        return Rk;
    }

    /**
     * Stationary distribution of a small irreducible generator: p Q = 0, sum p = 1.
     *
     * <p>The normalisation is APPENDED rather than substituted for the last
     * balance equation: overwriting a row discards a constraint and leaves Q'
     * singular to working precision from about |M_k| = 325 upwards, so the solve
     * returns NaN. The overdetermined system has full column rank whenever the
     * level chain is irreducible, and least squares through the normal equations
     * of the QR-equivalent system solves it stably.</p>
     */
    private static double[] solveStat(double[][] Q) {
        int n = Q.length;
        if (n == 1) {
            return new double[]{1.0};
        }
        // A = [Q'; ones], rhs = e_{n+1}; solve the least-squares problem
        double[][] A = new double[n + 1][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                A[j][i] = Q[i][j];
            }
        }
        for (int j = 0; j < n; j++) {
            A[n][j] = 1.0;
        }
        double[] rhs = new double[n + 1];
        rhs[n] = 1.0;
        double[] p = lstsq(A, rhs);
        double s = 0;
        for (int i = 0; i < n; i++) {
            if (p[i] < 0) {
                p[i] = 0;
            }
            s += p[i];
        }
        if (Double.isNaN(s) || Double.isInfinite(s) || s <= 0) {
            throw new RuntimeException("mdd_mcd: level CTMC of order " + n
                    + " admits no proper stationary distribution (the level generator is "
                    + "reducible or numerically degenerate).");
        }
        for (int i = 0; i < n; i++) {
            p[i] /= s;
        }
        return p;
    }

    /**
     * Least-squares solution of an overdetermined system by Householder QR.
     *
     * <p>QR is used rather than the normal equations because A'A squares the
     * condition number, which is exactly the failure the appended-normalisation
     * form is there to avoid.</p>
     */
    private static double[] lstsq(double[][] A, double[] b) {
        int m = A.length;
        int n = A[0].length;
        double[][] R = new double[m][n];
        for (int i = 0; i < m; i++) {
            System.arraycopy(A[i], 0, R[i], 0, n);
        }
        double[] y = b.clone();
        for (int k = 0; k < n; k++) {
            double norm = 0;
            for (int i = k; i < m; i++) {
                norm += R[i][k] * R[i][k];
            }
            norm = Math.sqrt(norm);
            if (norm == 0) {
                continue;
            }
            if (R[k][k] > 0) {
                norm = -norm;
            }
            double[] v = new double[m];
            for (int i = k; i < m; i++) {
                v[i] = R[i][k];
            }
            v[k] -= norm;
            double vtv = 0;
            for (int i = k; i < m; i++) {
                vtv += v[i] * v[i];
            }
            if (vtv == 0) {
                continue;
            }
            for (int j = k; j < n; j++) {
                double dot = 0;
                for (int i = k; i < m; i++) {
                    dot += v[i] * R[i][j];
                }
                double f = 2.0 * dot / vtv;
                for (int i = k; i < m; i++) {
                    R[i][j] -= f * v[i];
                }
            }
            double dot = 0;
            for (int i = k; i < m; i++) {
                dot += v[i] * y[i];
            }
            double f = 2.0 * dot / vtv;
            for (int i = k; i < m; i++) {
                y[i] -= f * v[i];
            }
        }
        double[] x = new double[n];
        for (int i = n - 1; i >= 0; i--) {
            double s = y[i];
            for (int j = i + 1; j < n; j++) {
                s -= R[i][j] * x[j];
            }
            x[i] = R[i][i] == 0 ? 0 : s / R[i][i];
        }
        return x;
    }
}
