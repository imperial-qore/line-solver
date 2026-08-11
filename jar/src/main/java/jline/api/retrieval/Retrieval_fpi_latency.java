/**
 * @file Retrieval_fpi_latency.java
 * @brief FPI-based delayed-hit count and expected latency for a delayed-hit cache.
 *
 * Port of matlab/src/api/retrieval/retrieval_fpi_latency.m (paper thm:di,
 * eq:moments, eq:latency tot). Per item i the fetch sojourn moments are obtained
 * from a reduced absorbing CTMC of its visits to the retrieval stations:
 *   E0[F^k] = k! * pe * (-D0)^{-k} * e,   d_i = phi_i*lambda_i*E0[F^2]/(2 E0[F]),
 *   Z = sum_i (phi_i + d_i) / sum_i lambda_i (phi_i + pi_{i,0}).
 *
 * Supported station types: IS, PS, SIRO, FCFS, LCFSPR. PS/SIRO/FCFS/LCFSPR use
 * the mean-field sharing slowdown 1/(1+phitilde); SIRO and FCFS require
 * exponential (single-phase) service with identical per-class rates.
 *
 * Representation: alpha[s][i] is the (1 x fsz_s) PH entry vector of item i at
 * station s; T[s][i] is its (fsz_s x fsz_s) PH subgenerator; R[i] is the
 * (S+1 x S+1) routing matrix of item i (index 0 = outside/cache, 1..S = stations).
 *
 * @since LINE 3.0
 */
package jline.api.retrieval;

public final class Retrieval_fpi_latency {
    private Retrieval_fpi_latency() {}

    /** Result: Z (latency), d (per-item delayed-hit count), phi, pi0 (each length n). */
    public static final class Result {
        public final double Z;
        public final double[] d;
        public final double[] phi;
        public final double[] pi0;
        Result(double Z, double[] d, double[] phi, double[] pi0) { this.Z = Z; this.d = d; this.phi = phi; this.pi0 = pi0; }
    }

    public static Result retrieval_fpi_latency(double[] m, double[] lambda, double[][] gamma,
                                               double[][][] alpha, double[][][][] T,
                                               double[][][] R, String[] stationType) {
        int n = lambda.length;
        int S = stationType.length;
        int[] fsz = new int[S];
        for (int s = 0; s < S; s++) fsz[s] = T[s][0].length;

        // --- supported types + SIRO/FCFS constraints ---
        for (int s = 0; s < S; s++) {
            String st = stationType[s];
            if (!st.equals("IS") && !st.equals("PS") && !st.equals("SIRO") && !st.equals("FCFS") && !st.equals("LCFSPR")) {
                throw new RuntimeException("retrieval_fpi_latency supports only IS, PS, SIRO, FCFS and LCFSPR stations; got " + st);
            }
            if ((st.equals("SIRO") || st.equals("FCFS")) && fsz[s] > 1) {
                throw new RuntimeException("retrieval_fpi_latency supports SIRO/FCFS only with exponential (single-phase) service; station " + s);
            }
            if (st.equals("SIRO") || st.equals("FCFS")) {
                double tau0 = phMean(alpha[s][0], T[s][0]);
                for (int i = 1; i < n; i++) {
                    if (Math.abs(phMean(alpha[s][i], T[s][i]) - tau0) > 1e-9 * Math.max(tau0, 1e-300)) {
                        throw new RuntimeException("retrieval_fpi_latency requires identical per-class rates at SIRO/FCFS station " + s);
                    }
                }
            }
        }

        // psIdx = PS/SIRO/FCFS/LCFSPR stations, isIdx = IS stations
        int[] psIdx = new int[S]; int r = 0;
        boolean[] isIS = new boolean[S];
        for (int s = 0; s < S; s++) {
            isIS[s] = stationType[s].equals("IS");
            if (!isIS[s]) psIdx[r++] = s;
        }

        // --- eta_fpi(i, 0)=sum IS, eta_fpi(i,1+p)=PS station psIdx[p] ---
        double[][] eta = new double[n][r + 1];
        for (int i = 0; i < n; i++) {
            double[] visits = visits(R[i], S);
            double[] tau = new double[S];
            for (int s = 0; s < S; s++) tau[s] = phMean(alpha[s][i], T[s][i]);
            for (int s = 0; s < S; s++) if (isIS[s]) eta[i][0] += visits[s] * tau[s];
            for (int p = 0; p < r; p++) eta[i][1 + p] = visits[psIdx[p]] * tau[psIdx[p]];
        }

        // --- step 1: FPI on the full system ---
        Retrieval_fpi.Result full = Retrieval_fpi.retrieval_fpi(m, lambda, eta, gamma);
        double[] pi0 = full.pmiss;
        double[] phi = new double[n];
        for (int i = 0; i < n; i++) { double sum = 0; for (double[] row : full.pdh) sum += row[i]; phi[i] = sum; }

        double[] d = new double[n];
        for (int i = 0; i < n; i++) {
            // --- step 2: FPI without item i -> per-station occupancy phitilde ---
            int nn = n - 1;
            double[] lambda_i = new double[nn];
            double[][] eta_i = new double[nn][];
            double[][] gamma_i = new double[nn][];
            int idx = 0;
            for (int k = 0; k < n; k++) { if (k == i) continue; lambda_i[idx] = lambda[k]; eta_i[idx] = eta[k]; gamma_i[idx] = gamma[k]; idx++; }
            Retrieval_fpi.Result without = Retrieval_fpi.retrieval_fpi(m, lambda_i, eta_i, gamma_i);
            double[] phitilde = new double[S];
            for (int p = 0; p < r; p++) {
                double sum = 0; for (int kk = 0; kk < nn; kk++) sum += without.pdh[1 + p][kk];
                phitilde[psIdx[p]] = sum;
            }

            // --- step 3: reduced absorbing CTMC for item i ---
            int Phi = 0; for (int s = 0; s < S; s++) Phi += fsz[s];
            int[] off = new int[S]; for (int s = 1; s < S; s++) off[s] = off[s - 1] + fsz[s - 1];
            double[][] D0 = new double[Phi][Phi];
            double[] pe = new double[Phi];
            double[][] Ri = R[i];
            for (int s = 0; s < S; s++) {
                double scale = isIS[s] ? 1.0 : 1.0 / (1.0 + phitilde[s]);
                double[][] blk = new double[fsz[s]][fsz[s]];
                for (int a = 0; a < fsz[s]; a++) for (int b = 0; b < fsz[s]; b++) blk[a][b] = scale * T[s][i][a][b];
                // diagonal block
                for (int a = 0; a < fsz[s]; a++) for (int b = 0; b < fsz[s]; b++) D0[off[s] + a][off[s] + b] += blk[a][b];
                // completion-rate vector compl = -blk*e
                double[] compl = new double[fsz[s]];
                for (int a = 0; a < fsz[s]; a++) { double rs = 0; for (int b = 0; b < fsz[s]; b++) rs += blk[a][b]; compl[a] = -rs; }
                // route to station sp with R(s+1,sp+1)*alpha_sp
                for (int sp = 0; sp < S; sp++) {
                    double rprob = Ri[s + 1][sp + 1];
                    if (rprob == 0) continue;
                    for (int a = 0; a < fsz[s]; a++)
                        for (int b = 0; b < fsz[sp]; b++)
                            D0[off[s] + a][off[sp] + b] += compl[a] * rprob * alpha[sp][i][b];
                }
                // entry distribution
                for (int b = 0; b < fsz[s]; b++) pe[off[s] + b] = Ri[0][s + 1] * alpha[s][i][b];
            }

            // --- step 4: moments via A=-D0, x1=A\e, x2=A\x1 ---
            double[][] A = new double[Phi][Phi];
            for (int a = 0; a < Phi; a++) for (int b = 0; b < Phi; b++) A[a][b] = -D0[a][b];
            double[] e = new double[Phi]; java.util.Arrays.fill(e, 1.0);
            double[] x1 = solve(A, e);
            double[] x2 = solve(A, x1);
            double M1 = 0, M2 = 0;
            for (int a = 0; a < Phi; a++) { M1 += pe[a] * x1[a]; M2 += pe[a] * x2[a]; }
            M2 *= 2;
            d[i] = phi[i] * lambda[i] * M2 / (2 * M1);
        }

        // --- step 5: Z ---
        double num = 0, den = 0;
        for (int i = 0; i < n; i++) { num += phi[i] + d[i]; den += lambda[i] * (phi[i] + pi0[i]); }
        return new Result(num / den, d, phi, pi0);
    }

    /** Mean of PH (alpha,T): -alpha*inv(T)*e = alpha*inv(-T)*e. */
    private static double phMean(double[] al, double[][] Tm) {
        int k = al.length;
        double[][] negT = new double[k][k];
        for (int a = 0; a < k; a++) for (int b = 0; b < k; b++) negT[a][b] = -Tm[a][b];
        double[] e = new double[k]; java.util.Arrays.fill(e, 1.0);
        double[] z = solve(negT, e);
        double mu = 0; for (int a = 0; a < k; a++) mu += al[a] * z[a];
        return mu;
    }

    /** visits = a * (I-P)^{-1} where a=R(0,1..S), P=R(1..S,1..S). */
    private static double[] visits(double[][] R, int S) {
        double[] a = new double[S];
        double[][] ImP = new double[S][S];
        for (int s = 0; s < S; s++) {
            a[s] = R[0][s + 1];
            for (int sp = 0; sp < S; sp++) ImP[s][sp] = (s == sp ? 1.0 : 0.0) - R[s + 1][sp + 1];
        }
        // solve visits*(I-P) = a  =>  (I-P)^T * visits^T = a^T
        double[][] At = new double[S][S];
        for (int s = 0; s < S; s++) for (int sp = 0; sp < S; sp++) At[s][sp] = ImP[sp][s];
        return solve(At, a);
    }

    /** Solve A x = b (dense, partial pivoting). */
    private static double[] solve(double[][] Ain, double[] bin) {
        int N = bin.length;
        double[][] A = new double[N][N];
        double[] b = bin.clone();
        for (int i = 0; i < N; i++) A[i] = Ain[i].clone();
        for (int col = 0; col < N; col++) {
            int piv = col; double best = Math.abs(A[col][col]);
            for (int rrow = col + 1; rrow < N; rrow++) if (Math.abs(A[rrow][col]) > best) { best = Math.abs(A[rrow][col]); piv = rrow; }
            if (piv != col) { double[] tr = A[piv]; A[piv] = A[col]; A[col] = tr; double tb = b[piv]; b[piv] = b[col]; b[col] = tb; }
            double d = A[col][col];
            for (int rrow = col + 1; rrow < N; rrow++) {
                double f = A[rrow][col] / d;
                if (f == 0) continue;
                for (int c = col; c < N; c++) A[rrow][c] -= f * A[col][c];
                b[rrow] -= f * b[col];
            }
        }
        double[] x = new double[N];
        for (int rrow = N - 1; rrow >= 0; rrow--) {
            double sum = b[rrow];
            for (int c = rrow + 1; c < N; c++) sum -= A[rrow][c] * x[c];
            x[rrow] = sum / A[rrow][rrow];
        }
        return x;
    }
}
