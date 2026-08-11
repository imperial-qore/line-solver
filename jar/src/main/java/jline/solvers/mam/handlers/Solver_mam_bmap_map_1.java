/**
 * BMAP/MAP/1 Queue Solver using M/G/1 type analysis with ETAQA.
 *
 * Solves a BMAP/MAP/1 queue where:
 * - Arrivals follow a Batch Markovian Arrival Process (BMAP)
 * - Service follows a Markovian Arrival Process (MAP) for single service
 * - Single server
 *
 * The queue is modeled as an M/G/1-type Markov chain because:
 * - BMAP arrivals can increase level by 1, 2, 3, ... (batch sizes)
 * - MAP service decreases level by exactly 1
 *
 * @since LINE 3.1.0
 */
package jline.solvers.mam.handlers;

import java.util.List;

import jline.api.mam.Map_prob;
import jline.lib.smc.MG1_ETAQA;
import jline.util.matrix.Matrix;

public final class Solver_mam_bmap_map_1 {
    private Solver_mam_bmap_map_1() {}

    /**
     * Solves a BMAP/MAP/1 queue using M/G/1 type matrix-analytic methods.
     */
    public static BMAPMAP1Result solver_mam_bmap_map_1(List<Matrix> D, Matrix S0, Matrix S1, int nMoments) {
        if (D.isEmpty()) {
            throw new IllegalArgumentException("BMAP must have at least D0 matrix");
        }
        if (D.size() < 2) {
            throw new IllegalArgumentException("BMAP must have at least D0 and D1 matrices");
        }

        int K = D.size() - 1; // Maximum batch size
        int ma = D.get(0).getNumRows(); // Number of BMAP phases
        int ms = S0.getNumRows();       // Number of MAP service phases
        int m = ma * ms;                // Combined phases per level

        // Validate BMAP matrices are square and same size
        for (int i = 0; i < D.size(); i++) {
            if (D.get(i).getNumRows() != ma || D.get(i).getNumCols() != ma) {
                throw new IllegalArgumentException("All BMAP matrices must be " + ma + "x" + ma);
            }
        }
        if (S0.getNumRows() != ms || S0.getNumCols() != ms) {
            throw new IllegalArgumentException("S0 must be " + ms + "x" + ms);
        }
        if (S1.getNumRows() != ms || S1.getNumCols() != ms) {
            throw new IllegalArgumentException("S1 must be " + ms + "x" + ms);
        }

        // Compute arrival rate from BMAP
        Matrix D0 = D.get(0);
        Matrix D1Total = Matrix.zeros(ma, ma);
        for (int k = 1; k <= K; k++) {
            D1Total = D1Total.add(D.get(k));
        }

        Matrix piBmap = Map_prob.map_prob(D0, D1Total);
        Matrix eMa = Matrix.ones(ma, 1);

        // Total customer arrival rate: sum_k (k * pi * Dk * e)
        double lambdaTotal = 0.0;
        double batchRate = 0.0;
        for (int k = 1; k <= K; k++) {
            double rateK = piBmap.mult(D.get(k)).mult(eMa).get(0, 0);
            lambdaTotal += k * rateK;
            batchRate += rateK;
        }

        // Mean batch size
        double meanBatchSize = (batchRate > 0) ? lambdaTotal / batchRate : 0.0;

        // Compute service rate from MAP
        Matrix piMap = Map_prob.map_prob(S0, S1);
        Matrix eMs = Matrix.ones(ms, 1);
        double mu = piMap.mult(S1).mult(eMs).get(0, 0);

        // Utilization
        double rho = lambdaTotal / mu;

        if (rho >= 1.0) {
            System.err.println("Warning: System is unstable (rho = " + rho + " >= 1). Results may be invalid.");
        }

        // Construct M/G/1-type matrices using Kronecker products
        Matrix I_ma = Matrix.eye(ma);
        Matrix I_ms = Matrix.eye(ms);

        // A = [A0, A1, A2, ..., A_{K+1}] for levels >= 1
        Matrix A = Matrix.zeros(m, m * (K + 2));

        // A0 = I_ma \otimes S1 (service completion)
        Matrix A0 = I_ma.kron(S1);
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < m; j++) {
                A.set(i, j, A0.get(i, j));
            }
        }

        // A1 = D0 \otimes I_ms + I_ma \otimes S0 (phase changes only)
        Matrix A1 = D0.kron(I_ms).add(I_ma.kron(S0));
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < m; j++) {
                A.set(i, m + j, A1.get(i, j));
            }
        }

        // A_{k+1} = D_k \otimes I_ms for k >= 1 (batch arrivals)
        for (int k = 1; k <= K; k++) {
            Matrix Ak = D.get(k).kron(I_ms);
            for (int i = 0; i < m; i++) {
                for (int j = 0; j < m; j++) {
                    A.set(i, (k + 1) * m + j, Ak.get(i, j));
                }
            }
        }

        // B = [B0, B1, B2, ..., BK] for level 0 (empty queue)
        Matrix B = Matrix.zeros(m, m * (K + 1));

        // B0: At level 0, no service (add S1 back as self-loop)
        Matrix B0 = D0.kron(I_ms).add(I_ma.kron(S0.add(S1)));
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < m; j++) {
                B.set(i, j, B0.get(i, j));
            }
        }

        // B_k = D_k \otimes I_ms for k >= 1 (batch arrivals from empty queue)
        for (int k = 1; k <= K; k++) {
            Matrix Bk = D.get(k).kron(I_ms);
            for (int i = 0; i < m; i++) {
                for (int j = 0; j < m; j++) {
                    B.set(i, k * m + j, Bk.get(i, j));
                }
            }
        }

        // Compute G matrix using ETAQA
        Matrix G = MG1_ETAQA.mg1_g_etaqa(A);

        // Compute stationary probabilities using ETAQA
        Matrix pi = MG1_ETAQA.mg1_pi_etaqa(B, A, G, null);

        // Compute queue length moments
        double[] qlenMoments = new double[nMoments];
        for (int n = 1; n <= nMoments; n++) {
            qlenMoments[n - 1] = MG1_ETAQA.mg1_qlen_etaqa(B, A, pi, n);
        }

        // Performance metrics
        double meanQueueLength = qlenMoments[0];
        double meanResponseTime = (lambdaTotal > 0) ? meanQueueLength / lambdaTotal : 0.0;

        return new BMAPMAP1Result(
                meanQueueLength,
                rho,
                meanResponseTime,
                lambdaTotal,
                pi,
                G,
                meanBatchSize
        );
    }

    public static BMAPMAP1Result solver_mam_bmap_map_1(List<Matrix> D, Matrix S0, Matrix S1) {
        return solver_mam_bmap_map_1(D, S0, S1, 3);
    }
}
