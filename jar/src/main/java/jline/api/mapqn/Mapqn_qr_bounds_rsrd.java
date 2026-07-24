/**
 * QR Bounds via RSRD Method.
 */
package jline.api.mapqn;

import com.quantego.josqp.CSCMatrix;
import com.quantego.josqp.OSQP;

import java.util.ArrayList;
import java.util.List;

import java.util.HashMap;
import java.util.Map;

public final class Mapqn_qr_bounds_rsrd {

    private Mapqn_qr_bounds_rsrd() {
    }

    public static Mapqn_solution solve(Mapqn_qr_bounds_rsrd_parameters params, int objectiveQueue) {
        return solve(params, objectiveQueue, "min");
    }

    public static Mapqn_solution solve(Mapqn_qr_bounds_rsrd_parameters params, int objectiveQueue, String sense) {
        params.validate();
        if (!(objectiveQueue >= 1 && objectiveQueue <= params.M)) {
            throw new IllegalArgumentException("Objective queue must be in range 1..M");
        }
        if (!("min".equals(sense) || "max".equals(sense))) {
            throw new IllegalArgumentException("Sense must be 'min' or 'max'");
        }

        int M = params.M;
        int N = params.N;
        int[] F = params.F;
        int[] K = params.K;
        Object r = params.r;
        Object[] mu = params.mu;
        Object[] v = params.v;
        double[][] alpha = params.alpha;

        // Compute transition rates q[i][j][ki][hi][n]
        double[][][][][] q = new double[M][M][][][];
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < M; j++) {
                q[i][j] = new double[K[i]][][];
                for (int ki = 0; ki < K[i]; ki++) {
                    q[i][j][ki] = new double[K[i]][];
                    for (int hi = 0; hi < K[i]; hi++) {
                        q[i][j][ki][hi] = new double[N + 1];
                        for (int n = 0; n <= N; n++) {
                            if (n == 0) {
                                q[i][j][ki][hi][n] = 0.0;
                            } else {
                                double alphaVal = (n <= alpha[i].length) ? alpha[i][n - 1] : 1.0;
                                if (j != i) {
                                    q[i][j][ki][hi][n] = ((jline.util.matrix.Matrix) r).get(i, j) *
                                            ((jline.util.matrix.Matrix) mu[i]).get(ki, hi) * alphaVal;
                                } else {
                                    q[i][j][ki][hi][n] = (((jline.util.matrix.Matrix) v[i]).get(ki, hi) +
                                            ((jline.util.matrix.Matrix) r).get(i, i) *
                                                    ((jline.util.matrix.Matrix) mu[i]).get(ki, hi)) * alphaVal;
                                }
                            }
                        }
                    }
                }
            }
        }

        Map<String, Integer> varIndex = new HashMap<String, Integer>();
        int varCount = 0;

        for (int j = 0; j < M; j++) {
            for (int nj = 0; nj <= F[j]; nj++) {
                for (int kj = 0; kj < K[j]; kj++) {
                    for (int i = 0; i < M; i++) {
                        for (int ni = 0; ni <= F[i]; ni++) {
                            for (int hi = 0; hi < K[i]; hi++) {
                                varIndex.put("p2_" + j + "_" + nj + "_" + kj + "_" + i + "_" + ni + "_" + hi, varCount++);
                            }
                        }
                    }
                }
            }
        }

        for (int i = 0; i < M; i++) {
            for (int ki = 0; ki < K[i]; ki++) {
                for (int ni = 1; ni <= F[i]; ni++) {
                    varIndex.put("U_" + i + "_" + ki + "_" + ni, varCount++);
                    varIndex.put("Ueff_" + i + "_" + ki + "_" + ni, varCount++);
                }
            }
        }

        // pb variables: blocking probability per queue
        // (qrboundsrsrd_skel.mod:51 `var pb{i in 1..M} >=0, <=1`)
        for (int i = 0; i < M; i++) {
            varIndex.put("pbvar_" + i, varCount++);
        }

        final int nVars = varCount;
        final Map<String, Integer> finalVarIndex = varIndex;

        // Helper functions implemented as inline lambdas
        IndexFn p2Idx = new IndexFn() {
            @Override
            public int idx(int j, int nj, int kj, int i, int ni, int hi) {
                Integer v = finalVarIndex.get("p2_" + j + "_" + nj + "_" + kj + "_" + i + "_" + ni + "_" + hi);
                return v != null ? v : -1;
            }
        };
        IndexFn3 uIdx = new IndexFn3() {
            @Override
            public int idx(int i, int ki, int ni) {
                Integer v = finalVarIndex.get("U_" + i + "_" + ki + "_" + ni);
                return v != null ? v : -1;
            }
        };
        IndexFn3 ueffIdx = new IndexFn3() {
            @Override
            public int idx(int i, int ki, int ni) {
                Integer v = finalVarIndex.get("Ueff_" + i + "_" + ki + "_" + ni);
                return v != null ? v : -1;
            }
        };

        double[] ub = new double[nVars];
        for (int i = 0; i < nVars; i++) ub[i] = Double.POSITIVE_INFINITY;

        for (int j = 0; j < M; j++) {
            for (int nj = 0; nj <= F[j]; nj++) {
                for (int kj = 0; kj < K[j]; kj++) {
                    for (int i = 0; i < M; i++) {
                        for (int ni = 0; ni <= F[i]; ni++) {
                            for (int hi = 0; hi < K[i]; hi++) {
                                int idx = p2Idx.idx(j, nj, kj, i, ni, hi);
                                if (idx < 0) continue;
                                if (i == j && nj == ni && hi != kj) ub[idx] = 0.0;
                                if (i == j && nj != ni) ub[idx] = 0.0;
                                if (i != j && nj + ni > N) ub[idx] = 0.0;
                                if (i != j) {
                                    int sumOtherF = 0;
                                    for (int y = 0; y < M; y++) {
                                        if (y != i && y != j) sumOtherF += F[y];
                                    }
                                    if (N - nj - ni > sumOtherF) ub[idx] = 0.0;
                                }
                            }
                        }
                    }
                }
                for (int kj = 0; kj < K[j]; kj++) {
                    int sumOtherF = 0;
                    for (int y = 0; y < M; y++) {
                        if (y != j) sumOtherF += F[y];
                    }
                    if (N - nj > sumOtherF) {
                        int idx = p2Idx.idx(j, nj, kj, j, nj, kj);
                        if (idx >= 0) ub[idx] = 0.0;
                    }
                }
            }
        }

        // see _kb/03-api-layer.md for rationale
        List<Row> rows = new ArrayList<Row>();
        int[] vars = new int[nVars];
        for (int idx = 0; idx < nVars; idx++) {
            vars[idx] = idx;
        }

        // see _kb/03-api-layer.md for rationale
        Row zeroRow = null;
        for (int idx = 0; idx < nVars; idx++) {
            if (ub[idx] == 0.0) {
                if (zeroRow == null) zeroRow = newRow(rows, nVars, "fixzero").level(0.0);
                zeroRow.set(idx, 1.0);
            }
        }

        for (int j = 0; j < M; j++) {
            Row expr = newRow(rows, nVars, "norm_" + j).level(1.0);
            for (int nj = 0; nj <= F[j]; nj++) {
                for (int kj = 0; kj < K[j]; kj++) {
                    int idx = p2Idx.idx(j, nj, kj, j, nj, kj);
                    if (idx >= 0) expr.set(vars[idx], 1.0);
                }
            }
        }

        int symCount = 0;
        for (int j = 0; j < M; j++) {
            for (int nj = 0; nj <= F[j]; nj++) {
                for (int kj = 0; kj < K[j]; kj++) {
                    for (int i = j + 1; i < M; i++) {
                        for (int ni = 0; ni <= F[i]; ni++) {
                            for (int hi = 0; hi < K[i]; hi++) {
                                int idx1 = p2Idx.idx(j, nj, kj, i, ni, hi);
                                int idx2 = p2Idx.idx(i, ni, hi, j, nj, kj);
                                if (idx1 >= 0 && idx2 >= 0 && idx1 != idx2 &&
                                        !(ub[idx1] == 0.0 && ub[idx2] == 0.0)) {
                                    Row expr = newRow(rows, nVars, "sym_" + (symCount++)).level(0.0);
                                    expr.set(vars[idx1], 1.0);
                                    expr.set(vars[idx2], -1.0);
                                }
                            }
                        }
                    }
                }
            }
        }

        int margCount = 0;
        for (int j = 0; j < M; j++) {
            for (int kj = 0; kj < K[j]; kj++) {
                for (int nj = 0; nj <= F[j]; nj++) {
                    for (int i = 0; i < M; i++) {
                        if (i == j) continue;
                        Row expr = newRow(rows, nVars, "marg_" + (margCount++)).level(0.0);
                        int idxDiag = p2Idx.idx(j, nj, kj, j, nj, kj);
                        if (idxDiag >= 0) expr.set(vars[idxDiag], 1.0);
                        for (int ni = 0; ni <= F[i]; ni++) {
                            for (int hi = 0; hi < K[i]; hi++) {
                                int idx = p2Idx.idx(j, nj, kj, i, ni, hi);
                                if (idx >= 0) expr.set(vars[idx], -1.0);
                            }
                        }
                    }
                }
            }
        }

        int uclCount = 0;
        for (int i = 0; i < M; i++) {
            for (int ki = 0; ki < K[i]; ki++) {
                for (int ni = 1; ni <= F[i]; ni++) {
                    Row expr = newRow(rows, nVars, "ucl_" + (uclCount++)).level(0.0);
                    int idxU = uIdx.idx(i, ki, ni);
                    int idxP2 = p2Idx.idx(i, ni, ki, i, ni, ki);
                    if (idxU >= 0) expr.set(vars[idxU], 1.0);
                    if (idxP2 >= 0) expr.set(vars[idxP2], -1.0);
                }
            }
        }

        int ueffCount = 0;
        for (int i = 0; i < M; i++) {
            for (int ki = 0; ki < K[i]; ki++) {
                for (int ni = 1; ni <= F[i]; ni++) {
                    Row expr = newRow(rows, nVars, "ueff_" + (ueffCount++)).level(0.0);
                    int idxUeff = ueffIdx.idx(i, ki, ni);
                    int idxP2Diag = p2Idx.idx(i, ni, ki, i, ni, ki);
                    if (idxUeff >= 0) expr.set(vars[idxUeff], 1.0);
                    if (idxP2Diag >= 0) expr.set(vars[idxP2Diag], -1.0);
                    for (int j = 0; j < M; j++) {
                        if (j != i && ((jline.util.matrix.Matrix) r).get(i, j) > 0) {
                            for (int hj = 0; hj < K[j]; hj++) {
                                int idxBlock = p2Idx.idx(i, ni, ki, j, F[j], hj);
                                if (idxBlock >= 0) expr.set(vars[idxBlock], ((jline.util.matrix.Matrix) r).get(i, j));
                            }
                        }
                    }
                }
            }
        }

        // THM2
        int thm2Count = 0;
        for (int i = 0; i < M; i++) {
            for (int ki = 0; ki < K[i]; ki++) {
                double[] coeffs = new double[nVars];
                for (int ni = 1; ni <= F[i]; ni++) {
                    for (int j = 0; j < M; j++) {
                        if (j == i) continue;
                        for (int hi = 0; hi < K[i]; hi++) {
                            if (hi == ki) continue;
                            double coef = q[i][j][ki][hi][ni];
                            int idxUeff = ueffIdx.idx(i, ki, ni);
                            if (idxUeff >= 0) coeffs[idxUeff] += coef;
                        }
                    }
                    for (int hi = 0; hi < K[i]; hi++) {
                        if (hi == ki) continue;
                        double coef = q[i][i][ki][hi][ni];
                        int idxP2 = p2Idx.idx(i, ni, ki, i, ni, ki);
                        if (idxP2 >= 0) coeffs[idxP2] += coef;
                    }
                }
                for (int ni = 1; ni <= F[i]; ni++) {
                    for (int j = 0; j < M; j++) {
                        if (j == i) continue;
                        for (int hi = 0; hi < K[i]; hi++) {
                            if (hi == ki) continue;
                            double coef = q[i][j][hi][ki][ni];
                            int idxUeff = ueffIdx.idx(i, hi, ni);
                            if (idxUeff >= 0) coeffs[idxUeff] -= coef;
                        }
                    }
                    for (int hi = 0; hi < K[i]; hi++) {
                        if (hi == ki) continue;
                        double coef = q[i][i][hi][ki][ni];
                        int idxP2 = p2Idx.idx(i, ni, hi, i, ni, hi);
                        if (idxP2 >= 0) coeffs[idxP2] -= coef;
                    }
                }
                Row expr = newRow(rows, nVars, "thm2_" + (thm2Count++)).level(0.0);
                emitEquilibrated(expr, vars, coeffs, nVars);
            }
        }

        // see _kb/03-api-layer.md for rationale
        int thm1Count = 0;
        for (int j = 0; j < M; j++) {
            for (int kj = 0; kj < K[j]; kj++) {
                double[] coeffs = new double[nVars];
                for (int nj = 1; nj <= F[j]; nj++) {
                    int idxDiag = p2Idx.idx(j, nj, kj, j, nj, kj);
                    if (idxDiag >= 0) coeffs[idxDiag] -= (double) N;
                    for (int i = 0; i < M; i++) {
                        for (int ni = 1; ni <= F[i]; ni++) {
                            for (int hi = 0; hi < K[i]; hi++) {
                                int idx = p2Idx.idx(j, nj, kj, i, ni, hi);
                                if (idx >= 0) coeffs[idx] += (double) ni;
                            }
                        }
                    }
                }
                Row expr = newRow(rows, nVars, "thm1_" + (thm1Count++)).level(0.0);
                emitEquilibrated(expr, vars, coeffs, nVars);
            }
        }

        // see _kb/03-api-layer.md for rationale
        int thm1cCount = 0;
        for (int j = 0; j < M; j++) {
            for (int kj = 0; kj < K[j]; kj++) {
                double[] coeffs = new double[nVars];
                int idxDiag = p2Idx.idx(j, 0, kj, j, 0, kj);
                if (idxDiag >= 0) coeffs[idxDiag] -= (double) N;
                for (int i = 0; i < M; i++) {
                    if (i == j) continue;
                    for (int ni = 1; ni <= F[i]; ni++) {
                        for (int hi = 0; hi < K[i]; hi++) {
                            int idx = p2Idx.idx(j, 0, kj, i, ni, hi);
                            if (idx >= 0) coeffs[idx] += (double) ni;
                        }
                    }
                }
                Row expr = newRow(rows, nVars, "thm1c_" + (thm1cCount++)).level(0.0);
                emitEquilibrated(expr, vars, coeffs, nVars);
            }
        }

        // PBLOCK: pb[i] = sum_{k, ni>=1} (U[i,k,ni] - Ueff[i,k,ni])
        // (qrboundsrsrd_skel.mod:52)
        for (int i = 0; i < M; i++) {
            Integer pbv = finalVarIndex.get("pbvar_" + i);
            if (pbv == null) continue;
            double[] coeffs = new double[nVars];
            coeffs[pbv] += 1.0;
            for (int ki = 0; ki < K[i]; ki++) {
                for (int ni = 1; ni <= F[i]; ni++) {
                    int idxU = uIdx.idx(i, ki, ni);
                    int idxUeff = ueffIdx.idx(i, ki, ni);
                    if (idxU >= 0) coeffs[idxU] -= 1.0;
                    if (idxUeff >= 0) coeffs[idxUeff] += 1.0;
                }
            }
            Row expr = newRow(rows, nVars, "pblock_" + i).level(0.0);
            emitEquilibrated(expr, vars, coeffs, nVars);
        }

        // PBB: pb[i] <= sum_{j!=i, r[i,j]>0} sum_h p2[j,F[j],h,j,F[j],h]
        // (qrboundsrsrd_skel.mod:53)
        for (int i = 0; i < M; i++) {
            Integer pbv = finalVarIndex.get("pbvar_" + i);
            if (pbv == null) continue;
            double[] coeffs = new double[nVars];
            coeffs[pbv] += 1.0;
            for (int j = 0; j < M; j++) {
                if (j == i || ((jline.util.matrix.Matrix) r).get(i, j) <= 0) continue;
                for (int hj = 0; hj < K[j]; hj++) {
                    int idx = p2Idx.idx(j, F[j], hj, j, F[j], hj);
                    if (idx >= 0) coeffs[idx] -= 1.0;
                }
            }
            Row expr = newRow(rows, nVars, "pbb_" + i).upper(0.0);
            emitEquilibrated(expr, vars, coeffs, nVars);
        }

        // THM3a
        int thm3aCount = 0;
        for (int i = 0; i < M; i++) {
            for (int ni = 1; ni < F[i]; ni++) {
                double[] coeffs = new double[nVars];
                for (int j = 0; j < M; j++) {
                    if (j == i) continue;
                    for (int kj = 0; kj < K[j]; kj++) {
                        for (int hj = 0; hj < K[j]; hj++) {
                            for (int ui = 0; ui < K[i]; ui++) {
                                for (int nj = 1; nj <= F[j]; nj++) {
                                    int idx = p2Idx.idx(j, nj, kj, i, ni, ui);
                                    if (idx >= 0) coeffs[idx] += q[j][i][kj][hj][nj];
                                }
                            }
                        }
                    }
                }
                for (int j = 0; j < M; j++) {
                    if (j == i) continue;
                    for (int ki = 0; ki < K[i]; ki++) {
                        for (int hi = 0; hi < K[i]; hi++) {
                            for (int uj = 0; uj < K[j]; uj++) {
                                for (int nj = 0; nj < F[j]; nj++) {
                                    int idx = p2Idx.idx(i, ni + 1, ki, j, nj, uj);
                                    if (idx >= 0 && ni + 1 <= F[i]) {
                                        coeffs[idx] -= q[i][j][ki][hi][ni + 1];
                                    }
                                }
                            }
                        }
                    }
                }
                Row expr = newRow(rows, nVars, "thm3a_" + (thm3aCount++)).level(0.0);
                emitEquilibrated(expr, vars, coeffs, nVars);
            }
        }

        // THM3b
        int thm3bCount = 0;
        for (int i = 0; i < M; i++) {
            for (int ui = 0; ui < K[i]; ui++) {
                double[] coeffs = new double[nVars];
                for (int j = 0; j < M; j++) {
                    if (j == i) continue;
                    for (int kj = 0; kj < K[j]; kj++) {
                        for (int hj = 0; hj < K[j]; hj++) {
                            for (int nj = 1; nj <= F[j]; nj++) {
                                int idx = p2Idx.idx(j, nj, kj, i, 0, ui);
                                if (idx >= 0) coeffs[idx] += q[j][i][kj][hj][nj];
                            }
                        }
                    }
                }
                for (int j = 0; j < M; j++) {
                    if (j == i) continue;
                    for (int ki = 0; ki < K[i]; ki++) {
                        for (int nj = 0; nj < F[j]; nj++) {
                            for (int hj = 0; hj < K[j]; hj++) {
                                int idx = p2Idx.idx(i, 1, ki, j, nj, hj);
                                if (idx >= 0) coeffs[idx] -= q[i][j][ki][ui][1];
                            }
                        }
                    }
                }
                Row expr = newRow(rows, nVars, "thm3b_" + (thm3bCount++)).level(0.0);
                emitEquilibrated(expr, vars, coeffs, nVars);
            }
        }

        // QBAL
        for (int i = 0; i < M; i++) {
            for (int ki = 0; ki < K[i]; ki++) {
                double[] coeffs = new double[nVars];
                for (int hi = 0; hi < K[i]; hi++) {
                    if (hi == ki) continue;
                    for (int j = 0; j < M; j++) {
                        if (j == i) continue;
                        for (int ni = 1; ni <= F[i]; ni++) {
                            for (int uj = 0; uj < K[j]; uj++) {
                                for (int nj = 0; nj < F[j]; nj++) {
                                    double coef = q[i][j][ki][hi][ni] * ni;
                                    int idx = p2Idx.idx(i, ni, ki, j, nj, uj);
                                    if (idx >= 0) coeffs[idx] += coef;
                                }
                            }
                        }
                    }
                }
                for (int hi = 0; hi < K[i]; hi++) {
                    if (hi == ki) continue;
                    for (int ni = 1; ni <= F[i]; ni++) {
                        double coef = q[i][i][ki][hi][ni] * ni;
                        int idx = p2Idx.idx(i, ni, ki, i, ni, ki);
                        if (idx >= 0) coeffs[idx] += coef;
                    }
                }
                for (int j = 0; j < M; j++) {
                    if (j == i) continue;
                    for (int hi = 0; hi < K[i]; hi++) {
                        for (int ni = 1; ni <= F[i]; ni++) {
                            for (int uj = 0; uj < K[j]; uj++) {
                                int maxNj = Math.min(F[j] - 1, N - ni);
                                for (int nj = 0; nj <= maxNj; nj++) {
                                    double coef = q[i][j][hi][ki][ni];
                                    int idx = p2Idx.idx(i, ni, hi, j, nj, uj);
                                    if (idx >= 0) coeffs[idx] += coef;
                                }
                            }
                        }
                    }
                }
                for (int j = 0; j < M; j++) {
                    if (j == i) continue;
                    for (int hj = 0; hj < K[j]; hj++) {
                        for (int ni = 0; ni < F[i]; ni++) {
                            for (int uj = 0; uj < K[j]; uj++) {
                                for (int nj = 1; nj <= F[j]; nj++) {
                                    double coef = q[j][i][hj][uj][nj];
                                    int idx = p2Idx.idx(i, ni, ki, j, nj, hj);
                                    if (idx >= 0) coeffs[idx] -= coef;
                                }
                            }
                        }
                    }
                }
                for (int hi = 0; hi < K[i]; hi++) {
                    if (hi == ki) continue;
                    for (int ni = 1; ni <= F[i]; ni++) {
                        double coef = q[i][i][hi][ki][ni] * ni;
                        int idx = p2Idx.idx(i, ni, hi, i, ni, hi);
                        if (idx >= 0) coeffs[idx] -= coef;
                    }
                }
                for (int hi = 0; hi < K[i]; hi++) {
                    if (hi == ki) continue;
                    for (int j = 0; j < M; j++) {
                        if (j == i) continue;
                        for (int ni = 1; ni <= F[i]; ni++) {
                            for (int uj = 0; uj < K[j]; uj++) {
                                for (int nj = 0; nj < F[j]; nj++) {
                                    double coef = q[i][j][hi][ki][ni] * ni;
                                    int idx = p2Idx.idx(i, ni, hi, j, nj, uj);
                                    if (idx >= 0) coeffs[idx] -= coef;
                                }
                            }
                        }
                    }
                }
                Row expr = newRow(rows, nVars, "qbal_" + i + "_" + ki).level(0.0);
                emitEquilibrated(expr, vars, coeffs, nVars);
            }
        }

        // THM4
        int thm4Count = 0;
        for (int j = 0; j < M; j++) {
            for (int kj = 0; kj < K[j]; kj++) {
                for (int i = 0; i < M; i++) {
                    double[] coeffs = new double[nVars];
                    for (int t = 0; t < M; t++) {
                        for (int ht = 0; ht < K[t]; ht++) {
                            for (int nj = 0; nj <= F[j]; nj++) {
                                for (int nt = 1; nt <= F[t]; nt++) {
                                    int idx = p2Idx.idx(j, nj, kj, t, nt, ht);
                                    if (idx >= 0) coeffs[idx] += (double) nt;
                                }
                            }
                        }
                    }
                    for (int hi = 0; hi < K[i]; hi++) {
                        for (int nj = 0; nj <= F[j]; nj++) {
                            for (int ni = 1; ni <= F[i]; ni++) {
                                int idx = p2Idx.idx(j, nj, kj, i, ni, hi);
                                if (idx >= 0) coeffs[idx] -= (double) N;
                            }
                        }
                    }
                    Row expr = newRow(rows, nVars, "thm4_" + (thm4Count++)).lower(0.0);
                    emitEquilibrated(expr, vars, coeffs, nVars);
                }
            }
        }

        double[] objectiveCoeffs = new double[nVars];
        int targetQueue = objectiveQueue - 1;
        for (int ki = 0; ki < K[targetQueue]; ki++) {
            for (int ni = 1; ni <= F[targetQueue]; ni++) {
                int idx = p2Idx.idx(targetQueue, ni, ki, targetQueue, ni, ki);
                if (idx >= 0) objectiveCoeffs[idx] = 1.0;
            }
        }

        // see _kb/03-api-layer.md for rationale
        String dumpPath = System.getProperty("jline.mapqn.dumpmps");
        if (dumpPath != null) {
            try {
                java.io.PrintWriter pw = new java.io.PrintWriter(dumpPath);
                pw.println("NVARS " + nVars);
                for (int idx = 0; idx < nVars; idx++) {
                    pw.println("BND " + idx + " 0.0 " + ((ub[idx] == 0.0) ? "0.0" : "1.0"));
                    if (objectiveCoeffs[idx] != 0.0) {
                        pw.println("OBJ " + idx + " " + objectiveCoeffs[idx]);
                    }
                }
                for (int r2 = 0; r2 < rows.size(); r2++) {
                    Row row = rows.get(r2);
                    pw.println("ROW " + row.name
                            + " " + (row.hasLower ? Double.toString(row.lower) : "NONE")
                            + " " + (row.hasUpper ? Double.toString(row.upper) : "NONE"));
                    for (int idx = 0; idx < nVars; idx++) {
                        if (row.coeffs[idx] != 0.0) {
                            pw.println("COEF " + idx + " " + row.coeffs[idx]);
                        }
                    }
                }
                pw.close();
                System.err.println("[mapqn_qr_bounds_rsrd] wrote LP dump to " + dumpPath);
            } catch (Exception dumpEx) {
                System.err.println("[mapqn_qr_bounds_rsrd] LP dump failed: " + dumpEx);
            }
        }

        // see _kb/03-api-layer.md for rationale
        int nRows = 0;
        for (int r2 = 0; r2 < rows.size(); r2++) {
            if (rows.get(r2).nonZeroCount() > 0) nRows++;
        }
        int mRows = nRows + nVars;

        List<List<int[]>> byCol = new ArrayList<List<int[]>>();
        List<List<Double>> valByCol = new ArrayList<List<Double>>();
        for (int j = 0; j < nVars; j++) {
            byCol.add(new ArrayList<int[]>());
            valByCol.add(new ArrayList<Double>());
        }
        double[] rowLower = new double[mRows];
        double[] rowUpper = new double[mRows];
        int ri = 0;
        for (int r2 = 0; r2 < rows.size(); r2++) {
            Row row = rows.get(r2);
            if (row.nonZeroCount() == 0) continue;
            rowLower[ri] = row.hasLower ? row.lower : -OSQP.OSQP_INFTY;
            rowUpper[ri] = row.hasUpper ? row.upper : OSQP.OSQP_INFTY;
            for (int idx = 0; idx < nVars; idx++) {
                if (row.coeffs[idx] != 0.0) {
                    byCol.get(idx).add(new int[]{ri});
                    valByCol.get(idx).add(Double.valueOf(row.coeffs[idx]));
                }
            }
            ri++;
        }
        for (int j = 0; j < nVars; j++) {
            byCol.get(j).add(new int[]{ri});
            valByCol.get(j).add(Double.valueOf(1.0));
            rowLower[ri] = 0.0;
            rowUpper[ri] = (ub[j] == 0.0) ? 0.0 : 1.0;
            ri++;
        }

        int nnz = 0;
        for (int j = 0; j < nVars; j++) nnz += byCol.get(j).size();
        int[] Ap = new int[nVars + 1];
        int[] Ai = new int[nnz];
        double[] Ax = new double[nnz];
        int pos = 0;
        for (int j = 0; j < nVars; j++) {
            Ap[j] = pos;
            List<int[]> colRows = byCol.get(j);
            List<Double> colVals = valByCol.get(j);
            for (int k = 0; k < colRows.size(); k++) {
                Ai[pos] = colRows.get(k)[0];
                Ax[pos] = colVals.get(k).doubleValue();
                pos++;
            }
        }
        Ap[nVars] = pos;
        CSCMatrix A = new CSCMatrix(mRows, nVars, nnz, Ap, Ai, Ax);
        CSCMatrix P = new CSCMatrix(nVars, nVars, 0, new int[nVars + 1], new int[0], new double[0]);

        // OSQP always minimises, so a maximisation is posed with a negated
        // linear cost and the objective is evaluated from the primal point.
        boolean maximise = "max".equals(sense);
        double[] linearCost = new double[nVars];
        for (int j = 0; j < nVars; j++) {
            linearCost[j] = maximise ? -objectiveCoeffs[j] : objectiveCoeffs[j];
        }

        try {
            OSQP.Settings settings = new OSQP.Settings();
            settings.max_iter = 200000;
            settings.eps_abs = OSQP_EPS;
            settings.eps_rel = OSQP_EPS;
            settings.polish = true;
            settings.verbose = false;

            OSQP solver = new OSQP(
                    new OSQP.Data(nVars, mRows, P, A, linearCost, rowLower, rowUpper, 0.0),
                    settings);
            OSQP.Status status = solver.solve();

            if (System.getProperty("jline.mapqn.debug") != null) {
                System.err.println("[mapqn_qr_bounds_rsrd] josqp status=" + status);
            }

            if (status != OSQP.Status.SOLVED && status != OSQP.Status.SOLVED_INACCURATE) {
                return new Mapqn_solution(Double.NaN, new HashMap<String, Double>());
            }

            double[] point = extractPrimal(solver, nVars);
            double objective = 0.0;
            for (int j = 0; j < nVars; j++) {
                objective += objectiveCoeffs[j] * point[j];
            }

            Map<String, Double> variables = new HashMap<String, Double>();
            for (int i = 0; i < M; i++) {
                double totalU = 0.0;
                double totalUeff = 0.0;
                for (int ki = 0; ki < K[i]; ki++) {
                    for (int ni = 1; ni <= F[i]; ni++) {
                        int idxU = uIdx.idx(i, ki, ni);
                        int idxUeff = ueffIdx.idx(i, ki, ni);
                        if (idxU >= 0) totalU += point[idxU];
                        if (idxUeff >= 0) totalUeff += point[idxUeff];
                    }
                }
                variables.put("U_" + (i + 1), totalU);
                variables.put("Ueff_" + (i + 1), totalUeff);
                variables.put("pb_" + (i + 1), totalU - totalUeff);
            }

            return new Mapqn_solution(objective, variables);
        } catch (Exception e) {
            if (System.getProperty("jline.mapqn.debug") != null) {
                System.err.println("[mapqn_qr_bounds_rsrd] exception: " + e);
            }
            return new Mapqn_solution(Double.NaN, new HashMap<String, Double>());
        }
    }

    /**
     * ADMM termination tolerance for the josqp backend, used for both eps_abs
     * and eps_rel.
     *
     * <p>TIGHTER IS NOT BETTER HERE, and the failure is not graceful. Measured
     * over eps in {1e-4 .. 1e-10} x rho in {1e-2, 1e-1, 1} x scaling in {10, 25}
     * on the two paper instances (objective error d against the AMPL optimum,
     * which glpsol reproduces from this method's own LP dump):
     *
     * <pre>
     *   eps    sec7.1 (0.7332453660)          M5N20 (0.8705798470)
     *   1e-4   SOLVED     d=3.6e-05  0.1 s    SOLVED     d=1.9e-03  0.3 s
     *   1e-5   SOLVED     d=2.1e-06  0.5 s    SOLVED     d=8.7e-05  0.4 s
     *   1e-6   MAX_ITER   (NaN)       23 s    SOLVED     d=6.4e-06  1.6 s
     *   1e-8   MAX_ITER   (NaN)       23 s    MAX_ITER   (NaN)       51 s
     * </pre>
     *
     * <p>1e-5 is the only value that converges on BOTH. 1e-6 is more accurate on
     * M5N20 but saturates 200000 iterations on the SMALLER sec7.1 and returns
     * nothing, so it cannot be the default. rho barely matters (its default is
     * already 1e-1) and scaling 10 vs 25 makes no measurable difference.
     *
     * <p>The residual ~1e-4 relative error is inherent to ADMM on this
     * degenerate LP; the polish step does not recover the exact vertex. MATLAB
     * (matlab/lib/qrf/qrf_rsrd.m) and native Python (api/mapqn/qr_bounds_rsrd.py)
     * agree with AMPL to ~2e-11, so the JAR bound is NOT interchangeable with
     * theirs beyond about four significant decimals. Re-measure both instances
     * before changing this.
     */
    private static final double OSQP_EPS = 1.0e-5;

    /**
     * Read the primal solution out of a solved josqp workspace.
     *
     * <p>josqp 0.6.5 exposes {@code Workspace.solution} and {@code Solution.x}
     * as package-private, with no accessor, so reflection is the only route to
     * the primal point from outside {@code com.quantego.josqp}.
     */
    private static double[] extractPrimal(OSQP solver, int nVars) throws Exception {
        Object workspace = solver.getWorkspace();
        java.lang.reflect.Field solutionField = workspace.getClass().getDeclaredField("solution");
        solutionField.setAccessible(true);
        Object solution = solutionField.get(workspace);
        java.lang.reflect.Field xField = solution.getClass().getDeclaredField("x");
        xField.setAccessible(true);
        double[] x = (double[]) xField.get(solution);
        double[] copy = new double[nVars];
        System.arraycopy(x, 0, copy, 0, nVars);
        return copy;
    }

    /**
     * Relative magnitude below which an accumulated coefficient is treated as a
     * structural zero rather than a genuine term. See {@link #emitEquilibrated}.
     */
    private static final double DUST_REL = 1.0e-12;

    /**
     * One assembled constraint row: a dense coefficient vector plus its bounds.
     *
     * <p>The builders below accumulate a row's left- and right-hand side terms into
     * a single coefficient array, adding the LHS blocks and subtracting the RHS
     * blocks, then hand it here. Collecting rows in this form keeps the assembly
     * independent of the LP backend and is what let the assembled program be dumped
     * and checked against glpsol.
     */
    private static final class Row {
        final String name;
        final double[] coeffs;
        double lower;
        double upper;
        boolean hasLower;
        boolean hasUpper;

        Row(String name, int nVars) {
            this.name = name;
            this.coeffs = new double[nVars];
        }

        Row level(double value) {
            lower = value;
            upper = value;
            hasLower = true;
            hasUpper = true;
            return this;
        }

        Row lower(double value) {
            lower = value;
            hasLower = true;
            return this;
        }

        Row upper(double value) {
            upper = value;
            hasUpper = true;
            return this;
        }

        void set(int idx, double value) {
            coeffs[idx] = value;
        }

        int nonZeroCount() {
            int count = 0;
            for (int idx = 0; idx < coeffs.length; idx++) {
                if (coeffs[idx] != 0.0) count++;
            }
            return count;
        }
    }

    private static Row newRow(List<Row> rows, int nVars, String name) {
        Row row = new Row(name, nVars);
        rows.add(row);
        return row;
    }

    /**
     * Emit an accumulated constraint row, row-equilibrated so that the
     * largest-magnitude coefficient becomes 1, with cancellation dust discarded.
     *
     * <p>Many terms appear on both sides of the balance equations and must cancel
     * exactly. Because the blocks are accumulated in different loop orders, and
     * because a block interleaves writes to an index a later block also touches,
     * floating-point addition is not associative and the cancellation leaves a
     * residue of a few ulp instead of a clean zero. On the paper's M=5, N=20
     * instance the QBAL rows carry 80 such residues at 1e-16..1e-18 while every
     * genuine coefficient in the same rows is at or above 1e-4 -- a gap of twelve
     * decades with nothing in between, which makes {@link #DUST_REL} unambiguous:
     * far above the roundoff, far below the smallest real term. The smallest real
     * term is bounded by the spread of the mu entries (1.016186 down to
     * 2.585708e-05, about 4e4).
     *
     * <p>Emitting the residue turns a structural zero of the constraint matrix into
     * a nonzero and inflates the coefficient ratio of the assembled program by
     * orders of magnitude. Discarding it restores the sparsity pattern the
     * formulation actually has; with it discarded, glpsol reports a ratio of
     * 3.9e4 and solves the M5N20 program to the AMPL optimum 0.8705798470 with no
     * instability warning.
     *
     * <p>Scaling a constraint by a positive constant leaves its solution set
     * unchanged PROVIDED the right-hand side is scaled identically. Every caller of
     * this helper emits a row whose RHS is ZERO -- {@code level(0.0)},
     * {@code upper(0.0)} or {@code lower(0.0)} -- for which the RHS is invariant
     * under scaling, so this is exact rather than an approximation. Do NOT use this
     * helper for a row with a non-zero RHS (such as the {@code norm_j}
     * normalisation at {@code level(1.0)}) without scaling that RHS too; that row
     * already has unit coefficients and needs no equilibration.
     */
    private static void emitEquilibrated(Row expr, int[] vars, double[] coeffs, int nVars) {
        double maxAbs = 0.0;
        for (int idx = 0; idx < nVars; idx++) {
            double a = Math.abs(coeffs[idx]);
            if (a > maxAbs) {
                maxAbs = a;
            }
        }
        double scale = (maxAbs > 0.0) ? (1.0 / maxAbs) : 1.0;
        double dust = maxAbs * DUST_REL;
        for (int idx = 0; idx < nVars; idx++) {
            if (Math.abs(coeffs[idx]) > dust) {
                expr.set(vars[idx], coeffs[idx] * scale);
            }
        }
    }

    private interface IndexFn {
        int idx(int j, int nj, int kj, int i, int ni, int hi);
    }

    private interface IndexFn3 {
        int idx(int i, int ki, int ni);
    }
}
