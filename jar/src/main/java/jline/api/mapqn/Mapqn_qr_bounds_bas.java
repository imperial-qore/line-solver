/**
 * QR Bounds via BAS Method.
 */
package jline.api.mapqn;

import com.quantego.josqp.CSCMatrix;
import com.quantego.josqp.OSQP;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Quadratic-reduction (QR) utilization bounds for closed MAP queueing networks
 * under the blocking-after-service (BAS) protocol.
 *
 * <p>The authoritative statement of this linear program is the AMPL/GMPL model
 * {@code qrboundsbas_skel.mod} (instantiated by {@code example_bas.mod}), and
 * every constraint block below names the AMPL constraint it implements. When
 * the two disagree, the AMPL model is right. On the paper's M=5, N=10 instance
 * glpsol solves that model to U1min = 0.4776625414 and U1max = 0.8077316296.
 */
public final class Mapqn_qr_bounds_bas {

    private Mapqn_qr_bounds_bas() {
    }

    public static Mapqn_solution solve(Mapqn_qr_bounds_bas_parameters params, int objectiveQueue) {
        return solve(params, objectiveQueue, "min");
    }

    public static Mapqn_solution solve(Mapqn_qr_bounds_bas_parameters params, int objectiveQueue, String sense) {
        params.validate();
        if (!(objectiveQueue >= 1 && objectiveQueue <= params.M)) {
            throw new IllegalArgumentException("Objective queue must be in range 1..M");
        }
        if (!("min".equals(sense) || "max".equals(sense))) {
            throw new IllegalArgumentException("Sense must be 'min' or 'max'");
        }

        int M = params.M;
        int N = params.N;
        int f = params.f - 1;
        int[] F = params.F;
        int[] K = params.K;
        int MR = params.MR;
        jline.util.matrix.Matrix BB = params.BB;
        jline.util.matrix.Matrix MM = params.MM;
        int[] ZZ = params.ZZ;
        int ZM = params.ZM;
        jline.util.matrix.Matrix MM1 = params.MM1;

        jline.util.matrix.Matrix[][] q = params.q();

        Map<String, Integer> varIndex = new HashMap<String, Integer>();
        int varCount = 0;

        for (int j = 0; j < M; j++) {
            for (int nj = 0; nj <= N; nj++) {
                for (int kj = 0; kj < K[j]; kj++) {
                    for (int i = 0; i < M; i++) {
                        for (int ni = 0; ni <= N; ni++) {
                            for (int hi = 0; hi < K[i]; hi++) {
                                for (int m = 0; m < MR; m++) {
                                    varIndex.put("p2_" + j + "_" + nj + "_" + kj + "_" + i + "_" + ni + "_" + hi + "_" + m, varCount++);
                                }
                            }
                        }
                    }
                }
            }
        }

        for (int i = 0; i < M; i++) {
            for (int ki = 0; ki < K[i]; ki++) {
                varIndex.put("e_" + i + "_" + ki, varCount++);
            }
        }

        final int nVars = varCount;
        final Map<String, Integer> finalVarIndex = varIndex;

        IndexFn7 p2Idx = new IndexFn7() {
            @Override
            public int idx(int j, int nj, int kj, int i, int ni, int hi, int m) {
                Integer val = finalVarIndex.get("p2_" + j + "_" + nj + "_" + kj + "_" + i + "_" + ni + "_" + hi + "_" + m);
                return val != null ? val : -1;
            }
        };
        IndexFn2 eIdx = new IndexFn2() {
            @Override
            public int idx(int i, int ki) {
                Integer val = finalVarIndex.get("e_" + i + "_" + ki);
                return val != null ? val : -1;
            }
        };

        double[] ub = new double[nVars];
        for (int i = 0; i < nVars; i++) ub[i] = 1.0;

        // see _kb/03-api-layer.md for rationale
        for (int j = 0; j < M; j++) {
            for (int nj = 0; nj <= N; nj++) {
                for (int kj = 0; kj < K[j]; kj++) {
                    for (int i = 0; i < M; i++) {
                        for (int ni = 0; ni <= N; ni++) {
                            for (int hi = 0; hi < K[i]; hi++) {
                                for (int m = 0; m < MR; m++) {
                                    int idx = p2Idx.idx(j, nj, kj, i, ni, hi, m);
                                    if (idx < 0) continue;

                                    if (i == j && nj == ni && hi != kj) ub[idx] = 0.0;
                                    if (i == j && nj != ni) ub[idx] = 0.0;
                                    if (i != j && nj + ni > N) ub[idx] = 0.0;
                                    // see _kb/03-api-layer.md for rationale
                                    if (nj > F[j] || ni > F[i]) ub[idx] = 0.0;
                                    if (m >= 1 && (int) BB.get(m, j) == 1 && nj == 0) ub[idx] = 0.0;
                                    if (m >= 1 && (int) BB.get(m, j) == 1 && i != j && i != f && ni + nj + F[f] > N) ub[idx] = 0.0;
                                    if (j == f && nj >= 1 && nj <= F[f] - 1 && m >= 1) ub[idx] = 0.0;
                                }
                            }
                        }
                    }
                }
            }
        }

        // ZERO4 (qrboundsbas_skel.mod:37): a sum of non-negative variables equal
        // to zero, hence each term is zero.
        for (int j = 0; j < M; j++) {
            if (j == f) continue;
            for (int nj = 0; nj <= N; nj++) {
                for (int kj = 0; kj < K[j]; kj++) {
                    for (int m = 1; m < MR; m++) {
                        for (int nf = 0; nf < F[f]; nf++) {
                            for (int hf = 0; hf < K[f]; hf++) {
                                int idx = p2Idx.idx(j, nj, kj, f, nf, hf, m);
                                if (idx >= 0) ub[idx] = 0.0;
                            }
                        }
                    }
                }
            }
        }

        // see _kb/03-api-layer.md for rationale
        int[] colMap = new int[nVars];
        int nCols = 0;
        for (int idx = 0; idx < nVars; idx++) {
            colMap[idx] = (ub[idx] == 0.0) ? -1 : nCols++;
        }

        List<Row> rows = new ArrayList<Row>();
        Acc acc = new Acc(colMap, nCols);

        // ONE (skel:33): normalisation of each queue's marginal.
        for (int j = 0; j < M; j++) {
            for (int nj = 0; nj <= N; nj++) {
                for (int kj = 0; kj < K[j]; kj++) {
                    for (int m = 0; m < MR; m++) {
                        acc.add(p2Idx.idx(j, nj, kj, j, nj, kj, m), 1.0);
                    }
                }
            }
            emit(rows, "one_" + j, acc, 1.0, 1.0, false);
        }

        // SIMMETRY (skel:42). Restricted to nj <= F[j], ni <= F[i] and nj+ni <= N;
        // outside that box both members of the pair are pinned to zero.
        int symCount = 0;
        for (int j = 0; j < M; j++) {
            for (int nj = 0; nj <= Math.min(N, F[j]); nj++) {
                for (int kj = 0; kj < K[j]; kj++) {
                    for (int i = j + 1; i < M; i++) {
                        for (int ni = 0; ni <= Math.min(N, F[i]); ni++) {
                            if (nj + ni > N) continue;
                            for (int hi = 0; hi < K[i]; hi++) {
                                for (int m = 0; m < MR; m++) {
                                    int idx1 = p2Idx.idx(j, nj, kj, i, ni, hi, m);
                                    int idx2 = p2Idx.idx(i, ni, hi, j, nj, kj, m);
                                    if (idx1 < 0 || idx2 < 0 || idx1 == idx2) continue;
                                    if (ub[idx1] == 0.0 && ub[idx2] == 0.0) continue;
                                    acc.add(idx1, 1.0);
                                    acc.add(idx2, -1.0);
                                    emit(rows, "sym_" + (symCount++), acc, 0.0, 0.0, false);
                                }
                            }
                        }
                    }
                }
            }
        }

        // see _kb/03-api-layer.md for rationale
        int margCount = 0;
        for (int j = 0; j < M; j++) {
            for (int kj = 0; kj < K[j]; kj++) {
                for (int nj = 0; nj <= Math.min(N, F[j]); nj++) {
                    for (int i = 0; i < M; i++) {
                        if (i == j) continue;
                        for (int m = 0; m < MR; m++) {
                            acc.add(p2Idx.idx(j, nj, kj, j, nj, kj, m), 1.0);
                            for (int ni = 0; ni <= Math.min(N - nj, F[i]); ni++) {
                                for (int hi = 0; hi < K[i]; hi++) {
                                    acc.add(p2Idx.idx(j, nj, kj, i, ni, hi, m), -1.0);
                                }
                            }
                            emit(rows, "marg_" + (margCount++), acc, 0.0, 0.0, false);
                        }
                    }
                }
            }
        }

        // see _kb/03-api-layer.md for rationale
        for (int j = 0; j < M; j++) {
            for (int i = 0; i < M; i++) {
                for (int ki = 0; ki < K[i]; ki++) {
                    acc.add(eIdx.idx(i, ki), -1.0);
                    for (int nj = 0; nj <= Math.min(N, F[j]); nj++) {
                        for (int kj = 0; kj < K[j]; kj++) {
                            for (int m = 0; m < MR; m++) {
                                if ((int) BB.get(m, i) != 0) continue;
                                for (int ni = 1; ni <= Math.min(N, F[i]); ni++) {
                                    acc.add(p2Idx.idx(j, nj, kj, i, ni, ki, m), 1.0);
                                }
                            }
                        }
                    }
                    emit(rows, "ueff_" + j + "_" + i + "_" + ki, acc, 0.0, 0.0, false);
                }
            }
        }

        // see _kb/03-api-layer.md for rationale
        for (int i = 0; i < M; i++) {
            for (int ki = 0; ki < K[i]; ki++) {
                for (int j = 0; j < M; j++) {
                    for (int hi = 0; hi < K[i]; hi++) {
                        if (j == i && hi == ki) continue;
                        acc.add(eIdx.idx(i, ki), q[i][j].get(ki, hi));
                        acc.add(eIdx.idx(i, hi), -q[i][j].get(hi, ki));
                    }
                }
                emit(rows, "thm1_" + i + "_" + ki, acc, 0.0, 0.0, true);
            }
        }

        // THM2 (skel:50): conditional population balance.
        for (int j = 0; j < M; j++) {
            for (int kj = 0; kj < K[j]; kj++) {
                for (int nj = 0; nj <= F[j]; nj++) {
                    for (int m = 0; m < MR; m++) {
                        acc.add(p2Idx.idx(j, nj, kj, j, nj, kj, m), -(double) N);
                        for (int i = 0; i < M; i++) {
                            for (int ni = 1; ni <= F[i]; ni++) {
                                for (int ki = 0; ki < K[i]; ki++) {
                                    acc.add(p2Idx.idx(j, nj, kj, i, ni, ki, m), (double) ni);
                                }
                            }
                        }
                        emit(rows, "thm2_" + j + "_" + kj + "_" + nj + "_" + m, acc, 0.0, 0.0, false);
                    }
                }
            }
        }

        // COR1 (skel:51): second moment of the joint population.
        for (int m = 0; m < MR; m++) {
            for (int i = 0; i < M; i++) {
                for (int j = 0; j < M; j++) {
                    for (int nj = 1; nj <= F[j]; nj++) {
                        for (int ni = 1; ni <= F[i]; ni++) {
                            for (int ki = 0; ki < K[i]; ki++) {
                                for (int kj = 0; kj < K[j]; kj++) {
                                    acc.add(p2Idx.idx(j, nj, kj, i, ni, ki, m), (double) (ni * nj));
                                }
                            }
                        }
                    }
                }
            }
        }
        emit(rows, "cor1", acc, (double) (N * N), (double) (N * N), false);

        // THM30 (skel:53): the ni = 0 case of THM3, kept disaggregated over the
        // phase u of queue i.
        for (int i = 0; i < M; i++) {
            if (i == f) continue;
            for (int ui = 0; ui < K[i]; ui++) {
                // LHS term 1: j <> i, j <> f, BB[m,j] == 0.
                for (int j = 0; j < M; j++) {
                    if (j == i || j == f) continue;
                    for (int nj = 1; nj <= F[j]; nj++) {
                        for (int kj = 0; kj < K[j]; kj++) {
                            for (int hj = 0; hj < K[j]; hj++) {
                                for (int m = 0; m < MR; m++) {
                                    if ((int) BB.get(m, j) != 0) continue;
                                    acc.add(p2Idx.idx(j, nj, kj, i, 0, ui, m), q[j][i].get(kj, hj));
                                }
                            }
                        }
                    }
                }
                // LHS term 2: j == f, MM[m,1] <> i.
                if (f != i) {
                    for (int nj = 1; nj <= F[f]; nj++) {
                        for (int kj = 0; kj < K[f]; kj++) {
                            for (int hj = 0; hj < K[f]; hj++) {
                                for (int m = 0; m < MR; m++) {
                                    if ((int) MM.get(m, 0) - 1 == i) continue;
                                    acc.add(p2Idx.idx(f, nj, kj, i, 0, ui, m), q[f][i].get(kj, hj));
                                }
                            }
                        }
                    }
                }
                // RHS term 1: j <> i, j <> f, BB[m,i] == 0.
                for (int j = 0; j < M; j++) {
                    if (j == i || j == f) continue;
                    for (int nj = 0; nj <= F[j]; nj++) {
                        for (int ki = 0; ki < K[i]; ki++) {
                            for (int hj = 0; hj < K[j]; hj++) {
                                for (int m = 0; m < MR; m++) {
                                    if ((int) BB.get(m, i) != 0) continue;
                                    acc.add(p2Idx.idx(j, nj, hj, i, 1, ki, m), -q[i][j].get(ki, ui));
                                }
                            }
                        }
                    }
                }
                // RHS term 2: j == f, BB[m,i] == 0, nj < F[f].
                if (f != i) {
                    for (int nj = 0; nj < F[f]; nj++) {
                        for (int ki = 0; ki < K[i]; ki++) {
                            for (int hj = 0; hj < K[f]; hj++) {
                                for (int m = 0; m < MR; m++) {
                                    if ((int) BB.get(m, i) != 0) continue;
                                    acc.add(p2Idx.idx(f, nj, hj, i, 1, ki, m), -q[i][f].get(ki, ui));
                                }
                            }
                        }
                    }
                }
                // see _kb/03-api-layer.md for rationale
                for (int m = 0; m < MR; m++) {
                    if ((int) BB.get(m, i) != 1 || (int) MM.get(m, 0) - 1 != i) continue;
                    for (int y = 0; y < K[f]; y++) {
                        double coef = 0.0;
                        for (int w = 0; w < M; w++) {
                            if (w == f || w == i) continue;
                            for (int p = 0; p < K[f]; p++) {
                                coef += q[f][w].get(y, p);
                            }
                        }
                        acc.add(p2Idx.idx(f, F[f], y, i, 1, ui, m), -coef);
                    }
                }
                emit(rows, "thm30_" + i + "_" + ui, acc, 0.0, 0.0, true);
            }
        }

        // THM3 (skel:54): marginal balance for ni in 0..F[i]-1, i <> f.
        for (int i = 0; i < M; i++) {
            if (i == f) continue;
            for (int ni = 0; ni < F[i]; ni++) {
                for (int j = 0; j < M; j++) {
                    if (j == i || j == f) continue;
                    for (int nj = 1; nj <= F[j]; nj++) {
                        for (int kj = 0; kj < K[j]; kj++) {
                            for (int hj = 0; hj < K[j]; hj++) {
                                for (int ui = 0; ui < K[i]; ui++) {
                                    for (int m = 0; m < MR; m++) {
                                        if ((int) BB.get(m, j) != 0) continue;
                                        acc.add(p2Idx.idx(j, nj, kj, i, ni, ui, m), q[j][i].get(kj, hj));
                                    }
                                }
                            }
                        }
                    }
                }
                if (f != i) {
                    for (int nj = 1; nj <= F[f]; nj++) {
                        for (int kj = 0; kj < K[f]; kj++) {
                            for (int hj = 0; hj < K[f]; hj++) {
                                for (int ui = 0; ui < K[i]; ui++) {
                                    for (int m = 0; m < MR; m++) {
                                        if ((int) MM.get(m, 0) - 1 == i) continue;
                                        acc.add(p2Idx.idx(f, nj, kj, i, ni, ui, m), q[f][i].get(kj, hj));
                                    }
                                }
                            }
                        }
                    }
                }
                for (int j = 0; j < M; j++) {
                    if (j == i || j == f) continue;
                    for (int nj = 0; nj <= F[j]; nj++) {
                        for (int ki = 0; ki < K[i]; ki++) {
                            for (int hi = 0; hi < K[i]; hi++) {
                                for (int uj = 0; uj < K[j]; uj++) {
                                    for (int m = 0; m < MR; m++) {
                                        if ((int) BB.get(m, i) != 0) continue;
                                        acc.add(p2Idx.idx(j, nj, uj, i, ni + 1, ki, m), -q[i][j].get(ki, hi));
                                    }
                                }
                            }
                        }
                    }
                }
                if (f != i) {
                    for (int nj = 0; nj < F[f]; nj++) {
                        for (int ki = 0; ki < K[i]; ki++) {
                            for (int hi = 0; hi < K[i]; hi++) {
                                for (int uj = 0; uj < K[f]; uj++) {
                                    for (int m = 0; m < MR; m++) {
                                        if ((int) BB.get(m, i) != 0) continue;
                                        acc.add(p2Idx.idx(f, nj, uj, i, ni + 1, ki, m), -q[i][f].get(ki, hi));
                                    }
                                }
                            }
                        }
                    }
                }
                for (int m = 0; m < MR; m++) {
                    if ((int) BB.get(m, i) != 1 || (int) MM.get(m, 0) - 1 != i) continue;
                    for (int ki = 0; ki < K[i]; ki++) {
                        for (int uj = 0; uj < K[f]; uj++) {
                            double coef = 0.0;
                            for (int w = 0; w < M; w++) {
                                if (w == f || w == i) continue;
                                for (int p = 0; p < K[f]; p++) {
                                    coef += q[f][w].get(uj, p);
                                }
                            }
                            acc.add(p2Idx.idx(f, F[f], uj, i, ni + 1, ki, m), -coef);
                        }
                    }
                }
                emit(rows, "thm3_" + i + "_" + ni, acc, 0.0, 0.0, true);
            }
        }

        // see _kb/03-api-layer.md for rationale
        for (int ni = 0; ni < F[f]; ni++) {
            for (int j = 0; j < M; j++) {
                if (j == f) continue;
                for (int nj = 1; nj <= F[j]; nj++) {
                    for (int kj = 0; kj < K[j]; kj++) {
                        for (int hj = 0; hj < K[j]; hj++) {
                            for (int uf = 0; uf < K[f]; uf++) {
                                for (int m = 0; m < MR; m++) {
                                    if ((int) BB.get(m, j) != 0) continue;
                                    acc.add(p2Idx.idx(j, nj, kj, f, ni, uf, m), q[j][f].get(kj, hj));
                                }
                            }
                        }
                    }
                }
            }
            for (int j = 0; j < M; j++) {
                if (j == f) continue;
                for (int nj = 0; nj <= F[j]; nj++) {
                    for (int kf = 0; kf < K[f]; kf++) {
                        for (int hf = 0; hf < K[f]; hf++) {
                            for (int uj = 0; uj < K[j]; uj++) {
                                acc.add(p2Idx.idx(j, nj, uj, f, ni + 1, kf, 0), -q[f][j].get(kf, hf));
                            }
                        }
                    }
                }
            }
            emit(rows, "thm3f_" + ni, acc, 0.0, 0.0, true);
        }

        // THM3I (skel:56): couples blocking configurations of successive depth.
        for (int z = 0; z < ZM; z++) {
            for (int j = 0; j < M; j++) {
                if (j == f) continue;
                for (int nj = 1; nj <= F[j]; nj++) {
                    for (int kj = 0; kj < K[j]; kj++) {
                        for (int hj = 0; hj < K[j]; hj++) {
                            for (int uf = 0; uf < K[f]; uf++) {
                                for (int m = 0; m < MR; m++) {
                                    if ((int) BB.get(m, j) != 0 || ZZ[m] != z) continue;
                                    acc.add(p2Idx.idx(j, nj, kj, f, F[f], uf, m), q[j][f].get(kj, hj));
                                }
                            }
                        }
                    }
                }
            }
            for (int j = 0; j < M; j++) {
                if (j == f) continue;
                for (int nj = 0; nj <= F[j]; nj++) {
                    for (int kf = 0; kf < K[f]; kf++) {
                        for (int hf = 0; hf < K[f]; hf++) {
                            for (int uj = 0; uj < K[j]; uj++) {
                                for (int m = 0; m < MR; m++) {
                                    if (ZZ[m] != z + 1) continue;
                                    acc.add(p2Idx.idx(j, nj, uj, f, F[f], kf, m), -q[f][j].get(kf, hf));
                                }
                            }
                        }
                    }
                }
            }
            emit(rows, "thm3i_" + z, acc, 0.0, 0.0, true);
        }

        // see _kb/03-api-layer.md for rationale
        for (int m = 0; m < MR; m++) {
            if (ZZ[m] != ZM - 1) continue;
            for (int j = 0; j < M; j++) {
                if (j == f || (int) BB.get(m, j) != 0 || (int) MM1.get(m, j) <= 0) continue;
                for (int nj = 1; nj <= F[j]; nj++) {
                    for (int kj = 0; kj < K[j]; kj++) {
                        for (int hj = 0; hj < K[j]; hj++) {
                            for (int uf = 0; uf < K[f]; uf++) {
                                acc.add(p2Idx.idx(j, nj, kj, f, F[f], uf, m), q[j][f].get(kj, hj));
                            }
                        }
                    }
                }
            }
            for (int j = 0; j < M; j++) {
                if (j == f || (int) BB.get(m, j) != 0 || (int) MM1.get(m, j) <= 0) continue;
                int mp = (int) MM1.get(m, j) - 1;
                if (mp < 0 || mp >= MR) continue;
                for (int kf = 0; kf < K[f]; kf++) {
                    for (int uf = 0; uf < K[f]; uf++) {
                        for (int w = 0; w < M; w++) {
                            if (w == f) continue;
                            acc.add(p2Idx.idx(f, F[f], kf, f, F[f], kf, mp), -q[f][w].get(kf, uf));
                        }
                    }
                }
            }
            emit(rows, "thm3l_" + m, acc, 0.0, 0.0, true);
        }

        // THM4 (skel:59): mean-population dominance inequality.
        for (int j = 0; j < M; j++) {
            for (int kj = 0; kj < K[j]; kj++) {
                for (int i = 0; i < M; i++) {
                    for (int m = 0; m < MR; m++) {
                        for (int t = 0; t < M; t++) {
                            for (int ht = 0; ht < K[t]; ht++) {
                                for (int njt = 0; njt <= F[j]; njt++) {
                                    for (int nt = 1; nt <= F[t]; nt++) {
                                        acc.add(p2Idx.idx(j, njt, kj, t, nt, ht, m), (double) nt);
                                    }
                                }
                            }
                        }
                        for (int hi = 0; hi < K[i]; hi++) {
                            for (int njt = 0; njt <= F[j]; njt++) {
                                for (int ni = 1; ni <= F[i]; ni++) {
                                    acc.add(p2Idx.idx(j, njt, kj, i, ni, hi, m), -(double) N);
                                }
                            }
                        }
                        emit(rows, "thm4_" + j + "_" + kj + "_" + i + "_" + m, acc,
                                0.0, Double.POSITIVE_INFINITY, false);
                    }
                }
            }
        }

        // Objective U1min / U1max (example_bas.mod:62).
        double[] objectiveCoeffs = new double[nCols];
        int targetQueue = objectiveQueue - 1;
        for (int m = 0; m < MR; m++) {
            for (int ki = 0; ki < K[targetQueue]; ki++) {
                for (int ni = 1; ni <= F[targetQueue]; ni++) {
                    int idx = p2Idx.idx(targetQueue, ni, ki, targetQueue, ni, ki, m);
                    if (idx >= 0 && colMap[idx] >= 0) objectiveCoeffs[colMap[idx]] = 1.0;
                }
            }
        }

        // see _kb/03-api-layer.md for rationale
        String dumpPath = System.getProperty("jline.mapqn.dumpmps");
        if (dumpPath != null) {
            try {
                java.io.PrintWriter pw = new java.io.PrintWriter(dumpPath);
                pw.println("NVARS " + nCols);
                for (int c = 0; c < nCols; c++) {
                    pw.println("BND " + c + " 0.0 1.0");
                    if (objectiveCoeffs[c] != 0.0) {
                        pw.println("OBJ " + c + " " + objectiveCoeffs[c]);
                    }
                }
                for (int r2 = 0; r2 < rows.size(); r2++) {
                    Row row = rows.get(r2);
                    pw.println("ROW " + row.name
                            + " " + (row.hasLower ? Double.toString(row.lower) : "NONE")
                            + " " + (row.hasUpper ? Double.toString(row.upper) : "NONE"));
                    for (int t = 0; t < row.idx.length; t++) {
                        pw.println("COEF " + row.idx[t] + " " + row.val[t]);
                    }
                }
                pw.close();
                System.err.println("[mapqn_qr_bounds_bas] wrote LP dump to " + dumpPath
                        + " (" + nCols + " cols, " + rows.size() + " rows)");
            } catch (Exception dumpEx) {
                System.err.println("[mapqn_qr_bounds_bas] LP dump failed: " + dumpEx);
            }
        }

        // see _kb/03-api-layer.md for rationale
        int nRows = 0;
        for (int r2 = 0; r2 < rows.size(); r2++) {
            if (rows.get(r2).idx.length > 0) nRows++;
        }
        int mRows = nRows + nCols;

        int[] colCount = new int[nCols];
        for (int r2 = 0; r2 < rows.size(); r2++) {
            Row row = rows.get(r2);
            if (row.idx.length == 0) continue;
            for (int t = 0; t < row.idx.length; t++) colCount[row.idx[t]]++;
        }
        int nnz = nCols;
        for (int c = 0; c < nCols; c++) nnz += colCount[c];

        int[] Ap = new int[nCols + 1];
        int pos = 0;
        for (int c = 0; c < nCols; c++) {
            Ap[c] = pos;
            pos += colCount[c] + 1;
        }
        Ap[nCols] = pos;

        int[] Ai = new int[nnz];
        double[] Ax = new double[nnz];
        int[] fill = new int[nCols];
        double[] rowLower = new double[mRows];
        double[] rowUpper = new double[mRows];
        int ri = 0;
        for (int r2 = 0; r2 < rows.size(); r2++) {
            Row row = rows.get(r2);
            if (row.idx.length == 0) continue;
            rowLower[ri] = row.hasLower ? row.lower : -OSQP.OSQP_INFTY;
            rowUpper[ri] = row.hasUpper ? row.upper : OSQP.OSQP_INFTY;
            for (int t = 0; t < row.idx.length; t++) {
                int c = row.idx[t];
                int w = Ap[c] + fill[c]++;
                Ai[w] = ri;
                Ax[w] = row.val[t];
            }
            ri++;
        }
        for (int c = 0; c < nCols; c++) {
            int w = Ap[c] + fill[c]++;
            Ai[w] = ri;
            Ax[w] = 1.0;
            rowLower[ri] = 0.0;
            rowUpper[ri] = 1.0;
            ri++;
        }

        CSCMatrix A = new CSCMatrix(mRows, nCols, nnz, Ap, Ai, Ax);
        CSCMatrix P = new CSCMatrix(nCols, nCols, 0, new int[nCols + 1], new int[0], new double[0]);

        boolean maximise = "max".equals(sense);
        double[] linearCost = new double[nCols];
        for (int c = 0; c < nCols; c++) {
            linearCost[c] = maximise ? -objectiveCoeffs[c] : objectiveCoeffs[c];
        }

        try {
            OSQP.Settings settings = new OSQP.Settings();
            settings.max_iter = OSQP_MAX_ITER;
            settings.eps_abs = OSQP_EPS;
            settings.eps_rel = OSQP_EPS;
            settings.rho = OSQP_RHO;
            settings.scaling = OSQP_SCALING;
            settings.polish = true;
            settings.verbose = false;

            OSQP solver = new OSQP(
                    new OSQP.Data(nCols, mRows, P, A, linearCost, rowLower, rowUpper, 0.0),
                    settings);
            OSQP.Status status = solver.solve();

            if (System.getProperty("jline.mapqn.debug") != null) {
                System.err.println("[mapqn_qr_bounds_bas] cols=" + nCols + " rows=" + mRows
                        + " nnz=" + nnz + " status=" + status);
            }

            if (status != OSQP.Status.SOLVED && status != OSQP.Status.SOLVED_INACCURATE) {
                return new Mapqn_solution(Double.NaN, new HashMap<String, Double>());
            }

            double[] point = extractPrimal(solver, nCols);
            double objective = 0.0;
            for (int c = 0; c < nCols; c++) {
                objective += objectiveCoeffs[c] * point[c];
            }

            Map<String, Double> variables = new HashMap<String, Double>();
            for (int i = 0; i < M; i++) {
                double totalU = 0.0;
                for (int m = 0; m < MR; m++) {
                    for (int ki = 0; ki < K[i]; ki++) {
                        for (int ni = 1; ni <= F[i]; ni++) {
                            totalU += value(point, colMap, p2Idx.idx(i, ni, ki, i, ni, ki, m));
                        }
                    }
                }
                variables.put("U_" + (i + 1), Double.valueOf(totalU));
            }
            for (int i = 0; i < M; i++) {
                double totalE = 0.0;
                for (int ki = 0; ki < K[i]; ki++) {
                    double ev = value(point, colMap, eIdx.idx(i, ki));
                    totalE += ev;
                    variables.put("e_" + (i + 1) + "_" + (ki + 1), Double.valueOf(ev));
                }
                variables.put("Ueff_" + (i + 1), Double.valueOf(totalE));
                variables.put("pb_" + (i + 1),
                        Double.valueOf(variables.get("U_" + (i + 1)).doubleValue() - totalE));
            }

            return new Mapqn_solution(objective, variables);
        } catch (Exception e) {
            if (System.getProperty("jline.mapqn.debug") != null) {
                System.err.println("[mapqn_qr_bounds_bas] exception: " + e);
            }
            return new Mapqn_solution(Double.NaN, new HashMap<String, Double>());
        }
    }

    private static double value(double[] point, int[] colMap, int idx) {
        if (idx < 0) return 0.0;
        int c = colMap[idx];
        return (c < 0) ? 0.0 : point[c];
    }

    /**
     * josqp settings. The three constants are tuned together; changing one
     * without re-measuring the others is not meaningful.
     *
     * <p>Measured on tomacs_qrf/example_bas.mod (M=5, N=10, MR=5), against the
     * AMPL optimum U1min = 0.4776625414, U1max = 0.8077316296, which glpsol also
     * reproduces from the LP this method dumps under
     * {@code -Djline.mapqn.dumpmps} (6707 and 9111 simplex iterations, ~3 s).
     * The column is the absolute deviation of U1min and the wall time:
     *
     * <pre>
     *   eps    rho=0.1 scaling=10      rho=1.0 scaling=25
     *   1e-4   6.1e-03    5 s          -
     *   1e-5   9.7e-04   13 s          -
     *   1e-6   7.5e-04   26 s          7.9e-05   33 s
     *   1e-7   9.8e-06  292 s          6.8e-06   20 s
     * </pre>
     *
     * <p>The deviation is NOT monotone in eps, which is characteristic of ADMM:
     * the iterate approaches the optimal face quickly and then crawls along it,
     * so a looser tolerance can happen to stop nearer the vertex. rho = 1.0 with
     * scaling = 25 is both the most accurate and, at eps = 1e-7, an order of
     * magnitude faster than the josqp default rho = 0.1, which spends 292 s to
     * reach the same accuracy. At the settings below the max direction lands
     * within 2.1e-07, and the two small instances asserted by MapqnAPITest
     * within 1e-16 (M=2, N=2) and 4.6e-07 (M=2, N=3).
     *
     * <p>ACCURACY LIMIT, and it is a real one. ADMM converges linearly and its
     * polish step does not recover the exact vertex of this degenerate LP, so
     * about five decimals is the ceiling. MATLAB (matlab/lib/qrf/qrf_bas.m) and
     * native Python (api/mapqn/qr_bounds_bas.py) agree with AMPL to ~3e-11, so a
     * JAR bound is NOT interchangeable with theirs beyond that.
     */
    private static final double OSQP_EPS = 1.0e-7;

    private static final double OSQP_RHO = 1.0;

    private static final int OSQP_SCALING = 25;

    private static final int OSQP_MAX_ITER = 200000;

    /**
     * Read the primal solution out of a solved josqp workspace.
     *
     * <p>josqp 0.6.5 exposes {@code Workspace.solution} and {@code Solution.x}
     * as package-private, with no accessor, so reflection is the only route to
     * the primal point from outside {@code com.quantego.josqp}.
     */
    private static double[] extractPrimal(OSQP solver, int nCols) throws Exception {
        Object workspace = solver.getWorkspace();
        java.lang.reflect.Field solutionField = workspace.getClass().getDeclaredField("solution");
        solutionField.setAccessible(true);
        Object solution = solutionField.get(workspace);
        java.lang.reflect.Field xField = solution.getClass().getDeclaredField("x");
        xField.setAccessible(true);
        double[] x = (double[]) xField.get(solution);
        double[] copy = new double[nCols];
        System.arraycopy(x, 0, copy, 0, nCols);
        return copy;
    }

    /**
     * Relative magnitude below which an accumulated coefficient is treated as a
     * structural zero rather than a genuine term. See {@link #emit}.
     */
    private static final double DUST_REL = 1.0e-12;

    /**
     * Accumulator for one constraint row over the compacted column space.
     *
     * <p>Each block below adds its left-hand-side terms and subtracts its
     * right-hand-side terms into this accumulator, which is then emitted as a
     * sparse row and reset. Only the columns actually touched are cleared, so
     * the cost is proportional to the row's number of nonzeros rather than to
     * the width of the program.
     */
    private static final class Acc {
        final int[] colMap;
        final double[] coeffs;
        final int[] touched;
        final boolean[] marked;
        int nTouched;

        Acc(int[] colMap, int nCols) {
            this.colMap = colMap;
            this.coeffs = new double[nCols];
            this.touched = new int[nCols];
            this.marked = new boolean[nCols];
        }

        /**
         * Accumulate {@code value} on variable {@code idx}. A variable that does
         * not exist, or that is pinned to zero by one of the ZERO families and
         * hence carries no column, contributes nothing and is dropped here.
         */
        void add(int idx, double value) {
            if (idx < 0) return;
            int c = colMap[idx];
            if (c < 0) return;
            if (!marked[c]) {
                marked[c] = true;
                touched[nTouched++] = c;
            }
            coeffs[c] += value;
        }

        void reset() {
            for (int t = 0; t < nTouched; t++) {
                coeffs[touched[t]] = 0.0;
                marked[touched[t]] = false;
            }
            nTouched = 0;
        }
    }

    /** One assembled constraint row in sparse form, plus its bounds. */
    private static final class Row {
        final String name;
        final int[] idx;
        final double[] val;
        final double lower;
        final double upper;
        final boolean hasLower;
        final boolean hasUpper;

        Row(String name, int[] idx, double[] val, double lower, double upper,
            boolean hasLower, boolean hasUpper) {
            this.name = name;
            this.idx = idx;
            this.val = val;
            this.lower = lower;
            this.upper = upper;
            this.hasLower = hasLower;
            this.hasUpper = hasUpper;
        }
    }

    /**
     * Emit the accumulated row and reset the accumulator.
     *
     * <p>When {@code equilibrate} is set the row is scaled so that its
     * largest-magnitude coefficient becomes 1 and cancellation dust below
     * {@link #DUST_REL} of that magnitude is discarded. Many terms of the
     * balance equations appear on both sides and must cancel exactly; because
     * the blocks accumulate in different loop orders, floating-point addition
     * leaves a residue of a few ulp instead of a clean zero, and emitting that
     * residue turns a structural zero of the constraint matrix into a nonzero
     * and inflates the coefficient ratio of the program by orders of magnitude.
     *
     * <p>Scaling a constraint by a positive constant leaves its solution set
     * unchanged only if the right-hand side is scaled identically, so
     * {@code equilibrate} may be set ONLY on a row whose bounds are zero (or
     * infinite). Every caller that sets it emits a row at level 0.
     */
    private static void emit(List<Row> rows, String name, Acc acc,
                             double lower, double upper, boolean equilibrate) {
        boolean hasLower = lower != Double.NEGATIVE_INFINITY;
        boolean hasUpper = upper != Double.POSITIVE_INFINITY;

        double maxAbs = 0.0;
        if (equilibrate) {
            for (int t = 0; t < acc.nTouched; t++) {
                double a = Math.abs(acc.coeffs[acc.touched[t]]);
                if (a > maxAbs) maxAbs = a;
            }
        }
        double scale = (equilibrate && maxAbs > 0.0) ? (1.0 / maxAbs) : 1.0;
        double dust = equilibrate ? maxAbs * DUST_REL : 0.0;

        int n = 0;
        for (int t = 0; t < acc.nTouched; t++) {
            double a = Math.abs(acc.coeffs[acc.touched[t]]);
            if (a > dust && a != 0.0) n++;
        }
        int[] idx = new int[n];
        double[] val = new double[n];
        int w = 0;
        for (int t = 0; t < acc.nTouched; t++) {
            int c = acc.touched[t];
            double a = Math.abs(acc.coeffs[c]);
            if (a > dust && a != 0.0) {
                idx[w] = c;
                val[w] = acc.coeffs[c] * scale;
                w++;
            }
        }
        acc.reset();
        rows.add(new Row(name, idx, val, lower, upper, hasLower, hasUpper));
    }

    private interface IndexFn7 {
        int idx(int j, int nj, int kj, int i, int ni, int hi, int m);
    }

    private interface IndexFn2 {
        int idx(int i, int ki);
    }
}
