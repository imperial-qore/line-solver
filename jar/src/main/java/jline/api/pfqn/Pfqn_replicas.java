/**
 * Station replica consolidation utilities for Product-Form Queueing Networks
 *
 * @since LINE 3.0
 */
package jline.api.pfqn;

import java.util.ArrayList;
import java.util.List;

import jline.util.Triple;

import jline.GlobalConstants;
import jline.util.matrix.Matrix;

public final class Pfqn_replicas {
    private Pfqn_replicas() {}

    /**
     * Consolidate replicated stations into unique stations with multiplicity.
     */
    public static PfqnUniqueResult pfqn_unique(Matrix L, Matrix mu, Matrix gamma) {
        int M = L.getNumRows();
        int R = L.getNumCols();
        double tol = GlobalConstants.FineTol;

        int fingerprintCols = R + ((mu != null) ? mu.getNumCols() : 0) + ((gamma != null) ? gamma.getNumCols() : 0);
        Matrix fingerprint = new Matrix(M, fingerprintCols);

        for (int i = 0; i < M; i++) {
            int col = 0;
            for (int j = 0; j < R; j++) {
                fingerprint.set(i, col++, L.get(i, j));
            }
            if (mu != null) {
                for (int j = 0; j < mu.getNumCols(); j++) {
                    fingerprint.set(i, col++, mu.get(i, j));
                }
            }
            if (gamma != null) {
                for (int j = 0; j < gamma.getNumCols(); j++) {
                    fingerprint.set(i, col++, gamma.get(i, j));
                }
            }
        }

        int[] mapping = new int[M];
        for (int i = 0; i < M; i++) mapping[i] = -1;
        List<Integer> uniqueIdx = new ArrayList<Integer>();
        List<Integer> miList = new ArrayList<Integer>();

        for (int i = 0; i < M; i++) {
            if (mapping[i] == -1) {
                uniqueIdx.add(i);
                int groupIdx = uniqueIdx.size() - 1;
                mapping[i] = groupIdx;
                int count = 1;

                for (int j = i + 1; j < M; j++) {
                    if (mapping[j] == -1) {
                        double maxDiff = 0.0;
                        for (int k = 0; k < fingerprintCols; k++) {
                            double diff = Math.abs(fingerprint.get(i, k) - fingerprint.get(j, k));
                            if (diff > maxDiff) maxDiff = diff;
                        }
                        if (maxDiff < tol) {
                            mapping[j] = groupIdx;
                            count++;
                        }
                    }
                }
                miList.add(count);
            }
        }

        int M_unique = uniqueIdx.size();
        Matrix L_unique = new Matrix(M_unique, R);
        for (int i = 0; i < M_unique; i++) {
            for (int j = 0; j < R; j++) {
                L_unique.set(i, j, L.get(uniqueIdx.get(i), j));
            }
        }

        Matrix mu_unique = null;
        if (mu != null) {
            mu_unique = new Matrix(M_unique, mu.getNumCols());
            for (int i = 0; i < M_unique; i++) {
                for (int j = 0; j < mu.getNumCols(); j++) {
                    mu_unique.set(i, j, mu.get(uniqueIdx.get(i), j));
                }
            }
        }

        Matrix gamma_unique = null;
        if (gamma != null) {
            gamma_unique = new Matrix(M_unique, gamma.getNumCols());
            for (int i = 0; i < M_unique; i++) {
                for (int j = 0; j < gamma.getNumCols(); j++) {
                    gamma_unique.set(i, j, gamma.get(uniqueIdx.get(i), j));
                }
            }
        }

        Matrix mi = new Matrix(1, M_unique);
        for (int i = 0; i < M_unique; i++) {
            mi.set(0, i, miList.get(i).doubleValue());
        }

        return new PfqnUniqueResult(L_unique, mu_unique, gamma_unique, mi, mapping);
    }

    public static PfqnUniqueResult pfqn_unique(Matrix L, Matrix mu) {
        return pfqn_unique(L, mu, null);
    }

    public static PfqnUniqueResult pfqn_unique(Matrix L) {
        return pfqn_unique(L, null, null);
    }

    /**
     * Expand per-station metrics from reduced model to original dimensions.
     */
    public static Triple<Matrix, Matrix, Matrix> pfqn_expand(Matrix QN, Matrix UN, Matrix CN, int[] mapping) {
        int M_original = mapping.length;
        int R = QN.getNumCols();

        Matrix QN_full = new Matrix(M_original, R);
        Matrix UN_full = new Matrix(M_original, R);
        Matrix CN_full = new Matrix(M_original, R);

        for (int i = 0; i < M_original; i++) {
            int uniqueIdx = mapping[i];
            for (int r = 0; r < R; r++) {
                QN_full.set(i, r, QN.get(uniqueIdx, r));
                UN_full.set(i, r, UN.get(uniqueIdx, r));
                CN_full.set(i, r, CN.get(uniqueIdx, r));
            }
        }

        return new Triple<Matrix, Matrix, Matrix>(QN_full, UN_full, CN_full);
    }

    /**
     * Combine user-provided multiplicity vector with detected replica multiplicity.
     */
    public static Matrix pfqn_combine_mi(Matrix mi, int[] mapping, int M_unique) {
        Matrix mi_combined = new Matrix(1, M_unique);
        for (int i = 0; i < mapping.length; i++) {
            mi_combined.set(0, mapping[i], mi_combined.get(0, mapping[i]) + mi.get(0, i));
        }
        return mi_combined;
    }
}
