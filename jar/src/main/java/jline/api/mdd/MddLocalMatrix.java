package jline.api.mdd;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * A local rate matrix W_k^e of the Kronecker descriptor, held row-compressed.
 *
 * <p>The aggregation only ever walks a row of W and needs its row sums (the
 * local enabling rates lambda_k^e), so the rows are stored as parallel
 * column/value arrays rather than as a general sparse matrix. Duplicate
 * triplets are summed when the matrix is built, so a caller may emit the same
 * (i,j) more than once.</p>
 */
public class MddLocalMatrix {

    /** Order of the (square) local matrix, i.e. the level domain. */
    public final int dim;
    /** cols[i] holds the column indices of the nonzeros of row i. */
    public final int[][] cols;
    /** vals[i] holds the values of the nonzeros of row i, aligned with cols[i]. */
    public final double[][] vals;
    /** rowSum[i] is the local enabling rate lambda[i]. */
    public final double[] rowSum;
    /** Total number of stored nonzeros. */
    public final int nnz;

    private MddLocalMatrix(int dim, int[][] cols, double[][] vals, double[] rowSum, int nnz) {
        this.dim = dim;
        this.cols = cols;
        this.vals = vals;
        this.rowSum = rowSum;
        this.nnz = nnz;
    }

    /** The identity of the given order, used for a level an event does not touch. */
    public static MddLocalMatrix identity(int dim) {
        int[][] cols = new int[dim][];
        double[][] vals = new double[dim][];
        double[] rowSum = new double[dim];
        for (int i = 0; i < dim; i++) {
            cols[i] = new int[]{i};
            vals[i] = new double[]{1.0};
            rowSum[i] = 1.0;
        }
        return new MddLocalMatrix(dim, cols, vals, rowSum, dim);
    }

    /** Incremental triplet builder; duplicate entries are accumulated. */
    public static class Builder {
        private final int dim;
        private final List<Map<Integer, Double>> rows;

        public Builder(int dim) {
            this.dim = dim;
            this.rows = new ArrayList<Map<Integer, Double>>(dim);
            for (int i = 0; i < dim; i++) {
                this.rows.add(new HashMap<Integer, Double>());
            }
        }

        /** Accumulate value into entry (i,j). A zero value is dropped. */
        public Builder add(int i, int j, double value) {
            if (value == 0.0) {
                return this;
            }
            Map<Integer, Double> row = rows.get(i);
            Integer key = Integer.valueOf(j);
            Double old = row.get(key);
            row.put(key, Double.valueOf(old == null ? value : old.doubleValue() + value));
            return this;
        }

        public MddLocalMatrix build() {
            int[][] cols = new int[dim][];
            double[][] vals = new double[dim][];
            double[] rowSum = new double[dim];
            int nnz = 0;
            for (int i = 0; i < dim; i++) {
                Map<Integer, Double> row = rows.get(i);
                int n = row.size();
                int[] c = new int[n];
                double[] v = new double[n];
                int p = 0;
                double s = 0.0;
                // ascending column order, so two builds of the same matrix agree
                int[] keys = new int[n];
                int q = 0;
                for (Integer k : row.keySet()) {
                    keys[q++] = k.intValue();
                }
                java.util.Arrays.sort(keys);
                for (int t = 0; t < n; t++) {
                    double value = row.get(Integer.valueOf(keys[t])).doubleValue();
                    c[p] = keys[t];
                    v[p] = value;
                    s += value;
                    p++;
                }
                cols[i] = c;
                vals[i] = v;
                rowSum[i] = s;
                nnz += n;
            }
            return new MddLocalMatrix(dim, cols, vals, rowSum, nnz);
        }
    }
}
