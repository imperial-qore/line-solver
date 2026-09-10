/**
 * @file Discrete-time Markov chain steady-state solver
 *
 * Computes the steady-state probability distribution for DTMCs by converting the
 * transition matrix problem (P-I)x = 0 into a CTMC-equivalent system and leveraging
 * the robust CTMC solver with automatic reducibility handling.
 *
 * @since LINE 3.0
 */
package jline.api.mc;

import java.util.Arrays;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.Map;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixEntry;

public final class Dtmc_solve {
    private Dtmc_solve() {}

    /**
     * Bounded memo over the exact content of P.
     *
     * dtmc_solve is a pure function of P, and the layered fixed point asks it the same
     * question over and over: a layer's routing topology does not change between LN
     * iterations -- only its rates do, and visits do not depend on rates -- so
     * snRefreshVisits re-solves one identical chain per layer per iteration.
     *
     * The key is the matrix's exact content, and a hit re-compares the whole key rather
     * than trusting a digest, so no collision can return another matrix's answer. Only
     * stored entries enter the key, which makes an explicit zero and an absent entry
     * different keys: that is a miss, never a wrong hit.
     */
    private static final int CACHE_MAX = 32;

    private static final Map<Key, Matrix> CACHE =
            new LinkedHashMap<Key, Matrix>(16, 0.75f, true) {
                private static final long serialVersionUID = 1L;
                @Override
                protected boolean removeEldestEntry(Map.Entry<Key, Matrix> eldest) {
                    return size() > CACHE_MAX;
                }
            };

    private static final class Key {
        private final int rows;
        private final int cols;
        private final long[] data;   // row, col, raw bits of the value, per stored entry
        private final int hash;

        Key(Matrix m) {
            this.rows = m.getNumRows();
            this.cols = m.getNumCols();
            int nz = m.getNonZeros();
            long[] d = new long[3 * nz];
            int k = 0;
            Iterator<MatrixEntry> it = m.nonZeroIterator();
            while (it.hasNext() && k + 3 <= d.length) {
                MatrixEntry e = it.next();
                d[k++] = e.row;
                d[k++] = e.col;
                d[k++] = Double.doubleToRawLongBits(e.value);
            }
            // The iterator may yield fewer entries than getNonZeros() reports; trim so a
            // trailing run of zeros cannot make two different matrices compare equal.
            this.data = (k == d.length) ? d : Arrays.copyOf(d, k);
            this.hash = 31 * (31 * rows + cols) + Arrays.hashCode(this.data);
        }

        @Override
        public boolean equals(Object o) {
            if (!(o instanceof Key)) return false;
            Key k = (Key) o;
            return rows == k.rows && cols == k.cols && Arrays.equals(data, k.data);
        }

        @Override
        public int hashCode() {
            return hash;
        }
    }

    /**
     * Returns the steady-state solution of a DTMC.
     *
     * @param P Transition matrix of the DTMC
     * @return Steady-state solution vector of the DTMC
     */
    public static Matrix dtmc_solve(Matrix P) {
        Key key = new Key(P);
        synchronized (CACHE) {
            Matrix hit = CACHE.get(key);
            // Copied out, so a caller that writes into the result cannot corrupt the entry.
            if (hit != null) return hit.copy();
        }

        Matrix Plocal = P.copy();
        // P - eye(size(P))
        for (int i = 0; i < Plocal.getNumRows(); i++) {
            Plocal.set(i, i, Plocal.get(i, i) - 1.0);
        }
        Matrix pi = Ctmc_solve.ctmc_solve(Plocal);

        synchronized (CACHE) {
            CACHE.put(key, pi.copy());
        }
        return pi;
    }
}
