/**
 * @file Order-independent (OI) load-dependent functional server.
 *
 * OI generalization of the load-dependent functional server (FNC) of Casale,
 * "On Single-Class Load-Dependent Normalizing Constant Equations", QEST 2006
 * (Theorem 3, Corollary 1). Given the balance function Phi of an existing OI
 * station and a queue-dependent target f(n), it builds an auxiliary OI station
 * whose balance function Psi satisfies (Psi * Phi)(n) = (1 + f(n)) Phi(n),
 * then inverts Psi to the FNC rate mu_f(n).
 *
 * Port of matlab/src/api/pfqn/pfqn_oi_fnc.m.
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.ld;

import jline.GlobalConstants;
import jline.io.Ret;
import jline.util.matrix.Matrix;

import java.util.function.ToDoubleFunction;

public final class Pfqn_oi_fnc {
    private Pfqn_oi_fnc() {}

    /** Default target f(n) = sum(n) (total occupancy). */
    public static Ret.pfqnOifnc pfqn_oi_fnc(Matrix Phi, int[] N) {
        return pfqn_oi_fnc(Phi, N, null);
    }

    /**
     * OI functional-server balance/rate for a target queue-dependent function.
     *
     * @param Phi balance function of the existing OI station over the lattice,
     *            given as a flat column-major vector of length prod(N+1).
     * @param N   (1 x R) closed population vector, finite.
     * @param f   target function f(n) with f(0)=0; null defaults to sum(n).
     * @return the FNC rate handle muf, and Psi/mu tabulated as flat
     *         column-major vectors over the lattice.
     */
    public static Ret.pfqnOifnc pfqn_oi_fnc(Matrix Phi, int[] N, ToDoubleFunction<int[]> f) {
        final ToDoubleFunction<int[]> fTarget;
        if (f == null) {
            fTarget = new ToDoubleFunction<int[]>() {
                @Override
                public double applyAsDouble(int[] n) {
                    double s = 0.0;
                    for (int v : n) {
                        s += v;
                    }
                    return s;
                }
            };
        } else {
            fTarget = f;
        }

        int R = N.length;
        final int[] shp = new int[R];
        int total = 1;
        for (int r = 0; r < R; r++) {
            shp[r] = N[r] + 1;
            total *= shp[r];
        }

        // Flatten Phi in column-major order.
        double[] Phiv = new double[total];
        int rows = Phi.getNumRows();
        int cols = Phi.getNumCols();
        if (rows * cols != total) {
            throw new RuntimeException("numel(Phi) must equal prod(N+1).");
        }
        int idxp = 0;
        for (int j = 0; j < cols; j++) {
            for (int i = 0; i < rows; i++) {
                Phiv[idxp++] = Phi.get(i, j);
            }
        }

        // Column-major strides and decoded subscripts.
        final int[] stride = new int[R];
        stride[0] = 1;
        for (int d = 1; d < R; d++) {
            stride[d] = stride[d - 1] * shp[d - 1];
        }
        int[][] subs = new int[total][R];
        for (int i = 0; i < total; i++) {
            int li = i;
            for (int d = 0; d < R; d++) {
                subs[i][d] = li % shp[d];
                li /= shp[d];
            }
        }

        // Step 1: deconvolve (Psi * Phi)(n) = (1 + f(n)) Phi(n) for Psi.
        double[] Psiv = new double[total];
        for (int i = 0; i < total; i++) {
            int[] n = subs[i];
            double acc = (1.0 + fTarget.applyAsDouble(n)) * Phiv[i];
            for (int j = 0; j < i; j++) {
                int[] k = subs[j];
                boolean sub = true;
                for (int d = 0; d < R; d++) {
                    if (k[d] > n[d]) {
                        sub = false;
                        break;
                    }
                }
                if (sub) {
                    int idx = 0;
                    for (int d = 0; d < R; d++) {
                        idx += (n[d] - k[d]) * stride[d];
                    }
                    acc -= Psiv[j] * Phiv[idx];
                }
            }
            Psiv[i] = acc;
        }

        // Step 2: balanced-fairness inversion of Psi to the FNC rate.
        final double[] muv = new double[total];
        for (int i = 0; i < total; i++) {
            int[] n = subs[i];
            boolean empty = true;
            for (int d = 0; d < R; d++) {
                if (n[d] != 0) {
                    empty = false;
                    break;
                }
            }
            if (empty) {
                muv[i] = 0.0;            // empty state
                continue;
            }
            double denom = Psiv[i];
            if (denom == 0) {
                muv[i] = GlobalConstants.Inf;   // non-physical / undefined rate
                continue;
            }
            double num = 0.0;
            for (int r = 0; r < R; r++) {
                if (n[r] > 0) {
                    num += Psiv[i - stride[r]];
                }
            }
            muv[i] = num / denom;
        }

        Matrix PsiM = new Matrix(total, 1);
        Matrix muM = new Matrix(total, 1);
        for (int i = 0; i < total; i++) {
            PsiM.set(i, 0, Psiv[i]);
            muM.set(i, 0, muv[i]);
        }

        ToDoubleFunction<int[]> muf = new ToDoubleFunction<int[]>() {
            @Override
            public double applyAsDouble(int[] n) {
                if (n.length != shp.length) {
                    return GlobalConstants.Inf;
                }
                int idx = 0;
                for (int d = 0; d < n.length; d++) {
                    if (n[d] < 0 || n[d] > shp[d] - 1) {
                        return GlobalConstants.Inf; // outside the tabulated lattice
                    }
                    idx += n[d] * stride[d];
                }
                return muv[idx];
            }
        };

        return new Ret.pfqnOifnc(muf, PsiM, muM, shp);
    }
}
