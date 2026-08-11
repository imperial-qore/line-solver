/**
 * @file Markovian Arrival Process Kronecker product composition
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import java.util.List;

import jline.util.matrix.Matrix;

public final class Map_kpc {
    private Map_kpc() {}

    /**
     * Convenience function for composing exactly two MAPs.
     */
    public static Matrix[] map_kpc(Matrix[] MAPa, Matrix[] MAPb) {
        if (MAPa.length != 2 || MAPb.length != 2) {
            throw new IllegalArgumentException("Each MAP must have exactly 2 matrices [D0, D1]");
        }

        Matrix D0 = MAPa[0].kron(MAPb[0]);
        D0.scale(-1.0);  // Negate the result

        Matrix D1 = MAPa[1].kron(MAPb[1]);

        return new Matrix[]{D0, D1};
    }

    /**
     * Convenience function for composing a list of MAPs.
     */
    public static Matrix[] map_kpc(List<Matrix[]> maps) {
        if (maps.size() < 2) {
            throw new IllegalArgumentException("Need at least 2 MAPs for composition");
        }

        Matrix[] result = map_kpc(maps.get(0), maps.get(1));

        for (int k = 2; k < maps.size(); k++) {
            result = map_kpc(result, maps.get(k));
        }

        return result;
    }
}
