/**
 * @file Markovian Arrival Process superposition operations
 *
 * Creates superposition of MAP processes using Kronecker product techniques.
 * Fundamental for modeling independent arrival stream combinations in queueing networks.
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import jline.util.matrix.MatrixCell;

public final class Map_super {
    private Map_super() {}

    /**
     * Creates a superposition of two Markovian Arrival Processes (MAPs) to form a new MAP.
     *
     * @param MAPa The first Markovian Arrival Process stored in a MatrixCell.
     * @param MAPb The second Markovian Arrival Process stored in a MatrixCell.
     * @return A MatrixCell representing the superposed MAP, formed by combining the input MAPs.
     */
    public static MatrixCell map_super(MatrixCell MAPa, MatrixCell MAPb) {
        MatrixCell sup = new MatrixCell();

        sup.set(0, MAPa.get(0).krons(MAPb.get(0)));
        sup.set(1, MAPa.get(1).krons(MAPb.get(1)));

        return Map_normalize.map_normalize(sup);
    }
}
