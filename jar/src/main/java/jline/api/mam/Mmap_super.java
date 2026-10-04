/**
 * @file Marked Markovian Arrival Process superposition operations
 *
 * Combines multiple MMAP processes into superposed multiclass arrival streams.
 * Fundamental for modeling aggregated traffic sources in complex queueing networks.
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Mmap_super {
    private Mmap_super() {}

    /**
     * Combines two MMAPs into one superposed MMAP.
     */
    public static MatrixCell mmap_super(MatrixCell MMAPa, MatrixCell MMAPb, String opt) {
        MatrixCell sup = new MatrixCell();
        if ("default".equals(opt)) {
            int K1 = MMAPa.size() - 2;
            int K2 = MMAPb.size() - 2;
            int n1 = MMAPa.get(0).length();
            int n2 = MMAPb.get(0).length();

            sup.set(0, MMAPa.get(0).krons(MMAPb.get(0)));
            sup.set(1, MMAPa.get(1).krons(MMAPb.get(1)));

            for (int i = 0; i < K1; i++) {
                sup.set(2 + i, MMAPa.get(2 + i).krons(new Matrix(n2, n2, 0)));
            }

            for (int j = 0; j < K2; j++) {
                Matrix a = new Matrix(n1, n1, 0);
                sup.set(2 + K1 + j, a.krons(MMAPb.get(2 + j)));
            }
        } else if ("match".equals(opt)) {
            int K1 = MMAPa.size();
            int K2 = MMAPb.size();

            if (K1 != K2) {
                throw new RuntimeException("class matching failed: MMAPs have different number of classes");
            }

            for (int i = 0; i < K1; i++) {
                // class c of both MMAPs maps to class c of the superposition
                // (MATLAB: krons(MMAPa{i}, MMAPb{i}); kron-summing MMAPa with
                // itself was a translation typo)
                sup.set(i, MMAPa.get(i).krons(MMAPb.get(i)));
            }
        } else {
            throw new RuntimeException("unrecognized option");
        }

        return Mmap_normalize.mmap_normalize(sup);
    }

    /**
     * Combines two MMAPs into one superposed MMAP using the default option.
     */
    public static MatrixCell mmap_super(MatrixCell MMAPa, MatrixCell MMAPb) {
        return mmap_super(MMAPa, MMAPb, "default");
    }

    /**
     * Combines a list of MMAPs into one superposed MMAP.
     */
    public static MatrixCell mmap_super(MatrixCell MMAPa) {
        MMAPa.removeNull();
        MatrixCell SUP = new MatrixCell();
        SUP.set(0, MMAPa.get(0));
        for (int i = 1; i < MMAPa.size(); i++) {
            MatrixCell MMAPb = new MatrixCell();
            MMAPb.set(0, MMAPa.get(i));
            SUP = mmap_super(SUP, MMAPb);
        }
        return SUP;
    }
}
