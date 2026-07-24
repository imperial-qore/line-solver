/**
 * @file M3A MMAP Compression
 *
 * @since LINE 3.0
 */
package jline.lib.m3a;

import jline.api.mam.Mamap2m_fit_gamma_fb_mmap;
import jline.api.mam.Mmap_isfeasible;
import jline.util.matrix.MatrixCell;

public final class M3aCompress {
    private M3aCompress() {}

    /**
     * Compresses a Marked Markovian Arrival Process (MMAP) using M3A fitting.
     */
    public static MatrixCell m3afit_compress(MatrixCell mmap, M3aCompressOptions options) {
        if (mmap.size() < 3) {
            throw new IllegalArgumentException(
                    "Input MMAP must have at least 3 matrices (D0, D1, and at least one class matrix)");
        }

        int numClasses = mmap.size() - 2;

        if (options.isVerbose()) {
            String mmapType = options.getNumStates() + "-state AMAP[" + numClasses + "]";
            System.out.println("Init: M3A will search for a " + mmapType);
        }

        MatrixCell result;
        switch (options.getMethod()) {
            case AMAP_2STATE:
                result = Mamap2m_fit_gamma_fb_mmap.mamap2m_fit_gamma_fb_mmap(mmap);
                break;
            default:
                throw new IllegalArgumentException("Unknown compression method: " + options.getMethod());
        }

        if (options.isVerbose()) {
            String resultType = result.get(0).getNumRows() + "-state M3PP[" + numClasses + "]";
            if (Mmap_isfeasible.mmap_isfeasible(result)) {
                System.out.println("Output: M3A found a valid " + resultType + ".");
            } else {
                System.out.println("Output: M3A could *not* obtain a valid MMAP.");
            }
        }

        return result;
    }

    public static MatrixCell m3afit_compress(MatrixCell mmap) {
        return m3afit_compress(mmap, new M3aCompressOptions());
    }
}
