/**
 * @file M3A automatic fitting functions
 *
 * @since LINE 3.0
 */
package jline.lib.m3a;

import java.util.HashSet;
import java.util.Set;

import jline.GlobalConstants;
import jline.VerboseLevel;
import jline.api.mam.Mamap22_fit_gamma_fs_trace;
import jline.api.mam.Mamap2m_fit_gamma_fb_trace;
import jline.api.mam.Mmap_isfeasible;
import jline.api.mam.m3pp.M3pp2m_fitc_trace;
import jline.api.mam.m3pp.M3pp_superpos_fitc_trace;
import jline.util.Pair;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class M3aFit {
    private M3aFit() {}

    private static boolean isVerbose() {
        return GlobalConstants.getVerbose() != VerboseLevel.SILENT;
    }

    /**
     * Prepares multiclass trace for M3A fitting.
     */
    public static MTrace m3afit_init(double[] S, int[] C) {
        Set<Integer> distinct = new HashSet<Integer>();
        for (int c : C) distinct.add(c);
        int numClasses = distinct.size();
        return new MTrace(S, C, numClasses);
    }

    /**
     * Prepares multiclass trace for M3A fitting from Matrix inputs.
     */
    public static MTrace m3afit_init(Matrix S, Matrix C) {
        double[] sArray = S.toArray1D();
        double[] cVals = C.toArray1D();
        int[] cArray = new int[C.getNumRows() * C.getNumCols()];
        for (int i = 0; i < cArray.length; i++) {
            cArray[i] = (int) cVals[i];
        }
        return m3afit_init(sArray, cArray);
    }

    /**
     * Automatic fitting of trace into a Marked Markovian Arrival Process.
     */
    public static MatrixCell m3afit_auto(MTrace mtrace, M3aFitOptions options) {
        double avg = 0.0;
        for (double s : mtrace.getS()) avg += s;
        avg /= mtrace.getS().length;
        double timescale = (options.getTimescale() != null) ? options.getTimescale() : 10 * avg;

        double sumS = 0.0;
        for (double s : mtrace.getS()) sumS += s;
        double timescaleAsy = (options.getTimescaleAsy() != null) ? options.getTimescaleAsy()
                : Math.max(10 * timescale, (sumS - mtrace.getS()[0]) / 100);

        String mmapType;
        MatrixCell mmap;

        if (mtrace.getNumClasses() == 2 && options.getNumStates() == 2 && options.getMethod() == 0) {
            mmapType = options.getNumStates() + "-state MAMAP[" + mtrace.getNumClasses() + "]";
            if (isVerbose()) System.out.println("Init: M3A algorithm will search for a " + mmapType);
            mmap = Mamap22_fit_gamma_fs_trace.mamap22_fit_gamma_fs_trace(mtrace.getS(), mtrace.getC());
        } else if (mtrace.getNumClasses() > 2 && options.getNumStates() >= 2 && options.getMethod() == 0) {
            mmapType = options.getNumStates() + "-state MAMAP[" + mtrace.getNumClasses() + "]";
            if (isVerbose()) System.out.println("Init: M3A algorithm will search for a " + mmapType);
            mmap = Mamap2m_fit_gamma_fb_trace.mamap2m_fit_gamma_fb_trace(mtrace.getS(), mtrace.getC());
        } else if (mtrace.getNumClasses() >= 2 && options.getNumStates() == 2 && options.getMethod() == 1) {
            mmapType = options.getNumStates() + "-state M3PP[" + mtrace.getNumClasses() + "]";
            if (isVerbose()) System.out.println("Init: M3A algorithm will search for a " + mmapType);
            Matrix[] result = M3pp2m_fitc_trace.m3pp2m_fitc_trace(
                    mtrace.getS(), mtrace.getC(), "approx_ag", timescale, timescaleAsy);
            mmap = arrayToMatrixCell(result);
        } else if (mtrace.getNumClasses() >= 2 && options.getNumStates() > 2 && options.getMethod() == 1) {
            mmapType = options.getNumStates() + "-state M3PP[" + mtrace.getNumClasses() + "]";
            if (isVerbose()) System.out.println("Init: M3A algorithm will search for a " + mmapType);
            Pair<Matrix[], ?> sup = M3pp_superpos_fitc_trace.m3pp_superpos_fitc_trace(
                    mtrace.getS(), mtrace.getC(), timescale, timescaleAsy);
            Matrix[] result = sup.getLeft();
            mmap = arrayToMatrixCell(result);
        } else {
            if (isVerbose()) System.out.println("Output: M3A algorithm could *not* obtain a valid MMAP.");
            return null;
        }

        String mmapResType = mmap.get(0).getNumRows() + "-state M3PP[" + (mmap.size() - 2) + "]";
        if (isVerbose()) {
            if (Mmap_isfeasible.mmap_isfeasible(mmap)) {
                System.out.println("Output: M3A algorithm found a valid " + mmapResType + ".");
            } else {
                System.out.println("Output: M3A algorithm could *not* obtain a valid MMAP.");
            }
        }

        return mmap;
    }

    /**
     * Automatic fitting with simple parameters.
     */
    public static MatrixCell m3afit_auto(double[] S, int[] C, int numStates, int method) {
        MTrace mtrace = m3afit_init(S, C);
        M3aFitOptions options = new M3aFitOptions(method, numStates, null, null);
        return m3afit_auto(mtrace, options);
    }

    public static MatrixCell m3afit_auto(double[] S, int[] C, int numStates) {
        return m3afit_auto(S, C, numStates, 1);
    }

    private static MatrixCell arrayToMatrixCell(Matrix[] array) {
        MatrixCell result = new MatrixCell(array.length);
        for (int i = 0; i < array.length; i++) {
            result.set(i, array[i]);
        }
        return result;
    }
}
