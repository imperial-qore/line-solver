/*
 * Ported from BUTools-family fluid tools (G. Horvath).
 */
package jline.lib.butools.mam;

import java.util.List;

import jline.util.matrix.Matrix;

/**
 * Matrix-exponential building blocks of a first/second-order level-dependent
 * fluid queue, as produced by SecondOrderLevelDependentFluidSolve and consumed
 * by LevelDependentFluidStationary.
 */
public final class LevelDependentFluidSolution {
    public final List<Matrix> masses; // K+1 point-mass row vectors
    public final List<Matrix> iniF, KF, cloF; // forward parameters per regime
    public final List<Matrix> iniB, KB, cloB; // backward parameters per regime

    public LevelDependentFluidSolution(List<Matrix> masses,
                                       List<Matrix> iniF, List<Matrix> KF, List<Matrix> cloF,
                                       List<Matrix> iniB, List<Matrix> KB, List<Matrix> cloB) {
        this.masses = masses;
        this.iniF = iniF; this.KF = KF; this.cloF = cloF;
        this.iniB = iniB; this.KB = KB; this.cloB = cloB;
    }
}
