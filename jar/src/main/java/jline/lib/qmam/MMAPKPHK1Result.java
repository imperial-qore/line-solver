package jline.lib.qmam;

import java.util.List;

import jline.util.matrix.Matrix;

/** Result of MMAP[K]/PH[K]/1 queue analysis. */
public final class MMAPKPHK1Result {
    public final List<Matrix> qlPerType;
    public final Matrix qlTotal;
    public final List<Matrix> sojAlpha;
    public final List<Matrix> waitAlpha;
    public final Matrix Smat; // nullable

    public MMAPKPHK1Result(List<Matrix> qlPerType, Matrix qlTotal,
                           List<Matrix> sojAlpha, List<Matrix> waitAlpha, Matrix Smat) {
        this.qlPerType = qlPerType;
        this.qlTotal = qlTotal;
        this.sojAlpha = sojAlpha;
        this.waitAlpha = waitAlpha;
        this.Smat = Smat;
    }

    public List<Matrix> getQlPerType() { return qlPerType; }
    public Matrix getQlTotal() { return qlTotal; }
    public List<Matrix> getSojAlpha() { return sojAlpha; }
    public List<Matrix> getWaitAlpha() { return waitAlpha; }
    public Matrix getSmat() { return Smat; }
}
