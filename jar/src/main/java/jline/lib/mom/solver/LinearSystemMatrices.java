package jline.lib.mom.solver;

import java.util.Objects;

import org.apache.commons.math3.linear.RealMatrix;

/**
 * Container for the linear system matrices used in MOM solver
 *
 * @property C Main coefficient matrix (block diagonal structure)
 * @property Cg Coupling matrix for previous populations
 * @property D Right-hand side matrix
 * @property Dr Special matrix for recursion
 */
public final class LinearSystemMatrices {
    public final RealMatrix C;
    public final RealMatrix Cg;
    public final RealMatrix D;
    public final RealMatrix Dr;

    public LinearSystemMatrices(RealMatrix C, RealMatrix Cg, RealMatrix D, RealMatrix Dr) {
        this.C = C;
        this.Cg = Cg;
        this.D = D;
        this.Dr = Dr;
    }

    public RealMatrix getC() {
        return C;
    }

    public RealMatrix getCg() {
        return Cg;
    }

    public RealMatrix getD() {
        return D;
    }

    public RealMatrix getDr() {
        return Dr;
    }

    public RealMatrix component1() {
        return C;
    }

    public RealMatrix component2() {
        return Cg;
    }

    public RealMatrix component3() {
        return D;
    }

    public RealMatrix component4() {
        return Dr;
    }

    public LinearSystemMatrices copy(RealMatrix C, RealMatrix Cg, RealMatrix D, RealMatrix Dr) {
        return new LinearSystemMatrices(C, Cg, D, Dr);
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof LinearSystemMatrices)) return false;
        LinearSystemMatrices that = (LinearSystemMatrices) o;
        return Objects.equals(C, that.C) && Objects.equals(Cg, that.Cg)
                && Objects.equals(D, that.D) && Objects.equals(Dr, that.Dr);
    }

    @Override
    public int hashCode() {
        return Objects.hash(C, Cg, D, Dr);
    }

    @Override
    public String toString() {
        return "LinearSystemMatrices(C=" + C + ", Cg=" + Cg + ", D=" + D + ", Dr=" + Dr + ")";
    }
}
