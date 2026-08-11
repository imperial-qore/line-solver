/**
 * @file CTMC stochastic complementarity analysis
 *
 * @since LINE 3.0
 */
package jline.api.mc;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;

import org.ejml.data.DMatrix;
import org.ejml.data.DMatrixRMaj;
import org.ejml.data.DMatrixSparseCSC;
import org.ejml.dense.row.CommonOps_DDRM;
import org.ejml.dense.row.factory.LinearSolverFactory_DDRM;
import org.ejml.interfaces.linsol.LinearSolverDense;
import org.ejml.ops.DConvertMatrixStruct;

import jline.solvers.ctmc.SolverCTMC;
import jline.util.matrix.Matrix;

public final class Ctmc_stochcomp {
    private Ctmc_stochcomp() {}

    /** Pivot-ratio below which the LU of the complement block is treated as singular. */
    private static final double LU_MIN_QUALITY = 1e-12;

    public static SolverCTMC.StochCompResult ctmc_stochcomp(Matrix Q, List<Double> I_list) {
        Matrix I = new Matrix(I_list);
        if (I_list.isEmpty()) {
            int halfSize = (int) Math.ceil((double) Q.getNumCols() / 2);
            List<Double> defaultList = new ArrayList<Double>();
            for (int idx = 0; idx < halfSize; idx++) {
                defaultList.add((double) idx);
            }
            I = new Matrix(defaultList);
        }

        HashSet<Integer> iSet = new HashSet<Integer>(I.getNumRows() * 2);
        for (int iRow = 0; iRow < I.getNumRows(); iRow++) {
            iSet.add((int) I.get(iRow));
        }
        List<Double> diff_values = new ArrayList<Double>();
        for (int checkValue = 0; checkValue < Q.getNumCols(); checkValue++) {
            if (!iSet.contains(checkValue)) {
                diff_values.add((double) checkValue);
            }
        }
        Matrix Ic = new Matrix(diff_values);

        Matrix Q11 = Q.getSubMatrix(I, I);
        Matrix Q12 = Q.getSubMatrix(I, Ic);
        Matrix Q21 = Q.getSubMatrix(Ic, I);
        Matrix Q22 = Q.getSubMatrix(Ic, Ic);

        int n = Q22.getNumRows();
        int nrhs = Q21.getNumCols();

        // see _kb/03-api-layer.md for rationale
        if (n > Ctmc_solve.GMRES_MIN_STATES) {
            Matrix negQ22 = Q22.neg();
            // A vanishing diagonal leaves the block singular; substituting one matches
            // what the dense path below does before factorizing.
            for (int i = 0; i < n; i++) {
                if (Math.abs(negQ22.get(i, i)) < 1e-10) negQ22.set(i, i, 1.0);
            }
            Matrix Tit = Ctmc_gmres.ctmc_gmres(negQ22, Q21, 0.0, 0, 0);
            if (Tit != null) {
                Matrix Tprod = Q12.mult(Tit);
                Matrix Sit = Q11.add(1.0, Tprod);
                return new SolverCTMC.StochCompResult(Sit, Q11, Q12, Q21, Q22, Tprod);
            }
        }

        DMatrixRMaj denseNegQ22 = new DMatrixRMaj(n, n);
        DConvertMatrixStruct.convert(Q22.toDMatrixSparseCSC(), denseNegQ22);
        CommonOps_DDRM.scale(-1.0, denseNegQ22);

        for (int i = 0; i < n; i++) {
            if (Math.abs(denseNegQ22.get(i, i)) < 1e-10) {
                denseNegQ22.set(i, i, 1.0);
            }
        }

        LinearSolverDense<DMatrixRMaj> luSolver = LinearSolverFactory_DDRM.lu(n);
        // see _kb/03-api-layer.md for rationale
        boolean luOk = luSolver.setA(denseNegQ22) && luSolver.quality() > LU_MIN_QUALITY;

        DMatrixRMaj denseQ21 = new DMatrixRMaj(n, nrhs);
        DConvertMatrixStruct.convert(Q21.toDMatrixSparseCSC(), denseQ21);
        DMatrixRMaj denseT = new DMatrixRMaj(n, nrhs);
        if (luOk) {
            luSolver.solve(denseQ21, denseT);
        } else {
            // see _kb/03-api-layer.md for rationale
            LinearSolverDense<DMatrixRMaj> svdSolver = LinearSolverFactory_DDRM.pseudoInverse(true);
            svdSolver.setA(denseNegQ22);
            svdSolver.solve(denseQ21, denseT);
        }
        Matrix T = new Matrix((DMatrix) DConvertMatrixStruct.convert(denseT, (DMatrixSparseCSC) null, 1e-15));
        T = Q12.mult(T);
        Matrix S = Q11.add(1.0, T);
        SolverCTMC.StochCompResult result = new SolverCTMC.StochCompResult(S, Q11, Q12, Q21, Q22, T);
        if (luOk) {
            result.denseLUSolver = luSolver;
        }
        return result;
    }
}
