/*
 * Ported from BUTools-family fluid tools (G. Horvath).
 *
 * Sojourn-time distribution of Markov-modulated fluid queues, returned as a
 * matrix-exponential (ME) or phase-type (PH) representation.
 *
 * Reference: Horvath G, Telek M, "Sojourn times in fluid queues with
 * independent and dependent input and output processes", Performance
 * Evaluation 79:160-181, 2014.
 */
package jline.lib.butools.queues;

import jline.lib.butools.mam.FluidTools;
import jline.lib.butools.mam.GeneralFluidSolution;
import jline.lib.butools.mam.GeneralFluidSolve;
import jline.lib.butools.mc.CTMCSolve;
import jline.lib.butools.reptrans.TransformToOnes;
import jline.util.matrix.Matrix;

public final class FluidQueueSTD {
    private FluidQueueSTD() {}

    /**
     * Sojourn-time distribution of a fluid queue with input rate matrix Rin and
     * output (service) rate matrix Rout, modulated by generator Q.
     *
     * @return {alpha, A}: an ME (transToPH=false) or PH (transToPH=true)
     *         representation of the sojourn time.
     */
    public static Matrix[] fluidQueueSTD(Matrix Q, Matrix Rin, Matrix Rout, Matrix Q0, boolean transToPH) {
        int N = Q.getNumRows();
        GeneralFluidSolution sol = GeneralFluidSolve.generalFluidSolve(Q, Rin.sub(Rout), Q0, 1e-14);
        Matrix mass0 = sol.getMass0(), ini = sol.getIni(), K = sol.getK(), clo = sol.getClo();
        int nk = K.getNumRows();

        Matrix iniKi = K.transpose().leftMatrixDivide(ini.transpose().scale(-1.0)).transpose(); // ini*inv(-K)
        double lambda = mass0.mult(Rin).add(iniKi.mult(clo).mult(Rin)).elementSum();

        Matrix alpha, A;
        if (transToPH) {
            Matrix Delta = FluidTools.diagFrom(iniKi.scale(1.0 / lambda));
            Matrix reshaped = clo.mult(Rin).columnMajorOrder().transpose(); // reshape(clo*Rin,1,N*nk)
            alpha = reshaped.mult(Matrix.eye(N).kron(Delta));
            A = Rout.kron(Delta.inv().mult(K.transpose()).mult(Delta)).add(Q.kron(Matrix.eye(nk)));
        } else {
            Matrix col = K.scale(-1.0).pinv().mult(clo).mult(Rin).columnMajorOrder(); // reshape(inv(-K)*clo*Rin,N*nk,1)
            Matrix B = TransformToOnes.transformToOnes(col);
            Matrix Bi = B.pinv();  // B is triangular with tiny determinant but full rank
            alpha = Matrix.ones(1, N).kron(ini.scale(1.0 / lambda)).mult(Bi);
            A = B.mult(Q.transpose().kron(Matrix.eye(nk)).add(Rout.kron(K))).mult(Bi);
        }
        return new Matrix[]{alpha, A};
    }

    /**
     * Sojourn-time distribution of a fluid queue in which both the arrival and
     * the service processes are Markov-modulated fluid flows. If srv0stop is
     * true the service stops while the server fluid level is zero.
     *
     * @return {alpha, A}: ME (transToPH=false) or PH (transToPH=true) sojourn time.
     */
    public static Matrix[] fluFluSTD(Matrix Qin, Matrix Rin, Matrix Qout, Matrix Rout, boolean srv0stop, boolean transToPH) {
        Matrix Iin = Matrix.eye(Qin.getNumRows());
        Matrix Iout = Matrix.eye(Qout.getNumRows());
        Matrix Rh = Rin.kron(Iout).sub(Iin.kron(Rout));
        Matrix Qh = Qin.kron(Rout).add(Rin.kron(Qout));
        GeneralFluidSolution sol = GeneralFluidSolve.generalFluidSolve(Qh, Rh, null, 1e-14);
        Matrix inih = sol.getIni(), Kh = sol.getK(), cloh = sol.getClo();

        double lambda = CTMCSolve.ctmcSolve(Qin).mult(Rin).elementSum();
        double mu = CTMCSolve.ctmcSolve(Qout).mult(Rout).elementSum();

        Matrix alpha, A;
        if (transToPH) {
            Matrix Delta = FluidTools.diagFrom(Kh.transpose().leftMatrixDivide(inih.transpose().scale(-1.0))); // diag(inih*inv(-Kh))
            A = Delta.inv().mult(Kh.transpose()).mult(Delta);
            Matrix M = srv0stop
                    ? Delta.mult(cloh).mult(Rin.kron(Rout)).scale(1.0 / lambda / mu)
                    : Delta.mult(cloh).mult(Rin.kron(Iout)).scale(1.0 / lambda);
            alpha = FluidTools.rowSums(M).transpose();
        } else {
            Matrix col = srv0stop
                    ? FluidTools.rowSums(cloh.mult(Rin.kron(Rout)).scale(1.0 / lambda / mu))
                    : FluidTools.rowSums(cloh.mult(Rin.kron(Iout)).scale(1.0 / lambda));
            Matrix B = TransformToOnes.transformToOnes(col);
            Matrix iB = B.pinv();  // B is triangular with tiny determinant but full rank
            A = B.mult(Kh).mult(iB);
            alpha = inih.mult(Kh.scale(-1.0).pinv()).mult(iB);
        }
        return new Matrix[]{alpha, A};
    }
}
