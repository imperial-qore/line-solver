/**
 * @file Knessl-Tier asymptotic ray method for normalizing constant computation
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.nc;

import org.apache.commons.math3.special.Gamma;
import org.apache.commons.math3.util.FastMath;

import jline.GlobalConstants;
import jline.api.pfqn.mva.Pfqn_aql;
import jline.api.pfqn.mva.Pfqn_bs;
import jline.io.Ret;
import jline.util.matrix.Matrix;

public final class Pfqn_kt {
    private Pfqn_kt() {}

    /**
     * Knessl-Tier asymptotic expansion of the normalizing constant using the ray method.
     */
    public static Ret.pfqnNc pfqn_kt(Matrix L, Matrix N, Matrix Z) {
        Matrix Llocal = L;
        Matrix Nlocal = N;
        Matrix Zlocal = Z;
        double Ntot0 = Nlocal.sumRows(0);
        Ret.pfqnNc out = new Ret.pfqnNc(0.0, 0.0);
        if (Llocal.isEmpty() || Nlocal.isEmpty() || Ntot0 <= GlobalConstants.Zero) {
            out.G = 1.0;
            return out;
        }
        if (Zlocal.isEmpty()) {
            Zlocal = new Matrix(1, Llocal.getNumCols());
        }

        // A class with no jobs contributes a factor of 1 to G, but its saddle point is
        // u_r -> 0, where N_r*log(u_r) and N_r/u_r^2 are indeterminate and lG comes back
        // NaN. Solve the reduced model, as the self-looping branch below already does.
        int nEmpty = 0;
        for (int r = 0; r < Llocal.getNumCols(); r++) {
            if (Nlocal.get(0, r) <= GlobalConstants.Zero) nEmpty++;
        }
        if (nEmpty > 0 && nEmpty < Llocal.getNumCols()) {
            int Rk = Llocal.getNumCols() - nEmpty;
            Matrix Lk = new Matrix(Llocal.getNumRows(), Rk);
            Matrix Nk = new Matrix(1, Rk);
            Matrix Zk = new Matrix(1, Rk);
            int c = 0;
            for (int r = 0; r < Llocal.getNumCols(); r++) {
                if (Nlocal.get(0, r) <= GlobalConstants.Zero) continue;
                for (int i = 0; i < Llocal.getNumRows(); i++) Lk.set(i, c, Llocal.get(i, r));
                Nk.set(0, c, Nlocal.get(0, r));
                Zk.set(0, c, Zlocal.get(0, r));
                c++;
            }
            return pfqn_kt(Lk, Nk, Zk);
        }

        int Morig = Llocal.getNumRows();
        int Rorig = Llocal.getNumCols();
        double slcDemandFactor = 0.0;

        if (Rorig > 1) {
            boolean[] isSLC = new boolean[Rorig];
            int[] slcStation = new int[Rorig];

            for (int r = 0; r < Rorig; r++) {
                int nnz = 0;
                int firstRow = -1;
                for (int k = 0; k < Morig; k++) {
                    if (Llocal.get(k, r) > GlobalConstants.Zero) {
                        nnz++;
                        firstRow = k;
                    }
                }
                if (nnz == 1 && Zlocal.get(0, r) == 0.0) {
                    isSLC[r] = true;
                    slcStation[r] = firstRow;
                }
            }

            // classes looping at the SAME station share one (1-V)^-(1+sum N) factor
            // and contribute the multinomial (sum N)!/prod N_r!
            boolean[] done = new boolean[Rorig];
            for (int r = 0; r < Rorig; r++) {
                if (!isSLC[r] || done[r]) continue;
                int ist = slcStation[r];
                double ntot = 0.0;
                for (int s = r; s < Rorig; s++) {
                    if (!isSLC[s] || slcStation[s] != ist) continue;
                    done[s] = true;
                    ntot += Nlocal.get(0, s);
                    slcDemandFactor += Nlocal.get(0, s) * FastMath.log(Llocal.get(ist, s))
                            - Gamma.logGamma(Nlocal.get(0, s) + 1.0);
                }
                slcDemandFactor += Gamma.logGamma(ntot + 1.0);
                int rep = (int) Math.round(ntot);
                if (rep > 0) {
                    Matrix block = new Matrix(rep, Rorig);
                    for (int i = 0; i < rep; i++) {
                        for (int c = 0; c < Rorig; c++) {
                            block.set(i, c, Llocal.get(ist, c));
                        }
                    }
                    Llocal = Matrix.concatRows(Llocal, block, null);
                }
            }

            int remaining = 0;
            for (boolean b : isSLC) if (!b) remaining++;

            Matrix newL = new Matrix(Llocal.getNumRows(), remaining);
            Matrix newN = new Matrix(1, remaining);
            Matrix newZ = new Matrix(1, remaining);

            int c = 0;
            for (int r = 0; r < Rorig; r++) {
                if (isSLC[r]) continue;
                for (int k = 0; k < Llocal.getNumRows(); k++) {
                    newL.set(k, c, Llocal.get(k, r));
                }
                newN.set(0, c, Nlocal.get(0, r));
                newZ.set(0, c, Zlocal.get(0, r));
                c++;
            }
            Llocal = newL;
            Nlocal = newN;
            Zlocal = newZ;
        }

        int M = Llocal.getNumRows();
        int R = Llocal.getNumCols();
        double Ntot = R == 0 ? 0.0 : Nlocal.sumRows(0);
        if (R == 0 || Ntot <= GlobalConstants.Zero) {
            // nothing left to expand: the demand factors are the exact answer
            out.lG = slcDemandFactor;
            out.G = FastMath.exp(slcDemandFactor);
            out.X = new Matrix(1, R);
            out.Q = new Matrix(M, R);
            return out;
        }

        Ret.pfqnAMVA XQ = (Ntot <= 4.0)
                ? Pfqn_bs.pfqn_bs(Llocal, Nlocal, Zlocal)
                : Pfqn_aql.pfqn_aql(Llocal, Nlocal, Zlocal);
        Matrix Xmat = XQ.X;

        // Solve the saddle-point equations by damped Newton, starting from X:
        //   g_r(u) = u_r*(Z_r + sum_k L_kr/(1-U_k)) - N_r = 0
        double[] u = new double[R];
        double[] Zc = new double[R];
        double[] Nc = new double[R];
        for (int r = 0; r < R; r++) {
            u[r] = Xmat.get(0, r);
            Zc[r] = Zlocal.get(0, r);
            Nc[r] = Nlocal.get(0, r);
        }
        double maxUk = Double.NEGATIVE_INFINITY;
        for (int k = 0; k < M; k++) {
            double uk = 0.0;
            for (int r = 0; r < R; r++) {
                uk += Llocal.get(k, r) * u[r];
            }
            maxUk = FastMath.max(maxUk, uk);
        }
        if (maxUk >= 1.0) {
            double scale = (1.0 - 1e-6) / maxUk;
            for (int r = 0; r < R; r++) {
                u[r] *= scale;
            }
        }
        boolean converged = false;
        for (int it = 0; it < 200; it++) {
            double[] D = new double[M];
            for (int k = 0; k < M; k++) {
                double uk = 0.0;
                for (int r = 0; r < R; r++) {
                    uk += Llocal.get(k, r) * u[r];
                }
                D[k] = 1.0 / (1.0 - uk);
            }
            double[] LtD = new double[R];
            for (int r = 0; r < R; r++) {
                double sum = 0.0;
                for (int k = 0; k < M; k++) {
                    sum += Llocal.get(k, r) * D[k];
                }
                LtD[r] = sum;
            }
            double[] g = new double[R];
            double gnorm = 0.0;
            for (int r = 0; r < R; r++) {
                g[r] = u[r] * (Zc[r] + LtD[r]) - Nc[r];
                gnorm += g[r] * g[r];
            }
            gnorm = FastMath.sqrt(gnorm);
            if (gnorm <= 1e-12 * Ntot) {
                converged = true;
                break;
            }
            // J = diag(Zc + L'*D) + (u*ones(1,R)) .* (L'*(D.^2 .* L))
            Matrix J = new Matrix(R, R);
            for (int r = 0; r < R; r++) {
                for (int s = 0; s < R; s++) {
                    double A = 0.0;
                    for (int k = 0; k < M; k++) {
                        A += Llocal.get(k, r) * (D[k] * D[k]) * Llocal.get(k, s);
                    }
                    double val = u[r] * A;
                    if (r == s) {
                        val += Zc[r] + LtD[r];
                    }
                    J.set(r, s, val);
                }
            }
            Matrix rhs = new Matrix(R, 1);
            for (int r = 0; r < R; r++) {
                rhs.set(r, 0, -g[r]);
            }
            Matrix duMat = new Matrix(R, 1);
            boolean ok = Matrix.solve(J, rhs, duMat);
            if (!ok) {
                break;
            }
            double[] du = new double[R];
            for (int r = 0; r < R; r++) {
                du[r] = duMat.get(r, 0);
            }
            double alpha = 1.0;
            while (true) {
                boolean bad = false;
                for (int r = 0; r < R; r++) {
                    if (u[r] + alpha * du[r] <= 0) {
                        bad = true;
                        break;
                    }
                }
                if (!bad) {
                    double mx = Double.NEGATIVE_INFINITY;
                    for (int k = 0; k < M; k++) {
                        double uk = 0.0;
                        for (int r = 0; r < R; r++) {
                            uk += Llocal.get(k, r) * (u[r] + alpha * du[r]);
                        }
                        mx = FastMath.max(mx, uk);
                    }
                    if (mx >= 1.0) {
                        bad = true;
                    }
                }
                if (!bad) {
                    break;
                }
                alpha = alpha / 2.0;
                if (alpha < 1e-12) {
                    break;
                }
            }
            if (alpha < 1e-12) {
                break;
            }
            for (int r = 0; r < R; r++) {
                u[r] += alpha * du[r];
            }
        }
        // Choose the evaluation point: the exact saddle if Newton converged,
        // otherwise the AQL/BS throughput (stationarity limits the damage).
        double[] Dfin = new double[M];
        for (int k = 0; k < M; k++) {
            double uk = 0.0;
            for (int r = 0; r < R; r++) {
                uk += Llocal.get(k, r) * u[r];
            }
            Dfin[k] = 1.0 / (1.0 - uk);
        }
        double resnorm = 0.0;
        for (int r = 0; r < R; r++) {
            double LtD = 0.0;
            for (int k = 0; k < M; k++) {
                LtD += Llocal.get(k, r) * Dfin[k];
            }
            double gr = u[r] * (Zc[r] + LtD) - Nc[r];
            resnorm += gr * gr;
        }
        resnorm = FastMath.sqrt(resnorm);
        double[] us = new double[R];
        if (converged && resnorm <= 1e-8 * Ntot) {
            for (int r = 0; r < R; r++) {
                us[r] = u[r];
            }
        } else {
            for (int r = 0; r < R; r++) {
                us[r] = Xmat.get(0, r);
            }
        }
        // Assemble the expansion at us
        double[] Uk = new double[M];
        double[] Dus = new double[M];
        for (int k = 0; k < M; k++) {
            double uk = 0.0;
            for (int r = 0; r < R; r++) {
                uk += Llocal.get(k, r) * us[r];
            }
            Uk[k] = uk;
            Dus[k] = 1.0 / FastMath.max(GlobalConstants.FineTol, 1.0 - uk);
        }
        // H = diag(Nc./us.^2) + L'*((D.^2).*L)
        Matrix H = new Matrix(R, R);
        for (int r = 0; r < R; r++) {
            for (int s = 0; s < R; s++) {
                double A = 0.0;
                for (int k = 0; k < M; k++) {
                    A += Llocal.get(k, r) * (Dus[k] * Dus[k]) * Llocal.get(k, s);
                }
                double val = A;
                if (r == s) {
                    val += Nc[r] / (us[r] * us[r]);
                }
                H.set(r, s, val);
            }
        }
        // F = Zc'*us - sum(log(max(FineTol,1-Uk))) - Nc'*log(us)
        double F = 0.0;
        for (int r = 0; r < R; r++) {
            F += Zc[r] * us[r];
        }
        for (int k = 0; k < M; k++) {
            F -= FastMath.log(FastMath.max(GlobalConstants.FineTol, 1.0 - Uk[k]));
        }
        for (int r = 0; r < R; r++) {
            F -= Nc[r] * FastMath.log(us[r]);
        }
        double sumLogUs = 0.0;
        for (int r = 0; r < R; r++) {
            sumLogUs += FastMath.log(us[r]);
        }
        double lG = F - sumLogUs - (R / 2.0) * FastMath.log(2.0 * FastMath.PI)
                - 0.5 * logDet(H) + slcDemandFactor;

        out.lG = lG;
        out.G = FastMath.exp(lG);
        // Expose the AQL/BS throughput and queue lengths as the 3rd/4th outputs,
        // matching the MATLAB [G,lG,X,Q]=pfqn_kt signature.
        out.X = new Matrix(XQ.X);
        out.Q = new Matrix(XQ.Q);
        return out;
    }

    /**
     * Logarithm of the determinant of a symmetric positive definite matrix, from its
     * Cholesky factor. det(H) of an R x R Hessian leaves double range well before its
     * logarithm does (it overflowed at R = 64, turning lG into -Infinity).
     */
    private static double logDet(Matrix H) {
        int n = H.getNumRows();
        double[][] a = new double[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) a[i][j] = H.get(i, j);
        }
        double acc = 0.0;
        for (int k = 0; k < n; k++) {
            double d = a[k][k];
            for (int j = 0; j < k; j++) d -= a[k][j] * a[k][j];
            if (d <= 0) { // not numerically positive definite
                return FastMath.log(H.det());
            }
            double lkk = Math.sqrt(d);
            a[k][k] = lkk;
            acc += 2 * FastMath.log(lkk);
            for (int i = k + 1; i < n; i++) {
                double sum = a[i][k];
                for (int j = 0; j < k; j++) sum -= a[i][j] * a[k][j];
                a[i][k] = sum / lkk;
            }
        }
        return acc;
    }
}
