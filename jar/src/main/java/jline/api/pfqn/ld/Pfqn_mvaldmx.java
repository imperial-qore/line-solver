/**
 * @file Load-dependent MVA for mixed open-closed queueing networks
 *
 * Implements Mean Value Analysis for mixed queueing networks with both open and closed classes
 * and load-dependent service rates.
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.ld;

import java.util.ArrayList;
import java.util.Arrays;

import org.apache.commons.math3.util.FastMath;

import jline.io.Ret;
import jline.util.Maths;
import jline.util.PopulationLattice;
import jline.util.Utils;
import jline.util.matrix.Matrix;

public final class Pfqn_mvaldmx {
    private Pfqn_mvaldmx() {}

    /**
     * MVA method for mixed queueing networks with load-dependent nodes.
     */
    public static Ret.pfqnMVALDMX pfqn_mvaldmx(Matrix lambda, Matrix D, Matrix N, Matrix Z,
                                                Matrix mu, Matrix S) {
        double NfiniteSum = 0.0;
        for (int i = 0; i < N.getNumRows(); i++) {
            for (int j = 0; j < N.getNumCols(); j++) {
                double num = N.get(i, j);
                if (Double.isFinite(num)) {
                    NfiniteSum += num;
                }
            }
        }
        if (mu == null && S == null) {
            mu = Matrix.ones(D.getNumRows(), (int) NfiniteSum);
            S = Matrix.ones(D.getNumRows(), 1);
        }
        if (mu.getNumCols() < NfiniteSum) {
            throw new RuntimeException("pfqn_mvaldmx: MVALDMX requires to specify the load-dependent rates "
                    + "with one job more than the maximum closed population.");
        }
        Matrix f = lambda.find();
        for (int i = 0; i < f.getNumRows(); i++) {
            double num = N.get((int) f.get(i));
            if (num > 0 && Double.isFinite(num)) {
                throw new RuntimeException("pfqn_mvaldmx: Arrival rate cannot be specified on closed classes.");
            }
        }
        int M = D.getNumRows();
        int R = D.getNumCols();
        ArrayList<Integer> openClasses = new ArrayList<Integer>();
        ArrayList<Integer> closedClasses = new ArrayList<Integer>();
        for (int i = 0; i < N.length(); i++) {
            if (Utils.isInf(N.get(i))) {
                openClasses.add(i);
            } else {
                closedClasses.add(i);
            }
        }
        Matrix XN = new Matrix(1, R);
        Matrix UN = new Matrix(M, R);
        Matrix CN = new Matrix(M, R);
        Matrix QN = new Matrix(M, R);
        double lGN = 0.0;
        Matrix newMu = new Matrix(mu.getNumRows(), mu.getNumCols() + 1);
        for (int i = 0; i < newMu.getNumRows(); i++) {
            for (int j = 0; j < newMu.getNumCols(); j++) {
                if (j < mu.getNumCols()) {
                    newMu.set(i, j, mu.get(i, j));
                } else {
                    newMu.set(i, j, mu.get(i, j - 1));
                }
            }
        }
        mu = newMu;
        Ret.pfqnLDMXEC ret1 = Pfqn_ldmx_ec.pfqn_ldmx_ec(lambda, D, new Matrix(mu));
        Matrix EC = ret1.EC;
        Matrix E = ret1.E;
        Matrix Eprime = ret1.Eprime;
        int C = closedClasses.size();
        if (C == 0) {
            // Purely open network: the closed lattice is the single empty
            // population, so P(station holds 0 closed jobs)=1 and the sums below
            // collapse to their n=0 term, which is the effective capacity. The
            // lattice recursion cannot run it -- pprod on an empty bound never
            // returns the -1 sentinel -- so it is evaluated in closed form here.
            Matrix Pc0 = Matrix.ones(M, 1);
            for (int r : openClasses) {
                XN.set(r, lambda.get(r));
                for (int ist = 0; ist < M; ist++) {
                    QN.set(ist, r, lambda.get(r) * D.get(ist, r) * EC.get(ist, 0));
                    CN.set(ist, r, D.get(ist, r) * EC.get(ist, 0));
                    UN.set(ist, r, lambda.get(r) * Eprime.get(ist, 1) / E.get(ist, 1));
                }
            }
            return new Ret.pfqnMVALDMX(XN, QN, UN, CN, lGN, Pc0);
        }
        Matrix Dc = new Matrix(D.getNumRows(), C);
        Matrix Nc = new Matrix(1, C);
        Matrix Zc = new Matrix(1, C);
        for (int i = 0; i < C; i++) {
            int c = closedClasses.get(i);
            for (int j = 0; j < Dc.getNumRows(); j++) {
                Dc.set(j, i, D.get(j, c));
            }
            Nc.set(0, i, N.get(c));
            Zc.set(0, i, Z.get(c));
        }
        Matrix prods = new Matrix(1, C);
        for (int r = 0; r < C; r++) {
            double prod = 1.0;
            for (int i = 0; i < r; i++) {
                prod *= (Nc.get(i) + 1);
            }
            prods.set(0, r, prod);
        }
        Matrix nvec = PopulationLattice.pprod(Nc);
        Matrix[] Pc = new Matrix[M];
        double ncProd = 1.0;
        for (int i = 0; i < C; i++) {
            ncProd *= (1 + Nc.get(i));
        }
        for (int ist = 0; ist < M; ist++) {
            Pc[ist] = new Matrix((int) (1 + Nc.elementSum()), (int) ncProd);
        }
        Matrix x = new Matrix(C, (int) ncProd);
        Matrix[] w = new Matrix[M];
        for (int ist = 0; ist < M; ist++) {
            w[ist] = new Matrix(C, (int) ncProd);
        }
        for (int ist = 0; ist < M; ist++) {
            Pc[ist].set(0, PopulationLattice.hashpop(nvec, Nc, C, prods), 1);
        }
        Matrix u = new Matrix(M, C);

        while (!(nvec.getNumRows() == 1 && nvec.getNumCols() == 1 && nvec.value() == -1.0)) {
            int hnvec = PopulationLattice.hashpop(nvec, Nc, C, prods);
            double nc = nvec.elementSum();
            for (int ist = 0; ist < M; ist++) {
                for (int c = 0; c < C; c++) {
                    if (nvec.get(c) > 0) {
                        int hnvec_c = PopulationLattice.hashpop(
                                Matrix.oner(nvec, new ArrayList<Integer>(Arrays.asList(Integer.valueOf(c)))),
                                Nc, C, prods);
                        int n = 0;
                        while (n < nc) {
                            w[ist].set(c, hnvec, w[ist].get(c, hnvec)
                                    + Dc.get(ist, c) * (n + 1) * EC.get(ist, n) * Pc[ist].get(n, hnvec_c));
                            n++;
                        }
                    }
                }
            }
            for (int c = 0; c < C; c++) {
                double sumw = 0.0;
                for (int ist = 0; ist < M; ist++) {
                    sumw += w[ist].get(c, hnvec);
                }
                x.set(c, hnvec, nvec.get(c) / (Zc.get(c) + sumw));
            }
            for (int ist = 0; ist < M; ist++) {
                int n = 0;
                while (n < nc) {
                    for (int c = 0; c < C; c++) {
                        if (nvec.get(c) > 0) {
                            int hnvec_c = PopulationLattice.hashpop(
                                    Matrix.oner(nvec, new ArrayList<Integer>(Arrays.asList(Integer.valueOf(c)))),
                                    Nc, C, prods);
                            Pc[ist].set(1 + n, hnvec, Pc[ist].get(1 + n, hnvec)
                                    + Dc.get(ist, c) * EC.get(ist, n) * x.get(c, hnvec) * Pc[ist].get(n, hnvec_c));
                        }
                    }
                    n++;
                }
                double sumpc = 0.0;
                int k = 0;
                while (k < nc) {
                    sumpc += Pc[ist].get(1 + k, hnvec);
                    k++;
                }
                Pc[ist].set(0, hnvec, Maths.max(Math.ulp(1.0), 1 - sumpc));
            }
            Matrix nvecFind = nvec.find();
            if (!nvecFind.isEmpty()) {
                int idx = nvecFind.getNumRows() - 1;
                while (idx >= 0 && nvec.get((int) nvecFind.get(idx)) <= 0) {
                    idx--;
                }
                int last_nnz = (int) nvecFind.get(idx);
                double sumnvec = 0.0;
                double sumnc = 0.0;
                double sumnvecp = 0.0;
                for (int i = 0; i < C; i++) {
                    if (i < last_nnz) {
                        sumnvec += nvec.get(i);
                        sumnc += Nc.get(i);
                    } else if (i > last_nnz) {
                        sumnvecp += nvec.get(i);
                    }
                }
                if (sumnvec == sumnc && sumnvecp == 0.0) {
                    // see _kb/03-api-layer.md for rationale
                    double xref = x.get(last_nnz, hnvec);
                    if (xref > 0) {
                        lGN -= FastMath.log(xref);
                    }
                }
            }
            nvec = PopulationLattice.pprod(nvec, Nc);
        }
        int hnvec = PopulationLattice.hashpop(Nc, Nc, C, prods);
        for (int c = 0; c < C; c++) {
            int hnvec_c = PopulationLattice.hashpop(
                    Matrix.oner(Nc, new ArrayList<Integer>(Arrays.asList(Integer.valueOf(c)))),
                    Nc, C, prods);
            for (int ist = 0; ist < M; ist++) {
                u.set(ist, c, 0);
                double sumNc = Nc.elementSum();
                int n = 0;
                while (n < sumNc) {
                    u.set(ist, c, u.get(ist, c)
                            + (Dc.get(ist, c) * x.get(c, hnvec) * Eprime.get(ist, n)) / E.get(ist, n)
                                    * Pc[ist].get(n, hnvec_c));
                    n++;
                }
            }
        }
        for (int i = 0; i < C; i++) {
            XN.set(closedClasses.get(i), x.get(i, hnvec));
            for (int j = 0; j < M; j++) {
                UN.set(j, closedClasses.get(i), u.get(j, i));
                CN.set(j, closedClasses.get(i), w[j].get(i, hnvec));
                QN.set(j, closedClasses.get(i),
                        XN.get(closedClasses.get(i)) * CN.get(j, closedClasses.get(i)));
            }
        }
        for (int r : openClasses) {
            XN.set(r, lambda.get(r));
            for (int ist = 0; ist < M; ist++) {
                QN.set(ist, r, 0);
                double sumNc = Nc.elementSum();
                {
                    int n = 0;
                    while (n <= sumNc) {
                        QN.set(ist, r, QN.get(ist, r)
                                + lambda.get(r) * D.get(ist, r) * (n + 1) * EC.get(ist, n) * Pc[ist].get(n, hnvec));
                        n++;
                    }
                }
                CN.set(ist, r, QN.get(ist, r) / lambda.get(r));
                UN.set(ist, r, 0);
                int n = -1;
                while (n < sumNc) {
                    UN.set(ist, r, UN.get(ist, r)
                            + lambda.get(r) * Eprime.get(ist, n + 2) / E.get(ist, n + 2) * Pc[ist].get(n + 1, hnvec));
                    n++;
                }
            }
        }
        Matrix newPc = new Matrix(M, (int) (1 + Nc.elementSum()));
        for (int ist = 0; ist < M; ist++) {
            for (int j = 0; j < newPc.getNumCols(); j++) {
                newPc.set(ist, j, Pc[ist].get(j, hnvec));
            }
        }
        return new Ret.pfqnMVALDMX(XN, QN, UN, CN, lGN, newPc);
    }
}
