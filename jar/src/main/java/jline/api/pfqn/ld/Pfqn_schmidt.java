package jline.api.pfqn.ld;

import java.util.ArrayList;
import java.util.List;

import jline.GlobalConstants;
import jline.io.Ret;
import jline.lang.constant.SchedStrategy;
import jline.util.PopulationLattice;
import jline.util.matrix.Matrix;
import org.apache.commons.math3.util.FastMath;

/**
 * Schmidt method for load-dependent MVA with multi-server stations.
 *
 * Reference: R. Schmidt, "An approximate MVA algorithm for exponential, class-dependent
 * multiple server stations," Performance Evaluation, vol. 29, no. 4, pp. 245-254, 1997.
 */
public final class Pfqn_schmidt {
    private Pfqn_schmidt() {}

    public static Ret.pfqnSchmidt pfqn_schmidt(Matrix D, Matrix N, Matrix S, List<SchedStrategy> sched) {
        long startTime = System.nanoTime();

        int M = D.getNumRows();
        int R = D.getNumCols();

        int C = R;
        Matrix Nc = new Matrix(N);

        Matrix XN = Matrix.zeros(1, R);
        Matrix UN = Matrix.zeros(M, R);
        Matrix CN = Matrix.zeros(M, R);
        Matrix QN = Matrix.zeros(M, R);
        List<Matrix> PN = new ArrayList<Matrix>();

        Matrix prods = Matrix.zeros(1, C);
        for (int r = 0; r < C; r++) {
            double prod = 1.0;
            for (int i = 0; i < r; i++) {
                prod *= (Nc.get(i) + 1);
            }
            prods.set(0, r, prod);
        }

        double totalStates = 1.0;
        for (int r = 0; r < R; r++) {
            totalStates *= (Nc.get(r) + 1);
        }
        int numStates = (int) totalStates;

        List<Matrix> L = new ArrayList<Matrix>();
        List<Matrix> Pc = new ArrayList<Matrix>();

        for (int ist = 0; ist < M; ist++) {
            L.add(Matrix.zeros(R, numStates));

            SchedStrategy s = sched.get(ist);
            if (s == SchedStrategy.INF) {
                Pc.add(null);
                PN.add(Matrix.zeros(1, 1));
            } else if (s == SchedStrategy.PS) {
                boolean isSingleServer = S.get(ist, 0) == 1.0;
                if (isSingleServer) {
                    Pc.add(null);
                    PN.add(Matrix.zeros(1, 1));
                } else {
                    int maxServers = (int) S.get(ist, 0);
                    Pc.add(Matrix.zeros(maxServers, numStates));
                    PN.add(Matrix.zeros(maxServers, numStates));
                }
            } else if (s == SchedStrategy.FCFS) {
                boolean classIndependent = true;
                for (int r = 1; r < R; r++) {
                    if (D.get(ist, r) != D.get(ist, 0)) {
                        classIndependent = false;
                        break;
                    }
                }
                boolean isSingleServer = S.get(ist, 0) == 1.0;
                if (classIndependent) {
                    if (isSingleServer) {
                        Pc.add(null);
                        PN.add(Matrix.zeros(1, 1));
                    } else {
                        int maxServers = (int) S.get(ist, 0);
                        Pc.add(Matrix.zeros(maxServers, numStates));
                        PN.add(Matrix.zeros(maxServers, numStates));
                    }
                } else {
                    Pc.add(Matrix.zeros(numStates, numStates));
                    PN.add(Matrix.zeros(numStates, numStates));
                }
            } else {
                Pc.add(null);
                PN.add(Matrix.zeros(1, 1));
            }
        }

        Matrix x = Matrix.zeros(C, numStates);
        Matrix[] w = new Matrix[M];
        for (int i = 0; i < M; i++) w[i] = Matrix.zeros(C, numStates);

        Matrix kvec = PopulationLattice.pprod(Nc);
        for (int ist = 0; ist < M; ist++) {
            Matrix pc = Pc.get(ist);
            if (pc != null) pc.set(0, hashpop(kvec, Nc, prods), 1.0);
        }

        kvec = PopulationLattice.pprod(kvec, Nc);

        while (allGE(kvec, Matrix.zeros(1, R)) && allLE(kvec, Nc)) {
            int hkvec = hashpop(kvec, Nc, prods);

            for (int ist = 0; ist < M; ist++) {
                for (int c = 0; c < C; c++) {
                    if (kvec.get(c) == 0.0) continue;

                    double ns = S.get(ist, 0);
                    Matrix kvec_c = oner(kvec, c);
                    int hkvec_c = hashpop(kvec_c, Nc, prods);

                    double serviceRate;
                    SchedStrategy s = sched.get(ist);
                    if (s == SchedStrategy.INF) {
                        serviceRate = (D.get(ist, c) > 0) ? 1.0 / D.get(ist, c) : 0.0;
                    } else if (s == SchedStrategy.PS) {
                        if (ns == 1.0) {
                            serviceRate = (D.get(ist, c) > 0) ? 1.0 / D.get(ist, c) : 0.0;
                        } else {
                            serviceRate = computeMultiServerPSRate(ist, c, kvec, (int) ns, D, Pc.get(ist));
                        }
                    } else if (s == SchedStrategy.FCFS) {
                        if (ns == 1.0) {
                            serviceRate = (D.get(ist, c) > 0) ? 1.0 / D.get(ist, c) : 0.0;
                        } else {
                            serviceRate = computeMultiServerFCFSRate(ist, c, kvec, (int) ns, D, Pc.get(ist));
                        }
                    } else {
                        serviceRate = (D.get(ist, c) > 0) ? 1.0 / D.get(ist, c) : 0.0;
                    }

                    double wij;
                    if (serviceRate > 0) {
                        wij = (1.0 + L.get(ist).get(c, hkvec_c)) / serviceRate;
                    } else {
                        wij = GlobalConstants.Inf;
                    }
                    w[ist].set(c, hkvec, wij);
                }
            }

            for (int c = 0; c < C; c++) {
                if (kvec.get(c) == 0.0) continue;
                double totalResponseTime = 0.0;
                for (int ist = 0; ist < M; ist++) {
                    totalResponseTime += w[ist].get(c, hkvec);
                }
                double xv = (totalResponseTime > 0) ? kvec.get(c) / totalResponseTime : 0.0;
                x.set(c, hkvec, xv);
            }

            for (int ist = 0; ist < M; ist++) {
                for (int c = 0; c < C; c++) {
                    L.get(ist).set(c, hkvec, x.get(c, hkvec) * w[ist].get(c, hkvec));
                }
            }

            for (int ist = 0; ist < M; ist++) {
                Matrix pc = Pc.get(ist);
                if (pc != null) {
                    updateStateProbabilities(ist, kvec, hkvec, sched.get(ist), (int) S.get(ist, 0), pc, L.get(ist));
                }
            }

            kvec = PopulationLattice.pprod(kvec, Nc);
        }

        int finalState = hashpop(Nc, Nc, prods);

        for (int c = 0; c < C; c++) {
            XN.set(0, c, x.get(c, finalState));
        }
        for (int ist = 0; ist < M; ist++) {
            for (int c = 0; c < C; c++) {
                CN.set(ist, c, w[ist].get(c, finalState));
                QN.set(ist, c, L.get(ist).get(c, finalState));
            }
        }
        for (int ist = 0; ist < M; ist++) {
            for (int c = 0; c < C; c++) {
                UN.set(ist, c, XN.get(0, c) * D.get(ist, c));
            }
        }

        Matrix RN = Matrix.zeros(M, R);
        for (int ist = 0; ist < M; ist++) {
            for (int c = 0; c < C; c++) {
                RN.set(ist, c, XN.get(0, c) > 0 ? QN.get(ist, c) / XN.get(0, c) : 0.0);
            }
        }

        Matrix TN = Matrix.zeros(M, R);
        for (int ist = 0; ist < M; ist++) {
            for (int c = 0; c < C; c++) {
                TN.set(ist, c, XN.get(0, c));
            }
        }

        for (int ist = 0; ist < M; ist++) {
            Matrix pc = Pc.get(ist);
            if (pc != null) {
                PN.set(ist, new Matrix(pc));
            }
        }

        double runtime = (System.nanoTime() - startTime) / 1e9;
        return new Ret.pfqnSchmidt(QN, UN, RN, TN, CN, XN, PN, "schmidt", numStates, runtime);
    }

    private static int hashpop(Matrix kvec, Matrix Nc, Matrix prods) {
        double hash = 0.0;
        for (int r = 0; r < kvec.length(); r++) {
            hash += kvec.get(r) * prods.get(0, r);
        }
        return (int) hash;
    }

    private static Matrix oner(Matrix kvec, int c) {
        Matrix result = new Matrix(kvec);
        if (result.get(c) > 0) {
            result.set(c, result.get(c) - 1);
        }
        return result;
    }

    private static boolean allGE(Matrix a, Matrix b) {
        for (int i = 0; i < a.length(); i++) {
            if (a.get(i) < b.get(i)) return false;
        }
        return true;
    }

    private static boolean allLE(Matrix a, Matrix b) {
        for (int i = 0; i < a.length(); i++) {
            if (a.get(i) > b.get(i)) return false;
        }
        return true;
    }

    private static double computeMultiServerPSRate(int stationIdx, int classIdx, Matrix kvec, int numServers,
                                                   Matrix D, Matrix Pc) {
        if (D.get(stationIdx, classIdx) == 0.0) return 0.0;
        int totalJobs = (int) kvec.elementSum();
        int activeServers = Math.min(totalJobs, numServers);
        if (activeServers > 0) {
            return (double) activeServers / D.get(stationIdx, classIdx);
        } else {
            return 1.0 / D.get(stationIdx, classIdx);
        }
    }

    private static double computeMultiServerFCFSRate(int stationIdx, int classIdx, Matrix kvec, int numServers,
                                                     Matrix D, Matrix Pc) {
        if (D.get(stationIdx, classIdx) == 0.0) return 0.0;
        int totalJobs = (int) kvec.elementSum();
        int busyServers = Math.min(totalJobs, numServers);
        if (busyServers > 0) {
            return (double) busyServers / D.get(stationIdx, classIdx);
        } else {
            return 1.0 / D.get(stationIdx, classIdx);
        }
    }

    private static void updateStateProbabilities(int stationIdx, Matrix kvec, int hkvec,
                                                 SchedStrategy sched, int numServers, Matrix Pc, Matrix L) {
        if (sched == SchedStrategy.PS || sched == SchedStrategy.FCFS) {
            if (numServers > 1) {
                int totalJobs = (int) kvec.elementSum();
                for (int j = 0; j < Math.min(totalJobs + 1, Pc.getNumRows()); j++) {
                    double meanJobs = L.sumRows().get(hkvec);
                    if (meanJobs > 0) {
                        Pc.set(j, hkvec, Math.exp(-meanJobs) * FastMath.pow(meanJobs, (double) j) / factorial(j));
                    }
                }
            }
        }
    }

    private static double factorial(int n) {
        if (n <= 1) return 1.0;
        double result = 1.0;
        for (int i = 2; i <= n; i++) result *= i;
        return result;
    }
}
