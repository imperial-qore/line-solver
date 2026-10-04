package jline.api.sn;

import jline.GlobalConstants;
import jline.io.Ret;
import jline.lang.NetworkStruct;
import jline.lang.nodes.Node;
import jline.lang.nodes.Sink;
import jline.util.Utils;
import jline.util.matrix.Matrix;

public final class SnDeaggregateChainResults {
    private SnDeaggregateChainResults() {}

    /**
     * Calculate class-based performance metrics for a queueing network based on performance measures of its chains.
     */
    public static Ret.snDeaggregateChainResults snDeaggregateChainResults(NetworkStruct sn,
                                                                          Matrix Lchain,
                                                                          Matrix ST,
                                                                          Matrix STchain,
                                                                          Matrix Vchain,
                                                                          Matrix alpha,
                                                                          Matrix Qchain,
                                                                          Matrix Uchain,
                                                                          Matrix Rchain,
                                                                          Matrix Tchain,
                                                                          Matrix Cchain,
                                                                          Matrix Xchain) {
        Matrix STlocal = ST;
        if (STlocal == null || STlocal.isEmpty()) {
            STlocal = new Matrix(0, 0);
            sn.rates.divide(1.0, STlocal, false);
            STlocal.removeNaN();
        }

        if (Cchain != null && !Cchain.isEmpty()) {
            throw new RuntimeException("Cchain input to snDeaggregateChainResults not yet supported");
        }

        int M = sn.nstations;
        int K = sn.nclasses;
        Matrix X = new Matrix(1, K);
        Matrix U = new Matrix(M, K);
        Matrix Q = new Matrix(M, K);
        Matrix T = new Matrix(M, K);
        Matrix R = new Matrix(M, K);
        Matrix C = new Matrix(1, K);

        int idxSink = 0;
        for (Node nodeIter : sn.nodes) {
            if (nodeIter instanceof Sink) idxSink = nodeIter.getNodeIndex();
        }

        Matrix Vsinktmp = new Matrix(sn.nodevisits.get(0).getNumRows(), sn.nodevisits.get(0).getNumCols());
        for (int i = 0; i < sn.nodevisits.size(); i++) {
            Vsinktmp = Vsinktmp.add(1.0, sn.nodevisits.get(i));
        }
        Matrix Vsink = Matrix.extractRows(Vsinktmp, idxSink, idxSink + 1, null);
        for (int c = 0; c < sn.nchains; c++) {
            Matrix inchain_c = sn.inchain.get(c);
            double sum = 0.0;
            for (int idx = 0; idx < inchain_c.getNumCols(); idx++) {
                sum += sn.njobs.get((int) inchain_c.get(idx));
            }
            for (int idx = 0; idx < inchain_c.getNumCols(); idx++) {
                int k = (int) inchain_c.get(0, idx);
                if (Utils.isInf(sum)) {
                    X.set(0, k, Xchain.get(0, c) * Vsink.get(0, k));
                } else {
                    X.set(0, k, Xchain.get(0, c) * alpha.get((int) sn.refstat.get(k, 0), k));
                }
                for (int i = 0; i < M; i++) {
                    if (Uchain == null || Uchain.isEmpty()) {
                        if (Utils.isInf(sn.nservers.get(i, 0))) {
                            U.set(i, k, STlocal.get(i, k) * (Xchain.get(0, c) * Vchain.get(i, c) / Vchain.get((int) sn.refstat.get(k, 0), c)) * alpha.get(i, k));
                        } else {
                            U.set(i, k, STlocal.get(i, k) * (Xchain.get(0, c) * Vchain.get(i, c) / Vchain.get((int) sn.refstat.get(k, 0), c)) * alpha.get(i, k) / sn.nservers.get(i, 0));
                        }
                    } else {
                        if (Utils.isInf(sn.nservers.get(i, 0))) {
                            U.set(i, k, STlocal.get(i, k) * (Xchain.get(0, c) * Vchain.get(i, c) / Vchain.get((int) sn.refstat.get(k, 0), c)) * alpha.get(i, k));
                        } else {
                            U.set(i, k, Uchain.get(i, c) * alpha.get(i, k));
                        }
                    }

                    if (Lchain.get(i, c) > 0) {
                        if (Qchain != null && !Qchain.isEmpty()) {
                            Q.set(i, k, Qchain.get(i, c) * alpha.get(i, k));
                        } else {
                            Q.set(i, k, Rchain.get(i, c) * STlocal.get(i, k) / STchain.get(i, c) * Xchain.get(0, c) * Vchain.get(i, c) / Vchain.get((int) sn.refstat.get(k, 0), c) * alpha.get(i, k));
                        }
                        T.set(i, k, Tchain.get(i, c) * alpha.get(i, k));
                        R.set(i, k, Q.get(i, k) / T.get(i, k));
                    } else {
                        T.set(i, k, 0.0);
                        R.set(i, k, 0.0);
                        Q.set(i, k, 0.0);
                    }
                }
                C.set(0, k, sn.njobs.get(0, k) / X.get(0, k));
            }
        }

        Q.absEq();
        R.absEq();
        X.absEq();
        U.absEq();
        T.absEq();
        C.absEq();
        Q.removeNaN();
        Q.apply(GlobalConstants.Inf, 0.0, "equal");
        R.removeNaN();
        R.apply(GlobalConstants.Inf, 0.0, "equal");
        X.removeNaN();
        X.apply(GlobalConstants.Inf, 0.0, "equal");
        U.removeNaN();
        U.apply(GlobalConstants.Inf, 0.0, "equal");
        T.removeNaN();
        T.apply(GlobalConstants.Inf, 0.0, "equal");
        C.removeNaN();
        C.apply(GlobalConstants.Inf, 0.0, "equal");

        return new Ret.snDeaggregateChainResults(Q, U, R, T, C, X);
    }
}
