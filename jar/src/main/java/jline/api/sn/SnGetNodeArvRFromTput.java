package jline.api.sn;

import jline.lang.NetworkStruct;
import jline.lang.constant.NodeType;
import jline.lang.nodeparam.CacheNodeParam;
import jline.solvers.AvgHandle;
import jline.util.matrix.Matrix;

public final class SnGetNodeArvRFromTput {
    private SnGetNodeArvRFromTput() {}

    public static Matrix snGetNodeArvRFromTput(NetworkStruct sn, Matrix TN, AvgHandle TH, Matrix AN) {
        int I = sn.nnodes;
        int C = sn.nchains;
        int M = sn.nstations;
        int R = sn.nclasses;
        Matrix ANn = new Matrix(I, R);

        // First, copy station arrival rates to station nodes
        if (AN != null && AN.getNumRows() > 0) {
            for (int ist = 0; ist < M; ist++) {
                int ind = (int) sn.stationToNode.get(ist);
                if (ind >= 0 && ind < I) {
                    for (int r = 0; r < R; r++) {
                        ANn.set(ind, r, AN.get(ist, r));
                    }
                }
            }
        }

        for (int ind = 0; ind < I; ind++) {
            int nodeToStation = (int) sn.nodeToStation.get(ind);
            if (nodeToStation >= 0) continue;

            for (int c = 0; c < C; c++) {
                Matrix inchain = sn.inchain.get(c);
                int refstat = (int) sn.refstat.get(c);
                for (int r = 0; r < inchain.length(); r++) {
                    int rIdx = (int) inchain.get(r);
                    if (sn.nodetype.get(ind) != NodeType.Source) {
                        if (sn.nodetype.get(ind) == NodeType.Cache) {
                            CacheNodeParam np = (CacheNodeParam) sn.nodeparam.get(sn.nodes.get(ind));
                            boolean isHitClass = !np.hitclass.findNumber((double) rIdx).isEmpty();
                            boolean isMissClass = !np.missclass.findNumber((double) rIdx).isEmpty();
                            if (!(isHitClass || isMissClass)) {
                                double num = sn.nodevisits.get(c).get(ind, rIdx);
                                double den = 0.0;
                                for (int r2 = 0; r2 < inchain.length(); r2++) {
                                    int rprime = (int) inchain.get(r2);
                                    den += sn.visits.get(c).get((int) sn.stationToStateful.get(refstat), rprime);
                                }
                                double coeff = 0.0;
                                for (int r2 = 0; r2 < inchain.length(); r2++) {
                                    int rprime = (int) inchain.get(r2);
                                    coeff += TN.get(refstat, rprime);
                                }
                                if (den > 0) {
                                    ANn.set(ind, rIdx, num / den * coeff);
                                } else {
                                    ANn.set(ind, rIdx, 0);
                                }
                            }
                        } else {
                            double num = sn.nodevisits.get(c).get(ind, rIdx);
                            double den = 0.0;
                            for (int r2 = 0; r2 < inchain.length(); r2++) {
                                int rprime = (int) inchain.get(r2);
                                den += sn.visits.get(c).get((int) sn.stationToStateful.get(refstat), rprime);
                            }
                            double coeff = 0.0;
                            for (int r2 = 0; r2 < inchain.length(); r2++) {
                                int rprime = (int) inchain.get(r2);
                                coeff += TN.get(refstat, rprime);
                            }
                            if (den > 0) {
                                ANn.set(ind, rIdx, num / den * coeff);
                            } else {
                                ANn.set(ind, rIdx, 0);
                            }
                        }
                    }
                }
            }
        }
        return ANn;
    }
}
