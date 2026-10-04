package jline.api.sn;

import jline.lang.NetworkStruct;
import jline.lang.constant.NodeType;
import jline.lang.nodeparam.CacheNodeParam;
import jline.solvers.AvgHandle;
import jline.util.matrix.Matrix;

public final class SnGetNodeTputFromTput {
    private SnGetNodeTputFromTput() {}

    public static Matrix snGetNodeTputFromTput(NetworkStruct sn, Matrix TN, AvgHandle TH, Matrix ANn) {
        int I = sn.nnodes;
        int C = sn.nchains;
        int M = sn.nstations;
        int R = sn.nclasses;
        Matrix TNn = new Matrix(I, R);

        // First, copy station throughput to station nodes
        for (int ist = 0; ist < M; ist++) {
            int ind = (int) sn.stationToNode.get(ist);
            if (ind >= 0 && ind < I) {
                for (int r = 0; r < R; r++) {
                    TNn.set(ind, r, TN.get(ist, r));
                }
            }
        }

        // Process Cache hit/miss classes using nodevisits formula
        for (int ind = 0; ind < I; ind++) {
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
                            if (isHitClass || isMissClass) {
                                double totalTput = 0.0;
                                for (int s = 0; s < inchain.length(); s++) {
                                    int sIdx = (int) inchain.get(s);
                                    totalTput = totalTput + TN.get(refstat, sIdx);
                                }
                                boolean haveActual = np.actualhitprob != null && !np.actualhitprob.isEmpty()
                                        && np.actualmissprob != null && !np.actualmissprob.isEmpty();
                                if (haveActual) {
                                    // see _kb/03-api-layer.md for rationale
                                    double acc = 0.0;
                                    boolean touched = false;
                                    for (int origClass = 0; origClass < np.hitclass.length(); origClass++) {
                                        boolean inChain = !inchain.findNumber((double) origClass).isEmpty();
                                        double arvTput = inChain ? TN.get(refstat, origClass) : 0.0;
                                        if (np.hitclass.get(origClass) == rIdx && !Double.isNaN(np.actualhitprob.get(origClass))) {
                                            // Hit-class throughput is (true hit + delayed hit) *
                                            // class arrival; delayed is zero for plain caches.
                                            double dh = 0.0;
                                            if (np.actualdelayedhitprob != null
                                                    && origClass < np.actualdelayedhitprob.length()
                                                    && !Double.isNaN(np.actualdelayedhitprob.get(origClass))) {
                                                dh = np.actualdelayedhitprob.get(origClass);
                                            }
                                            acc += arvTput * (np.actualhitprob.get(origClass) + dh);
                                            touched = true;
                                        } else if (np.missclass.get(origClass) == rIdx && !Double.isNaN(np.actualmissprob.get(origClass))) {
                                            acc += arvTput * np.actualmissprob.get(origClass);
                                            touched = true;
                                        }
                                    }
                                    if (touched) {
                                        TNn.set(ind, rIdx, acc);
                                    }
                                } else {
                                    // Fallback to nodevisits-based calculation
                                    double num = sn.nodevisits.get(c).get(ind, rIdx);
                                    double den = 0.0;
                                    for (int s = 0; s < inchain.length(); s++) {
                                        int sIdx = (int) inchain.get(s);
                                        den = den + sn.visits.get(c).get((int) sn.stationToStateful.get(refstat), sIdx);
                                    }
                                    if (den > 0) {
                                        TNn.set(ind, rIdx, num / den * totalTput);
                                    } else {
                                        TNn.set(ind, rIdx, 0);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        // Process non-station nodes using routing formula
        for (int ind = 0; ind < I; ind++) {
            int nodeToStation = (int) sn.nodeToStation.get(ind);
            if (nodeToStation >= 0) continue;

            for (int c = 0; c < C; c++) {
                Matrix inchain = sn.inchain.get(c);
                for (int r = 0; r < inchain.length(); r++) {
                    int rIdx = (int) inchain.get(r);
                    boolean anystateful = !sn.visits.get(c).getColumn(rIdx).isEmpty();
                    if (anystateful) {
                        if (sn.nodetype.get(ind) != NodeType.Sink && sn.nodetype.get(ind) != NodeType.Join) {
                            for (int s = 0; s < inchain.length(); s++) {
                                int sIdx = (int) inchain.get(s);
                                for (int jnd = 0; jnd < I; jnd++) {
                                    if (sn.nodetype.get(ind) == NodeType.Cache) {
                                        if (ind != jnd) {
                                            TNn.set(ind, sIdx, TNn.get(ind, sIdx) + ANn.get(ind, rIdx) * sn.rtnodes.get(ind * R + rIdx, jnd * R + sIdx));
                                        }
                                    } else {
                                        TNn.set(ind, sIdx, TNn.get(ind, sIdx) + ANn.get(ind, rIdx) * sn.rtnodes.get(ind * R + rIdx, jnd * R + sIdx));
                                    }
                                }
                            }
                        } else if (sn.nodetype.get(ind) == NodeType.Join) {
                            for (int s = 0; s < inchain.length(); s++) {
                                int sIdx = (int) inchain.get(s);
                                for (int jnd = 0; jnd < I; jnd++) {
                                    TNn.set(ind, sIdx, TNn.get(ind, sIdx) + ANn.get(ind, rIdx) * sn.rtnodes.get((ind) * R + rIdx, (jnd) * R + sIdx));
                                }
                            }
                        }
                    }
                }
            }
        }
        return TNn;
    }
}
