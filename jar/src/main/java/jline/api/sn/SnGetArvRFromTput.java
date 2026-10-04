/**
 * @file Compute arrival rates from network throughputs
 *
 * @since LINE 3.0
 */
package jline.api.sn;

import java.util.ArrayList;
import java.util.List;

import jline.lang.NetworkStruct;
import jline.lang.constant.NodeType;
import jline.lang.nodeparam.CacheNodeParam;
import jline.solvers.AvgHandle;
import jline.util.matrix.Matrix;

public final class SnGetArvRFromTput {
    private SnGetArvRFromTput() {}

    /**
     * Calculates the average arrival rates at each station from the network throughputs.
     */
    public static Matrix snGetArvRFromTput(NetworkStruct sn, Matrix TN, AvgHandle TH) {
        int M = sn.nstations;
        int R = sn.nclasses;
        Matrix AN = new Matrix(M, R);

        if (TH != null && TN != null) {
            List<Integer> statefulNodes = new ArrayList<Integer>();
            for (int i = 0; i < sn.isstateful.length(); i++) {
                if (sn.isstateful.get(i) == 1.0) {
                    statefulNodes.add(i);
                }
            }
            int nStateful = statefulNodes.size();

            Matrix TN_stateful = new Matrix(nStateful, R);

            for (int sf = 0; sf < nStateful; sf++) {
                int ind = statefulNodes.get(sf);
                int ist = (int) sn.nodeToStation.get(ind);
                if (ist >= 0) {
                    for (int r = 0; r < R; r++) {
                        TN_stateful.set(sf, r, TN.get(ist, r));
                    }
                }
            }

            for (int sf = 0; sf < nStateful; sf++) {
                int ind = statefulNodes.get(sf);
                int ist = (int) sn.nodeToStation.get(ind);
                if (ist < 0) {
                    if (sn.nodetype.get(ind) == NodeType.Cache) {
                        Object statefulNode = sn.stateful.get(sf);
                        Object np = sn.nodeparam.get(statefulNode);
                        CacheNodeParam cacheParam = (np instanceof CacheNodeParam) ? (CacheNodeParam) np : null;
                        if (cacheParam != null) {
                            Matrix hitclass = cacheParam.hitclass;
                            Matrix missclass = cacheParam.missclass;
                            Matrix actualHitProb = cacheParam.actualhitprob;
                            Matrix actualMissProb = cacheParam.actualmissprob;

                            if (actualHitProb != null && actualMissProb != null && !actualHitProb.isEmpty()) {
                                for (int c = 0; c < sn.nchains; c++) {
                                    Matrix inchain = sn.inchain.get(c);
                                    if (inchain == null) continue;
                                    int refstat = (int) sn.refstat.get(c);

                                    double totalTput = 0.0;
                                    for (int idx = 0; idx < inchain.length(); idx++) {
                                        int classIdx = (int) inchain.get(idx);
                                        if (classIdx >= 0 && classIdx < R) {
                                            totalTput += TN.get(refstat, classIdx);
                                        }
                                    }

                                    for (int k = 0; k < R; k++) {
                                        if (hitclass.length() > k && missclass.length() > k) {
                                            int hc = (int) hitclass.get(k);
                                            int mc = (int) missclass.get(k);
                                            if (hc >= 0 && hc < R && actualHitProb.length() > k) {
                                                double hitProbVal = actualHitProb.get(k);
                                                if (!Double.isNaN(hitProbVal)) {
                                                    TN_stateful.set(sf, hc, totalTput * hitProbVal);
                                                }
                                            }
                                            if (mc >= 0 && mc < R && actualMissProb.length() > k) {
                                                double missProbVal = actualMissProb.get(k);
                                                if (!Double.isNaN(missProbVal)) {
                                                    TN_stateful.set(sf, mc, totalTput * missProbVal);
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }

            // Non-station, non-Cache stateful nodes (e.g., Router) have no throughput of their own. Solve the flow balance
            // x_u = x_k*P_ku + x_u*P_uu exactly: a sweep in node order left a Router fed by a later-declared Router at zero.
            List<Integer> unk = new ArrayList<Integer>();
            boolean[] isUnkCol = new boolean[nStateful * R];
            for (int sf = 0; sf < nStateful; sf++) {
                int ind = statefulNodes.get(sf);
                int ist = (int) sn.nodeToStation.get(ind);
                if (ist < 0 && sn.nodetype.get(ind) != NodeType.Cache) {
                    unk.add(sf);
                    for (int k = 0; k < R; k++) {
                        isUnkCol[sf * R + k] = true;
                    }
                }
            }
            if (!unk.isEmpty()) {
                int nU = unk.size() * R;
                int[] ucols = new int[nU];
                for (int u = 0; u < unk.size(); u++) {
                    for (int k = 0; k < R; k++) {
                        ucols[u * R + k] = unk.get(u) * R + k;
                    }
                }
                // b = x_k * P_ku and A = I - P_uu, solved as A^T x_u^T = b^T
                Matrix b = new Matrix(nU, 1);
                Matrix At = Matrix.eye(nU);
                for (int j = 0; j < nU; j++) {
                    double acc = 0.0;
                    for (int col = 0; col < nStateful * R; col++) {
                        if (!isUnkCol[col]) {
                            acc += TN_stateful.get(col / R, col % R) * sn.rt.get(col, ucols[j]);
                        }
                    }
                    b.set(j, 0, acc);
                    for (int i = 0; i < nU; i++) {
                        At.set(j, i, At.get(j, i) - sn.rt.get(ucols[i], ucols[j]));
                    }
                }
                Matrix xu = At.leftMatrixDivide(b);
                for (int j = 0; j < nU; j++) {
                    TN_stateful.set(ucols[j] / R, ucols[j] % R, xu.get(j, 0));
                }
            }

            for (int ist = 0; ist < M; ist++) {
                int ind_ist = (int) sn.stationToNode.get(ist);
                if (sn.nodetype.get(ind_ist) == NodeType.Source) continue;

                int sf_ist = statefulNodes.indexOf(ind_ist);
                if (sf_ist < 0) continue;

                for (int k = 0; k < R; k++) {
                    double a = 0.0;
                    for (int sf_jst = 0; sf_jst < nStateful; sf_jst++) {
                        for (int r = 0; r < R; r++) {
                            a += TN_stateful.get(sf_jst, r) * sn.rt.get(sf_jst * R + r, sf_ist * R + k);
                        }
                    }
                    AN.set(ist, k, a);
                }
            }
        } else {
            AN = new Matrix(0, 0);
        }

        if (sn.fj.any()) {
            Matrix ANn = SnGetNodeArvRFromTput.snGetNodeArvRFromTput(sn, TN, TH, AN);
            for (int ist = 0; ist < M; ist++) {
                for (int r = 0; r < R; r++) {
                    AN.set(ist, r, ANn.get((int) sn.stationToNode.get(ist), r));
                }
            }
        }

        // see _kb/03-api-layer.md for rationale
        boolean hasPlace = false;
        for (int ind = 0; ind < sn.nnodes; ind++) {
            if (sn.nodetype.get(ind) == NodeType.Place) {
                hasPlace = true;
                break;
            }
        }
        if (!AN.isEmpty() && hasPlace) {
            SnPnFiringRates.Ret pnRet = SnPnFiringRates.snPnFiringRates(sn, TN, true);
            if (pnRet.rates != null) {
                for (int pp = 0; pp < pnRet.placeNodes.size(); pp++) {
                    int ist = (int) sn.nodeToStation.get(pnRet.placeNodes.get(pp));
                    if (ist < 0) {
                        continue;
                    }
                    for (int k = 0; k < R; k++) {
                        double arvr = 0.0;
                        for (int mm = 0; mm < pnRet.rates.getNumRows(); mm++) {
                            arvr += pnRet.produced[mm][pp][k] * pnRet.rates.get(mm, 0);
                        }
                        AN.set(ist, k, arvr);
                    }
                }
            }
        }

        return AN;
    }
}
