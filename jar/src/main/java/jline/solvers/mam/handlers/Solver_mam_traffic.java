package jline.solvers.mam.handlers;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import jline.GlobalConstants;
import jline.api.mam.Mmap_lambda;
import jline.api.mam.Mmap_normalize;
import jline.api.mam.Mmap_super;
import jline.api.mc.Dtmc_stochcomp;
import jline.api.npfqn.Npfqn_traffic_merge;
import jline.api.npfqn.Npfqn_traffic_split_cs;
import jline.lang.NetworkStruct;
import jline.lang.constant.NodeType;
import jline.solvers.SolverOptions;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Solver_mam_traffic {
    private Solver_mam_traffic() {}

    public static Map<Integer, MatrixCell> solver_mam_traffic(NetworkStruct sn,
                                                              Map<Integer, Map<Integer, MatrixCell>> DEP,
                                                              SolverOptions.Config config) {
        int I = sn.nnodes;
        int R = sn.nclasses;
        Matrix non_cs_classes = new Matrix(1, I * R, I * R);
        Matrix isNCS = new Matrix(1, I, I);
        Matrix nodeToNCS = new Matrix(1, I, I);
        int end = 0;
        for (int ind = 0; ind < I; ind++) {
            if (sn.nodetype.get(ind) != NodeType.ClassSwitch) {
                for (int i = 0; i < R; i++) {
                    non_cs_classes.set(i + end, (double) (ind * R + i));
                }
                end = end + R;
                isNCS.set(ind, 1.0);
                nodeToNCS.set(ind, isNCS.elementSum());
            } else {
                isNCS.set(ind, 0.0);
            }
        }

        List<Integer> non_cs_classes_list = new ArrayList<Integer>();
        for (int i = 0; i < non_cs_classes.length(); i++) {
            non_cs_classes_list.add((int) non_cs_classes.get(i));
        }

        Matrix rtncs = Dtmc_stochcomp.dtmc_stochcomp(sn.rtnodes, non_cs_classes_list);
        int Inc = I;
        for (int i = 0; i < sn.nodetype.size(); i++) {
            if (sn.nodetype.get(i) == NodeType.ClassSwitch) {
                Inc = Inc - 1;
            }
        }

        for (int ist = 0; ist < DEP.size(); ist++) {
            Map<Integer, MatrixCell> entry = DEP.get(ist);
            for (int r = 0; r < entry.size(); r++) {
                MatrixCell cell = entry.get(r);
                if (cell.isEmpty() || cell.get(0).hasNaN()) {
                    cell.set(0, new Matrix(1, 1, 0));
                    cell.set(1, new Matrix(1, 1, 0));
                    cell.set(2, new Matrix(1, 1, 0));
                } else if (cell.size() <= 2) {
                    // see _kb/06-solver-catalog.md for rationale
                    cell.set(2, cell.get(1));
                }
            }
        }

        Map<Integer, Map<Integer, MatrixCell>> DEP_new = new HashMap<Integer, Map<Integer, MatrixCell>>();
        Map<Integer, Map<Integer, MatrixCell>> LINKS = new HashMap<Integer, Map<Integer, MatrixCell>>();
        Map<Integer, MatrixCell> ARV = new HashMap<Integer, MatrixCell>();
        for (int ind = 0; ind < I; ind++) {
            if (isNCS.get(ind) == 1.0) {
                int inc = (int) nodeToNCS.get(ind);
                if (sn.nodetype.get(ind) == NodeType.Source
                        || sn.nodetype.get(ind) == NodeType.Delay
                        || sn.nodetype.get(ind) == NodeType.Queue) {
                    int ist = (int) sn.nodeToStation.get(ind);
                    DEP_new.put(inc, new HashMap<Integer, MatrixCell>());
                    // see _kb/06-solver-catalog.md for rationale
                    boolean markedSource = false;
                    if (sn.markidx != null && ist < sn.markidx.getNumRows()) {
                        for (int r = 0; r < R; r++) {
                            if (sn.markidx.get(ist, r) > 0) {
                                markedSource = true;
                                break;
                            }
                        }
                    }
                    if (markedSource) {
                        MatrixCell shared = null;
                        for (int r = 0; r < R; r++) {
                            if (sn.markidx.get(ist, r) > 0) {
                                shared = DEP.get(ist).get(r);
                                break;
                            }
                        }
                        int nph = shared.get(0).getNumRows();
                        MatrixCell expanded = new MatrixCell();
                        expanded.set(0, shared.get(0).copy());
                        expanded.set(1, shared.get(1).copy());
                        for (int r = 0; r < R; r++) {
                            int mk = (int) sn.markidx.get(ist, r);
                            if (mk > 0) {
                                expanded.set(2 + r, shared.get(1 + mk).copy());
                            } else {
                                // see _kb/06-solver-catalog.md for rationale
                                MatrixCell depr = DEP.get(ist).get(r);
                                if (depr != null && !depr.isEmpty() && depr.size() > 1
                                        && depr.get(1).elementSum() > GlobalConstants.FineTol
                                        && depr != shared) {
                                    throw new RuntimeException("SolverMAM: mixing a marked "
                                            + "(setMarkedArrival) source with independent per-class "
                                            + "arrivals at the same Source is not supported yet.");
                                }
                                expanded.set(2 + r, new Matrix(nph, nph));
                            }
                        }
                        DEP_new.get(inc).put(0, expanded);
                    } else if (R > 1) {
                        MatrixCell superposedMMAP = DEP.get(ist).get(0);
                        for (int r = 1; r < R; r++) {
                            superposedMMAP = Mmap_super.mmap_super(superposedMMAP, DEP.get(ist).get(r));
                        }
                        DEP_new.get(inc).put(0, superposedMMAP);
                    } else {
                        DEP_new.get(inc).put(0, DEP.get(ist).get(0));
                    }
                    Matrix Psplit = new Matrix(R, Inc * R, R * Inc * R);
                    for (int r = 0; r < R; r++) {
                        for (int jnd = 0; jnd < I; jnd++) {
                            if (isNCS.get(jnd) == 1.0) {
                                int jnc = (int) nodeToNCS.get(jnd);
                                for (int s = 0; s < R; s++) {
                                    Psplit.set(r, (jnc - 1) * R + s,
                                            rtncs.get((inc - 1) * R + r, (jnc - 1) * R + s));
                                }
                            }
                        }
                    }
                    @SuppressWarnings("unchecked")
                    Map<Integer, MatrixCell> Fsplit = (Map<Integer, MatrixCell>)
                            (Map<?, ?>) Npfqn_traffic_split_cs.npfqn_traffic_split_cs(
                                    DEP_new.get(inc).get(0), Psplit);
                    // see _kb/06-solver-catalog.md for rationale
                    LINKS.put(inc, new HashMap<Integer, MatrixCell>());
                    for (int jnc = 0; jnc < Inc; jnc++) {
                        LINKS.get(inc).put(jnc, Mmap_normalize.mmap_normalize(Fsplit.get(jnc)));
                    }
                }
            }
        }

        for (int ind = 0; ind < I; ind++) {
            Map<Integer, MatrixCell> FLOWS = new HashMap<Integer, MatrixCell>();
            if (isNCS.get(ind) == 1.0 && sn.nodetype.get(ind) != NodeType.Source) {
                int inc = (int) nodeToNCS.get(ind);
                // see _kb/06-solver-catalog.md for rationale
                for (int jnd = 1; jnd <= Inc; jnd++) {
                    Map<Integer, MatrixCell> outLinks = LINKS.get(jnd);
                    MatrixCell flow = (outLinks != null) ? outLinks.get(inc - 1) : null;
                    if (flow != null && !flow.isEmpty()
                            && Mmap_lambda.mmap_lambda(flow).elementSum() > GlobalConstants.FineTol) {
                        FLOWS.put(FLOWS.size(), flow);
                    }
                }
                if (FLOWS.size() > 1) {
                    ARV.put(ind, Npfqn_traffic_merge.npfqn_traffic_merge(FLOWS, config.merge, config.compress));
                } else if (FLOWS.size() == 1) {
                    ARV.put(ind, FLOWS.get(0));
                } else {
                    // all links are zeros, take one (MATLAB: LINKS{jnd,1} with
                    // the loop-final jnd)
                    MatrixCell fallback = null;
                    for (int j = 1; j <= Inc; j++) {
                        if (LINKS.get(j) != null && LINKS.get(j).get(0) != null) {
                            fallback = LINKS.get(j).get(0);
                        }
                    }
                    ARV.put(ind, fallback != null ? fallback : new MatrixCell());
                }
            } else {
                ARV.put(ind, new MatrixCell());
            }
        }

        return ARV;
    }
}
