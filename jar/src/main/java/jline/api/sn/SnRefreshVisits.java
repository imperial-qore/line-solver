/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.sn;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import jline.GlobalConstants;
import jline.api.mc.Dtmc_solve;
import jline.api.mc.Dtmc_solve_reducible;
import jline.lang.NetworkStruct;
import jline.lang.constant.NodeType;
import jline.util.Pair;
import jline.util.matrix.Matrix;

/**
 * Stochastic Network Visit Ratio Calculator.
 */
public final class SnRefreshVisits {
    private SnRefreshVisits() {}

    public static NetworkStruct snRefreshVisits(NetworkStruct sn, Matrix chains, Matrix rt, Matrix rtnodes) {
        int I = sn.nnodes;
        int M = sn.nstateful;
        int K = sn.nclasses;
        Matrix refstat = sn.refstat.copy();
        int nchains = sn.nchains;

        Map<Integer, Matrix> inchain = sn.inchain;
        for (int c = 0; c < nchains; c++) {
            Matrix inchain_c = inchain.get(c);
            double refstatValue = refstat.get((int) inchain_c.get(0, 0), 0);
            for (int col = 1; col < inchain_c.getNumCols(); col++) {
                int row = (int) inchain_c.get(0, col);
                if (refstatValue != refstat.get(row, 0)) refstat.set(row, 0, refstatValue);
            }
        }

        Map<Integer, List<Integer>> new_inchain = new HashMap<Integer, List<Integer>>();
        for (int c = 0; c < nchains; c++) {
            Matrix inchain_c = inchain.get(c);
            List<Integer> list = new ArrayList<Integer>();
            for (int i = 0; i < inchain_c.getNumCols(); i++) list.add((int) inchain_c.get(i));
            new_inchain.put(c, list);
        }

        boolean hasFork = false;
        for (int n = 0; n < sn.nodetype.size(); n++) {
            if (sn.nodetype.get(n) == NodeType.Fork) { hasFork = true; break; }
        }

        Map<Integer, Matrix> visits = new HashMap<Integer, Matrix>();
        for (int c = 0; c < nchains; c++) {
            List<Integer> inchain_c = new_inchain.get(c);
            List<Integer> cols = new ArrayList<Integer>();
            for (int ist = 0; ist < M; ist++) {
                for (int ik = 0; ik < inchain_c.size(); ik++) {
                    cols.add(ist * K + inchain_c.get(ik));
                }
            }

            Matrix Pchain = new Matrix(cols.size(), cols.size());
            for (int row = 0; row < cols.size(); row++) {
                for (int col = 0; col < cols.size(); col++) {
                    Pchain.set(row, col, rt.get(cols.get(row), cols.get(col)));
                }
            }

            replaceNaNWithEqualProb(Pchain);

            double[] rowSums = new double[Pchain.getNumRows()];
            for (int i = 0; i < rowSums.length; i++) rowSums[i] = 1.0;
            if (hasFork) {
                for (int row = 0; row < Pchain.getNumRows(); row++) {
                    double rs = Pchain.sumRows(row);
                    rowSums[row] = rs;
                    if (rs > 1e-8) {
                        for (int col = 0; col < Pchain.getNumCols(); col++) {
                            Pchain.set(row, col, Pchain.get(row, col) / rs);
                        }
                    }
                }
            }

            Matrix visited = new Matrix(Pchain.getNumRows(), 1);
            int countTrue = 0;
            for (int row = 0; row < Pchain.getNumRows(); row++) {
                if (Pchain.sumRows(row) > 0) {
                    countTrue++;
                    visited.set(row, 0, 1.0);
                }
            }

            Matrix input = new Matrix(countTrue, countTrue);
            int row_input = 0;
            for (int row = 0; row < visited.getNumRows(); row++) {
                if (visited.get(row, 0) > 0) {
                    int col_input = 0;
                    for (int col = 0; col < visited.getNumRows(); col++) {
                        if (visited.get(col, 0) > 0) {
                            input.set(row_input, col_input, Pchain.get(row, col));
                            col_input++;
                        }
                    }
                    row_input++;
                }
            }

            Matrix alpha_visited;
            try {
                alpha_visited = Dtmc_solve.dtmc_solve(input);
                if (alpha_visited.elementMax() == 0.0 || alpha_visited.hasNaN()) {
                    alpha_visited = Dtmc_solve_reducible.dtmc_solve_reducible(input).getLeft();
                }
            } catch (Exception e) {
                alpha_visited = Dtmc_solve_reducible.dtmc_solve_reducible(input).getLeft();
            }

            Matrix alpha = new Matrix(1, M * K);
            int idx = 0;
            for (int row = 0; row < visited.getNumRows(); row++) {
                if (visited.get(row, 0) > 0) alpha.set(0, row, alpha_visited.get(0, idx++));
            }

            if (hasFork) {
                boolean anyOver = false;
                for (double rs : rowSums) {
                    if (rs > 1.0 + 1e-8) { anyOver = true; break; }
                }
                if (anyOver) {
                    for (int col = 0; col < alpha.getNumCols(); col++) {
                        if (alpha.get(0, col) > 1e-8) alpha.set(0, col, 1.0);
                    }
                }
            }

            Matrix visits_c = new Matrix(M, K);
            for (int ist = 0; ist < M; ist++) {
                for (int k = 0; k < inchain_c.size(); k++) {
                    visits_c.set(ist, inchain_c.get(k), alpha.get(0, ist * inchain_c.size() + k));
                }
            }

            double sum = 0.0;
            int row = (int) sn.stationToStateful.get((int) refstat.get(inchain_c.get(0)));
            for (int i = 0; i < inchain_c.size(); i++) {
                sum += visits_c.get(row, inchain_c.get(i));
            }
            // see _kb/03-api-layer.md for rationale
            if (sum < GlobalConstants.Zero && countTrue == 0) {
                // see _kb/03-api-layer.md for rationale
                for (int i = 0; i < inchain_c.size(); i++) {
                    visits_c.set(row, inchain_c.get(i), 1.0);
                }
                sum = 1.0;
            }
            Matrix visits_c_divide;
            if (sum > GlobalConstants.FineTol) {
                visits_c_divide = new Matrix(0, 0);
                visits_c.divide(sum, visits_c_divide, true);
            } else {
                // Reducible/absorbing chain: leave visits unnormalized, matching
                // MATLAB's "if normSum > FineTol, visits = visits / normSum" guard.
                visits_c_divide = visits_c;
            }
            visits_c_divide.absEq();
            // see _kb/03-api-layer.md for rationale
            for (int i = 0; i < visits_c_divide.getNumRows(); i++) {
                for (int j = 0; j < visits_c_divide.getNumCols(); j++) {
                    if (visits_c_divide.get(i, j) < GlobalConstants.Zero) {
                        visits_c_divide.set(i, j, 0.0);
                    }
                }
            }
            visits.put(c, visits_c_divide);
        }

        // Node visits
        Map<Integer, Matrix> nodeVisits = new HashMap<Integer, Matrix>();
        for (int c = 0; c < nchains; c++) {
            List<Integer> inchain_c = new_inchain.get(c);
            List<Integer> nodes_cols = new ArrayList<Integer>();
            for (int ind = 0; ind < I; ind++) {
                for (int ik = 0; ik < inchain_c.size(); ik++) {
                    nodes_cols.add(ind * K + inchain_c.get(ik));
                }
            }

            Matrix nodes_Pchain = new Matrix(nodes_cols.size(), nodes_cols.size());
            for (int row = 0; row < nodes_cols.size(); row++) {
                for (int col = 0; col < nodes_cols.size(); col++) {
                    nodes_Pchain.set(row, col, rtnodes.get(nodes_cols.get(row), nodes_cols.get(col)));
                }
            }

            replaceNaNWithEqualProb(nodes_Pchain);

            double[] nodesRowSums = new double[nodes_Pchain.getNumRows()];
            for (int i = 0; i < nodesRowSums.length; i++) nodesRowSums[i] = 1.0;
            if (hasFork) {
                for (int row = 0; row < nodes_Pchain.getNumRows(); row++) {
                    double rs = nodes_Pchain.sumRows(row);
                    nodesRowSums[row] = rs;
                    if (rs > 1e-8) {
                        for (int col = 0; col < nodes_Pchain.getNumCols(); col++) {
                            nodes_Pchain.set(row, col, nodes_Pchain.get(row, col) / rs);
                        }
                    }
                }
            }

            Matrix nodes_visited = new Matrix(nodes_Pchain.getNumRows(), 1);
            int countTrue = 0;
            for (int row = 0; row < nodes_Pchain.getNumRows(); row++) {
                if (nodes_Pchain.sumRows(row) > 0) {
                    countTrue++;
                    nodes_visited.set(row, 0, 1.0);
                }
            }

            Matrix input = new Matrix(countTrue, countTrue);
            int row_input = 0;
            for (int row = 0; row < nodes_visited.getNumRows(); row++) {
                if (nodes_visited.get(row, 0) > 0) {
                    int col_input = 0;
                    for (int col = 0; col < nodes_visited.getNumRows(); col++) {
                        if (nodes_visited.get(col, 0) > 0) {
                            input.set(row_input, col_input, nodes_Pchain.get(row, col));
                            col_input++;
                        }
                    }
                    row_input++;
                }
            }

            Matrix nodes_alpha_visited;
            try {
                nodes_alpha_visited = Dtmc_solve.dtmc_solve(input);
                if (nodes_alpha_visited.elementMax() == 0.0 || nodes_alpha_visited.hasNaN()) {
                    nodes_alpha_visited = Dtmc_solve_reducible.dtmc_solve_reducible(input).getLeft();
                }
            } catch (Exception e) {
                nodes_alpha_visited = Dtmc_solve_reducible.dtmc_solve_reducible(input).getLeft();
            }

            Matrix nodes_alpha = new Matrix(1, I * K);
            int idx = 0;
            for (int row = 0; row < nodes_visited.getNumRows(); row++) {
                if (nodes_visited.get(row, 0) > 0) nodes_alpha.set(0, row, nodes_alpha_visited.get(0, idx++));
            }

            boolean anyOver = false;
            for (double rs : nodesRowSums) {
                if (rs > 1.0 + 1e-8) { anyOver = true; break; }
            }
            if (hasFork && anyOver) {
                int nIC = inchain_c.size();
                for (int col = 0; col < nodes_alpha.getNumCols(); col++) {
                    if (nodes_alpha.get(0, col) > 1e-8) {
                        int nd = col / nIC;
                        if (nd < I && sn.nodetype.get(nd) == NodeType.Join) {
                            int r = inchain_c.get(col % nIC);
                            int rtCol = nd * K + r;
                            int nSources = 0;
                            for (int row = 0; row < rtnodes.getNumRows(); row++) {
                                if (rtnodes.get(row, rtCol) > 1e-8) nSources++;
                            }
                            nodes_alpha.set(0, col, (double) nSources);
                        } else {
                            nodes_alpha.set(0, col, 1.0);
                        }
                    }
                }
            }

            Matrix node_visits_c = new Matrix(I, K);
            for (int ind = 0; ind < I; ind++) {
                for (int k = 0; k < inchain_c.size(); k++) {
                    node_visits_c.set(ind, inchain_c.get(k), nodes_alpha.get(0, ind * inchain_c.size() + k));
                }
            }

            int ref = (int) refstat.get(inchain_c.get(0));
            double sum = 0.0;
            for (int k = 0; k < inchain_c.size(); k++) {
                sum += node_visits_c.get((int) sn.statefulToNode.get(ref), inchain_c.get(k));
            }
            if (sum < GlobalConstants.Zero && countTrue == 0) {
                // see _kb/03-api-layer.md for rationale
                int refNode = (int) sn.statefulToNode.get(ref);
                for (int k = 0; k < inchain_c.size(); k++) {
                    node_visits_c.set(refNode, inchain_c.get(k), 1.0);
                }
                sum = 1.0;
            }
            Matrix node_visits_c_divide = new Matrix(0, 0);
            node_visits_c.divide(sum, node_visits_c_divide, true);

            for (int i = 0; i < node_visits_c_divide.getNumRows(); i++) {
                for (int j = 0; j < node_visits_c_divide.getNumCols(); j++) {
                    if (node_visits_c_divide.get(i, j) < 0) node_visits_c_divide.set(i, j, 0.0);
                }
            }
            for (int i = 0; i < node_visits_c_divide.getNumRows(); i++) {
                for (int j = 0; j < node_visits_c_divide.getNumCols(); j++) {
                    if (Double.isNaN(node_visits_c_divide.get(i, j))) node_visits_c_divide.set(i, j, 0.0);
                }
            }
            nodeVisits.put(c, node_visits_c_divide);
        }

        sn.visits = visits;
        sn.nodevisits = nodeVisits;
        sn.refstat = refstat;
        return sn;
    }

    /**
     * Replace NaN routing entries with equal probabilities, in place.
     *
     * A Cache leaves its hit/miss routing unknown (NaN) until the cache itself is
     * solved, but the visit ratios are needed first. Ports MATLAB
     * sn_refresh_visits.m: spread the probability mass the row is missing equally
     * over its NaN entries. Without this the NaN row sums to NaN, `sum > 0` is
     * false, the row is dropped from `visited`, and dropping it strips the only
     * outgoing edge of another state -- leaving an absorbing row that makes the
     * chain degenerate and its stationary distribution meaningless.
     */
    private static void replaceNaNWithEqualProb(Matrix Pchain) {
        for (int row = 0; row < Pchain.getNumRows(); row++) {
            int nNaN = 0;
            double nonNaNSum = 0.0;
            for (int col = 0; col < Pchain.getNumCols(); col++) {
                double v = Pchain.get(row, col);
                if (Double.isNaN(v)) {
                    nNaN++;
                } else {
                    nonNaNSum += v;
                }
            }
            if (nNaN == 0) continue;
            double remaining = Math.max(0.0, 1.0 - nonNaNSum);
            double share = (remaining > 0) ? remaining / nNaN : 0.0;
            for (int col = 0; col < Pchain.getNumCols(); col++) {
                if (Double.isNaN(Pchain.get(row, col))) Pchain.set(row, col, share);
            }
        }
    }

}
