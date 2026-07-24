/**
 * @file Routing Matrix Modification for NetworkStruct
 *
 * Provides functions to directly modify routing matrices (rt, rtnodes)
 * in a NetworkStruct without rebuilding the full Network object.
 *
 * @since LINE 3.0
 */
package jline.api.sn;

import java.util.ArrayList;
import java.util.List;

import jline.lang.NetworkStruct;
import jline.util.matrix.Matrix;

public final class SnSetRouting {
    private SnSetRouting() {}

    /**
     * Sets a routing probability between two stateful node-class pairs.
     */
    public static NetworkStruct snSetRoutingProb(NetworkStruct sn, int fromStateful, int fromClass,
                                                  int toStateful, int toClass, double probability,
                                                  ModifyMode mode, ValidationLevel validation,
                                                  boolean autoRefresh) {
        if (validation != ValidationLevel.NONE) {
            List<String> errors = new ArrayList<String>();

            if (fromStateful < 0 || fromStateful >= sn.nstateful) {
                errors.add("fromStateful=" + fromStateful + " is out of bounds [0, " + (sn.nstateful - 1) + "]");
            }
            if (toStateful < 0 || toStateful >= sn.nstateful) {
                errors.add("toStateful=" + toStateful + " is out of bounds [0, " + (sn.nstateful - 1) + "]");
            }
            String e1 = SnValidate.snValidateClassIndex(sn, fromClass, "fromClass");
            if (e1 != null) errors.add(e1);
            String e2 = SnValidate.snValidateClassIndex(sn, toClass, "toClass");
            if (e2 != null) errors.add(e2);

            if (validation == ValidationLevel.FULL) {
                if (Double.isNaN(probability)) {
                    errors.add("probability is NaN");
                } else if (probability < 0 || probability > 1) {
                    errors.add("probability=" + probability + " must be in [0, 1]");
                }
            }

            if (!errors.isEmpty()) {
                throw new SnValidationException("snSetRoutingProb validation failed", errors);
            }
        }

        NetworkStruct snWork = (mode == ModifyMode.COPY) ? sn.<NetworkStruct>copy() : sn;

        // Calculate indices in rt matrix
        int fromIdx = fromStateful * sn.nclasses + fromClass;
        int toIdx = toStateful * sn.nclasses + toClass;

        // Update rt matrix
        if (snWork.rt != null) {
            snWork.rt.set(fromIdx, toIdx, probability);
        }

        // Auto-refresh visit ratios if requested
        if (autoRefresh && snWork.chains != null && snWork.rt != null && snWork.rtnodes != null) {
            SnRefreshVisits.snRefreshVisits(snWork, snWork.chains, snWork.rt, snWork.rtnodes);
        }

        return snWork;
    }

    public static NetworkStruct snSetRoutingProb(NetworkStruct sn, int fromStateful, int fromClass,
                                                  int toStateful, int toClass, double probability,
                                                  ModifyMode mode, ValidationLevel validation) {
        return snSetRoutingProb(sn, fromStateful, fromClass, toStateful, toClass, probability,
                mode, validation, false);
    }

    public static NetworkStruct snSetRoutingProb(NetworkStruct sn, int fromStateful, int fromClass,
                                                  int toStateful, int toClass, double probability,
                                                  ModifyMode mode) {
        return snSetRoutingProb(sn, fromStateful, fromClass, toStateful, toClass, probability,
                mode, ValidationLevel.MINIMAL, false);
    }

    public static NetworkStruct snSetRoutingProb(NetworkStruct sn, int fromStateful, int fromClass,
                                                  int toStateful, int toClass, double probability) {
        return snSetRoutingProb(sn, fromStateful, fromClass, toStateful, toClass, probability,
                ModifyMode.IN_PLACE, ValidationLevel.MINIMAL, false);
    }

    /**
     * Sets the entire routing matrix for stateful nodes.
     */
    public static NetworkStruct snSetRoutingMatrix(NetworkStruct sn, Matrix rt, ModifyMode mode,
                                                    ValidationLevel validation, boolean autoRefresh) {
        int expectedSize = sn.nstateful * sn.nclasses;

        if (validation != ValidationLevel.NONE) {
            List<String> errors = new ArrayList<String>();

            if (rt.getNumRows() != expectedSize || rt.getNumCols() != expectedSize) {
                errors.add("rt matrix dimensions (" + rt.getNumRows() + "x" + rt.getNumCols()
                        + ") do not match expected (" + expectedSize + " x " + expectedSize + ")");
            }

            if (validation == ValidationLevel.FULL) {
                validateStochasticMatrix(rt, "rt", errors);
            }

            if (!errors.isEmpty()) {
                throw new SnValidationException("snSetRoutingMatrix validation failed", errors);
            }
        }

        NetworkStruct snWork = (mode == ModifyMode.COPY) ? sn.<NetworkStruct>copy() : sn;

        // Replace rt matrix
        snWork.rt = rt;

        // Auto-refresh visit ratios if requested
        if (autoRefresh && snWork.chains != null && snWork.rtnodes != null) {
            SnRefreshVisits.snRefreshVisits(snWork, snWork.chains, snWork.rt, snWork.rtnodes);
        }

        return snWork;
    }

    public static NetworkStruct snSetRoutingMatrix(NetworkStruct sn, Matrix rt, ModifyMode mode,
                                                    ValidationLevel validation) {
        return snSetRoutingMatrix(sn, rt, mode, validation, false);
    }

    public static NetworkStruct snSetRoutingMatrix(NetworkStruct sn, Matrix rt, ModifyMode mode) {
        return snSetRoutingMatrix(sn, rt, mode, ValidationLevel.MINIMAL, false);
    }

    public static NetworkStruct snSetRoutingMatrix(NetworkStruct sn, Matrix rt) {
        return snSetRoutingMatrix(sn, rt, ModifyMode.IN_PLACE, ValidationLevel.MINIMAL, false);
    }

    /**
     * Sets the entire routing matrix for all nodes.
     */
    public static NetworkStruct snSetRoutingNodesMatrix(NetworkStruct sn, Matrix rtnodes, ModifyMode mode,
                                                         ValidationLevel validation, boolean autoRefresh) {
        int expectedSize = sn.nnodes * sn.nclasses;

        if (validation != ValidationLevel.NONE) {
            List<String> errors = new ArrayList<String>();

            if (rtnodes.getNumRows() != expectedSize || rtnodes.getNumCols() != expectedSize) {
                errors.add("rtnodes matrix dimensions (" + rtnodes.getNumRows() + "x" + rtnodes.getNumCols()
                        + ") do not match expected (" + expectedSize + " x " + expectedSize + ")");
            }

            if (validation == ValidationLevel.FULL) {
                validateStochasticMatrix(rtnodes, "rtnodes", errors);
            }

            if (!errors.isEmpty()) {
                throw new SnValidationException("snSetRoutingNodesMatrix validation failed", errors);
            }
        }

        NetworkStruct snWork = (mode == ModifyMode.COPY) ? sn.<NetworkStruct>copy() : sn;

        // Replace rtnodes matrix
        snWork.rtnodes = rtnodes;

        // Auto-refresh visit ratios if requested
        if (autoRefresh && snWork.chains != null && snWork.rt != null) {
            SnRefreshVisits.snRefreshVisits(snWork, snWork.chains, snWork.rt, snWork.rtnodes);
        }

        return snWork;
    }

    public static NetworkStruct snSetRoutingNodesMatrix(NetworkStruct sn, Matrix rtnodes, ModifyMode mode,
                                                         ValidationLevel validation) {
        return snSetRoutingNodesMatrix(sn, rtnodes, mode, validation, false);
    }

    public static NetworkStruct snSetRoutingNodesMatrix(NetworkStruct sn, Matrix rtnodes, ModifyMode mode) {
        return snSetRoutingNodesMatrix(sn, rtnodes, mode, ValidationLevel.MINIMAL, false);
    }

    public static NetworkStruct snSetRoutingNodesMatrix(NetworkStruct sn, Matrix rtnodes) {
        return snSetRoutingNodesMatrix(sn, rtnodes, ModifyMode.IN_PLACE, ValidationLevel.MINIMAL, false);
    }

    /**
     * Validates that a matrix is stochastic (row sums approximately equal 1).
     */
    private static void validateStochasticMatrix(Matrix m, String name, List<String> errors) {
        double tolerance = 1e-6;

        for (int i = 0; i < m.getNumRows(); i++) {
            double rowSum = 0.0;
            boolean hasNonZero = false;
            boolean hasNegative = false;

            for (int j = 0; j < m.getNumCols(); j++) {
                double value = m.get(i, j);
                if (Double.isNaN(value)) {
                    errors.add(name + "[" + i + "," + j + "] is NaN");
                    continue;
                }
                if (value < 0) {
                    hasNegative = true;
                }
                if (value > 0) {
                    hasNonZero = true;
                }
                rowSum += value;
            }

            if (hasNegative) {
                errors.add(name + " row " + i + " contains negative values");
            }

            // Only check row sum if row has non-zero entries (active row)
            if (hasNonZero && Math.abs(rowSum - 1.0) > tolerance) {
                errors.add(name + " row " + i + " sum = " + rowSum + " (expected 1.0)");
            }
        }
    }
}
