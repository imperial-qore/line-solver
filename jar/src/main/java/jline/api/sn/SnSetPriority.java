/**
 * @file Class Priority Modification for NetworkStruct
 *
 * Provides functions to directly modify class priorities
 * in a NetworkStruct without rebuilding the full Network object.
 *
 * @since LINE 3.0
 */
package jline.api.sn;

import java.util.ArrayList;
import java.util.List;

import jline.lang.NetworkStruct;
import jline.util.matrix.Matrix;

public final class SnSetPriority {
    private SnSetPriority() {}

    /**
     * Sets the priority for a class.
     */
    public static NetworkStruct snSetPriority(NetworkStruct sn, int classIdx, double priority,
                                              ModifyMode mode, ValidationLevel validation) {
        if (validation != ValidationLevel.NONE) {
            List<String> errors = new ArrayList<String>();

            String validationErr = SnValidate.snValidateClassIndex(sn, classIdx, "classIdx");
            if (validationErr != null) {
                errors.add(validationErr);
            }

            if (validation == ValidationLevel.FULL) {
                if (Double.isNaN(priority)) {
                    errors.add("priority is NaN");
                }
            }

            if (!errors.isEmpty()) {
                throw new SnValidationException("snSetPriority validation failed", errors);
            }
        }

        NetworkStruct snWork = (mode == ModifyMode.COPY) ? sn.<NetworkStruct>copy() : sn;

        setPriorityValue(snWork.classprio, classIdx, priority);

        return snWork;
    }

    public static NetworkStruct snSetPriority(NetworkStruct sn, int classIdx, double priority, ModifyMode mode) {
        return snSetPriority(sn, classIdx, priority, mode, ValidationLevel.MINIMAL);
    }

    public static NetworkStruct snSetPriority(NetworkStruct sn, int classIdx, double priority) {
        return snSetPriority(sn, classIdx, priority, ModifyMode.IN_PLACE, ValidationLevel.MINIMAL);
    }

    /**
     * Sets priorities for multiple classes in a single operation.
     */
    public static NetworkStruct snSetPriorityBatch(NetworkStruct sn, Matrix priorities,
                                                   ModifyMode mode, ValidationLevel validation) {
        if (validation != ValidationLevel.NONE) {
            List<String> errors = new ArrayList<String>();

            int totalElements = priorities.getNumRows() * priorities.getNumCols();
            if (totalElements != sn.nclasses) {
                errors.add("priorities matrix total elements (" + totalElements
                        + ") do not match nclasses (" + sn.nclasses + ")");
            }

            if (!errors.isEmpty()) {
                throw new SnValidationException("snSetPriorityBatch validation failed", errors);
            }
        }

        NetworkStruct snWork = (mode == ModifyMode.COPY) ? sn.<NetworkStruct>copy() : sn;

        for (int k = 0; k < sn.nclasses; k++) {
            double p = getMatrixElement(priorities, k);
            if (!Double.isNaN(p)) {
                setPriorityValue(snWork.classprio, k, p);
            }
        }

        return snWork;
    }

    public static NetworkStruct snSetPriorityBatch(NetworkStruct sn, Matrix priorities, ModifyMode mode) {
        return snSetPriorityBatch(sn, priorities, mode, ValidationLevel.MINIMAL);
    }

    public static NetworkStruct snSetPriorityBatch(NetworkStruct sn, Matrix priorities) {
        return snSetPriorityBatch(sn, priorities, ModifyMode.IN_PLACE, ValidationLevel.MINIMAL);
    }

    /**
     * Sets priority value in classprio matrix.
     * Handles both (1 x K) and (K x 1) layouts.
     */
    private static void setPriorityValue(Matrix classprio, int classIdx, double value) {
        if (classprio == null) return;

        if (classprio.getNumRows() == 1) {
            classprio.set(0, classIdx, value);
        } else {
            classprio.set(classIdx, 0, value);
        }
    }

    /**
     * Gets element from a matrix that may be either (1 x K) or (K x 1).
     */
    private static double getMatrixElement(Matrix m, int idx) {
        if (m.getNumRows() == 1) {
            return m.get(0, idx);
        } else {
            return m.get(idx, 0);
        }
    }
}
