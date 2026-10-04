/**
 * @file Job Population Modification for NetworkStruct
 *
 * Provides functions to directly modify the number of jobs for closed classes
 * in a NetworkStruct without rebuilding the full Network object.
 *
 * Mirrors MATLAB implementation patterns for population modification.
 *
 * @since LINE 3.0
 */
package jline.api.sn;

import java.util.ArrayList;
import java.util.List;

import jline.lang.NetworkStruct;
import jline.util.matrix.Matrix;

public final class SnSetPopulation {
    private SnSetPopulation() {}

    /**
     * Sets the number of jobs for a closed class.
     */
    public static NetworkStruct snSetPopulation(NetworkStruct sn, int classIdx, double nJobs,
                                                ModifyMode mode, ValidationLevel validation,
                                                boolean autoRefresh) {
        // Validation
        if (validation != ValidationLevel.NONE) {
            List<String> errors = new ArrayList<String>();

            String validationErr = SnValidate.snValidateClassIndex(sn, classIdx, "classIdx");
            if (validationErr != null) {
                errors.add(validationErr);
            }

            if (validation == ValidationLevel.FULL) {
                if (Double.isNaN(nJobs)) {
                    errors.add("nJobs is NaN");
                } else if (nJobs < 0) {
                    errors.add("nJobs=" + nJobs + " must be non-negative");
                }
                // Check if this is a closed class (should not be Inf)
                double currentNJobs = getPopulationValue(sn.njobs, classIdx);
                if (Double.isInfinite(currentNJobs) && Double.isFinite(nJobs)) {
                    errors.add("Class " + classIdx + " is an open class (Inf jobs). Cannot set finite population.");
                }
            }

            if (!errors.isEmpty()) {
                throw new SnValidationException("snSetPopulation validation failed", errors);
            }
        }

        // Get working copy if needed
        NetworkStruct snWork = (mode == ModifyMode.COPY) ? sn.<NetworkStruct>copy() : sn;

        // Update njobs matrix
        setPopulationValue(snWork.njobs, classIdx, nJobs);

        // Recalculate nclosedjobs
        snWork.nclosedjobs = calculateTotalClosedJobs(snWork.njobs);

        // Auto-refresh visit ratios if requested
        if (autoRefresh && snWork.chains != null && snWork.rt != null && snWork.rtnodes != null) {
            SnRefreshVisits.snRefreshVisits(snWork, snWork.chains, snWork.rt, snWork.rtnodes);
        }

        return snWork;
    }

    public static NetworkStruct snSetPopulation(NetworkStruct sn, int classIdx, double nJobs,
                                                ModifyMode mode, ValidationLevel validation) {
        return snSetPopulation(sn, classIdx, nJobs, mode, validation, false);
    }

    public static NetworkStruct snSetPopulation(NetworkStruct sn, int classIdx, double nJobs, ModifyMode mode) {
        return snSetPopulation(sn, classIdx, nJobs, mode, ValidationLevel.MINIMAL, false);
    }

    public static NetworkStruct snSetPopulation(NetworkStruct sn, int classIdx, double nJobs) {
        return snSetPopulation(sn, classIdx, nJobs, ModifyMode.IN_PLACE, ValidationLevel.MINIMAL, false);
    }

    /**
     * Sets populations for multiple classes in a single operation.
     */
    public static NetworkStruct snSetPopulationBatch(NetworkStruct sn, Matrix njobs,
                                                     ModifyMode mode, ValidationLevel validation,
                                                     boolean autoRefresh) {
        // Validation
        if (validation != ValidationLevel.NONE) {
            List<String> errors = new ArrayList<String>();

            int totalElements = njobs.getNumRows() * njobs.getNumCols();
            if (totalElements != sn.nclasses) {
                errors.add("njobs matrix total elements (" + totalElements
                        + ") do not match nclasses (" + sn.nclasses + ")");
            }

            if (validation == ValidationLevel.FULL) {
                for (int k = 0; k < sn.nclasses; k++) {
                    double n = getMatrixElement(njobs, k);
                    if (!Double.isNaN(n) && n < 0) {
                        errors.add("njobs[" + k + "]=" + n + " must be non-negative");
                    }
                }
            }

            if (!errors.isEmpty()) {
                throw new SnValidationException("snSetPopulationBatch validation failed", errors);
            }
        }

        // Get working copy if needed
        NetworkStruct snWork = (mode == ModifyMode.COPY) ? sn.<NetworkStruct>copy() : sn;

        // Update njobs matrix
        for (int k = 0; k < sn.nclasses; k++) {
            double n = getMatrixElement(njobs, k);
            if (!Double.isNaN(n)) {
                setPopulationValue(snWork.njobs, k, n);
            }
        }

        // Recalculate nclosedjobs
        snWork.nclosedjobs = calculateTotalClosedJobs(snWork.njobs);

        // Auto-refresh visit ratios if requested
        if (autoRefresh && snWork.chains != null && snWork.rt != null && snWork.rtnodes != null) {
            SnRefreshVisits.snRefreshVisits(snWork, snWork.chains, snWork.rt, snWork.rtnodes);
        }

        return snWork;
    }

    public static NetworkStruct snSetPopulationBatch(NetworkStruct sn, Matrix njobs,
                                                     ModifyMode mode, ValidationLevel validation) {
        return snSetPopulationBatch(sn, njobs, mode, validation, false);
    }

    public static NetworkStruct snSetPopulationBatch(NetworkStruct sn, Matrix njobs, ModifyMode mode) {
        return snSetPopulationBatch(sn, njobs, mode, ValidationLevel.MINIMAL, false);
    }

    public static NetworkStruct snSetPopulationBatch(NetworkStruct sn, Matrix njobs) {
        return snSetPopulationBatch(sn, njobs, ModifyMode.IN_PLACE, ValidationLevel.MINIMAL, false);
    }

    /**
     * Gets population value from njobs matrix.
     * Handles both (1 x K) and (K x 1) layouts.
     */
    private static double getPopulationValue(Matrix njobs, int classIdx) {
        if (njobs == null) return Double.POSITIVE_INFINITY;

        if (njobs.getNumRows() == 1) {
            return njobs.get(0, classIdx);
        } else {
            return njobs.get(classIdx, 0);
        }
    }

    /**
     * Sets population value in njobs matrix.
     * Handles both (1 x K) and (K x 1) layouts.
     */
    private static void setPopulationValue(Matrix njobs, int classIdx, double value) {
        if (njobs == null) return;

        if (njobs.getNumRows() == 1) {
            njobs.set(0, classIdx, value);
        } else {
            njobs.set(classIdx, 0, value);
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

    /**
     * Calculates total number of closed jobs (sum of finite populations).
     */
    private static int calculateTotalClosedJobs(Matrix njobs) {
        if (njobs == null) return 0;

        double total = 0.0;
        int totalElements = njobs.getNumRows() * njobs.getNumCols();
        for (int i = 0; i < totalElements; i++) {
            double value;
            if (njobs.getNumRows() == 1) {
                value = njobs.get(0, i);
            } else {
                value = njobs.get(i, 0);
            }
            if (Double.isFinite(value)) {
                total += value;
            }
        }
        return (int) total;
    }
}
