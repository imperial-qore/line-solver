/**
 * @file Arrival Rate Modification for NetworkStruct
 *
 * Provides functions to directly modify arrival rates at the Source station
 * in a NetworkStruct without rebuilding the full Network object.
 * Delegates to snSetService with the Source station index.
 *
 * @since LINE 3.0
 */
package jline.api.sn;

import java.util.ArrayList;
import java.util.List;

import jline.lang.NetworkStruct;
import jline.lang.constant.NodeType;
import jline.util.matrix.Matrix;

public final class SnSetArrival {
    private SnSetArrival() {}

    /**
     * Sets the arrival rate for a class at the Source station.
     */
    public static NetworkStruct snSetArrival(NetworkStruct sn, int classIdx, double rate, double scv,
                                              ModifyMode mode, ValidationLevel validation, boolean autoRefresh) {
        // Find Source station index
        int sourceStationIdx = findSourceStationIndex(sn);

        if (sourceStationIdx < 0) {
            if (validation != ValidationLevel.NONE) {
                throw new SnValidationException("snSetArrival: No Source station found in network");
            }
            return sn;
        }

        // Validate that this is an open class (has Inf jobs)
        if (validation == ValidationLevel.FULL) {
            double njob = getPopulationValue(sn, classIdx);
            if (Double.isFinite(njob)) {
                throw new SnValidationException(
                        "snSetArrival: Class " + classIdx + " is a closed class (njobs=" + njob + "). " +
                                "Arrival rates can only be set for open classes.");
            }
        }

        // Delegate to snSetService
        return SnSetService.snSetService(sn, sourceStationIdx, classIdx, rate, scv, mode, validation, autoRefresh);
    }

    public static NetworkStruct snSetArrival(NetworkStruct sn, int classIdx, double rate, double scv,
                                              ModifyMode mode, ValidationLevel validation) {
        return snSetArrival(sn, classIdx, rate, scv, mode, validation, false);
    }

    public static NetworkStruct snSetArrival(NetworkStruct sn, int classIdx, double rate, double scv, ModifyMode mode) {
        return snSetArrival(sn, classIdx, rate, scv, mode, ValidationLevel.MINIMAL, false);
    }

    public static NetworkStruct snSetArrival(NetworkStruct sn, int classIdx, double rate, double scv) {
        return snSetArrival(sn, classIdx, rate, scv, ModifyMode.IN_PLACE, ValidationLevel.MINIMAL, false);
    }

    public static NetworkStruct snSetArrival(NetworkStruct sn, int classIdx, double rate) {
        return snSetArrival(sn, classIdx, rate, 1.0, ModifyMode.IN_PLACE, ValidationLevel.MINIMAL, false);
    }

    /**
     * Sets arrival rates for multiple classes in a single operation.
     */
    public static NetworkStruct snSetArrivalBatch(NetworkStruct sn, Matrix rates, Matrix scvs,
                                                   ModifyMode mode, ValidationLevel validation,
                                                   boolean autoRefresh) {
        int sourceStationIdx = findSourceStationIndex(sn);

        if (sourceStationIdx < 0) {
            if (validation != ValidationLevel.NONE) {
                throw new SnValidationException("snSetArrivalBatch: No Source station found in network");
            }
            return sn;
        }

        if (validation != ValidationLevel.NONE) {
            List<String> errors = new ArrayList<String>();
            int totalElements = rates.getNumRows() * rates.getNumCols();

            if (totalElements != sn.nclasses) {
                errors.add("rates matrix total elements (" + totalElements
                        + ") do not match nclasses (" + sn.nclasses + ")");
            }

            if (scvs != null) {
                int scvElements = scvs.getNumRows() * scvs.getNumCols();
                if (scvElements != sn.nclasses) {
                    errors.add("scvs matrix total elements (" + scvElements
                            + ") do not match nclasses (" + sn.nclasses + ")");
                }
            }

            if (validation == ValidationLevel.FULL) {
                for (int k = 0; k < sn.nclasses; k++) {
                    double r = getMatrixElement(rates, k);
                    if (!Double.isNaN(r) && r < 0) {
                        errors.add("rates[" + k + "]=" + r + " must be non-negative");
                    }
                    if (scvs != null) {
                        double s = getMatrixElement(scvs, k);
                        if (!Double.isNaN(s) && s < 0) {
                            errors.add("scvs[" + k + "]=" + s + " must be non-negative");
                        }
                    }
                }
            }

            if (!errors.isEmpty()) {
                throw new SnValidationException("snSetArrivalBatch validation failed", errors);
            }
        }

        NetworkStruct snWork = (mode == ModifyMode.COPY) ? sn.<NetworkStruct>copy() : sn;

        // Update rates at source station
        for (int k = 0; k < sn.nclasses; k++) {
            double r = getMatrixElement(rates, k);
            if (!Double.isNaN(r) && snWork.rates != null) {
                snWork.rates.set(sourceStationIdx, k, r);
            }
            if (scvs != null) {
                double s = getMatrixElement(scvs, k);
                if (!Double.isNaN(s) && snWork.scv != null) {
                    snWork.scv.set(sourceStationIdx, k, s);
                }
            }
        }

        // Auto-refresh if requested
        if (autoRefresh) {
            for (int k = 0; k < sn.nclasses; k++) {
                double r = getMatrixElement(rates, k);
                if (!Double.isNaN(r)) {
                    SnRefreshProcessFields.snRefreshProcessFields(snWork, sourceStationIdx, k);
                }
            }
        }

        return snWork;
    }

    public static NetworkStruct snSetArrivalBatch(NetworkStruct sn, Matrix rates, Matrix scvs, ModifyMode mode,
                                                   ValidationLevel validation) {
        return snSetArrivalBatch(sn, rates, scvs, mode, validation, false);
    }

    public static NetworkStruct snSetArrivalBatch(NetworkStruct sn, Matrix rates, Matrix scvs, ModifyMode mode) {
        return snSetArrivalBatch(sn, rates, scvs, mode, ValidationLevel.MINIMAL, false);
    }

    public static NetworkStruct snSetArrivalBatch(NetworkStruct sn, Matrix rates, Matrix scvs) {
        return snSetArrivalBatch(sn, rates, scvs, ModifyMode.IN_PLACE, ValidationLevel.MINIMAL, false);
    }

    public static NetworkStruct snSetArrivalBatch(NetworkStruct sn, Matrix rates) {
        return snSetArrivalBatch(sn, rates, null, ModifyMode.IN_PLACE, ValidationLevel.MINIMAL, false);
    }

    /**
     * Finds the index of the Source station in the network.
     */
    private static int findSourceStationIndex(NetworkStruct sn) {
        // First check nodetype for Source
        for (int i = 0; i < sn.nodetype.size(); i++) {
            if (sn.nodetype.get(i) == NodeType.Source) {
                int stationIdx = (sn.nodeToStation != null) ? (int) sn.nodeToStation.get(0, i) : -1;
                if (stationIdx >= 0) {
                    return stationIdx;
                }
            }
        }

        // Fallback: Source is typically the first station (index 0)
        if (sn.nstations > 0 && !sn.nodetype.isEmpty()) {
            NodeType firstNodeType = sn.nodetype.get(0);
            if (firstNodeType == NodeType.Source) {
                return 0;
            }
        }

        return -1;
    }

    /**
     * Gets population value for a class from njobs matrix.
     */
    private static double getPopulationValue(NetworkStruct sn, int classIdx) {
        if (sn.njobs == null) return Double.POSITIVE_INFINITY;

        if (sn.njobs.getNumRows() == 1) {
            return sn.njobs.get(0, classIdx);
        } else {
            return sn.njobs.get(classIdx, 0);
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
