/**
 * @file Server Count Modification for NetworkStruct
 *
 * @since LINE 3.0
 */
package jline.api.sn;

import java.util.ArrayList;
import java.util.List;

import jline.lang.NetworkStruct;
import jline.lang.constant.NodeType;
import jline.util.matrix.Matrix;

public final class SnSetServers {
    private SnSetServers() {}

    /**
     * Sets the number of servers at a station.
     */
    public static NetworkStruct snSetServers(NetworkStruct sn, int stationIdx, double nServers,
                                             ModifyMode mode, ValidationLevel validation) {
        if (validation != ValidationLevel.NONE) {
            List<String> errors = new ArrayList<String>();

            String validationErr = SnValidate.snValidateStationIndex(sn, stationIdx, "stationIdx");
            if (validationErr != null) {
                errors.add(validationErr);
            }

            if (validation == ValidationLevel.FULL) {
                if (Double.isNaN(nServers)) {
                    errors.add("nServers is NaN");
                } else if (nServers <= 0 && !Double.isInfinite(nServers)) {
                    int nodeIdx = (sn.stationToNode != null) ? (int) sn.stationToNode.get(stationIdx, 0) : stationIdx;
                    NodeType nodeType = (nodeIdx >= 0 && nodeIdx < sn.nodetype.size()) ? sn.nodetype.get(nodeIdx) : null;
                    if (nodeType != NodeType.Source && nodeType != NodeType.Sink) {
                        errors.add("nServers=" + nServers + " must be positive for station type " + nodeType);
                    }
                }
            }

            if (!errors.isEmpty()) {
                throw new SnValidationException("snSetServers validation failed", errors);
            }
        }

        NetworkStruct snWork = (mode == ModifyMode.COPY) ? sn.<NetworkStruct>copy() : sn;

        if (snWork.nservers != null) {
            snWork.nservers.set(stationIdx, 0, nServers);
        }

        return snWork;
    }

    public static NetworkStruct snSetServers(NetworkStruct sn, int stationIdx, double nServers, ModifyMode mode) {
        return snSetServers(sn, stationIdx, nServers, mode, ValidationLevel.MINIMAL);
    }

    public static NetworkStruct snSetServers(NetworkStruct sn, int stationIdx, double nServers) {
        return snSetServers(sn, stationIdx, nServers, ModifyMode.IN_PLACE, ValidationLevel.MINIMAL);
    }

    /**
     * Sets the number of servers for multiple stations in a single operation.
     */
    public static NetworkStruct snSetServersBatch(NetworkStruct sn, Matrix nServers,
                                                   ModifyMode mode, ValidationLevel validation) {
        if (validation != ValidationLevel.NONE) {
            List<String> errors = new ArrayList<String>();

            int totalElements = nServers.getNumRows() * nServers.getNumCols();
            if (totalElements != sn.nstations) {
                errors.add("nServers matrix total elements (" + totalElements
                        + ") do not match nstations (" + sn.nstations + ")");
            }

            if (validation == ValidationLevel.FULL) {
                for (int i = 0; i < sn.nstations; i++) {
                    double n = getMatrixElement(nServers, i);
                    if (!Double.isNaN(n) && n <= 0 && !Double.isInfinite(n)) {
                        int nodeIdx = (sn.stationToNode != null) ? (int) sn.stationToNode.get(i, 0) : i;
                        NodeType nodeType = (nodeIdx >= 0 && nodeIdx < sn.nodetype.size()) ? sn.nodetype.get(nodeIdx) : null;
                        if (nodeType != NodeType.Source && nodeType != NodeType.Sink) {
                            errors.add("nServers[" + i + "]=" + n + " must be positive");
                        }
                    }
                }
            }

            if (!errors.isEmpty()) {
                throw new SnValidationException("snSetServersBatch validation failed", errors);
            }
        }

        NetworkStruct snWork = (mode == ModifyMode.COPY) ? sn.<NetworkStruct>copy() : sn;

        for (int i = 0; i < sn.nstations; i++) {
            double n = getMatrixElement(nServers, i);
            if (!Double.isNaN(n)) {
                if (snWork.nservers != null) {
                    snWork.nservers.set(i, 0, n);
                }
            }
        }

        return snWork;
    }

    public static NetworkStruct snSetServersBatch(NetworkStruct sn, Matrix nServers, ModifyMode mode) {
        return snSetServersBatch(sn, nServers, mode, ValidationLevel.MINIMAL);
    }

    public static NetworkStruct snSetServersBatch(NetworkStruct sn, Matrix nServers) {
        return snSetServersBatch(sn, nServers, ModifyMode.IN_PLACE, ValidationLevel.MINIMAL);
    }

    /**
     * Gets element from a matrix that may be either (M x 1) or (1 x M).
     */
    private static double getMatrixElement(Matrix m, int idx) {
        if (m.getNumCols() == 1) {
            return m.get(idx, 0);
        } else {
            return m.get(0, idx);
        }
    }
}
