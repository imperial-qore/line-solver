/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
package jline.api.sn;

import java.util.ArrayList;
import java.util.List;

import jline.lang.NetworkStruct;
import jline.lang.constant.NodeType;
import jline.util.matrix.Matrix;

/**
 * NetworkStruct Validation Utilities.
 */
public final class SnValidate {
    private SnValidate() {}

    public static List<String> snValidate(NetworkStruct sn) {
        return snValidate(sn, ValidationLevel.FULL);
    }

    public static List<String> snValidate(NetworkStruct sn, ValidationLevel level) {
        if (level == ValidationLevel.NONE) {
            return new ArrayList<String>();
        }
        List<String> errors = new ArrayList<String>();
        errors.addAll(snValidateDimensions(sn));
        if (level == ValidationLevel.FULL) {
            errors.addAll(snValidateRates(sn));
            errors.addAll(snValidatePopulation(sn));
            errors.addAll(snValidateRouting(sn));
            errors.addAll(snValidateServers(sn));
        }
        return errors;
    }

    public static List<String> snValidateDimensions(NetworkStruct sn) {
        List<String> errors = new ArrayList<String>();
        int M = sn.nstations;
        int K = sn.nclasses;

        if (sn.rates != null) {
            if (sn.rates.getNumRows() != M || sn.rates.getNumCols() != K) {
                errors.add("rates matrix dimensions (" + sn.rates.getNumRows() + "x" + sn.rates.getNumCols()
                        + ") do not match (nstations=" + M + " x nclasses=" + K + ")");
            }
        }

        if (sn.scv != null) {
            if (sn.scv.getNumRows() != M || sn.scv.getNumCols() != K) {
                errors.add("scv matrix dimensions (" + sn.scv.getNumRows() + "x" + sn.scv.getNumCols()
                        + ") do not match (nstations=" + M + " x nclasses=" + K + ")");
            }
        }

        if (sn.nservers != null) {
            if (sn.nservers.getNumRows() != M || sn.nservers.getNumCols() != 1) {
                errors.add("nservers matrix dimensions (" + sn.nservers.getNumRows() + "x" + sn.nservers.getNumCols()
                        + ") do not match (nstations=" + M + " x 1)");
            }
        }

        if (sn.njobs != null) {
            int totalElements = sn.njobs.getNumRows() * sn.njobs.getNumCols();
            if (totalElements != K) {
                errors.add("njobs matrix total elements (" + totalElements + ") does not match nclasses=" + K);
            }
        }

        if (sn.classprio != null) {
            int totalElements = sn.classprio.getNumRows() * sn.classprio.getNumCols();
            if (totalElements != K) {
                errors.add("classprio matrix total elements (" + totalElements + ") does not match nclasses=" + K);
            }
        }

        if (sn.phases != null) {
            if (sn.phases.getNumRows() != M || sn.phases.getNumCols() != K) {
                errors.add("phases matrix dimensions (" + sn.phases.getNumRows() + "x" + sn.phases.getNumCols()
                        + ") do not match (nstations=" + M + " x nclasses=" + K + ")");
            }
        }
        return errors;
    }

    public static List<String> snValidateRates(NetworkStruct sn) {
        List<String> errors = new ArrayList<String>();
        if (sn.rates == null) return errors;

        for (int i = 0; i < sn.rates.getNumRows(); i++) {
            for (int j = 0; j < sn.rates.getNumCols(); j++) {
                double rate = sn.rates.get(i, j);
                if (Double.isNaN(rate)) continue;
                if (rate < 0) {
                    errors.add("rates[" + i + "," + j + "] = " + rate + " is negative");
                }
            }
        }

        if (sn.scv != null) {
            for (int i = 0; i < sn.scv.getNumRows(); i++) {
                for (int j = 0; j < sn.scv.getNumCols(); j++) {
                    double v = sn.scv.get(i, j);
                    if (Double.isNaN(v)) continue;
                    if (v < 0) {
                        errors.add("scv[" + i + "," + j + "] = " + v + " is negative (SCV must be >= 0)");
                    }
                }
            }
        }
        return errors;
    }

    public static List<String> snValidatePopulation(NetworkStruct sn) {
        List<String> errors = new ArrayList<String>();
        if (sn.njobs == null) return errors;
        for (int k = 0; k < sn.nclasses; k++) {
            double njob = getPopulation(sn.njobs, k);
            if (Double.isNaN(njob)) {
                errors.add("njobs for class " + k + " is NaN");
            } else if (Double.isFinite(njob) && njob < 0) {
                errors.add("njobs for class " + k + " = " + njob + " is negative");
            }
        }
        return errors;
    }

    public static List<String> snValidateRouting(NetworkStruct sn) {
        List<String> errors = new ArrayList<String>();
        if (sn.rt == null) return errors;
        double tolerance = 1e-6;

        for (int i = 0; i < sn.rt.getNumRows(); i++) {
            double rowSum = 0.0;
            boolean hasNonZero = false;
            for (int j = 0; j < sn.rt.getNumCols(); j++) {
                double v = sn.rt.get(i, j);
                if (Double.isNaN(v)) {
                    errors.add("rt[" + i + "," + j + "] is NaN");
                    continue;
                }
                if (v < 0) errors.add("rt[" + i + "," + j + "] = " + v + " is negative");
                if (v > 0) hasNonZero = true;
                rowSum += v;
            }
            if (hasNonZero && Math.abs(rowSum - 1.0) > tolerance) {
                errors.add("rt row " + i + " sum = " + rowSum + " (expected 1.0)");
            }
        }
        return errors;
    }

    public static List<String> snValidateServers(NetworkStruct sn) {
        List<String> errors = new ArrayList<String>();
        if (sn.nservers == null) return errors;

        for (int i = 0; i < sn.nservers.getNumRows(); i++) {
            double n = sn.nservers.get(i, 0);
            if (Double.isNaN(n)) {
                errors.add("nservers[" + i + "] is NaN");
            } else if (n <= 0 && !Double.isInfinite(n)) {
                NodeType nodeType = (i < sn.nodetype.size()) ? sn.nodetype.get(i) : null;
                if (nodeType != NodeType.Source && nodeType != NodeType.Sink) {
                    errors.add("nservers[" + i + "] = " + n + " must be positive");
                }
            }
        }
        return errors;
    }

    public static String snValidateStationIndex(NetworkStruct sn, int stationIdx) {
        return snValidateStationIndex(sn, stationIdx, "stationIdx");
    }

    public static String snValidateStationIndex(NetworkStruct sn, int stationIdx, String paramName) {
        if (stationIdx < 0 || stationIdx >= sn.nstations) {
            return paramName + "=" + stationIdx + " is out of bounds [0, " + (sn.nstations - 1) + "]";
        }
        return null;
    }

    public static String snValidateClassIndex(NetworkStruct sn, int classIdx) {
        return snValidateClassIndex(sn, classIdx, "classIdx");
    }

    public static String snValidateClassIndex(NetworkStruct sn, int classIdx, String paramName) {
        if (classIdx < 0 || classIdx >= sn.nclasses) {
            return paramName + "=" + classIdx + " is out of bounds [0, " + (sn.nclasses - 1) + "]";
        }
        return null;
    }

    public static String snValidateNodeIndex(NetworkStruct sn, int nodeIdx) {
        return snValidateNodeIndex(sn, nodeIdx, "nodeIdx");
    }

    public static String snValidateNodeIndex(NetworkStruct sn, int nodeIdx, String paramName) {
        if (nodeIdx < 0 || nodeIdx >= sn.nnodes) {
            return paramName + "=" + nodeIdx + " is out of bounds [0, " + (sn.nnodes - 1) + "]";
        }
        return null;
    }

    public static String snValidateNodeType(NetworkStruct sn, int nodeIdx, NodeType expectedType) {
        if (nodeIdx < 0 || nodeIdx >= sn.nodetype.size()) {
            return "nodeIdx=" + nodeIdx + " is out of bounds";
        }
        NodeType actualType = sn.nodetype.get(nodeIdx);
        if (actualType != expectedType) {
            return "Node " + nodeIdx + " is " + actualType + ", expected " + expectedType;
        }
        return null;
    }

    private static double getPopulation(Matrix njobs, int classIdx) {
        if (njobs.getNumRows() == 1) {
            return njobs.get(0, classIdx);
        } else {
            return njobs.get(classIdx, 0);
        }
    }
}
