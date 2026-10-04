/**
 * @file Service Rate Modification for NetworkStruct
 *
 * @since LINE 3.0
 */
package jline.api.sn;

import java.util.ArrayList;
import java.util.List;

import jline.lang.NetworkStruct;
import jline.util.Pair;
import jline.util.matrix.Matrix;

public final class SnSetService {
    private SnSetService() {}

    /**
     * Sets the service rate at a specific station and class.
     */
    public static NetworkStruct snSetService(NetworkStruct sn, int stationIdx, int classIdx,
                                             double rate, double scv, ModifyMode mode,
                                             ValidationLevel validation, boolean autoRefresh) {
        if (validation != ValidationLevel.NONE) {
            List<String> errors = new ArrayList<String>();

            String err1 = SnValidate.snValidateStationIndex(sn, stationIdx, "stationIdx");
            if (err1 != null) errors.add(err1);
            String err2 = SnValidate.snValidateClassIndex(sn, classIdx, "classIdx");
            if (err2 != null) errors.add(err2);

            if (validation == ValidationLevel.FULL) {
                if (Double.isNaN(rate)) {
                    errors.add("rate is NaN");
                } else if (rate < 0) {
                    errors.add("rate=" + rate + " must be non-negative");
                }
                if (Double.isNaN(scv)) {
                    errors.add("scv is NaN");
                } else if (scv < 0) {
                    errors.add("scv=" + scv + " must be non-negative");
                }
            }

            if (!errors.isEmpty()) {
                throw new SnValidationException("snSetService validation failed", errors);
            }
        }

        NetworkStruct snWork = (mode == ModifyMode.COPY) ? sn.<NetworkStruct>copy() : sn;

        if (snWork.rates != null) {
            snWork.rates.set(stationIdx, classIdx, rate);
        }

        if (snWork.scv != null) {
            snWork.scv.set(stationIdx, classIdx, scv);
        }

        if (autoRefresh) {
            SnRefreshProcessFields.snRefreshProcessFields(snWork, stationIdx, classIdx);
        }

        return snWork;
    }

    public static NetworkStruct snSetService(NetworkStruct sn, int stationIdx, int classIdx,
                                             double rate, double scv, ModifyMode mode,
                                             ValidationLevel validation) {
        return snSetService(sn, stationIdx, classIdx, rate, scv, mode, validation, false);
    }

    public static NetworkStruct snSetService(NetworkStruct sn, int stationIdx, int classIdx,
                                             double rate, double scv, ModifyMode mode) {
        return snSetService(sn, stationIdx, classIdx, rate, scv, mode, ValidationLevel.MINIMAL, false);
    }

    public static NetworkStruct snSetService(NetworkStruct sn, int stationIdx, int classIdx,
                                             double rate, double scv) {
        return snSetService(sn, stationIdx, classIdx, rate, scv, ModifyMode.IN_PLACE,
                ValidationLevel.MINIMAL, false);
    }

    public static NetworkStruct snSetService(NetworkStruct sn, int stationIdx, int classIdx, double rate) {
        return snSetService(sn, stationIdx, classIdx, rate, 1.0, ModifyMode.IN_PLACE,
                ValidationLevel.MINIMAL, false);
    }

    /**
     * Sets service rates for multiple station-class pairs in a single operation.
     */
    public static NetworkStruct snSetServiceBatch(NetworkStruct sn, Matrix rates, Matrix scvs,
                                                  ModifyMode mode, ValidationLevel validation,
                                                  boolean autoRefresh) {
        if (validation != ValidationLevel.NONE) {
            List<String> errors = new ArrayList<String>();

            if (rates.getNumRows() != sn.nstations || rates.getNumCols() != sn.nclasses) {
                errors.add("rates matrix dimensions (" + rates.getNumRows() + "x" + rates.getNumCols()
                        + ") do not match (" + sn.nstations + "x" + sn.nclasses + ")");
            }

            if (scvs != null) {
                if (scvs.getNumRows() != sn.nstations || scvs.getNumCols() != sn.nclasses) {
                    errors.add("scvs matrix dimensions (" + scvs.getNumRows() + "x" + scvs.getNumCols()
                            + ") do not match (" + sn.nstations + "x" + sn.nclasses + ")");
                }
            }

            if (validation == ValidationLevel.FULL) {
                for (int i = 0; i < rates.getNumRows(); i++) {
                    for (int j = 0; j < rates.getNumCols(); j++) {
                        double r = rates.get(i, j);
                        if (!Double.isNaN(r) && r < 0) {
                            errors.add("rates[" + i + "," + j + "]=" + r + " must be non-negative");
                        }
                        if (scvs != null) {
                            double s = scvs.get(i, j);
                            if (!Double.isNaN(s) && s < 0) {
                                errors.add("scvs[" + i + "," + j + "]=" + s + " must be non-negative");
                            }
                        }
                    }
                }
            }

            if (!errors.isEmpty()) {
                throw new SnValidationException("snSetServiceBatch validation failed", errors);
            }
        }

        NetworkStruct snWork = (mode == ModifyMode.COPY) ? sn.<NetworkStruct>copy() : sn;

        List<Pair<Integer, Integer>> updatedPairs = new ArrayList<Pair<Integer, Integer>>();

        for (int i = 0; i < rates.getNumRows(); i++) {
            for (int j = 0; j < rates.getNumCols(); j++) {
                double r = rates.get(i, j);
                if (!Double.isNaN(r)) {
                    if (snWork.rates != null) {
                        snWork.rates.set(i, j, r);
                    }
                    updatedPairs.add(new Pair<Integer, Integer>(i, j));
                }
            }
        }

        if (scvs != null && snWork.scv != null) {
            for (int i = 0; i < scvs.getNumRows(); i++) {
                for (int j = 0; j < scvs.getNumCols(); j++) {
                    double s = scvs.get(i, j);
                    if (!Double.isNaN(s)) {
                        snWork.scv.set(i, j, s);
                    }
                }
            }
        }

        if (autoRefresh) {
            for (Pair<Integer, Integer> p : updatedPairs) {
                SnRefreshProcessFields.snRefreshProcessFields(snWork, p.getLeft(), p.getRight());
            }
        }

        return snWork;
    }

    public static NetworkStruct snSetServiceBatch(NetworkStruct sn, Matrix rates, Matrix scvs,
                                                  ModifyMode mode, ValidationLevel validation) {
        return snSetServiceBatch(sn, rates, scvs, mode, validation, false);
    }

    public static NetworkStruct snSetServiceBatch(NetworkStruct sn, Matrix rates, Matrix scvs, ModifyMode mode) {
        return snSetServiceBatch(sn, rates, scvs, mode, ValidationLevel.MINIMAL, false);
    }

    public static NetworkStruct snSetServiceBatch(NetworkStruct sn, Matrix rates, Matrix scvs) {
        return snSetServiceBatch(sn, rates, scvs, ModifyMode.IN_PLACE, ValidationLevel.MINIMAL, false);
    }

    public static NetworkStruct snSetServiceBatch(NetworkStruct sn, Matrix rates) {
        return snSetServiceBatch(sn, rates, null, ModifyMode.IN_PLACE, ValidationLevel.MINIMAL, false);
    }
}
