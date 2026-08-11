/**
 * @file AoI parameter extraction
 *
 * @since LINE 3.0
 */
package jline.api.aoi;

import jline.lang.NetworkStruct;
import jline.lang.constant.SchedStrategy;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Aoi_extract_params {
    private Aoi_extract_params() {}

    public static AoiParams aoi_extract_params(NetworkStruct sn, AoiValidationResult aoiInfo) {
        return aoi_extract_params(sn, aoiInfo, Double.NaN);
    }

    /**
     * Extract parameters from LINE network for AoI analysis.
     */
    public static AoiParams aoi_extract_params(NetworkStruct sn, AoiValidationResult aoiInfo, double aoiPreemption) {
        if (!aoiInfo.isAoI()) {
            throw new IllegalArgumentException("Network is not a valid AoI topology: " + aoiInfo.getErrorMsg());
        }

        int sourceStation = aoiInfo.getSourceStation();
        int queueStation = aoiInfo.getQueueStation();
        int capacity = aoiInfo.getCapacity();

        // Find the open class
        int classIdx = -1;
        for (int r = 0; r < sn.nclasses; r++) {
            if (Double.isInfinite(sn.njobs.get(r))) {
                classIdx = r;
                break;
            }
        }
        if (classIdx < 0) {
            throw new IllegalArgumentException("No open class found");
        }

        Object queueStationObj = sn.stations.get(queueStation);
        Object jobClassObj = sn.jobclasses.get(classIdx);
        java.util.Map<?, MatrixCell> queueProcMap = sn.proc.get(queueStationObj);
        MatrixCell serviceProc = queueProcMap == null ? null : queueProcMap.get(jobClassObj);

        Matrix sigma;
        Matrix S;

        if (serviceProc == null || serviceProc.size() < 2 || serviceProc.get(0).hasNaN()) {
            double mu = sn.rates.get(queueStation, classIdx);
            sigma = new Matrix(1, 1);
            sigma.set(0, 0, 1.0);
            S = new Matrix(1, 1);
            S.set(0, 0, -mu);
        } else {
            Aoi_dist2phResult result = Aoi_dist2ph.aoi_dist2ph(serviceProc);
            sigma = result.getAlpha();
            S = result.getT();
        }

        SchedStrategy schedStrategy = (SchedStrategy) sn.sched.get(queueStationObj);

        if (capacity == 1) {
            // BUFFERLESS
            Object sourceStationObj = sn.stations.get(sourceStation);
            java.util.Map<?, MatrixCell> sourceProcMap = sn.proc.get(sourceStationObj);
            MatrixCell arrivalProc = sourceProcMap == null ? null : sourceProcMap.get(jobClassObj);

            Matrix tau;
            Matrix T;

            if (arrivalProc == null || arrivalProc.size() < 2 || arrivalProc.get(0).hasNaN()) {
                double lambda = sn.rates.get(sourceStation, classIdx);
                tau = new Matrix(1, 1);
                tau.set(0, 0, 1.0);
                T = new Matrix(1, 1);
                T.set(0, 0, -lambda);
            } else {
                Aoi_dist2phResult result = Aoi_dist2ph.aoi_dist2ph(arrivalProc);
                tau = result.getAlpha();
                T = result.getT();
            }

            double p;
            if (!Double.isNaN(aoiPreemption)) {
                p = aoiPreemption;
            } else {
                if (schedStrategy == SchedStrategy.FCFS) p = 0.0;
                else if (schedStrategy == SchedStrategy.LCFS) p = 0.0;
                else if (schedStrategy == SchedStrategy.LCFSPR) p = 1.0;
                else p = 0.0;
            }

            return new AoiParams(tau, T, sigma, S, p, Double.NaN, Double.NaN, "bufferless", "PH");
        } else {
            // SINGLE-BUFFER
            double lambda = sn.rates.get(sourceStation, classIdx);

            double r;
            if (!Double.isNaN(aoiPreemption)) {
                r = aoiPreemption;
            } else {
                if (schedStrategy == SchedStrategy.FCFS) r = 0.0;
                else if (schedStrategy == SchedStrategy.LCFS) r = 1.0;
                else if (schedStrategy == SchedStrategy.LCFSPR) r = 1.0;
                else r = 0.0;
            }

            return new AoiParams(null, null, sigma, S, Double.NaN, lambda, r, "singlebuffer", "M");
        }
    }
}
