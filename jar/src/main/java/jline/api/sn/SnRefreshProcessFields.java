/**
 * @file Process Fields Refresh for NetworkStruct
 *
 * Provides functions to refresh process-related fields (mu, phi, proc, pie, phases)
 * in a NetworkStruct based on the current rate and SCV values. This allows updating
 * derived process representations after modifying service rates directly.
 *
 * Mirrors MATLAB implementation patterns for process parameter computation.
 *
 * @since LINE 3.0
 */
package jline.api.sn;

import java.util.HashMap;

import jline.api.mam.Map_erlang;
import jline.api.mam.Map_exponential;
import jline.api.mam.Map_hyperexp;
import jline.api.mam.Map_pie;
import jline.lang.JobClass;
import jline.lang.NetworkStruct;
import jline.lang.nodes.Station;
import jline.lang.constant.ProcessType;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class SnRefreshProcessFields {
    private SnRefreshProcessFields() {}

    /**
     * Refreshes process fields for a specific station-class pair based on current rate and SCV.
     */
    public static NetworkStruct snRefreshProcessFields(NetworkStruct sn, int stationIdx, int classIdx) {
        // Get rate and SCV
        if (sn.rates == null) return sn;
        double rate = sn.rates.get(stationIdx, classIdx);
        double scv = (sn.scv != null) ? sn.scv.get(stationIdx, classIdx) : 1.0;

        // Skip if rate is invalid
        if (Double.isNaN(rate) || rate <= 0 || !Double.isFinite(rate)) {
            return sn;
        }

        // Get station and job class objects
        if (stationIdx < 0 || stationIdx >= sn.stations.size()) return sn;
        Station station = sn.stations.get(stationIdx);
        if (station == null) return sn;
        if (classIdx < 0 || classIdx >= sn.jobclasses.size()) return sn;
        JobClass jobClass = sn.jobclasses.get(classIdx);
        if (jobClass == null) return sn;

        double mean = 1.0 / rate;

        // Determine process type and create MAP based on SCV
        MatrixCell map;
        ProcessType processType;
        int nPhases;
        if (Double.isNaN(scv) || Math.abs(scv - 1.0) < 1e-10) {
            // Exponential
            map = Map_exponential.map_exponential(mean);
            processType = ProcessType.EXP;
            nPhases = 1;
        } else if (scv < 1.0) {
            // Erlang: k = ceil(1/scv)
            int k = Math.max(1, (int) Math.ceil(1.0 / scv));
            map = Map_erlang.map_erlang(mean, k);
            processType = ProcessType.ERLANG;
            nPhases = k;
        } else {
            // Hyperexponential (scv > 1)
            MatrixCell hyperMap = Map_hyperexp.map_hyperexp(mean, scv, 0.99);
            if (hyperMap != null) {
                map = hyperMap;
                processType = ProcessType.HYPEREXP;
                nPhases = 2;
            } else {
                // Fallback to exponential if hyperexp fails
                map = Map_exponential.map_exponential(mean);
                processType = ProcessType.EXP;
                nPhases = 1;
            }
        }

        // Update process fields
        updateProcessFields(sn, station, jobClass, stationIdx, classIdx, map, processType, nPhases);

        return sn;
    }

    /**
     * Refreshes process fields for all station-class pairs.
     */
    public static NetworkStruct snRefreshAllProcessFields(NetworkStruct sn) {
        for (int i = 0; i < sn.nstations; i++) {
            for (int j = 0; j < sn.nclasses; j++) {
                snRefreshProcessFields(sn, i, j);
            }
        }
        return sn;
    }

    /**
     * Updates all process-related fields in NetworkStruct for a given MAP.
     */
    private static void updateProcessFields(NetworkStruct sn, Station station, JobClass jobClass,
                                            int stationIdx, int classIdx, MatrixCell map,
                                            ProcessType processType, int nPhases) {
        Matrix d0 = map.get(0);
        Matrix d1 = map.get(1);

        // Initialize maps if null
        if (sn.proc == null) {
            sn.proc = new HashMap<jline.lang.nodes.Station, java.util.Map<jline.lang.JobClass, MatrixCell>>();
        }
        if (sn.procid == null) {
            sn.procid = new HashMap<jline.lang.nodes.Station, java.util.Map<jline.lang.JobClass, ProcessType>>();
        }
        if (sn.mu == null) {
            sn.mu = new HashMap<jline.lang.nodes.Station, java.util.Map<jline.lang.JobClass, Matrix>>();
        }
        if (sn.phi == null) {
            sn.phi = new HashMap<jline.lang.nodes.Station, java.util.Map<jline.lang.JobClass, Matrix>>();
        }
        if (sn.pie == null) {
            sn.pie = new HashMap<jline.lang.nodes.Station, java.util.Map<jline.lang.JobClass, Matrix>>();
        }

        // Initialize inner maps if null
        if (sn.proc.get(station) == null) {
            sn.proc.put(station, new HashMap<jline.lang.JobClass, MatrixCell>());
        }
        if (sn.procid.get(station) == null) {
            sn.procid.put(station, new HashMap<jline.lang.JobClass, ProcessType>());
        }
        if (sn.mu.get(station) == null) {
            sn.mu.put(station, new HashMap<jline.lang.JobClass, Matrix>());
        }
        if (sn.phi.get(station) == null) {
            sn.phi.put(station, new HashMap<jline.lang.JobClass, Matrix>());
        }
        if (sn.pie.get(station) == null) {
            sn.pie.put(station, new HashMap<jline.lang.JobClass, Matrix>());
        }

        // Update process representation
        sn.proc.get(station).put(jobClass, map);
        sn.procid.get(station).put(jobClass, processType);

        // Update phases
        if (sn.phases != null) {
            sn.phases.set(stationIdx, classIdx, (double) nPhases);
        }

        // Update phasessz
        if (sn.phasessz != null) {
            sn.phasessz.set(stationIdx, classIdx, (double) Math.max(nPhases, 1));
        }

        // Recompute phaseshift for this station (cumulative sum across classes)
        if (sn.phaseshift != null) {
            double cumSum = 0.0;
            sn.phaseshift.set(stationIdx, 0, 0.0);
            for (int c = 0; c < sn.nclasses; c++) {
                double sz = (sn.phasessz != null) ? sn.phasessz.get(stationIdx, c) : 1.0;
                cumSum += sz;
                if (c + 1 < sn.phaseshift.getNumCols()) {
                    sn.phaseshift.set(stationIdx, c + 1, cumSum);
                }
            }
        }

        // Update mu (rates from -diag(D0))
        Matrix muMatrix = new Matrix(nPhases, 1);
        for (int i = 0; i < nPhases; i++) {
            muMatrix.set(i, 0, -d0.get(i, i));
        }
        sn.mu.get(station).put(jobClass, muMatrix);

        // Update phi (completion probabilities: sum(D1,2) / -diag(D0))
        Matrix phiMatrix = new Matrix(nPhases, 1);
        for (int i = 0; i < nPhases; i++) {
            double d1RowSum = 0.0;
            for (int j = 0; j < d1.getNumCols(); j++) {
                d1RowSum += d1.get(i, j);
            }
            double d0Diag = -d0.get(i, i);
            phiMatrix.set(i, 0, (d0Diag != 0.0) ? d1RowSum / d0Diag : 0.0);
        }
        sn.phi.get(station).put(jobClass, phiMatrix);

        // Update pie (initial phase distribution)
        Matrix pieMatrix = Map_pie.map_pie(map);
        sn.pie.get(station).put(jobClass, pieMatrix);
    }
}
