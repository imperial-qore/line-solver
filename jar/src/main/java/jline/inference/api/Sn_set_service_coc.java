/*
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 */
package jline.inference.api;

import java.util.HashMap;
import java.util.Map;

import jline.api.mam.Map_erlang;
import jline.api.mam.Map_exponential;
import jline.api.mam.Map_hyperexp;
import jline.api.mam.Map_pie;
import jline.lang.JobClass;
import jline.lang.NetworkStruct;
import jline.lang.constant.ProcessType;
import jline.lang.nodes.Station;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Sn_set_service_coc {
    private Sn_set_service_coc() {}

    /**
     * Update service rate in NetworkStruct preserving Map-based format.
     */
    public static void sn_set_service_coc(NetworkStruct sn, int stationIdx, int classIdx,
                                          double rate, double scv) {
        Station station = sn.stations.get(stationIdx);
        JobClass jobclass = sn.jobclasses.get(classIdx);

        sn.rates.set(stationIdx, classIdx, rate);
        sn.scv.set(stationIdx, classIdx, scv);

        if (Double.isNaN(rate) || rate <= 0 || Double.isInfinite(rate)) return;

        double meanVal = 1.0 / rate;

        MatrixCell MAP;
        int nPhases;
        ProcessType procType;

        if (Double.isNaN(scv) || Math.abs(scv - 1.0) < 1e-10) {
            MAP = Map_exponential.map_exponential(meanVal);
            nPhases = 1;
            procType = ProcessType.EXP;
        } else if (scv < 1.0) {
            int k = Math.max(1, (int) Math.ceil(1.0 / scv));
            MAP = Map_erlang.map_erlang(meanVal, k);
            nPhases = k;
            procType = ProcessType.ERLANG;
        } else {
            MatrixCell hyperMAP = Map_hyperexp.map_hyperexp(meanVal, scv, 0.99);
            if (hyperMAP != null) {
                MAP = hyperMAP;
                nPhases = 2;
                procType = ProcessType.HYPEREXP;
            } else {
                MAP = Map_exponential.map_exponential(meanVal);
                nPhases = 1;
                procType = ProcessType.EXP;
            }
        }

        Matrix D0 = MAP.get(0);
        Matrix D1 = MAP.get(1);

        Map<JobClass, MatrixCell> stationProcMap = sn.proc.get(station);
        if (stationProcMap == null) {
            stationProcMap = new HashMap<JobClass, MatrixCell>();
            sn.proc.put(station, stationProcMap);
        }
        stationProcMap.put(jobclass, MAP);

        Map<JobClass, ProcessType> stationProcIdMap = sn.procid.get(station);
        if (stationProcIdMap == null) {
            stationProcIdMap = new HashMap<JobClass, ProcessType>();
            sn.procid.put(station, stationProcIdMap);
        }
        stationProcIdMap.put(jobclass, procType);

        sn.phases.set(stationIdx, classIdx, (double) nPhases);
        sn.phasessz.set(stationIdx, classIdx, (double) Math.max(nPhases, 1));

        double cumSum = 0.0;
        sn.phaseshift.set(stationIdx, 0, 0.0);
        for (int c = 0; c < sn.nclasses; c++) {
            cumSum += sn.phasessz.get(stationIdx, c);
            if (c + 1 < sn.phaseshift.getNumCols()) {
                sn.phaseshift.set(stationIdx, c + 1, cumSum);
            }
        }

        Matrix muVec = new Matrix(nPhases, 1);
        for (int i = 0; i < nPhases; i++) {
            muVec.set(i, 0, -D0.get(i, i));
        }
        Map<JobClass, Matrix> stationMuMap = sn.mu.get(station);
        if (stationMuMap == null) {
            stationMuMap = new HashMap<JobClass, Matrix>();
            sn.mu.put(station, stationMuMap);
        }
        stationMuMap.put(jobclass, muVec);

        Matrix phiVec = new Matrix(nPhases, 1);
        for (int i = 0; i < nPhases; i++) {
            double d1RowSum = 0.0;
            for (int j = 0; j < D1.getNumCols(); j++) {
                d1RowSum += D1.get(i, j);
            }
            double d0Diag = -D0.get(i, i);
            if (d0Diag != 0.0) {
                phiVec.set(i, 0, d1RowSum / d0Diag);
            }
        }
        Map<JobClass, Matrix> stationPhiMap = sn.phi.get(station);
        if (stationPhiMap == null) {
            stationPhiMap = new HashMap<JobClass, Matrix>();
            sn.phi.put(station, stationPhiMap);
        }
        stationPhiMap.put(jobclass, phiVec);

        Matrix pieVec = Map_pie.map_pie(MAP);
        Map<JobClass, Matrix> stationPieMap = sn.pie.get(station);
        if (stationPieMap == null) {
            stationPieMap = new HashMap<JobClass, Matrix>();
            sn.pie.put(station, stationPieMap);
        }
        stationPieMap.put(jobclass, pieVec);
    }

    public static void sn_set_service_coc(NetworkStruct sn, int stationIdx, int classIdx, double rate) {
        sn_set_service_coc(sn, stationIdx, classIdx, rate, 1.0);
    }
}
