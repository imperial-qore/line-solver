/**
 * @file Retrial queue detection
 *
 * @since LINE 3.1.0
 */
package jline.api.qsys;

import java.util.ArrayList;
import java.util.Map;

import jline.api.sn.SnIsOpenModel;
import jline.lang.NetworkStruct;
import jline.lang.constant.DropStrategy;
import jline.lang.constant.NodeType;
import jline.util.matrix.Matrix;

public final class Qsys_is_retrial {
    private Qsys_is_retrial() {}

    /**
     * Checks if network is a valid BMAP/PH/N/N bufferless retrial queue.
     */
    public static RetrialInfo qsys_is_retrial(NetworkStruct sn) {
        RetrialInfo retInfo = new RetrialInfo();

        if (!SnIsOpenModel.snIsOpenModel(sn)) {
            retInfo.setErrorMsg("BMAP/PH/N/N retrial solver requires open queueing model.");
            return retInfo;
        }

        if (sn.nclasses > 1) {
            retInfo.setErrorMsg("BMAP/PH/N/N retrial solver currently supports single class only.");
            return retInfo;
        }

        retInfo.setClassIdx(0);

        ArrayList<Integer> bufferlessStations = new ArrayList<Integer>();
        if (sn.stationToNode == null || sn.cap == null || sn.nservers == null) {
            retInfo.setErrorMsg("Missing network structure fields.");
            return retInfo;
        }
        for (int ist = 0; ist < sn.nstations; ist++) {
            if (ist >= sn.stationToNode.length()) continue;
            int nodeIdx = (int) sn.stationToNode.get(ist);
            if (nodeIdx >= 0 && nodeIdx < sn.nodetype.size() && sn.nodetype.get(nodeIdx) == NodeType.Queue) {
                if (ist >= sn.cap.length() || ist >= sn.nservers.length()) continue;
                // see _kb/03-api-layer.md for rationale
                boolean hasRetrialProc = false;
                if (sn.retrialProc != null && ist < sn.stations.size()) {
                    Map<?, ?> procs = sn.retrialProc.get(sn.stations.get(ist));
                    if (procs != null) {
                        for (Object proc : procs.values()) {
                            if (proc != null) {
                                hasRetrialProc = true;
                                break;
                            }
                        }
                    }
                }
                double cap = sn.cap.get(ist);
                double nservers = sn.nservers.get(ist);
                boolean isBufferless = Double.isFinite(cap) && Math.abs(cap - nservers) < 1e-10;
                if (hasRetrialProc || isBufferless) {
                    bufferlessStations.add(ist);
                }
            }
        }

        if (bufferlessStations.isEmpty()) {
            retInfo.setErrorMsg("No retrial queue found (configure setOrbit, or a bufferless station with setRetrial).");
            return retInfo;
        }

        int retrialStation = -1;
        for (int ist : bufferlessStations) {
            if (sn.droprule != null && ist < sn.stations.size()) {
                Object station = sn.stations.get(ist);
                Map<?, ?> classDropRules = (Map<?, ?>) sn.droprule.get(station);
                if (classDropRules != null) {
                    for (Object dropEntry : classDropRules.values()) {
                        if (dropEntry == DropStrategy.Retrial || dropEntry == DropStrategy.RetrialWithLimit) {
                            retrialStation = ist;
                            break;
                        }
                    }
                }
            }
            if (retrialStation >= 0) break;
        }

        if (retrialStation < 0) {
            retInfo.setErrorMsg("No retrial drop strategy configured on bufferless queue.");
            return retInfo;
        }

        retInfo.setStationIdx(retrialStation);
        retInfo.setNodeIdx((int) sn.stationToNode.get(retrialStation));
        retInfo.setN((int) sn.nservers.get(retrialStation));

        int sourceStation = -1;
        for (int ist = 0; ist < sn.nstations; ist++) {
            int nodeIdx = (int) sn.stationToNode.get(ist);
            if (nodeIdx >= 0 && nodeIdx < sn.nodetype.size() && sn.nodetype.get(nodeIdx) == NodeType.Source) {
                sourceStation = ist;
                break;
            }
        }

        if (sourceStation < 0) {
            retInfo.setErrorMsg("No Source node found.");
            return retInfo;
        }

        retInfo.setSourceIdx(sourceStation);
        retInfo.setAlpha(0.1);
        retInfo.setR(retInfo.getN() - 1);

        if (sn.nregions > 0 && sn.region != null) {
            for (int f = 0; f < sn.nregions; f++) {
                Matrix regionMatrix = sn.region.get(f);
                if (regionMatrix != null && regionMatrix.getNumRows() > retrialStation) {
                    int globalCapCol = regionMatrix.getNumCols() - 1;
                    double R_fcr = regionMatrix.get(retrialStation, globalCapCol);
                    if (R_fcr > 0 && R_fcr <= (double) (retInfo.getN() - 1)) {
                        retInfo.setR((int) R_fcr);
                        break;
                    }
                }
            }
        }

        retInfo.setRetrial(true);
        return retInfo;
    }
}
