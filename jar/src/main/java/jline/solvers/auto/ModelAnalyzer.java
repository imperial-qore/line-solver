package jline.solvers.auto;

import jline.api.sn.*;
import jline.lang.Network;
import jline.lang.NetworkStruct;
import jline.lang.constant.SchedStrategy;
import jline.util.matrix.Matrix;

/**
 * Helper class to analyze model characteristics for solver selection
 */
public class ModelAnalyzer {

    private final Network model;
    private final NetworkStruct sn;

    public ModelAnalyzer(Network model) {
        this.model = model;
        this.sn = model.getStruct(true);
    }

    /**
     * Get average jobs per chain
     */
    public double getAvgJobsPerChain() {
        int totalJobs = getTotalJobs();
        if (sn.nchains > 0) {
            return (double) totalJobs / sn.nchains;
        }
        return 0;
    }

    /**
     * Get number of chains
     */
    public int getNumChains() {
        return sn.nchains;
    }

    /**
     * Get number of classes
     */
    public int getNumClasses() {
        return sn.nclasses;
    }

    /**
     * Get number of stations (excluding reference stations)
     */
    public int getNumStations() {
        return sn.nstations;
    }

    /**
     * Get total number of jobs in closed classes
     */
    public int getTotalJobs() {
        int totalJobs = 0;
        Matrix N = sn.njobs;
        if (N != null) {
            for (int i = 0; i < N.getNumRows(); i++) {
                for (int j = 0; j < N.getNumCols(); j++) {
                    // Open classes carry njobs=Inf, which casts to MAX_VALUE
                    double nij = N.get(i, j);
                    if (!Double.isInfinite(nij) && !Double.isNaN(nij)) {
                        totalJobs += (int) nij;
                    }
                }
            }
        }
        return totalJobs;
    }

    /**
     * Check if model has class switching
     */
    public boolean hasClassSwitching() {
        return SnHasClassSwitching.snHasClassSwitching(sn);
    }

    /**
     * Check if model has closed classes
     */
    public boolean hasClosedClasses() {
        return SnHasClosedClasses.snHasClosedClasses(sn);
    }

    /**
     * Check if model has FCFS scheduling
     */
    public boolean hasFCFS() {
        return SnHasFCFS.snHasFCFS(sn);
    }

    /**
     * Check if model has fork-join
     */
    public boolean hasForkJoin() {
        return SnHasForkJoin.snHasForkJoin(sn);
    }

    /**
     * Check if model has load-dependent stations
     */
    public boolean hasLoadDependence() {
        return SnHasLoadDependence.snHasLoadDependence(sn);
    }

    /**
     * Check if model has multiple chains
     */
    public boolean hasMultiChain() {
        return sn.nchains > 1;
    }

    /**
     * Check if model has multi-server stations
     */
    public boolean hasMultiServer() {
        return SnHasMultiServer.snHasMultiServer(sn);
    }

    /**
     * Check if all stations are infinite servers
     */
    public boolean hasOnlyInfiniteServers() {
        for (int i = 0; i < sn.nstations; i++) {
            if (!sn.sched.get(this.model.getStations().get(i)).equals(SchedStrategy.INF)) {
                return false;
            }
        }
        return true;
    }

    /**
     * Check if model has open classes
     */
    public boolean hasOpenClasses() {
        return SnHasOpenClasses.snHasOpenClasses(sn);
    }

    /**
     * Check if model has PS or PSPRIO scheduling
     */
    public boolean hasPSorPSPRIO() {
        return SnHasPS.snHasPS(sn) || SnHasPSPRIO.snHasPSPRIO(sn);
    }

    /**
     * Check if model has priorities
     */
    public boolean hasPriorities() {
        return SnHasPriorities.snHasPriorities(sn);
    }

    /**
     * Check if model has product form
     */
    public boolean hasProductForm() {
        return SnHasProductForm.snHasProductForm(sn);
    }

    /**
     * Check if model has single chain
     */
    public boolean hasSingleChain() {
        return sn.nchains == 1;
    }

    /**
     * Check if model is closed (only closed classes)
     */
    public boolean isClosedModel() {
        return hasClosedClasses() && !hasOpenClasses();
    }

    /**
     * Check if model has homogeneous scheduling (all stations use the given strategy)
     */
    public boolean hasHomogeneousScheduling(SchedStrategy strategy) {
        return model.hasHomogeneousScheduling(strategy);
    }

    /**
     * Check if model is open (has open classes and no closed classes)
     */
    public boolean isOpenModel() {
        return hasOpenClasses() && !hasClosedClasses();
    }

    /**
     * Check if model is mixed (has both open and closed classes)
     */
    public boolean isMixedModel() {
        return hasOpenClasses() && hasClosedClasses();
    }
}