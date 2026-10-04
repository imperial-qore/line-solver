package jline.api.fes;

import jline.lang.Network;
import jline.lang.nodes.Queue;

/**
 * Result of Flow-Equivalent Server (FES) aggregation.
 *
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 */
public final class FESResult {
    /** New Network with FES replacing the subset. */
    public final Network fesModel;

    /** Reference to the FES Queue station. */
    public final Queue fesStation;

    /** Deaggregation information. */
    public final FESDeaggInfo deaggInfo;

    public FESResult(Network fesModel, Queue fesStation, FESDeaggInfo deaggInfo) {
        this.fesModel = fesModel;
        this.fesStation = fesStation;
        this.deaggInfo = deaggInfo;
    }

    public Network getFesModel() { return fesModel; }
    public Queue getFesStation() { return fesStation; }
    public FESDeaggInfo getDeaggInfo() { return deaggInfo; }

    @Override
    public boolean equals(Object other) {
        if (this == other) return true;
        if (!(other instanceof FESResult)) return false;
        FESResult o = (FESResult) other;
        return java.util.Objects.equals(fesModel, o.fesModel)
                && java.util.Objects.equals(fesStation, o.fesStation)
                && java.util.Objects.equals(deaggInfo, o.deaggInfo);
    }

    @Override
    public int hashCode() {
        return java.util.Objects.hash(fesModel, fesStation, deaggInfo);
    }

    @Override
    public String toString() {
        return "FESResult{fesModel=" + fesModel + ", fesStation=" + fesStation + ", deaggInfo=" + deaggInfo + "}";
    }
}
