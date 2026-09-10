package jline.api.fes;

import java.util.Arrays;
import java.util.List;
import java.util.Objects;

import jline.lang.Network;
import jline.lang.nodes.Station;
import jline.util.matrix.Matrix;

/**
 * Information needed to deaggregate FES results back to original model.
 */
public final class FESDeaggInfo {
    public final Network originalModel;
    public final List<Station> stationSubset;
    public final int[] subsetIndices;
    public final int[] complementIndices;
    public final List<Matrix> throughputTable;
    public final Matrix cutoffs;
    public final Matrix stochCompSubset;
    public final Matrix stochCompComplement;
    public final Network isolatedModel;
    public final int fesNodeIdx;

    public FESDeaggInfo(
            Network originalModel,
            List<Station> stationSubset,
            int[] subsetIndices,
            int[] complementIndices,
            List<Matrix> throughputTable,
            Matrix cutoffs,
            Matrix stochCompSubset,
            Matrix stochCompComplement,
            Network isolatedModel,
            int fesNodeIdx) {
        this.originalModel = originalModel;
        this.stationSubset = stationSubset;
        this.subsetIndices = subsetIndices;
        this.complementIndices = complementIndices;
        this.throughputTable = throughputTable;
        this.cutoffs = cutoffs;
        this.stochCompSubset = stochCompSubset;
        this.stochCompComplement = stochCompComplement;
        this.isolatedModel = isolatedModel;
        this.fesNodeIdx = fesNodeIdx;
    }

    public Network getOriginalModel() { return originalModel; }
    public List<Station> getStationSubset() { return stationSubset; }
    public int[] getSubsetIndices() { return subsetIndices; }
    public int[] getComplementIndices() { return complementIndices; }
    public List<Matrix> getThroughputTable() { return throughputTable; }
    public Matrix getCutoffs() { return cutoffs; }
    public Matrix getStochCompSubset() { return stochCompSubset; }
    public Matrix getStochCompComplement() { return stochCompComplement; }
    public Network getIsolatedModel() { return isolatedModel; }
    public int getFesNodeIdx() { return fesNodeIdx; }

    @Override
    public boolean equals(Object other) {
        if (this == other) return true;
        if (!(other instanceof FESDeaggInfo)) return false;
        FESDeaggInfo o = (FESDeaggInfo) other;
        return Objects.equals(originalModel, o.originalModel)
                && Objects.equals(stationSubset, o.stationSubset)
                && Arrays.equals(subsetIndices, o.subsetIndices)
                && Arrays.equals(complementIndices, o.complementIndices)
                && Objects.equals(throughputTable, o.throughputTable)
                && Objects.equals(cutoffs, o.cutoffs)
                && Objects.equals(stochCompSubset, o.stochCompSubset)
                && Objects.equals(stochCompComplement, o.stochCompComplement)
                && Objects.equals(isolatedModel, o.isolatedModel)
                && fesNodeIdx == o.fesNodeIdx;
    }

    @Override
    public int hashCode() {
        int result = originalModel != null ? originalModel.hashCode() : 0;
        result = 31 * result + (stationSubset != null ? stationSubset.hashCode() : 0);
        result = 31 * result + Arrays.hashCode(subsetIndices);
        result = 31 * result + Arrays.hashCode(complementIndices);
        result = 31 * result + (throughputTable != null ? throughputTable.hashCode() : 0);
        result = 31 * result + (cutoffs != null ? cutoffs.hashCode() : 0);
        result = 31 * result + (stochCompSubset != null ? stochCompSubset.hashCode() : 0);
        result = 31 * result + (stochCompComplement != null ? stochCompComplement.hashCode() : 0);
        result = 31 * result + (isolatedModel != null ? isolatedModel.hashCode() : 0);
        result = 31 * result + fesNodeIdx;
        return result;
    }
}
