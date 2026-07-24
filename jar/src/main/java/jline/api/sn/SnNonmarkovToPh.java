/**
 * @file Non-Markovian to Phase-Type Distribution Converter.
 *
 * @since LINE 3.0
 */
package jline.api.sn;

import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.function.DoubleUnaryOperator;

import jline.api.mam.Aph_bernstein;
import jline.api.mam.Map_erlang;
import jline.api.mam.Map_pie;
import jline.api.mam.Map_scale;
import jline.io.InputOutput;
import jline.lang.Event;
import jline.lang.JobClass;
import jline.lang.NetworkStruct;
import jline.lang.Sync;
import jline.lang.constant.EventType;
import jline.lang.constant.NodeType;
import jline.lang.constant.ProcessType;
import jline.lang.Mode;
import jline.lang.nodeparam.TransitionNodeParam;
import jline.lang.nodes.Station;
import jline.lang.nodes.StatefulNode;
import jline.lang.nodes.Transition;
import jline.lang.processes.Gamma;
import jline.lang.processes.Lognormal;
import jline.lang.processes.Weibull;
import jline.solvers.SolverOptions;
import jline.util.Pair;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class SnNonmarkovToPh {
    private SnNonmarkovToPh() {}

    private static final Set<ProcessType> SKIP_BERNSTEIN_CONVERSION = new HashSet<ProcessType>(Arrays.asList(
            ProcessType.EXP,
            ProcessType.ERLANG,
            ProcessType.HYPEREXP,
            ProcessType.PH,
            ProcessType.APH,
            ProcessType.MAP,
            ProcessType.DMAP,
            ProcessType.MMAP,
            ProcessType.BMAP,
            ProcessType.ME,
            ProcessType.RAP,
            ProcessType.COXIAN,
            ProcessType.COX2,
            ProcessType.MMPP2,
            ProcessType.IMMEDIATE,
            // see _kb/03-api-layer.md for rationale
            ProcessType.NHPP,
            ProcessType.DISABLED));

    public static NetworkStruct snNonmarkovToPh(NetworkStruct snInput, SolverOptions options) {
        // Default keeps the DET process-type tag after conversion so MAM can still
        // dispatch to the exact MAP/D/c (Crommelin) / D/M/c (Smith) solvers.
        return snNonmarkovToPh(snInput, options, true);
    }

    /**
     * @param keepDetTag when true, a converted DET keeps procid=DET (MAM needs this
     *   for exact deterministic dispatch); when false (CTMC/SSA), the converted DET
     *   is tagged MAP like every other converted distribution, so the state-space
     *   generator treats its Erlang-PH representation as a proper multi-phase process
     *   instead of a single-phase deterministic placeholder (fixes M/D/1 QLen=0).
     */
    public static NetworkStruct snNonmarkovToPh(NetworkStruct snInput, SolverOptions options, boolean keepDetTag) {
        String nonmkvMethod = (options.config != null && options.config.nonmkv != null) ? options.config.nonmkv : "bernstein";

        if (nonmkvMethod.equalsIgnoreCase("none")) {
            return snInput;
        }

        int nPhases = (options.config != null) ? options.config.nonmkvorder : 20;
        boolean preserveDet = (options.config != null) && options.config.preserveDet;

        for (int ist = 0; ist < snInput.nstations; ist++) {
            Station station = snInput.stations.get(ist);
            for (int r = 0; r < snInput.nclasses; r++) {
                JobClass jobClass = snInput.jobclasses.get(r);

                Map<JobClass, ProcessType> stProcid = snInput.procid.get(station);
                if (stProcid == null) continue;
                ProcessType procType = stProcid.get(jobClass);
                if (procType == null) continue;

                if (SKIP_BERNSTEIN_CONVERSION.contains(procType)) continue;
                if (preserveDet && procType == ProcessType.DET) continue;

                double rate = snInput.rates.get(ist, r);
                if (rate <= 0 || !Double.isFinite(rate)) continue;
                final double targetMean = 1.0 / rate;

                String distName = ProcessType.toText(procType);
                InputOutput.line_warning(
                        "snNonmarkovToPh",
                        "Distribution " + distName + " at station " + ist + " class " + r
                                + " is non-Markovian and will be converted to PH (" + nPhases + " phases).");

                Map<JobClass, MatrixCell> stProc = snInput.proc.get(station);
                if (stProc == null) continue;
                MatrixCell origProc = stProc.get(jobClass);
                if (origProc == null) continue;

                MatrixCell map;
                if (procType == ProcessType.GAMMA) {
                    final double shape = origProc.get(0).toDouble();
                    final double scale = origProc.get(1).toDouble();
                    DoubleUnaryOperator pdfFunc = new DoubleUnaryOperator() {
                        @Override
                        public double applyAsDouble(double x) {
                            if (x <= 0) return 0.0;
                            return new Gamma(shape, scale).evalPDF(x);
                        }
                    };
                    Pair<Matrix, Matrix> dd = Aph_bernstein.aph_bernstein(pdfFunc, nPhases);
                    map = Map_scale.map_scale(pairToMatrixCell(dd.getLeft(), dd.getRight()), targetMean);
                } else if (procType == ProcessType.WEIBULL) {
                    final double shapeParam = origProc.get(0).toDouble();
                    final double scaleParam = origProc.get(1).toDouble();
                    DoubleUnaryOperator pdfFunc = new DoubleUnaryOperator() {
                        @Override
                        public double applyAsDouble(double x) {
                            if (x <= 0) return 0.0;
                            return new Weibull(shapeParam, scaleParam).evalPDF(x);
                        }
                    };
                    Pair<Matrix, Matrix> dd = Aph_bernstein.aph_bernstein(pdfFunc, nPhases);
                    map = Map_scale.map_scale(pairToMatrixCell(dd.getLeft(), dd.getRight()), targetMean);
                } else if (procType == ProcessType.LOGNORMAL) {
                    final double mu = origProc.get(0).toDouble();
                    final double sigma = origProc.get(1).toDouble();
                    DoubleUnaryOperator pdfFunc = new DoubleUnaryOperator() {
                        @Override
                        public double applyAsDouble(double x) {
                            if (x <= 0) return 0.0;
                            return new Lognormal(mu, sigma).evalPDF(x);
                        }
                    };
                    Pair<Matrix, Matrix> dd = Aph_bernstein.aph_bernstein(pdfFunc, nPhases);
                    map = Map_scale.map_scale(pairToMatrixCell(dd.getLeft(), dd.getRight()), targetMean);
                } else if (procType == ProcessType.PARETO) {
                    final double shapeParam = origProc.get(0).toDouble();
                    final double scaleParam = origProc.get(1).toDouble();
                    DoubleUnaryOperator pdfFunc = new DoubleUnaryOperator() {
                        @Override
                        public double applyAsDouble(double x) {
                            if (x < scaleParam) return 0.0;
                            return shapeParam * Math.pow(scaleParam, shapeParam) / Math.pow(x, shapeParam + 1);
                        }
                    };
                    Pair<Matrix, Matrix> dd = Aph_bernstein.aph_bernstein(pdfFunc, nPhases);
                    map = Map_scale.map_scale(pairToMatrixCell(dd.getLeft(), dd.getRight()), targetMean);
                } else if (procType == ProcessType.UNIFORM) {
                    final double minVal = origProc.get(0).toDouble();
                    final double maxVal = origProc.get(1).toDouble();
                    DoubleUnaryOperator pdfFunc = new DoubleUnaryOperator() {
                        @Override
                        public double applyAsDouble(double x) {
                            if (x >= minVal && x <= maxVal) return 1.0 / (maxVal - minVal);
                            return 0.0;
                        }
                    };
                    Pair<Matrix, Matrix> dd = Aph_bernstein.aph_bernstein(pdfFunc, nPhases);
                    map = Map_scale.map_scale(pairToMatrixCell(dd.getLeft(), dd.getRight()), targetMean);
                } else if (procType == ProcessType.DET) {
                    map = Map_erlang.map_erlang(targetMean, nPhases);
                } else {
                    map = Map_erlang.map_erlang(targetMean, nPhases);
                }

                int actualPhases = map.get(0).getNumRows();
                updateSnForMAP(snInput, station, jobClass, ist, r, map, actualPhases, procType, keepDetTag);
            }
        }

        convertTransitionFiringDistributions(snInput, nPhases, preserveDet);

        return snInput;
    }

    /**
     * Walk every Transition node and convert non-Markovian firing distributions
     * (Det / Gamma / Weibull / Pareto / Uniform / Lognormal) to a phase-type
     * approximation. Mirrors the per-station loop above, and the SPN section of
     * sn_nonmarkov_toph.m:152-249. Without this an SPN firing time collapses to
     * the single exponential phase that Station.getServiceRates leaves in place,
     * so a Det or Gamma firing was solved as if it were exponential.
     */
    private static void convertTransitionFiringDistributions(NetworkStruct sn, int nPhases, boolean preserveDet) {
        if (sn.nodeparam == null || sn.nodeparam.isEmpty()) return;

        for (int ind = 0; ind < sn.nnodes; ind++) {
            if (sn.nodetype.get(ind) != NodeType.Transition) continue;
            Object rawParam = sn.nodeparam.get(sn.nodes.get(ind));
            if (!(rawParam instanceof TransitionNodeParam)) continue;
            TransitionNodeParam transParam = (TransitionNodeParam) rawParam;
            if (transParam.firingproc == null || transParam.firingprocid == null) continue;

            Transition transition = (Transition) sn.nodes.get(ind);
            for (int m = 0; m < transParam.nmodes; m++) {
                Mode modeObj = transition.getModes().get(m);

                // Markovian modes (Exp/Erlang/HyperExp/PH/APH/...) are populated with a
                // valid (D0,D1) PH and a finite firingphases by refreshPetriNetNodes.
                if (transParam.firingphases != null && transParam.firingphases.length() > m) {
                    double phasesM = transParam.firingphases.get(m);
                    if (!Double.isNaN(phasesM) && phasesM > 0) continue;
                }

                ProcessType procType = transParam.firingprocid.get(modeObj);
                if (procType == null) continue;
                if (SKIP_BERNSTEIN_CONVERSION.contains(procType)) continue;

                MatrixCell origProc = transParam.firingproc.get(modeObj);
                if (origProc == null || origProc.isEmpty()) continue;

                String distName = ProcessType.toText(procType);
                InputOutput.line_warning(
                        "snNonmarkovToPh",
                        "Firing distribution " + distName + " at Transition node " + ind + " mode " + m
                                + " is non-Markovian and will be converted to PH (" + nPhases + " phases).");

                final double targetMean;
                DoubleUnaryOperator pdfFunc;
                if (procType == ProcessType.GAMMA) {
                    final double shape = origProc.get(0).toDouble();
                    final double scale = origProc.get(1).toDouble();
                    targetMean = new Gamma(shape, scale).getMean();
                    pdfFunc = new DoubleUnaryOperator() {
                        @Override
                        public double applyAsDouble(double x) {
                            if (x <= 0) return 0.0;
                            return new Gamma(shape, scale).evalPDF(x);
                        }
                    };
                } else if (procType == ProcessType.WEIBULL) {
                    // Stored order matches Weibull.getProcess: {r (shape), alpha (scale)}
                    final double shapeParam = origProc.get(0).toDouble();
                    final double scaleParam = origProc.get(1).toDouble();
                    targetMean = new Weibull(shapeParam, scaleParam).getMean();
                    pdfFunc = new DoubleUnaryOperator() {
                        @Override
                        public double applyAsDouble(double x) {
                            if (x <= 0) return 0.0;
                            return new Weibull(shapeParam, scaleParam).evalPDF(x);
                        }
                    };
                } else if (procType == ProcessType.LOGNORMAL) {
                    final double mu = origProc.get(0).toDouble();
                    final double sigma = origProc.get(1).toDouble();
                    targetMean = new Lognormal(mu, sigma).getMean();
                    pdfFunc = new DoubleUnaryOperator() {
                        @Override
                        public double applyAsDouble(double x) {
                            if (x <= 0) return 0.0;
                            return new Lognormal(mu, sigma).evalPDF(x);
                        }
                    };
                } else if (procType == ProcessType.PARETO) {
                    final double shapeParam = origProc.get(0).toDouble();
                    final double scaleParam = origProc.get(1).toDouble();
                    // Undefined for alpha <= 1; MATLAB yields NaN and falls back below.
                    targetMean = (shapeParam > 1.0)
                            ? shapeParam * scaleParam / (shapeParam - 1.0)
                            : Double.NaN;
                    pdfFunc = new DoubleUnaryOperator() {
                        @Override
                        public double applyAsDouble(double x) {
                            if (x < scaleParam) return 0.0;
                            return shapeParam * Math.pow(scaleParam, shapeParam) / Math.pow(x, shapeParam + 1);
                        }
                    };
                } else if (procType == ProcessType.UNIFORM) {
                    final double minVal = origProc.get(0).toDouble();
                    final double maxVal = origProc.get(1).toDouble();
                    targetMean = (minVal + maxVal) / 2.0;
                    pdfFunc = new DoubleUnaryOperator() {
                        @Override
                        public double applyAsDouble(double x) {
                            if (x >= minVal && x <= maxVal) return 1.0 / (maxVal - minVal);
                            return 0.0;
                        }
                    };
                } else if (procType == ProcessType.DET) {
                    if (preserveDet) continue;
                    updateNodeparamForMAP(transParam, modeObj, m,
                            Map_erlang.map_erlang(origProc.get(0).toDouble(), nPhases));
                    continue;
                } else {
                    continue;
                }

                MatrixCell map;
                if (!isFinitePositive(targetMean)) {
                    map = Map_erlang.map_erlang(1.0, nPhases);
                } else {
                    try {
                        Pair<Matrix, Matrix> dd = Aph_bernstein.aph_bernstein(pdfFunc, nPhases);
                        map = Map_scale.map_scale(pairToMatrixCell(dd.getLeft(), dd.getRight()), targetMean);
                    } catch (RuntimeException e) {
                        map = Map_erlang.map_erlang(targetMean, nPhases);
                    }
                }
                updateNodeparamForMAP(transParam, modeObj, m, map);
            }
        }
    }

    private static boolean isFinitePositive(double v) {
        return !Double.isNaN(v) && !Double.isInfinite(v) && v > 0;
    }

    /**
     * Update a Transition mode's firing process to a phase-type representation.
     * Mirrors updateNodeparamForMAP in sn_nonmarkov_toph.m.
     */
    private static void updateNodeparamForMAP(TransitionNodeParam transParam, Mode modeObj, int m, MatrixCell map) {
        int actualPhases = map.get(0).getNumRows();
        transParam.firingproc.put(modeObj, map);
        transParam.firingpie.put(modeObj, Map_pie.map_pie(map));
        transParam.firingprocid.put(modeObj, ProcessType.MAP);
        if (transParam.firingphases != null && transParam.firingphases.length() > m) {
            transParam.firingphases.set(m, actualPhases);
        }
    }

    private static MatrixCell pairToMatrixCell(Matrix d0, Matrix d1) {
        MatrixCell cell = new MatrixCell();
        cell.set(0, d0);
        cell.set(1, d1);
        return cell;
    }

    private static void updateSnForMAP(NetworkStruct sn, Station station, JobClass jobClass,
                                       int ist, int r, MatrixCell map, int nPhases, ProcessType origProcType,
                                       boolean keepDetTag) {
        Matrix d0 = map.get(0);
        Matrix d1 = map.get(1);

        int oldPhases = (int) sn.phasessz.get(ist, r);

        Map<JobClass, MatrixCell> stProc = sn.proc.get(station);
        if (stProc != null) stProc.put(jobClass, map);
        Map<JobClass, ProcessType> stProcid = sn.procid.get(station);
        if (stProcid != null) {
            if (!keepDetTag) {
                // see _kb/03-api-layer.md for rationale
                stProcid.put(jobClass, ProcessType.APH);
            } else if (origProcType != ProcessType.DET) {
                // MAM path (matrix-analytic, uses the proc D0/D1 directly): tag MAP.
                // DET is left as DET for the exact MAP/D/c Crommelin dispatch.
                stProcid.put(jobClass, ProcessType.MAP);
            }
        }

        sn.phases.set(ist, r, (double) nPhases);
        sn.phasessz.set(ist, r, (double) Math.max(nPhases, 1));

        double cumSum = 0.0;
        sn.phaseshift.set(ist, 0, 0.0);
        for (int c = 0; c < sn.nclasses; c++) {
            cumSum += sn.phasessz.get(ist, c);
            if (c + 1 < sn.phaseshift.getNumCols()) {
                sn.phaseshift.set(ist, c + 1, cumSum);
            }
        }

        Matrix muMatrix = new Matrix(nPhases, 1);
        for (int i = 0; i < nPhases; i++) {
            muMatrix.set(i, 0, -d0.get(i, i));
        }
        if (sn.mu.get(station) == null) {
            sn.mu.put(station, new HashMap<JobClass, Matrix>());
        }
        sn.mu.get(station).put(jobClass, muMatrix);

        Matrix phiMatrix = new Matrix(nPhases, 1);
        for (int i = 0; i < nPhases; i++) {
            double d1RowSum = 0.0;
            for (int j = 0; j < d1.getNumCols(); j++) {
                d1RowSum += d1.get(i, j);
            }
            double d0Diag = -d0.get(i, i);
            phiMatrix.set(i, 0, d0Diag != 0.0 ? d1RowSum / d0Diag : 0.0);
        }
        if (sn.phi.get(station) == null) {
            sn.phi.put(station, new HashMap<JobClass, Matrix>());
        }
        sn.phi.get(station).put(jobClass, phiMatrix);

        Matrix pieMatrix = Map_pie.map_pie(map);
        if (sn.pie.get(station) == null) {
            sn.pie.put(station, new HashMap<JobClass, Matrix>());
        }
        sn.pie.get(station).put(jobClass, pieMatrix);

        int ind = (int) sn.stationToNode.get(ist);
        // see _kb/03-api-layer.md for rationale
        boolean addNvar = keepDetTag;
        if (addNvar) {
            sn.nvars.set(ind, r, sn.nvars.get(ind, r) + 1);
        }

        expandStateForMAP(sn, ind, r, oldPhases, nPhases, addNvar);

        if (nPhases > 1) {
            addPhaseSyncIfNeeded(sn, ind, r);
        }
    }

    private static void expandStateForMAP(NetworkStruct sn, int ind, int r, int oldPhases, int newPhases, boolean addNvar) {
        int isf = (int) sn.nodeToStateful.get(ind);
        if (isf < 0 || sn.state.isEmpty()) return;

        StatefulNode statefulNode = null;
        for (Map.Entry<StatefulNode, Matrix> entry : sn.state.entrySet()) {
            if (entry.getKey() != null && entry.getKey().getStatefulIndex() == isf) {
                statefulNode = entry.getKey();
                break;
            }
        }
        if (statefulNode == null) return;

        Matrix stateMatrix = sn.state.get(statefulNode);
        if (stateMatrix == null || stateMatrix.isEmpty()) return;

        int ist = (int) sn.nodeToStation.get(ind);
        int nRows = stateMatrix.getNumRows();

        int V = (int) sn.nvars.getRow(ind).elementSum();
        // see _kb/03-api-layer.md for rationale
        int V_old = addNvar ? V - 1 : V;

        int sumK_new = (int) sn.phasessz.getRow(ist).elementSum();
        int sumK_old = 0;
        for (int c = 0; c < sn.nclasses; c++) {
            sumK_old += (c == r) ? oldPhases : (int) sn.phasessz.get(ist, c);
        }

        int currentCols = stateMatrix.getNumCols();
        int bufSize = currentCols - sumK_old - V_old;
        if (bufSize < 0) bufSize = 0;

        Matrix spaceBuf = (bufSize > 0)
                ? Matrix.extract(stateMatrix, 0, nRows, 0, bufSize)
                : new Matrix(nRows, 0);

        Matrix spaceSrv = (sumK_old > 0)
                ? Matrix.extract(stateMatrix, 0, nRows, bufSize, bufSize + sumK_old)
                : new Matrix(nRows, 0);

        Matrix spaceVar = (V_old > 0)
                ? Matrix.extract(stateMatrix, 0, nRows, bufSize + sumK_old, currentCols)
                : new Matrix(nRows, 0);

        int phasesToAdd = newPhases - oldPhases;
        Matrix spaceSrvNew;
        if (phasesToAdd > 0 && spaceSrv.getNumCols() > 0) {
            int insertPos = 0;
            for (int c = 0; c < r; c++) {
                insertPos += (c == r) ? oldPhases : (int) sn.phasessz.get(ist, c);
            }
            insertPos += oldPhases;

            Matrix expanded = new Matrix(nRows, spaceSrv.getNumCols() + phasesToAdd);
            expanded.zero();
            for (int row = 0; row < nRows; row++) {
                for (int col = 0; col < insertPos; col++) {
                    if (col < spaceSrv.getNumCols()) {
                        expanded.set(row, col, spaceSrv.get(row, col));
                    }
                }
                for (int col = insertPos; col < spaceSrv.getNumCols(); col++) {
                    expanded.set(row, col + phasesToAdd, spaceSrv.get(row, col));
                }
            }
            spaceSrvNew = expanded;
        } else {
            spaceSrvNew = spaceSrv;
        }

        // MAP mode appends a phase-tracking nvar column (initialized to phase 1);
        // server-phase mode leaves the local-variable block unchanged.
        Matrix spaceVarNew;
        if (addNvar) {
            spaceVarNew = new Matrix(nRows, spaceVar.getNumCols() + 1);
            for (int row = 0; row < nRows; row++) {
                for (int col = 0; col < spaceVar.getNumCols(); col++) {
                    spaceVarNew.set(row, col, spaceVar.get(row, col));
                }
                spaceVarNew.set(row, spaceVar.getNumCols(), 1.0);
            }
        } else {
            spaceVarNew = spaceVar;
        }

        Matrix newState = Matrix.concatColumns(spaceBuf, spaceSrvNew, null);
        Matrix finalState = Matrix.concatColumns(newState, spaceVarNew, null);
        sn.state.put(statefulNode, finalState);
    }

    private static void addPhaseSyncIfNeeded(NetworkStruct sn, int ind, int r) {
        int local = sn.nnodes;

        for (Sync syncEntry : sn.sync.values()) {
            if (syncEntry != null && syncEntry.active != null && !syncEntry.active.isEmpty()) {
                Event activeEvent = syncEntry.active.get(0);
                if (activeEvent != null
                        && activeEvent.getEvent() == EventType.PHASE
                        && activeEvent.getNode() == ind
                        && activeEvent.getJobClass() == r) {
                    return;
                }
            }
        }

        Sync newSync = new Sync();
        Event activeEvent = new Event(EventType.PHASE, ind, r);
        Event passiveEvent = new Event(EventType.LOCAL, local, r);
        passiveEvent.setProb(1.0);
        newSync.active.put(0, activeEvent);
        newSync.passive.put(0, passiveEvent);

        int nextKey = sn.sync.size();
        sn.sync.put(nextKey, newSync);
    }
}
