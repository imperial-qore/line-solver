/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.solvers;

import static jline.GlobalConstants.Inf;

import jline.GlobalConstants;
import jline.VerboseLevel;
import jline.io.Ret.DistributionResult;
import jline.io.Ret.ProbabilityResult;
import jline.io.Ret.SampleResult;
import jline.lang.Network;
import jline.lang.NetworkStruct;
import jline.lang.constant.NodeType;
import jline.solvers.ctmc.SolverCTMC;
import jline.solvers.des.SolverDES;
import jline.solvers.fluid.SolverFluid;
import jline.solvers.jmt.SolverJMT;
import jline.solvers.mam.SolverMAM;
import jline.solvers.mva.SolverMVA;
import jline.solvers.nc.SolverNC;
import jline.solvers.ssa.SolverSSA;
import jline.util.matrix.Matrix;

import javax.xml.parsers.ParserConfigurationException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import static jline.api.sn.SnGetArvRFromTputKt.snGetArvRFromTput;
import static jline.api.sn.SnGetDemandsChainKt.snGetDemandsChain;
import static jline.api.sn.SnGetNodeArvRFromTputKt.snGetNodeArvRFromTput;
import static jline.api.sn.SnGetNodeTputFromTputKt.snGetNodeTputFromTput;
import static jline.api.sn.SnGetResidTFromRespTKt.snGetResidTFromRespT;
import static jline.io.InputOutputKt.*;
import static jline.util.Utils.isInf;

/**
 * Abstract base class for solvers applicable to queueing network models.
 * <p>
 * This class provides the core functionality for analyzing queueing networks using various
 * solution algorithms. It manages performance metrics computation including queue lengths,
 * utilizations, response times, throughputs, and arrival rates at both steady-state and
 * transient conditions.
 * <p>
 * The solver operates on {@link Network} models and produces results through various
 * table formats for different levels of aggregation (station-level, node-level, chain-level).
 *
 * @see Network
 * @see SolverOptions
 * @see SolverResult
 */
public abstract class NetworkSolver extends Solver {

    // ========== Class Fields ==========
    
    /**
     * The queueing network model to be solved
     */
    public Network model;

    /**
     * Internal data structure describing the network model
     */
    public NetworkStruct sn;

    /**
     * Handles for steady-state average performance metrics
     */
    public SolverAvgHandles avgHandles;

    /**
     * Handles for transient performance metrics
     */
    public SolverTranHandles tranHandles;

    /**
     * Constructs a NetworkSolver with the specified model, name, and options.
     *
     * @param model   the queueing network model to solve
     * @param name    the name identifier for this solver instance
     * @param options configuration options for the solver
     * @throws RuntimeException if the model is empty (has no nodes)
     */
    protected NetworkSolver(Network model, String name, SolverOptions options) {
        super(name, options);
        this.model = model;
        // Allow null model for LayeredNetwork-based solvers (e.g., SolverDES with LQN)
        if (model != null) {
            if (model.getNumberOfNodes() == 0) {
                throw new RuntimeException("The model supplied in input is empty.");
            }
            this.avgHandles = model.getAvgHandles();
            this.tranHandles = model.getTranHandles();
            this.sn = model.getStruct(true); // Force model to refresh
        }
    }

    /**
     * Constructs a NetworkSolver with the specified model and name using default options.
     *
     * @param model the queueing network model to solve
     * @param name  the name identifier for this solver instance
     */
    protected NetworkSolver(Network model, String name) {
        this(model, name, defaultOptions());
    }

    /**
     * Returns a list containing instances of all available network solvers for the given model.
     *
     * @param model the queueing network model
     * @return list of all available solver implementations
     */
    public static List<NetworkSolver> getAllSolvers(Network model) {
        SolverOptions options = new SolverOptions();
        List<NetworkSolver> solvers = new ArrayList<>();
        solvers.add(new SolverCTMC(model, options));
        solvers.add(new SolverDES(model, options));
        solvers.add(new SolverFluid(model, options));
        solvers.add(new SolverJMT(model, options));
        solvers.add(new SolverMAM(model, options));
        solvers.add(new SolverMVA(model, options));
        solvers.add(new SolverNC(model, options));
        solvers.add(new SolverSSA(model, options));
        return solvers;
    }

    // Basic solver results
    public SolverResult avg() {
        return getAvg();
    }

    public SolverResult avg(SolverAvgHandles avgHandles) {
        return getAvg(avgHandles);
    }

    public SolverResult avg(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A) {
        return getAvg(Q, U, R, W, T, A);
    }

    // Arrival rates
    public Matrix avgArvR() {
        return getAvgArvR();
    }

    public Matrix avgArvRChain() {
        return getAvgArvRChain();
    }

    public AvgHandle avgArvRHandles() {
        return getAvgArvRHandles();
    }

    // Chain-level results
    public SolverResult avgChain() {
        return getAvgChain();
    }

    public NetworkAvgChainTable avgChainTable() {
        return getAvgChainTable();
    }

    public NetworkAvgChainTable avgChainTable(SolverAvgHandles avgHandles) {
        return getAvgChainTable(avgHandles);
    }

    public NetworkAvgChainTable avgChainTable(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A) {
        return getAvgChainTable(Q, U, R, W, T, A);
    }

    public NetworkAvgChainTable avgChainTable(boolean keepDisabled) {
        return getAvgChainTable(keepDisabled);
    }

    public NetworkAvgChainTable avgChainTable(SolverAvgHandles avgHandles, boolean keepDisabled) {
        return getAvgChainTable(avgHandles, keepDisabled);
    }

    public NetworkAvgChainTable avgChainTable(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A, boolean keepDisabled) {
        return getAvgChainTable(Q, U, R, W, T, A, keepDisabled);
    }

    // avgChainTable -> avgChainT aliases
    public NetworkAvgChainTable avgChainT() {
        return avgChainTable();
    }

    public NetworkAvgChainTable avgChainT(SolverAvgHandles avgHandles) {
        return avgChainTable(avgHandles);
    }

    public NetworkAvgChainTable avgChainT(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A) {
        return avgChainTable(Q, U, R, W, T, A);
    }

    public NetworkAvgChainTable avgChainT(boolean keepDisabled) {
        return avgChainTable(keepDisabled);
    }

    public NetworkAvgChainTable avgChainT(SolverAvgHandles avgHandles, boolean keepDisabled) {
        return avgChainTable(avgHandles, keepDisabled);
    }

    public NetworkAvgChainTable avgChainT(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A, boolean keepDisabled) {
        return avgChainTable(Q, U, R, W, T, A, keepDisabled);
    }

    // Handles
    public SolverAvgHandles avgHandles() {
        return getAvgHandles();
    }

    // Node-level results
    public SolverResult avgNode() {
        return getAvgNode();
    }

    public Matrix avgNodeArvRChain() {
        return getAvgNodeArvRChain();
    }

    public SolverResult avgNodeChain() {
        return getAvgNodeChain();
    }

    public NetworkAvgNodeChainTable avgNodeChainTable() {
        return getAvgNodeChainTable();
    }

    public NetworkAvgNodeChainTable avgNodeChainTable(SolverAvgHandles avgHandles) {
        return getAvgNodeChainTable(avgHandles);
    }

    public NetworkAvgNodeChainTable avgNodeChainTable(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A) {
        return getAvgNodeChainTable(Q, U, R, W, T, A);
    }

    public NetworkAvgNodeChainTable avgNodeChainTable(boolean keepDisabled) {
        return getAvgNodeChainTable(keepDisabled);
    }

    public NetworkAvgNodeChainTable avgNodeChainTable(SolverAvgHandles avgHandles, boolean keepDisabled) {
        return getAvgNodeChainTable(avgHandles, keepDisabled);
    }

    public NetworkAvgNodeChainTable avgNodeChainTable(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A, boolean keepDisabled) {
        return getAvgNodeChainTable(Q, U, R, W, T, A, keepDisabled);
    }

    // avgNodeChainTable -> avgNodeChainT aliases
    public NetworkAvgNodeChainTable avgNodeChainT() {
        return avgNodeChainTable();
    }

    public NetworkAvgNodeChainTable avgNodeChainT(SolverAvgHandles avgHandles) {
        return avgNodeChainTable(avgHandles);
    }

    public NetworkAvgNodeChainTable avgNodeChainT(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A) {
        return avgNodeChainTable(Q, U, R, W, T, A);
    }

    public NetworkAvgNodeChainTable avgNodeChainT(boolean keepDisabled) {
        return avgNodeChainTable(keepDisabled);
    }

    public NetworkAvgNodeChainTable avgNodeChainT(SolverAvgHandles avgHandles, boolean keepDisabled) {
        return avgNodeChainTable(avgHandles, keepDisabled);
    }

    public NetworkAvgNodeChainTable avgNodeChainT(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A, boolean keepDisabled) {
        return avgNodeChainTable(Q, U, R, W, T, A, keepDisabled);
    }

    public Matrix avgNodeQLenChain() {
        return getAvgNodeQLenChain();
    }

    public Matrix avgNodeResidTChain() {
        return getAvgNodeResidTChain();
    }

    public Matrix avgNodeRespTChain() {
        return getAvgNodeRespTChain();
    }

    public NetworkAvgNodeTable avgNodeTable() {
        return getAvgNodeTable();
    }

    public NetworkAvgNodeTable avgNodeTable(SolverAvgHandles avgHandles) {
        return getAvgNodeTable(avgHandles);
    }

    public NetworkAvgNodeTable avgNodeTable(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A) {
        return getAvgNodeTable(Q, U, R, W, T, A);
    }

    public NetworkAvgNodeTable avgNodeTable(boolean keepDisabled) {
        return getAvgNodeTable(keepDisabled);
    }

    public NetworkAvgNodeTable avgNodeTable(SolverAvgHandles avgHandles, boolean keepDisabled) {
        return getAvgNodeTable(avgHandles, keepDisabled);
    }

    public NetworkAvgNodeTable avgNodeTable(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A, boolean keepDisabled) {
        return getAvgNodeTable(Q, U, R, W, T, A, keepDisabled);
    }

    // avgNodeTable -> avgNodeT aliases
    public NetworkAvgNodeTable avgNodeT() {
        return avgNodeTable();
    }

    public NetworkAvgNodeTable avgNodeT(SolverAvgHandles avgHandles) {
        return avgNodeTable(avgHandles);
    }

    public NetworkAvgNodeTable avgNodeT(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A) {
        return avgNodeTable(Q, U, R, W, T, A);
    }

    public NetworkAvgNodeTable avgNodeT(boolean keepDisabled) {
        return avgNodeTable(keepDisabled);
    }

    public NetworkAvgNodeTable avgNodeT(SolverAvgHandles avgHandles, boolean keepDisabled) {
        return avgNodeTable(avgHandles, keepDisabled);
    }

    public NetworkAvgNodeTable avgNodeT(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A, boolean keepDisabled) {
        return avgNodeTable(Q, U, R, W, T, A, keepDisabled);
    }

    public Matrix avgNodeTputChain() {
        return getAvgNodeTputChain();
    }

    public Matrix avgNodeUtilChain() {
        return getAvgNodeUtilChain();
    }

    // Queue lengths
    public Matrix avgQLen() {
        return getAvgQLen();
    }

    public Matrix avgQLenChain() {
        return getAvgQLenChain();
    }

    public AvgHandle avgQLenHandles() {
        return getAvgQLenHandles();
    }

    // Residence times
    public Matrix avgResidT() {
        return getAvgResidT();
    }

    public Matrix avgResidTChain() {
        return getAvgResidTChain();
    }

    public AvgHandle avgResidTHandles() {
        return getAvgResidTHandles();
    }

    // Response times
    public Matrix avgRespT() {
        return getAvgRespT();
    }

    public Matrix avgRespTChain() {
        return getAvgRespTChain();
    }

    public AvgHandle avgRespTHandles() {
        return getAvgRespTHandles();
    }

    // System-level results
    public void avgSys() {
        getAvgSys();
    }

    public Matrix avgSysRespT() {
        return getAvgSysRespT();
    }

    public NetworkAvgSysTable avgSysTable() {
        return getAvgSysTable();
    }

    public NetworkAvgSysTable avgSysTable(SolverAvgHandles avgHandles) {
        return getAvgSysTable(avgHandles);
    }

    public NetworkAvgSysTable avgSysTable(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A) {
        return getAvgSysTable(Q, U, R, W, T, A);
    }

    // avgSysTable(boolean keepDisabled) - method signature not available
    // avgSysTable with handles and boolean - method signature not available
    // avgSysTable with 6 handles and boolean - method signature not available

    // avgSysTable -> avgSysT aliases
    public NetworkAvgSysTable avgSysT() {
        return avgSysTable();
    }

    public NetworkAvgSysTable avgSysT(SolverAvgHandles avgHandles) {
        return avgSysTable(avgHandles);
    }

    public NetworkAvgSysTable avgSysT(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A) {
        return avgSysTable(Q, U, R, W, T, A);
    }

    public Matrix avgSysTput() {
        return getAvgSysTput();
    }

    // Tables
    public NetworkAvgTable avgTable() {
        return getAvgTable();
    }

    public NetworkAvgTable avgTable(SolverAvgHandles avgHandles) {
        return getAvgTable(avgHandles);
    }

    public NetworkAvgTable avgTable(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A) {
        return getAvgTable(Q, U, R, W, T, A);
    }

    public NetworkAvgTable avgTable(boolean keepDisabled) {
        return getAvgTable(keepDisabled);
    }

    public NetworkAvgTable avgTable(SolverAvgHandles avgHandles, boolean keepDisabled) {
        return getAvgTable(avgHandles, keepDisabled);
    }

    public NetworkAvgTable avgTable(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A, boolean keepDisabled) {
        return getAvgTable(Q, U, R, W, T, A, keepDisabled);
    }

    // Table -> T aliases
    public NetworkAvgTable avgT() {
        return avgTable();
    }

    public NetworkAvgTable avgT(SolverAvgHandles avgHandles) {
        return avgTable(avgHandles);
    }

    public NetworkAvgTable avgT(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A) {
        return avgTable(Q, U, R, W, T, A);
    }

    public NetworkAvgTable avgT(boolean keepDisabled) {
        return avgTable(keepDisabled);
    }

    public NetworkAvgTable avgT(SolverAvgHandles avgHandles, boolean keepDisabled) {
        return avgTable(avgHandles, keepDisabled);
    }

    public NetworkAvgTable avgT(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A, boolean keepDisabled) {
        return avgTable(Q, U, R, W, T, A, keepDisabled);
    }

    // aT aliases for getAvgTable
    public NetworkAvgTable aT() {
        return getAvgTable();
    }

    public NetworkAvgTable aT(SolverAvgHandles avgHandles) {
        return getAvgTable(avgHandles);
    }

    public NetworkAvgTable aT(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A) {
        return getAvgTable(Q, U, R, W, T, A);
    }

    public NetworkAvgTable aT(boolean keepDisabled) {
        return getAvgTable(keepDisabled);
    }

    public NetworkAvgTable aT(SolverAvgHandles avgHandles, boolean keepDisabled) {
        return getAvgTable(avgHandles, keepDisabled);
    }

    public NetworkAvgTable aT(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A, boolean keepDisabled) {
        return getAvgTable(Q, U, R, W, T, A, keepDisabled);
    }

    // Throughputs
    public Matrix avgTput() {
        return getAvgTput();
    }

    public Matrix avgTputChain() {
        return getAvgTputChain();
    }

    public AvgHandle avgTputHandles() {
        return getAvgTputHandles();
    }

    // Utilizations
    public Matrix avgUtil() {
        return getAvgUtil();
    }

    public Matrix avgUtilChain() {
        return getAvgUtilChain();
    }

    public AvgHandle avgUtilHandles() {
        return getAvgUtilHandles();
    }

    // Wait times
    public Matrix avgWaitT() {
        return getAvgWaitT();
    }

    // Distribution functions
    public DistributionResult cdfPassT() {
        return getCdfPassT();
    }

    public DistributionResult cdfPassT(AvgHandle R) {
        return getCdfPassT(R);
    }

    public DistributionResult cdfRespT() {
        return getCdfRespT();
    }

    public DistributionResult cdfRespT(AvgHandle R) {
        return getCdfRespT(R);
    }

    /**
     * Filters and processes performance metrics based on enabled handles and masks.
     * This internal method handles disabling metrics, applying zero masks, and
     * cleaning up numerical artifacts.
     *
     * @param avgHandle metric handle configuration
     * @param metric    raw metric matrix to filter
     * @param zeroMask  optional mask to zero out specific entries
     * @return filtered and processed metric matrix
     */
    private Matrix filterMetric(AvgHandle avgHandle, Matrix metric, Matrix zeroMask) {
        int M = sn.nstations;
        int K = sn.nclasses;
        Matrix outData = new Matrix(M, K);
        for (int k = 0; k < K; k++) {
            for (int i = 0; i < M; i++) {
                if (!avgHandle.get(this.model.getStations().get(i)).get(this.model.getClassByIndex(k)).isDisabled
                        && metric != null && !metric.isEmpty()) {
                    double value = metric.get(i, k);
                    if (value < GlobalConstants.FineTol || Double.isNaN(value)) {
                        value = 0; // Round to zero or disable the metric
                    }
                    outData.set(i, k, value);
                } else {
                    outData.set(i, k, Double.NaN); // Indicates that a metric is disabled
                }
            }
        }

        // NaN values indicate that a metric is disabled
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < K; j++) {
                if (Double.isNaN(outData.get(i, j))) {
                    outData.set(i, j, 0);
                }
            }
        }

        // Set to zero entries in the mask
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < K; j++) {
                if (zeroMask != null && !zeroMask.isEmpty() && (zeroMask.get(i, j) == 1)) {
                    outData.set(i, j, 0);
                }
            }
        }

        // Round to zero numerical perturbations
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < K; j++) {
                if (outData.get(i, j) < GlobalConstants.FineTol) {
                    outData.set(i, j, 0);
                }
            }
        }

        // set to zero metrics for classes that are unreachable
        // skip when no chains are defined (routing not specified)
        // skip for SPN models where Places don't have traditional visits
        boolean hasSPN = sn.nodetype.contains(NodeType.Place) || sn.nodetype.contains(NodeType.Transition);

        if (sn.nchains > 0 && !hasSPN) {
            for (int k = 0; k < K; k++) {
                Matrix chainCol = sn.chains.getColumn(k);
                if (chainCol.elementMax() > 0) { // Only check if class k belongs to a chain
                    int c = (int) chainCol.find().value();
                    for (int i = 0; i < M; i++) {
                        if (sn.visits.get(c).get(i, k) == 0) {
                            outData.set(i, k, 0);
                        }
                    }
                }
            }
        }

        return outData;
    }

    /**
     * Computes and returns average station metrics at steady-state.
     * This is the main method for obtaining steady-state performance metrics.
     *
     * @return solver result containing station-level average metrics
     * @throws RuntimeException if unable to compute results
     */
    public SolverResult getAvg() {

        if (this.avgHandles == null || this.avgHandles.Q == null || this.avgHandles.U == null || this.avgHandles.R == null ||
                this.avgHandles.T == null || this.avgHandles.A == null) {
            reset();
        }
        this.avgHandles = model.getAvgHandles();

        if (!isInf(options.timespan[1]) && !Double.isNaN(options.timespan[1])) {
            throw new RuntimeException(
                    "The getAvg method does not support the timespan option, use the getTranAvg method instead.");
        } else {
            this.options.timespan[0] = Inf;
            this.options.timespan[1] = Inf;
        }

        if (!this.hasAvgResults() || !this.options.cache) {
            try {
                runAnalyzer();
            } catch (IllegalAccessException e) {
                String errorMsg = "IllegalAccessException upon running runAnalyzer(): " + e.getMessage();
                line_error(mfilename(new Object() {
                }), errorMsg);
                throw new RuntimeException("Failed to execute analyzer due to access restrictions: " + e.getMessage(), e);
            } catch (ParserConfigurationException e) {
                String errorMsg = "ParserConfigurationException upon running runAnalyzer(): " + e.getMessage();
                line_error(mfilename(new Object() {
                }), errorMsg);
                throw new RuntimeException("Failed to configure XML parser: " + e.getMessage(), e);
            } catch (IOException e) {
                String errorMsg = "IOException upon running runAnalyzer(): " + e.getMessage();
                line_error(mfilename(new Object() {
                }), errorMsg);
                throw new RuntimeException("I/O error during analyzer execution: " + e.getMessage(), e);
            }
            if (!this.hasAvgResults()) {
                throw new RuntimeException("Line is unable to return results for this model. " +
                        "This may indicate that the solver " + this.getName() + " is not applicable to this model type, " +
                        "or the model parameters are invalid.");
            }
        } // else return cached value

        int M = sn.nstations;
        int K = sn.nclasses;

        int Vrows = sn.visits.get(0).getNumRows();
        int Vcols = sn.visits.get(0).getNumCols();
        int Vcells = sn.visits.size();
        Matrix V = new Matrix(Vrows, Vcols);
        for (int i = 0; i < Vrows; i++) {
            for (int j = 0; j < Vcols; j++) {
                double tmpSum = 0;
                for (int k = 0; k < Vcells; k++) {
                    tmpSum += sn.visits.get(k).get(i, j);
                }
                V.set(i, j, tmpSum);
            }
        }

        Matrix QNclass = new Matrix(0, 0);
        Matrix UNclass = new Matrix(0, 0);
        Matrix RNclass = new Matrix(0, 0);
        Matrix TNclass = new Matrix(0, 0);
        Matrix ANclass = new Matrix(0, 0);
        Matrix WNclass = new Matrix(0, 0);

        if (!this.avgHandles.R.isEmpty() && this.result.RN != null && !this.result.RN.isEmpty()) {
            RNclass = filterMetric(this.avgHandles.R, this.result.RN, null);
        }

        Matrix zeroMask = new Matrix(M, K);
        for (int i = 0; i < M; i++)
            for (int r = 0; r < K; r++)
                if (RNclass.get(i, r) < 10 * GlobalConstants.FineTol)
                    zeroMask.set(i, r, 1);

        // For queue lengths, exempt Place nodes from zeroMask since they don't have response times
        // (Place nodes track token counts, not job response times)
        Matrix zeroMaskQLen = zeroMask.copy();
        for (int ist = 0; ist < M; ist++) {
            if (sn.nodetype.get((int) sn.stationToNode.get(ist)) == NodeType.Place) {
                for (int r = 0; r < K; r++) {
                    zeroMaskQLen.set(ist, r, 0);  // Don't mask Place nodes
                }
            }
        }

        if (!this.avgHandles.Q.isEmpty() && this.result.QN != null && !this.result.QN.isEmpty()) {
            QNclass = filterMetric(this.avgHandles.Q, this.result.QN, zeroMaskQLen);
        }

        if (!this.avgHandles.U.isEmpty() && this.result.UN != null && !this.result.UN.isEmpty()) {
            UNclass = filterMetric(this.avgHandles.U, this.result.UN, zeroMask);
        }

        if (!this.avgHandles.T.isEmpty() && this.result.TN != null && !this.result.TN.isEmpty()) {
            TNclass = filterMetric(this.avgHandles.T, this.result.TN, null);
        }

        if (!this.avgHandles.A.isEmpty()) {
            // If result.AN is available, use it; otherwise compute from throughputs
            if (this.result.AN != null && !this.result.AN.isEmpty()) {
                Matrix zeroMaskSource = Matrix.createLike(zeroMask);
                for (int ist = 0; ist < M; ist++) {
                    if (sn.nodetype.get((int) sn.stationToNode.get(ist)) == NodeType.Source) {
                        for (int r = 0; r < K; r++) {
                            zeroMaskSource.set(ist, r, 1);
                        }
                    }
                }
                ANclass = filterMetric(this.avgHandles.A, this.result.AN, zeroMaskSource);
            } else if (!TNclass.isEmpty()) {
                // Compute arrival rates from throughputs when result.AN is not available
                Matrix computedAN = snGetArvRFromTput(sn, TNclass, this.avgHandles.T);
                Matrix zeroMaskSource = Matrix.createLike(zeroMask);
                for (int ist = 0; ist < M; ist++) {
                    if (sn.nodetype.get((int) sn.stationToNode.get(ist)) == NodeType.Source) {
                        for (int r = 0; r < K; r++) {
                            zeroMaskSource.set(ist, r, 1);
                        }
                    }
                }
                ANclass = filterMetric(this.avgHandles.A, computedAN, zeroMaskSource);
            }
        }

        if (!this.avgHandles.W.isEmpty()) {
            WNclass = snGetResidTFromRespT(sn, RNclass, this.avgHandles.W);
        }

        if (!UNclass.isEmpty()) {
            boolean unstableQueueFlag = false;
            for (int i = 0; i < M; i++) {
                double Uiopen = 0.0;
                for (int j = 0; j < K; j++) {
                    // sum util of open classes
                    if (isInf(sn.njobs.get(j))) {
                        Uiopen += UNclass.get(i, j);
                    }
                }
                if (Uiopen > 0.99 * sn.nservers.get(i, 0)) {
                    unstableQueueFlag = true;
                    break;
                }
            }

            boolean infJobsFlag = false;
            for (int i = 0; i < sn.njobs.length(); i++) {
                if (isInf(sn.njobs.get(0, i))) {
                    infJobsFlag = true;
                    break;
                }
            }

            if (unstableQueueFlag && infJobsFlag) {
                line_warning(mfilename(new Object() {
                        }),
                        "The model has unstable queues, performance metrics may grow unbounded.");
            }
        }
        double runtime = this.result.runtime;
        // Note: this.result.reset() is intentionally not called here as it would erase
        // transient measures (QNt, UNt, TNt, t, pi_t) that may be needed for analysis
        this.result.QN = QNclass;
        this.result.UN = UNclass;
        this.result.RN = RNclass;
        this.result.AN = ANclass;
        this.result.TN = TNclass;
        this.result.WN = WNclass;
        this.result.runtime = runtime;
        return this.result;
    }

    /**
     * Computes and returns average station metrics at steady-state using specified handles.
     *
     * @param avgHandles custom handles for performance metrics
     * @return solver result containing station-level average metrics
     * @throws RuntimeException if unable to compute results
     */
    public SolverResult getAvg(SolverAvgHandles avgHandles) {
        SolverAvgHandles originalHandles = this.avgHandles;
        this.avgHandles = avgHandles;
        try {
            return getAvg();
        } finally {
            this.avgHandles = originalHandles;
        }
    }

    /**
     * Computes and returns average station metrics at steady-state using individual handles.
     *
     * @param Q queue length handle
     * @param U utilization handle
     * @param R response time handle
     * @param W residence time handle
     * @param T throughput handle
     * @param A arrival rate handle
     * @return solver result containing station-level average metrics
     * @throws RuntimeException if unable to compute results
     */
    public SolverResult getAvg(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A) {
        SolverAvgHandles customHandles = new SolverAvgHandles(Q, U, R, W, T, A);
        return getAvg(customHandles);
    }

    /**
     * Computes and returns average arrival rates at steady-state.
     * If results are not available, triggers solver execution.
     *
     * @return matrix of average arrival rates [stations x classes]
     */
    public Matrix getAvgArvR() {
        if (!this.hasResults()) {
            this.getAvg();
        }
        return this.result.AN;
    }

    /**
     * Returns average arrival rates aggregated by job chains.
     *
     * @return matrix of arrival rates [stations x chains]
     */
    public Matrix getAvgArvRChain() {
        int C = sn.nchains;
        Matrix ANclass = getAvgArvR();
        Matrix AN = new Matrix(sn.nstations, C);
        for (int c = 0; c < C; c++) {
            Matrix inchain = sn.inchain.get(c);
            for (int i = 0; i < sn.nstations; i++) {
                for (int k1 = 0; k1 < inchain.length(); k1++) {
                    int k = (int) inchain.get(k1);
                    AN.set(i, c, AN.get(i, c) + ANclass.get(i, k));
                }
            }
        }
        return AN;
    }

    /**
     * Returns the average arrival rate metric handles.
     *
     * @return handles for arrival rate metrics
     */
    public AvgHandle getAvgArvRHandles() {
        return this.avgHandles.A;
    }

    /**
     * Returns average station metrics aggregated by job chains.
     *
     * @return solver result with metrics aggregated by chains
     */
    public SolverResult getAvgChain() {
        SolverResult res = new SolverResult();
        res.QN = this.getAvgQLenChain();
        res.UN = this.getAvgUtilChain();
        res.RN = this.getAvgRespTChain();
        res.TN = this.getAvgTputChain();
        res.AN = this.getAvgArvRChain();
        res.WN = this.getAvgResidTChain();
        // Preserve the method field from the original result
        if (this.result != null) {
            res.method = this.result.method;
            res.solver = this.result.solver;
            res.runtime = this.result.runtime;
            res.iter = this.result.iter;
        }
        // Note: Do NOT overwrite this.result here as it would replace class-level
        // results (M x K) with chain-level results (M x C), breaking subsequent
        // calls to getAvgSys() which expect class-level dimensions.
        return res;
    }

    /**
     * Returns a table of average station metrics aggregated by job chains.
     *
     * @return table containing station metrics organized by chains
     */
    public NetworkAvgChainTable getAvgChainTable() {

        this.sn = model.getStruct(true);

        boolean keepDisabled = false;
        this.avgHandles = model.getAvgHandles();

        int M = sn.nstations;
        int C = sn.nchains;

        SolverResult chainResult = null;
        try {
            if (!isInf(options.timespan[1]) && !Double.isNaN(options.timespan[1])) {
                this.tranHandles = model.getTranHandles();
                //getTranAvg();
            } else {
                chainResult = getAvgChain();
            }
        } catch (Exception e) {
            e.printStackTrace();
        }

        // Use the chain-aggregated results, not this.result which has class-level dimensions
        Matrix QN = chainResult != null ? chainResult.QN : null;
        Matrix UN = chainResult != null ? chainResult.UN : null;
        Matrix RN = chainResult != null ? chainResult.RN : null;
        Matrix TN = chainResult != null ? chainResult.TN : null;
        Matrix AN = chainResult != null ? chainResult.AN : null;
        Matrix WN = chainResult != null ? chainResult.WN : null;

        if (QN == null || QN.isEmpty()) {
            throw new RuntimeException(
                    "Unable to compute results and therefore unable to print AvgChainTable.");
        }

        if (!keepDisabled) {

            List<Double> Qval = new ArrayList<>();
            List<Double> Uval = new ArrayList<>();
            List<Double> Rval = new ArrayList<>();
            List<Double> Tval = new ArrayList<>();
            List<Double> ArvR = new ArrayList<>();
            List<Double> Residval = new ArrayList<>();
            List<String> className = new ArrayList<>();
            List<String> stationName = new ArrayList<>();

            for (int i = 0; i < M; i++) {
                for (int c = 0; c < C; c++) {
                    // if (QN.get(i, k) + UN.get(i, k) + RN.get(i, k) + TN.get(i, k) > 0) {
                    Qval.add(QN.get(i, c));
                    Uval.add(UN.get(i, c));
                    Rval.add(RN.get(i, c));
                    ArvR.add(AN.get(i, c));
                    Tval.add(TN.get(i, c));
                    Residval.add(WN.get(i, c));
                    int c1 = c + 1;
                    className.add("Chain" + c1);
                    stationName.add(this.model.getStations().get(i).getName());
                }
            }
            NetworkAvgChainTable avgChainTable = new NetworkAvgChainTable(Qval, Uval, Rval, Residval, ArvR, Tval);
            avgChainTable.setOptions(this.options);
            avgChainTable.setStationNames(stationName);

            java.util.List<String> chainNames = new ArrayList<>();
            java.util.List<String> inChainNames = new ArrayList<>();
            for (int i = 0; i < M; i++) {
                for (int c = 0; c < this.sn.nchains; c++) {
                    chainNames.add("Chain" + (c + 1));
                    Matrix inchain = sn.inchain.get(c);
                    String chainMembers = "(";
                    for (int j = 0; j < inchain.length(); j++) {
                        int r = (int) inchain.get(j);
                        if (j == 0) {
                            chainMembers = chainMembers.concat(sn.classnames.get(r));
                        } else {
                            chainMembers = chainMembers.concat(" " + sn.classnames.get(r));
                        }
                    }
                    inChainNames.add(chainMembers + ")");
                }
            }
            avgChainTable.setChainNames(chainNames);
            avgChainTable.setInChainNames(inChainNames);

            return avgChainTable;
        } else {
            // Keep all entries including disabled ones (keepDisabled == true)
            List<Double> Qval = new ArrayList<>();
            List<Double> Uval = new ArrayList<>();
            List<Double> Rval = new ArrayList<>();
            List<Double> Tval = new ArrayList<>();
            List<Double> ArvR = new ArrayList<>();
            List<Double> Residval = new ArrayList<>();
            List<String> className = new ArrayList<>();
            List<String> stationName = new ArrayList<>();

            for (int i = 0; i < M; i++) {
                for (int c = 0; c < C; c++) {
                    // Include all entries regardless of their values
                    Qval.add(QN.get(i, c));
                    Uval.add(UN.get(i, c));
                    Rval.add(RN.get(i, c));
                    ArvR.add(AN.get(i, c));
                    Tval.add(TN.get(i, c));
                    Residval.add(WN.get(i, c));
                    int c1 = c + 1;
                    className.add("Chain" + c1);
                    stationName.add(this.model.getStations().get(i).getName());
                }
            }
            NetworkAvgChainTable avgChainTable = new NetworkAvgChainTable(Qval, Uval, Rval, Residval, ArvR, Tval);
            avgChainTable.setOptions(this.options);
            avgChainTable.setStationNames(stationName);

            java.util.List<String> chainNames = new ArrayList<>();
            java.util.List<String> inChainNames = new ArrayList<>();
            for (int i = 0; i < M; i++) {
                for (int c = 0; c < this.sn.nchains; c++) {
                    chainNames.add("Chain" + (c + 1));
                    Matrix inchain = sn.inchain.get(c);
                    String chainMembers = "(";
                    for (int j = 0; j < inchain.length(); j++) {
                        int r = (int) inchain.get(j);
                        if (j == 0) {
                            chainMembers = chainMembers.concat(sn.classnames.get(r));
                        } else {
                            chainMembers = chainMembers.concat(" " + sn.classnames.get(r));
                        }
                    }
                    inChainNames.add(chainMembers + ")");
                }
            }
            avgChainTable.setChainNames(chainNames);
            avgChainTable.setInChainNames(inChainNames);

            return avgChainTable;
        }
    }

    /**
     * Returns a table of average station metrics aggregated by job chains using specified handles.
     *
     * @param avgHandles custom handles for performance metrics
     * @return table containing station metrics organized by chains
     */
    public NetworkAvgChainTable getAvgChainTable(SolverAvgHandles avgHandles) {
        SolverAvgHandles originalHandles = this.avgHandles;
        this.avgHandles = avgHandles;
        try {
            return getAvgChainTable();
        } finally {
            this.avgHandles = originalHandles;
        }
    }

    /**
     * Returns a table of average station metrics aggregated by job chains using individual handles.
     *
     * @param Q queue length handle
     * @param U utilization handle
     * @param R response time handle
     * @param W residence time handle
     * @param T throughput handle
     * @param A arrival rate handle
     * @return table containing station metrics organized by chains
     */
    public NetworkAvgChainTable getAvgChainTable(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A) {
        SolverAvgHandles customHandles = new SolverAvgHandles(Q, U, R, W, T, A);
        return getAvgChainTable(customHandles);
    }

    /**
     * Returns a table of average station metrics aggregated by job chains with keepDisabled option.
     *
     * @param keepDisabled whether to include disabled metrics in the table
     * @return table containing station metrics organized by chains
     */
    public NetworkAvgChainTable getAvgChainTable(boolean keepDisabled) {
        // Modify the implementation to support keepDisabled
        this.sn = model.getStruct(true);
        this.avgHandles = model.getAvgHandles();

        int M = sn.nstations;
        int C = sn.nchains;

        SolverResult chainResult = null;
        try {
            if (!isInf(options.timespan[1]) && !Double.isNaN(options.timespan[1])) {
                this.tranHandles = model.getTranHandles();
                //getTranAvg();
            } else {
                chainResult = getAvgChain();
            }
        } catch (Exception e) {
            e.printStackTrace();
        }

        // Use the chain-aggregated results, not this.result which has class-level dimensions
        Matrix QN = chainResult != null ? chainResult.QN : null;
        Matrix UN = chainResult != null ? chainResult.UN : null;
        Matrix RN = chainResult != null ? chainResult.RN : null;
        Matrix TN = chainResult != null ? chainResult.TN : null;
        Matrix AN = chainResult != null ? chainResult.AN : null;
        Matrix WN = chainResult != null ? chainResult.WN : null;

        if (QN == null || QN.isEmpty()) {
            throw new RuntimeException(
                    "Unable to compute results and therefore unable to print AvgChainTable.");
        }

        List<Double> Qval = new ArrayList<>();
        List<Double> Uval = new ArrayList<>();
        List<Double> Rval = new ArrayList<>();
        List<Double> Tval = new ArrayList<>();
        List<Double> ArvR = new ArrayList<>();
        List<Double> Residval = new ArrayList<>();
        List<String> className = new ArrayList<>();
        List<String> stationName = new ArrayList<>();

        for (int i = 0; i < M; i++) {
            for (int c = 0; c < C; c++) {
                // Include entry based on keepDisabled flag
                if (keepDisabled || (QN.get(i, c) + UN.get(i, c) + RN.get(i, c) + TN.get(i, c) + AN.get(i, c) + WN.get(i, c) > 0)) {
                    Qval.add(QN.get(i, c));
                    Uval.add(UN.get(i, c));
                    Rval.add(RN.get(i, c));
                    ArvR.add(AN.get(i, c));
                    Tval.add(TN.get(i, c));
                    Residval.add(WN.get(i, c));
                    int c1 = c + 1;
                    className.add("Chain" + c1);
                    stationName.add(this.model.getStations().get(i).getName());
                }
            }
        }
        NetworkAvgChainTable avgChainTable = new NetworkAvgChainTable(Qval, Uval, Rval, Residval, ArvR, Tval);
        avgChainTable.setOptions(this.options);
        avgChainTable.setStationNames(stationName);

        java.util.List<String> chainNames = new ArrayList<>();
        java.util.List<String> inChainNames = new ArrayList<>();
        for (int i = 0; i < M; i++) {
            for (int c = 0; c < this.sn.nchains; c++) {
                chainNames.add("Chain" + (c + 1));
                Matrix inchain = sn.inchain.get(c);
                String chainMembers = "(";
                for (int j = 0; j < inchain.length(); j++) {
                    int r = (int) inchain.get(j);
                    if (j == 0) {
                        chainMembers = chainMembers.concat(sn.classnames.get(r));
                    } else {
                        chainMembers = chainMembers.concat(" " + sn.classnames.get(r));
                    }
                }
                inChainNames.add(chainMembers + ")");
            }
        }
        avgChainTable.setChainNames(chainNames);
        avgChainTable.setInChainNames(inChainNames);

        return avgChainTable;
    }

    /**
     * Returns a table of average station metrics aggregated by job chains using specified handles and keepDisabled option.
     *
     * @param avgHandles custom handles for performance metrics
     * @param keepDisabled whether to include disabled metrics in the table
     * @return table containing station metrics organized by chains
     */
    public NetworkAvgChainTable getAvgChainTable(SolverAvgHandles avgHandles, boolean keepDisabled) {
        SolverAvgHandles originalHandles = this.avgHandles;
        this.avgHandles = avgHandles;
        try {
            return getAvgChainTable(keepDisabled);
        } finally {
            this.avgHandles = originalHandles;
        }
    }

    /**
     * Returns a table of average station metrics aggregated by job chains using individual handles and keepDisabled option.
     *
     * @param Q queue length handle
     * @param U utilization handle
     * @param R response time handle
     * @param W residence time handle
     * @param T throughput handle
     * @param A arrival rate handle
     * @param keepDisabled whether to include disabled metrics in the table
     * @return table containing station metrics organized by chains
     */
    public NetworkAvgChainTable getAvgChainTable(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A, boolean keepDisabled) {
        SolverAvgHandles customHandles = new SolverAvgHandles(Q, U, R, W, T, A);
        return getAvgChainTable(customHandles, keepDisabled);
    }

    /**
     * Alias for getAvgChainTable().
     *
     * @return table containing station metrics organized by chains
     */
    public NetworkAvgChainTable getChainAvgT() {
        return getAvgChainTable();
    }

    /**
     * Alias for getAvgChainTable(boolean keepDisabled).
     *
     * @param keepDisabled whether to include disabled metrics in the table
     * @return table containing station metrics organized by chains
     */
    public NetworkAvgChainTable getChainAvgT(boolean keepDisabled) {
        return getAvgChainTable(keepDisabled);
    }

    /**
     * Alias for getAvgChainTable(SolverAvgHandles avgHandles).
     *
     * @param avgHandles the average handles to use
     * @return table containing station metrics organized by chains
     */
    public NetworkAvgChainTable getChainAvgT(SolverAvgHandles avgHandles) {
        return getAvgChainTable(avgHandles);
    }

    /**
     * Alias for getAvgChainTable(SolverAvgHandles avgHandles, boolean keepDisabled).
     *
     * @param avgHandles the average handles to use
     * @param keepDisabled whether to include disabled metrics in the table
     * @return table containing station metrics organized by chains
     */
    public NetworkAvgChainTable getChainAvgT(SolverAvgHandles avgHandles, boolean keepDisabled) {
        return getAvgChainTable(avgHandles, keepDisabled);
    }

    /**
     * Alias for getAvgChainTable(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A).
     *
     * @param Q queue length handle
     * @param U utilization handle
     * @param R response time handle
     * @param W residence time handle
     * @param T throughput handle
     * @param A arrival rate handle
     * @return table containing station metrics organized by chains
     */
    public NetworkAvgChainTable getChainAvgT(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A) {
        return getAvgChainTable(Q, U, R, W, T, A);
    }

    /**
     * Alias for getAvgChainTable(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A, boolean keepDisabled).
     *
     * @param Q queue length handle
     * @param U utilization handle
     * @param R response time handle
     * @param W residence time handle
     * @param T throughput handle
     * @param A arrival rate handle
     * @param keepDisabled whether to include disabled metrics in the table
     * @return table containing station metrics organized by chains
     */
    public NetworkAvgChainTable getChainAvgT(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A, boolean keepDisabled) {
        return getAvgChainTable(Q, U, R, W, T, A, keepDisabled);
    }

    /**
     * Alias for getAvgChainTable().
     *
     * @return table containing station metrics organized by chains
     */
    public NetworkAvgChainTable chainAvgT() {
        return getAvgChainTable();
    }

    /**
     * Alias for getAvgChainTable(boolean keepDisabled).
     *
     * @param keepDisabled whether to include disabled metrics in the table
     * @return table containing station metrics organized by chains
     */
    public NetworkAvgChainTable chainAvgT(boolean keepDisabled) {
        return getAvgChainTable(keepDisabled);
    }

    /**
     * Alias for getAvgChainTable(SolverAvgHandles avgHandles).
     *
     * @param avgHandles the average handles to use
     * @return table containing station metrics organized by chains
     */
    public NetworkAvgChainTable chainAvgT(SolverAvgHandles avgHandles) {
        return getAvgChainTable(avgHandles);
    }

    /**
     * Alias for getAvgChainTable(SolverAvgHandles avgHandles, boolean keepDisabled).
     *
     * @param avgHandles the average handles to use
     * @param keepDisabled whether to include disabled metrics in the table
     * @return table containing station metrics organized by chains
     */
    public NetworkAvgChainTable chainAvgT(SolverAvgHandles avgHandles, boolean keepDisabled) {
        return getAvgChainTable(avgHandles, keepDisabled);
    }

    /**
     * Alias for getAvgChainTable(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A).
     *
     * @param Q queue length handle
     * @param U utilization handle
     * @param R response time handle
     * @param W residence time handle
     * @param T throughput handle
     * @param A arrival rate handle
     * @return table containing station metrics organized by chains
     */
    public NetworkAvgChainTable chainAvgT(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A) {
        return getAvgChainTable(Q, U, R, W, T, A);
    }

    /**
     * Alias for getAvgChainTable(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A, boolean keepDisabled).
     *
     * @param Q queue length handle
     * @param U utilization handle
     * @param R response time handle
     * @param W residence time handle
     * @param T throughput handle
     * @param A arrival rate handle
     * @param keepDisabled whether to include disabled metrics in the table
     * @return table containing station metrics organized by chains
     */
    public NetworkAvgChainTable chainAvgT(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A, boolean keepDisabled) {
        return getAvgChainTable(Q, U, R, W, T, A, keepDisabled);
    }

    // aCT aliases for getAvgChainTable
    /**
     * Alias for getAvgChainTable().
     *
     * @return table containing station metrics organized by chains
     */
    public NetworkAvgChainTable aCT() {
        return getAvgChainTable();
    }

    /**
     * Alias for getAvgChainTable(boolean keepDisabled).
     *
     * @param keepDisabled whether to include disabled metrics in the table
     * @return table containing station metrics organized by chains
     */
    public NetworkAvgChainTable aCT(boolean keepDisabled) {
        return getAvgChainTable(keepDisabled);
    }

    /**
     * Alias for getAvgChainTable(SolverAvgHandles avgHandles).
     *
     * @param avgHandles the average handles to use
     * @return table containing station metrics organized by chains
     */
    public NetworkAvgChainTable aCT(SolverAvgHandles avgHandles) {
        return getAvgChainTable(avgHandles);
    }

    /**
     * Alias for getAvgChainTable(SolverAvgHandles avgHandles, boolean keepDisabled).
     *
     * @param avgHandles the average handles to use
     * @param keepDisabled whether to include disabled metrics in the table
     * @return table containing station metrics organized by chains
     */
    public NetworkAvgChainTable aCT(SolverAvgHandles avgHandles, boolean keepDisabled) {
        return getAvgChainTable(avgHandles, keepDisabled);
    }

    /**
     * Alias for getAvgChainTable(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A).
     *
     * @param Q queue length handle
     * @param U utilization handle
     * @param R response time handle
     * @param W residence time handle
     * @param T throughput handle
     * @param A arrival rate handle
     * @return table containing station metrics organized by chains
     */
    public NetworkAvgChainTable aCT(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A) {
        return getAvgChainTable(Q, U, R, W, T, A);
    }

    /**
     * Alias for getAvgChainTable(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A, boolean keepDisabled).
     *
     * @param Q queue length handle
     * @param U utilization handle
     * @param R response time handle
     * @param W residence time handle
     * @param T throughput handle
     * @param A arrival rate handle
     * @param keepDisabled whether to include disabled metrics in the table
     * @return table containing station metrics organized by chains
     */
    public NetworkAvgChainTable aCT(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A, boolean keepDisabled) {
        return getAvgChainTable(Q, U, R, W, T, A, keepDisabled);
    }

    /**
     * Returns the steady-state average performance metric handles.
     *
     * @return the average handles object
     */
    public SolverAvgHandles getAvgHandles() {
        return this.avgHandles;
    }

    /**
     * Sets the steady-state average performance metric handles.
     *
     * @param handles the average handles to set
     */
    public void setAvgHandles(SolverAvgHandles handles) {
        this.avgHandles = handles;
    }


    // NOTE: the following LINE methods have not been migrated to JLINE
    // a) updateModel() - model is public and therefore no need for setter
    // b) all methods that are "not supported by this solver" - lack of existing method will suffice
    // rather than dedicated warning that method is not applicable

    /**
     * Computes average performance metrics at steady-state for all nodes.
     * This method aggregates station-level metrics to node-level metrics.
     *
     * @return solver result containing node-level average metrics
     */
    public SolverResult getAvgNode() {
        if (this.avgHandles == null || this.avgHandles.Q == null || this.avgHandles.U == null || this.avgHandles.R == null ||
                this.avgHandles.T == null || this.avgHandles.A == null) {
            reset();
        }
        this.avgHandles = model.getAvgHandles();
        SolverResult result = getAvg();

        SolverResult noderesult = new SolverResult();
        int I = sn.nnodes;  // Physical nodes only
        int M = sn.nstations;
        int R = sn.nclasses;
        int C = sn.nchains;
        int F = sn.nregions;  // Number of FCR virtual nodes

        // Total nodes includes physical nodes + FCR virtual nodes
        int totalNodes = I + F;

        Matrix QNn = Matrix.zeros(totalNodes, R);
        Matrix UNn = Matrix.zeros(totalNodes, R);
        Matrix RNn = Matrix.zeros(totalNodes, R);
        Matrix WNn = Matrix.zeros(totalNodes, R);
        for (int ist = 0; ist < M; ist++) {
            int ind = (int) sn.stationToNode.get(ist);
            for (int r = 0; r < R; r++) {
                QNn.set(ind, r, result.QN.get(ist, r));
                UNn.set(ind, r, result.UN.get(ist, r));
                RNn.set(ind, r, result.RN.get(ist, r));
                WNn.set(ind, r, result.WN.get(ist, r));
            }
        }

        AvgHandle T = getAvgTputHandles();
        Matrix ANn = snGetNodeArvRFromTput(sn, result.TN, T, result.AN);
        Matrix TNn = snGetNodeTputFromTput(sn, result.TN, T, ANn);

        // Extend ANn and TNn to include FCR rows
        if (F > 0) {
            Matrix ANnExtended = Matrix.zeros(totalNodes, R);
            Matrix TNnExtended = Matrix.zeros(totalNodes, R);
            for (int i = 0; i < I; i++) {
                for (int r = 0; r < R; r++) {
                    ANnExtended.set(i, r, ANn.get(i, r));
                    TNnExtended.set(i, r, TNn.get(i, r));
                }
            }
            ANn = ANnExtended;
            TNn = TNnExtended;
        }

        // Merge FCR metrics into node matrices at FCR node indices
        if (result.QNfcr != null && F > 0) {
            for (int f = 0; f < F; f++) {
                int fcrNodeIdx = I + f;  // FCR indices start after physical nodes
                for (int r = 0; r < R; r++) {
                    QNn.set(fcrNodeIdx, r, result.QNfcr.get(f, r));
                    UNn.set(fcrNodeIdx, r, result.UNfcr.get(f, r));
                    RNn.set(fcrNodeIdx, r, result.RNfcr.get(f, r));
                    WNn.set(fcrNodeIdx, r, result.WNfcr.get(f, r));
                    TNn.set(fcrNodeIdx, r, result.TNfcr.get(f, r));
                    ANn.set(fcrNodeIdx, r, result.ANfcr.get(f, r));
                }
            }
        }

        noderesult.QN = QNn;
        noderesult.UN = UNn;
        noderesult.RN = RNn;
        noderesult.TN = TNn;
        noderesult.AN = ANn;
        noderesult.WN = WNn;
        return noderesult;
    }

    // ========== State Probability Methods (Abstract - to be implemented by subclasses) ==========

    /**
     * Returns average node arrival rates aggregated by job chains.
     *
     * @return matrix of node arrival rates [nodes x chains]
     */
    public Matrix getAvgNodeArvRChain() {
        int C = sn.nchains;
        Matrix ANclass = getAvgNode().AN;
        Matrix AN = new Matrix(sn.nnodes, C);
        for (int c = 0; c < C; c++) {
            Matrix inchain = sn.inchain.get(c);
            for (int i = 0; i < sn.nnodes; i++) {
                for (int k1 = 0; k1 < inchain.length(); k1++) {
                    int k = (int) inchain.get(k1);
                    AN.set(i, c, AN.get(i, c) + ANclass.get(i, k));
                }
            }
        }
        return AN;
    }

    /**
     * Returns average node metrics aggregated by job chains.
     *
     * @return solver result with node metrics aggregated by chains
     */
    public SolverResult getAvgNodeChain() {
        SolverResult res = new SolverResult();
        res.QN = this.getAvgNodeQLenChain();
        res.UN = this.getAvgNodeUtilChain();
        res.RN = this.getAvgNodeRespTChain();
        res.TN = this.getAvgNodeTputChain();
        res.AN = this.getAvgNodeArvRChain();
        res.WN = this.getAvgNodeResidTChain();
        this.result = res;
        return res;
    }

    /**
     * Returns a table of average node metrics aggregated by chains.
     *
     * @return table containing node-level metrics organized by job chains
     */
    public NetworkAvgNodeChainTable getAvgNodeChainTable() {

        this.sn = model.getStruct(true);

        boolean keepDisabled = false;
        this.avgHandles = model.getAvgHandles();

        int I = sn.nnodes;
        int M = sn.nstations;
        int C = sn.nchains;

        try {
            if (!isInf(options.timespan[1]) && !Double.isNaN(options.timespan[1])) {
                this.tranHandles = model.getTranHandles();
                //getTranAvg();
            } else {
                getAvgNodeChain();
            }
        } catch (Exception e) {
            e.printStackTrace();
        }

        Matrix QN = this.result.QN;
        Matrix UN = this.result.UN;
        Matrix RN = this.result.RN;
        Matrix TN = this.result.TN;
        Matrix AN = this.result.AN;
        Matrix WN = this.result.WN;

        if (QN == null || QN.isEmpty()) {
            throw new RuntimeException(
                    "Unable to compute results and therefore unable to print AvgTable.");
        }

        if (!keepDisabled) {

            List<Double> Qval = new ArrayList<>();
            List<Double> Uval = new ArrayList<>();
            List<Double> Rval = new ArrayList<>();
            List<Double> Tval = new ArrayList<>();
            List<Double> ArvR = new ArrayList<>();
            List<Double> Residval = new ArrayList<>();
            List<String> className = new ArrayList<>();
            List<String> nodeName = new ArrayList<>();

            for (int i = 0; i < I; i++) {
                for (int c = 0; c < C; c++) {
                    //if (QN.get(i, c) + UN.get(i, c) + RN.get(i, c) + TN.get(i, c) > 0) {
                    Qval.add(QN.get(i, c));
                    Uval.add(UN.get(i, c));
                    Rval.add(RN.get(i, c));
                    ArvR.add(AN.get(i, c));
                    Tval.add(TN.get(i, c));
                    Residval.add(WN.get(i, c));
                    int c1 = c + 1;
                    className.add("Chain" + c1);
                    nodeName.add(this.model.getNodes().get(i).getName());
                    //}
                }
            }
            NetworkAvgNodeChainTable avgNodeChainTable = new NetworkAvgNodeChainTable(Qval, Uval, Rval, Residval, ArvR, Tval);
            avgNodeChainTable.setOptions(this.options);
            avgNodeChainTable.setNodeNames(nodeName);

            java.util.List<String> chainNames = new ArrayList<>();
            java.util.List<String> inChainNames = new ArrayList<>();
            for (int i = 0; i < I; i++) {
                for (int c = 0; c < this.sn.nchains; c++) {
                    chainNames.add("Chain" + (c + 1));
                    Matrix inchain = sn.inchain.get(c);
                    String chainMembers = "(";
                    for (int j = 0; j < inchain.length(); j++) {
                        int r = (int) inchain.get(j);
                        if (j == 0) {
                            chainMembers = chainMembers.concat(sn.classnames.get(r));
                        } else {
                            chainMembers = chainMembers.concat(" " + sn.classnames.get(r));
                        }
                    }
                    inChainNames.add(chainMembers + ")");
                }
            }
            avgNodeChainTable.setChainNames(chainNames);
            avgNodeChainTable.setInChainNames(inChainNames);

            return avgNodeChainTable;
        } else {
            // Keep all entries including disabled ones (keepDisabled == true)
            List<Double> Qval = new ArrayList<>();
            List<Double> Uval = new ArrayList<>();
            List<Double> Rval = new ArrayList<>();
            List<Double> Tval = new ArrayList<>();
            List<Double> ArvR = new ArrayList<>();
            List<Double> Residval = new ArrayList<>();
            List<String> className = new ArrayList<>();
            List<String> nodeName = new ArrayList<>();

            for (int i = 0; i < M; i++) {
                for (int c = 0; c < C; c++) {
                    // Include all entries regardless of their values
                    Qval.add(QN.get(i, c));
                    Uval.add(UN.get(i, c));
                    Rval.add(RN.get(i, c));
                    ArvR.add(AN.get(i, c));
                    Tval.add(TN.get(i, c));
                    Residval.add(WN.get(i, c));
                    int c1 = c + 1;
                    className.add("Chain" + c1);
                    nodeName.add(this.model.getNodes().get(i).getName());
                }
            }
            NetworkAvgNodeChainTable avgNodeChainTable = new NetworkAvgNodeChainTable(Qval, Uval, Rval, Residval, ArvR, Tval);
            avgNodeChainTable.setOptions(this.options);
            avgNodeChainTable.setNodeNames(nodeName);

            java.util.List<String> chainNames = new ArrayList<>();
            java.util.List<String> inChainNames = new ArrayList<>();
            for (int i = 0; i < M; i++) {
                for (int c = 0; c < this.sn.nchains; c++) {
                    chainNames.add("Chain" + (c + 1));
                    Matrix inchain = sn.inchain.get(c);
                    String chainMembers = "(";
                    for (int j = 0; j < inchain.length(); j++) {
                        int r = (int) inchain.get(j);
                        if (j == 0) {
                            chainMembers = chainMembers.concat(sn.classnames.get(r));
                        } else {
                            chainMembers = chainMembers.concat(" " + sn.classnames.get(r));
                        }
                    }
                    inChainNames.add(chainMembers + ")");
                }
            }
            avgNodeChainTable.setChainNames(chainNames);
            avgNodeChainTable.setInChainNames(inChainNames);

            return avgNodeChainTable;
        }
    }

    /**
     * Returns a table of average node metrics aggregated by chains using specified handles.
     *
     * @param avgHandles custom handles for performance metrics
     * @return table containing node-level metrics organized by job chains
     */
    public NetworkAvgNodeChainTable getAvgNodeChainTable(SolverAvgHandles avgHandles) {
        SolverAvgHandles originalHandles = this.avgHandles;
        this.avgHandles = avgHandles;
        try {
            return getAvgNodeChainTable();
        } finally {
            this.avgHandles = originalHandles;
        }
    }

    /**
     * Returns a table of average node metrics aggregated by chains using individual handles.
     *
     * @param Q queue length handle
     * @param U utilization handle
     * @param R response time handle
     * @param W residence time handle
     * @param T throughput handle
     * @param A arrival rate handle
     * @return table containing node-level metrics organized by job chains
     */
    public NetworkAvgNodeChainTable getAvgNodeChainTable(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A) {
        SolverAvgHandles customHandles = new SolverAvgHandles(Q, U, R, W, T, A);
        return getAvgNodeChainTable(customHandles);
    }

    /**
     * Returns a table of average node metrics aggregated by chains with keepDisabled option.
     *
     * @param keepDisabled whether to include disabled metrics in the table
     * @return table containing node-level metrics organized by job chains
     */
    public NetworkAvgNodeChainTable getAvgNodeChainTable(boolean keepDisabled) {
        // Implement keepDisabled functionality
        this.sn = model.getStruct(true);
        this.avgHandles = model.getAvgHandles();

        int M = sn.nnodes;
        int C = sn.nchains;

        try {
            if (!isInf(options.timespan[1]) && !Double.isNaN(options.timespan[1])) {
                this.tranHandles = model.getTranHandles();
                //getTranAvg();
            } else {
                getAvgNodeChain();
            }
        } catch (Exception e) {
            e.printStackTrace();
        }

        Matrix QN = this.result.QN;
        Matrix UN = this.result.UN;
        Matrix RN = this.result.RN;
        Matrix TN = this.result.TN;
        Matrix AN = this.result.AN;
        Matrix WN = this.result.WN;

        if (QN == null || QN.isEmpty()) {
            throw new RuntimeException(
                    "Unable to compute results and therefore unable to print AvgTable.");
        }

        List<Double> Qval = new ArrayList<>();
        List<Double> Uval = new ArrayList<>();
        List<Double> Rval = new ArrayList<>();
        List<Double> Tval = new ArrayList<>();
        List<Double> ArvR = new ArrayList<>();
        List<Double> Residval = new ArrayList<>();
        List<String> className = new ArrayList<>();
        List<String> nodeName = new ArrayList<>();

        for (int i = 0; i < M; i++) {
            for (int c = 0; c < C; c++) {
                // Include entry based on keepDisabled flag
                if (keepDisabled || (QN.get(i, c) + UN.get(i, c) + RN.get(i, c) + TN.get(i, c) + AN.get(i, c) + WN.get(i, c) > 0)) {
                    Qval.add(QN.get(i, c));
                    Uval.add(UN.get(i, c));
                    Rval.add(RN.get(i, c));
                    ArvR.add(AN.get(i, c));
                    Tval.add(TN.get(i, c));
                    Residval.add(WN.get(i, c));
                    int c1 = c + 1;
                    className.add("Chain" + c1);
                    nodeName.add(this.model.getNodes().get(i).getName());
                }
            }
        }
        NetworkAvgNodeChainTable avgNodeChainTable = new NetworkAvgNodeChainTable(Qval, Uval, Rval, Residval, ArvR, Tval);
        avgNodeChainTable.setOptions(this.options);
        avgNodeChainTable.setNodeNames(nodeName);

        java.util.List<String> chainNames = new ArrayList<>();
        java.util.List<String> inChainNames = new ArrayList<>();
        for (int i = 0; i < M; i++) {
            for (int c = 0; c < this.sn.nchains; c++) {
                chainNames.add("Chain" + (c + 1));
                Matrix inchain = sn.inchain.get(c);
                String chainMembers = "(";
                for (int j = 0; j < inchain.length(); j++) {
                    int r = (int) inchain.get(j);
                    if (j == 0) {
                        chainMembers = chainMembers.concat(sn.classnames.get(r));
                    } else {
                        chainMembers = chainMembers.concat(" " + sn.classnames.get(r));
                    }
                }
                inChainNames.add(chainMembers + ")");
            }
        }
        avgNodeChainTable.setChainNames(chainNames);
        avgNodeChainTable.setInChainNames(inChainNames);

        return avgNodeChainTable;
    }

    /**
     * Returns a table of average node metrics aggregated by chains using specified handles and keepDisabled option.
     *
     * @param avgHandles custom handles for performance metrics
     * @param keepDisabled whether to include disabled metrics in the table
     * @return table containing node-level metrics organized by job chains
     */
    public NetworkAvgNodeChainTable getAvgNodeChainTable(SolverAvgHandles avgHandles, boolean keepDisabled) {
        SolverAvgHandles originalHandles = this.avgHandles;
        this.avgHandles = avgHandles;
        try {
            return getAvgNodeChainTable(keepDisabled);
        } finally {
            this.avgHandles = originalHandles;
        }
    }

    /**
     * Returns a table of average node metrics aggregated by chains using individual handles and keepDisabled option.
     *
     * @param Q queue length handle
     * @param U utilization handle
     * @param R response time handle
     * @param W residence time handle
     * @param T throughput handle
     * @param A arrival rate handle
     * @param keepDisabled whether to include disabled metrics in the table
     * @return table containing node-level metrics organized by job chains
     */
    public NetworkAvgNodeChainTable getAvgNodeChainTable(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A, boolean keepDisabled) {
        SolverAvgHandles customHandles = new SolverAvgHandles(Q, U, R, W, T, A);
        return getAvgNodeChainTable(customHandles, keepDisabled);
    }

    /**
     * Alias for getAvgNodeChainTable().
     *
     * @return table containing node-level metrics organized by job chains
     */
    public NetworkAvgNodeChainTable getNodeChainAvgT() {
        return getAvgNodeChainTable();
    }

    /**
     * Alias for getAvgNodeChainTable(boolean keepDisabled).
     *
     * @param keepDisabled whether to include disabled metrics in the table
     * @return table containing node-level metrics organized by job chains
     */
    public NetworkAvgNodeChainTable getNodeChainAvgT(boolean keepDisabled) {
        return getAvgNodeChainTable(keepDisabled);
    }

    /**
     * Alias for getAvgNodeChainTable(SolverAvgHandles avgHandles).
     *
     * @param avgHandles the average handles to use
     * @return table containing node-level metrics organized by job chains
     */
    public NetworkAvgNodeChainTable getNodeChainAvgT(SolverAvgHandles avgHandles) {
        return getAvgNodeChainTable(avgHandles);
    }

    /**
     * Alias for getAvgNodeChainTable(SolverAvgHandles avgHandles, boolean keepDisabled).
     *
     * @param avgHandles the average handles to use
     * @param keepDisabled whether to include disabled metrics in the table
     * @return table containing node-level metrics organized by job chains
     */
    public NetworkAvgNodeChainTable getNodeChainAvgT(SolverAvgHandles avgHandles, boolean keepDisabled) {
        return getAvgNodeChainTable(avgHandles, keepDisabled);
    }

    /**
     * Alias for getAvgNodeChainTable(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A).
     *
     * @param Q queue length handle
     * @param U utilization handle
     * @param R response time handle
     * @param W residence time handle
     * @param T throughput handle
     * @param A arrival rate handle
     * @return table containing node-level metrics organized by job chains
     */
    public NetworkAvgNodeChainTable getNodeChainAvgT(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A) {
        return getAvgNodeChainTable(Q, U, R, W, T, A);
    }

    /**
     * Alias for getAvgNodeChainTable(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A, boolean keepDisabled).
     *
     * @param Q queue length handle
     * @param U utilization handle
     * @param R response time handle
     * @param W residence time handle
     * @param T throughput handle
     * @param A arrival rate handle
     * @param keepDisabled whether to include disabled metrics in the table
     * @return table containing node-level metrics organized by job chains
     */
    public NetworkAvgNodeChainTable getNodeChainAvgT(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A, boolean keepDisabled) {
        return getAvgNodeChainTable(Q, U, R, W, T, A, keepDisabled);
    }

    /**
     * Alias for getAvgNodeChainTable().
     *
     * @return table containing node-level metrics organized by job chains
     */
    public NetworkAvgNodeChainTable nodeChainAvgT() {
        return getAvgNodeChainTable();
    }

    /**
     * Alias for getAvgNodeChainTable(boolean keepDisabled).
     *
     * @param keepDisabled whether to include disabled metrics in the table
     * @return table containing node-level metrics organized by job chains
     */
    public NetworkAvgNodeChainTable nodeChainAvgT(boolean keepDisabled) {
        return getAvgNodeChainTable(keepDisabled);
    }

    /**
     * Alias for getAvgNodeChainTable(SolverAvgHandles avgHandles).
     *
     * @param avgHandles the average handles to use
     * @return table containing node-level metrics organized by job chains
     */
    public NetworkAvgNodeChainTable nodeChainAvgT(SolverAvgHandles avgHandles) {
        return getAvgNodeChainTable(avgHandles);
    }

    /**
     * Alias for getAvgNodeChainTable(SolverAvgHandles avgHandles, boolean keepDisabled).
     *
     * @param avgHandles the average handles to use
     * @param keepDisabled whether to include disabled metrics in the table
     * @return table containing node-level metrics organized by job chains
     */
    public NetworkAvgNodeChainTable nodeChainAvgT(SolverAvgHandles avgHandles, boolean keepDisabled) {
        return getAvgNodeChainTable(avgHandles, keepDisabled);
    }

    /**
     * Alias for getAvgNodeChainTable(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A).
     *
     * @param Q queue length handle
     * @param U utilization handle
     * @param R response time handle
     * @param W residence time handle
     * @param T throughput handle
     * @param A arrival rate handle
     * @return table containing node-level metrics organized by job chains
     */
    public NetworkAvgNodeChainTable nodeChainAvgT(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A) {
        return getAvgNodeChainTable(Q, U, R, W, T, A);
    }

    /**
     * Alias for getAvgNodeChainTable(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A, boolean keepDisabled).
     *
     * @param Q queue length handle
     * @param U utilization handle
     * @param R response time handle
     * @param W residence time handle
     * @param T throughput handle
     * @param A arrival rate handle
     * @param keepDisabled whether to include disabled metrics in the table
     * @return table containing node-level metrics organized by job chains
     */
    public NetworkAvgNodeChainTable nodeChainAvgT(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A, boolean keepDisabled) {
        return getAvgNodeChainTable(Q, U, R, W, T, A, keepDisabled);
    }

    /**
     * Short alias for getAvgNodeChainTable().
     *
     * @return table containing node-level metrics organized by job chains
     */
    public NetworkAvgNodeChainTable aNCT() {
        return getAvgNodeChainTable();
    }

    /**
     * Short alias for getAvgNodeChainTable(boolean keepDisabled).
     *
     * @param keepDisabled whether to include disabled metrics in the table
     * @return table containing node-level metrics organized by job chains
     */
    public NetworkAvgNodeChainTable aNCT(boolean keepDisabled) {
        return getAvgNodeChainTable(keepDisabled);
    }

    /**
     * Short alias for getAvgNodeChainTable(SolverAvgHandles avgHandles).
     *
     * @param avgHandles the average handles to use
     * @return table containing node-level metrics organized by job chains
     */
    public NetworkAvgNodeChainTable aNCT(SolverAvgHandles avgHandles) {
        return getAvgNodeChainTable(avgHandles);
    }

    /**
     * Short alias for getAvgNodeChainTable(SolverAvgHandles avgHandles, boolean keepDisabled).
     *
     * @param avgHandles the average handles to use
     * @param keepDisabled whether to include disabled metrics in the table
     * @return table containing node-level metrics organized by job chains
     */
    public NetworkAvgNodeChainTable aNCT(SolverAvgHandles avgHandles, boolean keepDisabled) {
        return getAvgNodeChainTable(avgHandles, keepDisabled);
    }

    /**
     * Short alias for getAvgNodeChainTable(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A).
     *
     * @param Q queue length handle
     * @param U utilization handle
     * @param R response time handle
     * @param W residence time handle
     * @param T throughput handle
     * @param A arrival rate handle
     * @return table containing node-level metrics organized by job chains
     */
    public NetworkAvgNodeChainTable aNCT(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A) {
        return getAvgNodeChainTable(Q, U, R, W, T, A);
    }

    /**
     * Short alias for getAvgNodeChainTable(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A, boolean keepDisabled).
     *
     * @param Q queue length handle
     * @param U utilization handle
     * @param R response time handle
     * @param W residence time handle
     * @param T throughput handle
     * @param A arrival rate handle
     * @param keepDisabled whether to include disabled metrics in the table
     * @return table containing node-level metrics organized by job chains
     */
    public NetworkAvgNodeChainTable aNCT(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A, boolean keepDisabled) {
        return getAvgNodeChainTable(Q, U, R, W, T, A, keepDisabled);
    }

    /**
     * Returns average node queue lengths aggregated by job chains.
     *
     * @return matrix of node queue lengths [nodes x chains]
     */
    public Matrix getAvgNodeQLenChain() {
        int C = sn.nchains;
        Matrix QNclass = getAvgNode().QN;
        Matrix QN = new Matrix(sn.nnodes, C);
        for (int c = 0; c < C; c++) {
            Matrix inchain = sn.inchain.get(c);
            for (int i = 0; i < sn.nnodes; i++) {
                for (int k1 = 0; k1 < inchain.length(); k1++) {
                    int k = (int) inchain.get(k1);
                    QN.set(i, c, QN.get(i, c) + QNclass.get(i, k));
                }
            }
        }
        return QN;
    }

    // ========== Sampling Methods (Abstract - to be implemented by subclasses) ==========

    /**
     * Returns average node residence times aggregated by job chains.
     *
     * @return matrix of node residence times [nodes x chains]
     */
    public Matrix getAvgNodeResidTChain() {
        int C = sn.nchains;
        Matrix WNclass = getAvgNode().WN;
        Matrix WN = new Matrix(sn.nnodes, C);
        Matrix alpha = snGetDemandsChain(sn).alpha;
        for (int c = 0; c < C; c++) {
            Matrix inchain = sn.inchain.get(c);
            for (int i = 0; i < sn.nnodes; i++) {
                for (int k1 = 0; k1 < inchain.length(); k1++) {
                    int k = (int) inchain.get(k1);
                    if (sn.isstation.get(i) == 1)
                        WN.set(i, c, WN.get(i, c) + WNclass.get(i, k) * alpha.get(i, k));
                }
            }
        }
        return WN;
    }

    /**
     * Returns average node response times aggregated by job chains.
     *
     * @return matrix of node response times [nodes x chains]
     */
    public Matrix getAvgNodeRespTChain() {
        int C = sn.nchains;
        Matrix RNclass = getAvgNode().RN;
        Matrix RN = new Matrix(sn.nnodes, C);
        Matrix alpha = snGetDemandsChain(sn).alpha;
        for (int c = 0; c < C; c++) {
            Matrix inchain = sn.inchain.get(c);
            for (int i = 0; i < sn.nnodes; i++) {
                for (int k1 = 0; k1 < inchain.length(); k1++) {
                    int k = (int) inchain.get(k1);
                    if (sn.isstation.get(i) == 1)
                        RN.set(i, c, RN.get(i, c) + RNclass.get(i, k) * alpha.get(i, k));
                }
            }
        }
        return RN;
    }

    /**
     * Returns a table of average node metrics organized by job classes.
     *
     * @return table containing node-level metrics for each class
     */
    public NetworkAvgNodeTable getAvgNodeTable() {

        this.sn = model.getStruct(true);

        boolean keepDisabled = false;
        this.avgHandles = model.getAvgHandles();

        // Total nodes includes physical nodes + FCR virtual nodes
        int I = sn.nnodes + sn.nregions;
        int K = sn.nclasses;

        SolverResult noderesult = null;
        try {
            if (!isInf(options.timespan[1]) && !Double.isNaN(options.timespan[1])) {
                this.tranHandles = model.getTranHandles();
                //getTranAvg();
            } else {
                noderesult = getAvgNode();
            }
        } catch (Exception e) {
            e.printStackTrace();
        }

        // Check if noderesult is null before accessing its fields
        if (noderesult == null) {
            throw new RuntimeException("Unable to compute results - solver execution failed or did not complete.");
        }

        Matrix QN = noderesult.QN;
        Matrix UN = noderesult.UN;
        Matrix RN = noderesult.RN;
        Matrix WN = noderesult.WN;
        Matrix TN = noderesult.TN;
        Matrix AN = noderesult.AN;

        if (QN == null || QN.isEmpty()) {
            throw new RuntimeException(
                    "Unable to compute results and therefore unable to print AvgTable.");
        }

        if (!keepDisabled) {
            Matrix V = new Matrix(sn.nnodes, K);
            for (int i = 0; i < sn.nodevisits.size(); i++) {
                V = V.add(1, sn.nodevisits.get(i));
            }
            if (V.isEmpty()) { // SSA
                // Implementation for SSA: sum all chain visits across all chains
                for (Integer chainIndex : sn.nodevisits.keySet()) {
                    Matrix chainVisits = sn.nodevisits.get(chainIndex);
                    if (chainVisits != null) {
                        V = V.add(1.0, chainVisits);
                    }
                }
            }

            List<Double> Qval = new ArrayList<>();
            List<Double> Uval = new ArrayList<>();
            List<Double> Rval = new ArrayList<>();
            List<Double> Tval = new ArrayList<>();
            List<Double> ArvR = new ArrayList<>();
            List<Double> Residval = new ArrayList<>();
            List<String> className = new ArrayList<>();
            List<String> nodeName = new ArrayList<>();
            for (int i = 0; i < I; i++) {
                for (int k = 0; k < K; k++) {
                    int c = -1;
                    for (int row = 0; row < sn.chains.getNumRows(); row++) {
                        if (sn.chains.get(row, k) > 0) {
                            c = row;
                            break;
                        }
                    }
                    Qval.add(QN.get(i, k));
                    Uval.add(UN.get(i, k));
                    Rval.add(RN.get(i, k));
                    ArvR.add(AN.get(i, k));
                    Tval.add(TN.get(i, k));
                    className.add(model.getClasses().get(k).getName());
                    nodeName.add(sn.nodenames.get(i));
                    Residval.add(WN.get(i, k));
                }
            }
            NetworkAvgNodeTable avgTable = new NetworkAvgNodeTable(Qval, Uval, Rval, Residval, ArvR, Tval);
            avgTable.setOptions(this.options);
            avgTable.setClassNames(className);
            avgTable.setNodeNames(nodeName);

            return avgTable;
        } else {
            // Keep all entries including disabled ones (keepDisabled == true)
            Matrix V = new Matrix(sn.nnodes, K);
            for (int i = 0; i < sn.nodevisits.size(); i++) {
                V = V.add(1, sn.nodevisits.get(i));
            }
            if (V.isEmpty()) { // SSA
                // Implementation for SSA: sum all chain visits across all chains
                for (Integer chainIndex : sn.nodevisits.keySet()) {
                    Matrix chainVisits = sn.nodevisits.get(chainIndex);
                    if (chainVisits != null) {
                        V = V.add(1.0, chainVisits);
                    }
                }
            }

            List<Double> Qval = new ArrayList<>();
            List<Double> Uval = new ArrayList<>();
            List<Double> Rval = new ArrayList<>();
            List<Double> Tval = new ArrayList<>();
            List<Double> ArvR = new ArrayList<>();
            List<Double> Residval = new ArrayList<>();
            List<String> className = new ArrayList<>();
            List<String> nodeName = new ArrayList<>();
            for (int i = 0; i < I; i++) {
                for (int k = 0; k < K; k++) {
                    // Include all entries regardless of their values
                    Qval.add(QN.get(i, k));
                    Uval.add(UN.get(i, k));
                    Rval.add(RN.get(i, k));
                    ArvR.add(AN.get(i, k));
                    Tval.add(TN.get(i, k));
                    className.add(model.getClasses().get(k).getName());
                    nodeName.add(sn.nodenames.get(i));
                    Residval.add(WN.get(i, k));
                }
            }
            NetworkAvgNodeTable avgTable = new NetworkAvgNodeTable(Qval, Uval, Rval, Residval, ArvR, Tval);
            avgTable.setOptions(this.options);
            avgTable.setClassNames(className);
            avgTable.setNodeNames(nodeName);

            return avgTable;
        }
    }

    /**
     * Returns a table of average node metrics organized by job classes using specified handles.
     *
     * @param avgHandles custom handles for performance metrics
     * @return table containing node-level metrics for each class
     */
    public NetworkAvgNodeTable getAvgNodeTable(SolverAvgHandles avgHandles) {
        SolverAvgHandles originalHandles = this.avgHandles;
        this.avgHandles = avgHandles;
        try {
            return getAvgNodeTable();
        } finally {
            this.avgHandles = originalHandles;
        }
    }

    // ========== Distribution Methods (Abstract - to be implemented by subclasses) ==========

    /**
     * Returns a table of average node metrics organized by job classes using individual handles.
     *
     * @param Q queue length handle
     * @param U utilization handle
     * @param R response time handle
     * @param W residence time handle
     * @param T throughput handle
     * @param A arrival rate handle
     * @return table containing node-level metrics for each class
     */
    public NetworkAvgNodeTable getAvgNodeTable(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A) {
        SolverAvgHandles customHandles = new SolverAvgHandles(Q, U, R, W, T, A);
        return getAvgNodeTable(customHandles);
    }

    /**
     * Returns a table of average node metrics organized by job classes with keepDisabled option.
     *
     * @param keepDisabled whether to include disabled metrics in the table
     * @return table containing node-level metrics for each class
     */
    public NetworkAvgNodeTable getAvgNodeTable(boolean keepDisabled) {
        // Implement keepDisabled functionality
        this.sn = model.getStruct(true);
        this.avgHandles = model.getAvgHandles();

        // Total nodes includes physical nodes + FCR virtual nodes
        int I = sn.nnodes + sn.nregions;
        int K = sn.nclasses;

        SolverResult noderesult = null;
        try {
            if (!isInf(options.timespan[1]) && !Double.isNaN(options.timespan[1])) {
                this.tranHandles = model.getTranHandles();
                //getTranAvg();
            } else {
                noderesult = getAvgNode();
            }
        } catch (Exception e) {
            e.printStackTrace();
        }

        if (noderesult == null) {
            throw new RuntimeException("Unable to compute results - solver execution failed or did not complete.");
        }

        Matrix QN = noderesult.QN;
        Matrix UN = noderesult.UN;
        Matrix RN = noderesult.RN;
        Matrix WN = noderesult.WN;
        Matrix TN = noderesult.TN;
        Matrix AN = noderesult.AN;

        if (QN == null || QN.isEmpty()) {
            throw new RuntimeException(
                    "Unable to compute results and therefore unable to print AvgTable.");
        }

        Matrix V = new Matrix(sn.nnodes, K);
        for (int i = 0; i < sn.nodevisits.size(); i++) {
            V = V.add(1, sn.nodevisits.get(i));
        }
        if (V.isEmpty()) { // SSA
            // Implementation for SSA: sum all chain visits across all chains
            for (Integer chainIndex : sn.nodevisits.keySet()) {
                Matrix chainVisits = sn.nodevisits.get(chainIndex);
                if (chainVisits != null) {
                    V = V.add(1.0, chainVisits);
                }
            }
        }

        List<Double> Qval = new ArrayList<>();
        List<Double> Uval = new ArrayList<>();
        List<Double> Rval = new ArrayList<>();
        List<Double> Tval = new ArrayList<>();
        List<Double> ArvR = new ArrayList<>();
        List<Double> Residval = new ArrayList<>();
        List<String> className = new ArrayList<>();
        List<String> nodeName = new ArrayList<>();
        for (int i = 0; i < I; i++) {
            for (int k = 0; k < K; k++) {
                // Include entry based on keepDisabled flag
                if (keepDisabled || (QN.get(i, k) + UN.get(i, k) + RN.get(i, k) + TN.get(i, k) + AN.get(i, k) + WN.get(i, k) > 0)) {
                    Qval.add(QN.get(i, k));
                    Uval.add(UN.get(i, k));
                    Rval.add(RN.get(i, k));
                    ArvR.add(AN.get(i, k));
                    Tval.add(TN.get(i, k));
                    className.add(model.getClasses().get(k).getName());
                    nodeName.add(sn.nodenames.get(i));
                    Residval.add(WN.get(i, k));
                }
            }
        }
        NetworkAvgNodeTable avgTable = new NetworkAvgNodeTable(Qval, Uval, Rval, Residval, ArvR, Tval);
        avgTable.setOptions(this.options);
        avgTable.setClassNames(className);
        avgTable.setNodeNames(nodeName);

        return avgTable;
    }

    /**
     * Returns a table of average node metrics organized by job classes using specified handles and keepDisabled option.
     *
     * @param avgHandles custom handles for performance metrics
     * @param keepDisabled whether to include disabled metrics in the table
     * @return table containing node-level metrics for each class
     */
    public NetworkAvgNodeTable getAvgNodeTable(SolverAvgHandles avgHandles, boolean keepDisabled) {
        SolverAvgHandles originalHandles = this.avgHandles;
        this.avgHandles = avgHandles;
        try {
            return getAvgNodeTable(keepDisabled);
        } finally {
            this.avgHandles = originalHandles;
        }
    }

    /**
     * Returns a table of average node metrics organized by job classes using individual handles and keepDisabled option.
     *
     * @param Q queue length handle
     * @param U utilization handle
     * @param R response time handle
     * @param W residence time handle
     * @param T throughput handle
     * @param A arrival rate handle
     * @param keepDisabled whether to include disabled metrics in the table
     * @return table containing node-level metrics for each class
     */
    public NetworkAvgNodeTable getAvgNodeTable(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A, boolean keepDisabled) {
        SolverAvgHandles customHandles = new SolverAvgHandles(Q, U, R, W, T, A);
        return getAvgNodeTable(customHandles, keepDisabled);
    }

    /**
     * Alias for getAvgNodeTable().
     *
     * @return table containing node-level metrics for each class
     */
    public NetworkAvgNodeTable getNodeAvgT() {
        return getAvgNodeTable();
    }

    /**
     * Alias for getAvgNodeTable(boolean keepDisabled).
     *
     * @param keepDisabled whether to include disabled metrics in the table
     * @return table containing node-level metrics for each class
     */
    public NetworkAvgNodeTable getNodeAvgT(boolean keepDisabled) {
        return getAvgNodeTable(keepDisabled);
    }

    /**
     * Alias for getAvgNodeTable(SolverAvgHandles avgHandles).
     *
     * @param avgHandles the average handles to use
     * @return table containing node-level metrics for each class
     */
    public NetworkAvgNodeTable getNodeAvgT(SolverAvgHandles avgHandles) {
        return getAvgNodeTable(avgHandles);
    }

    /**
     * Alias for getAvgNodeTable(SolverAvgHandles avgHandles, boolean keepDisabled).
     *
     * @param avgHandles the average handles to use
     * @param keepDisabled whether to include disabled metrics in the table
     * @return table containing node-level metrics for each class
     */
    public NetworkAvgNodeTable getNodeAvgT(SolverAvgHandles avgHandles, boolean keepDisabled) {
        return getAvgNodeTable(avgHandles, keepDisabled);
    }

    /**
     * Alias for getAvgNodeTable(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A).
     *
     * @param Q queue length handle
     * @param U utilization handle
     * @param R response time handle
     * @param W residence time handle
     * @param T throughput handle
     * @param A arrival rate handle
     * @return table containing node-level metrics for each class
     */
    public NetworkAvgNodeTable getNodeAvgT(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A) {
        return getAvgNodeTable(Q, U, R, W, T, A);
    }

    /**
     * Alias for getAvgNodeTable(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A, boolean keepDisabled).
     *
     * @param Q queue length handle
     * @param U utilization handle
     * @param R response time handle
     * @param W residence time handle
     * @param T throughput handle
     * @param A arrival rate handle
     * @param keepDisabled whether to include disabled metrics in the table
     * @return table containing node-level metrics for each class
     */
    public NetworkAvgNodeTable getNodeAvgT(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A, boolean keepDisabled) {
        return getAvgNodeTable(Q, U, R, W, T, A, keepDisabled);
    }

    /**
     * Alias for getAvgNodeTable().
     *
     * @return table containing node-level metrics for each class
     */
    public NetworkAvgNodeTable nodeAvgT() {
        return getAvgNodeTable();
    }

    /**
     * Alias for getAvgNodeTable(boolean keepDisabled).
     *
     * @param keepDisabled whether to include disabled metrics in the table
     * @return table containing node-level metrics for each class
     */
    public NetworkAvgNodeTable nodeAvgT(boolean keepDisabled) {
        return getAvgNodeTable(keepDisabled);
    }

    /**
     * Alias for getAvgNodeTable(SolverAvgHandles avgHandles).
     *
     * @param avgHandles the average handles to use
     * @return table containing node-level metrics for each class
     */
    public NetworkAvgNodeTable nodeAvgT(SolverAvgHandles avgHandles) {
        return getAvgNodeTable(avgHandles);
    }

    /**
     * Alias for getAvgNodeTable(SolverAvgHandles avgHandles, boolean keepDisabled).
     *
     * @param avgHandles the average handles to use
     * @param keepDisabled whether to include disabled metrics in the table
     * @return table containing node-level metrics for each class
     */
    public NetworkAvgNodeTable nodeAvgT(SolverAvgHandles avgHandles, boolean keepDisabled) {
        return getAvgNodeTable(avgHandles, keepDisabled);
    }

    /**
     * Alias for getAvgNodeTable(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A).
     *
     * @param Q queue length handle
     * @param U utilization handle
     * @param R response time handle
     * @param W residence time handle
     * @param T throughput handle
     * @param A arrival rate handle
     * @return table containing node-level metrics for each class
     */
    public NetworkAvgNodeTable nodeAvgT(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A) {
        return getAvgNodeTable(Q, U, R, W, T, A);
    }

    /**
     * Alias for getAvgNodeTable(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A, boolean keepDisabled).
     *
     * @param Q queue length handle
     * @param U utilization handle
     * @param R response time handle
     * @param W residence time handle
     * @param T throughput handle
     * @param A arrival rate handle
     * @param keepDisabled whether to include disabled metrics in the table
     * @return table containing node-level metrics for each class
     */
    public NetworkAvgNodeTable nodeAvgT(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A, boolean keepDisabled) {
        return getAvgNodeTable(Q, U, R, W, T, A, keepDisabled);
    }

    // aNT aliases for getAvgNodeTable
    /**
     * Alias for getAvgNodeTable().
     *
     * @return table containing node-level metrics for each class
     */
    public NetworkAvgNodeTable aNT() {
        return getAvgNodeTable();
    }

    /**
     * Alias for getAvgNodeTable(boolean keepDisabled).
     *
     * @param keepDisabled whether to include disabled metrics in the table
     * @return table containing node-level metrics for each class
     */
    public NetworkAvgNodeTable aNT(boolean keepDisabled) {
        return getAvgNodeTable(keepDisabled);
    }

    /**
     * Alias for getAvgNodeTable(SolverAvgHandles avgHandles).
     *
     * @param avgHandles the average handles to use
     * @return table containing node-level metrics for each class
     */
    public NetworkAvgNodeTable aNT(SolverAvgHandles avgHandles) {
        return getAvgNodeTable(avgHandles);
    }

    /**
     * Alias for getAvgNodeTable(SolverAvgHandles avgHandles, boolean keepDisabled).
     *
     * @param avgHandles the average handles to use
     * @param keepDisabled whether to include disabled metrics in the table
     * @return table containing node-level metrics for each class
     */
    public NetworkAvgNodeTable aNT(SolverAvgHandles avgHandles, boolean keepDisabled) {
        return getAvgNodeTable(avgHandles, keepDisabled);
    }

    /**
     * Alias for getAvgNodeTable(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A).
     *
     * @param Q queue length handle
     * @param U utilization handle
     * @param R response time handle
     * @param W residence time handle
     * @param T throughput handle
     * @param A arrival rate handle
     * @return table containing node-level metrics for each class
     */
    public NetworkAvgNodeTable aNT(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A) {
        return getAvgNodeTable(Q, U, R, W, T, A);
    }

    /**
     * Alias for getAvgNodeTable(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A, boolean keepDisabled).
     *
     * @param Q queue length handle
     * @param U utilization handle
     * @param R response time handle
     * @param W residence time handle
     * @param T throughput handle
     * @param A arrival rate handle
     * @param keepDisabled whether to include disabled metrics in the table
     * @return table containing node-level metrics for each class
     */
    public NetworkAvgNodeTable aNT(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A, boolean keepDisabled) {
        return getAvgNodeTable(Q, U, R, W, T, A, keepDisabled);
    }

    /**
     * Returns average node throughputs aggregated by job chains.
     *
     * @return matrix of node throughputs [nodes x chains]
     */
    public Matrix getAvgNodeTputChain() {
        int C = sn.nchains;
        Matrix TNclass = getAvgNode().TN;
        Matrix TN = new Matrix(sn.nnodes, C);
        for (int c = 0; c < C; c++) {
            Matrix inchain = sn.inchain.get(c);
            for (int i = 0; i < sn.nnodes; i++) {
                for (int k1 = 0; k1 < inchain.length(); k1++) {
                    int k = (int) inchain.get(k1);
                    TN.set(i, c, TN.get(i, c) + TNclass.get(i, k));
                }
            }
        }
        return TN;
    }

    /**
     * Returns average node utilizations aggregated by job chains.
     *
     * @return matrix of node utilizations [nodes x chains]
     */
    public Matrix getAvgNodeUtilChain() {
        int C = sn.nchains;
        Matrix UNclass = getAvgNode().UN;
        Matrix UN = new Matrix(sn.nnodes, C);
        for (int c = 0; c < C; c++) {
            Matrix inchain = sn.inchain.get(c);
            for (int i = 0; i < sn.nnodes; i++) {
                for (int k1 = 0; k1 < inchain.length(); k1++) {
                    int k = (int) inchain.get(k1);
                    UN.set(i, c, UN.get(i, c) + UNclass.get(i, k));
                }
            }
        }
        return UN;
    }

    /**
     * Computes and returns average queue lengths at steady-state.
     * If results are not available, triggers solver execution.
     *
     * @return matrix of average queue lengths [stations x classes]
     */
    public Matrix getAvgQLen() {
        if (!this.hasResults()) {
            this.getAvg();
        }
        return this.result.QN;
    }

    /**
     * Returns average queue lengths aggregated by job chains.
     *
     * @return matrix of queue lengths [stations x chains]
     */
    public Matrix getAvgQLenChain() {
        int C = sn.nchains;
        Matrix QNclass = getAvgQLen();
        Matrix QN = new Matrix(sn.nstations, C);
        for (int c = 0; c < C; c++) {
            Matrix inchain = sn.inchain.get(c);
            for (int i = 0; i < sn.nstations; i++) {
                for (int k1 = 0; k1 < inchain.length(); k1++) {
                    int k = (int) inchain.get(k1);
                    QN.set(i, c, QN.get(i, c) + QNclass.get(i, k));
                }
            }
        }
        return QN;
    }

    /**
     * Returns the average queue length metric handles.
     *
     * @return handles for queue length metrics
     */
    public AvgHandle getAvgQLenHandles() {
        return this.avgHandles.Q;
    }

    // ========== Kotlin-style Alias Methods ==========
    // Aliases for get* methods following Kotlin naming conventions

    /**
     * Computes and returns average residence times in queue (including service).
     * If results are not available, triggers solver execution.
     *
     * @return matrix of average residence times [stations x classes]
     */
    public Matrix getAvgResidT() {
        if (!this.hasResults()) {
            this.getAvg();
        }
        return this.result.WN;
    }

    /**
     * Returns average residence times aggregated by job chains.
     *
     * @return matrix of residence times [stations x chains]
     */
    public Matrix getAvgResidTChain() {
        int C = sn.nchains;
        Matrix WNclass = getAvgResidT();
        Matrix WN = new Matrix(sn.nstations, C);
        for (int c = 0; c < C; c++) {
            Matrix inchain = sn.inchain.get(c);
            for (int i = 0; i < sn.nstations; i++) {
                double sum = 0.0;
                for (int k1 = 0; k1 < inchain.length(); k1++) {
                    int k = (int) inchain.get(k1);
                    sum += WNclass.get(i, k);
                }
                WN.set(i, c, sum);
            }
        }
        return WN;
    }

    /**
     * Returns the average residence time metric handles.
     *
     * @return handles for residence time metrics
     */
    public AvgHandle getAvgResidTHandles() {
        return this.avgHandles.W;
    }

    /**
     * Computes and returns average response times at steady-state.
     * If results are not available, triggers solver execution.
     *
     * @return matrix of average response times [stations x classes]
     */
    public Matrix getAvgRespT() {
        if (!this.hasResults()) {
            this.getAvg();
        }
        return this.result.RN;
    }

    /**
     * Returns average response times aggregated by job chains.
     *
     * @return matrix of response times [stations x chains]
     */
    public Matrix getAvgRespTChain() {
        int C = sn.nchains;
        Matrix RNclass = getAvgRespT();
        Matrix RN = new Matrix(sn.nstations, C);
        Matrix alpha = snGetDemandsChain(sn).alpha;
        for (int c = 0; c < C; c++) {
            Matrix inchain = sn.inchain.get(c);
            for (int i = 0; i < sn.nstations; i++) {
                for (int k1 = 0; k1 < inchain.length(); k1++) {
                    int k = (int) inchain.get(k1);
                    RN.set(i, c, RN.get(i, c) + RNclass.get(i, k) * alpha.get(i, k));
                }
            }
        }
        return RN;
    }

    /**
     * Returns the average response time metric handles.
     *
     * @return handles for response time metrics
     */
    public AvgHandle getAvgRespTHandles() {
        return this.avgHandles.R;
    }

    /**
     * Computes average system-level metrics at steady state.
     * This includes system response times and throughputs aggregated across all chains.
     */
    public void getAvgSys() {

        this.sn = model.getStruct(true);

        this.getAvg();

        // Fork-join support for system response time computation
        if (this.model.hasFork() || this.model.hasJoin()) {
            if (this.model.hasOpenClasses()) {
                line_error(mfilename(new Object() {
                        }),
                        "System response time computation not yet supported with open classes in the presence of fork/join nodes.");
                this.result.RN.fill(Double.NaN);
            } else {
                // For closed networks with fork-join, ensure response times are properly adjusted
                // Fork-join networks require special handling of synchronization delays
                // The join node response time should include synchronization delay
                if (sn.fj != null && !sn.fj.isEmpty()) {
                    // Apply fork-join transformations to response times
                    // This is a simplified implementation - full fork-join support would require
                    // detailed path analysis and synchronization delay computation
                    for (int i = 0; i < sn.nnodes; i++) {
                        if (sn.nodetype.get(i) == NodeType.Join) {
                            // Join nodes get synchronization delay - placeholder implementation
                            // In full implementation, this would be computed using order statistics
                            // based on parallel path response times
                            for (int r = 0; r < sn.nclasses; r++) {
                                if (this.result.RN != null && !this.result.RN.isEmpty()) {
                                    int stationIndex = (int) sn.nodeToStation.get(i);
                                    double currentRT = this.result.RN.get(stationIndex, r);
                                    if (!Double.isNaN(currentRT) && currentRT > 0) {
                                        // Apply basic synchronization delay factor
                                        this.result.RN.set(stationIndex, r, currentRT * 1.2);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        boolean[] completes = new boolean[sn.nclasses];
        for (int idx = 0; idx < sn.nclasses; idx++) {
            completes[idx] = model.getClasses().get(idx).getCompletes();
        }

        // Optimization: Check if the model has any open classes (infinite jobs)
        // This could be optimised by computing the statistics only for open chains
        boolean hasOpenClasses = false;
        for (int k = 0; k < sn.nclasses; k++) {
            if (isInf(sn.njobs.get(k))) {
                hasOpenClasses = true;
                break;
            }
        }
        // Note: Optimization for open classes not yet implemented - would skip closed chain calculations

        // Compute chain visits
        Matrix alpha = new Matrix(sn.nstations, sn.nclasses);
        Matrix CNclass = new Matrix(1, sn.nclasses);
        if (!this.model.hasJoin() && !this.model.hasFork()) {
            for (int c = 0; c < sn.nchains; c++) {
                Matrix inchain = sn.inchain.get(c);
                for (int i = 0; i < inchain.length(); i++) {
                    int r = (int) inchain.get(i);
                    for (int j = 0; j < sn.nstations; j++) {
                        // Not empty and not a source
                        if (this.result.RN != null && !this.result.RN.isEmpty() &&
                                (!(isInf(this.sn.njobs.get(r)) && j == sn.refstat.get(r)))) {
                            CNclass.set(
                                    0,
                                    r,
                                    CNclass.get(0, r)
                                            + sn.visits.get(c).get((int) sn.stationToStateful.get(j), r)
                                            * this.result.RN.get(j, r)
                                            / sn.visits
                                            .get(c)
                                            .get((int) sn.stationToStateful.get((int) sn.refstat.get(r)), r));
                        }
                    }
                }
            }
        }


        for (int c = 0; c < sn.nchains; c++) {
            Matrix inchain = sn.inchain.get(c);
            Matrix completingClasses = Matrix.extractRows(sn.chains, c, c + 1, null);
            for (int i = 0; i < completingClasses.length(); i++) {
                if (!completes[i]) {
                    completingClasses.set(0, i, Double.NaN);
                }
            }

            for (int i = 0; i < sn.nstations; i++) {
                if (sn.refclass.get(c) >= 0) {
                    // For all classes within the chain (a class belongs to a single chain, the reference
                    // station must be identical for all classes within a chain)
                    List<Double> intersection = Matrix.intersect(sn.refclass.findNonNegative(), inchain);
                    for (double value : intersection) {
                        int k = (int) value;
                        double sumVisits = 0.0;
                        for (int idx = 0; idx < completingClasses.length(); idx++) {
                            if (completingClasses.get(idx) == 1) {
                                sumVisits +=
                                        sn.visits
                                                .get(c)
                                                .get(
                                                        (int) sn.stationToStateful.get((int) sn.refstat.get(k)),
                                                        idx);
                            }
                        }
                        alpha.set(
                                i,
                                k,
                                alpha.get(i, k)
                                        + sn.visits.get(c).get((int) sn.stationToStateful.get(i), k) / sumVisits);
                    }
                } else {
                    // For all classes within the chain (a class belongs to a single chain, the reference
                    // station must be identical for all classes within a chain)
                    for (int j = 0; j < inchain.length(); j++) {
                        int k = (int) inchain.get(j);
                        double sumVisits = 0.0;
                        for (int idx = 0; idx < completingClasses.length(); idx++) {
                            if (completingClasses.get(idx) == 1) {
                                sumVisits += sn.visits.get(c).get((int) sn.stationToStateful.get((int) sn.refstat.get(k)), idx);
                            }
                        }
                        alpha.set(i, k, alpha.get(i, k) + sn.visits.get(c).get((int) sn.stationToStateful.get(i), k) / sumVisits);
                    }
                }
            }
        }
        for (int i = 0; i < sn.nstations; i++) {
            for (int k = 0; k < sn.nclasses; k++) {
                if (isInf(alpha.get(i, k)) || Double.isNaN(alpha.get(i, k))) {
                    alpha.set(i, k, 0.0);
                }
            }
        }

        // Compute average chain metrics
        this.result.CN = new Matrix(1, sn.nchains);
        this.result.XN = new Matrix(1, sn.nchains);

        for (int c = 0; c < sn.nchains; c++) {
            Matrix inchain = sn.inchain.get(c);
            Matrix completingClasses = Matrix.extractRows(sn.chains, c, c + 1, null).find();
            completingClasses = completingClasses.transpose();
            for (int i = 0; i < inchain.length(); i++) {
                int classIndex = (int) inchain.get(i);
                if (!this.model.getClasses().get(classIndex).getCompletes()) {
                    completingClasses.set(0, i, Double.NaN);
                }
            }

            if (!result.TN.isEmpty()) {
                // All classes in same chain must share the same refstation, so we use the first one
                int ref = (int) sn.refstat.get((int) inchain.get(0));
                // We now compute the incoming system throughput to the reference station from completing
                // classes
                for (int i = 0; i < sn.nstations; i++) {
                    for (int j = 0; j < completingClasses.length(); j++) {
                        int r = (int) completingClasses.get(j);
                        if (completingClasses.get(j) >= 0) {
                            List<Double> intersection = Matrix.intersect(sn.refclass.findNonNegative(), inchain);
                            for (double value : intersection) {
                                int s = (int) value;
                                if (!Double.isNaN(this.result.TN.get(i, r))) {
                                    this.result.XN.set(
                                            0,
                                            c,
                                            this.result.XN.get(0, c)
                                                    + sn.rt.get(i * sn.nclasses + r, ref * sn.nclasses + s)
                                                    * this.result.TN.get(i, r));
                                }
                            }
                            for (int k = 0; k < inchain.length(); k++) {
                                int s = (int) inchain.get(k);
                                if (!Double.isNaN(this.result.TN.get(i, r))) {
                                    this.result.XN.set(
                                            0,
                                            c,
                                            this.result.XN.get(0, c)
                                                    + sn.rt.get(i * sn.nclasses + r, ref * sn.nclasses + s)
                                                    * this.result.TN.get(i, r));
                                }
                            }
                        }
                    }
                }
            }

            // If this is a closed chain we simply apply Little's law
            int nJobsChain = 0;
            for (int i = 0; i < sn.chains.getNumCols(); i++) {
                if (sn.chains.get(c, i) > 0) {
                    nJobsChain += sn.njobs.get(i);
                }
            }

            if (this.model.hasFork() && this.model.hasJoin()) {
                // In this case, CN is unreliable as it sums the contribution across all stations,
                // which would include also forked tasks, we use Little's law instead
                this.result.CN.set(0, c, nJobsChain / this.result.XN.get(0, c));
            } else {
                // Standard chain response time computation for non-fork-join networks
                if (isInf(nJobsChain)) {
                    if (inchain.length() != completingClasses.length()) {
                        throw new RuntimeException(
                                "Edge-based chain definition not yet supported for open queueing networks.");
                    }
                }
                double sumFinite = 0;
                for (int i = 0; i < inchain.length(); i++) {
                    double value = alpha.get((int) sn.refstat.get((int) inchain.get(0)), (int) inchain.get(i))
                            * CNclass.get((int) inchain.get(i));
                    if (!isInf(value) && !Double.isNaN(value))
                        sumFinite += value;
                }

                this.result.CN.set(0, c, sumFinite);
            }
        }
    }

    /**
     * Computes average system-level metrics at steady state using specified handles.
     * This includes system response times and throughputs aggregated across all chains.
     *
     * @param avgHandles custom handles for performance metrics
     */
    public void getAvgSys(SolverAvgHandles avgHandles) {
        SolverAvgHandles originalHandles = this.avgHandles;
        this.avgHandles = avgHandles;
        try {
            getAvgSys();
        } finally {
            this.avgHandles = originalHandles;
        }
    }

    /**
     * Computes average system-level metrics at steady state using individual handles.
     * This includes system response times and throughputs aggregated across all chains.
     *
     * @param Q queue length handle
     * @param U utilization handle
     * @param R response time handle
     * @param W residence time handle
     * @param T throughput handle
     * @param A arrival rate handle
     */
    public void getAvgSys(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A) {
        SolverAvgHandles customHandles = new SolverAvgHandles(Q, U, R, W, T, A);
        getAvgSys(customHandles);
    }

    /**
     * Returns average system response times at steady state.
     *
     * @return matrix of system response times by chain
     */
    public Matrix getAvgSysRespT() {
        if (!this.hasResults()) {
            this.getAvgSys();
        }
        return this.result.CN;
    }

    /**
     * Returns a table of average system-level metrics.
     *
     * @return table containing system response times and throughputs by chain
     */
    public NetworkAvgSysTable getAvgSysTable() {

        this.getAvgSys();

        NetworkAvgSysTable avgSysTable = new NetworkAvgSysTable(this.result.CN.toList1D(), this.result.XN.toList1D(), this.options);


        java.util.List<String> chainNames = new ArrayList<>();
        java.util.List<String> inChainNames = new ArrayList<>();
        for (int c = 0; c < this.sn.nchains; c++) {
            chainNames.add("Chain" + (c + 1));
            Matrix inchain = sn.inchain.get(c);
            String chainMembers = "(";
            for (int i = 0; i < inchain.length(); i++) {
                int r = (int) inchain.get(i);
                if (i == 0) {
                    chainMembers = chainMembers.concat(sn.classnames.get(r));
                } else {
                    chainMembers = chainMembers.concat(" " + sn.classnames.get(r));
                }
            }
            inChainNames.add(chainMembers + ")");
        }

        avgSysTable.setChainNames(chainNames);
        avgSysTable.setInChainNames(inChainNames);

        return avgSysTable;
    }

    /**
     * Returns a table of average system-level metrics using specified handles.
     *
     * @param avgHandles custom handles for performance metrics
     * @return table containing system response times and throughputs by chain
     */
    public NetworkAvgSysTable getAvgSysTable(SolverAvgHandles avgHandles) {
        SolverAvgHandles originalHandles = this.avgHandles;
        this.avgHandles = avgHandles;
        try {
            return getAvgSysTable();
        } finally {
            this.avgHandles = originalHandles;
        }
    }

    /**
     * Returns a table of average system-level metrics using individual handles.
     *
     * @param Q queue length handle
     * @param U utilization handle
     * @param R response time handle
     * @param W residence time handle
     * @param T throughput handle
     * @param A arrival rate handle
     * @return table containing system response times and throughputs by chain
     */
    public NetworkAvgSysTable getAvgSysTable(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A) {
        SolverAvgHandles customHandles = new SolverAvgHandles(Q, U, R, W, T, A);
        return getAvgSysTable(customHandles);
    }

    /**
     * Alias for getAvgSysTable().
     *
     * @return table containing system response times and throughputs by chain
     */
    public NetworkAvgSysTable getSysAvgT() {
        return getAvgSysTable();
    }

    /**
     * Alias for getAvgSysTable(SolverAvgHandles avgHandles).
     *
     * @param avgHandles the average handles to use
     * @return table containing system response times and throughputs by chain
     */
    public NetworkAvgSysTable getSysAvgT(SolverAvgHandles avgHandles) {
        return getAvgSysTable(avgHandles);
    }

    /**
     * Alias for getAvgSysTable(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A).
     *
     * @param Q queue length handle
     * @param U utilization handle
     * @param R response time handle
     * @param W residence time handle
     * @param T throughput handle
     * @param A arrival rate handle
     * @return table containing system response times and throughputs by chain
     */
    public NetworkAvgSysTable getSysAvgT(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A) {
        return getAvgSysTable(Q, U, R, W, T, A);
    }

    /**
     * Alias for getAvgSysTable().
     *
     * @return table containing system response times and throughputs by chain
     */
    public NetworkAvgSysTable sysAvgT() {
        return getAvgSysTable();
    }

    /**
     * Alias for getAvgSysTable(SolverAvgHandles avgHandles).
     *
     * @param avgHandles the average handles to use
     * @return table containing system response times and throughputs by chain
     */
    public NetworkAvgSysTable sysAvgT(SolverAvgHandles avgHandles) {
        return getAvgSysTable(avgHandles);
    }

    /**
     * Alias for getAvgSysTable(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A).
     *
     * @param Q queue length handle
     * @param U utilization handle
     * @param R response time handle
     * @param W residence time handle
     * @param T throughput handle
     * @param A arrival rate handle
     * @return table containing system response times and throughputs by chain
     */
    public NetworkAvgSysTable sysAvgT(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A) {
        return getAvgSysTable(Q, U, R, W, T, A);
    }

    /**
     * Short alias for getAvgSysTable().
     *
     * @return table containing system response times and throughputs by chain
     */
    public NetworkAvgSysTable aST() {
        return getAvgSysTable();
    }

    /**
     * Short alias for getAvgSysTable(SolverAvgHandles avgHandles).
     *
     * @param avgHandles the average handles to use
     * @return table containing system response times and throughputs by chain
     */
    public NetworkAvgSysTable aST(SolverAvgHandles avgHandles) {
        return getAvgSysTable(avgHandles);
    }

    /**
     * Short alias for getAvgSysTable(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A).
     *
     * @param Q queue length handle
     * @param U utilization handle
     * @param R response time handle
     * @param W residence time handle
     * @param T throughput handle
     * @param A arrival rate handle
     * @return table containing system response times and throughputs by chain
     */
    public NetworkAvgSysTable aST(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A) {
        return getAvgSysTable(Q, U, R, W, T, A);
    }

    /**
     * Returns average system throughputs at steady state.
     *
     * @return matrix of system throughputs by chain
     */
    public Matrix getAvgSysTput() {
        if (!this.hasResults()) {
            this.getAvgSys();
        }
        return this.result.XN;
    }

    /**
     * Returns a table of average station metrics organized by job classes.
     *
     * @return table containing station-level metrics for each class
     */
    public NetworkAvgTable getAvgTable() {

        this.sn = model.getStruct(true);

        boolean keepDisabled = false;
        this.avgHandles = model.getAvgHandles();

        int M = sn.nstations;
        int K = sn.nclasses;

        try {
            if (Double.isFinite(options.timespan[1])) {
                this.tranHandles = model.getTranHandles();
                getTranAvg();
            } else {
                getAvg();
            }
        } catch (Exception e) {
            e.printStackTrace();
            throw new RuntimeException(
                    "Unable to compute results and therefore unable to print AvgTable.", e);
        }

        // Check if result is null or incomplete before accessing fields
        if (this.result == null) {
            throw new RuntimeException(
                    "Unable to compute results and therefore unable to print AvgTable.");
        }

        Matrix QN = this.result.QN;
        Matrix UN = this.result.UN;
        Matrix RN = this.result.RN;
        Matrix WN = this.result.WN;
        Matrix TN = this.result.TN;
        Matrix AN = this.result.AN;

        if (QN == null || QN.isEmpty()) {
            throw new RuntimeException(
                    "Unable to compute results and therefore unable to print AvgTable.");
        }

        if (!keepDisabled) {
            Matrix V = new Matrix(sn.nstateful, K);
            for (int i = 0; i < sn.visits.size(); i++) {
                V = V.add(1, sn.visits.get(i));
            }
            if (V.isEmpty()) { // SSA
                // Implementation for SSA: sum all chain visits across all chains
                for (Integer chainIndex : sn.visits.keySet()) {
                    Matrix chainVisits = sn.visits.get(chainIndex);
                    if (chainVisits != null) {
                        V = V.add(1.0, chainVisits);
                    }
                }
            }

            List<Double> Qval = new ArrayList<>();
            List<Double> Uval = new ArrayList<>();
            List<Double> Rval = new ArrayList<>();
            List<Double> Tval = new ArrayList<>();
            List<Double> ArvR = new ArrayList<>();
            List<Double> Residval = new ArrayList<>();
            List<String> className = new ArrayList<>();
            List<String> stationName = new ArrayList<>();
            for (int i = 0; i < M; i++) {
                for (int k = 0; k < K; k++) {
                    if (QN.get(i, k) + UN.get(i, k) + RN.get(i, k) + TN.get(i, k) + AN.get(i, k) + WN.get(i, k) > 0) {
                        int c = -1;
                        for (int row = 0; row < sn.chains.getNumRows(); row++) {
                            if (sn.chains.get(row, k) > 0) {
                                c = row;
                                break;
                            }
                        }
                        Qval.add(QN.get(i, k));
                        Uval.add(UN.get(i, k));
                        Rval.add(RN.get(i, k));
                        ArvR.add(AN.get(i, k));
                        Tval.add(TN.get(i, k));
                        className.add(model.getClasses().get(k).getName());
                        stationName.add(this.model.getStations().get(i).getName());
                        Residval.add(WN.get(i, k));
                    }
                }
            }
            NetworkAvgTable avgTable = new NetworkAvgTable(Qval, Uval, Rval, Residval, ArvR, Tval);
            avgTable.setOptions(this.options);
            avgTable.setClassNames(className);
            avgTable.setStationNames(stationName);
            return avgTable;
        } else {
            // Keep all entries including disabled ones (keepDisabled == true)
            Matrix V = new Matrix(sn.nstateful, K);
            for (int i = 0; i < sn.visits.size(); i++) {
                V = V.add(1, sn.visits.get(i));
            }
            if (V.isEmpty()) { // SSA
                // Implementation for SSA: sum all chain visits across all chains
                for (Integer chainIndex : sn.visits.keySet()) {
                    Matrix chainVisits = sn.visits.get(chainIndex);
                    if (chainVisits != null) {
                        V = V.add(1.0, chainVisits);
                    }
                }
            }

            List<Double> Qval = new ArrayList<>();
            List<Double> Uval = new ArrayList<>();
            List<Double> Rval = new ArrayList<>();
            List<Double> Tval = new ArrayList<>();
            List<Double> ArvR = new ArrayList<>();
            List<Double> Residval = new ArrayList<>();
            List<String> className = new ArrayList<>();
            List<String> stationName = new ArrayList<>();
            for (int i = 0; i < M; i++) {
                for (int k = 0; k < K; k++) {
                    // Include all entries regardless of their values
                    Qval.add(QN.get(i, k));
                    Uval.add(UN.get(i, k));
                    Rval.add(RN.get(i, k));
                    ArvR.add(AN.get(i, k));
                    Tval.add(TN.get(i, k));
                    className.add(model.getClasses().get(k).getName());
                    stationName.add(this.model.getStations().get(i).getName());
                    Residval.add(WN.get(i, k));
                }
            }
            NetworkAvgTable avgTable = new NetworkAvgTable(Qval, Uval, Rval, Residval, ArvR, Tval);
            avgTable.setOptions(this.options);
            avgTable.setClassNames(className);
            avgTable.setStationNames(stationName);
            return avgTable;
        }
    }

    /**
     * Returns a table of average station metrics organized by job classes using specified handles.
     *
     * @param avgHandles custom handles for performance metrics
     * @return table containing station-level metrics for each class
     */
    public NetworkAvgTable getAvgTable(SolverAvgHandles avgHandles) {
        SolverAvgHandles originalHandles = this.avgHandles;
        this.avgHandles = avgHandles;
        try {
            return getAvgTable();
        } finally {
            this.avgHandles = originalHandles;
        }
    }

    /**
     * Returns a table of average station metrics organized by job classes using individual handles.
     *
     * @param Q queue length handle
     * @param U utilization handle
     * @param R response time handle
     * @param W residence time handle
     * @param T throughput handle
     * @param A arrival rate handle
     * @return table containing station-level metrics for each class
     */
    public NetworkAvgTable getAvgTable(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A) {
        SolverAvgHandles customHandles = new SolverAvgHandles(Q, U, R, W, T, A);
        return getAvgTable(customHandles);
    }

    /**
     * Returns a table of average station metrics organized by job classes with keepDisabled option.
     *
     * @param keepDisabled whether to include disabled metrics in the table
     * @return table containing station-level metrics for each class
     */
    public NetworkAvgTable getAvgTable(boolean keepDisabled) {
        // Implement keepDisabled functionality
        this.sn = model.getStruct(true);
        this.avgHandles = model.getAvgHandles();

        int M = sn.nstations;
        int K = sn.nclasses;

        try {
            if (!isInf(options.timespan[1]) && !Double.isNaN(options.timespan[1])) {
                this.tranHandles = model.getTranHandles();
                //getTranAvg();
            } else {
                getAvg();
            }
        } catch (Exception e) {
            e.printStackTrace();
        }

        Matrix QN = this.result.QN;
        Matrix UN = this.result.UN;
        Matrix RN = this.result.RN;
        Matrix TN = this.result.TN;
        Matrix AN = this.result.AN;
        Matrix WN = this.result.WN;

        if (QN == null || QN.isEmpty()) {
            throw new RuntimeException(
                    "Unable to compute results and therefore unable to print AvgTable.");
        }

        Matrix V = new Matrix(sn.nstateful, K);
        for (int i = 0; i < sn.visits.size(); i++) {
            V = V.add(1, sn.visits.get(i));
        }
        if (V.isEmpty()) { // SSA
            // Implementation for SSA: sum all chain visits across all chains
            for (Integer chainIndex : sn.visits.keySet()) {
                Matrix chainVisits = sn.visits.get(chainIndex);
                if (chainVisits != null) {
                    V = V.add(1.0, chainVisits);
                }
            }
        }

        List<Double> Qval = new ArrayList<>();
        List<Double> Uval = new ArrayList<>();
        List<Double> Rval = new ArrayList<>();
        List<Double> Tval = new ArrayList<>();
        List<Double> ArvR = new ArrayList<>();
        List<Double> Residval = new ArrayList<>();
        List<String> className = new ArrayList<>();
        List<String> stationName = new ArrayList<>();
        for (int i = 0; i < M; i++) {
            for (int k = 0; k < K; k++) {
                // Include entry based on keepDisabled flag
                if (keepDisabled || (QN.get(i, k) + UN.get(i, k) + RN.get(i, k) + TN.get(i, k) + AN.get(i, k) + WN.get(i, k) > 0)) {
                    Qval.add(QN.get(i, k));
                    Uval.add(UN.get(i, k));
                    Rval.add(RN.get(i, k));
                    ArvR.add(AN.get(i, k));
                    Tval.add(TN.get(i, k));
                    className.add(model.getClasses().get(k).getName());
                    stationName.add(this.model.getStations().get(i).getName());
                    Residval.add(WN.get(i, k));
                }
            }
        }
        NetworkAvgTable avgTable = new NetworkAvgTable(Qval, Uval, Rval, Residval, ArvR, Tval);
        avgTable.setOptions(this.options);
        avgTable.setClassNames(className);
        avgTable.setStationNames(stationName);
        return avgTable;
    }

    /**
     * Returns a table of average station metrics organized by job classes using specified handles and keepDisabled option.
     *
     * @param avgHandles custom handles for performance metrics
     * @param keepDisabled whether to include disabled metrics in the table
     * @return table containing station-level metrics for each class
     */
    public NetworkAvgTable getAvgTable(SolverAvgHandles avgHandles, boolean keepDisabled) {
        SolverAvgHandles originalHandles = this.avgHandles;
        this.avgHandles = avgHandles;
        try {
            return getAvgTable(keepDisabled);
        } finally {
            this.avgHandles = originalHandles;
        }
    }

    /**
     * Returns a table of average station metrics organized by job classes using individual handles and keepDisabled option.
     *
     * @param Q queue length handle
     * @param U utilization handle
     * @param R response time handle
     * @param W residence time handle
     * @param T throughput handle
     * @param A arrival rate handle
     * @param keepDisabled whether to include disabled metrics in the table
     * @return table containing station-level metrics for each class
     */
    public NetworkAvgTable getAvgTable(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A, boolean keepDisabled) {
        SolverAvgHandles customHandles = new SolverAvgHandles(Q, U, R, W, T, A);
        return getAvgTable(customHandles, keepDisabled);
    }

    /**
     * Alias for getAvgTable().
     *
     * @return table containing station-level metrics for each class
     */
    public NetworkAvgTable getAvgT() {
        return getAvgTable();
    }

    /**
     * Alias for getAvgTable(boolean keepDisabled).
     *
     * @param keepDisabled whether to include disabled metrics in the table
     * @return table containing station-level metrics for each class
     */
    public NetworkAvgTable getAvgT(boolean keepDisabled) {
        return getAvgTable(keepDisabled);
    }

    /**
     * Alias for getAvgTable(SolverAvgHandles avgHandles).
     *
     * @param avgHandles the average handles to use
     * @return table containing station-level metrics for each class
     */
    public NetworkAvgTable getAvgT(SolverAvgHandles avgHandles) {
        return getAvgTable(avgHandles);
    }

    /**
     * Alias for getAvgTable(SolverAvgHandles avgHandles, boolean keepDisabled).
     *
     * @param avgHandles the average handles to use
     * @param keepDisabled whether to include disabled metrics in the table
     * @return table containing station-level metrics for each class
     */
    public NetworkAvgTable getAvgT(SolverAvgHandles avgHandles, boolean keepDisabled) {
        return getAvgTable(avgHandles, keepDisabled);
    }

    /**
     * Alias for getAvgTable(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A).
     *
     * @param Q queue length handle
     * @param U utilization handle
     * @param R response time handle
     * @param W residence time handle
     * @param T throughput handle
     * @param A arrival rate handle
     * @return table containing station-level metrics for each class
     */
    public NetworkAvgTable getAvgT(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A) {
        return getAvgTable(Q, U, R, W, T, A);
    }

    /**
     * Alias for getAvgTable(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A, boolean keepDisabled).
     *
     * @param Q queue length handle
     * @param U utilization handle
     * @param R response time handle
     * @param W residence time handle
     * @param T throughput handle
     * @param A arrival rate handle
     * @param keepDisabled whether to include disabled metrics in the table
     * @return table containing station-level metrics for each class
     */
    public NetworkAvgTable getAvgT(AvgHandle Q, AvgHandle U, AvgHandle R, AvgHandle W, AvgHandle T, AvgHandle A, boolean keepDisabled) {
        return getAvgTable(Q, U, R, W, T, A, keepDisabled);
    }

    /**
     * Returns a table of deadline-related metrics (response time and tardiness) organized by station and job class.
     *
     * @return table containing response time and tardiness metrics, or null if tardiness data is not available
     */
    public NetworkAvgTable getDeadlineTable() {
        this.sn = model.getStruct(true);

        int M = sn.nstations;
        int K = sn.nclasses;

        try {
            if (Double.isFinite(options.timespan[1])) {
                this.tranHandles = model.getTranHandles();
                getTranAvg();
            } else {
                getAvg();
            }
        } catch (Exception e) {
            e.printStackTrace();
            throw new RuntimeException(
                    "Unable to compute results and therefore unable to get deadline table.");
        }

        if (this.result == null || this.result.TardN == null || this.result.SysTardN == null) {
            return null;
        }

        Matrix RN = this.result.RN;
        Matrix TardN = this.result.TardN;
        Matrix SysTardN = this.result.SysTardN;

        if (TardN.isEmpty() || SysTardN.isEmpty()) {
            return null;
        }

        List<Double> Rval = new ArrayList<>();
        List<Double> Tardval = new ArrayList<>();
        List<Double> SysTardval = new ArrayList<>();
        List<String> className = new ArrayList<>();
        List<String> stationName = new ArrayList<>();

        for (int i = 0; i < M; i++) {
            for (int k = 0; k < K; k++) {
                if (RN.get(i, k) > 0 || TardN.get(i, k) > 0 || SysTardN.get(0, k) > 0) {
                    Rval.add(RN.get(i, k));
                    Tardval.add(TardN.get(i, k));
                    SysTardval.add(SysTardN.get(0, k));
                    className.add(model.getClasses().get(k).getName());
                    stationName.add(this.model.getStations().get(i).getName());
                }
            }
        }

        NetworkAvgTable deadlineTable = new NetworkAvgTable(
                new ArrayList<>(), new ArrayList<>(), Rval,
                new ArrayList<>(), new ArrayList<>(), new ArrayList<>(),
                Tardval, SysTardval);
        deadlineTable.setOptions(this.options);
        deadlineTable.setClassNames(className);
        deadlineTable.setStationNames(stationName);
        return deadlineTable;
    }

    /**
     * Computes and returns average throughputs at steady-state.
     * If results are not available, triggers solver execution.
     *
     * @return matrix of average throughputs [stations x classes]
     */
    public Matrix getAvgTput() {
        if (!this.hasResults()) {
            this.getAvg();
        }
        return this.result.TN;
    }

    /**
     * Returns average throughputs aggregated by job chains.
     *
     * @return matrix of throughputs [stations x chains]
     */
    public Matrix getAvgTputChain() {
        int C = sn.nchains;
        Matrix TNclass = getAvgTput();
        Matrix TN = new Matrix(sn.nstations, C);
        for (int c = 0; c < C; c++) {
            Matrix inchain = sn.inchain.get(c);
            for (int i = 0; i < sn.nstations; i++) {
                for (int k1 = 0; k1 < inchain.length(); k1++) {
                    int k = (int) inchain.get(k1);
                    TN.set(i, c, TN.get(i, c) + TNclass.get(i, k));
                }
            }
        }
        return TN;
    }

    /**
     * Returns the average throughput metric handles.
     *
     * @return handles for throughput metrics
     */
    public AvgHandle getAvgTputHandles() {
        return this.avgHandles.T;
    }

    /**
     * Computes and returns average server utilizations at steady-state.
     * If results are not available, triggers solver execution.
     *
     * @return matrix of average utilizations [stations x classes]
     */
    public Matrix getAvgUtil() {
        if (!this.hasResults()) {
            this.getAvg();
        }
        return this.result.UN;
    }

    /**
     * Returns average server utilizations aggregated by job chains.
     *
     * @return matrix of utilizations [stations x chains]
     */
    public Matrix getAvgUtilChain() {
        int C = sn.nchains;
        Matrix UNclass = getAvgUtil();
        Matrix UN = new Matrix(sn.nstations, C);
        for (int c = 0; c < C; c++) {
            Matrix inchain = sn.inchain.get(c);
            for (int i = 0; i < sn.nstations; i++) {
                for (int k1 = 0; k1 < inchain.length(); k1++) {
                    int k = (int) inchain.get(k1);
                    UN.set(i, c, UN.get(i, c) + UNclass.get(i, k));
                }
            }
        }
        return UN;
    }

    /**
     * Returns the average utilization metric handles.
     *
     * @return handles for utilization metrics
     */
    public AvgHandle getAvgUtilHandles() {
        return this.avgHandles.U;
    }

    /**
     * Computes and returns average waiting times in queue excluding service time.
     * Waiting time = Response time - Service time (1/rate)
     *
     * @return matrix of average waiting times [stations x classes]
     */
    public Matrix getAvgWaitT() {
        if (!this.hasResults()) {
            this.getAvg();
        }

        Matrix RN = this.result.RN;
        if (RN == null || RN.isEmpty()) {
            return new Matrix(0, 0);
        }

        // Get service rates for waiting time calculation
        this.sn = model.getStruct(true);
        Matrix WT = RN.copy();

        // Calculate waiting time = response time - service time (1/rate)
        for (int i = 0; i < sn.nstations; i++) {
            for (int k = 0; k < sn.nclasses; k++) {
                double serviceTime = 1.0 / sn.rates.get(i, k);
                WT.set(i, k, RN.get(i, k) - serviceTime);

                // Set waiting time to 0 for source nodes (they don't have queues)
                if (sn.nodetype.get((int) sn.stationToNode.get(i)) == NodeType.Source) {
                    WT.set(i, k, 0.0);
                }

                // Ensure non-negative waiting times (numerical safety)
                if (WT.get(i, k) < 0) {
                    WT.set(i, k, 0.0);
                }
            }
        }

        return WT;
    }

    /**
     * Returns cumulative distribution functions of passage times at steady-state.
     * This is an abstract method that must be implemented by concrete solver subclasses.
     *
     * @param R response time handles (optional)
     * @return result containing CDFs for passage times
     * @throws RuntimeException if not implemented by the concrete solver
     */
    public DistributionResult getCdfPassT(AvgHandle R) {
        throw new RuntimeException("getCdfPassT is not supported by " + this.getClass().getSimpleName());
    }

    /**
     * Returns cumulative distribution functions of passage times at steady-state.
     * Uses default response time handles.
     *
     * @return result containing CDFs for passage times
     * @throws RuntimeException if not implemented by the concrete solver
     */
    public DistributionResult getCdfPassT() {
        return getCdfPassT(getAvgRespTHandles());
    }

    /**
     * Returns cumulative distribution functions of response times at steady-state.
     * This is an abstract method that must be implemented by concrete solver subclasses.
     *
     * @param R response time handles (optional)
     * @return result containing CDFs for response times [stations x classes]
     * @throws RuntimeException if not implemented by the concrete solver
     */
    public DistributionResult getCdfRespT(AvgHandle R) {
        throw new RuntimeException("getCdfRespT is not supported by " + this.getClass().getSimpleName());
    }

    /**
     * Returns cumulative distribution functions of response times at steady-state.
     * Uses default response time handles.
     *
     * @return result containing CDFs for response times [stations x classes]
     * @throws RuntimeException if not implemented by the concrete solver
     */
    public DistributionResult getCdfRespT() {
        return getCdfRespT(getAvgRespTHandles());
    }

    /**
     * Returns the queueing network model being solved.
     *
     * @return the network model
     */
    public Network getModel() {
        return model;
    }

    /**
     * Sets the queueing network model to be solved.
     *
     * @param model the network model to set
     */
    public void setModel(Network model) {
        this.model = model;
    }

    /**
     * Returns marginal state probabilities for a specific node and state.
     * This is an abstract method that must be implemented by concrete solver subclasses.
     *
     * @param node  the node index for which to compute probabilities
     * @param state the state vector to query (optional, null for all states)
     * @return result containing marginal state probabilities
     * @throws RuntimeException if not implemented by the concrete solver
     */
    public ProbabilityResult getProb(int node, Matrix state) {
        throw new RuntimeException("getProb is not supported by " + this.getClass().getSimpleName());
    }

    /**
     * Returns marginal state probabilities for a specific node (all states).
     * This is an abstract method that must be implemented by concrete solver subclasses.
     *
     * @param node the node index for which to compute probabilities
     * @return result containing marginal state probabilities
     * @throws RuntimeException if not implemented by the concrete solver
     */
    public ProbabilityResult getProb(int node) {
        return getProb(node, null);
    }

    /**
     * Probability of a SPECIFIC per-class job distribution at a station.
     * Returns P(n1 jobs of class 1, n2 jobs of class 2, ...) for given state.
     *
     * <p>Compare with {@link #getProbMarg}: returns queue-length distribution for a
     * single class, i.e., P(n jobs of class r) for n=0,1,...,N(r).</p>
     *
     * @param node    the node index for which to compute probabilities
     * @param state_a per-class job counts, e.g., [2,1] = 2 class-1, 1 class-2
     * @return scalar probability in [0,1]
     * @throws RuntimeException if not implemented by the concrete solver
     */
    public ProbabilityResult getProbAggr(int node, Matrix state_a) {
        throw new RuntimeException("getProbAggr is not supported by " + this.getClass().getSimpleName());
    }

    /**
     * Probability of a SPECIFIC per-class job distribution at a station (current state).
     * Returns P(n1 jobs of class 1, n2 jobs of class 2, ...).
     *
     * <p>Compare with {@link #getProbMarg}: returns queue-length distribution for a
     * single class, i.e., P(n jobs of class r) for n=0,1,...,N(r).</p>
     *
     * @param node the node index for which to compute probabilities
     * @return scalar probability in [0,1]
     * @throws RuntimeException if not implemented by the concrete solver
     */
    public ProbabilityResult getProbAggr(int node) {
        return getProbAggr(node, null);
    }

    /**
     * Probability distribution for queue length of a SINGLE class at a station.
     * Returns P(n jobs of class r) for n=0,1,...,N(r).
     *
     * <p>Compare with {@link #getProbAggr}: returns probability of a specific per-class
     * distribution, e.g., P(2 class-1, 1 class-2) as a scalar.</p>
     *
     * @param node     the node index for which to compute probabilities
     * @param jobclass the job class index for marginalization
     * @param state_m  specific states to query, or null for all
     * @return vector where element n+1 = P(n jobs of this class)
     * @throws RuntimeException if not implemented by the concrete solver
     */
    public ProbabilityResult getProbMarg(int node, int jobclass, Matrix state_m) {
        throw new RuntimeException("getProbMarg is not supported by " + this.getClass().getSimpleName());
    }

    /**
     * Probability distribution for queue length of a SINGLE class at a station (all states).
     * Returns P(n jobs of class r) for n=0,1,...,N(r).
     *
     * <p>Compare with {@link #getProbAggr}: returns probability of a specific per-class
     * distribution, e.g., P(2 class-1, 1 class-2) as a scalar.</p>
     *
     * @param node     the node index for which to compute probabilities
     * @param jobclass the job class index for marginalization
     * @return vector where element n+1 = P(n jobs of this class)
     * @throws RuntimeException if not implemented by the concrete solver
     */
    public ProbabilityResult getProbMarg(int node, int jobclass) {
        return getProbMarg(node, jobclass, null);
    }

    /**
     * Returns the logarithm of the normalizing constant of state probabilities.
     * This is an abstract method that must be implemented by concrete solver subclasses.
     *
     * @return result containing the log normalizing constant
     * @throws RuntimeException if not implemented by the concrete solver
     */
    public ProbabilityResult getProbNormConstAggr() {
        throw new RuntimeException("getProbNormConstAggr is not supported by " + this.getClass().getSimpleName());
    }

    /**
     * Returns joint state probabilities for the entire system.
     * This is an abstract method that must be implemented by concrete solver subclasses.
     *
     * @return result containing joint state probabilities
     * @throws RuntimeException if not implemented by the concrete solver
     */
    public ProbabilityResult getProbSys() {
        throw new RuntimeException("getProbSys is not supported by " + this.getClass().getSimpleName());
    }

    /**
     * Returns aggregated joint state probabilities for the entire system.
     * This is an abstract method that must be implemented by concrete solver subclasses.
     *
     * @return result containing aggregated joint state probabilities
     * @throws RuntimeException if not implemented by the concrete solver
     */
    public ProbabilityResult getProbSysAggr() {
        throw new RuntimeException("getProbSysAggr is not supported by " + this.getClass().getSimpleName());
    }

    /**
     * Returns a table of average stage metrics organized by job classes.
     * Note: This method is not yet fully implemented and will throw an exception.
     *
     * @return table containing stage-level metrics for each class
     * @throws UnsupportedOperationException indicating this method is not yet implemented
     */
    public NetworkAvgTable getStageTable() {
        throw new UnsupportedOperationException("getStageTable is not yet implemented in NetworkSolver");
    }

    /**
     * Returns a table of average stage metrics organized by job classes with keepDisabled option.
     * Note: This method is not yet fully implemented and will throw an exception.
     *
     * @param keepDisabled whether to include disabled metrics in the table
     * @return table containing stage-level metrics for each class
     * @throws UnsupportedOperationException indicating this method is not yet implemented
     */
    public NetworkAvgTable getStageTable(boolean keepDisabled) {
        throw new UnsupportedOperationException("getStageTable is not yet implemented in NetworkSolver");
    }

    /**
     * Computes transient average station metrics over the specified time interval.
     * The timespan is defined in the solver options.
     */
    public void getTranAvg() {

        this.tranHandles = model.getTranHandles();

        // NOTE: This was in LINE, but I believe is legacy as 'matrix' method can provide tran results
    /*    if (!Objects.equals(options.method, "default")) {
      System.err.println(
          "getTranAvg is not offered by the specified method. Setting the solution method to \"closing\".");
      resetResults();
    }
    options.method = "closing";*/

        sn = model.getStruct(true);
        double minRate = sn.rates.elementMin();
        if (!hasTranResults()) {
            if (isInf(options.timespan[0]) && isInf(options.timespan[1])) {
                options.timespan[0] = 0;
                options.timespan[1] = 30 / minRate;
                System.out.format(
                        "Timespan of transient analysis unspecified, setting the timespan option to [0, %f].\n",
                        options.timespan[1]);
            } else if (isInf(options.timespan[0])) {
                options.timespan[0] = 0;
                System.out.format(
                        "Start time of transient analysis unspecified, setting the timespan option to [0, %f].\n",
                        options.timespan[1]);
            } else if (isInf(options.timespan[1])) {
                options.timespan[1] = 30 / minRate;
                System.out.format(
                        "End time of transient analysis unspecified, setting the timespan option to [%f, %f].\n",
                        options.timespan[0], options.timespan[1]);
            }
            try {
                runAnalyzer();
            } catch (IllegalAccessException e) {
                line_error(mfilename(new Object[]{}), "IllegalAccessException upon running runAnalyzer()");
            } catch (ParserConfigurationException e) {
                line_error(mfilename(new Object[]{}), "ParserConfigurationException upon running runAnalyzer()");
            } catch (IOException e) {
                line_error(mfilename(new Object[]{}), "IOException upon running runAnalyzer()");
            }
        }

        // Note: Metrics storage functionality not yet implemented
        // This would involve storing computed transient metrics for later retrieval
        // Implementation depends on the specific metrics storage system to be defined
    }

    /**
     * Returns cumulative distribution functions of passage times during transient analysis.
     * This is an abstract method that must be implemented by concrete solver subclasses.
     *
     * @param R response time handles (optional)
     * @return result containing transient CDFs for passage times
     * @throws RuntimeException if not implemented by the concrete solver
     */
    public DistributionResult getTranCdfPassT(AvgHandle R) {
        throw new RuntimeException("getTranCdfPassT is not supported by " + this.getClass().getSimpleName());
    }

    /**
     * Returns cumulative distribution functions of passage times during transient analysis.
     * Uses default response time handles.
     *
     * @return result containing transient CDFs for passage times
     * @throws RuntimeException if not implemented by the concrete solver
     */
    public DistributionResult getTranCdfPassT() {
        return getTranCdfPassT(getAvgRespTHandles());
    }

    /**
     * Returns cumulative distribution functions of response times during transient analysis.
     * This is an abstract method that must be implemented by concrete solver subclasses.
     *
     * @param R response time handles (optional)
     * @return result containing transient CDFs for response times
     * @throws RuntimeException if not implemented by the concrete solver
     */
    public DistributionResult getTranCdfRespT(AvgHandle R) {
        throw new RuntimeException("getTranCdfRespT is not supported by " + this.getClass().getSimpleName());
    }

    /**
     * Returns cumulative distribution functions of response times during transient analysis.
     * Uses default response time handles.
     *
     * @return result containing transient CDFs for response times
     * @throws RuntimeException if not implemented by the concrete solver
     */
    public DistributionResult getTranCdfRespT() {
        return getTranCdfRespT(getAvgRespTHandles());
    }

    /**
     * Returns the transient performance metric handles.
     *
     * @return the transient handles object
     */
    public SolverTranHandles getTranHandles() {
        return this.tranHandles;
    }

    /**
     * Sets the transient performance metric handles.
     *
     * @param handles the transient handles to set
     */
    public void setTranHandles(SolverTranHandles handles) {
        this.tranHandles = handles;
    }

    /**
     * Checks if the solver has computed steady-state average metrics.
     *
     * @return true if steady-state results are available, false otherwise
     */
    protected boolean hasAvgResults() {
        return !(result == null) && (!((result.QN == null || result.QN.isEmpty()) &&
                (result.UN == null || result.UN.isEmpty()) &&
                (result.RN == null || result.RN.isEmpty()) &&
                (result.TN == null || result.TN.isEmpty()) &&
                (result.CN == null || result.CN.isEmpty()) &&
                (result.XN == null || result.XN.isEmpty())));
    }

    /**
     * Checks if the solver has computed steady-state distribution metrics.
     *
     * @return true if distribution results are available
     */
    public boolean hasDistribResults() {
        // Check if any distribution methods have been implemented and can return results
        try {
            DistributionResult cdfResult = getCdfRespT();
            return cdfResult != null && cdfResult.cdfData != null && !cdfResult.cdfData.isEmpty();
        } catch (RuntimeException e) {
            // Distribution methods not implemented for this solver
            return false;
        }
    }

    /**
     * Helper method to check if a matrix array contains NaN values.
     *
     * @param matrices 2D array of matrices to check
     * @return true if any matrix contains NaN values, false otherwise
     */
    private boolean hasNaNValues(Matrix[][] matrices) {
        if (matrices == null) {
            return false;
        }
        for (Matrix[] matrixRow : matrices) {
            if (matrixRow == null) {
                continue;
            }
            for (Matrix matrix : matrixRow) {
                if (matrix != null && matrix.hasNaN()) {
                    return true;
                }
            }
        }
        return false;
    }

    /**
     * Checks if the solver has computed transient average metrics.
     *
     * @return true if transient results are available, false otherwise
     */
    public boolean hasTranResults() {
        if (this.hasResults()) {
            return result.QNt.length > 0 && result.QNt[0].length > 0 && !result.QNt[0][0].isEmpty();
        }
        return false;
    }

    /**
     * Initializes performance metric handles from the model.
     * This method retrieves and sets both average and transient handles.
     */
    protected void initHandles() {
        this.avgHandles = model.getAvgHandles();
        this.tranHandles = model.getTranHandles();
        // Force model to refresh struct if needed
        this.sn = model.getStruct(true);
    }

    // Model
    public Network model() {
        return getModel();
    }

    /**
     * Detailed print function that displays all contents of the NetworkSolver.
     * This function prints actual values, not reference addresses, and includes ALL fields
     * from both NetworkSolver and its parent Solver class for comprehensive comparison.
     */
    public void print() {
        // NetworkSolver fields
        System.out.println("name: " + (this.name != null ? "\"" + this.name + "\"" : "null"));
        System.out.println("model: " + (this.model != null ? "\"" + this.model.getName() + "\"" : "null"));
        System.out.println("enableChecks: " + this.enableChecks);
        System.out.println("options: " + (this.options != null ? this.options.toString() : "null"));
        System.out.println("result: " + (this.result != null ? "SolverResult[" + this.result.getClass().getSimpleName() + "]" : "null"));
        System.out.println("avgHandles: " + (this.avgHandles != null ? "SolverAvgHandles[" + this.avgHandles.getClass().getSimpleName() + "]" : "null"));
        System.out.println("tranHandles: " + (this.tranHandles != null ? "SolverTranHandles[" + this.tranHandles.getClass().getSimpleName() + "]" : "null"));

        // NetworkStruct (sn) information
        if (this.sn == null) {
            System.out.println("sn: null");
        } else {
            this.sn.print();
        }
    }

    // Probabilities
    public ProbabilityResult prob(int node) {
        return getProb(node);
    }

    public ProbabilityResult prob(int node, Matrix state_a) {
        return getProb(node, state_a);
    }

    public ProbabilityResult probAggr(int node) {
        return getProbAggr(node);
    }

    public ProbabilityResult probAggr(int node, Matrix state_a) {
        return getProbAggr(node, state_a);
    }

    public ProbabilityResult probMarg(int node, int jobclass) {
        return getProbMarg(node, jobclass);
    }

    public ProbabilityResult probMarg(int node, int jobclass, Matrix state_m) {
        return getProbMarg(node, jobclass, state_m);
    }

    public ProbabilityResult probNormConstAggr() {
        return getProbNormConstAggr();
    }

    public ProbabilityResult probSys() {
        return getProbSys();
    }

    public ProbabilityResult probSysAggr() {
        return getProbSysAggr();
    }

    /**
     * Validates model compatibility and method support before analysis.
     *
     * @param options solver options containing method specification
     * @throws RuntimeException if model contains unsupported features or method is invalid
     */
    public void runAnalyzerChecks(SolverOptions options) {
        // Basic model validation - check for empty model
        if (model == null) {
            throw new RuntimeException("Model cannot be null");
        }

        if (model.getNumberOfNodes() == 0) {
            throw new RuntimeException("Model must contain at least one node");
        }

        // Basic options validation
        if (options == null) {
            throw new RuntimeException("SolverOptions cannot be null");
        }

        // Additional model structure validation
        NetworkStruct sn = model.getStruct(false);
        if (sn == null) {
            throw new RuntimeException("Unable to obtain model structure");
        }

        // Check for finite timespan values if applicable
        if (options.timespan != null && options.timespan.length >= 2) {
            if (options.timespan[0] < 0 || (Double.isFinite(options.timespan[1]) && options.timespan[1] <= options.timespan[0])) {
                throw new RuntimeException("Invalid timespan configuration: start time must be non-negative and end time must be greater than start time");
            }
        }
    }

    /**
     * Samples state trajectories for a specific node.
     * This is an abstract method that must be implemented by concrete solver subclasses.
     *
     * @param node      the node index to sample from
     * @param numEvents the number of events to sample
     * @return result containing sampled state trajectories
     * @throws RuntimeException if not implemented by the concrete solver
     */
    public SampleResult sample(int node, int numEvents) {
        throw new RuntimeException("sample is not supported by " + this.getClass().getSimpleName());
    }

    /**
     * Samples aggregated state trajectories for a specific node.
     * This is an abstract method that must be implemented by concrete solver subclasses.
     *
     * @param node      the node index to sample from
     * @param numEvents the number of events to sample
     * @return result containing sampled aggregated state trajectories
     * @throws RuntimeException if not implemented by the concrete solver
     */
    public SampleResult sampleAggr(int node, int numEvents) {
        throw new RuntimeException("sampleAggr is not supported by " + this.getClass().getSimpleName());
    }

    /**
     * Samples joint system state trajectories.
     * This is an abstract method that must be implemented by concrete solver subclasses.
     *
     * @param numEvents the number of events to sample
     * @return result containing sampled joint system state trajectories
     * @throws RuntimeException if not implemented by the concrete solver
     */
    public SampleResult sampleSys(int numEvents) {
        throw new RuntimeException("sampleSys is not supported by " + this.getClass().getSimpleName());
    }

    /**
     * Samples aggregated joint system state trajectories.
     * This is an abstract method that must be implemented by concrete solver subclasses.
     *
     * @param numEvents the number of events to sample
     * @return result containing sampled aggregated joint system state trajectories
     * @throws RuntimeException if not implemented by the concrete solver
     */
    public SampleResult sampleSysAggr(int numEvents) {
        throw new RuntimeException("sampleSysAggr is not supported by " + this.getClass().getSimpleName());
    }

    /**
     * Stores computed average metrics at steady-state in the solver result.
     *
     * @param Q       queue length matrix [stations x classes]
     * @param U       utilization matrix [stations x classes]
     * @param R       response time matrix [stations x classes]
     * @param T       throughput matrix [stations x classes]
     * @param A       arrival rate matrix [stations x classes]
     * @param W       residence time matrix [stations x classes]
     * @param C       system response time vector [chains]
     * @param X       system throughput vector [chains]
     * @param runtime computation time in seconds
     * @param method  solution algorithm used
     * @param iter    number of iterations performed
     */
    public void setAvgResults(Matrix Q, Matrix U, Matrix R, Matrix T, Matrix A, Matrix W, Matrix C, Matrix X,
                              double runtime, String method, int iter) {
        this.result.solver = this.getName();
        this.result.QN = Q.copy();
        this.result.UN = U.copy();
        this.result.RN = R.copy();
        this.result.TN = T.copy();
        this.result.AN = A.copy();
        this.result.WN = W.copy();
        this.result.CN = C.copy();
        this.result.XN = X.copy();
        this.result.runtime = runtime;
        this.result.method = method;

        if (this.options.verbose != VerboseLevel.SILENT) {
            if (iter <= 1) {
                System.out.printf(
                        "%s analysis [method: %s, lang: %s, env: %s] completed in %fs.\n",
                        this.name.replaceFirst("^Solver", ""),   // solver name with prefix stripped
                        this.result.method,                      // algorithm/method
                        "java",                                  // language label
                        System.getProperty("java.version"),      // actual JVM version in use
                        this.result.runtime                      // elapsed time in seconds
                );
            } else {
                System.out.printf(
                        "%s analysis [method: %s, lang: %s, env: %s] completed in %fs. Iterations: %d.\n",
                        this.name.replaceFirst("^Solver", ""),   // solver name with prefix stripped
                        this.result.method,                      // algorithm/method
                        "java",                                  // language label
                        System.getProperty("java.version"),      // actual JVM version in use
                        this.result.runtime,                     // elapsed time in seconds
                        iter                                     // iteration count
                );
            }
            System.out.flush();
        }
    }

    /**
     * Stores distribution metrics at steady-state.
     *
     * @param RD      distribution data [stations x classes] containing CDFs
     * @param runtime computation time for distribution analysis
     */
    protected final void setDistribResults(Matrix RD, double runtime) {
        this.result.solver = this.getName();
        //this.result.RD = RD.copy();
        this.result.runtime = runtime;
    }

    /**
     * Sets the language for solver execution.
     * This method configures the solver to use either MATLAB or Java backend.
     */
    protected void setLang() {
        // Java implementation always uses Java backend
        // This method is provided for consistency with MATLAB interface
    }

    /**
     * Stores computed transient average metrics.
     *
     * @param Qt       transient queue length matrices [time][stations x classes]
     * @param Ut       transient utilization matrices [time][stations x classes]
     * @param Rt       transient response time matrices [time][stations x classes]
     * @param Tt       transient throughput matrices [time][stations x classes]
     * @param Ct       transient system response time matrices [time][chains]
     * @param Xt       transient system throughput matrices [time][chains]
     * @param runtimet computation time for transient analysis
     */
    public void setTranAvgResults(Matrix[][] Qt, Matrix[][] Ut, Matrix[][] Rt, Matrix[][] Tt, Matrix[][]
            Ct, Matrix[][] Xt, double runtimet) {
        this.result.solver = getName();
        this.result.method = getOptions().method;
        // NaN values in transient metrics are normal (e.g., throughput at time 0
        // when no jobs exist at a station), so store the arrays as-is
        this.result.QNt = Qt != null ? Qt.clone() : new Matrix[0][0];
        this.result.UNt = Ut != null ? Ut.clone() : new Matrix[0][0];
        this.result.RNt = Rt != null ? Rt.clone() : new Matrix[0][0];
        this.result.TNt = Tt != null ? Tt.clone() : new Matrix[0][0];
        this.result.CNt = Ct != null ? Ct.clone() : new Matrix[0][0];
        this.result.XNt = Xt != null ? Xt.clone() : new Matrix[0][0];
        this.result.runtime = runtimet;
    }

    /**
     * Stores transient probability distributions.
     *
     * @param t        time points vector
     * @param pi_t     transient probability distributions
     * @param SS       steady-state distribution
     * @param runtimet computation time for transient analysis
     */
    protected final void setTranProb(Matrix t, Matrix pi_t, Matrix SS, double runtimet) {
        this.result.solver = getName();
        this.result.method = getOptions().method;
        this.result.t = t.copy();
        this.result.pi_t = pi_t.copy();
        this.result.SS = SS.copy();
        this.result.runtime = runtimet;
    }

    // Stage table
    public Object stageTable() {
        return getStageTable();
    }

    // stageTable -> stageT aliases
    public Object stageT() {
        return stageTable();
    }

    public NetworkAvgTable getStageT() {
        return getStageTable();
    }

    public NetworkAvgTable getStageT(boolean keepDisabled) {
        return getStageTable(keepDisabled);
    }

    // Transient analysis
    public void tranAvg() {
        getTranAvg();
    }

    public DistributionResult tranCdfPassT() {
        return getTranCdfPassT();
    }

    public DistributionResult tranCdfPassT(AvgHandle R) {
        return getTranCdfPassT(R);
    }

    public DistributionResult tranCdfRespT() {
        return getTranCdfRespT();
    }

    public DistributionResult tranCdfRespT(AvgHandle R) {
        return getTranCdfRespT(R);
    }

    public SolverTranHandles tranHandles() {
        return getTranHandles();
    }

}
