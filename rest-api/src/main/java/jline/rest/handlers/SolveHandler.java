/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.rest.handlers;

import jline.lang.Model;
import jline.lang.Network;
import jline.lang.layered.LayeredNetwork;
import jline.rest.model.SolveRequest;
import jline.rest.model.SolveResponse;
import jline.rest.model.SolveResponse.AvgTableResult;
import jline.rest.model.SolveResponse.AvgSysTableResult;
import jline.rest.model.SolveResponse.ResultTables;
import jline.rest.model.SolverOptionsDTO;
import jline.rest.util.ModelParser;
import jline.rest.util.ModelParser.ModelParseException;
import jline.solvers.*;
import jline.solvers.ctmc.SolverCTMC;
import jline.solvers.des.SolverDES;
import jline.solvers.fluid.SolverFluid;
import jline.solvers.jmt.SolverJMT;
import jline.solvers.ln.SolverLN;
import jline.solvers.lqns.SolverLQNS;
import jline.solvers.mva.SolverMVA;
import jline.solvers.nc.SolverNC;
import jline.solvers.ssa.SolverSSA;
import jline.util.matrix.Matrix;
import jline.lang.JobClass;
import jline.lang.nodes.StatefulNode;

import java.util.HashMap;
import java.util.Map;

/**
 * Handler for solving queueing network models.
 * Extracts and encapsulates the solve logic from LineCLI for REST API use.
 */
public class SolveHandler {

    private final ModelParser modelParser;

    /**
     * Create a new solve handler.
     */
    public SolveHandler() {
        this.modelParser = new ModelParser();
    }

    /**
     * Solve a model based on the request.
     *
     * @param request The solve request
     * @return The solve response
     */
    public SolveResponse solve(SolveRequest request) {
        // Validate request
        String validationError = request.validate();
        if (validationError != null) {
            return SolveResponse.failure(request.getSolver(), validationError);
        }

        String solver = request.getSolver();
        String analysis = request.getAnalysis();
        if (analysis == null || analysis.isEmpty()) {
            analysis = "all";
        }

        long startTime = System.currentTimeMillis();

        try {
            // Parse model
            Model model = modelParser.parse(request.getModel());

            // Get solver options
            SolverOptions options = getSolverOptions(request);

            // Create and run solver
            Solver solverObj = createSolver(model, solver, options);
            if (solverObj == null) {
                return SolveResponse.failure(solver, "Failed to create solver");
            }

            // Determine what to compute
            boolean wantAvgTable = "all".equals(analysis) || "avg".equals(analysis);
            boolean wantSysTable = "all".equals(analysis) || "sys".equals(analysis);

            // For LQN models, sys table is not available
            if (model instanceof LayeredNetwork) {
                wantSysTable = false;
            }

            // Get results
            ResultTables results = new ResultTables();

            if (wantAvgTable) {
                AvgTable avgTable = getAvgTable(solverObj, model);
                if (avgTable != null) {
                    results.setAvgTable(convertAvgTable(avgTable, model));
                }
            }

            if (wantSysTable && solverObj instanceof NetworkSolver) {
                NetworkSolver netSolver = (NetworkSolver) solverObj;
                try {
                    NetworkAvgSysTable sysTable = netSolver.getAvgSysTable();
                    if (sysTable != null) {
                        results.setSysTable(convertSysTable(sysTable, (Network) model));
                    }
                } catch (Exception e) {
                    // Sys table may not be available for all solvers
                }
            }

            double runtime = (System.currentTimeMillis() - startTime) / 1000.0;
            return SolveResponse.success(solver, runtime, results);

        } catch (ModelParseException e) {
            return SolveResponse.failure(solver, "Model parse error: " + e.getMessage());
        } catch (Exception e) {
            return SolveResponse.failure(solver, "Solver error: " + e.getMessage());
        }
    }

    /**
     * Create a solver instance based on the solver type.
     */
    private Solver createSolver(Model model, String solverType, SolverOptions options) {
        if (model instanceof Network) {
            Network network = (Network) model;
            switch (solverType) {
                case "mva":
                    return new SolverMVA(network, options);
                case "ctmc":
                    options.force(true);
                    return new SolverCTMC(network, options);
                case "fluid":
                    return new SolverFluid(network, options);
                case "jmt":
                    return new SolverJMT(network, options);
                case "nc":
                    return new SolverNC(network, options);
                case "ssa":
                    return new SolverSSA(network, options);
                case "des":
                    return new SolverDES(network, options);
                default:
                    return null;
            }
        } else if (model instanceof LayeredNetwork) {
            LayeredNetwork lqn = (LayeredNetwork) model;
            switch (solverType) {
                case "lqns":
                    return new SolverLQNS(lqn, options);
                case "ln":
                case "mva":
                    return new SolverLN(lqn, options);
                case "nc":
                    return new SolverLN(lqn,
                        (net) -> new SolverNC(net, options), options);
                default:
                    return null;
            }
        }
        return null;
    }

    /**
     * Get the average table from the solver.
     */
    private AvgTable getAvgTable(Solver solver, Model model) {
        try {
            if (solver instanceof NetworkSolver) {
                return ((NetworkSolver) solver).getAvgTable();
            } else if (solver instanceof SolverLQNS) {
                return ((SolverLQNS) solver).getAvgTable();
            } else if (solver instanceof SolverLN) {
                return ((SolverLN) solver).getAvgTable();
            }
        } catch (Exception e) {
            // Table computation may fail
        }
        return null;
    }

    /**
     * Convert solver options from DTO.
     */
    private SolverOptions getSolverOptions(SolveRequest request) {
        SolverOptionsDTO dto = request.getOptions();
        if (dto != null) {
            return dto.toSolverOptions();
        }
        return new SolverOptions();
    }

    /**
     * Convert NetworkAvgTable to REST API format.
     */
    private AvgTableResult convertAvgTable(AvgTable avgTable, Model model) {
        AvgTableResult result = new AvgTableResult();

        if (!(avgTable instanceof NetworkAvgTable) || !(model instanceof Network)) {
            return result;
        }

        NetworkAvgTable netAvgTable = (NetworkAvgTable) avgTable;
        Network network = (Network) model;

        // Get station and class names
        java.util.List<StatefulNode> statefulNodes = network.getStatefulNodes();
        java.util.List<JobClass> jobClasses = network.getClasses();
        int numStations = statefulNodes.size();
        int numClasses = jobClasses.size();

        String[] stations = new String[numStations];
        for (int i = 0; i < numStations; i++) {
            stations[i] = statefulNodes.get(i).getName();
        }
        result.setStations(stations);

        String[] classes = new String[numClasses];
        for (int k = 0; k < numClasses; k++) {
            classes[k] = jobClasses.get(k).getName();
        }
        result.setClasses(classes);

        // Get metrics
        Map<String, double[][]> metrics = new HashMap<String, double[][]>();

        // Metric columns: QLen, Util, RespT, ResidT, ArvR, Tput
        String[] metricNames = {"QLen", "Util", "RespT", "ResidT", "ArvR", "Tput"};
        for (int m = 0; m < metricNames.length; m++) {
            java.util.List<Double> col = netAvgTable.get(m);
            if (col != null && !col.isEmpty()) {
                double[][] metricData = new double[numStations][numClasses];
                for (int i = 0; i < numStations; i++) {
                    for (int k = 0; k < numClasses; k++) {
                        int row = i * numClasses + k;
                        if (row < col.size()) {
                            metricData[i][k] = col.get(row);
                        }
                    }
                }
                metrics.put(metricNames[m], metricData);
            }
        }

        result.setMetrics(metrics);
        return result;
    }

    /**
     * Convert NetworkAvgSysTable to REST API format.
     */
    private AvgSysTableResult convertSysTable(NetworkAvgSysTable sysTable, Network network) {
        AvgSysTableResult result = new AvgSysTableResult();

        // Get chain names from the table if available, otherwise from network
        java.util.List<String> chainNamesList = sysTable.getChainNames();
        if (chainNamesList != null && !chainNamesList.isEmpty()) {
            result.setChains(chainNamesList.toArray(new String[0]));
        } else {
            int numChains = network.getNumberOfChains();
            String[] chains = new String[numChains];
            for (int c = 0; c < numChains; c++) {
                chains[c] = "Chain" + (c + 1);
            }
            result.setChains(chains);
        }

        // Get metrics
        Map<String, double[]> metrics = new HashMap<String, double[]>();

        // System metrics: SysRespT, SysTput
        String[] metricNames = {"SysRespT", "SysTput"};
        for (int m = 0; m < metricNames.length; m++) {
            java.util.List<Double> col = sysTable.get(m);
            if (col != null && !col.isEmpty()) {
                double[] metricData = new double[col.size()];
                for (int c = 0; c < col.size(); c++) {
                    metricData[c] = col.get(c);
                }
                metrics.put(metricNames[m], metricData);
            }
        }

        result.setMetrics(metrics);
        return result;
    }
}
