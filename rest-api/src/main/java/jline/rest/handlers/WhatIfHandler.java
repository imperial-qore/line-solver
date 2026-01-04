/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.rest.handlers;

import jline.lang.Model;
import jline.lang.Network;
import jline.lang.JobClass;
import jline.lang.nodes.*;
import jline.lang.processes.Exp;
import jline.rest.model.*;
import jline.rest.model.WhatIfResponse.ParameterPoint;
import jline.rest.model.SensitivityResponse.ParameterSensitivity;
import jline.rest.model.BottleneckResponse.Bottleneck;
import jline.rest.model.BottleneckResponse.StationMetrics;
import jline.rest.util.ModelParser;
import jline.rest.util.ModelParser.ModelParseException;
import jline.solvers.*;

import java.util.*;
import java.util.concurrent.*;

/**
 * Handler for what-if analysis, sensitivity analysis, and bottleneck detection.
 */
public class WhatIfHandler {

    private final ModelParser modelParser;
    private final SolveHandler solveHandler;
    private final ExecutorService executor;

    public WhatIfHandler() {
        this.modelParser = new ModelParser();
        this.solveHandler = new SolveHandler();
        this.executor = Executors.newFixedThreadPool(
            Runtime.getRuntime().availableProcessors()
        );
    }

    /**
     * Execute a what-if parameter sweep analysis.
     */
    public WhatIfResponse whatIf(WhatIfRequest request) {
        String validationError = request.validate();
        if (validationError != null) {
            WhatIfResponse response = new WhatIfResponse();
            response.setStatus("failed");
            response.setError(validationError);
            return response;
        }

        long startTime = System.currentTimeMillis();
        List<ParameterSweep> parameters = request.getParameters();

        // Generate all parameter combinations
        List<Map<String, Double>> combinations = generateCombinations(parameters);
        int totalPoints = combinations.size();

        WhatIfResponse response = new WhatIfResponse();
        response.setSolver(request.getSolver());
        response.setTotalPoints(totalPoints);

        List<ParameterPoint> results = new ArrayList<ParameterPoint>();
        int completed = 0;
        int failed = 0;

        // Execute each combination
        boolean parallel = request.getParallel() != null && request.getParallel();

        if (parallel && totalPoints > 1) {
            // Parallel execution
            List<Future<ParameterPoint>> futures = new ArrayList<Future<ParameterPoint>>();
            for (final Map<String, Double> combo : combinations) {
                futures.add(executor.submit(new Callable<ParameterPoint>() {
                    @Override
                    public ParameterPoint call() {
                        return evaluatePoint(request, combo);
                    }
                }));
            }

            for (Future<ParameterPoint> future : futures) {
                try {
                    ParameterPoint point = future.get(300, TimeUnit.SECONDS);
                    results.add(point);
                    if (point.getError() == null) {
                        completed++;
                    } else {
                        failed++;
                    }
                } catch (Exception e) {
                    ParameterPoint errorPoint = new ParameterPoint();
                    errorPoint.setError("Execution error: " + e.getMessage());
                    results.add(errorPoint);
                    failed++;
                }
            }
        } else {
            // Sequential execution
            for (Map<String, Double> combo : combinations) {
                ParameterPoint point = evaluatePoint(request, combo);
                results.add(point);
                if (point.getError() == null) {
                    completed++;
                } else {
                    failed++;
                }
            }
        }

        response.setResults(results);
        response.setCompletedPoints(completed);
        response.setFailedPoints(failed);
        response.setRuntime((System.currentTimeMillis() - startTime) / 1000.0);
        response.setStatus(failed == totalPoints ? "failed" : "completed");

        return response;
    }

    /**
     * Execute sensitivity analysis.
     */
    public SensitivityResponse sensitivity(SensitivityRequest request) {
        String validationError = request.validate();
        if (validationError != null) {
            SensitivityResponse response = new SensitivityResponse();
            response.setStatus("failed");
            response.setError(validationError);
            return response;
        }

        long startTime = System.currentTimeMillis();
        double delta = request.getEffectiveDelta();
        boolean normalized = request.isNormalized();

        SensitivityResponse response = new SensitivityResponse();
        response.setSolver(request.getSolver());
        response.setDelta(delta);
        response.setNormalized(normalized);

        try {
            // Get baseline metrics
            SolveRequest baselineRequest = createSolveRequest(request.getModel(),
                request.getSolver(), request.getAnalysis(), request.getOptions());
            SolveResponse baselineResponse = solveHandler.solve(baselineRequest);

            if (!"completed".equals(baselineResponse.getStatus())) {
                response.setStatus("failed");
                response.setError("Baseline solve failed: " + baselineResponse.getError());
                return response;
            }

            Map<String, Object> baselineMetrics = extractMetricsMap(baselineResponse);
            response.setBaselineMetrics(baselineMetrics);

            // Compute sensitivities for each parameter
            List<ParameterSensitivity> sensitivities = new ArrayList<ParameterSensitivity>();

            for (ParameterSweep param : request.getParameters()) {
                ParameterSensitivity sens = computeSensitivity(
                    request, param, baselineMetrics, delta, normalized);
                sensitivities.add(sens);
            }

            response.setSensitivities(sensitivities);
            response.setStatus("completed");

        } catch (Exception e) {
            response.setStatus("failed");
            response.setError("Sensitivity analysis error: " + e.getMessage());
        }

        response.setRuntime((System.currentTimeMillis() - startTime) / 1000.0);
        return response;
    }

    /**
     * Execute bottleneck detection.
     */
    public BottleneckResponse bottleneck(BottleneckRequest request) {
        String validationError = request.validate();
        if (validationError != null) {
            BottleneckResponse response = new BottleneckResponse();
            response.setStatus("failed");
            response.setError(validationError);
            return response;
        }

        long startTime = System.currentTimeMillis();
        double threshold = request.getEffectiveThreshold();

        BottleneckResponse response = new BottleneckResponse();
        response.setSolver(request.getSolver());
        response.setThreshold(threshold);

        try {
            // Solve the model
            SolveRequest solveRequest = createSolveRequest(request.getModel(),
                request.getSolver(), "all", request.getOptions());
            SolveResponse solveResponse = solveHandler.solve(solveRequest);

            if (!"completed".equals(solveResponse.getStatus())) {
                response.setStatus("failed");
                response.setError("Solve failed: " + solveResponse.getError());
                return response;
            }

            // Extract station metrics and identify bottlenecks
            List<Bottleneck> bottlenecks = new ArrayList<Bottleneck>();
            List<StationMetrics> allStations = new ArrayList<StationMetrics>();

            SolveResponse.ResultTables results = solveResponse.getResults();
            if (results != null && results.getAvgTable() != null) {
                SolveResponse.AvgTableResult avgTable = results.getAvgTable();
                String[] stations = avgTable.getStations();
                String[] classes = avgTable.getClasses();
                Map<String, double[][]> metrics = avgTable.getMetrics();

                double[][] util = metrics.get("Util");
                double[][] qlen = metrics.get("QLen");
                double[][] respT = metrics.get("RespT");
                double[][] tput = metrics.get("Tput");

                if (stations != null && util != null) {
                    for (int i = 0; i < stations.length; i++) {
                        // Aggregate utilization across classes
                        double totalUtil = 0;
                        double totalQLen = 0;
                        double totalRespT = 0;
                        double totalTput = 0;
                        Map<String, Double> classBreakdown = new HashMap<String, Double>();

                        for (int k = 0; k < classes.length; k++) {
                            if (util != null && i < util.length && k < util[i].length) {
                                double u = util[i][k];
                                if (!Double.isNaN(u)) {
                                    totalUtil += u;
                                    classBreakdown.put(classes[k], u);
                                }
                            }
                            if (qlen != null && i < qlen.length && k < qlen[i].length) {
                                double q = qlen[i][k];
                                if (!Double.isNaN(q)) totalQLen += q;
                            }
                            if (respT != null && i < respT.length && k < respT[i].length) {
                                double r = respT[i][k];
                                if (!Double.isNaN(r)) totalRespT += r;
                            }
                            if (tput != null && i < tput.length && k < tput[i].length) {
                                double t = tput[i][k];
                                if (!Double.isNaN(t)) totalTput += t;
                            }
                        }

                        StationMetrics sm = new StationMetrics();
                        sm.setStation(stations[i]);
                        sm.setUtilization(totalUtil);
                        sm.setQueueLength(totalQLen);
                        sm.setResponseTime(totalRespT);
                        sm.setThroughput(totalTput);
                        sm.setBottleneck(totalUtil >= threshold);
                        allStations.add(sm);

                        // Check if bottleneck
                        if (totalUtil >= threshold) {
                            Bottleneck bn = new Bottleneck();
                            bn.setStation(stations[i]);
                            bn.setUtilization(totalUtil);
                            bn.setQueueLength(totalQLen);
                            bn.setResponseTime(totalRespT);
                            bn.setClassBreakdown(classBreakdown);

                            // Determine severity
                            if (totalUtil >= 0.95) {
                                bn.setSeverity("critical");
                            } else if (totalUtil >= 0.9) {
                                bn.setSeverity("high");
                            } else {
                                bn.setSeverity("moderate");
                            }

                            bottlenecks.add(bn);
                        }
                    }
                }

                // Sort bottlenecks by utilization (descending)
                Collections.sort(bottlenecks, new Comparator<Bottleneck>() {
                    @Override
                    public int compare(Bottleneck a, Bottleneck b) {
                        return Double.compare(b.getUtilization(), a.getUtilization());
                    }
                });
            }

            response.setBottlenecks(bottlenecks);
            if (request.shouldIncludeAll()) {
                response.setAllStations(allStations);
            }

            // Extract system metrics
            if (results != null && results.getSysTable() != null) {
                Map<String, Object> sysMetrics = new HashMap<String, Object>();
                SolveResponse.AvgSysTableResult sysTable = results.getSysTable();
                if (sysTable.getMetrics() != null) {
                    for (Map.Entry<String, double[]> entry : sysTable.getMetrics().entrySet()) {
                        sysMetrics.put(entry.getKey(), entry.getValue());
                    }
                }
                response.setSystemMetrics(sysMetrics);
            }

            response.setStatus("completed");

        } catch (Exception e) {
            response.setStatus("failed");
            response.setError("Bottleneck analysis error: " + e.getMessage());
        }

        response.setRuntime((System.currentTimeMillis() - startTime) / 1000.0);
        return response;
    }

    /**
     * Generate all combinations of parameter values.
     */
    private List<Map<String, Double>> generateCombinations(List<ParameterSweep> parameters) {
        List<Map<String, Double>> result = new ArrayList<Map<String, Double>>();
        result.add(new HashMap<String, Double>());

        for (ParameterSweep param : parameters) {
            List<Double> values = param.getEffectiveValues();
            String key = param.getDescription();

            List<Map<String, Double>> newResult = new ArrayList<Map<String, Double>>();
            for (Map<String, Double> existing : result) {
                for (Double value : values) {
                    Map<String, Double> combo = new HashMap<String, Double>(existing);
                    combo.put(key, value);
                    newResult.add(combo);
                }
            }
            result = newResult;
        }

        return result;
    }

    /**
     * Evaluate a single parameter point by modifying the model and solving.
     */
    private ParameterPoint evaluatePoint(WhatIfRequest request, Map<String, Double> parameters) {
        ParameterPoint point = new ParameterPoint();
        point.setParameters(parameters);

        try {
            // Parse the base model
            Model model = modelParser.parse(request.getModel());

            if (!(model instanceof Network)) {
                point.setError("What-if analysis requires a queueing network model");
                return point;
            }

            Network network = (Network) model;

            // Apply parameter modifications
            for (ParameterSweep param : request.getParameters()) {
                String key = param.getDescription();
                Double value = parameters.get(key);
                if (value != null) {
                    applyParameterToNetwork(network, param, value);
                }
            }

            // Solve the modified model
            SolveRequest solveRequest = new SolveRequest();
            solveRequest.setSolver(request.getSolver());
            solveRequest.setAnalysis(request.getAnalysis());
            solveRequest.setOptions(request.getOptions());

            // We need to solve the already-parsed network directly
            SolveResponse solveResponse = solveNetworkDirectly(network, request.getSolver(),
                request.getAnalysis(), request.getOptions());

            if ("completed".equals(solveResponse.getStatus())) {
                point.setMetrics(extractMetricsMap(solveResponse));
            } else {
                point.setError(solveResponse.getError());
            }
        } catch (Exception e) {
            point.setError("Evaluation error: " + e.getMessage());
        }

        return point;
    }

    /**
     * Apply a parameter value to a network.
     */
    private void applyParameterToNetwork(Network network, ParameterSweep param, double value) {
        String property = param.getProperty();
        String stationName = param.getStation();
        String className = param.getClassName();

        // Find the station
        Node node = network.getNodeByName(stationName);
        if (node == null) {
            throw new RuntimeException("Station not found: " + stationName);
        }

        // Find the job class (if specified)
        JobClass jobClass = null;
        if (className != null && !className.isEmpty()) {
            jobClass = network.getClassByName(className);
            if (jobClass == null) {
                throw new RuntimeException("Class not found: " + className);
            }
        }

        // Apply the parameter based on property type
        switch (property.toLowerCase()) {
            case "servicerate":
                applyServiceRate(node, jobClass, value);
                break;
            case "arrivalrate":
                applyArrivalRate(node, jobClass, value);
                break;
            case "thinktime":
                applyThinkTime(node, jobClass, value);
                break;
            case "servers":
                applyNumberOfServers(node, (int) value);
                break;
            default:
                throw new RuntimeException("Unsupported property: " + property);
        }

        // Reset the network structure to force rebuild with new parameters
        network.resetStruct();
    }

    /**
     * Apply service rate change to a station.
     */
    private void applyServiceRate(Node node, JobClass jobClass, double rate) {
        if (node instanceof jline.lang.nodes.Queue) {
            jline.lang.nodes.Queue queue = (jline.lang.nodes.Queue) node;
            // Create exponential distribution with the given rate
            Exp exp = new Exp(rate);
            if (jobClass != null) {
                queue.setService(jobClass, exp);
            } else {
                // Apply to all classes
                for (JobClass jc : queue.getModel().getClasses()) {
                    if (queue.containsJobClass(jc)) {
                        queue.setService(jc, exp);
                    }
                }
            }
        } else if (node instanceof Delay) {
            Delay delay = (Delay) node;
            Exp exp = new Exp(rate);
            if (jobClass != null) {
                delay.setService(jobClass, exp);
            } else {
                for (JobClass jc : delay.getModel().getClasses()) {
                    delay.setService(jc, exp);
                }
            }
        } else {
            throw new RuntimeException("Cannot set service rate on node type: " + node.getClass().getSimpleName());
        }
    }

    /**
     * Apply arrival rate change to a source.
     */
    private void applyArrivalRate(Node node, JobClass jobClass, double rate) {
        if (node instanceof Source) {
            Source source = (Source) node;
            Exp exp = new Exp(rate);
            if (jobClass != null) {
                source.setArrival(jobClass, exp);
            } else {
                for (JobClass jc : source.getModel().getClasses()) {
                    source.setArrival(jc, exp);
                }
            }
        } else {
            throw new RuntimeException("Cannot set arrival rate on node type: " + node.getClass().getSimpleName());
        }
    }

    /**
     * Apply think time change (as mean service time at a delay station).
     */
    private void applyThinkTime(Node node, JobClass jobClass, double thinkTime) {
        if (node instanceof Delay) {
            Delay delay = (Delay) node;
            // Think time is the mean, so rate = 1/thinkTime
            double rate = 1.0 / thinkTime;
            Exp exp = new Exp(rate);
            if (jobClass != null) {
                delay.setService(jobClass, exp);
            } else {
                for (JobClass jc : delay.getModel().getClasses()) {
                    delay.setService(jc, exp);
                }
            }
        } else {
            throw new RuntimeException("Cannot set think time on node type: " + node.getClass().getSimpleName());
        }
    }

    /**
     * Apply number of servers change.
     */
    private void applyNumberOfServers(Node node, int servers) {
        if (node instanceof Station) {
            ((Station) node).setNumberOfServers(servers);
        } else {
            throw new RuntimeException("Cannot set number of servers on node type: " + node.getClass().getSimpleName());
        }
    }

    /**
     * Solve a network directly (already parsed).
     */
    private SolveResponse solveNetworkDirectly(Network network, String solverType,
            String analysis, SolverOptionsDTO optionsDto) {
        if (analysis == null || analysis.isEmpty()) {
            analysis = "all";
        }

        long startTime = System.currentTimeMillis();

        try {
            SolverOptions options = optionsDto != null ? optionsDto.toSolverOptions() : new SolverOptions();
            Solver solver = createSolver(network, solverType, options);

            if (solver == null) {
                return SolveResponse.failure(solverType, "Failed to create solver");
            }

            boolean wantAvgTable = "all".equals(analysis) || "avg".equals(analysis);
            boolean wantSysTable = "all".equals(analysis) || "sys".equals(analysis);

            SolveResponse.ResultTables results = new SolveResponse.ResultTables();

            if (wantAvgTable) {
                AvgTable avgTable = getAvgTable(solver, network);
                if (avgTable != null) {
                    results.setAvgTable(convertAvgTable(avgTable, network));
                }
            }

            if (wantSysTable && solver instanceof NetworkSolver) {
                NetworkSolver netSolver = (NetworkSolver) solver;
                try {
                    NetworkAvgSysTable sysTable = netSolver.getAvgSysTable();
                    if (sysTable != null) {
                        results.setSysTable(convertSysTable(sysTable, network));
                    }
                } catch (Exception e) {
                    // Sys table may not be available
                }
            }

            double runtime = (System.currentTimeMillis() - startTime) / 1000.0;
            return SolveResponse.success(solverType, runtime, results);

        } catch (Exception e) {
            return SolveResponse.failure(solverType, "Solver error: " + e.getMessage());
        }
    }

    /**
     * Create a solver for the network.
     */
    private Solver createSolver(Network network, String solverType, SolverOptions options) {
        switch (solverType.toLowerCase()) {
            case "mva":
                return new jline.solvers.mva.SolverMVA(network, options);
            case "ctmc":
                options.force(true);
                return new jline.solvers.ctmc.SolverCTMC(network, options);
            case "fluid":
                return new jline.solvers.fluid.SolverFluid(network, options);
            case "jmt":
                return new jline.solvers.jmt.SolverJMT(network, options);
            case "nc":
                return new jline.solvers.nc.SolverNC(network, options);
            case "ssa":
                return new jline.solvers.ssa.SolverSSA(network, options);
            default:
                return null;
        }
    }

    /**
     * Get average table from solver.
     */
    private AvgTable getAvgTable(Solver solver, Network network) {
        try {
            if (solver instanceof NetworkSolver) {
                return ((NetworkSolver) solver).getAvgTable();
            }
        } catch (Exception e) {
            // Table computation may fail
        }
        return null;
    }

    /**
     * Convert NetworkAvgTable to REST API format.
     */
    private SolveResponse.AvgTableResult convertAvgTable(AvgTable avgTable, Network network) {
        SolveResponse.AvgTableResult result = new SolveResponse.AvgTableResult();

        if (!(avgTable instanceof NetworkAvgTable)) {
            return result;
        }

        NetworkAvgTable netAvgTable = (NetworkAvgTable) avgTable;

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

        Map<String, double[][]> metrics = new HashMap<String, double[][]>();
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
    private SolveResponse.AvgSysTableResult convertSysTable(NetworkAvgSysTable sysTable, Network network) {
        SolveResponse.AvgSysTableResult result = new SolveResponse.AvgSysTableResult();

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

        Map<String, double[]> metrics = new HashMap<String, double[]>();
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

    /**
     * Create a solve request from analysis request components.
     */
    private SolveRequest createSolveRequest(ModelInput model, String solver,
            String analysis, SolverOptionsDTO options) {
        SolveRequest request = new SolveRequest();
        request.setModel(model);
        request.setSolver(solver);
        request.setAnalysis(analysis);
        request.setOptions(options);
        return request;
    }

    /**
     * Extract metrics from solve response as a flat map.
     */
    private Map<String, Object> extractMetricsMap(SolveResponse response) {
        Map<String, Object> result = new HashMap<String, Object>();

        if (response.getResults() != null) {
            SolveResponse.ResultTables results = response.getResults();

            if (results.getAvgTable() != null) {
                result.put("avgTable", results.getAvgTable());
            }
            if (results.getSysTable() != null) {
                result.put("sysTable", results.getSysTable());
            }
        }

        return result;
    }

    /**
     * Compute sensitivity of metrics with respect to a parameter using finite differences.
     */
    private ParameterSensitivity computeSensitivity(SensitivityRequest request,
            ParameterSweep param, Map<String, Object> baselineMetrics,
            double delta, boolean normalized) {

        ParameterSensitivity sens = new ParameterSensitivity();
        sens.setParameter(param.getDescription());

        try {
            // Get the base parameter value from the model
            Model baseModel = modelParser.parse(request.getModel());
            if (!(baseModel instanceof Network)) {
                sens.setDerivatives(new HashMap<String, Object>());
                sens.setBaseValue(0.0);
                return sens;
            }

            Network baseNetwork = (Network) baseModel;
            double baseValue = getParameterValue(baseNetwork, param);
            sens.setBaseValue(baseValue);

            // Compute perturbed values
            double perturbedValuePlus = baseValue * (1.0 + delta);
            double perturbedValueMinus = baseValue * (1.0 - delta);

            // Solve with perturbed values
            Map<String, Object> metricsPlus = solveWithParameterValue(request, param, perturbedValuePlus);
            Map<String, Object> metricsMinus = solveWithParameterValue(request, param, perturbedValueMinus);

            // Compute numerical derivatives using central differences
            Map<String, Object> derivatives = computeNumericalDerivatives(
                baselineMetrics, metricsPlus, metricsMinus, baseValue, delta, normalized);

            sens.setDerivatives(derivatives);

        } catch (Exception e) {
            sens.setDerivatives(new HashMap<String, Object>());
            sens.setBaseValue(0.0);
        }

        return sens;
    }

    /**
     * Get the current value of a parameter from the network.
     */
    private double getParameterValue(Network network, ParameterSweep param) {
        String property = param.getProperty();
        String stationName = param.getStation();
        String className = param.getClassName();

        Node node = network.getNodeByName(stationName);
        if (node == null) {
            return 1.0; // Default
        }

        JobClass jobClass = null;
        if (className != null && !className.isEmpty()) {
            jobClass = network.getClassByName(className);
        }

        switch (property.toLowerCase()) {
            case "servicerate":
                if (node instanceof jline.lang.nodes.Queue) {
                    jline.lang.nodes.Queue queue = (jline.lang.nodes.Queue) node;
                    if (jobClass != null && queue.containsJobClass(jobClass)) {
                        return queue.getService(jobClass).getRate();
                    }
                } else if (node instanceof Delay) {
                    Delay delay = (Delay) node;
                    if (jobClass != null) {
                        return delay.getService(jobClass).getRate();
                    }
                }
                return 1.0;
            case "servers":
                if (node instanceof Station) {
                    return ((Station) node).getNumberOfServers();
                }
                return 1.0;
            default:
                return 1.0;
        }
    }

    /**
     * Solve the model with a specific parameter value.
     */
    private Map<String, Object> solveWithParameterValue(SensitivityRequest request,
            ParameterSweep param, double value) throws ModelParseException {

        Model model = modelParser.parse(request.getModel());
        if (!(model instanceof Network)) {
            return new HashMap<String, Object>();
        }

        Network network = (Network) model;
        applyParameterToNetwork(network, param, value);

        SolveResponse response = solveNetworkDirectly(network, request.getSolver(),
            request.getAnalysis(), request.getOptions());

        return extractMetricsMap(response);
    }

    /**
     * Compute numerical derivatives using central differences.
     */
    private Map<String, Object> computeNumericalDerivatives(
            Map<String, Object> baselineMetrics,
            Map<String, Object> metricsPlus,
            Map<String, Object> metricsMinus,
            double baseValue,
            double delta,
            boolean normalized) {

        Map<String, Object> derivatives = new HashMap<String, Object>();
        double h = 2.0 * baseValue * delta; // Total perturbation for central differences

        // Extract avgTable derivatives
        SolveResponse.AvgTableResult baseAvg = (SolveResponse.AvgTableResult) baselineMetrics.get("avgTable");
        SolveResponse.AvgTableResult plusAvg = (SolveResponse.AvgTableResult) metricsPlus.get("avgTable");
        SolveResponse.AvgTableResult minusAvg = (SolveResponse.AvgTableResult) metricsMinus.get("avgTable");

        if (baseAvg != null && plusAvg != null && minusAvg != null) {
            Map<String, double[][]> avgDerivatives = new HashMap<String, double[][]>();

            for (String metricName : baseAvg.getMetrics().keySet()) {
                double[][] baseData = baseAvg.getMetrics().get(metricName);
                double[][] plusData = plusAvg.getMetrics().get(metricName);
                double[][] minusData = minusAvg.getMetrics().get(metricName);

                if (baseData != null && plusData != null && minusData != null) {
                    int rows = baseData.length;
                    int cols = baseData[0].length;
                    double[][] deriv = new double[rows][cols];

                    for (int i = 0; i < rows; i++) {
                        for (int j = 0; j < cols; j++) {
                            // Central difference: (f(x+h) - f(x-h)) / (2h)
                            double derivative = (plusData[i][j] - minusData[i][j]) / h;

                            // Normalize if requested (elasticity)
                            if (normalized && !Double.isNaN(baseData[i][j]) && baseData[i][j] != 0) {
                                derivative = derivative * baseValue / baseData[i][j];
                            }

                            deriv[i][j] = derivative;
                        }
                    }
                    avgDerivatives.put(metricName, deriv);
                }
            }

            derivatives.put("avgTable", avgDerivatives);
        }

        return derivatives;
    }

    /**
     * Shutdown the handler's executor.
     */
    public void shutdown() {
        executor.shutdown();
        try {
            if (!executor.awaitTermination(5, TimeUnit.SECONDS)) {
                executor.shutdownNow();
            }
        } catch (InterruptedException e) {
            executor.shutdownNow();
            Thread.currentThread().interrupt();
        }
    }
}
