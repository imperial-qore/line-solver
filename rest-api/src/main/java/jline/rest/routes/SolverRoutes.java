/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.rest.routes;

import jline.rest.util.JsonTransformer;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import static spark.Spark.get;

/**
 * Solver listing endpoint for the REST API.
 * Provides information about available solvers and their capabilities.
 */
public class SolverRoutes {

    /**
     * Register solver routes.
     *
     * @param basePath The API base path (e.g., "/api/v1")
     * @param json The JSON transformer
     */
    public static void register(String basePath, JsonTransformer json) {
        // List all available solvers
        get(basePath + "/solvers", (req, res) -> {
            List<Map<String, Object>> solvers = new ArrayList<Map<String, Object>>();

            // MVA - Mean Value Analysis
            Map<String, Object> mva = new HashMap<String, Object>();
            mva.put("id", "mva");
            mva.put("name", "Mean Value Analysis");
            mva.put("description", "Analytical solver for queueing networks using MVA algorithm");
            mva.put("type", "analytical");
            mva.put("formats", new String[]{"jsim", "jsimg", "jsimw", "lqnx", "xml"});
            mva.put("features", new String[]{"closed", "open", "mixed", "multiclass", "multiserver"});
            mva.put("async", false);
            solvers.add(mva);

            // CTMC - Continuous-Time Markov Chain
            Map<String, Object> ctmc = new HashMap<String, Object>();
            ctmc.put("id", "ctmc");
            ctmc.put("name", "Continuous-Time Markov Chain");
            ctmc.put("description", "Exact state-space analysis using CTMC");
            ctmc.put("type", "analytical");
            ctmc.put("formats", new String[]{"jsim", "jsimg", "jsimw"});
            ctmc.put("features", new String[]{"closed", "open", "mixed", "multiclass", "exact"});
            ctmc.put("async", false);
            ctmc.put("notes", "State space explosion for large models");
            solvers.add(ctmc);

            // Fluid
            Map<String, Object> fluid = new HashMap<String, Object>();
            fluid.put("id", "fluid");
            fluid.put("name", "Fluid Solver");
            fluid.put("description", "Mean-field ODE solver for queueing networks");
            fluid.put("type", "analytical");
            fluid.put("formats", new String[]{"jsim", "jsimg", "jsimw"});
            fluid.put("features", new String[]{"closed", "open", "mixed", "multiclass", "transient"});
            fluid.put("async", false);
            solvers.add(fluid);

            // JMT - Java Modelling Tools
            Map<String, Object> jmt = new HashMap<String, Object>();
            jmt.put("id", "jmt");
            jmt.put("name", "Java Modelling Tools");
            jmt.put("description", "Discrete event simulation using JMT");
            jmt.put("type", "simulation");
            jmt.put("formats", new String[]{"jsim", "jsimg", "jsimw"});
            jmt.put("features", new String[]{"closed", "open", "mixed", "multiclass", "multiserver", "general_distributions"});
            jmt.put("async", true);
            jmt.put("notes", "Requires JMT.jar");
            solvers.add(jmt);

            // NC - Normalizing Constant
            Map<String, Object> nc = new HashMap<String, Object>();
            nc.put("id", "nc");
            nc.put("name", "Normalizing Constant");
            nc.put("description", "Analytical solver using normalizing constant methods");
            nc.put("type", "analytical");
            nc.put("formats", new String[]{"jsim", "jsimg", "jsimw", "lqnx", "xml"});
            nc.put("features", new String[]{"closed", "multiclass", "product_form"});
            nc.put("async", false);
            solvers.add(nc);

            // SSA - Stochastic Simulation Algorithm
            Map<String, Object> ssa = new HashMap<String, Object>();
            ssa.put("id", "ssa");
            ssa.put("name", "Stochastic Simulation Algorithm");
            ssa.put("description", "Event-driven stochastic simulation");
            ssa.put("type", "simulation");
            ssa.put("formats", new String[]{"jsim", "jsimg", "jsimw"});
            ssa.put("features", new String[]{"closed", "open", "mixed", "multiclass", "exact_sampling"});
            ssa.put("async", true);
            solvers.add(ssa);

            // DES - Discrete Event Simulation
            Map<String, Object> des = new HashMap<String, Object>();
            des.put("id", "des");
            des.put("name", "Discrete Event Simulation");
            des.put("description", "Native discrete event simulation using SSJ library");
            des.put("type", "simulation");
            des.put("formats", new String[]{"jsim", "jsimg", "jsimw"});
            des.put("features", new String[]{"open", "multiclass", "multiserver", "general_distributions"});
            des.put("async", true);
            solvers.add(des);

            // LN - Layered Network
            Map<String, Object> ln = new HashMap<String, Object>();
            ln.put("id", "ln");
            ln.put("name", "Layered Network Solver");
            ln.put("description", "Solver for layered queueing networks");
            ln.put("type", "analytical");
            ln.put("formats", new String[]{"lqnx", "xml"});
            ln.put("features", new String[]{"layered", "multiclass", "multiprocessor"});
            ln.put("async", false);
            solvers.add(ln);

            // LQNS - LQN Solver
            Map<String, Object> lqns = new HashMap<String, Object>();
            lqns.put("id", "lqns");
            lqns.put("name", "LQN Solver");
            lqns.put("description", "Layered queueing network solver using LQNS");
            lqns.put("type", "analytical");
            lqns.put("formats", new String[]{"lqnx", "xml"});
            lqns.put("features", new String[]{"layered", "multiclass", "multiprocessor"});
            lqns.put("async", false);
            lqns.put("notes", "External LQNS installation may be required");
            solvers.add(lqns);

            Map<String, Object> response = new HashMap<String, Object>();
            response.put("solvers", solvers);
            response.put("count", solvers.size());

            return response;
        }, json);

        // Get details for a specific solver
        get(basePath + "/solvers/:id", (req, res) -> {
            String solverId = req.params(":id");

            Map<String, Object> solver = getSolverInfo(solverId);
            if (solver == null) {
                res.status(404);
                Map<String, Object> error = new HashMap<String, Object>();
                error.put("error", "Solver not found");
                error.put("solverId", solverId);
                return error;
            }

            return solver;
        }, json);
    }

    /**
     * Get detailed information for a specific solver.
     */
    private static Map<String, Object> getSolverInfo(String solverId) {
        Map<String, Object> solver = new HashMap<String, Object>();

        switch (solverId) {
            case "mva":
                solver.put("id", "mva");
                solver.put("name", "Mean Value Analysis");
                solver.put("description", "Analytical solver for queueing networks using Mean Value Analysis algorithm. " +
                        "Computes average performance metrics without enumerating the state space.");
                solver.put("type", "analytical");
                solver.put("formats", new String[]{"jsim", "jsimg", "jsimw", "lqnx", "xml"});
                solver.put("features", new String[]{"closed", "open", "mixed", "multiclass", "multiserver"});
                solver.put("async", false);
                solver.put("options", getMvaOptions());
                break;
            case "ctmc":
                solver.put("id", "ctmc");
                solver.put("name", "Continuous-Time Markov Chain");
                solver.put("description", "Exact state-space analysis by constructing the underlying CTMC. " +
                        "Provides exact results but may be limited by state space size.");
                solver.put("type", "analytical");
                solver.put("formats", new String[]{"jsim", "jsimg", "jsimw"});
                solver.put("features", new String[]{"closed", "open", "mixed", "multiclass", "exact"});
                solver.put("async", false);
                break;
            case "fluid":
                solver.put("id", "fluid");
                solver.put("name", "Fluid Solver");
                solver.put("description", "Approximates the queueing network using a deterministic fluid model " +
                        "based on ordinary differential equations.");
                solver.put("type", "analytical");
                solver.put("formats", new String[]{"jsim", "jsimg", "jsimw"});
                solver.put("features", new String[]{"closed", "open", "mixed", "multiclass", "transient"});
                solver.put("async", false);
                break;
            case "jmt":
                solver.put("id", "jmt");
                solver.put("name", "Java Modelling Tools");
                solver.put("description", "Discrete event simulation using the JMT simulation engine. " +
                        "Supports general service time distributions.");
                solver.put("type", "simulation");
                solver.put("formats", new String[]{"jsim", "jsimg", "jsimw"});
                solver.put("features", new String[]{"closed", "open", "mixed", "multiclass", "multiserver", "general_distributions"});
                solver.put("async", true);
                break;
            case "nc":
                solver.put("id", "nc");
                solver.put("name", "Normalizing Constant");
                solver.put("description", "Computes performance metrics using normalizing constant methods " +
                        "for product-form queueing networks.");
                solver.put("type", "analytical");
                solver.put("formats", new String[]{"jsim", "jsimg", "jsimw", "lqnx", "xml"});
                solver.put("features", new String[]{"closed", "multiclass", "product_form"});
                solver.put("async", false);
                break;
            case "ssa":
                solver.put("id", "ssa");
                solver.put("name", "Stochastic Simulation Algorithm");
                solver.put("description", "Event-driven stochastic simulation with exact sampling capabilities.");
                solver.put("type", "simulation");
                solver.put("formats", new String[]{"jsim", "jsimg", "jsimw"});
                solver.put("features", new String[]{"closed", "open", "mixed", "multiclass", "exact_sampling"});
                solver.put("async", true);
                break;
            case "des":
                solver.put("id", "des");
                solver.put("name", "Discrete Event Simulation");
                solver.put("description", "Native discrete event simulation using the SSJ library. " +
                        "Supports general distributions and open networks.");
                solver.put("type", "simulation");
                solver.put("formats", new String[]{"jsim", "jsimg", "jsimw"});
                solver.put("features", new String[]{"open", "multiclass", "multiserver", "general_distributions"});
                solver.put("async", true);
                break;
            case "ln":
                solver.put("id", "ln");
                solver.put("name", "Layered Network Solver");
                solver.put("description", "Analytical solver for layered queueing networks using fixed-point iteration.");
                solver.put("type", "analytical");
                solver.put("formats", new String[]{"lqnx", "xml"});
                solver.put("features", new String[]{"layered", "multiclass", "multiprocessor"});
                solver.put("async", false);
                break;
            case "lqns":
                solver.put("id", "lqns");
                solver.put("name", "LQN Solver");
                solver.put("description", "Solves layered queueing networks using the LQNS analytical solver.");
                solver.put("type", "analytical");
                solver.put("formats", new String[]{"lqnx", "xml"});
                solver.put("features", new String[]{"layered", "multiclass", "multiprocessor"});
                solver.put("async", false);
                break;
            default:
                return null;
        }

        return solver;
    }

    /**
     * Get MVA solver options.
     */
    private static List<Map<String, Object>> getMvaOptions() {
        List<Map<String, Object>> options = new ArrayList<Map<String, Object>>();

        Map<String, Object> method = new HashMap<String, Object>();
        method.put("name", "method");
        method.put("type", "string");
        method.put("description", "MVA variant to use");
        method.put("values", new String[]{"default", "exact", "amva", "qd"});
        method.put("default", "default");
        options.add(method);

        Map<String, Object> tolerance = new HashMap<String, Object>();
        tolerance.put("name", "tolerance");
        tolerance.put("type", "number");
        tolerance.put("description", "Convergence tolerance for iterative methods");
        tolerance.put("default", 1e-6);
        options.add(tolerance);

        Map<String, Object> maxIter = new HashMap<String, Object>();
        maxIter.put("name", "maxIterations");
        maxIter.put("type", "integer");
        maxIter.put("description", "Maximum number of iterations");
        maxIter.put("default", 1000);
        options.add(maxIter);

        return options;
    }
}
