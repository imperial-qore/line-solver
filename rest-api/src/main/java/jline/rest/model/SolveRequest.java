/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.rest.model;

/**
 * Request payload for the /models/solve endpoint.
 */
public class SolveRequest {

    /**
     * The model to solve.
     */
    private ModelInput model;

    /**
     * The solver to use.
     * Supported values: mva, ctmc, fluid, jmt, nc, ssa, ln, lqns, des
     */
    private String solver;

    /**
     * The type of analysis to perform.
     * Supported values: all, avg, sys
     */
    private String analysis;

    /**
     * Optional solver-specific options.
     */
    private SolverOptionsDTO options;

    /**
     * Default constructor for JSON deserialization.
     */
    public SolveRequest() {
    }

    public ModelInput getModel() {
        return model;
    }

    public void setModel(ModelInput model) {
        this.model = model;
    }

    public String getSolver() {
        return solver;
    }

    public void setSolver(String solver) {
        this.solver = solver;
    }

    public String getAnalysis() {
        return analysis;
    }

    public void setAnalysis(String analysis) {
        this.analysis = analysis;
    }

    public SolverOptionsDTO getOptions() {
        return options;
    }

    public void setOptions(SolverOptionsDTO options) {
        this.options = options;
    }

    /**
     * Validate the solve request.
     * @return null if valid, error message if invalid
     */
    public String validate() {
        if (model == null) {
            return "Model is required";
        }
        String modelError = model.validate();
        if (modelError != null) {
            return modelError;
        }

        // Validate solver
        if (solver == null || solver.trim().isEmpty()) {
            return "Solver is required";
        }
        String[] validSolvers = {"mva", "ctmc", "fluid", "jmt", "nc", "ssa", "ln", "lqns", "des"};
        boolean validSolver = false;
        for (String valid : validSolvers) {
            if (valid.equals(solver)) {
                validSolver = true;
                break;
            }
        }
        if (!validSolver) {
            return "Invalid solver: " + solver + ". Valid solvers: mva, ctmc, fluid, jmt, nc, ssa, ln, lqns, des";
        }

        // Validate analysis type (optional, defaults to "all")
        if (analysis != null && !analysis.trim().isEmpty()) {
            String[] validAnalysis = {"all", "avg", "sys"};
            boolean valid = false;
            for (String va : validAnalysis) {
                if (va.equals(analysis)) {
                    valid = true;
                    break;
                }
            }
            if (!valid) {
                return "Invalid analysis type: " + analysis + ". Valid types: all, avg, sys";
            }
        }

        // Validate solver compatibility with format
        String format = model.getFormat();
        if ("lqnx".equals(format) || "xml".equals(format)) {
            String[] lqnSolvers = {"ln", "lqns", "mva", "nc"};
            boolean compatible = false;
            for (String s : lqnSolvers) {
                if (s.equals(solver)) {
                    compatible = true;
                    break;
                }
            }
            if (!compatible) {
                return "Solver '" + solver + "' is not compatible with format '" + format + "'. " +
                       "Valid solvers for LQN/XML: ln, lqns, mva, nc";
            }
        } else {
            String[] jsimSolvers = {"mva", "ctmc", "fluid", "jmt", "nc", "ssa", "des"};
            boolean compatible = false;
            for (String s : jsimSolvers) {
                if (s.equals(solver)) {
                    compatible = true;
                    break;
                }
            }
            if (!compatible) {
                return "Solver '" + solver + "' is not compatible with format '" + format + "'. " +
                       "Valid solvers for JSIM: mva, ctmc, fluid, jmt, nc, ssa, des";
            }
        }

        return null;
    }
}
