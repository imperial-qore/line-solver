package jline.solvers.auto;

import jline.lang.Network;
import jline.solvers.NetworkSolver;
import jline.solvers.Solver;
import jline.solvers.SolverOptions;
import jline.solvers.ctmc.SolverCTMC;
import jline.solvers.fluid.SolverFluid;
import jline.solvers.jmt.SolverJMT;
import jline.solvers.mam.SolverMAM;
import jline.solvers.mva.SolverMVA;
import jline.solvers.nc.SolverNC;
import jline.solvers.ssa.SolverSSA;

/**
 * LINE solver - Java/Kotlin implementation equivalent to MATLAB LINE.m
 * <p>
 * This class serves as an alias for SolverAUTO and provides a static factory method
 * for creating solver instances based on the chosen method, matching the MATLAB interface.
 * <p>
 * Copyright (c) 2012-2026, Imperial College London
 * All rights reserved.
 */
public class LINE extends SolverAUTO {

    /**
     * Constructor - alias for SolverAUTO
     */
    public LINE(Network model, Object... varargin) {
        super(model, varargin);
    }

    /**
     * Constructor with options
     */
    public LINE(Network model, SolverOptions options) {
        super(model, options);
    }

    /**
     * Constructor with method string
     */
    public LINE(Network model, String method) {
        super(model, method);
    }

    /**
     * Convenience method - create LINE solver with default settings
     */
    public static LINE create(Network model) {
        return new LINE(model);
    }

    /**
     * Convenience method - create LINE solver with method
     */
    public static LINE create(Network model, String method) {
        return new LINE(model, method);
    }

    /**
     * Static factory method to create and configure a solver based on the chosen method
     * This matches the MATLAB LINE.load() static method functionality
     *
     * @param chosenMethod The solver method to use
     * @param model        The network model
     * @param varargin     Variable arguments for solver options
     * @return Configured solver instance
     */
    public static NetworkSolver load(String chosenMethod, Network model, Object... varargin) {
        // Parse options - use default options if none provided
        SolverOptions options = Solver.parseOptions(Solver.defaultOptions(), varargin);
        options.method = chosenMethod;

        // Create solver based on method
        switch (options.method.toLowerCase()) {
            // AUTO/DEFAULT methods
            case "default":
            case "auto":
                if ("auto".equals(options.method)) {
                    options.method = "default";
                }
                return new LINE(model, options);

            // CTMC methods
            case "ctmc":
            case "ctmc.gpu":
            case "gpu":
                if ("ctmc".equals(options.method)) {
                    options.method = "default";
                }
                options.method = options.method.replace("ctmc.", "");
                return new SolverCTMC(model, options);

            // MVA methods
            case "mva":
            case "mva.exact":
            case "amva":
            case "mva.amva":
            case "qna":
            case "mva.qna":
            case "sqrt":
            case "mva.sqrt":
            case "amva.qd":
            case "mva.amva.qd":
            case "amva.bs":
            case "mva.amva.bs":
            case "amva.qli":
            case "mva.amva.qli":
            case "amva.fli":
            case "mva.amva.fli":
            case "amva.lin":
            case "mva.amva.lin":
            case "mm1":
            case "mmk":
            case "mg1":
            case "gm1":
            case "gig1":
            case "gim1":
            case "gig1.kingman":
            case "gigk":
            case "gigk.kingman_approx":
            case "gig1.gelenbe":
            case "gig1.heyman":
            case "gig1.kimura":
            case "gig1.allen":
            case "gig1.kobayashi":
            case "gig1.klb":
            case "gig1.marchal":
            case "aba.upper":
            case "aba.lower":
            case "bjb.upper":
            case "bjb.lower":
            case "gb.upper":
            case "gb.lower":
            case "pb.upper":
            case "pb.lower":
            case "sb.upper":
            case "sb.lower":
                if ("mva".equals(options.method)) {
                    options.method = "default";
                }
                options.method = options.method.replace("mva.", "");
                return new SolverMVA(model, options);

            // SSA methods
            case "nrm":
            case "ssa":
            case "ssa.serial":
            case "ssa.parallel":
            case "serial":
            case "parallel":
                if ("ssa".equals(options.method)) {
                    options.method = "default";
                }
                options.method = options.method.replace("ssa.", "");
                return new SolverSSA(model, options);

            // JMT methods
            case "jmt":
            case "jsim":
            case "jmva":
            case "jmva.mva":
            case "jmva.recal":
            case "jmva.comom":
            case "jmva.chow":
            case "jmva.bs":
            case "jmva.aql":
            case "jmva.lin":
            case "jmva.dmlin":
            case "jmt.jsim":
            case "jmt.jmva":
            case "jmt.jmva.mva":
            case "jmt.jmva.amva":
            case "jmva.amva":
            case "jmt.jmva.recal":
            case "jmt.jmva.comom":
            case "jmt.jmva.chow":
            case "jmt.jmva.bs":
            case "jmt.jmva.aql":
            case "jmt.jmva.lin":
            case "jmt.jmva.dmlin":
                if ("jmt".equals(options.method)) {
                    options.method = "default";
                }
                options.method = options.method.replace("jmt.", "");
                return new SolverJMT(model, options);

            // Fluid methods
            case "fluid":
            case "fluid.softmin":
            case "fluid.statedep":
            case "fluid.closing":
                if ("fluid".equals(options.method)) {
                    options.method = "default";
                }
                options.method = options.method.replace("fluid.", "");
                return new SolverFluid(model, options);

            // NC methods
            case "nc":
            case "nc.exact":
            case "nc.imci":
            case "nc.ls":
            case "comomrm":
            case "comomld":
            case "cub":
            case "ls":
            case "nc.le":
            case "le":
            case "mmint2":
            case "nc.panacea":
            case "nc.pana":
            case "nc.mmint2":
            case "nc.kt":
            case "nc.deterministic":
            case "nc.sampling":
            case "nc.propfair":
            case "nc.comom":
            case "nc.comomld":
            case "nc.mom":
            case "nc.cub":
            case "nc.brute":
            case "nc.rd":
            case "nc.nr.probit":
            case "nc.nr.logit":
            case "nc.gm":
                if ("nc".equals(options.method)) {
                    options.method = "default";
                }
                options.method = options.method.replace("nc.", "");
                return new SolverNC(model, options);

            // MAM methods
            case "mam":
            case "mam.dec.source":
            case "mam.dec.mmap":
            case "mam.dec.poisson":
                if ("mam".equals(options.method)) {
                    options.method = "default";
                }
                options.method = options.method.replace("mam.", "");
                return new SolverMAM(model, options);

            // Default fallback
            default:
                if ("auto".equals(options.method)) {
                    options.method = "default";
                }
                return new LINE(model, options);
        }
    }

    /**
     * Loads a Network model from a file.
     * Detects the file format and loads the appropriate model type.
     *
     * @param filename The path to the file to load
     * @return The loaded Network or LayeredNetwork model
     * @throws RuntimeException if the file format is unsupported or loading fails
     */
    public static Network load(String filename) {
        return load(filename, false);
    }

    /**
     * Loads a Network model from a file with optional verbose output.
     * Detects the file format and loads the appropriate model type.
     *
     * @param filename The path to the file to load
     * @param verbose Whether to print verbose loading information
     * @return The loaded Network or LayeredNetwork model
     * @throws RuntimeException if the file format is unsupported or loading fails
     */
    public static Network load(String filename, boolean verbose) {
        // Get file extension
        String lowerFilename = filename.toLowerCase();

        if (lowerFilename.endsWith(".xml") || lowerFilename.endsWith(".lqn") || lowerFilename.endsWith(".lqnx")) {
            // Load as LayeredNetwork (LQN model)
            try {
                jline.lang.layered.LayeredNetwork lqnModel = jline.lang.layered.LayeredNetwork.load(filename, verbose);
                if (verbose) {
                    System.out.println("Loaded LQN model from: " + filename);
                }
                return lqnModel.getModel(0);
            } catch (Exception e) {
                throw new RuntimeException("Failed to load LQN model from: " + filename, e);
            }
        } else if (lowerFilename.endsWith(".jsim") || lowerFilename.endsWith(".jsimg") || lowerFilename.endsWith(".jsimw")) {
            // JSIM/JMT format - would need JSIM2LINE implementation
            throw new UnsupportedOperationException("JSIM/JMT format loading not yet implemented in Java. Use MATLAB LINE.load() for JSIM files.");
        } else if (lowerFilename.endsWith(".jmva")) {
            // JMVA format - would need JMVA2LINE implementation
            throw new UnsupportedOperationException("JMVA format loading not yet implemented in Java. Use MATLAB LINE.load() for JMVA files.");
        } else {
            throw new IllegalArgumentException("Unsupported file format: " + filename +
                ". Supported formats: .xml, .lqn, .lqnx (LQN models)");
        }
    }
}