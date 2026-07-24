/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.cli;

import jline.lang.Network;
import jline.lang.NetworkStruct;
import jline.solvers.ldes.SolverLDES;
import jline.solvers.ldes.LDESOptions;
import jline.solvers.ldes.LDESResult;
import jline.io.LineModelIO;
import jline.io.LDESResultIO;

/**
 * CLI for the standalone ldes.jar (LINE Discrete Event Simulator).
 *
 * <p>Usage: {@code java -jar ldes.jar solve model.json -o result.json [OPTIONS]}
 */
public class LdesCLI {

    /**
     * Prints help for the solve subcommand.
     */
    static void printHelp() {
        System.out.println("LDES - LINE Discrete Event Simulator");
        System.out.println();
        System.out.println("USAGE:");
        System.out.println("  java -jar ldes.jar solve <model.json> -o <result.json> [OPTIONS]");
        System.out.println();
        System.out.println("OPTIONS:");
        System.out.println("  -o <path>              Output result JSON file (required)");
        System.out.println("  -s, --samples <N>      Number of service completion events (default: 200000)");
        System.out.println("  -e, --maxevents <N>    Maximum number of simulation events (default: no limit)");
        System.out.println("      --maxtime <S>      Wall-clock time budget in seconds (default: no limit)");
        System.out.println("  --seed <N>             Random seed (default: auto)");
        System.out.println("  --method <M>           Solution method (default: default)");
        System.out.println("  --cnvgon               Enable convergence-based stopping");
        System.out.println("  --cnvgtol <V>          Convergence tolerance (default: 0.05)");
        System.out.println("  --tranfilter <F>       Transient filter: mser5, fixed, none (default: mser5)");
        System.out.println("  --warmupfrac <V>       Warmup fraction for fixed filter (default: 0.2)");
        System.out.println("  --cimethod <M>         CI method: obm, bm, spectral, none (default: obm)");
        System.out.println("  --mserbatch <N>        MSER batch size, used when tranfilter=mser5 (default: 5)");
        System.out.println("  --obmoverlap <V>       OBM overlap fraction, used when cimethod=obm (default: 0.5)");
        System.out.println("  --ciminbatch <N>       Minimum batch size for CI computation (default: 10)");
        System.out.println("  --ciminobs <N>         Minimum post-warmup observations for a CI (default: 100)");
        System.out.println("  --spectrallowfreqfrac <V>  Low-frequency fraction for spectral CI (default: 0.25)");
        System.out.println("  --cnvgbatch <N>        Minimum batches before the first convergence check (default: 20)");
        System.out.println("  --cnvgchk <N>          Events between convergence checks, 0=auto (default: 0)");
        System.out.println("  --replications <N>     Number of independent replications (default: 1)");
        System.out.println("  --numthreads <N>       Parallel threads for replications (default: auto)");
        System.out.println("  --timespan <T0,T1>     Transient analysis time span");
        System.out.println("  --initsol <v1,v2,...>  Initial placement, station-major [st0_cl0, st0_cl1, ...]");
        System.out.println("  --trajectory           Request trajectory data");
        System.out.println("  --respt-samples        Export per-job response time samples (empirical CDF input)");
        System.out.println("  --slotted              Discrete-time mode: all intervals must fall on the slot lattice");
        System.out.println("  --slotlength <V>       Slot length in model time units, implies --slotted (default: 1)");
        System.out.println("  -h, --help             Show this help");
        System.out.println();
        System.out.println("EXAMPLES:");
        System.out.println("  java -jar ldes.jar solve model.json -o result.json");
        System.out.println("  java -jar ldes.jar solve model.json -o result.json -e 500000 --seed 42");
        System.out.println("  java -jar ldes.jar solve model.json -o result.json --cnvgon --cnvgtol 0.01");
        System.out.println("  java -jar ldes.jar solve model.json -o result.json --timespan 0,100 --trajectory");
        System.out.println("  java -jar ldes.jar solve model.json -o result.json --replications 5 --numthreads 4");
    }

    /**
     * Handles the 'solve' command for LDES simulation.
     *
     * @param args command-line arguments after "solve"
     * @return 0 on success, non-zero on error
     */
    public static int handleSolveCommand(String[] args) {
        if (args.length == 0) {
            printHelp();
            return 1;
        }

        String modelPath = null;
        String outputPath = null;
        LDESOptions opts = new LDESOptions();
        int seed = -1;
        double[] timespan = null;
        boolean trajectory = false;
        boolean resptSamples = false;
        boolean hasSamples = false;
        boolean hasMaxEvents = false;

        for (int i = 0; i < args.length; i++) {
            switch (args[i]) {
                case "-h":
                case "--help":
                    printHelp();
                    return 0;
                case "-o":
                    if (i + 1 >= args.length) {
                        System.err.println("Error: -o requires a file path.");
                        return 1;
                    }
                    outputPath = args[++i];
                    break;
                case "-s":
                case "--samples":
                    if (i + 1 >= args.length) {
                        System.err.println("Error: -s/--samples requires a number.");
                        return 1;
                    }
                    try {
                        opts.samples = Integer.parseInt(args[++i]);
                        hasSamples = true;
                    } catch (NumberFormatException e) {
                        System.err.println("Error: Invalid samples value.");
                        return 1;
                    }
                    break;
                case "-e":
                case "--maxevents":
                    if (i + 1 >= args.length) {
                        System.err.println("Error: -e/--maxevents requires a number.");
                        return 1;
                    }
                    try {
                        opts.maxSimEvents = Integer.parseInt(args[++i]);
                        hasMaxEvents = true;
                    } catch (NumberFormatException e) {
                        System.err.println("Error: Invalid maxevents value.");
                        return 1;
                    }
                    break;
                case "--maxtime":
                    if (i + 1 >= args.length) {
                        System.err.println("Error: --maxtime requires a number (seconds).");
                        return 1;
                    }
                    try {
                        opts.maxTime = Double.parseDouble(args[++i]);
                    } catch (NumberFormatException e) {
                        System.err.println("Error: Invalid maxtime value.");
                        return 1;
                    }
                    break;
                case "--seed":
                    if (i + 1 >= args.length) {
                        System.err.println("Error: --seed requires a number.");
                        return 1;
                    }
                    try {
                        seed = Integer.parseInt(args[++i]);
                    } catch (NumberFormatException e) {
                        System.err.println("Error: Invalid seed value.");
                        return 1;
                    }
                    break;
                case "--method":
                    if (i + 1 >= args.length) {
                        System.err.println("Error: --method requires a value.");
                        return 1;
                    }
                    opts.method(args[++i]);
                    break;
                case "--cnvgon":
                    opts.cnvgon = true;
                    break;
                case "--slotted":
                    opts.slotted = true;
                    break;
                case "--slotlength":
                    if (i + 1 >= args.length) {
                        System.err.println("Error: --slotlength requires a value.");
                        return 1;
                    }
                    try {
                        opts.setSlotLength(Double.parseDouble(args[++i]));
                    } catch (NumberFormatException e) {
                        System.err.println("Error: Invalid slotlength value.");
                        return 1;
                    } catch (IllegalArgumentException e) {
                        System.err.println("Error: " + e.getMessage());
                        return 1;
                    }
                    break;
                case "--cnvgtol":
                    if (i + 1 >= args.length) {
                        System.err.println("Error: --cnvgtol requires a value.");
                        return 1;
                    }
                    try {
                        opts.cnvgtol = Double.parseDouble(args[++i]);
                    } catch (NumberFormatException e) {
                        System.err.println("Error: Invalid cnvgtol value.");
                        return 1;
                    }
                    break;
                case "--tranfilter":
                    if (i + 1 >= args.length) {
                        System.err.println("Error: --tranfilter requires a value.");
                        return 1;
                    }
                    opts.tranfilter = args[++i];
                    break;
                case "--warmupfrac":
                    if (i + 1 >= args.length) {
                        System.err.println("Error: --warmupfrac requires a value.");
                        return 1;
                    }
                    try {
                        opts.warmupfrac = Double.parseDouble(args[++i]);
                    } catch (NumberFormatException e) {
                        System.err.println("Error: Invalid warmupfrac value.");
                        return 1;
                    }
                    break;
                case "--cimethod":
                    if (i + 1 >= args.length) {
                        System.err.println("Error: --cimethod requires a value.");
                        return 1;
                    }
                    opts.cimethod = args[++i];
                    break;
                case "--mserbatch":
                    if (i + 1 >= args.length) {
                        System.err.println("Error: --mserbatch requires a number.");
                        return 1;
                    }
                    try {
                        opts.mserbatch = Integer.parseInt(args[++i]);
                    } catch (NumberFormatException e) {
                        System.err.println("Error: Invalid mserbatch value.");
                        return 1;
                    }
                    break;
                case "--obmoverlap":
                    if (i + 1 >= args.length) {
                        System.err.println("Error: --obmoverlap requires a value.");
                        return 1;
                    }
                    try {
                        opts.obmoverlap = Double.parseDouble(args[++i]);
                    } catch (NumberFormatException e) {
                        System.err.println("Error: Invalid obmoverlap value.");
                        return 1;
                    }
                    break;
                case "--ciminbatch":
                    if (i + 1 >= args.length) {
                        System.err.println("Error: --ciminbatch requires a number.");
                        return 1;
                    }
                    try {
                        opts.ciminbatch = Integer.parseInt(args[++i]);
                    } catch (NumberFormatException e) {
                        System.err.println("Error: Invalid ciminbatch value.");
                        return 1;
                    }
                    break;
                case "--ciminobs":
                    if (i + 1 >= args.length) {
                        System.err.println("Error: --ciminobs requires a number.");
                        return 1;
                    }
                    try {
                        opts.ciminobs = Integer.parseInt(args[++i]);
                    } catch (NumberFormatException e) {
                        System.err.println("Error: Invalid ciminobs value.");
                        return 1;
                    }
                    break;
                case "--spectrallowfreqfrac":
                    if (i + 1 >= args.length) {
                        System.err.println("Error: --spectrallowfreqfrac requires a value.");
                        return 1;
                    }
                    try {
                        opts.spectralLowFreqFrac = Double.parseDouble(args[++i]);
                    } catch (NumberFormatException e) {
                        System.err.println("Error: Invalid spectrallowfreqfrac value.");
                        return 1;
                    }
                    break;
                case "--cnvgbatch":
                    if (i + 1 >= args.length) {
                        System.err.println("Error: --cnvgbatch requires a number.");
                        return 1;
                    }
                    try {
                        opts.cnvgbatch = Integer.parseInt(args[++i]);
                    } catch (NumberFormatException e) {
                        System.err.println("Error: Invalid cnvgbatch value.");
                        return 1;
                    }
                    break;
                case "--cnvgchk":
                    if (i + 1 >= args.length) {
                        System.err.println("Error: --cnvgchk requires a number.");
                        return 1;
                    }
                    try {
                        opts.cnvgchk = Integer.parseInt(args[++i]);
                    } catch (NumberFormatException e) {
                        System.err.println("Error: Invalid cnvgchk value.");
                        return 1;
                    }
                    break;
                case "--replications":
                    if (i + 1 >= args.length) {
                        System.err.println("Error: --replications requires a number.");
                        return 1;
                    }
                    try {
                        opts.replications = Integer.parseInt(args[++i]);
                    } catch (NumberFormatException e) {
                        System.err.println("Error: Invalid replications value.");
                        return 1;
                    }
                    break;
                case "--numthreads":
                    if (i + 1 >= args.length) {
                        System.err.println("Error: --numthreads requires a number.");
                        return 1;
                    }
                    try {
                        opts.numThreads = Integer.parseInt(args[++i]);
                    } catch (NumberFormatException e) {
                        System.err.println("Error: Invalid numthreads value.");
                        return 1;
                    }
                    break;
                case "--timespan":
                    if (i + 1 >= args.length) {
                        System.err.println("Error: --timespan requires T0,T1 values.");
                        return 1;
                    }
                    try {
                        String[] parts = args[++i].split(",");
                        if (parts.length != 2) {
                            System.err.println("Error: --timespan requires exactly two comma-separated values (T0,T1).");
                            return 1;
                        }
                        timespan = new double[]{Double.parseDouble(parts[0]), Double.parseDouble(parts[1])};
                    } catch (NumberFormatException e) {
                        System.err.println("Error: Invalid timespan values.");
                        return 1;
                    }
                    break;
                case "--respt-samples":
                    resptSamples = true;
                    break;
                case "--trajectory":
                    trajectory = true;
                    break;
                case "--export-histogram":
                    opts.exportStateHistogram = true;
                    break;
                case "--initsol":
                    if (i + 1 >= args.length) {
                        System.err.println("Error: --initsol requires comma-separated values.");
                        return 1;
                    }
                    try {
                        String[] vals = args[++i].split(",");
                        jline.util.matrix.Matrix initSol = new jline.util.matrix.Matrix(1, vals.length);
                        for (int v = 0; v < vals.length; v++) {
                            initSol.set(0, v, Double.parseDouble(vals[v]));
                        }
                        opts.init_sol = initSol;
                    } catch (NumberFormatException e) {
                        System.err.println("Error: Invalid initsol values.");
                        return 1;
                    }
                    break;
                default:
                    if (args[i].startsWith("-")) {
                        System.err.println("Error: Unknown option: " + args[i]);
                        printHelp();
                        return 1;
                    }
                    if (modelPath == null) {
                        modelPath = args[i];
                    } else {
                        System.err.println("Error: Unexpected argument: " + args[i]);
                        return 1;
                    }
                    break;
            }
        }

        if (modelPath == null) {
            System.err.println("Error: Model file path is required.");
            printHelp();
            return 1;
        }
        if (outputPath == null) {
            System.err.println("Error: Output file path (-o) is required.");
            printHelp();
            return 1;
        }

        if (hasSamples && hasMaxEvents) {
            System.err.println("Warning: Both -s/--samples and -e/--maxevents specified. "
                    + "Simulation will stop when either limit is reached first.");
        }

        if (seed >= 0) {
            opts.seed(seed);
        }

        try {
            Object loaded = LineModelIO.load(modelPath);
            if (!(loaded instanceof Network)) {
                System.err.println("Error: LDES currently supports Network models only (not LayeredNetwork).");
                return 1;
            }
            Network model = (Network) loaded;

            // Set the transient horizon on `opts` BEFORE constructing the solver:
            // the SolverLDES constructor snapshots options, so a timespan assigned
            // afterwards would not reach solver.options. Ensemble averaging of the
            // transient is driven by --replications (opts.replications), parsed
            // above; replications > 1 routes through the parallel analyzer, which
            // now averages the per-bucket transient series across replications.
            if (timespan != null) {
                opts.timespan = timespan;
            }
            SolverLDES solver = new SolverLDES(model, opts);
            if (timespan != null) {
                solver.getTranAvg();
            } else {
                solver.getAvg();
            }

            LDESResult ldesResult = (LDESResult) solver.result;
            NetworkStruct sn = model.getStruct(false);
            LDESResultIO.save(ldesResult, sn, outputPath, trajectory, resptSamples);

            return 0;
        } catch (Exception e) {
            System.err.println("Error: " + e.getMessage());
            e.printStackTrace(System.err);
            return 1;
        }
    }

    /**
     * Main entry point for the LDES standalone JAR.
     * @param args Command line arguments
     */
    public static void main(String[] args) {
        // Accept both "solve model.json ..." and "model.json ..." (implicit solve)
        if (args.length > 0 && "solve".equals(args[0])) {
            String[] solveArgs = new String[args.length - 1];
            System.arraycopy(args, 1, solveArgs, 0, solveArgs.length);
            System.exit(handleSolveCommand(solveArgs));
        } else {
            System.exit(handleSolveCommand(args));
        }
    }
}
