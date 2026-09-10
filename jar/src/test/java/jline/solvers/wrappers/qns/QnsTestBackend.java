package jline.solvers.wrappers.qns;

import jline.VerboseLevel;
import jline.solvers.SolverOptions;
import jline.solvers.wrappers.qns.analyzers.Solver_qns_analyzer;

/**
 * Decides, once per test run, whether the {@link SolverQNS} tests can reach
 * qnsolver, and configures each solver accordingly.
 *
 * <p>The only backend is a {@code qnsolver} binary on the PATH. qnsolver ships
 * with LQNS, whose licence is an evaluation agreement forbidding
 * redistribution, so neither LINE nor its test suite pulls or runs a container
 * image holding it. To test against a containerised build, run
 * {@code run-tests.sh --lqns-docker <image>}: that puts {@code lqns},
 * {@code lqsim} and {@code qnsolver} shims on the PATH which exec the
 * container, so everything below still sees a plain local binary and needs no
 * Docker knowledge of its own.
 *
 * <p>When no binary is found the QNS tests are skipped, not failed.
 */
final class QnsTestBackend {
    private QnsTestBackend() {}

    private static boolean initialized;
    private static boolean enabled;

    /** Whether the QNS tests can run at all (used by {@code @EnabledIf}). */
    static synchronized boolean isEnabled() {
        init();
        return enabled;
    }

    /** Fresh QNS options for a test. */
    static synchronized SolverOptions options() {
        init();
        SolverOptions options = SolverQNS.defaultOptions();
        options.verbose = VerboseLevel.SILENT;
        return options;
    }

    private static void init() {
        if (initialized) return;
        initialized = true;
        enabled = Solver_qns_analyzer.hasNativeQNSolver();
        if (!enabled) {
            System.out.println("[LINE-test] No qnsolver on the PATH; QNS tests skipped. "
                    + "Install LQNS, or use run-tests.sh --lqns-docker <image>.");
        }
    }
}
