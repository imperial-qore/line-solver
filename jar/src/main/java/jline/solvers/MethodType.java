package jline.solvers;

import java.util.HashMap;
import java.util.Map;

/**
 * Classification of a solution method, as printed in the solver banner:
 * <code>&lt;accuracy&gt;, &lt;randomness&gt;</code> with accuracy in
 * {exact, approximate, bound} and randomness in {deterministic, randomized}.
 *
 * <p>Conventions, applied uniformly across the four codebases (MATLAB
 * <code>line_method_type.m</code>, python <code>solvers/base.py</code>, C++
 * <code>method_type.h</code>):</p>
 * <ul>
 *   <li><b>exact</b>: the algorithm targets the metric with no modeling
 *       approximation. Numerical truncation and floating-point error do not
 *       make a method approximate, so an integral representation or a
 *       transform inversion is exact while an asymptotic expansion is not.</li>
 *   <li><b>approximate</b>: the algorithm introduces a heuristic, an
 *       asymptotic expansion, a decomposition, or a statistical estimate.</li>
 *   <li><b>bound</b>: the algorithm returns a formal one-sided bound on the
 *       metric, not a point estimate: the value is guaranteed to lie on the
 *       stated side of the exact one, and the two sides of a family bracket
 *       it. The side is read off the method label and printed with it
 *       ("gb.upper" -&gt; "upper bound").</li>
 *   <li><b>randomized</b>: the algorithm consumes pseudo-random numbers, so
 *       two runs agree only if the seed does.</li>
 * </ul>
 *
 * <p>Perfect sampling (cftp) is classified by the law it samples from, which is
 * the stationary one, hence exact; the ordinary simulators are approximate
 * because a finite horizon leaves warm-up bias on top of the sampling error.</p>
 *
 * <p>Keep this registry in step with the citation registry
 * (<code>line_citations.m</code>) and with the three twins named above.</p>
 */
public class MethodType {

    public static final String EXACT_DET = "exact, deterministic";
    public static final String EXACT_RND = "exact, randomized";
    public static final String APPROX_DET = "approximate, deterministic";
    public static final String APPROX_RND = "approximate, randomized";
    public static final String BOUND_DET = "bound, deterministic";

    private static final Map<String, String> REGISTRY = new HashMap<String, String>();

    private static void reg(String label, String... methodNames) {
        for (int i = 0; i < methodNames.length; i++) {
            REGISTRY.put(methodNames[i], label);
        }
    }

    static {
        // product-form evaluation and the exact single-queue closed forms
        reg(EXACT_DET, "exact", "mva", "mvac", "recal", "conv", "ca", "comom", "comomld",
                "rd", "nrp", "nrl", "nre", "clw", "gleint", "mmint2", "lcfsqn.ca", "rgf", "dnc", "ger",
                "divdiff",
                "nintmva", "nc.oi.exact", "sdr", "sdr.mva");
        reg(EXACT_DET, "mm1", "mmk", "mxm1", "mm1k", "mg1", "gm1", "mapm1ps", "pas");
        reg(EXACT_DET, "mg1.prio", "mg1.fb", "mg1.srpt", "mg1.psjf", "mg1.setf", "mg1.lrpt",
                "mm1.dps");
        // mixed open/closed limited load dependence: the Bruell-Balbo-Afshari
        // effective-capacity MVA evaluates the product form itself, no expansion
        reg(EXACT_DET, "ncldmx");
        // state-space enumeration: every CTMC path but the perfect samplers
        reg(EXACT_DET, "ctmc", "sync", "flat", "gpu", "fd", "uniformization");
        reg(EXACT_DET, "exact.mapmap1");
        reg(EXACT_DET, "jmva", "jmva.mva", "jmva.recal", "jmva.comom");
        reg(EXACT_DET, "lossn.exact");
        // MDD-rec: the normalising constant of a product form, summed EXACTLY
        // over the reachable set by one memoised walk of the decision diagram
        // that holds it. On a Petri net (Solver_nc_spn_analyzer) and on a loss
        // network (Lossn_rec) alike.
        reg(EXACT_DET, "rec", "lossn.rec", "mdd.rec");
        // discrete-time (slotted) product form: the Bernoulli server of
        // chapter 2 and the closed cycle of chapter 3 in Daduna (2001)
        reg(EXACT_DET, "dt.bernoulli1", "dt.cycle", "dt.cycleld");

        // coupling from the past samples the stationary law itself
        reg(EXACT_RND, "cftp", "ctmc.cftp");

        // approximate MVA and its variants
        reg(APPROX_DET, "amva", "bs", "aql", "qsa", "lin", "gflin", "egflin", "dmlin", "qd",
                "qdlin", "qdaql", "qli", "fli", "ab", "schmidt", "schmidt-ext", "schmidtext",
                "sqni", "sum", "esum", "cl", "chandy-lakshmi", "shadow", "seidmann",
                "linearizerms", "conway", "rolia", "zhou", "suri", "reiser.ms", "chow",
                "lcp", "pamb", "pami", "pamt", "clust",
                "marie", "sqd", "mapqn", "interp", "highvar", "balanced", "tay", "scat");
        // open-network decomposition and general-service closed forms
        reg(APPROX_DET, "qna", "rqna", "gig1", "gigk", "klb", "kraemer", "mg1k",
                "mm1k.approx");
        // asymptotic expansions, entropy and fixed-point methods
        reg(APPROX_DET, "le", "ble", "aghq", "cub", "kt", "bkt", "lekt", "pana", "panald", "mem", "mem.blocking",
                "gm", "erlangfp", "propfair", "fpi", "spm", "ttl");
        // balanced fairness aggregation is not closed under composition
        reg(APPROX_DET, "oi", "balancedfairness", "stationtime");
        // mean-field, diffusion and ODE methods
        reg(APPROX_DET, "fld", "fluid", "matrix", "closing", "statedep", "softmin", "pnorm",
                "mfq", "rmf", "tbi", "diffusion", "kp");
        // phase-type network decomposition
        reg(APPROX_DET, "mam", "dec", "mna", "inap", "inapplus", "inapinf", "ldqbd", "qbd",
                "qiu", "cdf", "reneging", "retrial");
        // layered and environment decomposition
        reg(APPROX_DET, "ln", "ln.mva", "layers", "ln.dec", "enhanced", "ln.fluid", "moment3", "lqns",
                "srvn", "lqnsdefault", "exactmva", "srvn.exactmva", "qns", "env", "env.blend",
                "blend", "dec.avg");
        reg(APPROX_DET, "jmva.amva", "jmva.chow", "jmva.bs", "jmva.aql", "jmva.lin",
                "jmva.dmlin");
        // solver selection is a meta-method; the selected solver's banner carries the
        // real classification
        reg(APPROX_DET, "auto", "tree", "auto.tree", "forest", "cart");

        // formal one-sided bounds; "cub.upper", "qrf.mem" and "qrf.bas.mem" are spelt
        // out because their tails "cub" and "mem" are NC method names
        reg(BOUND_DET, "ba", "aba", "bjb", "mbjb", "gb", "pb", "sb", "mwba", "pbh", "pbk",
                "bjbk", "cbh", "ssd", "sib", "scb", "ldbcmp", "qr", "lr", "qrf", "harel",
                "cub.upper", "qrf.mem", "qrf.bas.mem", "auto.upper", "auto.lower", "looping", "bpt", "bgt");

        // discrete-event simulation
        reg(APPROX_RND, "ssa", "ldes", "serial", "para", "parallel", "nrm", "jsim",
                "replication", "jmt", "lqsim", "sim", "uq");
        // Monte Carlo and sampling-based normalizing constants
        reg(APPROX_RND, "mci", "imci", "amci", "lhsmci", "ls", "is", "sampling", "lossn.mci");
        // Markov chain Monte Carlo on the regularized network (Chen-O'Cinneide)
        reg(APPROX_RND, "mcmc", "nc.mcmc");
        // the approximate sampler stops before coalescence
        reg(APPROX_RND, "cftp.approx");

        // per-solver default for a method with no entry of its own; '#' keeps the
        // solver name out of the method namespace, since 'mva' is also a method
        reg(EXACT_DET, "#ctmc");
        reg(APPROX_RND, "#ssa", "#ldes", "#jmt");
        reg(BOUND_DET, "#ba");
        reg(APPROX_DET, "#mva", "#nc", "#fld", "#mam", "#ln", "#env", "#lqns",
                "#qns", "#auto", "#uq");
    }

    /**
     * Banner classification of a solution method.
     *
     * Lookup order: <code>solver.method</code>, <code>method</code>, the tail after each
     * dot of the method (longest suffix first), the head before its first dot,
     * <code>#solver</code>, then "approximate, deterministic". The unknown-method default
     * is the conservative one: claiming exactness a method does not have is the costlier
     * error.
     *
     * @param solvername banner solver name, with or without the "Solver" prefix
     * @param method     resolved method label, possibly carrying the "default/" prefix
     * @return the classification string, e.g. "exact, deterministic"
     */
    public static String of(String solvername, String method) {
        String solver = solvername == null ? "" : solvername.trim().toLowerCase();
        if (solver.startsWith("solver")) {
            solver = solver.substring("solver".length());
        }
        String m = method == null ? "" : method.trim().toLowerCase();
        // the banner label is 'default/<resolved>' once a default has been resolved
        int slash = m.lastIndexOf('/');
        if (slash >= 0) {
            m = m.substring(slash + 1);
        }
        return boundSide(lookup(solver, m), m);
    }

    /**
     * A bound is reported with the side it lies on, which the method label carries as
     * its last component ("gb.upper" -&gt; "upper bound"). A family with no sided variant
     * (the qrf reductions) stays the unqualified "bound".
     */
    private static String boundSide(String label, String method) {
        if (!label.startsWith("bound")) {
            return label;
        }
        if (method.endsWith(".upper")) {
            return "upper " + label;
        }
        if (method.endsWith(".lower")) {
            return "lower " + label;
        }
        return label;
    }

    private static String lookup(String solver, String m) {
        if (!solver.isEmpty() && !m.isEmpty()) {
            String hit = REGISTRY.get(solver + "." + m);
            if (hit != null) {
                return hit;
            }
        }
        if (!m.isEmpty()) {
            String hit = REGISTRY.get(m);
            if (hit != null) {
                return hit;
            }
            int dot = m.indexOf('.');
            int from = dot;
            while (from >= 0) {     // 'a.b.c' -> 'b.c', 'c'
                hit = REGISTRY.get(m.substring(from + 1));
                if (hit != null) {
                    return hit;
                }
                from = m.indexOf('.', from + 1);
            }
            if (dot >= 0) {
                hit = REGISTRY.get(m.substring(0, dot));
                if (hit != null) {
                    return hit;
                }
            }
        }
        if (!solver.isEmpty()) {
            String hit = REGISTRY.get("#" + solver);
            if (hit != null) {
                return hit;
            }
        }
        return APPROX_DET;
    }
}
