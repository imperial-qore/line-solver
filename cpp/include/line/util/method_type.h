#ifndef LINE_UTIL_METHOD_TYPE_H
#define LINE_UTIL_METHOD_TYPE_H

/**
 * Classification of a solution method, as printed in the solver banner:
 * "<accuracy>, <randomness>" with accuracy in {exact, approximate, bound} and
 * randomness in {deterministic, randomized}.
 *
 * Conventions, applied uniformly across the four codebases (MATLAB
 * line_method_type.m, JAR MethodType.java, python solvers/base.py):
 *   - exact:       the algorithm targets the metric with no modeling
 *                  approximation. Numerical truncation and floating-point
 *                  error do not make a method approximate, so an integral
 *                  representation or a transform inversion is exact while an
 *                  asymptotic expansion is not.
 *   - approximate: the algorithm introduces a heuristic, an asymptotic
 *                  expansion, a decomposition, or a statistical estimate.
 *   - bound:       the algorithm returns a formal one-sided bound on the
 *                  metric, not a point estimate: the value is guaranteed to
 *                  lie on the stated side of the exact one, and the two sides
 *                  of a family bracket it. The side is read off the method
 *                  label and printed with it ('gb.upper' -> 'upper bound').
 *   - randomized:  the algorithm consumes pseudo-random numbers, so two runs
 *                  agree only if the seed does.
 *
 * Perfect sampling (cftp) is classified by the law it samples from, which is
 * the stationary one, hence exact; the ordinary simulators are approximate
 * because a finite horizon leaves warm-up bias on top of the sampling error.
 *
 * The C++ CLI keeps its own key=value banner grammar, so this classification is
 * printed there as `type=exact,deterministic` rather than inside the bracketed
 * field list the other three codebases use.
 */

#include <map>
#include <string>

namespace line {
namespace util {

namespace detail {

inline const std::map<std::string, std::string>& method_type_registry() {
    static const std::map<std::string, std::string> reg = [] {
        std::map<std::string, std::string> r;
        struct Group {
            const char* label;
            const char* tokens;   // space separated
        };
        // product-form evaluation and the exact single-queue closed forms
        const Group groups[] = {
            {"exact,deterministic",
             "exact mva mvac recal conv ca comom comomld rd nrp nrl nre clw gleint mmint2 rgf "
             "ger "
             "lcfsqn.ca nc.oi.exact mm1 mmk mxm1 mm1k mg1 gm1 mapm1ps pas mg1.prio mg1.fb "
             "mg1.srpt mg1.psjf mg1.setf mg1.lrpt mm1.dps ctmc sync flat gpu fd "
             "uniformization exact.mapmap1 jmva jmva.mva jmva.recal jmva.comom lossn.exact "
             // Krzesinski state-dependent routing: eq. (16) by enumeration, and
             // its Section 4 MVA and convolution
             "sdr sdr.mva "
             // MDD-rec: the normalising constant of a product form, summed
             // EXACTLY over the reachable set by one memoised walk of the
             // decision diagram that holds it. On a Petri net
             // (solver_nc_spn_analyzer) and on a loss network (lossn_rec) alike
             "rec lossn.rec mdd.rec "
             // discrete-time (slotted) product form, Daduna (2001)
             "dt.bernoulli1 dt.cycle dt.cycleld "
             // mixed open/closed limited load dependence: the Bruell-Balbo-Afshari
             // effective-capacity MVA evaluates the product form itself, no expansion
             "ncldmx"},
            // coupling from the past samples the stationary law itself
            {"exact,randomized", "cftp ctmc.cftp"},
            // approximate MVA, open-network decomposition, expansions, mean-field,
            // phase-type decomposition, layered and environment decomposition
            {"approximate,deterministic",
             "amva bs aql qsa lin gflin egflin dmlin qd qdlin qdaql qli fli ab schmidt "
             "schmidt-ext schmidtext sqni sum esum cl chandy-lakshmi shadow seidmann "
             "linearizerms conway rolia zhou suri reiser.ms chow marie sqd mapqn interp highvar "
             "balanced qna rqna rqt gig1 gigk klb kraemer mg1k mm1k.approx le ble dir aghq cub kt bkt lekt pana "
             "panald mem mem.blocking gm erlangfp propfair fpi spm ttl oi "
             "balancedfairness stationtime fld fluid matrix closing statedep softmin pnorm "
             "mfq rmf tbi diffusion minnormal refined kp "
             "mam dec mna inap inapplus inapinf ldqbd qbd qiu cdf "
             "reneging retrial ln layers ln.dec enhanced ln.fluid moment3 lqns srvn "
             "lqnsdefault exactmva srvn.exactmva ln.mva qns env env.blend blend dec.avg jmva.amva "
             "jmva.chow jmva.bs jmva.aql jmva.lin jmva.dmlin auto tree auto.tree forest cart"},
            // formal one-sided bounds; 'cub.upper', 'qrf.mem' and 'qrf.bas.mem' are
            // spelt out because their tails 'cub' and 'mem' are NC method names
            {"bound,deterministic",
             "ba aba bjb mbjb gb pb sb mwba mwba.upper mwba.lower pbh pbk bjbk cbh ssd sib scb "
             "ldbcmp qr lr qrf harel bpt bgt cub.upper qrf.mem qrf.bas.mem spnlp"},
            // discrete-event simulation, Monte Carlo normalizing constants, and the
            // sampler that stops before coalescence
            {"approximate,randomized",
             "ssa ldes serial para parallel nrm jsim replication jmt lqsim sim uq mci imci "
             "ls is sampling lossn.mci cftp.approx "
             // Markov chain Monte Carlo on the regularized network (Chen-O'Cinneide)
             "mcmc nc.mcmc"},
            // per-solver default for a method with no entry of its own; '#' keeps the
            // solver name out of the method namespace, since 'mva' is also a method
            {"exact,deterministic", "#ctmc"},
            {"approximate,randomized", "#ssa #ldes #jmt"},
            {"bound,deterministic", "#ba"},
            {"approximate,deterministic",
             "#mva #nc #fld #fluid #mam #ln #env #lqns #qns #auto #uq"},
        };
        for (std::size_t g = 0; g < sizeof(groups) / sizeof(groups[0]); ++g) {
            const std::string toks(groups[g].tokens);
            std::size_t i = 0;
            while (i < toks.size()) {
                const std::size_t j = toks.find(' ', i);
                const std::string tok = toks.substr(i, j == std::string::npos ? j : j - i);
                if (!tok.empty()) r[tok] = groups[g].label;
                if (j == std::string::npos) break;
                i = j + 1;
            }
        }
        return r;
    }();
    return reg;
}

inline std::string lower_trim(const std::string& s) {
    std::size_t b = s.find_first_not_of(" \t");
    if (b == std::string::npos) return std::string();
    std::size_t e = s.find_last_not_of(" \t");
    std::string out = s.substr(b, e - b + 1);
    for (std::size_t i = 0; i < out.size(); ++i) {
        if (out[i] >= 'A' && out[i] <= 'Z') out[i] = char(out[i] - 'A' + 'a');
    }
    return out;
}

// A bound is reported with the side it lies on, which the method label carries
// as its last component ('gb.upper' -> 'upper bound'). A family with no sided
// variant (the qrf reductions) stays the unqualified 'bound'.
inline std::string bound_side(const std::string& label, const std::string& method) {
    if (label.compare(0, 5, "bound") != 0) return label;
    const std::size_t n = method.size();
    if (n >= 6 && method.compare(n - 6, 6, ".upper") == 0) return "upper " + label;
    if (n >= 6 && method.compare(n - 6, 6, ".lower") == 0) return "lower " + label;
    return label;
}

inline std::string method_type_lookup(const std::string& solvername, const std::string& method) {
    const std::map<std::string, std::string>& reg = detail::method_type_registry();
    std::string solver = detail::lower_trim(solvername);
    if (solver.compare(0, 6, "solver") == 0) solver = solver.substr(6);
    std::string m = detail::lower_trim(method);
    // the banner label is 'default/<resolved>' once a default has been resolved
    const std::size_t slash = m.find_last_of('/');
    if (slash != std::string::npos) m = m.substr(slash + 1);

    std::map<std::string, std::string>::const_iterator it;
    if (!solver.empty() && !m.empty()) {
        it = reg.find(solver + "." + m);
        if (it != reg.end()) return it->second;
    }
    if (!m.empty()) {
        it = reg.find(m);
        if (it != reg.end()) return it->second;
        const std::size_t dot = m.find('.');
        std::size_t from = dot;
        while (from != std::string::npos) {   // 'a.b.c' -> 'b.c', 'c'
            it = reg.find(m.substr(from + 1));
            if (it != reg.end()) return it->second;
            from = m.find('.', from + 1);
        }
        if (dot != std::string::npos) {
            it = reg.find(m.substr(0, dot));
            if (it != reg.end()) return it->second;
        }
    }
    if (!solver.empty()) {
        it = reg.find("#" + solver);
        if (it != reg.end()) return it->second;
    }
    return "approximate,deterministic";
}

}  // namespace detail

/**
 * Banner classification of a solution method.
 *
 * Lookup order: "<solver>.<method>", "<method>", the tail after each dot of the
 * method (longest suffix first), the head before its first dot, "#<solver>",
 * then "approximate,deterministic". The unknown-method default is the
 * conservative one: claiming exactness a method does not have is the costlier
 * error.
 *
 * @param solvername banner solver name, with or without the "Solver" prefix
 * @param method     resolved method label, possibly carrying a "default/" prefix
 */
inline std::string method_type(const std::string& solvername, const std::string& method) {
    std::string m = detail::lower_trim(method);
    const std::size_t slash0 = m.find_last_of('/');
    if (slash0 != std::string::npos) m = m.substr(slash0 + 1);
    return detail::bound_side(detail::method_type_lookup(solvername, method), m);
}

}  // namespace util
}  // namespace line

#endif  // LINE_UTIL_METHOD_TYPE_H
