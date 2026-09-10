/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_FLUID_FLUID_EXPORT_ODES_H
#define LINE_SOLVERS_FLUID_FLUID_EXPORT_ODES_H

/**
 * `@@SolverFLD/exportODEs.m`: the fluid ODE system as a standalone LaTeX
 * document, in a form meant to be read by a person AND parsed by a program.
 *
 * THE MACHINE-READABLE PART is the comment header, before the preamble: one
 * `% STATE s station=... class=... phase=...` line per state variable and, in
 * the J form, one `% EVENT e var=... type=... coeff=...` line per event. A tool
 * that wants the system and not the typesetting reads those and stops at
 * `\documentclass`. The comments are the contract; the body below them is
 * presentation, and the two are generated from the same `FluidSymSystem`.
 *
 * TWO NOTATIONS. `scalar` writes one expanded equation per state variable,
 * which is what one reads to understand a small model. `matrix` writes the
 * compact form -- dx/dt = W' theta(x) + lambda, or dx/dt = J r(x) -- with the
 * numeric matrices printed once, which is the only readable option once the
 * state space passes a few dozen entries.
 *
 * WHAT IS NOT REPRODUCED, and is called out in the exported Remarks exactly as
 * the reference calls it out: the regularization of vanishing denominators, the
 * FCFS non-exponential re-fitting loop, and any immediate-transition
 * elimination. The exported system is the nominal one.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <string>
#include <vector>

#include "line/lang/qn/network_struct.h"
#include "line/solvers/fluid/fluid_symodes.h"
#include "line/solvers/fluid/solver_fluid.h"
#include "line/util/error.h"

namespace line {
namespace fluid {

namespace detail {

/** `fmtnum`: integers plain, everything else at %.8g, infinity as `\infty`. */
inline std::string sym_fmtnum(double v) {
    if (std::isinf(v)) return v > 0 ? "\\infty" : "-\\infty";
    char buf[64];
    if (v == std::floor(v) && std::fabs(v) < 1e15) {
        std::snprintf(buf, sizeof(buf), "%lld", static_cast<long long>(v));
        return std::string(buf);
    }
    std::snprintf(buf, sizeof(buf), "%.8g", v);
    return std::string(buf);
}

/** `texesc`: escape the identifier characters LaTeX would otherwise eat. */
inline std::string sym_texesc(const std::string& s) {
    std::string out;
    for (char ch : s) {
        if (ch == '_' || ch == '%' || ch == '&' || ch == '#') out.push_back('\\');
        out.push_back(ch);
    }
    return out;
}

inline std::string sym_int(std::size_t v) {
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%llu", static_cast<unsigned long long>(v));
    return std::string(buf);
}

/** `factor_tex`: the state-dependent factor of a term driven by variable v. */
inline std::string sym_factor_tex(std::size_t v1, const SymFactor& f) {
    const std::string xv = "x_{" + sym_int(v1) + "}";
    if (f.type == "lin") return xv;
    if (f.type == "dpspw")
        return xv + "\\,g_{" + sym_int(f.station) + "," + sym_int(f.cls) + "}(\\mathbf{x})";
    if (f.type == "ext1") {
        if (f.others.empty()) return "";  // a single-phase source class is unit mass
        std::string s = "\\bigl(1 - ";
        for (std::size_t i = 0; i < f.others.size(); ++i) {
            if (i) s += " - ";
            s += "x_{" + sym_int(f.others[i]) + "}";
        }
        return s + "\\bigr)";
    }
    return xv + "\\,g_{" + sym_int(f.station) + "}(\\mathbf{x})";
}

/** `term_tex`: a term with positive coefficient c and factor fstr. */
inline std::string sym_term_tex(double c, const std::string& fstr) {
    if (fstr.empty()) return sym_fmtnum(c);
    if (c == 1.0) return fstr;
    return sym_fmtnum(c) + "\\," + fstr;
}

inline std::string sym_num_matrix(const Matrix<double>& A) {
    std::string body;
    for (std::size_t i = 0; i < A.rows(); ++i) {
        if (i) body += " \\\\ ";
        for (std::size_t j = 0; j < A.cols(); ++j) {
            if (j) body += " & ";
            body += sym_fmtnum(A(i, j));
        }
    }
    const std::size_t mx = std::max(A.rows(), A.cols());
    if (mx > 12) return "{\\scriptsize\\begin{bmatrix} " + body + " \\end{bmatrix}}";
    return "\\begin{bmatrix} " + body + " \\end{bmatrix}";
}

inline std::string sym_num_vector(const std::vector<double>& v) {
    std::string body;
    for (std::size_t i = 0; i < v.size(); ++i) {
        if (i) body += " & ";
        body += sym_fmtnum(v[i]);
    }
    return "\\begin{pmatrix} " + body + " \\end{pmatrix}";
}

/**
 * `build_terms`: T(s,v) is the constant coefficient of the term driven by
 * variable v in the equation of state s, `varFactor` the factor of v, and
 * `constTerm` the additive constant of state s.
 */
struct SymTerms {
    Matrix<double> T;
    std::vector<SymFactor> var_factor;
    std::vector<bool> has_factor;
    std::vector<double> const_term;
};
inline SymTerms sym_build_terms(const FluidSymSystem& sys) {
    const std::size_t n = sys.nstates;
    SymTerms t;
    t.T = Matrix<double>(n, n, 0.0);
    t.var_factor.assign(n, SymFactor());
    t.has_factor.assign(n, false);
    t.const_term.assign(n, 0.0);
    if (sys.form == "W") {
        for (std::size_t s = 0; s < n; ++s)
            for (std::size_t v = 0; v < n; ++v)
                t.T(s, v) = sys.is_source[v] ? 0.0 : sys.W(v, s);  // W', Source theta is zero
        for (std::size_t v = 0; v < n; ++v) {
            if (sys.is_source[v]) continue;
            const std::size_t i = sys.state_station[v];
            SymFactor f;
            f.station = i;
            f.cls = sys.state_class[v];
            f.type = std::isinf(sys.S[i - 1]) ? "lin" : sys.smoothing;
            t.var_factor[v] = f;
            t.has_factor[v] = true;
        }
        t.const_term = sys.alambda;
        return t;
    }
    for (std::size_t e = 0; e < sys.nevents; ++e) {
        const std::size_t v = sys.event_var[e];
        for (std::size_t s = 0; s < n; ++s) t.T(s, v) += sys.J(s, e) * sys.coeff[e];
        if (!t.has_factor[v]) {
            t.var_factor[v] = sys.factor[e];
            t.has_factor[v] = true;
        }
    }
    return t;
}

/** `build_defs`: the station-level auxiliary quantities the factors refer to. */
inline std::vector<std::string> sym_build_defs(const FluidSymSystem& sys, const SymTerms& t) {
    const std::size_t n = sys.nstates, M = sys.station_names.size(),
                      K = sys.class_names.size();
    std::vector<bool> need_n(M, false), need_nt(M, false), need_nh(M, false);
    std::vector<std::string> gdef(M);
    for (std::size_t v = 0; v < n; ++v) {
        if (!t.has_factor[v]) continue;
        const SymFactor& f = t.var_factor[v];
        const std::size_t i = f.station, i0 = i - 1;
        if (f.type == "lin") continue;
        if (f.type == "min") {
            need_n[i0] = true;
            gdef[i0] = "g_{" + sym_int(i) + "}(\\mathbf{x}) &= \\frac{\\min(n_{" + sym_int(i) +
                       "}(\\mathbf{x}),\\, " + sym_fmtnum(sys.S[i0]) + ")}{n_{" + sym_int(i) +
                       "}(\\mathbf{x})}";
        } else if (f.type == "pnorm") {
            need_n[i0] = true;
            const std::string ps = sym_fmtnum(i0 < sys.pstar.size() ? sys.pstar[i0] : 0.0);
            gdef[i0] = "g_{" + sym_int(i) + "}(\\mathbf{x}) &= \\Bigl(1 + \\bigl(n_{" +
                       sym_int(i) + "}(\\mathbf{x})/" + sym_fmtnum(sys.S[i0]) + "\\bigr)^{" + ps +
                       "}\\Bigr)^{-1/" + ps + "}";
        } else if (f.type == "dpsmin") {
            need_n[i0] = true;
            need_nt[i0] = true;
            gdef[i0] = "g_{" + sym_int(i) + "}(\\mathbf{x}) &= \\frac{\\min(n_{" + sym_int(i) +
                       "}(\\mathbf{x}),\\, " + sym_fmtnum(sys.S[i0]) + ")}{\\tilde{n}_{" +
                       sym_int(i) + "}(\\mathbf{x})}";
        } else if (f.type == "dpspw") {
            need_n[i0] = true;
            need_nt[i0] = true;
        } else if (f.type == "fcfsw" || f.type == "fcfsws") {
            need_n[i0] = true;
            need_nh[i0] = true;
            const std::string inner = (f.type == "fcfsw")
                                          ? "\\min(n_{" + sym_int(i) + "}(\\mathbf{x}),\\, " +
                                                sym_fmtnum(sys.S[i0]) + ")"
                                          : "\\mathrm{softmin}\\bigl(n_{" + sym_int(i) +
                                                "}(\\mathbf{x}),\\, " + sym_fmtnum(sys.S[i0]) +
                                                "\\bigr)";
            gdef[i0] = "g_{" + sym_int(i) + "}(\\mathbf{x}) &= \\frac{" + inner + "}{\\hat{n}_{" +
                       sym_int(i) + "}(\\mathbf{x})}";
        }
    }
    std::vector<std::string> defs;
    for (std::size_t i0 = 0; i0 < M; ++i0) {
        if (!need_n[i0]) continue;
        std::string s = "n_{" + sym_int(i0 + 1) + "}(\\mathbf{x}) &= ";
        bool first = true;
        for (std::size_t v = 0; v < n; ++v)
            if (sys.state_station[v] == i0 + 1) {
                if (!first) s += " + ";
                s += "x_{" + sym_int(v + 1) + "}";
                first = false;
            }
        defs.push_back(s + "\\\\");
    }
    for (std::size_t i0 = 0; i0 < M; ++i0) {
        if (!need_nt[i0]) continue;
        std::string s = "\\tilde{n}_{" + sym_int(i0 + 1) + "}(\\mathbf{x}) &= ";
        bool firstpart = true;
        for (std::size_t r = 0; r < K; ++r) {
            std::string inner;
            bool first = true;
            for (std::size_t v = 0; v < n; ++v)
                if (sys.state_station[v] == i0 + 1 && sys.state_class[v] == r + 1) {
                    if (!first) inner += " + ";
                    inner += "x_{" + sym_int(v + 1) + "}";
                    first = false;
                }
            if (inner.empty()) continue;
            if (!firstpart) s += " + ";
            s += sym_fmtnum(sys.dpsw(i0, r)) + "\\,(" + inner + ")";
            firstpart = false;
        }
        defs.push_back(s + "\\\\");
    }
    for (std::size_t i0 = 0; i0 < M; ++i0) {
        if (!need_nh[i0]) continue;
        std::string s = "\\hat{n}_{" + sym_int(i0 + 1) + "}(\\mathbf{x}) &= ";
        bool first = true;
        for (std::size_t v = 0; v < n; ++v)
            if (sys.state_station[v] == i0 + 1) {
                if (!first) s += " + ";
                s += sym_fmtnum(sys.fcfs_phase_w[v]) + "\\,x_{" + sym_int(v + 1) + "}";
                first = false;
            }
        defs.push_back(s + "\\\\");
    }
    for (std::size_t i0 = 0; i0 < M; ++i0)
        if (!gdef[i0].empty()) defs.push_back(gdef[i0] + "\\\\");
    for (std::size_t v = 0; v < n; ++v) {
        if (!t.has_factor[v] || t.var_factor[v].type != "dpspw") continue;
        const std::size_t i = t.var_factor[v].station, r = t.var_factor[v].cls;
        const std::string Si = sym_fmtnum(sys.S[i - 1]);
        const std::string d = "g_{" + sym_int(i) + "," + sym_int(r) +
                              "}(\\mathbf{x}) &= \\begin{cases} 1 & n_{" + sym_int(i) +
                              "}(\\mathbf{x}) \\le " + Si + "\\\\ \\dfrac{" +
                              sym_fmtnum(sys.S[i - 1] * sys.dpsw(i - 1, r - 1)) + "}{\\tilde{n}_{" +
                              sym_int(i) + "}(\\mathbf{x})} & n_{" + sym_int(i) +
                              "}(\\mathbf{x}) > " + Si + " \\end{cases}\\\\";
        if (std::find(defs.begin(), defs.end(), d) == defs.end()) defs.push_back(d);
    }
    if (sys.alpha > 0.0) {
        bool any_ws = false;
        for (std::size_t v = 0; v < n; ++v)
            if (t.has_factor[v] && t.var_factor[v].type == "fcfsws") any_ws = true;
        if (any_ws)
            defs.push_back(
                "\\mathrm{softmin}(a,b) &= \\frac{a\\,e^{-\\alpha a} + b\\,e^{-\\alpha "
                "b}}{e^{-\\alpha a} + e^{-\\alpha b}}, \\qquad \\alpha = " +
                sym_fmtnum(sys.alpha) + "\\\\");
    }
    if (!defs.empty()) {
        std::string& last = defs.back();
        if (last.size() >= 2 && last.compare(last.size() - 2, 2, "\\\\") == 0)
            last.erase(last.size() - 2);
    }
    return defs;
}

/** `render_equation`: one align row of the scalar notation. */
inline std::string sym_render_equation(const FluidSymSystem& sys, std::size_t s, const SymTerms& t,
                                       bool is_last) {
    std::vector<std::pair<double, std::string>> terms;
    for (std::size_t v = 0; v < sys.nstates; ++v) {
        const double c = t.T(s, v);
        if (c == 0.0) continue;
        terms.push_back(std::make_pair(c, sym_factor_tex(v + 1, t.var_factor[v])));
    }
    if (t.const_term[s] != 0.0)
        terms.push_back(std::make_pair(t.const_term[s], std::string()));
    std::string rhs;
    if (terms.empty()) {
        rhs = "0";
    } else {
        for (std::size_t k = 0; k < terms.size(); ++k) {
            const double c = terms[k].first;
            const std::string body = sym_term_tex(std::fabs(c), terms[k].second);
            if (k == 0)
                rhs += (c < 0.0) ? ("-" + body) : body;
            else
                rhs += (c < 0.0) ? (" - " + body) : (" + " + body);
            if ((k + 1) % 4 == 0 && k + 1 < terms.size()) rhs += "\\nonumber\\\\\n&\\quad ";
        }
    }
    return "\\frac{\\mathrm{d}x_{" + sym_int(s + 1) + "}}{\\mathrm{d}t} &= " + rhs +
           (is_last ? "" : "\\\\");
}

}  // namespace detail

/**
 * Render the fluid ODE system of `sn` as a LaTeX document.
 *
 * @param notation "scalar" (one expanded ODE per state) or "matrix"
 * @param model_name the name printed in the header, `model.getName` in the reference
 * @param sn the refreshed network struct
 * @param opt fluid options, which fix the method whose ODEs are exported
 */
template <class T>
std::string fluid_export_odes(const qn::NetworkStruct<T>& sn, const FluidOptions& opt,
                              const std::string& notation = "scalar",
                              const std::string& model_name = "model") {
    using namespace detail;
    if (notation != "scalar" && notation != "matrix")
        throw InputError("fluid_export_odes: unknown notation '" + notation +
                         "'; valid notations are scalar and matrix");

    std::string m = opt.method;
    if (m.compare(0, 6, "fluid.") == 0) m = m.substr(6);
    const bool use_pnorm = (m == "pnorm");
    std::vector<double> init = opt.init_sol;
    if (init.empty()) init = detail::fluid_default_initsol(sn, fluid_layout(sn));
    const FluidSymSystem sys =
        fluid_symodes(sn, opt.method, use_pnorm ? opt.pstar : 0.0, init);
    const std::size_t n = sys.nstates;
    const SymTerms terms = sym_build_terms(sys);

    std::vector<std::string> L;
    // ---- machine-readable header -----------------------------------------
    L.push_back("% Mean-field fluid ODE system exported by LINE SolverFLD");
    L.push_back("% model: " + model_name);
    L.push_back("% method: " + sys.method);
    L.push_back(sys.form == "W" ? "% form: dx/dt = W^T*theta(x) + lambda"
                                : "% form: dx/dt = J*r(x)");
    L.push_back("% notation: " + notation);
    L.push_back("% nstates: " + sym_int(n));
    if (sys.form == "J") L.push_back("% nevents: " + sym_int(sys.nevents));
    for (std::size_t s = 0; s < n; ++s)
        L.push_back("% STATE " + sym_int(s + 1) + " station=" +
                    sys.station_names[sys.state_station[s] - 1] + " class=" +
                    sys.class_names[sys.state_class[s] - 1] + " phase=" +
                    sym_int(sys.state_phase[s]));
    if (sys.form == "J")
        for (std::size_t e = 0; e < sys.nevents; ++e) {
            char buf[64];
            std::snprintf(buf, sizeof(buf), "%.15g", sys.coeff[e]);
            L.push_back("% EVENT " + sym_int(e + 1) + " var=" + sym_int(sys.event_var[e] + 1) +
                        " type=" + sys.factor[e].type + " coeff=" + std::string(buf));
        }

    // ---- preamble ---------------------------------------------------------
    L.push_back("\\documentclass{article}");
    L.push_back("\\usepackage{amsmath}");
    L.push_back("\\usepackage[margin=2.5cm]{geometry}");
    L.push_back("\\allowdisplaybreaks");
    L.push_back("\\setcounter{MaxMatrixCols}{500}");
    L.push_back("\\begin{document}");
    L.push_back("\\section*{Mean-field fluid ODE system}");
    L.push_back("\\noindent Model: \\texttt{" + sym_texesc(model_name) +
                "}. Solver: \\texttt{SolverFLD}, method \\texttt{" + sym_texesc(sys.method) +
                "}, " + notation + " notation.");
    if (sys.form == "W")
        L.push_back("The system has " + sym_int(n) +
                    " state variables and reads $\\frac{\\mathrm{d}\\mathbf{x}}{\\mathrm{d}t} = "
                    "W^{\\top}\\theta(\\mathbf{x}) + \\boldsymbol{\\lambda}$.");
    else
        L.push_back("The system has " + sym_int(n) + " state variables and " +
                    sym_int(sys.nevents) +
                    " events and reads $\\frac{\\mathrm{d}\\mathbf{x}}{\\mathrm{d}t} = "
                    "J\\,r(\\mathbf{x})$.");

    // ---- state legend -----------------------------------------------------
    L.push_back("\\subsection*{State variables}");
    L.push_back(
        "Each state variable $x_{s}$ is the mean number of jobs of a class in a service phase at "
        "a station:");
    L.push_back("\\begin{center}");
    const std::size_t chunk = 48;
    for (std::size_t s0 = 0; s0 < n; s0 += chunk) {
        const std::size_t s1 = std::min(n, s0 + chunk);
        L.push_back("\\begin{tabular}{rlll}");
        L.push_back("\\hline");
        L.push_back("$s$ & station & class & phase\\\\");
        L.push_back("\\hline");
        for (std::size_t s = s0; s < s1; ++s)
            L.push_back(sym_int(s + 1) + " & \\texttt{" +
                        sym_texesc(sys.station_names[sys.state_station[s] - 1]) + "} & \\texttt{" +
                        sym_texesc(sys.class_names[sys.state_class[s] - 1]) + "} & " +
                        sym_int(sys.state_phase[s]) + "\\\\");
        L.push_back("\\hline");
        L.push_back("\\end{tabular}");
        if (s1 < n) L.push_back("\\par\\medskip");
    }
    L.push_back("\\end{center}");

    std::vector<std::size_t> used;
    for (std::size_t s = 0; s < n; ++s)
        if (std::find(used.begin(), used.end(), sys.state_station[s]) == used.end())
            used.push_back(sys.state_station[s]);
    std::sort(used.begin(), used.end());
    L.push_back("\\begin{center}");
    L.push_back("\\begin{tabular}{rlll}");
    L.push_back("\\hline");
    L.push_back("$i$ & station & scheduling & $S_{i}$\\\\");
    L.push_back("\\hline");
    for (std::size_t i : used)
        L.push_back(sym_int(i) + " & \\texttt{" + sym_texesc(sys.station_names[i - 1]) + "} & " +
                    sym_texesc(sys.sched_names[i - 1]) + " & $" + sym_fmtnum(sys.S[i - 1]) +
                    "$\\\\");
    L.push_back("\\hline");
    L.push_back("\\end{tabular}");
    L.push_back("\\end{center}");

    // ---- definitions ------------------------------------------------------
    const std::vector<std::string> defs = sym_build_defs(sys, terms);
    if (!defs.empty()) {
        L.push_back("\\subsection*{Definitions}");
        L.push_back("\\begin{align*}");
        for (const std::string& d : defs) L.push_back(d);
        L.push_back("\\end{align*}");
    }

    // ---- the system itself ------------------------------------------------
    if (notation == "scalar") {
        L.push_back("\\subsection*{ODE system (scalar notation)}");
        L.push_back("\\begin{align}");
        for (std::size_t s = 0; s < n; ++s)
            L.push_back(sym_render_equation(sys, s, terms, s + 1 == n));
        L.push_back("\\end{align}");
    } else {
        L.push_back("\\subsection*{ODE system (matrix notation)}");
        if (sys.form == "W") {
            bool have_lambda = false;
            for (double v : sys.alambda)
                if (v != 0.0) have_lambda = true;
            L.push_back("\\begin{equation}");
            L.push_back(have_lambda
                            ? "\\frac{\\mathrm{d}\\mathbf{x}}{\\mathrm{d}t} = "
                              "W^{\\top}\\,\\theta(\\mathbf{x}) + \\boldsymbol{\\lambda}"
                            : "\\frac{\\mathrm{d}\\mathbf{x}}{\\mathrm{d}t} = "
                              "W^{\\top}\\,\\theta(\\mathbf{x})");
            L.push_back("\\end{equation}");
            L.push_back("with $\\theta_{s}(\\mathbf{x})$ given componentwise by");
            L.push_back("\\begin{equation*}");
            {
                std::string rows;
                for (std::size_t v = 0; v < n; ++v) {
                    if (v) rows += " \\\\ ";
                    rows += sys.is_source[v] ? std::string("0")
                                             : sym_factor_tex(v + 1, terms.var_factor[v]);
                }
                L.push_back("\\theta(\\mathbf{x}) = \\begin{bmatrix} " + rows +
                            " \\end{bmatrix}");
            }
            L.push_back("\\end{equation*}");
            L.push_back("and");
            L.push_back("\\begin{equation*}");
            {
                Matrix<double> Wt(n, n, 0.0);
                for (std::size_t i = 0; i < n; ++i)
                    for (std::size_t j = 0; j < n; ++j) Wt(i, j) = sys.W(j, i);
                L.push_back("W^{\\top} = " + sym_num_matrix(Wt));
            }
            L.push_back("\\end{equation*}");
            if (have_lambda) {
                L.push_back("\\begin{equation*}");
                L.push_back("\\boldsymbol{\\lambda} = " + sym_num_vector(sys.alambda) +
                            "^{\\top}");
                L.push_back("\\end{equation*}");
            }
        } else {
            L.push_back("\\begin{equation}");
            L.push_back("\\frac{\\mathrm{d}\\mathbf{x}}{\\mathrm{d}t} = J\\,r(\\mathbf{x})");
            L.push_back("\\end{equation}");
            L.push_back("with stoichiometry matrix");
            L.push_back("\\begin{equation*}");
            L.push_back("J = " + sym_num_matrix(sys.J));
            L.push_back("\\end{equation*}");
            L.push_back("and event rate functions");
            L.push_back("\\begin{align*}");
            for (std::size_t e = 0; e < sys.nevents; ++e)
                L.push_back("r_{" + sym_int(e + 1) + "}(\\mathbf{x}) &= " +
                            sym_term_tex(sys.coeff[e],
                                         sym_factor_tex(sys.event_var[e] + 1, sys.factor[e])) +
                            (e + 1 < sys.nevents ? "\\\\" : ""));
            L.push_back("\\end{align*}");
        }
    }

    // ---- initial condition and remarks ------------------------------------
    if (!sys.x0.empty()) {
        L.push_back("\\subsection*{Initial condition}");
        L.push_back("\\begin{equation*}");
        L.push_back("\\mathbf{x}(0) = " + sym_num_vector(sys.x0) + "^{\\top}");
        L.push_back("\\end{equation*}");
    }
    L.push_back("\\subsection*{Remarks}");
    L.push_back("\\begin{itemize}");
    L.push_back(
        "\\item For each station $i$, $n_{i}(\\mathbf{x})$ denotes the total mass at the station "
        "and $S_{i}$ the number of servers (infinite-server stations use the closed job "
        "population, $\\infty$ denotes infinity).");
    L.push_back(
        "\\item The numerical solver regularizes vanishing denominators with a small positive "
        "constant; these regularizations are omitted here.");
    if (sys.form == "J") {
        bool fcfs_fac = false, dps_fac = false;
        for (const SymFactor& f : sys.factor) {
            if (f.type == "fcfsw" || f.type == "fcfsws") fcfs_fac = true;
            if (f.type == "dpsmin") dps_fac = true;
        }
        if (fcfs_fac)
            L.push_back(
                "\\item At FCFS stations, the mean phase residence times $w_{u} = "
                "-1/[D_{0}]_{kk}$ weight the backlog $\\hat{n}_{i}$; the factors $w_{u}$ of the "
                "departing phases are folded into the rate coefficients.");
        if (dps_fac)
            L.push_back(
                "\\item At DPS stations, weights are normalized to sum to one and the weight "
                "$w_{ir}$ of the departing class is folded into the rate coefficient; the class "
                "shares $w_{ir}x/\\tilde{n}_{i}$ divide the station capacity $\\min(n_{i},S_{i})$, "
                "so they sum to one whenever the station is busy.");
    }
    {
        bool any_fcfs = false;
        for (std::size_t s = 0; s < n; ++s)
            if (sys.sched[sys.state_station[s] - 1] == lang::SchedStrategy::FCFS) any_fcfs = true;
        if (any_fcfs && (sys.method == "matrix" || sys.method == "closing"))
            L.push_back(
                "\\item For FCFS stations with non-exponential service, the solver may "
                "iteratively re-fit the service distributions (non-exponential approximation); "
                "the exported system uses the nominal model parameters.");
    }
    L.push_back("\\end{itemize}");
    L.push_back("\\end{document}");

    std::string tex;
    for (std::size_t i = 0; i < L.size(); ++i) {
        tex += L[i];
        tex += "\n";
    }
    return tex;
}

}  // namespace fluid
}  // namespace line

#endif  // LINE_SOLVERS_FLUID_FLUID_EXPORT_ODES_H
