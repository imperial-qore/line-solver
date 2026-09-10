/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_IO_JSIM_READER_H
#define LINE_IO_JSIM_READER_H

/**
 * Read a JMT `.jsim` / `.jsimg` / `.jsimw` model into a `qn::Network`.
 *
 * Port of `matlab/src/io/JSIM2LINE.m` (and `jline.io.M2M.JSIM2LINE`), the
 * inverse of `jmt_writer.h`. It is what lets `line-cli -i jsimg` solve a model
 * a user drew in JMT's GUI, which until now only the MATLAB and JAR front ends
 * could do.
 *
 * THE MODEL IS UNDER `<sim>` AND NOTHING ELSE IS READ. A `.jsimg` also carries
 * a sibling `<jmodel>` block holding GUI data -- station coordinates, class
 * colours -- and it repeats every `<userClass name=...>`. Reading the document
 * root instead of `<sim>` therefore doubles every class, silently, on exactly
 * the files a user is most likely to import.
 *
 * THE SECTION ORDER IS THE INTERFACE. A JMT node is three `<section>` elements
 * -- input, service, output -- and every parameter inside one is found by
 * POSITION, not by name: `saveGetStrategy` writes the put strategy fourth
 * because the reader counts to four. This reader therefore addresses parameters
 * positionally where the reference does, and by `name` attribute where a name
 * exists, which is what makes it survive the two JMT XML dialects (the
 * retrial-carrying 1.2.0 form shifts the put strategy from parameter 4 to 5,
 * and the reference tells them apart by reading parameter 3's name).
 *
 * WHAT IS NOT READ IS REFUSED BY NAME rather than dropped. A `.jsimg` can
 * describe blocking regions, cache sections and measures this reader does not
 * reconstruct; importing such a file and solving the remainder would answer a
 * question about a different model. The refusals name the section.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <limits>
#include <fstream>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "line/lang/lang_types.h"
#include "line/lang/qn/network_builder.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"
#include "line/util/xml.h"

#include <unistd.h>

namespace line {
namespace io {

namespace jsim_detail {

using xml::Element;

inline double to_num(const std::string& s) { return std::atof(s.c_str()); }

/** The text of an element, trimmed; empty when the element is absent. */
inline std::string text_of(const Element* e) {
    if (!e) return std::string();
    std::string t = e->text;
    std::size_t a = 0, b = t.size();
    while (a < b && (t[a] == ' ' || t[a] == '\n' || t[a] == '\t' || t[a] == '\r')) ++a;
    while (b > a && (t[b - 1] == ' ' || t[b - 1] == '\n' || t[b - 1] == '\t' || t[b - 1] == '\r'))
        --b;
    return t.substr(a, b - a);
}

/**
 * `<value>` of a parameter, which JMT writes as a child element rather than as
 * text on the parameter itself.
 */
inline std::string value_of(const Element* e) {
    if (!e) return std::string();
    const std::vector<const Element*> v = e->child_tags("value");
    if (!v.empty()) return text_of(v[0]);
    return text_of(e);
}

/** Direct `<parameter>` children, in document order. */
inline std::vector<const Element*> params(const Element* e) {
    return e ? e->child_tags("parameter") : std::vector<const Element*>();
}

/** Direct `<subParameter>` children, in document order. */
inline std::vector<const Element*> subs(const Element* e) {
    return e ? e->child_tags("subParameter") : std::vector<const Element*>();
}

/** The first direct child `<parameter>`/`<subParameter>` carrying `name`. */
inline const Element* named(const std::vector<const Element*>& v, const std::string& name) {
    for (std::size_t i = 0; i < v.size(); ++i)
        if (v[i]->attr("name") == name) return v[i];
    return nullptr;
}

/** Element i of a list, or null -- positional access that cannot read past the end. */
inline const Element* at(const std::vector<const Element*>& v, std::size_t i) {
    return i < v.size() ? v[i] : nullptr;
}

/** `\` and `/` become `_`, as the reference renames a node on import. */
inline std::string sanitize(const std::string& s) {
    std::string out = s;
    for (std::size_t i = 0; i < out.size(); ++i)
        if (out[i] == '/' || out[i] == '\\') out[i] = '_';
    return out;
}

/**
 * A JMT distribution block -- the `<subParameter name="...">` pair naming the
 * law and its parameters -- lowered to a `Distrib`.
 *
 * THE INVERSE OF `jmt_append_distribution`, and of the JMT GUI's own writer,
 * which are the same format. `dist` is the element whose `name` attribute is
 * the display name ("Exponential", "Burst (MAP)", ...) and `par` the
 * `distrPar` block beside it.
 *
 * A LAW THIS READER DOES NOT KNOW IS REFUSED BY NAME. The reference falls back
 * to matching two moments with an APH, which silently replaces the user's
 * distribution with a different one of the same mean and SCV; that is a
 * reasonable default in an interactive session and a wrong answer in a CLI, so
 * the name is reported instead.
 */
template <class T>
lang::Distrib<T> read_distribution(const Element* dist, const Element* par,
                                   const std::string& who) {
    if (!dist) return lang::Distrib<T>::disabled_dist();
    const std::string name = dist->attr("name");
    const std::string cp = dist->attr("classPath");
    if (name == "DisabledServiceTimeStrategy" || cp.find("DisabledServiceTime") != std::string::npos)
        return lang::Distrib<T>::disabled_dist();
    if (name == "ZeroServiceTimeStrategy" || cp.find("ZeroServiceTime") != std::string::npos)
        return lang::Distrib<T>::immediate();

    const std::vector<const Element*> p = subs(par);
    auto num = [&](std::size_t i) -> T {
        return num_traits<T>::from_double(to_num(value_of(at(p, i))));
    };
    auto bykey = [&](const char* key, std::size_t fallback) -> T {
        const Element* e = named(p, key);
        return e ? num_traits<T>::from_double(to_num(value_of(e))) : num(fallback);
    };

    if (name == "Exponential") return lang::Distrib<T>::exp_rate(bykey("lambda", 0));
    if (name == "Deterministic") return lang::Distrib<T>::det(bykey("t", 0));
    if (name == "Erlang") {
        // JMT's `alpha` is the PHASE rate and `r` the phase count, which is
        // exactly `Distrib::erlang`'s pair -- no conversion, and in particular
        // NOT the mean.
        const T a = bykey("alpha", 0);
        const Element* re = named(p, "r");
        const double r = re ? to_num(value_of(re)) : to_num(value_of(at(p, 1)));
        if (!(r >= 1.0))
            throw InputError(who + ": an Erlang needs at least one phase, got r=" +
                             std::to_string(r));
        return lang::Distrib<T>::erlang(a, static_cast<std::size_t>(r));
    }
    if (name == "Hyperexponential")
        return lang::Distrib<T>::hyperexp(bykey("p", 0), bykey("lambda1", 1), bykey("lambda2", 2));
    if (name == "Coxian") {
        // `Coxian([mu1 mu2], [phi1 1])` of the reference: the third parameter is
        // the completion probability of phase 1, and phase 2 always completes.
        std::vector<T> mu, phi;
        mu.push_back(bykey("lambda0", 0));
        mu.push_back(bykey("lambda1", 1));
        phi.push_back(bykey("phi0", 2));
        phi.push_back(num_traits<T>::from_int(1));
        return lang::Distrib<T>::coxian(mu, phi);
    }
    if (name == "Pareto") return lang::Distrib<T>::pareto(bykey("alpha", 0), bykey("k", 1));
    if (name == "Gamma") return lang::Distrib<T>::gamma_dist(bykey("alpha", 0), bykey("beta", 1));
    if (name == "Uniform") return lang::Distrib<T>::uniform(bykey("min", 0), bykey("max", 1));
    if (name == "Weibull") {
        // The constructor takes (scale, shape) and JMT writes alpha=scale,
        // r=shape; the reference's own comment records that the two are
        // inverted relative to the JMT parameter order.
        return lang::Distrib<T>::weibull(bykey("alpha", 0), bykey("r", 1));
    }
    if (name == "Lognormal")
        return lang::Distrib<T>::lognormal(bykey("mu", 0), bykey("sigma", 1));
    if (name == "Replayer" || name == "Trace") {
        const Element* fe = named(p, "fileName");
        const std::string path = fe ? value_of(fe) : value_of(at(p, 0));
        if (path.empty())
            throw InputError(who + ": a Replayer names no trace file");
        std::vector<T> samples;
        {
            std::ifstream tr(path.c_str());
            double v = 0.0;
            while (tr >> v) samples.push_back(num_traits<T>::from_double(v));
        }
        if (samples.empty())
            throw InputError(who + ": the Replayer trace '" + path +
                             "' is missing or empty; the model cannot be reconstructed without "
                             "the samples it replays");
        lang::Distrib<T> d = lang::Distrib<T>::replayer(samples);
        d.trace_file = path;
        return d;
    }
    if (name == "Burst (MMPP2)") {
        const double l0 = to_num(value_of(at(p, 0))), l1 = to_num(value_of(at(p, 1)));
        const double s0 = to_num(value_of(at(p, 2))), s1 = to_num(value_of(at(p, 3)));
        Matrix<T> D0(2, 2, num_traits<T>::from_int(0)), D1(2, 2, num_traits<T>::from_int(0));
        D1(0, 0) = num_traits<T>::from_double(l0);
        D1(1, 1) = num_traits<T>::from_double(l1);
        D0(0, 0) = num_traits<T>::from_double(-(l0 + s0));
        D0(0, 1) = num_traits<T>::from_double(s0);
        D0(1, 0) = num_traits<T>::from_double(s1);
        D0(1, 1) = num_traits<T>::from_double(-(l1 + s1));
        return lang::Distrib<T>::map_dist(D0, D1, lang::ProcessType::MMPP2);
    }
    // The two matrix forms: a square block of rows, each row an array of
    // `<value>` scalars. `jmt_object_array` writes the same nesting for both.
    auto read_matrix = [&](const Element* blk) {
        std::vector<std::vector<double> > rows;
        const std::vector<const Element*> rr = subs(blk);
        for (std::size_t i = 0; i < rr.size(); ++i) {
            std::vector<double> row;
            const std::vector<const Element*> cc = subs(rr[i]);
            for (std::size_t j = 0; j < cc.size(); ++j) row.push_back(to_num(value_of(cc[j])));
            if (!row.empty()) rows.push_back(row);
        }
        Matrix<T> M(rows.size(), rows.empty() ? 0 : rows[0].size(), num_traits<T>::from_int(0));
        for (std::size_t i = 0; i < rows.size(); ++i)
            for (std::size_t j = 0; j < rows[i].size(); ++j)
                M(i, j) = num_traits<T>::from_double(rows[i][j]);
        return M;
    };
    if (name == "Burst (MAP)") {
        const Element* d0 = named(p, "D0") ? named(p, "D0") : at(p, 0);
        const Element* d1 = named(p, "D1") ? named(p, "D1") : at(p, 1);
        return lang::Distrib<T>::map_dist(read_matrix(d0), read_matrix(d1),
                                          lang::ProcessType::MAP);
    }
    if (name == "Phase-Type") {
        const Element* av = named(p, "alpha") ? named(p, "alpha") : at(p, 0);
        const Element* Tb = named(p, "T") ? named(p, "T") : at(p, 1);
        // `alpha` is an array wrapping a `vector` array of scalars, one level
        // deeper than `T`'s rows.
        std::vector<T> alpha;
        const std::vector<const Element*> a1 = subs(av);
        const std::vector<const Element*> a2 = a1.empty() ? a1 : subs(a1[0]);
        for (std::size_t i = 0; i < a2.size(); ++i)
            alpha.push_back(num_traits<T>::from_double(to_num(value_of(a2[i]))));
        const Matrix<T> A = read_matrix(Tb);
        // APH when the sub-diagonal is empty, general PH otherwise -- the
        // reference's own test, and the flag decides which representation the
        // struct records rather than which numbers it holds.
        bool acyclic = true;
        for (std::size_t i = 0; i < A.rows(); ++i)
            for (std::size_t j = 0; j < i; ++j)
                if (num_traits<T>::to_double(A(i, j)) > 0.0) acyclic = false;
        return lang::Distrib<T>::phase_type(alpha, A, acyclic);
    }
    throw UnsupportedError(who + ": JMT distribution '" + name +
                           "' has no counterpart in this reader. The MATLAB importer matches its "
                           "first two moments with an APH instead; that substitutes a different "
                           "law for the user's, which a solve would then report as the model's "
                           "answer");
}

/** The (dist, par) pair inside a strategy wrapper, as the writer emits it. */
template <class T>
lang::Distrib<T> read_strategy_dist(const Element* wrapper, const std::string& who) {
    if (!wrapper) return lang::Distrib<T>::disabled_dist();
    const std::string nm = wrapper->attr("name");
    if (nm == "DisabledServiceTimeStrategy") return lang::Distrib<T>::disabled_dist();
    if (nm == "ZeroServiceTimeStrategy") return lang::Distrib<T>::immediate();
    const std::vector<const Element*> ss = subs(wrapper);
    if (ss.empty()) return lang::Distrib<T>::disabled_dist();
    return read_distribution<T>(at(ss, 0), at(ss, 1), who);
}

}  // namespace jsim_detail

/**
 * Read a JSIM document into a Network.
 *
 * @param path  the `.jsim` / `.jsimg` / `.jsimw` file
 * @param name  the model name; the document's own when empty
 */
template <class T>
qn::Network<T> read_jsim(const std::string& path, const std::string& name = std::string()) {
    using namespace jsim_detail;
    std::unique_ptr<Element> doc = xml::parse_file(path);
    if (!doc) throw InputError("read_jsim: cannot parse '" + path + "' as XML");

    // `<sim>` is the model; a `.jsimg` wraps it beside the GUI's `<jmodel>`.
    const Element* sim = doc->name == "sim" ? doc.get() : nullptr;
    if (!sim) {
        const std::vector<const Element*> s = doc->by_tag("sim");
        if (s.empty())
            throw InputError("read_jsim: '" + path +
                             "' carries no <sim> element, so it is not a JSIM model");
        sim = s[0];
    }

    qn::Network<T> net(name.empty() ? sim->attr("name") : name);
    const std::vector<const Element*> xnodes = sim->child_tags("node");
    const std::vector<const Element*> xclasses = sim->child_tags("userClass");
    if (xnodes.empty()) throw InputError("read_jsim: the model declares no node");
    if (xclasses.empty()) throw InputError("read_jsim: the model declares no job class");
    const std::size_t N = xnodes.size(), K = xclasses.size();

    std::vector<std::string> orig_name(N), node_name(N);
    for (std::size_t i = 0; i < N; ++i) {
        orig_name[i] = xnodes[i]->attr("name");
        node_name[i] = sanitize(orig_name[i]);
    }
    // The three sections of each node, and the class names of each section.
    std::vector<std::vector<const Element*> > sec(N);
    for (std::size_t i = 0; i < N; ++i) sec[i] = xnodes[i]->child_tags("section");

    auto secclass = [&](std::size_t i, std::size_t k) -> std::string {
        const Element* e = at(sec[i], k);
        return e ? e->attr("className") : std::string();
    };

    // ---- pass 1: create the nodes -----------------------------------------
    std::vector<std::size_t> nidx(N, 0);            ///< 1-based node index in the Network
    std::vector<lang::SchedStrategy> sched(N, lang::SchedStrategy::FCFS);
    std::vector<bool> is_station(N, false), is_source(N, false), is_delay(N, false);
    std::vector<std::vector<T> > schedparam(N);
    std::size_t forknode = 0;

    for (std::size_t i = 0; i < N; ++i) {
        const std::string in = secclass(i, 0), svc = secclass(i, 1), out = secclass(i, 2);
        if (in == "JobSink") {
            nidx[i] = net.add_sink(node_name[i]);
            continue;
        }
        if (in == "RandomSource") {
            nidx[i] = net.add_source(node_name[i]);
            is_station[i] = is_source[i] = true;
            continue;
        }
        if (in == "Join") {
            if (!forknode)
                throw UnsupportedError(
                    "read_jsim: a Join appears before any Fork; the importer supports at most one "
                    "fork-join pair and reads them in document order, as the reference does");
            nidx[i] = net.add_join(node_name[i], forknode);
            continue;
        }
        if (in == "Storage") {
            nidx[i] = net.add_place(node_name[i]);
            is_station[i] = true;
            continue;
        }
        if (in == "Enabling") {
            // Built in pass 3, once the mode structure has been read: the
            // constructor takes the whole TransitionParam and there is no
            // setter to fill it afterwards.
            continue;
        }
        if (in != "Queue" && in != "Buffer")
            throw UnsupportedError("read_jsim: node '" + orig_name[i] +
                                   "' has input section '" + in +
                                   "', which this reader does not build");

        if (out == "Fork") {
            const std::vector<const Element*> op = params(at(sec[i], 2));
            const double tpl = op.empty() ? 1.0 : to_num(value_of(op[0]));
            nidx[i] = net.add_fork(node_name[i], tpl);
            forknode = nidx[i];
            continue;
        }
        if (svc == "ServiceTunnel") {
            nidx[i] = net.add_router(node_name[i]);
            continue;
        }
        if (svc == "ClassSwitch" || svc == "StatelessClassSwitcher") {
            // The node is created HERE, in file order, and its matrix installed
            // in pass 3 once the classes exist: the matrix is indexed by class
            // and a `.jsimg` lists its nodes first. Creating the node later
            // instead would renumber every node after it, which the reference's
            // own commented-out attempt to "create the cs elements last" was
            // uncertain about for exactly that reason.
            nidx[i] = net.add_class_switch(node_name[i]);
            continue;
        }

        // ---- a queueing station: read the put strategy for the discipline --
        const std::vector<const Element*> ip = params(at(sec[i], 0));
        // The 1.2.0 dialect inserts `retrialDistributions` at position 3, which
        // pushes the put strategy from parameter 4 to parameter 5. The
        // reference tells the two apart by reading parameter 3's NAME, and so
        // does this: counting positions without the test reads the retrial
        // block as a scheduling strategy.
        const std::size_t putpos =
            (ip.size() > 2 && ip[2]->attr("name") == "retrialDistributions") ? 4 : 3;
        const std::vector<const Element*> put = subs(at(ip, putpos));
        const std::string putname = put.empty() ? std::string() : put[0]->attr("name");
        if (putname == "TailStrategy") sched[i] = lang::SchedStrategy::FCFS;
        else if (putname == "TailStrategyPriority") sched[i] = lang::SchedStrategy::HOL;
        else if (putname == "HeadStrategy") sched[i] = lang::SchedStrategy::LCFS;
        else if (putname == "RandStrategy") sched[i] = lang::SchedStrategy::SIRO;
        else if (putname == "SJFStrategy") sched[i] = lang::SchedStrategy::SJF;
        else if (putname == "SEPTStrategy") sched[i] = lang::SchedStrategy::SEPT;
        else if (putname == "LJFStrategy") sched[i] = lang::SchedStrategy::LJF;
        else if (putname == "LEPTStrategy") sched[i] = lang::SchedStrategy::LEPT;

        if (svc == "Delay" || svc == "InfiniteServer") {
            nidx[i] = net.add_delay(node_name[i]);
            is_station[i] = is_delay[i] = true;
        } else if (svc == "PSServer" || svc == "SharedServer") {
            // The PS family names its discipline in the PREEMPTIVE strategy
            // block of the service section, not in the buffer's put strategy.
            const std::vector<const Element*> sp = params(at(sec[i], 1));
            std::string ps = "EPSStrategy";
            if (sp.size() > 3) {
                const std::vector<const Element*> ss = subs(sp[3]);
                if (!ss.empty()) ps = ss[0]->attr("name");
            }
            if (ps == "EPSStrategy") sched[i] = lang::SchedStrategy::PS;
            else if (ps == "DPSStrategy") sched[i] = lang::SchedStrategy::DPS;
            else if (ps == "GPSStrategy") sched[i] = lang::SchedStrategy::GPS;
            else if (ps == "EPSStrategyPriority") sched[i] = lang::SchedStrategy::PSPRIO;
            else if (ps == "DPSStrategyPriority") sched[i] = lang::SchedStrategy::DPSPRIO;
            else if (ps == "GPSStrategyPriority") sched[i] = lang::SchedStrategy::GPSPRIO;
            // The per-class weights sit in the parameter after it.
            if (sp.size() > 4)
                for (const Element* w : subs(sp[4]))
                    schedparam[i].push_back(num_traits<T>::from_double(to_num(value_of(w))));
            nidx[i] = net.add_queue(node_name[i], sched[i]);
            is_station[i] = true;
        } else if (svc == "Server" || svc == "PreemptiveServer") {
            nidx[i] = net.add_queue(node_name[i], sched[i]);
            is_station[i] = true;
        } else {
            throw UnsupportedError("read_jsim: node '" + orig_name[i] + "' has service section '" +
                                   svc + "', which this reader does not build");
        }

        // Buffer size and, for a Server, the number of servers.
        if (!ip.empty()) {
            const double cap = to_num(value_of(ip[0]));
            // JMT writes -1 for "unbounded", which is Inf here and NOT a
            // capacity of -1: a negative buffer would refuse every arrival.
            net.set_capacity(nidx[i], cap < 0 ? std::numeric_limits<double>::infinity() : cap);
        }
        if (svc != "Delay" && svc != "InfiniteServer") {
            const std::vector<const Element*> sp = params(at(sec[i], 1));
            const Element* ns = named(sp, "maxJobs");
            if (!ns && !sp.empty()) ns = sp[0];
            if (ns) {
                const double c = to_num(value_of(ns));
                net.set_number_of_servers(
                    nidx[i], c < 0 ? std::numeric_limits<double>::infinity() : c);
            }
        }
    }

    // ---- pass 2: the job classes ------------------------------------------
    // JMT reads a HIGHER priority value as MORE important and LINE reads a
    // LOWER one as more important, so the values are reflected about the
    // maximum rather than copied. Copying them would invert every priority
    // ordering in the imported model, which no metric would flag as wrong.
    int maxprio = 0;
    for (std::size_t r = 0; r < K; ++r)
        maxprio = std::max(maxprio, static_cast<int>(to_num(xclasses[r]->attr("priority"))));

    std::map<std::string, std::size_t> node_by_name;
    for (std::size_t i = 0; i < N; ++i)
        if (nidx[i]) {
            node_by_name[node_name[i]] = nidx[i];
            node_by_name[orig_name[i]] = nidx[i];
        }

    std::vector<std::size_t> cidx(K, 0);
    std::vector<bool> cs_referenced(K, false);
    for (std::size_t r = 0; r < K; ++r) {
        const std::string cname = xclasses[r]->attr("name");
        const std::string type = xclasses[r]->attr("type");
        const int prio = maxprio - static_cast<int>(to_num(xclasses[r]->attr("priority")));
        const std::string refsrc = xclasses[r]->attr("referenceSource");
        if (type == "closed") {
            const std::map<std::string, std::size_t>::const_iterator it =
                node_by_name.find(refsrc);
            if (it == node_by_name.end())
                throw InputError("read_jsim: class '" + cname + "' names reference source '" +
                                 refsrc + "', which is not a node of the model");
            cidx[r] = net.add_closed_class(cname, to_num(xclasses[r]->attr("customers")),
                                           it->second, prio);
        } else {
            cidx[r] = net.add_open_class(cname, prio);
            // 'ClassSwitch' / 'StatelessClassSwitcher' is JMT's marker for an
            // open class that enters by class switching only; its Source
            // arrival is DISABLED rather than absent, which is a different
            // model from one whose arrival was simply not written.
            if (refsrc == "ClassSwitch" || refsrc == "StatelessClassSwitcher")
                cs_referenced[r] = true;
        }
    }

    // ---- pass 3: what needs the classes: ClassSwitch matrices, Petri nets --
    for (std::size_t i = 0; i < N; ++i) {
        const std::string svc3 = secclass(i, 1);
        if (nidx[i] && (svc3 == "ClassSwitch" || svc3 == "StatelessClassSwitcher")) {
            const std::vector<const Element*> sp = params(at(sec[i], 1));
            Matrix<T> C(K, K, num_traits<T>::from_int(0));
            if (!sp.empty()) {
                // `<parameter name="matrix">` holds one `<subParameter name="row">`
                // per class, each holding one `<subParameter name="cell">` per
                // class. TWO levels, not three: descending once more reads the
                // FIRST row's cells as the whole matrix, which leaves every
                // other row zero -- and a zero row means "this class switches
                // into nothing", so the classes it fed become unreachable and
                // the model is rejected several frames later for a reason that
                // names neither the matrix nor this node.
                const std::vector<const Element*> rows = subs(sp[0]);
                for (std::size_t r = 0; r < rows.size() && r < K; ++r) {
                    const std::vector<const Element*> cols = subs(rows[r]);
                    for (std::size_t c = 0; c < cols.size() && c < K; ++c)
                        C(r, c) = num_traits<T>::from_double(to_num(value_of(cols[c])));
                }
            }
            net.set_class_switch_matrix(nidx[i], C);
        }
        const std::string in = secclass(i, 0);
        if (in == "Storage") {
            const std::vector<const Element*> ip = params(at(sec[i], 0));
            if (!ip.empty()) {
                const double cap = to_num(value_of(ip[0]));
                net.set_capacity(nidx[i],
                                 cap < 0 ? std::numeric_limits<double>::infinity() : cap);
            }
            if (ip.size() > 1) {
                const std::vector<const Element*> pc = subs(ip[1]);
                for (std::size_t c = 0; c < pc.size() && c < K; ++c) {
                    const double v = to_num(value_of(pc[c]));
                    net.set_class_capacity(nidx[i], cidx[c],
                                           v < 0 ? std::numeric_limits<double>::infinity() : v);
                }
            }
            if (ip.size() > 2) {
                const std::vector<const Element*> dr = subs(ip[2]);
                for (std::size_t c = 0; c < dr.size() && c < K; ++c) {
                    const std::string rule = value_of(dr[c]);
                    if (rule == "BAS blocking")
                        net.set_drop_rule(nidx[i], cidx[c], lang::DropStrategy::BAS);
                    else if (rule == "drop")
                        net.set_drop_rule(nidx[i], cidx[c], lang::DropStrategy::DROP);
                    else if (rule == "waiting queue")
                        net.set_drop_rule(nidx[i], cidx[c], lang::DropStrategy::WAITQ);
                }
            }
            continue;
        }
        if (in != "Enabling") continue;

        // ---- a Transition: enabling, timing and firing sections ------------
        const std::vector<const Element*> ep = params(at(sec[i], 0));
        const std::vector<const Element*> tp = params(at(sec[i], 1));
        const std::vector<const Element*> fp = params(at(sec[i], 2));
        if (ep.size() < 2 || tp.size() < 5 || fp.empty())
            throw InputError("read_jsim: transition '" + orig_name[i] +
                             "' is missing one of the enabling, timing or firing parameters");
        const std::vector<const Element*> enmodes = subs(ep[0]);
        const std::vector<const Element*> inmodes = subs(ep[1]);
        const std::vector<const Element*> names = subs(tp[0]);
        const std::vector<const Element*> nserv = subs(tp[1]);
        const std::vector<const Element*> timing = subs(tp[2]);
        const std::vector<const Element*> fprio = subs(tp[3]);
        const std::vector<const Element*> fweight = subs(tp[4]);
        const std::vector<const Element*> fmodes = subs(fp[0]);
        const std::size_t M = enmodes.size();

        qn::TransitionParam<T> par;
        par.nmodes = M;
        // The arcs are indexed by PLACE, which is the whole node index space
        // here: a transition names its neighbours by node name, and a name that
        // is not a Place would silently become an arc to a queue.
        const std::size_t inf_marker = 0;
        (void)inf_marker;
        const double dinf = std::numeric_limits<double>::infinity();
        par.enabling.assign(M, Matrix<T>(N, K, num_traits<T>::from_int(0)));
        par.inhibiting.assign(M, Matrix<T>(N, K, num_traits<T>::from_double(dinf)));
        par.firing.assign(M, Matrix<T>(N, K, num_traits<T>::from_int(0)));

        // ONE VALUE PER CLASS in the file, and one per (place, class) in the
        // struct: a COLOURED net is what JSIM writes and what this reads back,
        // so the per-class values are kept apart instead of being collapsed onto
        // the place. `cidx` maps the document's class order onto the struct's.
        auto read_arcs = [&](const Element* mode, Matrix<T>& row, bool inhibiting) {
            const std::vector<const Element*> lvl1 = subs(mode);
            if (lvl1.empty()) return;
            for (const Element* arc : subs(lvl1[0])) {
                const std::vector<const Element*> ap = subs(arc);
                if (ap.size() < 2) continue;
                const std::string target = value_of(ap[0]);
                const std::map<std::string, std::size_t>::const_iterator it =
                    node_by_name.find(target);
                if (it == node_by_name.end()) continue;
                const std::vector<const Element*> per = subs(ap[1]);
                for (std::size_t c = 0; c < per.size() && c < K; ++c) {
                    const double w = to_num(value_of(per[c]));
                    const std::size_t rr = cidx[c] ? cidx[c] - 1 : c;
                    if (inhibiting) {
                        // JMT writes 0 or -1 for "no inhibitor arc"; a threshold
                        // of 0 would inhibit at zero tokens, i.e. always, and
                        // deadlock the transition -- the reference's own note.
                        row(it->second - 1, rr) = (w <= 0.0)
                                                      ? num_traits<T>::from_double(dinf)
                                                      : num_traits<T>::from_double(w);
                    } else {
                        row(it->second - 1, rr) = (w < 0.0)
                                                      ? num_traits<T>::from_double(dinf)
                                                      : num_traits<T>::from_double(w);
                    }
                }
            }
        };

        for (std::size_t m = 0; m < M; ++m) {
            par.modenames.push_back(m < names.size() ? value_of(names[m])
                                                     : "mode" + std::to_string(m + 1));
            read_arcs(enmodes[m], par.enabling[m], false);
            if (m < inmodes.size()) read_arcs(inmodes[m], par.inhibiting[m], true);
            if (m < fmodes.size()) read_arcs(fmodes[m], par.firing[m], false);

            const double ns = m < nserv.size() ? to_num(value_of(nserv[m])) : 1.0;
            par.nmodeservers.push_back(ns < 0 ? dinf : ns);
            par.firingprio.push_back(m < fprio.size() ? to_num(value_of(fprio[m])) : 0.0);
            par.fireweight.push_back(num_traits<T>::from_double(
                m < fweight.size() ? to_num(value_of(fweight[m])) : 1.0));

            const Element* tm = at(timing, m);
            const std::string cp = tm ? tm->attr("classPath") : std::string();
            if (cp.find("ZeroServiceTimeStrategy") != std::string::npos) {
                par.timing.push_back(lang::TimingStrategy::IMMEDIATE);
                par.firingproc.push_back(lang::Distrib<T>::immediate());
                par.firingphases.push_back(1);
            } else {
                par.timing.push_back(lang::TimingStrategy::TIMED);
                const std::vector<const Element*> ts = subs(tm);
                const lang::Distrib<T> d = read_distribution<T>(
                    at(ts, 0), at(ts, 1), "read_jsim (transition '" + orig_name[i] + "')");
                par.firingproc.push_back(d);
                par.firingphases.push_back(d.phases());
            }
        }
        par.firingdep.assign(M, std::function<T(const std::vector<T>&)>());
        nidx[i] = net.add_transition(node_name[i], par);
        node_by_name[node_name[i]] = nidx[i];
        node_by_name[orig_name[i]] = nidx[i];
    }

    // ---- pass 4: arrival and service processes -----------------------------
    for (std::size_t i = 0; i < N; ++i) {
        if (!nidx[i] || !is_station[i]) continue;
        if (is_source[i]) {
            const std::vector<const Element*> sp = params(at(sec[i], 0));
            if (sp.empty()) continue;
            const std::vector<const Element*> per = subs(sp[0]);
            for (std::size_t r = 0; r < K; ++r) {
                if (cs_referenced[r]) {
                    net.set_arrival(nidx[i], cidx[r], lang::Distrib<T>::disabled_dist());
                    continue;
                }
                net.set_arrival(nidx[i], cidx[r],
                                read_strategy_dist<T>(at(per, r),
                                                      "read_jsim (arrival at '" + orig_name[i] +
                                                          "')"));
            }
            continue;
        }
        const std::string svc = secclass(i, 1);
        if (svc == "ClassSwitch" || svc == "StatelessClassSwitcher") continue;
        // The service strategy is the LAST parameter block whose subParameters
        // are per class: parameter 0 for a Delay (which has no server count),
        // parameter 2 for a Server. Located by name rather than by position so
        // an extra block -- server visits, heterogeneous policy -- does not
        // shift it.
        const std::vector<const Element*> sp = params(at(sec[i], 1));
        const Element* strat = named(sp, "ServiceStrategy");
        if (!strat) strat = is_delay[i] ? at(sp, 0) : at(sp, 2);
        if (!strat) continue;
        const std::vector<const Element*> per = subs(strat);
        for (std::size_t r = 0; r < K; ++r) {
            net.set_service(nidx[i], cidx[r],
                            read_strategy_dist<T>(at(per, r), "read_jsim (service at '" +
                                                                  orig_name[i] + "')"));
            if (r < schedparam[i].size()) net.set_sched_param(nidx[i], cidx[r], schedparam[i][r]);
        }
    }

    // ---- pass 5: links and routing ----------------------------------------
    std::vector<std::vector<bool> > conn(N, std::vector<bool>(N, false));
    std::map<std::string, std::size_t> pos_by_name;
    for (std::size_t i = 0; i < N; ++i) {
        pos_by_name[orig_name[i]] = i;
        pos_by_name[node_name[i]] = i;
    }
    for (const Element* c : sim->child_tags("connection")) {
        const std::map<std::string, std::size_t>::const_iterator a =
            pos_by_name.find(c->attr("source"));
        const std::map<std::string, std::size_t>::const_iterator b =
            pos_by_name.find(c->attr("target"));
        if (a == pos_by_name.end() || b == pos_by_name.end()) continue;
        conn[a->second][b->second] = true;
    }

    qn::RoutingMatrix<T> P;
    const T one = num_traits<T>::from_int(1);
    for (std::size_t from = 0; from < N; ++from) {
        if (!nidx[from]) continue;
        const std::string in = secclass(from, 0);
        if (in == "JobSink" || in == "Storage" || in == "Enabling") continue;
        const std::vector<const Element*> op = params(at(sec[from], 2));
        const std::vector<const Element*> per = op.empty() ? std::vector<const Element*>()
                                                           : subs(op[0]);
        std::vector<std::size_t> targets;
        for (std::size_t j = 0; j < N; ++j)
            if (conn[from][j] && nidx[j]) targets.push_back(j);

        for (std::size_t r = 0; r < K; ++r) {
            const Element* st = at(per, r);
            const std::string rs = st ? st->attr("name") : std::string("Random");
            if (rs == "Disabled") {
                net.set_routing(nidx[from], cidx[r], lang::RoutingStrategy::DISABLED);
                continue;
            }
            if (rs == "Probabilities" || rs == "Weighted Round Robin") {
                const bool wrr = rs != "Probabilities";
                net.set_routing(nidx[from], cidx[r],
                                wrr ? lang::RoutingStrategy::WRROBIN
                                    : lang::RoutingStrategy::PROB);
                // `<subParameter name="EmpiricalEntryArray">` holds one entry
                // per destination, each a (name, value) pair.
                const std::vector<const Element*> arr = subs(st);
                const std::vector<const Element*> entries =
                    arr.empty() ? arr : subs(arr[0]);
                std::map<std::size_t, double> w;
                for (const Element* e : entries) {
                    const std::vector<const Element*> kv = subs(e);
                    if (kv.size() < 2) continue;
                    const std::map<std::string, std::size_t>::const_iterator it =
                        pos_by_name.find(value_of(kv[0]));
                    if (it == pos_by_name.end() || !nidx[it->second]) continue;
                    w[it->second] = to_num(value_of(kv[1]));
                }
                if (wrr) {
                    std::map<std::size_t, double> wt;
                    for (std::map<std::size_t, double>::const_iterator it = w.begin();
                         it != w.end(); ++it)
                        wt[nidx[it->first]] = it->second;
                    net.set_routing_weights(nidx[from], cidx[r], wt);
                    // A weighted round robin still has to reach its
                    // destinations, so the links are declared with equal shares;
                    // the WEIGHTS above are what the refresh reads.
                    for (std::size_t j = 0; j < targets.size(); ++j)
                        P.set(cidx[r], cidx[r], nidx[from], nidx[targets[j]],
                              num_traits<T>::from_double(1.0 /
                                                         static_cast<double>(targets.size())));
                } else {
                    for (std::map<std::size_t, double>::const_iterator it = w.begin();
                         it != w.end(); ++it)
                        P.set(cidx[r], cidx[r], nidx[from], nidx[it->first],
                              num_traits<T>::from_double(it->second));
                }
                continue;
            }
            // The state-dependent strategies declare only WHERE a job may go;
            // the choice is made at run time, so every reachable link carries an
            // equal share and the strategy is what a solver reads.
            lang::RoutingStrategy strat = lang::RoutingStrategy::RAND;
            if (rs == "Round Robin") strat = lang::RoutingStrategy::RROBIN;
            else if (rs == "Join the Shortest Queue (JSQ)") strat = lang::RoutingStrategy::JSQ;
            else if (rs == "Power of k") strat = lang::RoutingStrategy::SQ;
            else if (rs != "Random")
                throw UnsupportedError("read_jsim: node '" + orig_name[from] +
                                       "' uses routing strategy '" + rs +
                                       "', which this reader does not build");
            net.set_routing(nidx[from], cidx[r], strat);
            if (strat == lang::RoutingStrategy::SQ) {
                int kk = 2;  // JMT's default when <k> is absent
                for (const Element* pj : subs(st)) {
                    if (pj->attr("name") == "k") kk = static_cast<int>(to_num(value_of(pj)));
                    if (pj->attr("name") == "withMemory") {
                        const std::string wm = value_of(pj);
                        if (wm == "true" || wm == "1")
                            throw UnsupportedError(
                                "read_jsim: node '" + orig_name[from] +
                                "' selects Power-of-k WITH MEMORY, which is Anselmi & Dufour's "
                                "SQ(d,N) and not the memoryless SQ(d) this port implements; "
                                "importing it as SQ(d) would answer about a different policy");
                    }
                }
                net.set_routing_param(nidx[from], cidx[r], kk);
            }
            for (std::size_t j = 0; j < targets.size(); ++j)
                P.set(cidx[r], cidx[r], nidx[from], nidx[targets[j]],
                      targets.size() == 1
                          ? one
                          : num_traits<T>::from_double(1.0 /
                                                       static_cast<double>(targets.size())));
        }
    }
    net.link(P);

    // ---- pass 6: the preload marking --------------------------------------
    //
    // THE REFERENCE CALLS `initFromMarginal` AND THIS PORT HAS NO SUCH ENTRY.
    // Its stateful nodes derive their initial state from the class populations
    // and reference stations (`default_init_state`), and the one declared
    // marking it carries is a Place's token count -- which is the case that
    // MATTERS, because an SPN with no tokens is a dead net and its answer is
    // "nothing ever fires", not an approximation of the modelled one.
    //
    // A preload that places CLOSED jobs anywhere but their reference station is
    // therefore refused rather than dropped: dropping it solves a model whose
    // population sits somewhere else, which no metric would flag.
    const std::vector<const Element*> pre = sim->child_tags("preload");
    if (!pre.empty()) {
        for (const Element* sp : pre[0]->child_tags("stationPopulations")) {
            const std::map<std::string, std::size_t>::const_iterator it =
                pos_by_name.find(sp->attr("stationName"));
            if (it == pos_by_name.end() || !nidx[it->second]) continue;
            const std::size_t pos = it->second;
            const bool is_place = secclass(pos, 0) == "Storage";
            std::vector<T> tokens(K, num_traits<T>::from_int(0));
            for (const Element* cp : sp->child_tags("classPopulation")) {
                const std::string cn = cp->attr("refClass");
                const double pop = to_num(cp->attr("population"));
                for (std::size_t r = 0; r < K; ++r) {
                    if (xclasses[r]->attr("name") != cn) continue;
                    tokens[r] = num_traits<T>::from_double(pop);
                    if (is_place || pop == 0.0) continue;
                    // A closed class whose preload sits at its own reference
                    // station is exactly what the default derivation produces,
                    // so it is not a divergence and needs no marking.
                    const std::string type = xclasses[r]->attr("type");
                    if (type == "closed" && xclasses[r]->attr("referenceSource") == orig_name[pos])
                        continue;
                    throw UnsupportedError(
                        "read_jsim: the <preload> block places " + std::to_string(pop) +
                        " job(s) of class '" + cn + "' at '" + orig_name[pos] +
                        "', which is not that class's reference station. This port derives a "
                        "queueing station's initial state from the class populations and "
                        "reference stations and has no initFromMarginal to override it, so the "
                        "marking is refused rather than dropped -- dropping it would solve a "
                        "model whose population starts somewhere else");
                }
            }
            if (is_place) net.set_initial_marking(nidx[pos], tokens);
        }
    }
    return net;
}

/**
 * Write a piped XML model document to a temporary file and return its path.
 *
 * Shared by the JMT and PNML paths: both readers walk a DOM built by
 * `xml::parse_file`, which takes a path, and a document arriving on stdin has
 * none. Staging it is the whole of the difference; the caller removes the file.
 * The messages name no reader, since either may be the caller.
 */
inline std::string jsim_stage_stdin(const std::string& text) {
    if (text.empty())
        throw InputError("stdin carried no model document");
    const char* tmpdir = std::getenv("TMPDIR");
    std::string path = std::string(tmpdir && *tmpdir ? tmpdir : "/tmp") + "/line-cli-jsim-XXXXXX";
    std::vector<char> buf(path.begin(), path.end());
    buf.push_back('\0');
    const int fd = ::mkstemp(&buf[0]);
    if (fd < 0) throw InputError("cannot create a temporary file for the piped model");
    ::close(fd);
    path.assign(&buf[0]);
    std::ofstream out(path.c_str());
    out << text;
    out.close();
    return path;
}

}  // namespace io
}  // namespace line

#endif  // LINE_IO_JSIM_READER_H
