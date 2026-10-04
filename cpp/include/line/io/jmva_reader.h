/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_IO_JMVA_READER_H
#define LINE_IO_JMVA_READER_H

/**
 * @file
 * @ingroup line_io
 * Port of `matlab/src/io/JMVA2LINE.m` (and `jline.io.M2M.JMVA2LINE`): a JMT
 * `.jmva` product-form model read into a `qn::Network`, the inverse of
 * `jmva_writer.h`.
 *
 * WHAT A JMVA DOCUMENT CAN SAY IS LESS THAN A NETWORK. It carries service
 * times and visit ratios per station and class, not routing, so the model is
 * rebuilt the way the reference rebuilds it: each station's service becomes an
 * exponential of mean `servicetime * visits` (the DEMAND), and each class is
 * routed serially through the stations it visits, in document-type order
 * (delay stations, then load-independent, then load-dependent), open classes
 * from a `Source` to a `Sink`. Every station-level mean of the product-form
 * solution depends on the demands alone, so this is exact for what the format
 * carries; per-visit response times are those of a visit ratio of one.
 *
 * STATION ORDER IS BY TYPE, NOT BY DOCUMENT POSITION, as in the reference:
 * `<delaystation>` elements first, then `<listation>`, then `<ldstation>`.
 *
 * LOAD-DEPENDENT STATIONS follow MATLAB, not the JAR (which refuses them). The
 * per-class `<servicetimes customerclass=...>` body is the semicolon (or comma)
 * separated list d(1);d(2);... of per-population service times, with a leading
 * d(0) = 0 dropped when present. The station is served at d(1) * visits and
 * load-scaled by alpha(n) = d(1) / d(n), which is how `jmva_writer.h` encodes a
 * multiserver queue. The `servers` attribute becomes the server count. A
 * `<servicetime>` scalar under an `<ldstation>`, the spelling `lqn2ps -Ojmva`
 * emits for a multiserver station, is read as d(1) alone: `servers` servers and
 * no load dependence.
 *
 * TWO DELIBERATE READINGS where the reference indexes by position:
 *   - servicetime and visit entries are matched to a class by their
 *     `customerclass` attribute, not by position. The reference takes the r-th
 *     `<visit>` for the class named by the r-th `<servicetime>`, which agrees
 *     wherever both lists are in class order and is wrong otherwise.
 *   - a numeric attribute or body that does not parse (a SPEX variable such as
 *     `population="$N"`) is refused by name. The reference's `Str2Num` leaves
 *     it a string and then fails later in arithmetic with no mention of where.
 *
 * REFERENCE STATION. The reference reads `<ReferenceStation>` and then ignores
 * it, giving every closed class the FIRST station created; the JAR does the
 * same. This port reproduces that, so system-level metrics agree across the
 * codebases on the same file.
 */

#include <cstdlib>
#include <memory>
#include <string>
#include <vector>

#include "line/lang/lang_types.h"
#include "line/lang/qn/network_builder.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/xml.h"

namespace line {
namespace io {

namespace jmva_detail {

using xml::Element;

/** Strict decimal parse: the whole (trimmed) text must be a number. */
inline double to_double_strict(const std::string& s, const std::string& what) {
    std::size_t b = 0, e = s.size();
    while (b < e && (s[b] == ' ' || s[b] == '\t' || s[b] == '\n' || s[b] == '\r')) ++b;
    while (e > b && (s[e - 1] == ' ' || s[e - 1] == '\t' || s[e - 1] == '\n' || s[e - 1] == '\r'))
        --e;
    const std::string t = s.substr(b, e - b);
    if (t.empty()) throw InputError("jmva2line: " + what + " is empty");
    char* end = nullptr;
    const double v = std::strtod(t.c_str(), &end);
    if (end == t.c_str() || *end != '\0')
        throw InputError("jmva2line: " + what + " is '" + t +
                         "', which is not a number (a SPEX template variable must be "
                         "substituted before import)");
    return v;
}

/** The first direct child with the given tag, or null. */
inline const Element* child(const Element* e, const std::string& tag) {
    if (!e) return nullptr;
    const std::vector<const Element*> v = e->child_tags(tag);
    return v.empty() ? nullptr : v[0];
}

/** The entry of `tag` under `parent` whose `customerclass` is `cls`, or null. */
inline const Element* by_class(const Element* parent, const std::string& tag,
                               const std::string& cls) {
    if (!parent) return nullptr;
    for (const Element* e : parent->child_tags(tag))
        if (e->attr("customerclass") == cls) return e;
    return nullptr;
}

/** `str2double(strsplit(s, ';'))` (or ','), NaN entries dropped as the reference does. */
inline std::vector<double> split_demands(const std::string& s) {
    const char sep = s.find(';') != std::string::npos ? ';' : ',';
    std::vector<double> out;
    std::size_t start = 0;
    while (start <= s.size()) {
        std::size_t end = s.find(sep, start);
        if (end == std::string::npos) end = s.size();
        const std::string tok = s.substr(start, end - start);
        char* stop = nullptr;
        const double v = std::strtod(tok.c_str(), &stop);
        bool blank = true;
        for (char c : tok)
            if (!(c == ' ' || c == '\t' || c == '\n' || c == '\r')) blank = false;
        if (!blank) {
            while (*stop == ' ' || *stop == '\t' || *stop == '\n' || *stop == '\r') ++stop;
            if (stop == tok.c_str() || *stop != '\0')
                throw InputError("jmva2line: load-dependent service time list '" + s +
                                 "' holds the non-numeric entry '" + tok + "'");
            out.push_back(v);
        }
        start = end + 1;
    }
    return out;
}

/** Basename without directory or extension, MATLAB's `[~,fname] = fileparts`. */
inline std::string file_stem(const std::string& path) {
    const std::size_t slash = path.find_last_of("/\\");
    std::string f = slash == std::string::npos ? path : path.substr(slash + 1);
    const std::size_t dot = f.find_last_of('.');
    if (dot != std::string::npos && dot > 0) f = f.substr(0, dot);
    return f;
}

}  // namespace jmva_detail

/**
 * Port of `JMVA2LINE(filename, modelName)`.
 *
 * @param path  the `.jmva` file
 * @param name  the model name; empty takes the file's basename, as the reference does
 * @return the network, linked and ready for `get_struct()`
 */
template <class T>
qn::Network<T> jmva2line(const std::string& path, const std::string& name = std::string()) {
    using namespace jmva_detail;
    typedef lang::Distrib<T> D;
    std::unique_ptr<Element> doc = xml::parse_file(path);
    if (!doc) throw InputError("jmva2line: cannot parse '" + path + "' as XML");

    // `xDoc = xDoc.sim` when present: the parameters may sit under <sim>.
    const Element* root = doc.get();
    if (const Element* sim = child(root, "sim")) root = sim;
    const Element* params = child(root, "parameters");
    if (!params) throw InputError("jmva2line: '" + path + "' carries no <parameters> element");
    const Element* xstations = child(params, "stations");
    const Element* xclasses = child(params, "classes");
    if (!xstations) throw InputError("jmva2line: '" + path + "' carries no <stations> element");
    if (!xclasses) throw InputError("jmva2line: '" + path + "' carries no <classes> element");

    qn::Network<T> model(name.empty() ? file_stem(path) : name);

    const std::vector<const Element*> xdelay = xstations->child_tags("delaystation");
    const std::vector<const Element*> xli = xstations->child_tags("listation");
    const std::vector<const Element*> xld = xstations->child_tags("ldstation");

    // ---- stations, by type ----------------------------------------------
    std::vector<std::size_t> nodes;             // node index of each station, creation order
    std::vector<const Element*> xnode;          // its XML element
    for (const Element* e : xdelay) {
        nodes.push_back(model.add_delay(e->attr("name")));
        xnode.push_back(e);
    }
    for (const Element* e : xli) {
        nodes.push_back(model.add_queue(e->attr("name"), lang::SchedStrategy::PS));
        xnode.push_back(e);
    }
    for (const Element* e : xld) {
        const std::size_t nd = model.add_queue(e->attr("name"), lang::SchedStrategy::PS);
        double nservers = 1.0;
        if (e->has_attr("servers"))
            nservers = to_double_strict(e->attr("servers"),
                                        "the 'servers' of ldstation '" + e->attr("name") + "'");
        model.set_number_of_servers(nd, nservers);
        nodes.push_back(nd);
        xnode.push_back(e);
    }
    const std::size_t nInf = xdelay.size(), nLI = xli.size();
    const std::size_t M = nodes.size();

    // ---- classes ----------------------------------------------------------
    std::vector<std::string> cname;
    std::vector<std::size_t> cls;
    const std::vector<const Element*> xopen = xclasses->child_tags("openclass");
    const std::vector<const Element*> xclosed = xclasses->child_tags("closedclass");
    std::size_t source = 0, sink = 0;
    if (!xopen.empty()) {
        source = model.add_source("Source");
        sink = model.add_sink("Sink");
        for (const Element* e : xopen) {
            const std::string nm = e->attr("name");
            const std::size_t r = model.add_open_class(nm, 0);
            const double rate = to_double_strict(e->attr("rate"), "the rate of open class '" + nm + "'");
            model.set_arrival(source, r, D::exp_rate(num_traits<T>::from_double(rate)));
            cname.push_back(nm);
            cls.push_back(r);
        }
    }
    if (!xclosed.empty() && M == 0)
        throw InputError("jmva2line: '" + path + "' declares closed classes but no station");
    for (const Element* e : xclosed) {
        const std::string nm = e->attr("name");
        const double pop =
            to_double_strict(e->attr("population"), "the population of closed class '" + nm + "'");
        // The reference's refstat: the first station created (see the file comment).
        const std::size_t r = model.add_closed_class(nm, pop, nodes[0], 0);
        cname.push_back(nm);
        cls.push_back(r);
    }
    const std::size_t K = cls.size();
    const std::size_t nOpen = xopen.size();

    // ---- service ----------------------------------------------------------
    std::vector<std::vector<bool>> visited(M, std::vector<bool>(K, false));
    for (std::size_t i = 0; i < M; ++i) {
        const Element* st = xnode[i];
        const std::string sname = st->attr("name");
        const Element* xst = child(st, "servicetimes");
        const Element* xvis = child(st, "visits");
        const bool is_ld = i >= nInf + nLI;
        for (std::size_t r = 0; r < K; ++r) {
            const Element* v = by_class(xvis, "visit", cname[r]);
            if (!v)
                throw InputError("jmva2line: station '" + sname + "' gives no visit for class '" +
                                 cname[r] + "'");
            const double visits =
                to_double_strict(v->text, "the visits of class '" + cname[r] + "' at '" + sname + "'");
            if (!is_ld) {
                const Element* s = by_class(xst, "servicetime", cname[r]);
                if (!s)
                    throw InputError("jmva2line: station '" + sname +
                                     "' gives no service time for class '" + cname[r] + "'");
                const double stime = to_double_strict(
                    s->text, "the service time of class '" + cname[r] + "' at '" + sname + "'");
                if (visits > 0.0) {
                    model.set_service(nodes[i], cls[r],
                                      D::exp_mean(num_traits<T>::from_double(stime * visits)));
                    visited[i][r] = true;
                } else {
                    model.set_service(nodes[i], cls[r], D::disabled_dist());
                }
                continue;
            }
            // Load-dependent: the per-population list, or the scalar lqn2ps form.
            std::vector<double> demands;
            if (const Element* s = by_class(xst, "servicetimes", cname[r])) {
                demands = split_demands(s->text);
                if (demands.empty())
                    for (const Element* c : s->child_tags("servicetime"))
                        demands.push_back(to_double_strict(
                            c->text, "a service time of class '" + cname[r] + "' at '" + sname + "'"));
            } else if (const Element* s1 = by_class(xst, "servicetime", cname[r])) {
                demands.push_back(to_double_strict(
                    s1->text, "the service time of class '" + cname[r] + "' at '" + sname + "'"));
            } else {
                throw InputError("jmva2line: ldstation '" + sname +
                                 "' gives no service times for class '" + cname[r] + "'");
            }
            if (demands.size() > 1 && demands[0] == 0.0) demands.erase(demands.begin());
            if (visits > 0.0 && !demands.empty() && demands[0] > 0.0) {
                const double base = demands[0] * visits;
                model.set_service(nodes[i], cls[r], D::exp_mean(num_traits<T>::from_double(base)));
                visited[i][r] = true;
                if (demands.size() > 1) {
                    std::vector<T> alpha(demands.size());
                    for (std::size_t n = 0; n < demands.size(); ++n)
                        alpha[n] = num_traits<T>::from_double(
                            demands[n] > 0.0 ? demands[0] / demands[n] : 0.0);
                    // As the reference: setLoadDependence is per STATION, so the
                    // last class with a list is the one that sets it.
                    model.set_load_dependence(nodes[i], alpha);
                }
            } else {
                model.set_service(nodes[i], cls[r], D::disabled_dist());
            }
        }
    }

    // ---- routing: serial through the visited stations -------------------
    qn::RoutingMatrix<T> P = model.init_routing_matrix();
    for (std::size_t r = 0; r < K; ++r) {
        std::vector<std::size_t> path_nodes;
        if (r < nOpen) path_nodes.push_back(source);
        for (std::size_t i = 0; i < M; ++i)
            if (visited[i][r]) path_nodes.push_back(nodes[i]);
        if (r < nOpen) path_nodes.push_back(sink);
        if (path_nodes.empty()) continue;
        P.set(cls[r], model.serial_routing(path_nodes));
    }
    model.link(P);
    return model;
}

}  // namespace io
}  // namespace line

#endif  // LINE_IO_JMVA_READER_H
