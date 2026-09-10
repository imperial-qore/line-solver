/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_IO_PNML_H
#define LINE_IO_PNML_H

/**
 * PNML (ISO/IEC 15909-2) place/transition nets, read and written.
 *
 * Port of matlab/src/io/pnml_save.m and pnml_load.m,
 * jar/src/main/java/jline/io/PnmlIO.java and
 * python/line_solver/io/pnml_io.py. The grammar written is
 * http://www.pnml.org/version-2009/grammar/ptnet, so that a LINE net can be
 * read by the tools built around that corpus (GreatSPN, TINA, the Model
 * Checking Contest harnesses) and a net from that corpus can be analysed here.
 *
 * THE P/T GRAMMAR IS UNCOLOURED, so what it can carry is narrower than what
 * LINE can express, and the difference is REFUSED rather than approximated:
 * more than one job class, an open class or a Source/Sink, a queueing place, a
 * firing-rate dependence, and any distribution outside the scalar-parameter
 * families listed in dist_param_names() below.
 *
 * TIMING RIDES IN A TOOLSPECIFIC BLOCK, which is where the grammar puts what it
 * does not define. Each LINE MODE becomes one PNML transition, so that the arcs
 * of a mode are the arcs of a transition as the grammar requires; the block
 * records which LINE transition and mode the PNML transition came from, so the
 * reader regroups the modes the writer split. A reader that ignores the block
 * still sees a correct untimed P/T net, and a P/T net with no such block is read
 * with every transition TIMED and EXPONENTIAL AT RATE 1, the convention of the
 * stochastic Petri net literature and GreatSPN's own default.
 *
 * PARAMETER NAMES ARE MATLAB'S, not this port's argument names. `Distrib::params`
 * is documented to hold the constructor arguments in MATLAB getParam order, so
 * the names below are attached positionally to that order; they are what makes a
 * file written by any of the four codebases readable by the other three.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <fstream>
#include <limits>
#include <map>
#include <string>
#include <vector>

#include "line/lang/lang_types.h"
#include "line/lang/qn/network_builder.h"
#include "line/lang/qn/network_struct.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"
#include "line/util/xml.h"

namespace line {
namespace io {

namespace pnml_detail {

/**
 * The parameter names of the scalar-parameter families, in `Distrib::params`
 * order. An empty vector means the family has no PNML form and the writer
 * refuses it by name.
 */
inline std::vector<std::string> dist_param_names(lang::ProcessType t) {
    std::vector<std::string> v;
    switch (t) {
        case lang::ProcessType::EXP: v.push_back("lambda"); break;
        case lang::ProcessType::DET: v.push_back("t"); break;
        case lang::ProcessType::ERLANG:
            v.push_back("alpha");
            v.push_back("r");
            break;
        case lang::ProcessType::HYPEREXP:
            v.push_back("p");
            v.push_back("lambda1");
            v.push_back("lambda2");
            break;
        case lang::ProcessType::UNIFORM:
            v.push_back("min");
            v.push_back("max");
            break;
        case lang::ProcessType::GAMMA:
            v.push_back("alpha");
            v.push_back("beta");
            break;
        case lang::ProcessType::PARETO:
            v.push_back("alpha");
            v.push_back("k");
            break;
        case lang::ProcessType::WEIBULL:
            v.push_back("alpha");
            v.push_back("r");
            break;
        case lang::ProcessType::LOGNORMAL:
            v.push_back("mu");
            v.push_back("sigma");
            break;
        default: break;
    }
    return v;
}

/** Rebuild a distribution from its name and named parameters. */
template <class T>
lang::Distrib<T> dist_from(const std::string& name, const std::map<std::string, double>& p) {
    typedef lang::Distrib<T> D;
    struct Need {
        static double get(const std::map<std::string, double>& m, const std::string& key,
                          const std::string& who) {
            std::map<std::string, double>::const_iterator it = m.find(key);
            if (it == m.end())
                throw InputError("pnml: distribution " + who + " is missing the parameter \"" + key + "\"");
            return it->second;
        }
    };
    if (name == "Immediate") return D::immediate();
    if (name == "Disabled") return D::disabled_dist();
    if (name == "Exp") return D::exp_rate(num_traits<T>::from_double(Need::get(p, "lambda", name)));
    if (name == "Det") return D::det(num_traits<T>::from_double(Need::get(p, "t", name)));
    if (name == "Erlang")
        return D::erlang(num_traits<T>::from_double(Need::get(p, "alpha", name)),
                         static_cast<std::size_t>(std::llround(Need::get(p, "r", name))));
    if (name == "HyperExp")
        return D::hyperexp(num_traits<T>::from_double(Need::get(p, "p", name)),
                           num_traits<T>::from_double(Need::get(p, "lambda1", name)),
                           num_traits<T>::from_double(Need::get(p, "lambda2", name)));
    if (name == "Uniform")
        return D::uniform(num_traits<T>::from_double(Need::get(p, "min", name)),
                          num_traits<T>::from_double(Need::get(p, "max", name)));
    if (name == "Gamma")
        return D::gamma_dist(num_traits<T>::from_double(Need::get(p, "alpha", name)),
                             num_traits<T>::from_double(Need::get(p, "beta", name)));
    if (name == "Pareto")
        return D::pareto(num_traits<T>::from_double(Need::get(p, "alpha", name)),
                         num_traits<T>::from_double(Need::get(p, "k", name)));
    if (name == "Weibull")
        // stored (alpha = scale, r = shape); the factory takes (scale, shape),
        // so naming both is what keeps a round trip from transposing them.
        return D::weibull(num_traits<T>::from_double(Need::get(p, "alpha", name)),
                          num_traits<T>::from_double(Need::get(p, "r", name)));
    if (name == "Lognormal")
        return D::lognormal(num_traits<T>::from_double(Need::get(p, "mu", name)),
                            num_traits<T>::from_double(Need::get(p, "sigma", name)));
    throw InputError("pnml: the timing block names distribution \"" + name +
                     "\", which the reader does not construct. The families it reads are Exp, Det, "
                     "Erlang, HyperExp, Uniform, Gamma, Pareto, Weibull, Lognormal, Immediate and "
                     "Disabled, which are the ones the writer writes");
}

/**
 * Shortest form that reads back exactly, so an integral count does not acquire a
 * decimal point and an infinite server count keeps the spelling the reader
 * expects.
 */
/** The five XML entities a name or an attribute value may need. */
inline std::string escape(const std::string& in) {
    std::string out;
    out.reserve(in.size());
    for (std::size_t i = 0; i < in.size(); ++i) {
        switch (in[i]) {
            case '&': out += "&amp;"; break;
            case '<': out += "&lt;"; break;
            case '>': out += "&gt;"; break;
            case '"': out += "&quot;"; break;
            default: out += in[i]; break;
        }
    }
    return out;
}

inline std::string num_text(double v) {
    if (std::isinf(v)) return v > 0 ? "Inf" : "-Inf";
    if (v == std::floor(v) && std::fabs(v) < 9.007199254740992e15) {
        char buf[32];
        std::snprintf(buf, sizeof buf, "%lld", static_cast<long long>(v));
        return std::string(buf);
    }
    char buf[64];
    std::snprintf(buf, sizeof buf, "%.17g", v);
    return std::string(buf);
}

inline double text_number(const xml::Element& e, const std::string& label, double fallback) {
    const std::vector<const xml::Element*> labels = e.child_tags(label);
    if (labels.empty()) return fallback;
    std::string txt;
    const std::vector<const xml::Element*> texts = labels[0]->child_tags("text");
    txt = texts.empty() ? labels[0]->text : texts[0]->text;
    // trim
    std::size_t a = txt.find_first_not_of(" \t\r\n");
    std::size_t b = txt.find_last_not_of(" \t\r\n");
    if (a == std::string::npos) return fallback;
    txt = txt.substr(a, b - a + 1);
    // Some tools write the marking as "3" and some as "1`3" (a coloured multiset
    // of one colour); the plain integer is the one the grammar defines.
    const std::size_t tick = txt.rfind('`');
    if (tick != std::string::npos) txt = txt.substr(tick + 1);
    try {
        return std::stod(txt);
    } catch (const std::exception&) {
        throw InputError("pnml: label <" + label + "> holds \"" + txt + "\", which is not a number");
    }
}

inline double attr_number(const xml::Element& e, const std::string& key, double fallback) {
    std::string txt = e.attr(key);
    const std::size_t a = txt.find_first_not_of(" \t\r\n");
    if (a == std::string::npos) return fallback;
    const std::size_t b = txt.find_last_not_of(" \t\r\n");
    txt = txt.substr(a, b - a + 1);
    if (txt == "Inf" || txt == "inf") return std::numeric_limits<double>::infinity();
    if (txt == "-Inf" || txt == "-inf") return -std::numeric_limits<double>::infinity();
    try {
        return std::stod(txt);
    } catch (const std::exception&) {
        throw InputError("pnml: attribute " + key + " holds \"" + txt + "\", which is not a number");
    }
}

inline std::string element_id(const xml::Element& e) {
    std::string id = e.attr("id");
    if (id.empty()) {
        const std::vector<const xml::Element*> names = e.child_tags("name");
        if (!names.empty()) {
            const std::vector<const xml::Element*> texts = names[0]->child_tags("text");
            id = texts.empty() ? names[0]->text : texts[0]->text;
        }
    }
    if (id.empty()) throw InputError("pnml: a place or transition carries neither an id nor a name");
    return id;
}

inline bool is_inhibitor(const xml::Element& arc) {
    const std::vector<const xml::Element*> types = arc.child_tags("type");
    for (std::size_t i = 0; i < types.size(); ++i) {
        std::string v = types[i]->attr("value");
        std::transform(v.begin(), v.end(), v.begin(), ::tolower);
        if (v == "inhibitor") return true;
    }
    std::string v = arc.attr("type");
    std::transform(v.begin(), v.end(), v.begin(), ::tolower);
    return v == "inhibitor";
}

/** The <mode> of a transition's LINE toolspecific block, or null. */
inline const xml::Element* line_toolspecific(const xml::Element& tr) {
    const std::vector<const xml::Element*> blocks = tr.child_tags("toolspecific");
    for (std::size_t i = 0; i < blocks.size(); ++i) {
        std::string tool = blocks[i]->attr("tool");
        std::transform(tool.begin(), tool.end(), tool.begin(), ::toupper);
        if (tool != "LINE") continue;
        const std::vector<const xml::Element*> modes = blocks[i]->child_tags("mode");
        if (!modes.empty()) return modes[0];
    }
    return 0;
}

}  // namespace pnml_detail

/**
 * Write the Petri net of a refreshed NetworkStruct to a PNML place/transition
 * file.
 *
 * @param sn   struct of a net holding only places and transitions, one closed class
 * @param path output path
 */
template <class T>
void pnml_save(const qn::NetworkStruct<T>& sn, const std::string& path) {
    const double inf = std::numeric_limits<double>::infinity();
    if (sn.classes.size() != 1)
        throw InputError("pnml_save: the PNML place/transition grammar is UNCOLOURED, so it cannot carry a "
                         "net with " + std::to_string(sn.classes.size()) +
                         " job classes: its tokens are indistinguishable. Export a single-class net, or "
                         "use the JSON writer for the full model");
    if (sn.classes[0].type != qn::JobClassType::CLOSED)
        throw InputError("pnml_save: the PNML place/transition grammar has no unbounded token source, so an "
                         "open class cannot be represented. Close the class, or use the JSON writer");

    std::vector<std::size_t> places, transitions;  // 1-based node indices
    for (std::size_t i = 0; i < sn.nodes.size(); ++i) {
        const lang::NodeType nt = sn.nodes[i].nodetype;
        if (nt == lang::NodeType::Place) {
            places.push_back(i + 1);
        } else if (nt == lang::NodeType::Transition) {
            transitions.push_back(i + 1);
        } else {
            throw InputError("pnml_save: node " + sn.nodes[i].name + " is a " +
                             lang::node_type_to_text(nt) +
                             ". A PNML place/transition net holds only places and transitions; a Source, a "
                             "Sink or a queueing station has no counterpart in the grammar");
        }
    }
    if (places.empty())
        throw InputError("pnml_save: the model holds no Place, so there is no Petri net to write");

    // Token count of each place. The recorded marking is authoritative; a place
    // with none holds the class population when it is the reference station,
    // which is the default LINE itself applies.
    std::vector<long long> marking(places.size(), 0);
    for (std::size_t p = 0; p < places.size(); ++p) {
        const typename std::map<std::size_t, std::vector<T> >::const_iterator it =
            sn.initmarking.find(places[p]);
        if (it != sn.initmarking.end() && !it->second.empty()) {
            marking[p] = std::llround(num_traits<T>::to_double(it->second[0]));
        } else if (sn.classes[0].refstat != 0 &&
                   sn.station_to_node[sn.classes[0].refstat - 1] == places[p]) {
            marking[p] = std::llround(sn.classes[0].population);
        }
    }

    // THE DOCUMENT IS BUILT AS TEXT, not through xml::Element, and the reason is
    // byte identity: the shared serializer puts every element on its own line,
    // while MATLAB, the JAR and python all write `<name><text>x</text></name>`
    // inline. Emitting the same bytes as the other three is what lets a parity
    // check diff the FILES rather than reparse them. Reading still goes through
    // xml::parse_file, where the layout does not matter.
    const std::string netname = sn.name.empty() ? std::string("net") : sn.name;
    std::vector<std::string> sb;
    sb.push_back("<?xml version=\"1.0\" encoding=\"UTF-8\"?>");
    sb.push_back("<pnml xmlns=\"http://www.pnml.org/version-2009/grammar/pnml\">");
    sb.push_back("  <net id=\"" + pnml_detail::escape(netname) +
                 "\" type=\"http://www.pnml.org/version-2009/grammar/ptnet\">");
    sb.push_back("    <name><text>" + pnml_detail::escape(netname) + "</text></name>");
    sb.push_back("    <page id=\"page0\">");

    for (std::size_t p = 0; p < places.size(); ++p) {
        const std::string nm = pnml_detail::escape(sn.nodes[places[p] - 1].name);
        sb.push_back("      <place id=\"" + nm + "\">");
        sb.push_back("        <name><text>" + nm + "</text></name>");
        sb.push_back("        <initialMarking><text>" + std::to_string(marking[p]) +
                     "</text></initialMarking>");
        sb.push_back("      </place>");
    }

    // Arcs are collected and appended after every transition, so the document
    // reads places, transitions, arcs, as the other three codebases write it.
    struct ArcRec {
        std::string source, target;
        long long weight;
        bool inhibitor;
    };
    std::vector<ArcRec> arcs;

    for (std::size_t t = 0; t < transitions.size(); ++t) {
        const std::size_t nd = transitions[t];
        const typename std::map<std::size_t, qn::TransitionParam<T> >::const_iterator it =
            sn.transparam.find(nd);
        if (it == sn.transparam.end())
            throw InputError("pnml_save: transition " + sn.nodes[nd - 1].name +
                             " declares no mode, so it has no firing behaviour to write");
        const qn::TransitionParam<T>& tp = it->second;
        for (std::size_t m = 0; m < tp.nmodes; ++m) {
            if (m < tp.firingdep.size() && tp.firingdep[m])
                throw InputError("pnml_save: transition " + sn.nodes[nd - 1].name + " mode " +
                                 std::to_string(m + 1) +
                                 " declares a marking-dependent firing rate, which no PNML element can carry");
            const std::string mname =
                m < tp.modenames.size() && !tp.modenames[m].empty() ? tp.modenames[m]
                                                                    : "Mode" + std::to_string(m + 1);
            const std::string tid =
                tp.nmodes == 1 ? sn.nodes[nd - 1].name : sn.nodes[nd - 1].name + "." + mname;
            const bool immediate =
                m < tp.timing.size() && tp.timing[m] == lang::TimingStrategy::IMMEDIATE;

            const std::string etid = pnml_detail::escape(tid);
            sb.push_back("      <transition id=\"" + etid + "\">");
            sb.push_back("        <name><text>" + etid + "</text></name>");
            sb.push_back("        <toolspecific tool=\"LINE\" version=\"3.0\">");
            sb.push_back(
                "          <mode transition=\"" + pnml_detail::escape(sn.nodes[nd - 1].name) +
                "\" name=\"" + pnml_detail::escape(mname) + "\" timing=\"" +
                (immediate ? "immediate" : "timed") + "\" servers=\"" +
                pnml_detail::num_text(m < tp.nmodeservers.size() ? tp.nmodeservers[m] : 1.0) +
                "\" priority=\"" +
                pnml_detail::num_text(m < tp.firingprio.size() ? tp.firingprio[m] : 1.0) +
                "\" weight=\"" +
                pnml_detail::num_text(m < tp.fireweight.size()
                                          ? num_traits<T>::to_double(tp.fireweight[m])
                                          : 1.0) +
                "\">");
            if (!immediate) {
                if (m >= tp.firingproc.size() || tp.firingproc[m].disabled)
                    throw InputError("pnml_save: transition " + sn.nodes[nd - 1].name + " mode " +
                                     std::to_string(m + 1) + " is timed but carries no distribution");
                const lang::Distrib<T>& d = tp.firingproc[m];
                const std::string dname = lang::process_to_text(d.type);
                const std::vector<std::string> pnames = pnml_detail::dist_param_names(d.type);
                if (pnames.empty())
                    throw InputError(
                        "pnml_save: transition " + sn.nodes[nd - 1].name + " mode " +
                        std::to_string(m + 1) + " holds a " + dname +
                        ", whose parameters are not scalars. The PNML timing block carries scalar "
                        "parameters only; a matrix-parameterized law (PH, APH, MAP, MMPP2, ME, RAP) has "
                        "no PNML representation and writing its mean rate instead would read back as a "
                        "different model");
                if (d.params.size() < pnames.size())
                    throw InputError("pnml_save: transition " + sn.nodes[nd - 1].name + " mode " +
                                     std::to_string(m + 1) + " holds a " + dname + " with " +
                                     std::to_string(d.params.size()) + " parameters, expected " +
                                     std::to_string(pnames.size()));
                sb.push_back("            <distribution name=\"" + pnml_detail::escape(dname) +
                             "\">");
                for (std::size_t k = 0; k < pnames.size(); ++k) {
                    sb.push_back("              <parameter name=\"" +
                                 pnml_detail::escape(pnames[k]) + "\" value=\"" +
                                 pnml_detail::num_text(num_traits<T>::to_double(d.params[k])) +
                                 "\"/>");
                }
                sb.push_back("            </distribution>");
            }
            sb.push_back("          </mode>");
            sb.push_back("        </toolspecific>");
            sb.push_back("      </transition>");

            for (std::size_t p = 0; p < places.size(); ++p) {
                const std::size_t row = places[p] - 1;
                const double w = m < tp.enabling.size() && row < tp.enabling[m].rows()
                                     ? num_traits<T>::to_double(tp.enabling[m](row, 0))
                                     : 0.0;
                if (w > 0) {
                    ArcRec a;
                    a.source = sn.nodes[places[p] - 1].name;
                    a.target = tid;
                    a.weight = std::llround(w);
                    a.inhibitor = false;
                    arcs.push_back(a);
                }
                const double h = m < tp.inhibiting.size() && row < tp.inhibiting[m].rows()
                                     ? num_traits<T>::to_double(tp.inhibiting[m](row, 0))
                                     : inf;
                if (std::isfinite(h)) {
                    ArcRec a;
                    a.source = sn.nodes[places[p] - 1].name;
                    a.target = tid;
                    a.weight = std::llround(h);
                    a.inhibitor = true;
                    arcs.push_back(a);
                }
                const double f = m < tp.firing.size() && row < tp.firing[m].rows()
                                     ? num_traits<T>::to_double(tp.firing[m](row, 0))
                                     : 0.0;
                if (f > 0) {
                    ArcRec a;
                    a.source = tid;
                    a.target = sn.nodes[places[p] - 1].name;
                    a.weight = std::llround(f);
                    a.inhibitor = false;
                    arcs.push_back(a);
                }
            }
        }
    }

    for (std::size_t a = 0; a < arcs.size(); ++a) {
        sb.push_back("      <arc id=\"a" + std::to_string(a + 1) + "\" source=\"" +
                     pnml_detail::escape(arcs[a].source) + "\" target=\"" +
                     pnml_detail::escape(arcs[a].target) + "\">");
        if (arcs[a].inhibitor) {
            // An inhibitor arc is not in the P/T grammar itself; <type
            // value="inhibitor"/> is the extension GreatSPN, TINA and PIPE all
            // read, so it is the one written here.
            sb.push_back("        <type value=\"inhibitor\"/>");
        }
        sb.push_back("        <inscription><text>" + std::to_string(arcs[a].weight) +
                     "</text></inscription>");
        sb.push_back("      </arc>");
    }

    sb.push_back("    </page>");
    sb.push_back("  </net>");
    sb.push_back("</pnml>");

    std::ofstream out(path.c_str());
    if (!out) throw InputError("pnml_save: cannot open " + path + " for writing");
    for (std::size_t i = 0; i < sb.size(); ++i) out << sb[i] << "\n";
}

/**
 * Read one net of a PNML place/transition document into a Network.
 *
 * @param path   input path
 * @param net_id id of the net to read; empty selects the first
 */
template <class T>
qn::Network<T> pnml_load(const std::string& path, const std::string& net_id = std::string()) {
    const double inf = std::numeric_limits<double>::infinity();
    std::unique_ptr<xml::Element> root = xml::parse_file(path);
    const std::vector<const xml::Element*> nets = root->child_tags("net");
    if (nets.empty()) throw InputError("pnml_load: " + path + " holds no <net> element");
    const xml::Element* net = 0;
    for (std::size_t i = 0; i < nets.size(); ++i)
        if (net_id.empty() || nets[i]->attr("id") == net_id) {
            net = nets[i];
            break;
        }
    if (net == 0) throw InputError("pnml_load: " + path + " holds no net with id \"" + net_id + "\"");
    const std::string nettype = net->attr("type");
    if (!nettype.empty() && nettype.find("ptnet") == std::string::npos)
        throw InputError("pnml_load: net \"" + net->attr("id") + "\" declares type " + nettype +
                         ". Only the place/transition grammar "
                         "(http://www.pnml.org/version-2009/grammar/ptnet) is read: a coloured or a "
                         "symmetric net carries token colours that a single-class LINE net cannot hold");

    // Places, transitions and arcs may sit directly under <net> or under any
    // <page>; the grammar allows both and tools differ, so the whole subtree is
    // searched rather than one level of it.
    const std::vector<const xml::Element*> place_elems = net->by_tag("place");
    const std::vector<const xml::Element*> trans_elems = net->by_tag("transition");
    const std::vector<const xml::Element*> arc_elems = net->by_tag("arc");
    if (place_elems.empty())
        throw InputError("pnml_load: net \"" + net->attr("id") + "\" holds no place");

    std::vector<std::string> place_names;
    std::vector<long long> place_marking;
    long long total = 0;
    for (std::size_t i = 0; i < place_elems.size(); ++i) {
        place_names.push_back(pnml_detail::element_id(*place_elems[i]));
        const long long mk =
            std::llround(pnml_detail::text_number(*place_elems[i], "initialMarking", 0.0));
        place_marking.push_back(mk);
        total += mk;
    }
    if (total == 0)
        throw InputError("pnml_load: the initial marking of this net is empty. A LINE closed class needs "
                         "tokens to hold, and a net with none has no reachable behaviour to analyse");

    // Each PNML transition is one LINE MODE; the toolspecific block says which
    // LINE transition it belongs to, so the modes the writer split regroup here.
    std::vector<std::string> trans_ids, mode_owner, mode_name;
    std::vector<lang::TimingStrategy> mode_timing;
    std::vector<double> mode_servers, mode_prio, mode_weight;
    std::vector<lang::Distrib<T> > mode_dist;
    for (std::size_t i = 0; i < trans_elems.size(); ++i) {
        const std::string id = pnml_detail::element_id(*trans_elems[i]);
        trans_ids.push_back(id);
        std::string owner = id, name = "Mode1";
        lang::TimingStrategy timing = lang::TimingStrategy::TIMED;
        double servers = 1.0, prio = 1.0, weight = 1.0;
        lang::Distrib<T> dist = lang::Distrib<T>::exp_rate(num_traits<T>::from_int(1));
        const xml::Element* spec = pnml_detail::line_toolspecific(*trans_elems[i]);
        if (spec != 0) {
            if (!spec->attr("transition").empty()) owner = spec->attr("transition");
            if (!spec->attr("name").empty()) name = spec->attr("name");
            std::string tm = spec->attr("timing");
            std::transform(tm.begin(), tm.end(), tm.begin(), ::tolower);
            if (tm == "immediate") timing = lang::TimingStrategy::IMMEDIATE;
            servers = pnml_detail::attr_number(*spec, "servers", 1.0);
            prio = pnml_detail::attr_number(*spec, "priority", 1.0);
            weight = pnml_detail::attr_number(*spec, "weight", 1.0);
            if (timing != lang::TimingStrategy::IMMEDIATE) {
                const std::vector<const xml::Element*> ds = spec->child_tags("distribution");
                if (!ds.empty()) {
                    std::map<std::string, double> params;
                    const std::vector<const xml::Element*> ps = ds[0]->child_tags("parameter");
                    for (std::size_t k = 0; k < ps.size(); ++k)
                        params[ps[k]->attr("name")] =
                            pnml_detail::attr_number(*ps[k], "value",
                                                     std::numeric_limits<double>::quiet_NaN());
                    dist = pnml_detail::dist_from<T>(ds[0]->attr("name"), params);
                }
            } else {
                dist = lang::Distrib<T>::immediate();
            }
        }
        mode_owner.push_back(owner);
        mode_name.push_back(name);
        mode_timing.push_back(timing);
        mode_servers.push_back(servers);
        mode_prio.push_back(prio);
        mode_weight.push_back(weight);
        mode_dist.push_back(dist);
    }

    std::vector<std::string> owner_names;
    std::vector<std::size_t> owner_of(trans_ids.size(), 0);
    for (std::size_t i = 0; i < trans_ids.size(); ++i) {
        std::size_t k = owner_names.size();
        for (std::size_t j = 0; j < owner_names.size(); ++j)
            if (owner_names[j] == mode_owner[i]) {
                k = j;
                break;
            }
        if (k == owner_names.size()) owner_names.push_back(mode_owner[i]);
        owner_of[i] = k;
    }

    std::string name = net_id.empty() ? net->attr("id") : net_id;
    if (name.empty()) name = "pnml";
    qn::Network<T> model(name);

    std::vector<std::size_t> place_node;
    for (std::size_t i = 0; i < place_names.size(); ++i)
        place_node.push_back(model.add_place(place_names[i]));

    // The reference station is the first place holding tokens, so the class
    // starts where the marking says it does.
    std::size_t ref = 0;
    for (std::size_t i = 0; i < place_marking.size(); ++i)
        if (place_marking[i] > 0) {
            ref = i;
            break;
        }
    const std::size_t cls =
        model.add_closed_class("Class1", static_cast<double>(total), place_node[ref], 0);

    // A transition is itself a node, so the arc matrices must be sized against
    // the FINISHED node count, which is known before any of them is added.
    const std::size_t nnodes = place_names.size() + owner_names.size();

    // The arcs are read first, because a mode's matrices are built whole and
    // handed to add_transition rather than set afterwards.
    std::vector<qn::TransitionParam<T> > params(owner_names.size());
    std::vector<std::size_t> mode_index(trans_ids.size(), 0);
    for (std::size_t i = 0; i < trans_ids.size(); ++i) {
        qn::TransitionParam<T>& tp = params[owner_of[i]];
        mode_index[i] = tp.nmodes;
        tp.nmodes++;
        tp.modenames.push_back(mode_name[i]);
        tp.timing.push_back(mode_timing[i]);
        tp.firingproc.push_back(mode_dist[i]);
        tp.firingphases.push_back(
            mode_timing[i] == lang::TimingStrategy::IMMEDIATE || mode_dist[i].disabled
                ? 0
                : lang::dist_to_map(mode_dist[i]).order());
        tp.nmodeservers.push_back(mode_servers[i]);
        tp.firingprio.push_back(mode_prio[i]);
        tp.fireweight.push_back(num_traits<T>::from_double(mode_weight[i]));
        tp.enabling.push_back(Matrix<T>(nnodes, 1, num_traits<T>::from_int(0)));
        tp.inhibiting.push_back(Matrix<T>(nnodes, 1, num_traits<T>::from_double(inf)));
        tp.firing.push_back(Matrix<T>(nnodes, 1, num_traits<T>::from_int(0)));
        tp.firingdep.push_back(std::function<T(const std::vector<T>&)>());
    }

    // Node index of a LINE transition, in the order add_transition will assign:
    // the places come first, then one node per owner.
    std::vector<std::size_t> trans_node(owner_names.size(), 0);
    for (std::size_t i = 0; i < owner_names.size(); ++i)
        trans_node[i] = place_names.size() + i + 1;

    qn::RoutingMatrix<T> R;
    for (std::size_t a = 0; a < arc_elems.size(); ++a) {
        const std::string src = arc_elems[a]->attr("source");
        const std::string tgt = arc_elems[a]->attr("target");
        const double w = pnml_detail::text_number(*arc_elems[a], "inscription", 1.0);
        if (!(w > 0))
            throw InputError("pnml_load: arc " + src + " -> " + tgt +
                             " carries a non-positive inscription");
        std::size_t ip = place_names.size(), it = trans_ids.size();
        for (std::size_t k = 0; k < place_names.size(); ++k)
            if (place_names[k] == src) ip = k;
        for (std::size_t k = 0; k < trans_ids.size(); ++k)
            if (trans_ids[k] == tgt) it = k;
        if (ip < place_names.size() && it < trans_ids.size()) {
            qn::TransitionParam<T>& tp = params[owner_of[it]];
            const std::size_t row = place_node[ip] - 1;
            if (pnml_detail::is_inhibitor(*arc_elems[a]))
                tp.inhibiting[mode_index[it]](row, 0) = num_traits<T>::from_double(w);
            else
                tp.enabling[mode_index[it]](row, 0) = num_traits<T>::from_double(w);
            R.set(cls, cls, place_node[ip], trans_node[owner_of[it]], num_traits<T>::from_int(1));
            continue;
        }
        ip = place_names.size();
        it = trans_ids.size();
        for (std::size_t k = 0; k < trans_ids.size(); ++k)
            if (trans_ids[k] == src) it = k;
        for (std::size_t k = 0; k < place_names.size(); ++k)
            if (place_names[k] == tgt) ip = k;
        if (ip < place_names.size() && it < trans_ids.size()) {
            qn::TransitionParam<T>& tp = params[owner_of[it]];
            tp.firing[mode_index[it]](place_node[ip] - 1, 0) = num_traits<T>::from_double(w);
            R.set(cls, cls, trans_node[owner_of[it]], place_node[ip], num_traits<T>::from_int(1));
            continue;
        }
        throw InputError("pnml_load: arc " + src + " -> " + tgt +
                         " connects two places or two transitions, which the place/transition grammar "
                         "does not allow");
    }

    for (std::size_t i = 0; i < owner_names.size(); ++i) {
        const std::size_t nd = model.add_transition(owner_names[i], params[i]);
        if (nd != trans_node[i])
            throw InputError("pnml_load: internal node numbering disagreed with the arc matrices");
    }

    model.link(R);
    for (std::size_t i = 0; i < place_node.size(); ++i)
        model.set_initial_marking(
            place_node[i], std::vector<T>(1, num_traits<T>::from_double(
                                                 static_cast<double>(place_marking[i]))));
    return model;
}

}  // namespace io
}  // namespace line

#endif  // LINE_IO_PNML_H
