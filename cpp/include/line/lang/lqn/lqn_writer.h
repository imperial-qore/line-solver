/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_LANG_LQN_LQN_WRITER_H
#define LINE_LANG_LQN_LQN_WRITER_H

/**
 * LqnModel -> .lqnx, a port of matlab/src/lang/layered/@@LayeredNetwork/writeXML.m.
 *
 * WHY IT WRITES THE INTERMEDIATE MODEL AND NOT THE STRUCT. getStruct flattens a
 * precedence block into edges of `graph`, so an AND-fork and two independent
 * sequences leave the same trace there, and a POST_LOOP loses its counts to the
 * branch shares. lqns rejects a document whose activity graph does not name its
 * blocks, so a writer working from the struct would have to guess them; the
 * reference writes from the handle graph for the same reason, and `LqnModel` is
 * this port's stand-in for it (lqn_reader.h).
 *
 * WHAT THE SCHEMA CANNOT CARRY, and what this does about it:
 *
 * - A think time on a NON-reference task. lqns rejects `think-time` there
 *   outright ('Task "X" is not a reference task'), so the attribute is dropped
 *   and the caller is told through `LqnWriteReport::dropped`, never silently.
 *   See _kb, "lqnx cannot carry non-ref think time".
 * - PRE_OR branch shares. An OR-JOIN takes whichever branch arrives, so the
 *   schema puts no `prob` on a `pre-OR` activity; writeXML.m omits them too.
 *   The reader accepts them when present, so a document that carries them
 *   round-trips through THIS port and not through the reference.
 * - Cache tasks, item entries and admission constraints. They reach this port
 *   through the JSON interchange or the builder, and the LQN schema has no
 *   element for any of them. A model that declares one is REFUSED by name
 *   rather than written as a plain task, because lqns would answer the
 *   resulting document and the answer would describe a different model.
 *   A SETUP TASK IS NOT IN THAT LIST: `<setup>`/`<delay-off>` are a LINE
 *   extension that every codebase here writes and reads, so the model survives
 *   the round trip; lqns ignores the two elements and answers the model
 *   without the cold start, which is what it would do with them absent too.
 *
 * REPLIES. lqns requires every synchronously-called entry of a non-reference
 * task to name its reply activity. A model built in code, or read from a
 * document that left them implicit, has none declared, so the implicit rule of
 * getStruct.m:641-671 is reproduced here: a leaf activity of the task (no
 * successor within the same task) replies to the entry reached by walking the
 * graph backwards, and an entry that declares any reply keeps its own.
 */

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "line/lang/lqn/lqn_reader.h"
#include "line/util/error.h"
#include "line/util/xml.h"

namespace line {
namespace lqn {

/** What the schema could not carry, one human-readable line per loss. */
struct LqnWriteReport {
    std::vector<std::string> dropped;
};

namespace detail {

/**
 * Shortest decimal that reads back as the same double.
 *
 * num2str, which the reference uses, keeps five significant digits, so a demand
 * of 1/3 reaches lqns as 0.33333 and the answer differs in the fourth digit
 * from the one this port computes in-process. The file is a wire format for a
 * solver, not a display, so it carries the value.
 */
inline std::string lqnx_num(double v) {
    if (std::isinf(v)) return v > 0 ? "inf" : "-inf";
    char buf[40];
    for (int prec = 15; prec <= 17; ++prec) {
        std::snprintf(buf, sizeof(buf), "%.*g", prec, v);
        if (std::strtod(buf, nullptr) == v) return std::string(buf);
    }
    return std::string(buf);
}

template <class T>
std::string lqnx_num_of(const T& v) {
    return lqnx_num(num_traits<T>::to_double(v));
}

/** The `pre` / `post` element name of a precedence kind. */
inline const char* precedence_tag(PrecedenceType t) {
    switch (t) {
        case PrecedenceType::PRE_SEQ: return "pre";
        case PrecedenceType::PRE_AND: return "pre-AND";
        case PrecedenceType::PRE_OR: return "pre-OR";
        case PrecedenceType::POST_SEQ: return "post";
        case PrecedenceType::POST_AND: return "post-AND";
        case PrecedenceType::POST_OR: return "post-OR";
        case PrecedenceType::POST_LOOP: return "post-LOOP";
        // `post-CACHE` IS A LINE EXTENSION OF THE SCHEMA, and all three reference
        // codebases already write and read it: ActivityPrecedenceType.toText
        // names it, writeXML.m:306 and layered.py:3287 emit it, parseXML.m:449
        // and layered.py:4179 take it back. Refusing it here left this port the
        // only one that could not round-trip a layered cache-queueing model
        // through .lqnx, and its own reader threw on the element.
        case PrecedenceType::POST_CACHE: return "post-CACHE";
        default:
            throw UnsupportedError("lqn writer: precedence with no schema element");
    }
}

/**
 * A call group's `strategy` name, spelled as the JSON interchange spells a
 * routing strategy: the enum CONSTANT, upper case, not the lower-case name
 * `lang::routing_to_text` prints in a struct dump.
 *
 * Only the two strategies a call group can be built with are named: WRROBIN
 * would need per-target weights the group API does not take, and the remaining
 * strategies are not dispatch policies at all, so an unnamed one is an error
 * rather than a silent PROB.
 */
inline const char* callgroup_to_lqnx(lang::RoutingStrategy r) {
    if (r == lang::RoutingStrategy::RROBIN) return "RROBIN";
    if (r == lang::RoutingStrategy::JSQ) return "JSQ";
    throw UnsupportedError("lqn writer: call groups carry RROBIN or JSQ; routing strategy '" +
                           std::string(lang::routing_to_text(r)) + "' cannot be written to .lqnx");
}

/**
 * The reply activities of each entry, declared where declared and inferred
 * where not, keyed by entry name.
 *
 * Port of the implicit-reply pass of getStruct.m: a LEAF activity of a task
 * (one with no successor among that task's own activities) replies, and the
 * entry it replies to is found by walking `graph` backwards through first
 * ancestors until an ENTRY is reached. An entry with an explicit reply keeps
 * it, which is what makes a phase-2 activity possible.
 */
template <class T>
std::map<std::string, std::vector<std::string>> reply_activities(const LqnModel<T>& m,
                                                                 const LqnStruct<T>& sn) {
    std::map<std::string, std::vector<std::string>> out;
    std::set<std::string> explicit_reply;
    for (std::size_t e = 0; e < m.entries.size(); ++e)
        if (!m.entries[e].reply_activities.empty()) {
            out[m.entries[e].name] = m.entries[e].reply_activities;
            explicit_reply.insert(m.entries[e].name);
        }

    for (std::size_t t = 1; t <= sn.ntasks; ++t) {
        const std::size_t tidx = sn.tshift + t;
        for (std::size_t aidx : sn.actsof[tidx]) {
            bool is_reply = true;
            const std::vector<std::size_t> post = sn.graph.succ(aidx);
            for (std::size_t p : post)
                if (std::find(sn.actsof[tidx].begin(), sn.actsof[tidx].end(), p) !=
                    sn.actsof[tidx].end())
                    is_reply = false;
            if (!is_reply) continue;
            // A leaf: walk back to the entry it belongs to, first ancestor only,
            // exactly as the reference does.
            std::size_t parent = aidx;
            std::size_t hops = 0;
            while (sn.type[parent] != LqnElement::ENTRY && hops++ <= sn.nidx) {
                const std::vector<std::size_t> anc = sn.graph.pred(parent);
                if (anc.empty()) break;
                parent = anc[0];
            }
            if (sn.type[parent] != LqnElement::ENTRY) continue;
            const std::string& ename = sn.names[parent];
            if (explicit_reply.count(ename)) continue;
            out[ename].push_back(sn.names[aidx]);
        }
    }
    return out;
}

}  // namespace detail

/**
 * Write a layered model as a .lqnx document.
 *
 * @param m                    the intermediate model, from the builder or from
 *                             read_lqnx_model
 * @param path                 file to create
 * @param model_name           the `name` attribute of `<lqn-model>`
 * @param use_abstract_names   rename elements P1/T1/E1/A1, as writeXML's third
 *                             argument does, for a document that is compared
 *                             rather than read
 * @return                     what the schema could not carry
 */
template <class T>
LqnWriteReport write_lqnx(const LqnModel<T>& m, const std::string& path,
                          const std::string& model_name = std::string("LQN"),
                          bool use_abstract_names = false) {
    LqnWriteReport report;
    const LqnStruct<T> sn = lqn_finalize(m);

    // ---- constructs with no element in the schema are refused by name -----
    for (std::size_t t = 0; t < m.tasks.size(); ++t) {
        const detail::RawTask<T>& tk = m.tasks[t];
        if (tk.nitems > 0)
            throw UnsupportedError("write_lqnx: task '" + tk.name +
                                   "' is a CacheTask, which the LQN XML schema cannot express; "
                                   "solve it with SolverLN, which models the cache directly");
        if (!tk.linconrows.empty() || tk.lincon_A.rows() > 0)
            throw UnsupportedError("write_lqnx: task '" + tk.name +
                                   "' declares an admission constraint, which the LQN XML schema "
                                   "cannot express; it travels as JSON `admissionConstraints`");
    }
    for (std::size_t e = 0; e < m.entries.size(); ++e)
        if (m.entries[e].cardinality > 0)
            throw UnsupportedError("write_lqnx: entry '" + m.entries[e].name +
                                   "' is an ItemEntry, which the LQN XML schema cannot express");
    if (!m.proc_linconrows.empty() || !m.proc_lincon.empty())
        throw UnsupportedError(
            "write_lqnx: a processor declares an admission constraint, which the LQN XML schema "
            "cannot express; it travels as JSON `admissionConstraints`");

    // ---- name map, either identity or the abstract P1/T1/E1/A1 -----------
    //
    // ONE MAP PER KIND, not the reference's single nodeHashMap. An LQN
    // routinely names a processor, its task and that task's entry alike (`c0`
    // throughout the lqngen corpus), and every reference in the document is to
    // a KNOWN kind: `dest` is an entry, `bound-to-entry` an entry, a precedence
    // operand an activity. With one map the last declaration wins and the
    // abstract-name mode would emit `E1` where the processor should be, writing
    // a document that names elements which do not exist.
    std::map<std::string, std::string> nm_host, nm_task, nm_entry, nm_act;
    {
        std::size_t tctr = 0, ectr = 0, actr = 0;
        for (std::size_t p = 0; p < m.procs.size(); ++p) {
            char buf[32];
            std::snprintf(buf, sizeof(buf), "P%zu", p + 1);
            nm_host[m.procs[p].name] = use_abstract_names ? buf : m.procs[p].name;
            for (std::size_t t = 0; t < m.tasks.size(); ++t) {
                if (m.tasks[t].proc_slot != p) continue;
                std::snprintf(buf, sizeof(buf), "T%zu", ++tctr);
                nm_task[m.tasks[t].name] = use_abstract_names ? buf : m.tasks[t].name;
                for (std::size_t e = 0; e < m.entries.size(); ++e) {
                    if (m.entries[e].task_slot != t) continue;
                    std::snprintf(buf, sizeof(buf), "E%zu", ++ectr);
                    nm_entry[m.entries[e].name] = use_abstract_names ? buf : m.entries[e].name;
                }
                for (std::size_t a = 0; a < m.acts.size(); ++a) {
                    if (m.acts[a].task_slot != t) continue;
                    std::snprintf(buf, sizeof(buf), "A%zu", ++actr);
                    nm_act[m.acts[a].name] = use_abstract_names ? buf : m.acts[a].name;
                }
            }
        }
    }
    // A name the map has never seen is a dangling reference, not a name to
    // invent: the document would name an element that does not exist and lqns
    // would refuse the file with no indication of which model built it.
    auto lookup = [](const std::map<std::string, std::string>& nm, const std::string& raw,
                     const char* kind) -> const std::string& {
        const std::map<std::string, std::string>::const_iterator it = nm.find(raw);
        if (it == nm.end())
            throw InputError("write_lqnx: '" + raw + "' is referenced as a " + kind +
                             " but no " + kind + " declares it");
        return it->second;
    };
    auto host_name = [&](const std::string& r) -> const std::string& {
        return lookup(nm_host, r, "processor");
    };
    auto task_name = [&](const std::string& r) -> const std::string& {
        return lookup(nm_task, r, "task");
    };
    auto entry_name = [&](const std::string& r) -> const std::string& {
        return lookup(nm_entry, r, "entry");
    };
    auto act_name = [&](const std::string& r) -> const std::string& {
        return lookup(nm_act, r, "activity");
    };

    const std::map<std::string, std::vector<std::string>> replies =
        detail::reply_activities(m, sn);

    // Which entries lqns needs a reply for: called synchronously, or the target
    // of a forwarding chain. An entry reached only by an asynchronous call has
    // nobody to reply TO, and a reply-entry declared for it is an error there.
    std::set<std::string> needs_reply;
    for (std::size_t e = 1; e <= sn.nentries; ++e) {
        const std::size_t eidx = sn.eshift + e;
        bool wanted = sn.issynccaller.any_col(eidx);
        for (std::size_t c = 1; !wanted && c <= sn.ncalls; ++c)
            if (sn.calltype[c] == CallType::FWD && sn.callpair_dst[c] == eidx) wanted = true;
        if (wanted) needs_reply.insert(sn.names[eidx]);
    }

    xml::Element root;
    root.name = "lqn-model";
    root.set_attr("xmlns:xsi", "http://www.w3.org/2001/XMLSchema-instance");
    root.set_attr("xsi:noNamespaceSchemaLocation", "lqn.xsd");
    root.set_attr("name", model_name);

    for (std::size_t p = 0; p < m.procs.size(); ++p) {
        const detail::RawProc& pr = m.procs[p];
        xml::Element& pe = root.add_child("processor");
        pe.set_attr("name", host_name(pr.name));
        pe.set_attr("scheduling", lang::sched_to_lqnx(pr.sched));
        if (pr.repl > 1.0) pe.set_attr("replication", detail::lqnx_num(pr.repl));
        if (pr.sched != SchedStrategy::INF) {
            // An infinite multiplicity on a finite discipline has no spelling;
            // the reference writes 1, which is what a single server means.
            const double mult = std::isinf(pr.mult) ? 1.0 : pr.mult;
            pe.set_attr("multiplicity", detail::lqnx_num(mult));
        }
        if (pr.sched == SchedStrategy::PS && pr.quantum > 0.0)
            pe.set_attr("quantum", detail::lqnx_num(pr.quantum));
        pe.set_attr("speed-factor", detail::lqnx_num(pr.speed_factor));

        for (std::size_t t = 0; t < m.tasks.size(); ++t) {
            if (m.tasks[t].proc_slot != p) continue;
            const detail::RawTask<T>& tk = m.tasks[t];
            xml::Element& te = pe.add_child("task");
            te.set_attr("name", task_name(tk.name));
            te.set_attr("scheduling", lang::sched_to_lqnx(tk.sched));
            if (tk.repl > 1.0) te.set_attr("replication", detail::lqnx_num(tk.repl));
            if (tk.sched != SchedStrategy::INF)
                te.set_attr("multiplicity",
                            detail::lqnx_num(std::isinf(tk.mult) ? 1.0 : tk.mult));
            const double think = num_traits<T>::to_double(tk.thinktime.mean);
            if (tk.sched == SchedStrategy::REF) {
                te.set_attr("think-time", detail::lqnx_num(tk.thinktime.disabled ? 0.0 : think));
            } else if (!tk.thinktime.disabled && think > 0.0) {
                report.dropped.push_back("task '" + tk.name + "' has a think time of " +
                                         detail::lqnx_num(think) +
                                         ", which the schema accepts on reference tasks only");
            }
            // <setup>/<delay-off> are a LINE extension the reference writes
            // ahead of fan-out (writeXML.m:170-183), and every reader in the
            // project takes them; refusing them here made the C++ row the only
            // one that could not round-trip a SetupTask.
            if (!tk.setuptime.disabled &&
                num_traits<T>::to_double(tk.setuptime.mean) > lang::GlobalConstants::FineTol) {
                xml::Element& se = te.add_child("setup");
                se.set_attr("mean", detail::lqnx_num(num_traits<T>::to_double(tk.setuptime.mean)));
                se.set_attr("scv", detail::lqnx_num(num_traits<T>::to_double(tk.setuptime.scv)));
            }
            if (!tk.delayofftime.disabled &&
                num_traits<T>::to_double(tk.delayofftime.mean) > lang::GlobalConstants::FineTol) {
                xml::Element& de = te.add_child("delay-off");
                de.set_attr("mean",
                            detail::lqnx_num(num_traits<T>::to_double(tk.delayofftime.mean)));
                de.set_attr("scv",
                            detail::lqnx_num(num_traits<T>::to_double(tk.delayofftime.scv)));
            }
            // lqn-core.xsd (TaskType) places fan-out and fan-in before the entries.
            for (std::size_t f = 0; f < tk.fanout.size(); ++f) {
                xml::Element& fe = te.add_child("fan-out");
                fe.set_attr("dest", task_name(tk.fanout[f].first));
                fe.set_attr("value", detail::lqnx_num(tk.fanout[f].second));
            }
            for (std::size_t f = 0; f < tk.fanin.size(); ++f) {
                xml::Element& fe = te.add_child("fan-in");
                fe.set_attr("source", task_name(tk.fanin[f].first));
                fe.set_attr("value", detail::lqnx_num(tk.fanin[f].second));
            }

            for (std::size_t e = 0; e < m.entries.size(); ++e) {
                if (m.entries[e].task_slot != t) continue;
                const detail::RawEntry<T>& en = m.entries[e];
                xml::Element& ee = te.add_child("entry");
                ee.set_attr("name", entry_name(en.name));
                ee.set_attr("type", "NONE");
                if (en.has_arrival && !en.arrival.disabled) {
                    const double mean = num_traits<T>::to_double(en.arrival.mean);
                    if (std::isfinite(mean) && mean > lang::GlobalConstants::FineTol)
                        ee.set_attr("open-arrival-rate", detail::lqnx_num(1.0 / mean));
                }
                for (std::size_t f = 0; f < en.fwd_dest.size(); ++f) {
                    xml::Element& fe = ee.add_child("forwarding");
                    fe.set_attr("dest", entry_name(en.fwd_dest[f]));
                    fe.set_attr("prob", detail::lqnx_num_of<T>(en.fwd_prob[f]));
                }
            }

            xml::Element& ta = te.add_child("task-activities");
            for (std::size_t a = 0; a < m.acts.size(); ++a) {
                if (m.acts[a].task_slot != t) continue;
                const detail::RawActivity<T>& ac = m.acts[a];
                xml::Element& ae = ta.add_child("activity");
                ae.set_attr("host-demand-mean",
                            detail::lqnx_num(ac.hostdem.disabled
                                                 ? 0.0
                                                 : num_traits<T>::to_double(ac.hostdem.mean)));
                ae.set_attr("host-demand-cvsq",
                            detail::lqnx_num(ac.hostdem.disabled
                                                 ? 1.0
                                                 : num_traits<T>::to_double(ac.hostdem.scv)));
                if (!ac.bound_to_entry.empty())
                    ae.set_attr("bound-to-entry", entry_name(ac.bound_to_entry));
                ae.set_attr("call-order", "STOCHASTIC");
                ae.set_attr("name", act_name(ac.name));
                const double athink = num_traits<T>::to_double(ac.thinktime.mean);
                if (!ac.thinktime.disabled && athink > lang::GlobalConstants::FineTol)
                    ae.set_attr("think-time", detail::lqnx_num(athink));
                for (std::size_t c = 0; c < ac.sync_calls.size(); ++c) {
                    xml::Element& ce = ae.add_child("synch-call");
                    ce.set_attr("dest", entry_name(ac.sync_calls[c].dest));
                    ce.set_attr("calls-mean", detail::lqnx_num_of<T>(ac.sync_calls[c].mean));
                }
                for (std::size_t c = 0; c < ac.async_calls.size(); ++c) {
                    xml::Element& ce = ae.add_child("asynch-call");
                    ce.set_attr("dest", entry_name(ac.async_calls[c].dest));
                    ce.set_attr("calls-mean", detail::lqnx_num_of<T>(ac.async_calls[c].mean));
                }
                // LINE dialect: which of the synch-calls above one dispatcher
                // issues, and under which strategy. The member calls stay
                // ordinary synch-calls, so a reader that ignores this element
                // still sees the same aggregate call means -- which is what
                // lqns and lqsim, having no dispatcher, should see.
                for (std::size_t g = 0; g < ac.call_groups.size(); ++g) {
                    xml::Element& ge = ae.add_child("call-group");
                    ge.set_attr("strategy", detail::callgroup_to_lqnx(ac.call_groups[g].first));
                    for (std::size_t d = 0; d < ac.call_groups[g].second.size(); ++d)
                        ge.add_child("dest").set_attr("name",
                                                      entry_name(ac.call_groups[g].second[d]));
                }
            }

            for (std::size_t q = 0; q < tk.precedences.size(); ++q) {
                const detail::RawPrecedence<T>& pc = tk.precedences[q];
                xml::Element& pce = ta.add_child("precedence");

                xml::Element& pre = pce.add_child(detail::precedence_tag(pc.pretype));
                if (pc.pretype == PrecedenceType::PRE_AND && pc.has_quorum)
                    pre.set_attr("quorum", detail::lqnx_num(static_cast<double>(pc.quorum)));
                for (std::size_t i = 0; i < pc.preacts.size(); ++i)
                    pre.add_child("activity").set_attr("name", act_name(pc.preacts[i]));

                xml::Element& post = pce.add_child(detail::precedence_tag(pc.posttype));
                if (pc.posttype == PrecedenceType::POST_OR) {
                    for (std::size_t i = 0; i < pc.postacts.size(); ++i) {
                        xml::Element& ae = post.add_child("activity");
                        ae.set_attr("name", act_name(pc.postacts[i]));
                        if (i < pc.postparams.size())
                            ae.set_attr("prob", detail::lqnx_num_of<T>(pc.postparams[i]));
                    }
                } else if (pc.posttype == PrecedenceType::POST_LOOP) {
                    // The LAST post activity is the loop exit and is named by
                    // the `end` attribute, not by an <activity> of its own --
                    // which is also how the reader takes it apart again.
                    if (pc.postacts.empty())
                        throw InputError("write_lqnx: a post-LOOP names no activity");
                    for (std::size_t i = 0; i + 1 < pc.postacts.size(); ++i) {
                        xml::Element& ae = post.add_child("activity");
                        ae.set_attr("name", act_name(pc.postacts[i]));
                        if (i < pc.postparams.size())
                            ae.set_attr("count", detail::lqnx_num_of<T>(pc.postparams[i]));
                    }
                    post.set_attr("end", act_name(pc.postacts.back()));
                } else if (pc.posttype == PrecedenceType::POST_CACHE) {
                    // NAME THE BRANCH, do not leave it to position: hit first and
                    // miss second is the builder's order, but a reader that sorts
                    // or a writer that reorders would otherwise swap them
                    // silently. writeXML.m:305-317 sets the same attribute.
                    static const char* const kResult[2] = {"hit", "miss"};
                    for (std::size_t i = 0; i < pc.postacts.size(); ++i) {
                        xml::Element& ae = post.add_child("activity");
                        ae.set_attr("name", act_name(pc.postacts[i]));
                        if (i < 2) ae.set_attr("cache-result", kResult[i]);
                    }
                } else {
                    for (std::size_t i = 0; i < pc.postacts.size(); ++i)
                        post.add_child("activity").set_attr("name", act_name(pc.postacts[i]));
                }
            }

            if (tk.sched != SchedStrategy::REF) {
                for (std::size_t e = 0; e < m.entries.size(); ++e) {
                    if (m.entries[e].task_slot != t) continue;
                    const std::string& ename = m.entries[e].name;
                    if (!needs_reply.count(ename)) continue;
                    const std::map<std::string, std::vector<std::string>>::const_iterator it =
                        replies.find(ename);
                    if (it == replies.end() || it->second.empty())
                        throw InputError(
                            "write_lqnx: entry '" + ename +
                            "' is called synchronously but no activity replies to it, and none "
                            "can be inferred; declare one with replies_to()");
                    xml::Element& re = ta.add_child("reply-entry");
                    re.set_attr("name", entry_name(ename));
                    for (std::size_t r = 0; r < it->second.size(); ++r)
                        re.add_child("reply-activity").set_attr("name", act_name(it->second[r]));
                }
            }
        }
    }

    xml::write_file(path, root);
    return report;
}

}  // namespace lqn
}  // namespace line

#endif  // LINE_LANG_LQN_LQN_WRITER_H
