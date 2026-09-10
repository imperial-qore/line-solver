/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_LANG_QN_FEATURE_SET_H
#define LINE_LANG_QN_FEATURE_SET_H

/**
 * The language-feature gate: what a MODEL uses against what a SOLVER declares.
 *
 * SCOPE. This is the port of three MATLAB pieces that together form the single
 * universal check every `runAnalyzer` performs before it computes anything:
 *
 *   matlab/src/solvers/SolverFeatureSet.m          the canonical name registry
 *                                                  and the SUPPORTS comparison
 *   matlab/src/lang/@@MNetwork/getUsedLangFeatures.m   the USED side
 *   matlab/src/solvers/@@NetworkSolver/NetworkSolver.m:169-216
 *                                                  runAnalyzerChecks, which
 *                                                  raises line_error naming the
 *                                                  offending feature
 *
 * WHY IT EXISTS. Without it, every refusal is ad hoc and per solver, and the
 * constructs nobody wrote a refusal for are neither handled nor rejected: the
 * solver returns numbers for a model the user did not describe. A finite
 * capacity region silently ignored is not a small error, it is a different
 * model. A generic "unsupported model" is nearly as bad, because the user
 * cannot tell which construct to remove, so `supports` returns the offending
 * FEATURES and a reason naming them.
 *
 * THE ASYMMETRY THAT MATTERS. A feature DECLARED but not implemented is worse
 * than one left undeclared: the undeclared one produces a clean refusal, the
 * over-declared one produces a wrong number. So a C++ solver's declared set is
 * what its C++ code actually handles, never a transcription of the MATLAB set.
 *
 * WHAT THE REGISTRY IS FOR. `Feature` lists every name in
 * SolverFeatureSet.fields, including names this port can never emit. That is
 * deliberate and is the reference's own convention: SUPPORTS iterates the
 * registry, so a capability name MISSING from the registry is invisible to the
 * gate and passes as if the capability were absent. Keeping the registry
 * complete keeps that failure mode out.
 *
 * WHAT NetworkStruct CANNOT EXPRESS. `used_lang_features` derives each feature
 * from a field that exists. Where MATLAB has a feature and NetworkStruct has no
 * field for it, the feature is NEVER emitted, and that is recorded by name at
 * the point where it would have been derived rather than faked from a proxy.
 * The list is SetupDelayOff, Reneging, Balking, QueueingPlace and
 * CacheRetrieval's MATLAB-only cousins; see the comments in
 * used_lang_features. Breakdown LEFT that list on 2026-08-15, when
 * `NetworkStruct::breakdownparam` gave it a field to be derived from.
 */

#include <cmath>
#include <cstddef>
#include <algorithm>
#include <initializer_list>
#include <sstream>
#include <string>
#include <vector>

#include "line/lang/lang_types.h"
#include "line/lang/qn/network_struct.h"
#include "line/util/error.h"

namespace line {
namespace qn {

// ---------------------------------------------------------------------------
// The canonical feature registry
// ---------------------------------------------------------------------------

/**
 * SolverFeatureSet.fields, in its order, as an X-list.
 *
 * The order is the reference's and is kept so that the two lists can be
 * diffed entry by entry; `supports` reports in this order too, which makes the
 * reason string reproducible.
 */
#define LINE_QN_FEATURE_LIST(X)     \
    X(ClassSwitch)                  \
    X(Cache)                        \
    X(Delay)                        \
    X(DelayStation)                 \
    X(Fork)                         \
    X(Join)                         \
    X(Logger)                       \
    X(Place)                        \
    X(QueueingPlace)                \
    X(Queue)                        \
    X(JobSink)                      \
    X(Sink)                         \
    X(Source)                       \
    X(Router)                       \
    X(Transition)                   \
    X(Coxian)                       \
    X(Cox2)                         \
    X(APH)                          \
    X(Det)                          \
    X(Disabled)                     \
    X(Erlang)                       \
    X(Exp)                          \
    X(Gamma)                        \
    X(HyperExp)                     \
    X(Immediate)                    \
    X(Lognormal)                    \
    X(MAP)                          \
    X(DMAP)                         \
    X(MMAP)                         \
    X(BMAP)                         \
    X(MMPP2)                        \
    X(NHPP)                         \
    X(MAPt)                         \
    X(PHt)                          \
    X(EmpiricalCdf)                 \
    X(Expolynomial)                 \
    X(Normal)                       \
    X(Pareto)                       \
    X(PH)                           \
    X(ME)                           \
    X(RAP)                          \
    X(Replayer)                     \
    X(Trace)                        \
    X(Uniform)                      \
    X(Weibull)                      \
    X(Bernoulli)                    \
    X(Binomial)                     \
    X(Geometric)                    \
    X(Poisson)                      \
    X(DiscreteSampler)              \
    X(DiscreteUniform)              \
    X(Empirical)                    \
    X(GMM)                          \
    X(MMDP)                         \
    X(MMDP2)                        \
    X(MMPP)                         \
    X(MultivariateNormal)           \
    X(NegBinomial)                  \
    X(Prior)                        \
    X(Zipf)                         \
    X(StatelessClassSwitcher)       \
    X(CacheClassSwitcher)           \
    X(CacheRetrieval)               \
    X(CacheItemSize)                \
    X(InfiniteServer)               \
    X(Forker)                       \
    X(Joiner)                       \
    X(LogTunnel)                    \
    X(SharedServer)                 \
    X(Buffer)                       \
    X(Region)                       \
    X(Linkage)                      \
    X(Enabling)                     \
    X(Inhibiting)                   \
    X(Timing)                       \
    X(Firing)                       \
    X(Storage)                      \
    X(RandomSource)                 \
    X(Dispatcher)                   \
    X(Server)                       \
    X(ServiceTunnel)                \
    X(RoutingStrategy_PROB)         \
    X(RoutingStrategy_RAND)         \
    X(RoutingStrategy_RROBIN)       \
    X(RoutingStrategy_WRROBIN)      \
    X(RoutingStrategy_JSQ)          \
    X(RoutingStrategy_SQ)           \
    X(RoutingStrategy_SDR)          \
    X(SchedStrategy_INF)            \
    X(SchedStrategy_FCFS)           \
    X(SchedStrategy_FCFSPR)         \
    X(SchedStrategy_FCFSPI)         \
    X(SchedStrategy_FCFSPRIO)       \
    X(SchedStrategy_FCFSPRPRIO)     \
    X(SchedStrategy_FCFSPIPRIO)     \
    X(SchedStrategy_LCFS)           \
    X(SchedStrategy_LCFSPR)         \
    X(SchedStrategy_LCFSPI)         \
    X(SchedStrategy_LCFSPRIO)       \
    X(SchedStrategy_LCFSPRPRIO)     \
    X(SchedStrategy_LCFSPIPRIO)     \
    X(SchedStrategy_SEPT)           \
    X(SchedStrategy_LEPT)           \
    X(SchedStrategy_SJF)            \
    X(SchedStrategy_LJF)            \
    X(SchedStrategy_SRPT)           \
    X(SchedStrategy_SRPTPRIO)       \
    X(SchedStrategy_PSJF)           \
    X(SchedStrategy_FB)             \
    X(SchedStrategy_LRPT)           \
    X(SchedStrategy_SETF)           \
    X(SchedStrategy_FSP)            \
    X(SchedStrategy_PAS)            \
    X(SchedStrategy_OI)             \
    X(SchedStrategy_PS)             \
    X(SchedStrategy_DPS)            \
    X(SchedStrategy_GPS)            \
    X(SchedStrategy_PSPRIO)         \
    X(SchedStrategy_DPSPRIO)        \
    X(SchedStrategy_GPSPRIO)        \
    X(SchedStrategy_SIRO)           \
    X(SchedStrategy_HOL)            \
    X(SchedStrategy_EXT)            \
    X(SchedStrategy_POLLING)        \
    X(SchedStrategy_EDD)            \
    X(SchedStrategy_EDF)            \
    X(SchedStrategy_LPS)            \
    X(ReplacementStrategy_RR)       \
    X(ReplacementStrategy_FIFO)     \
    X(ReplacementStrategy_SFIFO)    \
    X(ReplacementStrategy_LRU)      \
    X(ReplacementStrategy_HLRU)     \
    X(ReplacementStrategy_CLIMB)    \
    X(ReplacementStrategy_QLRU)     \
    X(ClosedClass)                  \
    X(OpenClass)                    \
    X(SelfLoopingClass)             \
    X(OpenSignal)                   \
    X(ClosedSignal)                 \
    X(SignalType_NEGATIVE)          \
    X(SignalType_REPLY)             \
    X(SignalType_CATASTROPHE)       \
    X(SignalBatchRemoval)           \
    X(SignalRemovalPolicy)          \
    X(BatchArrival)                 \
    X(LoadDependence)               \
    X(ClassDependence)              \
    X(JointDependence)              \
    X(GlobalDependence)             \
    X(SetupDelayOff)                \
    X(ServerParallelism)            \
    X(Retrial)                      \
    X(Balking)                      \
    X(Reneging)                     \
    X(Breakdown)                    \
    X(HeteroServers)                \
    X(DepartureDiscipline)          \
    X(JoinPartial)                  \
    X(Host)                         \
    X(Processor)                    \
    X(Task)                         \
    X(Entry)                        \
    X(Activity)                     \
    X(SyncCall)                     \
    X(AsyncCall)                    \
    X(ActivityPrecedence_PRE_SEQ)   \
    X(ActivityPrecedence_POST_SEQ)  \
    X(ActivityPrecedence_PRE_AND)   \
    X(ActivityPrecedence_POST_AND)  \
    X(ActivityPrecedence_PRE_OR)    \
    X(ActivityPrecedence_POST_OR)   \
    X(SchedStrategy_REF)            \
    /* Layered cache-queueing models: a CacheTask holds a segmented cache whose
       items are ItemEntry entries, and a read is a call to one of them whose
       bound activity branches on a POST_CACHE precedence into hit and miss.
       Registered in all four codebases because SolverLDES.getLNFeatureSet
       DECLARES all three, and an unregistered name makes setTrue line_error --
       so that whole feature set threw before it could compare anything. */ \
    X(CacheTask)                    \
    X(ItemEntry)                    \
    X(ActivityPrecedence_POST_CACHE) \
    /* Variable forking levels. A Fork emits tasksPerLink jobs on every
       outgoing link; these three name the ways that degree stops being one
       number. ForkFanoutVector: the count differs by destination or by class
       (sn.nodeparam{f}.fanOutLink). ForkFanoutRandom: the count is a draw from
       a DiscreteSampler, redrawn per link and per forked job
       (sn.nodeparam{f}.fanOutDist). ForkBranchProbability: a branch fires only
       with probability p, so the SIBLING COUNT is random even when each link
       carries a fixed number (sn.nodeparam{f}.fanOutProb). Appended at the tail
       so every earlier index is unchanged. */ \
    X(ForkFanoutVector)             \
    X(ForkFanoutRandom)             \
    X(ForkBranchProbability)        \
    /* Ported from MATLAB SolverFeatureSet.m (2026-09-05), in this order, so
       that every earlier enumerator is unchanged. Both are "having something"
       properties the registry could not name, so every rule about them lived
       only in structural predicates and was invisible to the gate.
       MultiServer: a finite-server station serving more than one job at once,
       i.e. a station whose nservers is finite and > 1 (a Delay / INF station is
       not one), which `NetworkStruct::has_multi_server` already tests.
       FiniteCapacity: a station or per-class buffer that can BIND, exactly the
       condition `has_binding_capacity` tests (station cap / classcap below the
       population that can reach it; an open class always binds; a Cache model
       is exempt), which is also what `check_binding_capacity` refuses on. */ \
    X(MultiServer)                  \
    X(FiniteCapacity)

/** One language feature. COUNT is the registry size and is not a feature. */
enum class Feature : int {
#define LINE_QN_FEATURE_ENUM(id) id,
    LINE_QN_FEATURE_LIST(LINE_QN_FEATURE_ENUM)
#undef LINE_QN_FEATURE_ENUM
        COUNT
};

/** The number of registered features, MATLAB `numel(SolverFeatureSet.fields)`. */
inline constexpr std::size_t feature_count() { return static_cast<std::size_t>(Feature::COUNT); }

/** The canonical registry name, byte for byte the MATLAB field name. */
inline const char* feature_name(Feature f) {
    switch (f) {
#define LINE_QN_FEATURE_NAME(id) \
    case Feature::id:            \
        return #id;
        LINE_QN_FEATURE_LIST(LINE_QN_FEATURE_NAME)
#undef LINE_QN_FEATURE_NAME
        default:
            return "Unknown";
    }
}

/**
 * A human-readable phrase for the reason string, e.g. "a Finite Capacity
 * Region". Structured names are derived from their suffix so that the whole
 * SchedStrategy_ / RoutingStrategy_ / ReplacementStrategy_ / SignalType_ /
 * ActivityPrecedence_ families need no per-entry table to stay in step.
 */
inline std::string feature_phrase(Feature f) {
    switch (f) {
        case Feature::Region: return "a Finite Capacity Region";
        case Feature::Fork: return "a Fork node";
        case Feature::Join: return "a Join node";
        case Feature::Forker: return "a fork output section";
        case Feature::ForkFanoutVector:
            return "a fork whose tasks per link differ by destination or by class";
        case Feature::ForkFanoutRandom:
            return "a fork whose tasks per link are drawn from a distribution";
        case Feature::ForkBranchProbability:
            return "a fork whose branches fire only with a probability";
        case Feature::Joiner: return "a join input section";
        case Feature::Cache: return "a Cache node";
        case Feature::CacheRetrieval: return "a delayed-hit cache retrieval system";
        case Feature::CacheItemSize: return "per-item storage costs with per-list cost caps";
        case Feature::ClassSwitch: return "a ClassSwitch node";
        case Feature::Router: return "a Router node";
        case Feature::Logger: return "a Logger node";
        case Feature::Source: return "a Source node";
        case Feature::Sink: return "a Sink node";
        case Feature::Place: return "an SPN Place";
        case Feature::QueueingPlace: return "a queueing SPN Place";
        case Feature::Transition: return "an SPN Transition";
        case Feature::Inhibiting: return "an SPN inhibitor arc";
        case Feature::Retrial: return "a retrial orbit";
        case Feature::Balking: return "job balking";
        case Feature::Reneging: return "job reneging";
        case Feature::Breakdown: return "server breakdowns";
        case Feature::SetupDelayOff: return "server setup and delay-off times";
        case Feature::LoadDependence: return "load-dependent service rates";
        case Feature::ClassDependence: return "class-dependent service rates";
        case Feature::JointDependence: return "joint-dependent service rates";
        case Feature::GlobalDependence: return "globally state-dependent (Whittle) service rates";
        case Feature::BatchArrival: return "batch arrivals";
        case Feature::OpenSignal: return "an open G-network signal class";
        case Feature::ClosedSignal: return "a closed G-network signal class";
        case Feature::SignalBatchRemoval: return "a signal batch-removal distribution";
        case Feature::SignalRemovalPolicy: return "a non-random signal removal policy";
        case Feature::OpenClass: return "an open job class";
        case Feature::ClosedClass: return "a closed job class";
        case Feature::SelfLoopingClass: return "a self-looping job class";
        case Feature::Queue: return "a Queue station";
        case Feature::Delay:
        case Feature::DelayStation: return "a Delay station";
        // "multi-server" is the vocabulary every other refusal in this tree uses
        // (jmt_method_refusal, ba, ssa nrm), and the gate reports must agree with them.
        case Feature::MultiServer: return "a multi-server station";
        case Feature::FiniteCapacity: return "a finite station or per-class buffer that binds";
        default: break;
    }
    const std::string nm = feature_name(f);
    const std::size_t us = nm.find('_');
    if (us != std::string::npos) {
        const std::string head = nm.substr(0, us), tail = nm.substr(us + 1);
        if (head == "SchedStrategy") return "the " + tail + " scheduling discipline";
        if (head == "RoutingStrategy") return "the " + tail + " routing strategy";
        if (head == "ReplacementStrategy") return "the " + tail + " cache replacement policy";
        if (head == "SignalType") return tail + " signals";
        if (head == "ActivityPrecedence") return "the " + tail + " activity precedence";
    }
    return "the " + nm + " feature";
}

// ---------------------------------------------------------------------------
// The set
// ---------------------------------------------------------------------------

/** A subset of the registry: MATLAB's SolverFeatureSet, whose `list` is a flag per field. */
class FeatureSet {
public:
    FeatureSet() : bits_(feature_count(), false) {}

    FeatureSet& set(Feature f) {
        if (f != Feature::COUNT) bits_[static_cast<std::size_t>(f)] = true;
        return *this;
    }
    FeatureSet& set(std::initializer_list<Feature> fs) {
        for (Feature f : fs) set(f);
        return *this;
    }
    FeatureSet& unset(Feature f) {
        if (f != Feature::COUNT) bits_[static_cast<std::size_t>(f)] = false;
        return *this;
    }
    bool has(Feature f) const {
        return f != Feature::COUNT && bits_[static_cast<std::size_t>(f)];
    }
    /** The features held, in registry order. */
    std::vector<Feature> list() const {
        std::vector<Feature> out;
        for (std::size_t i = 0; i < bits_.size(); ++i)
            if (bits_[i]) out.push_back(static_cast<Feature>(i));
        return out;
    }
    bool empty() const {
        for (std::size_t i = 0; i < bits_.size(); ++i)
            if (bits_[i]) return false;
        return true;
    }

private:
    std::vector<bool> bits_;
};

/**
 * The registry name a specialization falls back to when it is not declared.
 *
 * Some registry names denote a SPECIAL CASE of another name rather than a
 * capability of their own: a Cox2 is a Coxian restricted to two phases, and a
 * Trace is a Replayer under another class name. Recording only the general name
 * left the specific entry unreachable -- dead registry surface that no model
 * could ever set (see _kb/06-solver-catalog.md, "A registered name nothing
 * emits gates nothing"). Recording the specific name
 * instead would silently REJECT those models at every solver that declares only
 * the general one, which is every solver that accepts them today.
 *
 * So the recorder emits the most specific name it can, and the gate resolves an
 * undeclared specific name against its generalization here. A solver that
 * genuinely supports only the special case (two-phase Coxian, say) keeps the
 * option of declaring the specialization alone: the fallback is consulted only
 * when the specific name is missing, never in the other direction.
 *
 * COUNT means "no generalization"; the feature stands on its own.
 */
inline Feature feature_generalization(Feature f) {
    switch (f) {
        case Feature::Cox2: return Feature::Coxian;
        // Trace is emitted by MATLAB, the JAR and python, which each carry a
        // Trace class distinct from Replayer. This port has no such type -- the
        // JSON wire writes both as "Replayer" -- so the entry stays inert here
        // and the fallback is kept only so the four tables read alike.
        case Feature::Trace: return Feature::Replayer;
        default: return Feature::COUNT;
    }
}

/** The outcome of the gate: the verdict, the offending features, the message. */
struct SupportResult {
    bool ok = true;
    std::vector<Feature> missing;  ///< used but not declared, in registry order
    std::string reason;            ///< empty when ok
};

/**
 * SolverFeatureSet.supports: is every feature the model uses declared?
 *
 * `solver` names the solver in the reason so the message is actionable on its
 * own; MATLAB gets that name from the mfilename of the raising runAnalyzer.
 */
inline SupportResult feature_set_supports(const std::string& solver, const FeatureSet& declared,
                                          const FeatureSet& used) {
    SupportResult r;
    for (std::size_t i = 0; i < feature_count(); ++i) {
        const Feature f = static_cast<Feature>(i);
        if (!used.has(f) || declared.has(f)) continue;
        // A specialization the solver did not name is covered by the general
        // capability when that one IS declared; see feature_generalization.
        const Feature g = feature_generalization(f);
        if (g != Feature::COUNT && declared.has(g)) continue;
        r.missing.push_back(f);
    }
    if (r.missing.empty()) return r;
    r.ok = false;
    if (r.missing.size() == 1) {
        r.reason = solver + ": this model uses " + feature_phrase(r.missing[0]) +
                   ", which this solver does not support";
    } else {
        r.reason = solver + ": this model uses features which this solver does not support: ";
        for (std::size_t i = 0; i < r.missing.size(); ++i) {
            if (i) r.reason += ", ";
            r.reason += feature_phrase(r.missing[i]);
        }
    }
    return r;
}

// ---------------------------------------------------------------------------
// Enum to feature
// ---------------------------------------------------------------------------

/**
 * The feature name of a distribution, MATLAB `serviceProcess{r}{3}.name`.
 *
 * Immediate and Disabled return COUNT: they are internal placeholders, not
 * user-facing distributions, and getUsedLangFeatures excludes them by name.
 * Expolynomial and Normal have no ProcessType in MATLAB either, so those two
 * registry entries can never be emitted by this port.
 */
inline Feature feature_of_process(ProcessType p) {
    switch (p) {
        case ProcessType::EXP: return Feature::Exp;
        case ProcessType::ERLANG: return Feature::Erlang;
        case ProcessType::HYPEREXP: return Feature::HyperExp;
        case ProcessType::PH: return Feature::PH;
        case ProcessType::APH: return Feature::APH;
        case ProcessType::MAP: return Feature::MAP;
        case ProcessType::UNIFORM: return Feature::Uniform;
        case ProcessType::DET: return Feature::Det;
        case ProcessType::COXIAN: return Feature::Coxian;
        case ProcessType::GAMMA: return Feature::Gamma;
        case ProcessType::PARETO: return Feature::Pareto;
        case ProcessType::MMPP2: return Feature::MMPP2;
        case ProcessType::REPLAYER: return Feature::Replayer;
        case ProcessType::COX2: return Feature::Cox2;
        case ProcessType::WEIBULL: return Feature::Weibull;
        case ProcessType::LOGNORMAL: return Feature::Lognormal;
        case ProcessType::DUNIFORM: return Feature::DiscreteUniform;
        case ProcessType::BERNOULLI: return Feature::Bernoulli;
        case ProcessType::BINOMIAL: return Feature::Binomial;
        case ProcessType::POISSON: return Feature::Poisson;
        case ProcessType::GEOMETRIC: return Feature::Geometric;
        case ProcessType::BMAP: return Feature::BMAP;
        case ProcessType::ME: return Feature::ME;
        case ProcessType::RAP: return Feature::RAP;
        case ProcessType::DISCRETESAMPLER: return Feature::DiscreteSampler;
        case ProcessType::ZIPF: return Feature::Zipf;
        case ProcessType::DMAP: return Feature::DMAP;
        case ProcessType::MMAP: return Feature::MMAP;
        case ProcessType::EMPIRICALCDF: return Feature::EmpiricalCdf;
        // A Prior IS user-facing, unlike Immediate and Disabled, and emitting it
        // is what makes every solver but SolverUQ refuse the model instead of
        // lowering the epistemic mixture to a rate and answering with it.
        case ProcessType::PRIOR: return Feature::Prior;
        case ProcessType::NHPP: return Feature::NHPP;
        case ProcessType::MAPT: return Feature::MAPt;
        case ProcessType::PHT: return Feature::PHt;
        default: return Feature::COUNT;
    }
}

/**
 * SchedStrategy.toFeature. FORK and NONE are internal markers with no registry
 * name; MATLAB aliases FCFSPRIO to HOL, so HOL covers both.
 */
inline Feature feature_of_sched(SchedStrategy s) {
    switch (s) {
        case SchedStrategy::INF: return Feature::SchedStrategy_INF;
        case SchedStrategy::FCFS: return Feature::SchedStrategy_FCFS;
        case SchedStrategy::FCFSPR: return Feature::SchedStrategy_FCFSPR;
        case SchedStrategy::FCFSPI: return Feature::SchedStrategy_FCFSPI;
        case SchedStrategy::FCFSPRPRIO: return Feature::SchedStrategy_FCFSPRPRIO;
        case SchedStrategy::FCFSPIPRIO: return Feature::SchedStrategy_FCFSPIPRIO;
        case SchedStrategy::LCFS: return Feature::SchedStrategy_LCFS;
        case SchedStrategy::LCFSPR: return Feature::SchedStrategy_LCFSPR;
        case SchedStrategy::LCFSPI: return Feature::SchedStrategy_LCFSPI;
        case SchedStrategy::LCFSPRIO: return Feature::SchedStrategy_LCFSPRIO;
        case SchedStrategy::LCFSPRPRIO: return Feature::SchedStrategy_LCFSPRPRIO;
        case SchedStrategy::LCFSPIPRIO: return Feature::SchedStrategy_LCFSPIPRIO;
        case SchedStrategy::SEPT: return Feature::SchedStrategy_SEPT;
        case SchedStrategy::LEPT: return Feature::SchedStrategy_LEPT;
        case SchedStrategy::SJF: return Feature::SchedStrategy_SJF;
        case SchedStrategy::LJF: return Feature::SchedStrategy_LJF;
        case SchedStrategy::SRPT: return Feature::SchedStrategy_SRPT;
        case SchedStrategy::SRPTPRIO: return Feature::SchedStrategy_SRPTPRIO;
        case SchedStrategy::PSJF: return Feature::SchedStrategy_PSJF;
        case SchedStrategy::FB: return Feature::SchedStrategy_FB;
        case SchedStrategy::LRPT: return Feature::SchedStrategy_LRPT;
        case SchedStrategy::SETF: return Feature::SchedStrategy_SETF;
        case SchedStrategy::FSP: return Feature::SchedStrategy_FSP;
        case SchedStrategy::PAS: return Feature::SchedStrategy_PAS;
        case SchedStrategy::OI: return Feature::SchedStrategy_OI;
        case SchedStrategy::PS: return Feature::SchedStrategy_PS;
        case SchedStrategy::DPS: return Feature::SchedStrategy_DPS;
        case SchedStrategy::GPS: return Feature::SchedStrategy_GPS;
        case SchedStrategy::PSPRIO: return Feature::SchedStrategy_PSPRIO;
        case SchedStrategy::DPSPRIO: return Feature::SchedStrategy_DPSPRIO;
        case SchedStrategy::GPSPRIO: return Feature::SchedStrategy_GPSPRIO;
        case SchedStrategy::SIRO: return Feature::SchedStrategy_SIRO;
        case SchedStrategy::HOL: return Feature::SchedStrategy_HOL;
        case SchedStrategy::EXT: return Feature::SchedStrategy_EXT;
        case SchedStrategy::POLLING: return Feature::SchedStrategy_POLLING;
        case SchedStrategy::EDD: return Feature::SchedStrategy_EDD;
        case SchedStrategy::EDF: return Feature::SchedStrategy_EDF;
        case SchedStrategy::LPS: return Feature::SchedStrategy_LPS;
        case SchedStrategy::REF: return Feature::SchedStrategy_REF;
        default: return Feature::COUNT;
    }
}

/**
 * RoutingStrategy.toFeature. FIRING and DISABLED are internal markers, not
 * user-selectable capabilities, and map to no registry name (RoutingStrategy.m:132-139).
 */
inline Feature feature_of_routing(RoutingStrategy r) {
    switch (r) {
        case RoutingStrategy::PROB: return Feature::RoutingStrategy_PROB;
        case RoutingStrategy::RAND: return Feature::RoutingStrategy_RAND;
        case RoutingStrategy::RROBIN: return Feature::RoutingStrategy_RROBIN;
        case RoutingStrategy::WRROBIN: return Feature::RoutingStrategy_WRROBIN;
        case RoutingStrategy::JSQ: return Feature::RoutingStrategy_JSQ;
        case RoutingStrategy::SQ: return Feature::RoutingStrategy_SQ;
        case RoutingStrategy::SDR: return Feature::RoutingStrategy_SDR;
        default: return Feature::COUNT;
    }
}

/** ReplacementStrategy.toFeature. */
inline Feature feature_of_replacement(ReplacementStrategy s) {
    switch (s) {
        case ReplacementStrategy::RR: return Feature::ReplacementStrategy_RR;
        case ReplacementStrategy::FIFO: return Feature::ReplacementStrategy_FIFO;
        case ReplacementStrategy::SFIFO: return Feature::ReplacementStrategy_SFIFO;
        case ReplacementStrategy::LRU: return Feature::ReplacementStrategy_LRU;
        case ReplacementStrategy::HLRU: return Feature::ReplacementStrategy_HLRU;
        case ReplacementStrategy::CLIMB: return Feature::ReplacementStrategy_CLIMB;
        case ReplacementStrategy::QLRU: return Feature::ReplacementStrategy_QLRU;
        default: return Feature::COUNT;
    }
}

/**
 * The SERVER section a Queue gets for a discipline, MATLAB Queue.m:79-108.
 *
 * The section class name is itself a feature: MNetwork.addNode registers
 * class(node.server) alongside class(node), which is where SharedServer,
 * InfiniteServer and Server enter the used set. PreemptiveServer and
 * PollingServer are absent from the registry, so those disciplines contribute
 * only their SchedStrategy_ name, exactly as in the reference.
 */
inline Feature feature_of_server_section(SchedStrategy s) {
    switch (s) {
        case SchedStrategy::PS:
        case SchedStrategy::DPS:
        case SchedStrategy::GPS:
        case SchedStrategy::PSPRIO:
        case SchedStrategy::DPSPRIO:
        case SchedStrategy::GPSPRIO:
        case SchedStrategy::LPS: return Feature::SharedServer;
        case SchedStrategy::INF: return Feature::InfiniteServer;
        case SchedStrategy::LCFSPR:
        case SchedStrategy::LCFSPRPRIO:
        case SchedStrategy::LCFSPI:
        case SchedStrategy::LCFSPIPRIO:
        case SchedStrategy::FCFSPR:
        case SchedStrategy::FCFSPRPRIO:
        case SchedStrategy::FCFSPI:
        case SchedStrategy::FCFSPIPRIO:
        case SchedStrategy::EDF: return Feature::COUNT;  // PreemptiveServer, unregistered
        case SchedStrategy::POLLING: return Feature::COUNT;  // PollingServer, unregistered
        default: return Feature::Server;
    }
}

// ---------------------------------------------------------------------------
// The used side
// ---------------------------------------------------------------------------

/**
 * getUsedLangFeatures: the features the MODEL uses.
 *
 * MATLAB derives them from the Network OBJECT, where a node still knows its
 * MATLAB class and its section classes; this port derives them from the
 * refreshed struct, which is the only thing a C++ solver runner holds. The two
 * agree because the section classes are a pure function of the node kind and
 * the discipline (Queue.m:48-108, Source.m:66-68, Cache.m:54-76, Fork.m:106-109,
 * Join.m:33-35, Place.m:27-31, Transition.m:29-36), which is what
 * feature_of_server_section and the node switch below encode.
 *
 * FEATURES THIS PORT CANNOT EMIT, because NetworkStruct carries no field for
 * them. Each is a real MATLAB feature, and a model using one reaches a C++
 * solver unflagged; the fix is a struct field, not a proxy here.
 *   SetupDelayOff    MATLAB Queue.setupTime / delayoffTime, carried by
 *                    NetworkStruct::setupparam and solved by solver_mam_basic.
 *   Reneging         MATLAB Queue.impatienceTypes. state_events.h takes the
 *                    reneging classes as a CALL PARAMETER, so the struct does
 *                    not record them and they cannot be derived here.
 *   Balking          MATLAB Queue.balkingStrategies. No struct field.
 *   QueueingPlace    MATLAB Place.queueing. add_place() builds every Place as
 *                    an INF station with no queueing flag.
 * Breakdown was on this list until 2026-08-15, when
 * `NetworkStruct::breakdownparam` was added for the LDES engine; it is derived
 * below like any other field-backed feature. `state_events.h` still takes the
 * breakdown nodes as a call parameter, which is why the CTMC declaration of it
 * stays inert -- the two are independent.
 *
 * FEATURES EMITTED HERE THAT MATLAB DOES NOT EMIT. LoadDependence is in the
 * registry and is declared by solvers, but getUsedLangFeatures never set it, so
 * in MATLAB it was unenforced. Deriving it from `lldscaling` is a strict
 * tightening: the solvers that declare it are unaffected, and the ones that do
 * not now refuse instead of silently dropping the scaling. Retrial was in this
 * paragraph until 2026-09-05, when MATLAB's recorder gained the same arm.
 */
template <class T>
bool has_binding_capacity(const NetworkStruct<T>& sn);

template <class T>
FeatureSet used_lang_features(const NetworkStruct<T>& sn) {
    FeatureSet u;
    const std::size_t R = sn.classes.size();

    for (std::size_t r = 0; r < R; ++r) {
        u.set(sn.classes[r].type == JobClassType::CLOSED ? Feature::ClosedClass
                                                         : Feature::OpenClass);
        // A self-looping class is a closed class whose routing returns to its own
        // reference station. It was declared by eight solvers and never marked, so
        // a solver that omitted the name (the bounds, RCAT, QNS) was never asked;
        // it is now refused unless it declares the name.
        if (sn.classes[r].self_looping) u.set(Feature::SelfLoopingClass);
        if (r < sn.issignal.size() && sn.issignal[r]) {
            // Open vs closed signal follows the class type; MATLAB resolves the
            // unresolved Signal placeholder by the presence of a Source, which
            // is the same decision the refresh has already made here.
            u.set(sn.classes[r].type == JobClassType::CLOSED ? Feature::ClosedSignal
                                                             : Feature::OpenSignal);
            if (r < sn.signaltype.size()) {
                switch (sn.signaltype[r]) {
                    case lang::SignalType::NEGATIVE: u.set(Feature::SignalType_NEGATIVE); break;
                    case lang::SignalType::REPLY: u.set(Feature::SignalType_REPLY); break;
                    case lang::SignalType::CATASTROPHE:
                        u.set(Feature::SignalType_CATASTROPHE);
                        break;
                }
            }
            if (r < sn.signalremdist.size() && !sn.signalremdist[r].empty())
                u.set(Feature::SignalBatchRemoval);
            if (r < sn.signalrempolicy.size() &&
                sn.signalrempolicy[r] != lang::RemovalPolicy::RANDOM)
                u.set(Feature::SignalRemovalPolicy);
        }
    }

    // Globally state-dependent scaling phi(n) over the full network state, the
    // Whittle primitive. Only SolverCTMC plumbs it, so every other solver must
    // reject the model rather than solve it unscaled.
    if (static_cast<bool>(sn.gdscaling)) u.set(Feature::GlobalDependence);

    if (!sn.regions.empty()) u.set(Feature::Region);
    if (sn.sourceIdx != 0) u.set(Feature::Source);
    if (sn.sinkNode != 0) u.set(Feature::Sink);
    if (!sn.fj.empty()) u.set({Feature::Fork, Feature::Forker, Feature::Join, Feature::Joiner});

    for (std::size_t i = 0; i < sn.nodes.size(); ++i) {
        const NodeDef& nd = sn.nodes[i];
        const std::size_t ind = i + 1;
        const std::size_t ist = nd.station;  // 1-based, 0 when the node serves no jobs

        switch (nd.nodetype) {
            case NodeType::Queue:
            case NodeType::Delay: {
                u.set({nd.nodetype == NodeType::Delay ? Feature::Delay : Feature::Queue,
                       Feature::Buffer, Feature::Dispatcher});
                if (ist == 0) break;
                const Station<T>& st = sn.stations[ist - 1];
                u.set(feature_of_server_section(st.sched));
                bool any_class = false;
                for (std::size_t r = 0; r < R; ++r) {
                    // A disabled pair is MATLAB's empty serviceProcess: it emits
                    // no distribution, no discipline and no routing.
                    if (ist - 1 < sn.disabled.size() && r < sn.disabled[ist - 1].size() &&
                        sn.disabled[ist - 1][r])
                        continue;
                    any_class = true;
                    if (ist - 1 < sn.service.size() && r < sn.service[ist - 1].size())
                        u.set(feature_of_process(sn.service[ist - 1][r].type));
                    if (r < nd.routing.size()) u.set(feature_of_routing(nd.routing[r]));
                }
                if (any_class) u.set(feature_of_sched(st.sched));
                if (!st.lldscaling.empty()) u.set(Feature::LoadDependence);
                if (static_cast<bool>(st.cdscaling)) u.set(Feature::ClassDependence);
                if (static_cast<bool>(st.jdscaling)) u.set(Feature::JointDependence);
                // Impatience, balking and heterogeneous pools. MATLAB's
                // getUsedLangFeatures emits Balking and Reneging from the same
                // declarations; HeteroServers has no MATLAB counterpart and is
                // emitted here for the same reason Retrial and LoadDependence
                // are -- a solver that ignores the pools would answer for the
                // homogeneous station and say nothing about the difference.
                for (std::size_t r = 0; r < R; ++r) {
                    if (r < st.impatience.size() && st.impatience[r] != lang::ImpatienceType::NONE)
                        u.set(st.impatience[r] == lang::ImpatienceType::BALKING ? Feature::Balking
                                                                               : Feature::Reneging);
                    if (r < st.orbit_impatience.size() && !st.orbit_impatience[r].disabled)
                        u.set(Feature::Reneging);
                    if (r < st.balking.size() &&
                        st.balking[r].strategy != lang::BalkingStrategy::NONE)
                        u.set(Feature::Balking);
                }
                if (!st.server_types.empty()) u.set(Feature::HeteroServers);
                for (std::size_t n : st.server_parallelism)
                    if (n > 1) { u.set(Feature::ServerParallelism); break; }
                break;
            }
            case NodeType::Source: {
                u.set({Feature::Source, Feature::RandomSource, Feature::ServiceTunnel,
                       Feature::Dispatcher});
                // The Source branch emits the arrival distributions only: MATLAB
                // registers no discipline and no routing for it.
                if (ist == 0 || ist - 1 >= sn.service.size()) break;
                for (std::size_t r = 0; r < R && r < sn.service[ist - 1].size(); ++r)
                    u.set(feature_of_process(sn.service[ist - 1][r].type));
                // A batch arrival releases several jobs per epoch, which no
                // single-arrival event handler produces.
                for (const Distrib<T>& b : sn.stations[ist - 1].arrival_batch)
                    if (!b.disabled) { u.set(Feature::BatchArrival); break; }
                break;
            }
            case NodeType::Sink:
                u.set(Feature::Sink);
                break;
            case NodeType::Router: {
                // A ROUTER EMITS ITS ROUTING STRATEGY AND NOTHING ELSE. All
                // three reference codebases agree and this port was the only
                // one that did not: getUsedLangFeatures.m:77-78 (case 'Router'),
                // Network.java:3184-3190 and network.py:1983-1987 each record
                // RoutingStrategy.toFeature(...) alone. This branch additionally
                // emitted Router, Buffer, Dispatcher and ServiceTunnel, and
                // since neither mva_feature_set nor nc_feature_set declares
                // Router -- correctly, SolverMVA.m:231-257 does not either --
                // every model with a dispatcher was refused by the two solvers
                // MATLAB and python solve it with. The whole cluster family
                // (Network.cluster/clusterClosed/clusterMixed) is exactly that
                // shape, and the Router is eliminated by the stochastic
                // complement long before any solver sees it.
                for (std::size_t r = 0; r < R && r < nd.routing.size(); ++r)
                    u.set(feature_of_routing(nd.routing[r]));
                break;
            }
            case NodeType::ClassSwitch:
                u.set({Feature::ClassSwitch, Feature::StatelessClassSwitcher, Feature::Buffer,
                       Feature::Dispatcher});
                break;
            case NodeType::Cache: {
                u.set({Feature::Cache, Feature::CacheClassSwitcher, Feature::Buffer,
                       Feature::Dispatcher});
                auto it = sn.nodeparam.find(ind);
                if (it == sn.nodeparam.end()) break;
                u.set(feature_of_replacement(it->second.replacestrat));
                if (it->second.retrieval_capacity > 0 || !it->second.retrieval_classes.empty())
                    u.set(Feature::CacheRetrieval);
                if (!it->second.costcap.empty()) u.set(Feature::CacheItemSize);
                break;
            }
            case NodeType::Fork: {
                u.set({Feature::Fork, Feature::Forker, Feature::Buffer, Feature::ServiceTunnel});
                // Variable forking levels. A name declared but never marked is
                // a name no solver can ever refuse, so all three are marked
                // here. A plain fork has no override block at all, so this
                // costs one map lookup in the classic case.
                const qn::ForkParam<T>* fp = sn.fork_param_of(i);
                if (fp != 0) {
                    const double tpl = sn.nodes[i - 1].tasks_per_link;
                    bool vector_valued = false, random_valued = false, probabilistic = false;
                    for (std::size_t k = 0; k < fp->fan_out_link.rows(); ++k)
                        for (std::size_t r = 0; r < fp->fan_out_link.cols(); ++r) {
                            const double p = num_traits<T>::to_double(fp->fan_out_prob(k, r));
                            if (p == 0.0) continue;  // link this class does not take
                            if (p != 1.0) probabilistic = true;
                            if (num_traits<T>::to_double(fp->fan_out_link(k, r)) != tpl)
                                vector_valued = true;
                            if (!fp->fan_out_dist[k][r].disabled) random_valued = true;
                        }
                    if (vector_valued) u.set(Feature::ForkFanoutVector);
                    if (random_valued) u.set(Feature::ForkFanoutRandom);
                    if (probabilistic) u.set(Feature::ForkBranchProbability);
                }
                break;
            }
            case NodeType::Join:
                u.set({Feature::Join, Feature::Joiner, Feature::Dispatcher,
                       Feature::ServiceTunnel});
                break;
            case NodeType::Logger:
                u.set({Feature::Logger, Feature::LogTunnel, Feature::Buffer, Feature::Dispatcher});
                break;
            case NodeType::Place:
                u.set({Feature::Place, Feature::Storage, Feature::Linkage,
                       Feature::ServiceTunnel});
                break;
            case NodeType::Transition: {
                u.set({Feature::Transition, Feature::Enabling, Feature::Timing, Feature::Firing});
                auto it = sn.transparam.find(ind);
                if (it == sn.transparam.end()) break;
                // Inhibiting only for a FINITE threshold: Inf is "never blocks",
                // and a plain SPN must not be gated out of solvers lacking it.
                for (const auto& row : it->second.inhibiting) {
                    bool finite = false;
                    for (std::size_t p = 0; p < row.rows() && !finite; ++p)
                        for (std::size_t r = 0; r < row.cols(); ++r)
                            if (std::isfinite(num_traits<T>::to_double(row(p, r)))) {
                                finite = true;
                                break;
                            }
                    if (finite) { u.set(Feature::Inhibiting); break; }
                }
                break;
            }
            case NodeType::Region:
                u.set(Feature::Region);
                break;
        }
    }

    // A ClassSwitch is never materialised as a node by this port (see the
    // NetworkStruct header), so the switch matrix is the only trace it leaves.
    if (!sn.csmatrix.empty())
        u.set({Feature::ClassSwitch, Feature::StatelessClassSwitcher, Feature::Buffer,
               Feature::Dispatcher});

    if (!sn.retrialparam.empty()) u.set(Feature::Retrial);
    // Breakdown became a STRUCT FIELD on 2026-08-15, so it is derivable here
    // where balking and reneging still are not: a solver that does not declare
    // it is now refused on a model that carries it, instead of running as
    // though the server never failed.
    if (!sn.breakdownparam.empty()) u.set(Feature::Breakdown);
    if (!sn.setupparam.empty()) u.set(Feature::SetupDelayOff);
    // A PARTIAL join fires on a quorum rather than on every sibling, and a FIFO
    // depository releases a served token only after the earlier ones: both
    // change which transitions exist, so a solver that has neither must refuse.
    for (const auto& kv : sn.joindecl)
        if (kv.second.strategy == lang::JoinStrategy::PARTIAL) u.set(Feature::JoinPartial);
    for (const Station<T>& st : sn.stations)
        for (lang::DepartureDiscipline d : st.departure_discipline)
            if (d != lang::DepartureDiscipline::NORMAL) u.set(Feature::DepartureDiscipline);
    for (const auto& row : sn.droprule)
        for (DropStrategy d : row)
            if (d == DropStrategy::RETRIAL || d == DropStrategy::RETRIAL_WITH_LIMIT)
                u.set(Feature::Retrial);

    // A finite-server station serving several jobs at once. A Delay carries an
    // infinite nservers and is NOT one, which is what `has_multi_server` already
    // tests: the single-server recursions answered a c-server station as one
    // server of the same rate, and only a structural predicate could say so.
    if (sn.has_multi_server()) u.set(Feature::MultiServer);
    // A station or per-class buffer that can BIND: the one predicate
    // `check_binding_capacity` refuses on (station caps against the class
    // populations, open classes always bind, Cache models exempt), asked here so
    // the refusal has a registry name and a solver method that does not declare
    // it is gated on exactly the models the structural gate refuses.
    if (has_binding_capacity(sn)) u.set(Feature::FiniteCapacity);

    return u;
}

// ---------------------------------------------------------------------------
// The gate
// ---------------------------------------------------------------------------

/**
 * runAnalyzerChecks: refuse a model the solver does not declare, by name.
 *
 * Throws UnsupportedError, which is what every other by-name refusal in this
 * port throws and what the CLI reports without a stack.
 *
 * THE TWO MESSAGE FORMS, NetworkSolver.m:184-190. The reference emits a
 * different message when `resolveMethod` changed the method under the user's
 * feet than when the user named it:
 *
 *   method == options.method  "features not supported by the solver"
 *   otherwise                 "features not supported by the solver's '%s' method"
 *
 * That distinction is not cosmetic. The only resolution today is MVA's
 * `default` -> `rqna` upgrade, and rqna's set WITHDRAWS ClosedClass and
 * SelfLoopingClass, so a user who asked for nothing but `default` can be
 * refused for having a closed class. Without the method in the message the
 * refusal cannot be understood, because nothing the user typed mentions rqna.
 *
 * Both method arguments default to empty, which selects the plain form and
 * leaves every call site that does not resolve a method unchanged.
 */
template <class T>
void feature_gate(const std::string& solver, const FeatureSet& declared,
                  const NetworkStruct<T>& sn, const std::string& requested_method = "",
                  const std::string& resolved_method = "") {
    const bool upgraded = !resolved_method.empty() && resolved_method != requested_method;
    const std::string label = upgraded ? solver + "'s '" + resolved_method + "' method" : solver;
    const SupportResult r = feature_set_supports(label, declared, used_lang_features(sn));
    if (!r.ok) throw UnsupportedError(r.reason);
}

/** `%g` for a capacity, so the message reads 2 and not 2.000000. */
inline std::string capacity_str(double v) {
    std::ostringstream os;
    os << v;
    return os.str();
}

/**
 * NetworkSolver.checkBindingCapacity (NetworkSolver.m:1172-1228): the shared
 * structural gate for finite station capacity (setCapacity) and finite
 * per-class buffers (classCap), used by the product-form solvers.
 *
 * There is no LINE_QN_FEATURE_LIST enumerator for plain capacity -- it is a
 * NUMBER on a station, not a construct -- so `feature_gate` above cannot see
 * it, and a product-form solver has no representation of a finite buffer.
 * Without this check the solve returns the UNCONSTRAINED answer (QLen=4 where
 * the M/M/1/2 value is 0.8525), which is a wrong number rather than a refusal.
 *
 * The test reads the STATION-level cap/classcap the user set, never the derived
 * `sn.cap`/`sn.classcap`: refresh_capacity gives every closed model a finite
 * derived classcap (the chain population), so an sn-level test would reject
 * every closed model.
 *
 * Only a capacity that can actually BIND is refused. A closed model whose
 * station capacity is at least the total population can never block a job, so
 * the declaration is a no-op and the product-form answer stays exact -- a
 * common idiom is setCapacity(N) on a station of an N-job closed model. The
 * population is infinite for an open class, so any finite capacity an open
 * class reaches binds.
 *
 * Cache models are exempt: a Cache sets classcap=1 on the retrieval queues it
 * builds, and MVA/NC solve those through their cache analyzers rather than as a
 * buffer constraint.
 */
/**
 * Which solvers to point at when a finite capacity is refused.
 *
 * The two lists differ, and naming the wrong one sends the user to a solver that
 * also refuses. An OPEN refused arrival is LOST, which SolverJMT reproduces (its
 * queue section carries the drop rule directly). A CLOSED one BLOCKS: LINE
 * disables the upstream departure and holds the job where it is, and no JMT drop
 * strategy expresses that -- 'waiting queue' does not enforce the size at all and
 * 'BAS blocking' completes the service before blocking, a different queueing
 * model. The JMT writer refuses the closed case by name
 * (`assert_station_cap_exportable` in `jmt_writer.h`, BUG-81), so it must not be
 * advertised here for it.
 */
inline std::string capacity_fallback_advice(bool is_open_class) {
    return is_open_class ? "Use SolverCTMC, SolverJMT or SolverLDES"
                         : "Use SolverCTMC, SolverSSA or SolverLDES";
}

/**
 * The refusal `check_binding_capacity` raises, as a string, or empty when no
 * buffer binds. Split out so that a caller which must DECIDE on the same rule --
 * `fluid_resolve_method`, which sends a blocked model to the one method that
 * carries the constraint -- asks the gate's own question rather than a second
 * copy of it that could drift from it.
 */
template <class T>
std::string binding_capacity_reason(const std::string& solver, const NetworkStruct<T>& sn) {
    for (const NodeDef& nd : sn.nodes)
        if (nd.nodetype == NodeType::Cache) return std::string();

    double totaljobs = 0.0;  // infinite as soon as one class is open
    bool any_open = false;
    for (const JobClass& c : sn.classes) {
        totaljobs += c.population;
        if (std::isinf(c.population)) any_open = true;
    }

    for (const Station<T>& st : sn.stations) {
        if (st.nodetype == NodeType::Source || st.nodetype == NodeType::Sink) continue;
        if (std::isfinite(st.cap) && st.cap >= 0.0 && st.cap < totaljobs)
            return solver + ": finite station capacity (setCapacity=" + capacity_str(st.cap) +
                   ") at station '" + st.name +
                   "' is not supported by this solver, which has no "
                   "representation of a finite buffer and would return the "
                   "unconstrained answer. " +
                   capacity_fallback_advice(any_open);
        const std::size_t n = std::min(st.classcap.size(), sn.classes.size());
        for (std::size_t r = 0; r < n; ++r)
            if (std::isfinite(st.classcap[r]) && st.classcap[r] > 0.0 &&
                st.classcap[r] < sn.classes[r].population)
                return solver + ": finite per-class capacity (classCap=" +
                       capacity_str(st.classcap[r]) + " for class " + std::to_string(r + 1) +
                       ") at station '" + st.name +
                       "' is not supported by this solver, which has no representation of a "
                       "finite buffer and would return the unconstrained answer. " +
                       capacity_fallback_advice(std::isinf(sn.classes[r].population));
    }
    return std::string();
}

/** Does any station buffer BIND? The gate's own question, asked without raising. */
template <class T>
bool has_binding_capacity(const NetworkStruct<T>& sn) {
    return !binding_capacity_reason(std::string("solver"), sn).empty();
}

template <class T>
void check_binding_capacity(const std::string& solver, const NetworkStruct<T>& sn) {
    const std::string reason = binding_capacity_reason(solver, sn);
    if (!reason.empty()) throw UnsupportedError(reason);
}

}  // namespace qn
}  // namespace line

#endif  // LINE_LANG_QN_FEATURE_SET_H
