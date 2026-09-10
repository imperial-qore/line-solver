/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_LANG_LANG_TYPES_H
#define LINE_LANG_LANG_TYPES_H

/**
 * Enumerations and the minimal distribution descriptor shared by the model
 * layer of the C++ port.
 *
 * The numeric values are the MATLAB ones (matlab/src/lang/constant), not a
 * fresh numbering, because they cross the JSON boundary and appear in the
 * dumps used as parity oracles. _kb/11 records that MATLAB and Python disagree
 * on ProcessType numbering (MATLAB starts at EXP=0, Python at EXP=1); this port
 * follows MATLAB, the reference implementation, so a numeric comparison against
 * a MATLAB dump is meaningful and one against a Python dump is not -- compare
 * by name there.
 *
 * SCOPE: the model layer added here exists to run SolverLN over a layered
 * queueing network whose layers are solved by SolverMVA. It carries the
 * scheduling disciplines, node kinds and precedence types that path reaches
 * and REFUSES the rest by name rather than silently mapping them onto a
 * neighbour, because a discipline that is silently treated as FCFS returns a
 * plausible number that is wrong.
 */

#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace lang {

/** Solver output metrics, with the numeric values of MATLAB `MetricType`. */
enum class MetricType {
    ResidT = 0,
    RespT = 1,
    DropRate = 2,
    QLen = 3,
    QueueT = 4,
    FCRWeight = 5,
    FCRMemOcc = 6,
    FJQLen = 7,
    FJRespT = 8,
    RespTSink = 9,
    SysQLen = 10,
    SysRespT = 11,
    SysTput = 12,
    Tput = 13,
    ArvR = 14,
    TputSink = 15,
    Util = 16,
    TranQLen = 17,
    TranUtil = 18,
    TranTput = 19,
    TranRespT = 20,
    Tard = 21,
    SysTard = 22
};

/** Port of `MetricType.toText`. */
inline const char* metric_to_text(MetricType metric) {
    switch (metric) {
        case MetricType::ResidT: return "Residence Time";
        case MetricType::RespT: return "Response Time";
        case MetricType::DropRate: return "Drop Rate";
        case MetricType::QLen: return "Number of Customers";
        case MetricType::QueueT: return "Queue Time";
        case MetricType::FCRWeight: return "FCR Total Weight";
        case MetricType::FCRMemOcc: return "FCR Memory Occupation";
        case MetricType::FJQLen: return "Fork Join Number of Customers";
        case MetricType::FJRespT: return "Fork Join Response Time";
        case MetricType::RespTSink: return "Response Time per Sink";
        case MetricType::SysQLen: return "System Number of Customers";
        case MetricType::SysRespT: return "System Response Time";
        case MetricType::SysTput: return "System Throughput";
        case MetricType::Tput: return "Throughput";
        case MetricType::ArvR: return "Arrival Rate";
        case MetricType::TputSink: return "Throughput per Sink";
        case MetricType::Util: return "Utilization";
        case MetricType::TranQLen: return "Tran Number of Customers";
        case MetricType::TranUtil: return "Tran Utilization";
        case MetricType::TranTput: return "Tran Throughput";
        case MetricType::TranRespT: return "Tran Response Time";
        case MetricType::Tard: return "Tardiness";
        case MetricType::SysTard: return "System Tardiness";
        default: return "Unknown Metric";
    }
}

/**
 * The events a state can undergo, with the values of MATLAB EventType.
 *
 * An event is ACTIVE at the node that schedules it and PASSIVE at the node
 * that receives it: a DEP at one station is the ARV at the next, and only the
 * active half carries a rate. The passive half is marked with a rate of -1,
 * which the generator assembly replaces with the active rate -- a convention
 * that only reads as a sentinel because a rate can never be negative.
 */
enum class EventType {
    INIT = -1,     ///< the model is initialized, t = 0
    LOCAL = 0,     ///< dummy event, no state change outside the node
    ARV = 1,       ///< a job arrives
    DEP = 2,       ///< a job departs
    PHASE = 3,     ///< service advances a phase WITHOUT departing
    READ = 4,      ///< a cache item is read
    STAGE = 5,     ///< a random environment changes stage
    ENABLE = 6,    ///< an SPN mode becomes enabled
    FIRE = 7,      ///< an SPN mode fires
    PRE = 8,       ///< consume from a place or queue buffer, no server effect
    POST = 9,      ///< produce to a place or queue buffer
    RENEGE = 10,   ///< a waiting job abandons the queue (impatience)
    RETRY = 11,    ///< an orbiting job retries entry at a retrial station
    SWITCH = 12,   ///< a polling server advances its switchover timer
    FAILURE = 13,  ///< the server breaks down, going from up to down
    REPAIR = 14,   ///< the server is repaired, going from down to up, resuming
                   ///< the held job, which is why it emits no START
    START = 15,    ///< a job begins or resumes holding a server
    PREEMPT = 16   ///< a job holding a server is pushed back into the buffer
};
// START and PREEMPT are instantaneous tags on the arc of the ARV or DEP that
// causes them, never the active half of an sn.sync entry: no clock, no state,
// no change to any numerical result. PREEMPT is spelled in full because PRE
// already names the Petri-net pre-arc.

inline const char* event_to_text(EventType e) {
    switch (e) {
        case EventType::ARV: return "ARV";
        case EventType::DEP: return "DEP";
        case EventType::PHASE: return "PHASE";
        case EventType::READ: return "READ";
        case EventType::LOCAL: return "LOCAL";
        case EventType::STAGE: return "STAGE";
        case EventType::ENABLE: return "ENABLE";
        case EventType::FIRE: return "FIRE";
        case EventType::PRE: return "PRE";
        case EventType::POST: return "POST";
        case EventType::RENEGE: return "RENEGE";
        case EventType::RETRY: return "RETRY";
        case EventType::SWITCH: return "SWITCH";
        case EventType::FAILURE: return "FAILURE";
        case EventType::REPAIR: return "REPAIR";
        case EventType::START: return "START";
        case EventType::PREEMPT: return "PREEMPT";
        default: return "INIT";
    }
}

/**
 * G-network signal classes, with the values of MATLAB SignalType.
 *
 * A signal is not a job: it never joins a station, it acts on the jobs already
 * there and is annihilated. REPLY is the odd one out -- it completes a
 * synchronous call and then joins as an ordinary job.
 */
enum class SignalType {
    REPLY = 0,       ///< completes a synchronous call, releasing a held server
    NEGATIVE = 1,    ///< removes a batch of jobs (Gelenbe's negative customer)
    CATASTROPHE = 2  ///< removes EVERY job at the station
};

/** Which job a negative signal removes, with the values of MATLAB RemovalPolicy. */
enum class RemovalPolicy {
    RANDOM = 0,  ///< uniform over waiting AND in-service jobs
    FCFS = 1,    ///< the oldest waiting job; servers only once nobody waits
    LCFS = 2     ///< the newest waiting job; servers only once nobody waits
};

/** Scheduling disciplines, with the values of MATLAB SchedStrategy. */
enum class SchedStrategy {
    INF = 0,
    FCFS = 1,
    LCFS = 2,
    SIRO = 3,
    SJF = 4,
    LJF = 5,
    PS = 6,
    DPS = 7,
    GPS = 8,
    SEPT = 9,
    LEPT = 10,
    HOL = 11,
    FORK = 12,
    EXT = 13,
    REF = 14,
    LCFSPR = 15,
    POLLING = 16,
    // The preemptive and priority families. PR resumes an interrupted job in
    // the phase it held; PI restarts it from the entry phase, which is why the
    // two cannot share an encoding: PR must carry the phase of EVERY preempted
    // job, PI need not. FCFSPRIO is MATLAB's alias for HOL, so it is not a
    // distinct enumerator here.
    PSPRIO = 17,
    DPSPRIO = 18,
    GPSPRIO = 19,
    LCFSPI = 20,
    LCFSPRIO = 21,
    LCFSPRPRIO = 22,
    LCFSPIPRIO = 23,
    FCFSPR = 24,
    FCFSPI = 25,
    FCFSPRPRIO = 26,
    FCFSPIPRIO = 27,
    SRPT = 28,
    SRPTPRIO = 29,
    EDD = 30,
    EDF = 31,
    LPS = 32,
    PSJF = 33,
    FB = 34,
    LRPT = 35,
    SETF = 36,
    FSP = 37,
    PAS = 38,
    OI = 39,
    NONE = -1
};

inline const char* sched_to_text(SchedStrategy s) {
    switch (s) {
        case SchedStrategy::INF: return "inf";
        case SchedStrategy::FCFS: return "fcfs";
        case SchedStrategy::LCFS: return "lcfs";
        case SchedStrategy::SIRO: return "siro";
        case SchedStrategy::SJF: return "sjf";
        case SchedStrategy::LJF: return "ljf";
        case SchedStrategy::PS: return "ps";
        case SchedStrategy::DPS: return "dps";
        case SchedStrategy::GPS: return "gps";
        case SchedStrategy::SEPT: return "sept";
        case SchedStrategy::LEPT: return "lept";
        case SchedStrategy::HOL: return "hol";
        case SchedStrategy::FORK: return "fork";
        case SchedStrategy::EXT: return "ext";
        case SchedStrategy::REF: return "ref";
        case SchedStrategy::LCFSPR: return "lcfspr";
        case SchedStrategy::POLLING: return "polling";
        case SchedStrategy::SRPT: return "srpt";
        case SchedStrategy::LPS: return "lps";
        case SchedStrategy::PSJF: return "psjf";
        case SchedStrategy::FB: return "fb";
        case SchedStrategy::LRPT: return "lrpt";
        case SchedStrategy::SETF: return "setf";
        case SchedStrategy::FSP: return "fsp";
        case SchedStrategy::PSPRIO: return "psprio";
        case SchedStrategy::DPSPRIO: return "dpsprio";
        case SchedStrategy::GPSPRIO: return "gpsprio";
        case SchedStrategy::LCFSPI: return "lcfspi";
        case SchedStrategy::LCFSPRIO: return "lcfsprio";
        case SchedStrategy::LCFSPRPRIO: return "lcfsprprio";
        case SchedStrategy::LCFSPIPRIO: return "lcfspiprio";
        case SchedStrategy::FCFSPR: return "fcfspr";
        case SchedStrategy::FCFSPI: return "fcfspi";
        case SchedStrategy::FCFSPRPRIO: return "fcfsprprio";
        case SchedStrategy::FCFSPIPRIO: return "fcfspiprio";
        case SchedStrategy::SRPTPRIO: return "srptprio";
        case SchedStrategy::EDD: return "edd";
        case SchedStrategy::EDF: return "edf";
        case SchedStrategy::PAS: return "pas";
        case SchedStrategy::OI: return "oi";
        default: return "none";
    }
}

/** Parse the `scheduling` attribute of an .lqnx processor or task. */
inline SchedStrategy sched_from_lqnx(const std::string& s) {
    if (s == "inf" || s == "INF") return SchedStrategy::INF;
    if (s == "fcfs" || s == "FCFS") return SchedStrategy::FCFS;
    if (s == "ps" || s == "PS") return SchedStrategy::PS;
    if (s == "ref" || s == "REF") return SchedStrategy::REF;
    if (s == "hol" || s == "HOL") return SchedStrategy::HOL;
    // lqns spells preemptive priority resume (SCHEDULE_PPR) "pri"; "pp" is the
    // stale lqn-core.xsd spelling, absent from the lqns 6.2.31 sources. It is
    // preemptive, so it is FCFSPRPRIO and not the non-preemptive HOL
    if (s == "pri" || s == "PRI" || s == "pp") return SchedStrategy::FCFSPRPRIO;
    if (s == "rand" || s == "siro") return SchedStrategy::SIRO;
    if (s == "sjf") return SchedStrategy::SJF;
    if (s == "ljf") return SchedStrategy::LJF;
    if (s == "lcfs") return SchedStrategy::LCFS;
    if (s == "burst" || s == "poll") return SchedStrategy::FCFS;
    throw UnsupportedError("lqn reader: unsupported scheduling discipline '" + s + "'");
}

/**
 * The `scheduling` attribute an .lqnx processor or task carries for a strategy.
 *
 * The inverse of sched_from_lqnx over the disciplines the schema spells, and a
 * refusal by name for every other one. It refuses rather than falling back on
 * fcfs because the file is handed to an external solver: a task written as fcfs
 * when the model says lcfspr is answered, not rejected, and the discipline
 * would be lost inside a number that looks ordinary.
 */
inline std::string sched_to_lqnx(SchedStrategy s) {
    switch (s) {
        case SchedStrategy::INF: return "inf";
        case SchedStrategy::FCFS: return "fcfs";
        case SchedStrategy::PS: return "ps";
        case SchedStrategy::REF: return "ref";
        case SchedStrategy::HOL: return "hol";
        case SchedStrategy::FCFSPRPRIO: return "pri";
        case SchedStrategy::SIRO: return "rand";
        case SchedStrategy::SJF: return "sjf";
        case SchedStrategy::LJF: return "ljf";
        case SchedStrategy::LCFS: return "lcfs";
        default:
            throw UnsupportedError(
                std::string("the LQN XML schema has no spelling for scheduling discipline '") +
                sched_to_text(s) + "'; it carries inf, fcfs, ps, ref, hol, pri, rand, sjf, ljf and lcfs");
    }
}

/** Node kinds, with the values of MATLAB NodeType. */
enum class NodeType {
    Queue = 0,
    Source = 1,
    Delay = 2,
    ClassSwitch = 3,
    Logger = 4,
    Cache = 5,
    Router = 6,
    Fork = 7,
    Place = 8,
    Transition = 9,
    Region = 10,
    Join = 11,
    Sink = 12
};

/** Name of a node kind, for diagnostics. The JSON spelling lives in the writer. */
inline const char* node_type_to_text(NodeType t) {
    switch (t) {
        case NodeType::Queue: return "Queue";
        case NodeType::Source: return "Source";
        case NodeType::Delay: return "Delay";
        case NodeType::ClassSwitch: return "ClassSwitch";
        case NodeType::Logger: return "Logger";
        case NodeType::Cache: return "Cache";
        case NodeType::Router: return "Router";
        case NodeType::Fork: return "Fork";
        case NodeType::Place: return "Place";
        case NodeType::Transition: return "Transition";
        case NodeType::Region: return "Region";
        case NodeType::Join: return "Join";
        case NodeType::Sink: return "Sink";
    }
    return "Unknown";
}

/** SPN transition timing, with the values of MATLAB TimingStrategy. */
enum class TimingStrategy {
    TIMED = 0,      ///< fires after its firing distribution elapses
    IMMEDIATE = 1   ///< fires with zero delay, resolved by weight and priority
};

/** Job class kinds, with the values of MATLAB JobClassType. */
enum class JobClassType { OPEN = 0, CLOSED = 1 };

/** Polling service disciplines, with the values of MATLAB PollingType. */
enum class PollingType {
    GATED = 0,        ///< serve exactly the jobs present at the polling instant
    EXHAUSTIVE = 1,   ///< serve until the queue empties
    KLIMITED = 2,     ///< serve at most K per visit (K in pollingPar)
    DECREMENTING = 3  ///< serve until the queue is one shorter than at arrival
};

/** Cache replacement policies, with the values of MATLAB ReplacementStrategy. */
enum class ReplacementStrategy {
    RR = 0,     ///< random replacement
    FIFO = 1,   ///< first in, first out
    SFIFO = 2,  ///< strict FIFO
    LRU = 3,    ///< least recently used
    HLRU = 4,   ///< h-LRU / LRU(m): h lists, promote i -> i+1 on a hit
    CLIMB = 5,  ///< move up one position on a hit (transposition rule)
    QLRU = 6    ///< q-LRU: LRU with probabilistic admission on a miss
};

/** Routing strategies, with the values of MATLAB RoutingStrategy. */
enum class RoutingStrategy {
    RAND = 0,
    PROB = 1,
    RROBIN = 2,
    WRROBIN = 3,
    JSQ = 4,
    FIRING = 5,
    SQ = 6,
    /** Krzesinski (1987) product-form state-dependent routing. */
    SDR = 7,
    DISABLED = -1
};

inline const char* routing_to_text(RoutingStrategy r) {
    switch (r) {
        case RoutingStrategy::RAND: return "rand";
        case RoutingStrategy::PROB: return "prob";
        case RoutingStrategy::RROBIN: return "rrobin";
        case RoutingStrategy::WRROBIN: return "wrrobin";
        case RoutingStrategy::JSQ: return "jsq";
        case RoutingStrategy::FIRING: return "firing";
        case RoutingStrategy::SQ: return "sq";
        case RoutingStrategy::SDR: return "sdr";
        default: return "disabled";
    }
}

/**
 * Blocking and loss rules, with the values of MATLAB DropStrategy.
 *
 * WAITQ is -1 and is also the marker `refreshCapacity` writes where the rule is
 * never consulted (an unbounded station, or a closed class), so it means two
 * different things depending on the station's capacity; see the comment in
 * refresh_capacity().
 */
enum class DropStrategy {
    WAITQ = -1,
    DROP = 1,
    BAS = 2,
    BBS = 3,
    RSRD = 4,
    RETRIAL = 5,
    RETRIAL_WITH_LIMIT = 6
};

/**
 * Impatience kinds, with the values of MATLAB ImpatienceType.
 *
 * RENEGING is a timer a job started at a queue; BALKING is a decision taken
 * BEFORE joining, on the state of the queue, and is parameterized by
 * `Station::balking` rather than by a distribution; RETRIAL sends the job to an
 * orbit and is parameterized by `RetrialParam`.
 */
enum class ImpatienceType { NONE = 0, RENEGING = 1, BALKING = 2, RETRIAL = 3 };

/** Balking rules, with the values of MATLAB BalkingStrategy. */
enum class BalkingStrategy { NONE = 0, QUEUE_LENGTH = 1, EXPECTED_WAIT = 2, COMBINED = 3 };

/**
 * How a heterogeneous station picks among its server types, MATLAB
 * HeteroSchedPolicy. ORDER is the default: the declared order of the types.
 */
enum class HeteroSchedPolicy { ORDER = 0, ALIS = 1, ALFS = 2, FAIRNESS = 3, FSF = 4, RAIS = 5 };

/**
 * When a Place releases a served token, MATLAB DepartureDiscipline. NORMAL is
 * the standard queueing-Petri-net rule (available on completion); FIFO holds
 * it until every earlier arrival to the depository has been released.
 */
enum class DepartureDiscipline { NORMAL = 0, FIFO = 1 };

/** Join rules, with the values of MATLAB JoinStrategy. */
enum class JoinStrategy { STD = 1, PARTIAL = 2 };

/** LQN element kinds, with the values of MATLAB LayeredNetworkElement. */
enum class LqnElement { HOST = 0, TASK = 1, ENTRY = 2, ACTIVITY = 3, CALL = 4 };

/** Call kinds, with the values of MATLAB CallType. */
enum class CallType { NONE = 0, SYNC = 1, ASYNC = 2, FWD = 3 };

/** Activity precedence kinds, with the values of MATLAB ActivityPrecedenceType. */
enum class PrecedenceType {
    NONE = 0,
    PRE_SEQ = 1,
    PRE_AND = 2,
    PRE_OR = 3,
    POST_SEQ = 11,
    POST_AND = 12,
    POST_OR = 13,
    POST_LOOP = 14,
    POST_CACHE = 15
};

/** Distribution kinds, with the values of MATLAB ProcessType. */
enum class ProcessType {
    EXP = 0,
    ERLANG = 1,
    HYPEREXP = 2,
    PH = 3,
    APH = 4,
    MAP = 5,
    UNIFORM = 6,
    DET = 7,
    COXIAN = 8,
    GAMMA = 9,
    PARETO = 10,
    MMPP2 = 11,
    REPLAYER = 12,
    IMMEDIATE = 13,
    DISABLED = 14,
    COX2 = 15,
    WEIBULL = 16,
    LOGNORMAL = 17,
    DUNIFORM = 18,
    BERNOULLI = 19,
    /**
     * A `Prior`: a weighted set of ALTERNATIVE distributions, or a density over
     * a scalar parameter plus a factory from it. It is not a mixture -- each
     * alternative is a separate model realization -- and only SolverUQ consumes
     * it; every other solver refuses it through Feature::Prior.
     */
    PRIOR = 20,
    BINOMIAL = 21,
    POISSON = 22,
    GEOMETRIC = 23,
    BMAP = 24,
    ME = 25,
    RAP = 26,
    DISCRETESAMPLER = 27,
    ZIPF = 28,
    DMAP = 29,
    MMAP = 31,
    EMPIRICALCDF = 32,
    /**
     * The time-INHOMOGENEOUS families of Ko and Pender (ORL 45, 2017): an
     * NHPP is a rate schedule lambda(t), a MAPt a (D0(t), D1(t)) schedule and a
     * PHt an (alpha(t), S(t)) one, all piecewise constant on one breakpoint
     * vector and optionally cyclic. The numeric values are MATLAB's
     * (`ProcessType.m:41-43`).
     *
     * THEY CARRY A NOMINAL PAIR TOO. `Distrib::D0`/`D1` hold the width-weighted
     * time average of the schedule, which is what `sn_schedule_nominal` returns
     * as its first two outputs and what every consumer that has no notion of
     * time -- the phase count, the rate, the fluid layout -- reads. The schedule
     * itself lives in `sched_bp`/`sched_D0`/`sched_D1` beside it, and only a
     * solver that integrates in time looks at it.
     */
    NHPP = 33,
    MAPT = 34,
    PHT = 35,
    /**
     * A Gaussian, and the ONE family whose value is not MATLAB's, because
     * MATLAB has none to copy: `ProcessType.m` stops at 35 and `Normal.m` is a
     * `ContinuousDistribution` with no id at all, exactly as `Normal.java` and
     * the python `Normal` have none. That is not an oversight in the reference
     * -- a Gaussian has mass below zero, so it can never be a service or
     * interarrival process and can never appear in `sn.proc`. It reaches this
     * port only as the PARAMETER density of a continuous `Prior`, where it is
     * read through `dist_cdf` and `dist_quantile` and never through
     * `dist_to_map`.
     *
     * The value is deliberately far outside 0..35 so that it can never collide
     * with an id MATLAB assigns later; `sn.procid` must never carry it, and
     * `dist_to_map` refuses it by name rather than handing back the Erlang fit
     * its default arm would otherwise produce for a law with negative support.
     */
    NORMAL = 100,
    NONE = -1
};

/** The MATLAB ProcessType name, as `sn.procid` prints it. */
inline const char* process_to_text(ProcessType p) {
    switch (p) {
        case ProcessType::EXP: return "Exp";
        case ProcessType::ERLANG: return "Erlang";
        case ProcessType::HYPEREXP: return "HyperExp";
        case ProcessType::PH: return "PH";
        case ProcessType::APH: return "APH";
        case ProcessType::MAP: return "MAP";
        case ProcessType::UNIFORM: return "Uniform";
        case ProcessType::DET: return "Det";
        case ProcessType::COXIAN: return "Coxian";
        case ProcessType::GAMMA: return "Gamma";
        case ProcessType::PARETO: return "Pareto";
        case ProcessType::MMPP2: return "MMPP2";
        case ProcessType::REPLAYER: return "Replayer";
        case ProcessType::IMMEDIATE: return "Immediate";
        case ProcessType::DISABLED: return "Disabled";
        case ProcessType::COX2: return "Cox2";
        case ProcessType::WEIBULL: return "Weibull";
        case ProcessType::LOGNORMAL: return "Lognormal";
        case ProcessType::DUNIFORM: return "DiscreteUniform";
        case ProcessType::BERNOULLI: return "Bernoulli";
        case ProcessType::BINOMIAL: return "Binomial";
        case ProcessType::POISSON: return "Poisson";
        case ProcessType::GEOMETRIC: return "Geometric";
        case ProcessType::BMAP: return "BMAP";
        case ProcessType::ME: return "ME";
        case ProcessType::RAP: return "RAP";
        case ProcessType::DISCRETESAMPLER: return "DiscreteSampler";
        case ProcessType::ZIPF: return "Zipf";
        case ProcessType::DMAP: return "DMAP";
        case ProcessType::MMAP: return "MMAP";
        case ProcessType::EMPIRICALCDF: return "EmpiricalCDF";
        case ProcessType::PRIOR: return "Prior";
        case ProcessType::NHPP: return "NHPP";
        case ProcessType::MAPT: return "MAPt";
        case ProcessType::PHT: return "PHt";
        case ProcessType::NORMAL: return "Normal";
        default: return "none";
    }
}

/**
 * `ProcessType.isMarkovian`: true when `sn.proc` carries an exact matrix
 * representation of the law, rather than the Erlang fit `convertToMAP` leaves
 * there for the parameter-only families. ME and RAP count -- their
 * representation is the matrix-exponential analogue, not a generator -- and the
 * discrete families do not, so a solver reading `sn.proc` as the law must gate
 * on this exactly as MATLAB `ProcessType.m` does.
 */
inline bool process_is_markovian(ProcessType p) {
    switch (p) {
        case ProcessType::EXP:
        case ProcessType::ERLANG:
        case ProcessType::HYPEREXP:
        case ProcessType::PH:
        case ProcessType::APH:
        case ProcessType::MAP:
        case ProcessType::COXIAN:
        case ProcessType::COX2:
        case ProcessType::MMPP2:
        case ProcessType::ME:
        case ProcessType::RAP:
        case ProcessType::DMAP:
        case ProcessType::BMAP:
        case ProcessType::MMAP: return true;
        default: return false;
    }
}

/**
 * A class-dependent scaling map, `sn.cdscaling`.
 *
 * It takes the per-class population vector at one station and returns the
 * per-class rate multipliers, which is the signature `pfqn_cdfun` consumes; the
 * alias resolves to the same std::function type as `pfqn::CdScaling`, so a map
 * built here is passed straight through to the api layer.
 */
template <class T>
using CdScaling = std::function<std::vector<T>(const std::vector<T>&)>;

/**
 * A globally state-dependent scaling, `sn.gdscaling`.
 *
 * Unlike CdScaling it is declared on the NETWORK, not on a station: the argument
 * is the FULL population matrix, given row-major as nstations rows of nclasses
 * entries, and the result is either one scalar, one entry per station, or one
 * entry per (station, class) in the same row-major order. This is the Whittle
 * primitive -- a rate that reads the whole state -- and no per-station scaling
 * can express it when one route holds several resources at once.
 */
template <class T>
using GdScaling = std::function<std::vector<T>(const std::vector<T>&)>;

// ---------------------------------------------------------------------------
// Global constants
// ---------------------------------------------------------------------------

/**
 * The MATLAB GlobalConstants, as reported by lineStart at its defaults.
 *
 * These are doubles on purpose even in the exact instantiation: they are
 * tolerances and sentinels of the ALGORITHM, not quantities of the model, and
 * `lineStart` prints exactly these values. Converting them through
 * num_traits<T>::from_double keeps the exact backend reproducing the same
 * branch decisions as the reference rather than a mathematically cleaner set.
 */
struct GlobalConstants {
    static constexpr double FineTol = 1e-8;
    static constexpr double CoarseTol = 1e-3;
    static constexpr double Zero = 1e-14;
    /** Below this an off-diagonal entry is NO ARC of the phase / state graph. */
    static constexpr double ArcTol = 1e-12;
    /** Rate of an Immediate distribution; its mean is 1/Immediate = 1e-8. */
    static constexpr double Immediate = 1e8;
    /**
     * Stand-in for an unbounded COUNT, MATLAB `GlobalConstants.MaxInt`. Used
     * where a state row must hold a server count and Inf is not a count.
     */
    static constexpr double MaxInt = 2147483647.0;
};

// ---------------------------------------------------------------------------
// Distribution descriptor
// ---------------------------------------------------------------------------

/**
 * A LINE Distribution, as the model layer and `sn` carry it.
 *
 * WHAT EACH CONSUMER READS, which is why all of it is here:
 *   sn.rates, sn.scv     the first two moments -- every AMVA path
 *   sn.procid            the type tag -- the qsys and QNA dispatch
 *   sn.proc, sn.pie      the (D0,D1) pair -- QNA, RQNA, cache, polling
 *   sn.mu, sn.phi        the phase rates and completion probabilities
 *   sn.phases            the order of that representation
 *   sn.lst               the Laplace-Stieltjes transform -- M/G/1 analyzers
 *
 * `params` holds the constructor arguments in MATLAB's getParam order, so a
 * dump can be compared parameter by parameter rather than through the moments,
 * which two different distributions can share.
 *
 * Two values are special and must not be confused, because they enter the
 * struct as opposite extremes:
 *   Immediate  mean = 1/GlobalConstants.Immediate = 1e-8, rate = 1e8
 *   Disabled   rate = NaN, which marks a (station, class) pair the class never
 *              visits; the chain and visit machinery keys on it.
 *
 * ARITHMETIC. The Markovian families are rational in their parameters and are
 * built exactly. Gamma, Weibull and Lognormal are not -- their moments call
 * tgamma or exp -- so their factories refuse by name under exact arithmetic
 * rather than returning a rounded rational that would look exact.
 */
template <class T>
struct PriorSpec;

template <class T>
struct Distrib {
    ProcessType type = ProcessType::DISABLED;
    /**
     * The law as DECLARED, when this one is a surrogate fitted over it.
     *
     * `sn_nonmarkov_toph` installs a fitted (D0,D1) over a Gamma or a Lognormal
     * and retags `type` PH or ME, after which nothing names or evaluates the law
     * the user wrote. MMAP[K]/G[K]/1 reads that law's TRANSFORM rather than the
     * fit, so it needs the original; every other consumer wants the surrogate
     * and reads this struct as before. Null when no substitution has happened.
     * A POINTER, and shared, for the reason `prior` below is one: an inline
     * member would make the type self-embedding.
     */
    std::shared_ptr<Distrib<T>> declared;
    T mean = num_traits<T>::from_int(0);
    T scv = num_traits<T>::from_int(1);
    bool disabled = true;
    /** Constructor arguments, in MATLAB getParam order. */
    std::vector<T> params;
    /** Replayer / Trace samples; empty for every other type. */
    std::vector<T> trace;
    /**
     * The trace FILE a Replayer was read from, when there was one.
     *
     * The samples above are what every solver in this port uses, so the path is
     * carried only for the exporters: `saveServiceStrategy` hands JMT a
     * `ReplayerPar` naming a file, and a Replayer exported without it is a JMT
     * model that reads nothing. Empty when the samples were supplied directly.
     * `network_writer.h` emits it for the same reason: the samples have no wire
     * form, so without the path a Replayer is written back as the moments and
     * reloads as a different law.
     */
    std::string trace_file;
    /**
     * The (D0,D1) pair when the type carries one directly.
     *
     * EMPTY for Det, Uniform, Pareto, Gamma, Weibull, Lognormal and Replayer:
     * MATLAB's getProcess returns their raw PARAMETERS there, and
     * refreshProcessRepresentations replaces them with an Erlang approximation
     * (`convertToMAP`) on the way into sn.proc. That conversion is a property
     * of the refresh, not of the distribution, so it is not done here; see
     * dist_to_map() in lang/distribution.h.
     */
    Matrix<T> D0, D1;
    /** MMAP per-class D1 blocks / BMAP per-batch-size blocks; empty otherwise. */
    std::vector<Matrix<T>> Dmark;
    /**
     * The alternatives of a `Prior`, set only when `type == PRIOR`.
     *
     * A POINTER, and shared: `PriorSpec` holds `Distrib<T>` values, so an
     * inline member would make the type self-embedding, and the spec is
     * immutable once built, so copying a service table copies a pointer rather
     * than a design. `mean` and `scv` beside it are the MIXTURE moments, as
     * MATLAB's `Prior.getMean`/`getSCV` return: a Prior that reaches
     * `refresh_rates` therefore lowers to a rate rather than to a NaN. That is
     * for honesty of the struct dump only -- the featset gate refuses the model
     * before any solver reads those rates, and `dist_to_map`, `dist_lst` and
     * `dist_moment` refuse a Prior by name.
     */
    std::shared_ptr<PriorSpec<T>> prior;
    bool is_prior() const { return type == ProcessType::PRIOR; }

    /**
     * `sn.proc{i}{r} = {breakpoints, A, B, cyclic}` of a MAPt / PHt / NHPP.
     *
     * `sched_bp` has one more entry than there are segments -- it is the
     * BOUNDARY vector, so segment k is in force on [bp(k), bp(k+1)) -- and
     * `sched_D0[k]`, `sched_D1[k]` are the pair of segment k, ALREADY in MAP
     * form. A PHt is stored converted, D0 = S and D1 = (-S e) alpha, because
     * `sn_schedule_nominal` converts it on every read and keeping the raw
     * (alpha, S) here would make every consumer repeat that conversion and one
     * of them eventually forget. The raw form is not needed: the conversion is
     * lossless and nothing downstream asks for alpha again.
     *
     * EMPTY FOR EVERY OTHER TYPE. `has_schedule()` is the test, and a solver
     * with no notion of time simply never calls it -- the nominal pair in
     * `D0`/`D1` is a complete, time-averaged answer for such a solver.
     */
    std::vector<T> sched_bp;
    std::vector<Matrix<T>> sched_D0, sched_D1;
    bool sched_cyclic = false;
    bool has_schedule() const { return !sched_D0.empty(); }

    static Distrib exp_mean(const T& m) {
        const T one = num_traits<T>::from_int(1);
        Distrib d;
        d.type = ProcessType::EXP;
        d.mean = m;
        d.scv = one;
        d.disabled = false;
        const T lambda = m > num_traits<T>::from_int(0) ? T(one / m)
                                                       : num_traits<T>::from_double(
                                                             GlobalConstants::Immediate);
        d.params.push_back(lambda);
        d.D0 = Matrix<T>(1, 1, T(-lambda));
        d.D1 = Matrix<T>(1, 1, lambda);
        return d;
    }
    static Distrib exp_rate(const T& r) {
        const T zero = num_traits<T>::from_int(0);
        if (r <= zero) {
            // Exp.fitRate(0) rationale: see _kb/04-networkstruct.md (cpp port notes)
            return exp_mean(num_traits<T>::from_double(1.0 / GlobalConstants::Zero));
        }
        return exp_mean(T(num_traits<T>::from_int(1) / r));
    }
    /**
     * The Immediate singleton. Its MEAN is zero and its RATE is 1e8, and the
     * two are deliberately not reciprocal: MATLAB's Immediate.getMean() returns
     * 0 while Immediate.getRate() returns GlobalConstants.Immediate, and both
     * are read, by different callers. SolverLN reads the mean (a task with an
     * Immediate think time contributes no think time); Network.refreshRates
     * reads the rate (the station serves the class in 1e-8 time units, not in
     * zero, which would be an infinite service rate the MVA recursion cannot
     * carry). Collapsing them onto 1/mean or 1/rate breaks one caller or the
     * other, so the type tag decides.
     */
    /**
     * The point mass at zero.
     *
     * Its SCV is 1, NOT the 0 of a degenerate distribution. The reference makes
     * this explicit -- `Immediate.getSCV` returns 1 in MATLAB, in the JAR and in
     * Python -- because Immediate is realised downstream as an exponential of
     * rate GlobalConstants.Immediate rather than as a Dirac: `rate()` below
     * returns 1e8, and the SCV has to be the one that goes with it. Declaring 0
     * here is invisible on a layer of PS or infinite-server stations, where the
     * AMVA correction does not read the SCV at all, and shows up only once a
     * layer holds an FCFS or multiserver station -- so it survives a model like
     * lqn_ofbiz and breaks a model like lqn_basic.
     */
    static Distrib immediate() {
        Distrib d;
        d.type = ProcessType::IMMEDIATE;
        d.mean = num_traits<T>::from_int(0);
        d.scv = num_traits<T>::from_int(1);
        d.disabled = false;
        const T imm = num_traits<T>::from_double(GlobalConstants::Immediate);
        d.D0 = Matrix<T>(1, 1, T(-imm));
        d.D1 = Matrix<T>(1, 1, imm);
        return d;
    }
    static Distrib disabled_dist() { return Distrib(); }
    static Distrib det(const T& m) {
        Distrib d;
        d.type = ProcessType::DET;
        d.mean = m;
        d.scv = num_traits<T>::from_int(0);
        d.disabled = false;
        d.params.push_back(m);
        return d;
    }

    // -----------------------------------------------------------------------
    // The Markovian families: (D0,D1) is built here, exactly
    // -----------------------------------------------------------------------

    /** Erlang(alpha, r): r phases of rate alpha, as MATLAB's Erlang(phaseRate, nphases). */
    static Distrib erlang(const T& phase_rate, std::size_t r) {
        const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
        if (r == 0) throw InputError("Erlang: the number of phases must be positive");
        if (!(phase_rate > zero)) throw InputError("Erlang: the phase rate must be positive");
        Distrib d;
        d.type = ProcessType::ERLANG;
        d.disabled = false;
        d.params.push_back(phase_rate);
        d.params.push_back(num_traits<T>::from_int(static_cast<long>(r)));
        d.mean = T(num_traits<T>::from_int(static_cast<long>(r)) / phase_rate);
        d.scv = T(one / num_traits<T>::from_int(static_cast<long>(r)));
        d.D0 = Matrix<T>(r, r, zero);
        d.D1 = Matrix<T>(r, r, zero);
        for (std::size_t i = 0; i < r; ++i) {
            d.D0(i, i) = T(-phase_rate);
            if (i + 1 < r) d.D0(i, i + 1) = phase_rate;
        }
        d.D1(r - 1, 0) = phase_rate;
        return d;
    }

    /**
     * Erlang fitted to a mean and an SCV, as MATLAB's Erlang.fitMeanAndSCV.
     *
     * AN SCV ABOVE ONE IS REFUSED, not answered. An Erlang of order r has
     * SCV = 1/r, so the family reaches 1 and no higher; ceil(1/c2) is 1 for
     * every c2 > 1, and returning that means handing back an EXPONENTIAL under
     * the name of the distribution the caller asked for. MATLAB errors here and
     * this port now does too, so a mis-specified SCV is a diagnostic rather
     * than a silently different service process.
     */
    static Distrib erlang_fit(const T& m, const T& c2) {
        const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
        if (!(c2 > zero)) throw InputError("Erlang.fitMeanAndSCV: the SCV must be positive");
        if (c2 > one)
            throw InputError(
                "Erlang.fitMeanAndSCV: the Erlang distribution requires a squared coefficient "
                "of variation <= 1; use HyperExp, Coxian or APH above 1");
        // MATLAB: r = ceil(1/scv); alpha = r/mean
        const double r_d = std::ceil(1.0 / num_traits<T>::to_double(c2));
        const std::size_t r = static_cast<std::size_t>(r_d < 1.0 ? 1.0 : r_d);
        return erlang(T(num_traits<T>::from_int(static_cast<long>(r)) / m), r);
    }

    /**
     * HyperExp(p, lambda1, lambda2): phase i chosen with probability p_i.
     *
     * D1(i,j) = mu(i) p(j) -- the OUTER product. Associating the other way
     * gives D1(i,j) = mu(i) p(i) replicated across the row, whose rows no
     * longer sum with D0 to zero; _kb records that trap in the MATLAB class.
     */
    /**
     * The same for any number of branches, which is what MATLAB `HyperExp`
     * accepts and what the writers emit as vector `p` and `lambda`. The
     * two-branch entry point stays because it carries MATLAB's getParam order
     * (p, lambda1, lambda2), which a parameter-by-parameter dump compares
     * against.
     */
    static Distrib hyperexp_n(const std::vector<T>& p, const std::vector<T>& lambda) {
        const T zero = num_traits<T>::from_int(0), two = num_traits<T>::from_int(2);
        const std::size_t n = p.size();
        if (n == 0 || lambda.size() != n)
            throw InputError("HyperExp: p and lambda must be non-empty and of equal length");
        Distrib d;
        d.type = ProcessType::HYPEREXP;
        d.disabled = false;
        for (const T& v : p) d.params.push_back(v);
        for (const T& v : lambda) d.params.push_back(v);
        d.D0 = Matrix<T>(n, n, zero);
        d.D1 = Matrix<T>(n, n, zero);
        T m1 = zero, m2 = zero;
        for (std::size_t i = 0; i < n; ++i) {
            if (!(lambda[i] > zero)) throw InputError("HyperExp: the phase rates must be positive");
            d.D0(i, i) = T(-lambda[i]);
            for (std::size_t j = 0; j < n; ++j) d.D1(i, j) = T(lambda[i] * p[j]);
            m1 += T(p[i] / lambda[i]);
            m2 += T(two * p[i] / (lambda[i] * lambda[i]));
        }
        d.mean = m1;
        d.scv = T((m2 - m1 * m1) / (m1 * m1));
        return d;
    }

    static Distrib hyperexp(const T& p, const T& lambda1, const T& lambda2) {
        const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
        if (!(lambda1 > zero) || !(lambda2 > zero))
            throw InputError("HyperExp: the phase rates must be positive");
        if (p < zero || p > one) throw InputError("HyperExp: p is not a probability");
        Distrib d;
        d.type = ProcessType::HYPEREXP;
        d.disabled = false;
        d.params.push_back(p);
        d.params.push_back(lambda1);
        d.params.push_back(lambda2);
        const T q = T(one - p);
        d.mean = T(p / lambda1 + q / lambda2);
        const T m2 = T(num_traits<T>::from_int(2) *
                       (p / (lambda1 * lambda1) + q / (lambda2 * lambda2)));
        d.scv = T((m2 - d.mean * d.mean) / (d.mean * d.mean));
        d.D0 = Matrix<T>(2, 2, zero);
        d.D1 = Matrix<T>(2, 2, zero);
        d.D0(0, 0) = T(-lambda1);
        d.D0(1, 1) = T(-lambda2);
        d.D1(0, 0) = T(lambda1 * p);
        d.D1(0, 1) = T(lambda1 * q);
        d.D1(1, 0) = T(lambda2 * p);
        d.D1(1, 1) = T(lambda2 * q);
        return d;
    }

    /**
     * Coxian(mu, phi): phase i completes with probability phi(i) and otherwise
     * moves to phase i+1. The last phase always completes.
     */
    static Distrib coxian(const std::vector<T>& mu, const std::vector<T>& phi) {
        const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
        const std::size_t n = mu.size();
        if (n == 0 || phi.size() != n)
            throw InputError("Coxian: mu and phi must be non-empty and of equal length");
        Distrib d;
        d.type = n == 2 ? ProcessType::COX2 : ProcessType::COXIAN;
        d.disabled = false;
        for (const T& v : mu) d.params.push_back(v);
        for (const T& v : phi) d.params.push_back(v);
        d.D0 = Matrix<T>(n, n, zero);
        d.D1 = Matrix<T>(n, n, zero);
        for (std::size_t i = 0; i < n; ++i) {
            if (!(mu[i] > zero)) throw InputError("Coxian: the phase rates must be positive");
            d.D0(i, i) = T(-mu[i]);
            if (i + 1 < n) d.D0(i, i + 1) = T(mu[i] * (one - phi[i]));
            d.D1(i, 0) = T(mu[i] * phi[i]);
        }
        // moments of the absorbing chain started in phase 1: m_k = k! e_1 (-D0)^-k e
        d.mean = ph_moment_from(d.D0, 1, 1);
        const T m2 = ph_moment_from(d.D0, 1, 2);
        d.scv = T((m2 - d.mean * d.mean) / (d.mean * d.mean));
        return d;
    }

    /** Cox2(mu1, mu2, phi1), MATLAB's two-phase Coxian constructor. */
    static Distrib cox2(const T& mu1, const T& mu2, const T& phi1) {
        std::vector<T> mu, phi;
        mu.push_back(mu1);
        mu.push_back(mu2);
        phi.push_back(phi1);
        phi.push_back(num_traits<T>::from_int(1));
        return coxian(mu, phi);
    }

    /**
     * PH / APH given by (alpha, A): D0 = A and D1 = (-A e) alpha.
     *
     * `acyclic` selects the type tag only; the representation is the same, and
     * no consumer of sn.proc distinguishes them.
     */
    static Distrib phase_type(const std::vector<T>& alpha, const Matrix<T>& A, bool acyclic) {
        const T zero = num_traits<T>::from_int(0);
        const std::size_t n = alpha.size();
        if (n == 0 || A.rows() != n || A.cols() != n)
            throw InputError("PH: alpha and the subgenerator have inconsistent sizes");
        Distrib d;
        d.type = acyclic ? ProcessType::APH : ProcessType::PH;
        d.disabled = false;
        for (const T& v : alpha) d.params.push_back(v);
        d.D0 = A;
        d.D1 = Matrix<T>(n, n, zero);
        for (std::size_t i = 0; i < n; ++i) {
            T out = zero;
            for (std::size_t j = 0; j < n; ++j) out += A(i, j);
            for (std::size_t j = 0; j < n; ++j) d.D1(i, j) = T(-out * alpha[j]);
        }
        d.mean = ph_moment(alpha, A, 1);
        const T m2 = ph_moment(alpha, A, 2);
        d.scv = T((m2 - d.mean * d.mean) / (d.mean * d.mean));
        return d;
    }

    /** A MAP given by its two matrices; the moments are those of its stationary phase. */
    static Distrib map_dist(const Matrix<T>& D0, const Matrix<T>& D1, ProcessType tag) {
        if (D0.rows() != D0.cols() || D1.rows() != D1.cols() || D0.rows() != D1.rows())
            throw InputError("MAP: D0 and D1 must be square and of the same order");
        Distrib d;
        d.type = tag;
        d.disabled = false;
        d.D0 = D0;
        d.D1 = D1;
        // MAP moment rationale: see _kb/04-networkstruct.md (cpp port notes)
        d.mean = num_traits<T>::from_int(0);
        d.scv = num_traits<T>::from_int(1);
        return d;
    }

    // -----------------------------------------------------------------------
    // The time-inhomogeneous families (Ko-Pender)
    // -----------------------------------------------------------------------

    /** The shared constructor of the three schedule families. */
    static Distrib sched_dist(const std::vector<T>& breakpoints,
                              const std::vector<Matrix<T>>& D0segs,
                              const std::vector<Matrix<T>>& D1segs, bool cyclic, ProcessType tag) {
        const T zero = num_traits<T>::from_int(0);
        const std::string who = process_to_text(tag);
        const std::size_t n = D0segs.size();
        if (n == 0 || D1segs.size() != n)
            throw InputError(who + ": the two segment lists must be non-empty and equally long");
        if (breakpoints.size() != n + 1)
            throw InputError(who +
                             ": breakpoints is the BOUNDARY vector and must hold one more entry "
                             "than there are segments");
        for (std::size_t k = 0; k + 1 < breakpoints.size(); ++k)
            if (!(breakpoints[k + 1] > breakpoints[k]))
                throw InputError(who + ": breakpoints must be strictly increasing");
        const std::size_t order = D0segs[0].rows();
        for (std::size_t k = 0; k < n; ++k)
            if (D0segs[k].rows() != order || D0segs[k].cols() != order ||
                D1segs[k].rows() != order || D1segs[k].cols() != order)
                throw InputError(who +
                                 ": every segment must have the same order; the schedule "
                                 "modulates one phase structure and does not switch between them");

        Distrib d;
        d.type = tag;
        d.disabled = false;
        d.sched_bp = breakpoints;
        d.sched_D0 = D0segs;
        d.sched_D1 = D1segs;
        d.sched_cyclic = cyclic;
        // The nominal pair: the width-weighted time average, i.e. the first two
        // outputs of `sn_schedule_nominal`.
        T total = zero;
        for (std::size_t k = 0; k < n; ++k) total += T(breakpoints[k + 1] - breakpoints[k]);
        d.D0 = Matrix<T>(order, order, zero);
        d.D1 = Matrix<T>(order, order, zero);
        for (std::size_t k = 0; k < n; ++k) {
            const T w = T(T(breakpoints[k + 1] - breakpoints[k]) / total);
            for (std::size_t a = 0; a < order; ++a)
                for (std::size_t b = 0; b < order; ++b) {
                    d.D0(a, b) += w * D0segs[k](a, b);
                    d.D1(a, b) += w * D1segs[k](a, b);
                }
        }
        // Same convention as `map_dist`: the moments of a MAP are computed by
        // `dist_moment`, not stored.
        d.mean = zero;
        d.scv = num_traits<T>::from_int(1);
        return d;
    }

    /**
     * MAPt(breakpoints, {D0_k}, {D1_k}, cyclic): a piecewise-constant
     * (D0(t), D1(t)).
     *
     * `breakpoints` is the boundary vector, so it holds one more entry than
     * there are segments and must be strictly increasing. Every segment must
     * have the SAME ORDER: the schedule modulates one phase structure, it does
     * not switch between structures, and a solver that integrated across a
     * change of order would have no way to map the phase occupancy across the
     * boundary. That is the reference's own constructor requirement.
     *
     * `D0` / `D1` are set to the WIDTH-WEIGHTED TIME AVERAGE of the segments,
     * which is `sn_schedule_nominal`'s nominal pair: it is the stationary
     * carrier of the phase structure the schedule modulates, so the phase count
     * and the mean rate a time-blind consumer reads are the ones the model
     * actually has.
     */
    static Distrib mapt(const std::vector<T>& breakpoints, const std::vector<Matrix<T>>& D0segs,
                        const std::vector<Matrix<T>>& D1segs, bool cyclic) {
        return sched_dist(breakpoints, D0segs, D1segs, cyclic, ProcessType::MAPT);
    }

    /**
     * PHt(breakpoints, {alpha_k}, {S_k}, cyclic), stored as its equivalent MAP
     * schedule: D0 = S and D1 = (-S e) alpha, the pair `sn_schedule_nominal`
     * builds from a PHt slot.
     */
    static Distrib pht(const std::vector<T>& breakpoints, const std::vector<std::vector<T>>& alphas,
                       const std::vector<Matrix<T>>& Ssegs, bool cyclic) {
        const T zero = num_traits<T>::from_int(0);
        if (alphas.size() != Ssegs.size() || alphas.empty())
            throw InputError("PHt: alpha and S must have the same, non-zero number of segments");
        std::vector<Matrix<T>> D0segs, D1segs;
        for (std::size_t k = 0; k < Ssegs.size(); ++k) {
            const Matrix<T>& S = Ssegs[k];
            const std::vector<T>& a = alphas[k];
            if (S.rows() != S.cols() || S.rows() != a.size())
                throw InputError("PHt: alpha and the sub-generator disagree in order");
            Matrix<T> D1(S.rows(), S.cols(), zero);
            for (std::size_t i = 0; i < S.rows(); ++i) {
                T s = zero;
                for (std::size_t j = 0; j < S.cols(); ++j) s += S(i, j);
                for (std::size_t j = 0; j < S.cols(); ++j) D1(i, j) = T(-s * a[j]);
            }
            D0segs.push_back(S);
            D1segs.push_back(D1);
        }
        Distrib d = sched_dist(breakpoints, D0segs, D1segs, cyclic, ProcessType::PHT);
        return d;
    }

    /**
     * NHPP(breakpoints, rates, cyclic): a MAPt of ORDER ONE, which is what an
     * inhomogeneous Poisson process is. Building it through the same path is
     * what makes every schedule consumer see one representation.
     */
    static Distrib nhpp(const std::vector<T>& breakpoints, const std::vector<T>& rates,
                        bool cyclic) {
        std::vector<Matrix<T>> D0segs, D1segs;
        for (std::size_t k = 0; k < rates.size(); ++k) {
            Matrix<T> a(1, 1, T(-rates[k])), b(1, 1, rates[k]);
            D0segs.push_back(a);
            D1segs.push_back(b);
        }
        return sched_dist(breakpoints, D0segs, D1segs, cyclic, ProcessType::NHPP);
    }

    // -----------------------------------------------------------------------
    // The non-Markovian families: parameters only, as MATLAB's getProcess
    // -----------------------------------------------------------------------

    /** Uniform(a, b). */
    static Distrib uniform(const T& a, const T& b) {
        if (!(b > a)) throw InputError("Uniform: the upper bound must exceed the lower bound");
        Distrib d;
        d.type = ProcessType::UNIFORM;
        d.disabled = false;
        d.params.push_back(a);
        d.params.push_back(b);
        d.mean = T((a + b) / num_traits<T>::from_int(2));
        const T w = T(b - a);
        d.scv = T((w * w / num_traits<T>::from_int(12)) / (d.mean * d.mean));
        return d;
    }

    /** Pareto(shape, scale), with the MATLAB parameter order (alpha, k). */
    static Distrib pareto(const T& shape, const T& scale) {
        const T one = num_traits<T>::from_int(1), two = num_traits<T>::from_int(2);
        if (!(shape > two))
            throw InputError("Pareto: the shape must exceed 2 for a finite variance");
        Distrib d;
        d.type = ProcessType::PARETO;
        d.disabled = false;
        d.params.push_back(shape);
        d.params.push_back(scale);
        d.mean = T(shape * scale / (shape - one));
        const T var = T(scale * scale * shape / ((shape - one) * (shape - one)) / (shape - two));
        d.scv = T(var / (d.mean * d.mean));
        return d;
    }

    /**
     * Gamma(shape, scale), Weibull(scale, shape) and Lognormal(mu, sigma).
     *
     * Their moments are values of the gamma function or of exp, so they exist
     * only where the arithmetic has transcendentals. Under exact arithmetic the
     * factory REFUSES rather than storing a rounded rational: a rational that
     * came out of tgamma is not the exact moment of the distribution, and every
     * downstream claim of exactness would be false.
     */
    static Distrib gamma_dist(const T& shape, const T& scale) {
        const T one = num_traits<T>::from_int(1);
        Distrib d;
        d.type = ProcessType::GAMMA;
        d.disabled = false;
        d.params.push_back(shape);
        d.params.push_back(scale);
        d.mean = T(shape * scale);
        d.scv = T(one / shape);
        return d;
    }

    static Distrib weibull(const T& scale, const T& shape) {
        if constexpr (!num_traits<T>::has_transcendental) {
            throw UnsupportedError(
                "Weibull: its moments are values of the gamma function, which exact arithmetic "
                "has no representation for; use the double or real backend");
        } else {
            const double a = num_traits<T>::to_double(scale);
            const double r = num_traits<T>::to_double(shape);
            if (!(a > 0.0) || !(r > 0.0))
                throw InputError("Weibull: the scale and the shape must be positive");
            const double g1 = std::tgamma(1.0 + 1.0 / r);
            const double g2 = std::tgamma(1.0 + 2.0 / r);
            Distrib d;
            d.type = ProcessType::WEIBULL;
            d.disabled = false;
            d.params.push_back(scale);
            d.params.push_back(shape);
            d.mean = num_traits<T>::from_double(a * g1);
            d.scv = num_traits<T>::from_double((g2 - g1 * g1) / (g1 * g1));
            return d;
        }
    }

    static Distrib lognormal(const T& logmean, const T& logsigma) {
        if constexpr (!num_traits<T>::has_transcendental) {
            throw UnsupportedError(
                "Lognormal: its moments are values of exp, which exact arithmetic has no "
                "representation for; use the double or real backend");
        } else {
            const double mu = num_traits<T>::to_double(logmean);
            const double sg = num_traits<T>::to_double(logsigma);
            if (!(sg > 0.0)) throw InputError("Lognormal: sigma must be positive");
            Distrib d;
            d.type = ProcessType::LOGNORMAL;
            d.disabled = false;
            d.params.push_back(logmean);
            d.params.push_back(logsigma);
            d.mean = num_traits<T>::from_double(std::exp(mu + sg * sg / 2.0));
            d.scv = num_traits<T>::from_double(std::exp(sg * sg) - 1.0);
            return d;
        }
    }

    /**
     * `Normal(mu, sigma)`: the Gaussian, for use as a continuous `Prior`'s
     * parameter density.
     *
     * `scv` is Inf at mu = 0, which is `Normal.m:52-63`'s own answer and not a
     * degradation: the SCV of a zero-mean law is not defined, and the reference
     * says so rather than dividing.
     */
    static Distrib normal(const T& mu, const T& sigma) {
        if (!(num_traits<T>::to_double(sigma) > 0.0))
            throw InputError("Normal: sigma must be positive");
        Distrib d;
        d.type = ProcessType::NORMAL;
        d.disabled = false;
        d.params.push_back(mu);
        d.params.push_back(sigma);
        d.mean = mu;
        const double m = num_traits<T>::to_double(mu), s = num_traits<T>::to_double(sigma);
        d.scv = std::abs(m) < GlobalConstants::FineTol
                    ? num_traits<T>::from_double(std::numeric_limits<double>::infinity())
                    : num_traits<T>::from_double((s * s) / (m * m));
        return d;
    }

    /**
     * Replayer / Trace read FROM A FILE, which keeps the path beside the samples.
     *
     * Every solver in this port replays `trace`; the path is what the EXPORTERS
     * need. JMT is handed a `ReplayerPar` naming a file and has no way to take
     * samples inline, so a Replayer exported without it is a JMT model that
     * reads nothing -- `oqn_trace_driven` was refused outright for exactly that.
     */
    static Distrib replayer_from(const std::vector<T>& samples, const std::string& path) {
        Distrib d = replayer(samples);
        d.trace_file = path;
        return d;
    }

    /** Replayer / Trace: the samples, with their empirical first two moments. */
    static Distrib replayer(const std::vector<T>& samples) {
        if (samples.empty()) throw InputError("Replayer: the trace is empty");
        Distrib d;
        d.type = ProcessType::REPLAYER;
        d.disabled = false;
        d.trace = samples;
        const T n = num_traits<T>::from_int(static_cast<long>(samples.size()));
        T s1 = num_traits<T>::from_int(0), s2 = num_traits<T>::from_int(0);
        for (const T& x : samples) {
            s1 += x;
            s2 += T(x * x);
        }
        d.mean = T(s1 / n);
        d.scv = T((s2 / n - d.mean * d.mean) / (d.mean * d.mean));
        return d;
    }

    // -----------------------------------------------------------------------
    // The discrete families
    //
    // MATLAB's getProcess for each of these returns [mean, scv] and nothing
    // else, so they carry NO (D0,D1): refreshProcessRepresentations replaces
    // them with the Erlang fit convertToMAP, which `dist_to_map` reproduces
    // from the two moments alone. Storing an invented representation here
    // would make sn.proc disagree with the reference for the same model.
    // -----------------------------------------------------------------------

    /** DiscreteUniform(a, b) over the integers a..b inclusive. */
    static Distrib discrete_uniform(const T& a, const T& b) {
        const T two = num_traits<T>::from_int(2), one = num_traits<T>::from_int(1);
        if (!(b >= a)) throw InputError("DiscreteUniform: the upper bound must not be below the lower");
        Distrib d;
        d.type = ProcessType::DUNIFORM;
        d.disabled = false;
        d.params.push_back(a);
        d.params.push_back(b);
        d.mean = T((a + b) / two);
        const T w = T(b - a + one);
        const T var = T((w * w - one) / num_traits<T>::from_int(12));
        d.scv = T(var / (d.mean * d.mean));
        return d;
    }

    /** Bernoulli(p): one trial, mean p and variance p(1-p). */
    static Distrib bernoulli(const T& p) {
        const T one = num_traits<T>::from_int(1);
        Distrib d;
        d.type = ProcessType::BERNOULLI;
        d.disabled = false;
        d.params.push_back(p);
        d.mean = p;
        d.scv = T((one - p) / p);
        return d;
    }

    /** Binomial(n, p). */
    static Distrib binomial(const T& n, const T& p) {
        const T one = num_traits<T>::from_int(1);
        Distrib d;
        d.type = ProcessType::BINOMIAL;
        d.disabled = false;
        d.params.push_back(n);
        d.params.push_back(p);
        d.mean = T(n * p);
        d.scv = T((one - p) / (n * p));
        return d;
    }

    /**
     * Poisson(lambda), whose SCV is 1/lambda -- the count's variance is lambda
     * and its mean is lambda, so this is NOT the exponential's SCV of 1.
     */
    static Distrib poisson(const T& lambda) {
        const T one = num_traits<T>::from_int(1);
        Distrib d;
        d.type = ProcessType::POISSON;
        d.disabled = false;
        d.params.push_back(lambda);
        d.mean = lambda;
        d.scv = T(one / lambda);
        return d;
    }

    /**
     * Geometric(p) on the MATLAB convention: the NUMBER OF TRIALS to the first
     * success, support {1, 2, ...}, so the mean is 1/p and the SCV is 1-p.
     */
    static Distrib geometric(const T& p) {
        const T one = num_traits<T>::from_int(1);
        Distrib d;
        d.type = ProcessType::GEOMETRIC;
        d.disabled = false;
        d.params.push_back(p);
        d.mean = T(one / p);
        d.scv = T(one - p);
        return d;
    }

    /**
     * Zipf(s, n) over the ranks 1..n, with the generalized harmonic moments
     * H(s-1,n)/H(s,n) and H(s-2,n)/H(s,n) MATLAB `Zipf.m` uses.
     *
     * The harmonic sums call pow for a non-integer shape, so the factory is
     * gated on transcendental arithmetic exactly as Weibull and Lognormal are.
     */
    static Distrib zipf(const T& s, std::size_t n) {
        if constexpr (!num_traits<T>::has_transcendental) {
            (void)s;
            (void)n;
            throw UnsupportedError(
                "Zipf: its moments are generalized harmonic sums of a real exponent, which exact "
                "arithmetic has no representation for; use the double or real backend");
        } else {
            if (n == 0) throw InputError("Zipf: the item count must be positive");
            const double sv = num_traits<T>::to_double(s);
            auto harmonic = [n](double e) {
                double acc = 0.0;
                for (std::size_t k = 1; k <= n; ++k) acc += std::pow(double(k), -e);
                return acc;
            };
            const double h0 = harmonic(sv), h1 = harmonic(sv - 1.0), h2 = harmonic(sv - 2.0);
            Distrib d;
            d.type = ProcessType::ZIPF;
            d.disabled = false;
            d.params.push_back(s);
            d.params.push_back(num_traits<T>::from_int(static_cast<long>(n)));
            const double m1 = h1 / h0;
            d.mean = num_traits<T>::from_double(m1);
            d.scv = num_traits<T>::from_double((h2 / h0 - m1 * m1) / (m1 * m1));
            return d;
        }
    }

    /**
     * DiscreteSampler(p, x): the pmf p over the points x.
     *
     * THE MOMENTS ARE TAKEN OVER x. This port, and MATLAB
     * `DiscreteSampler.getMean` with it, used to weight by the RANKS 1..n
     * instead; the two agree on the default x = 1:n, which is the form the
     * cache popularity vectors are written in, so the rank form survived
     * unnoticed until a fork's jobs-per-link distribution arrived on a shifted
     * support. The JAR and native Python already weighted by x.
     */
    static Distrib discrete_sampler(const std::vector<T>& p, const std::vector<T>& x) {
        if (p.empty()) throw InputError("DiscreteSampler: the probability vector is empty");
        if (!x.empty() && x.size() != p.size())
            throw InputError("DiscreteSampler: p and x must have the same length");
        Distrib d;
        d.type = ProcessType::DISCRETESAMPLER;
        d.disabled = false;
        d.params = p;
        d.trace = x;
        T m1 = num_traits<T>::from_int(0), m2 = num_traits<T>::from_int(0);
        T tot = num_traits<T>::from_int(0);
        for (std::size_t k = 0; k < p.size(); ++k) {
            // an absent x is the default support 1..n
            const T pt = x.empty() ? num_traits<T>::from_int(static_cast<long>(k + 1)) : x[k];
            m1 += T(p[k] * pt);
            m2 += T(p[k] * pt * pt);
            tot += p[k];
        }
        m1 = T(m1 / tot);
        m2 = T(m2 / tot);
        d.mean = m1;
        d.scv = T((m2 - m1 * m1) / (m1 * m1));
        return d;
    }

    /**
     * EmpiricalCDF(x, F): the moments of the MIDPOINT rule over the CDF bins,
     * which is what MATLAB `EmpiricalCDF.getMoments` integrates -- each bin
     * contributes its midpoint raised to the moment order, weighted by the CDF
     * increment. The rows are the (F, x) pairs in the order they arrive.
     */
    static Distrib empirical_cdf(const std::vector<T>& x, const std::vector<T>& F) {
        if (x.size() != F.size() || x.size() < 2)
            throw InputError("EmpiricalCDF: x and F must be equally long and hold at least two points");
        const T two = num_traits<T>::from_int(2);
        Distrib d;
        d.type = ProcessType::EMPIRICALCDF;
        d.disabled = false;
        d.trace = x;
        d.params = F;
        T m1 = num_traits<T>::from_int(0), m2 = num_traits<T>::from_int(0);
        for (std::size_t i = 0; i + 1 < x.size(); ++i) {
            const T mid = T((x[i + 1] - x[i]) / two + x[i]);
            const T w = T(F[i + 1] - F[i]);
            m1 += T(mid * w);
            m2 += T(mid * mid * w);
        }
        d.mean = m1;
        d.scv = T(m2 / (m1 * m1) - num_traits<T>::from_int(1));
        return d;
    }

    // -----------------------------------------------------------------------
    // The matrix-exponential and discrete-time Markovian families
    // -----------------------------------------------------------------------

    /**
     * ME(alpha, A): the matrix-exponential distribution, whose moments are the
     * phase-type ones -- k! alpha (-A)^-k e -- evaluated by DEFINITION rather
     * than through the stationary vector of A + (-Ae)alpha. MATLAB `ME.getMean`
     * makes the same choice and says why: the stationary solve is a
     * probabilistic object that a non-Markovian A degrades badly (a CME of
     * order 101 lost 2.6e-4 in the SCV that way, against 1e-13 by definition).
     */
    static Distrib me(const std::vector<T>& alpha, const Matrix<T>& A) {
        Distrib d = phase_type(alpha, A, false);
        d.type = ProcessType::ME;
        return d;
    }

    /** RAP(H0, H1): a rational arrival process, whose moments are the MAP ones. */
    static Distrib rap(const Matrix<T>& H0, const Matrix<T>& H1) {
        return map_dist(H0, H1, ProcessType::RAP);
    }

    /**
     * DMAP(D0, D1): a DISCRETE-time MAP, where D0 + D1 is stochastic rather
     * than a generator. Its moments cannot come from `dist_refresh_moments`'s
     * continuous formulas; `dmap_moments` in lang/distribution.h fills them.
     */
    static Distrib dmap(const Matrix<T>& D0, const Matrix<T>& D1) {
        return map_dist(D0, D1, ProcessType::DMAP);
    }

    /**
     * MMAP: D0 plus one D1 block per mark. `D1` is their sum, the aggregate
     * arrival matrix every unmarked consumer reads, and the blocks stay in
     * `Dmark` for the ones that distinguish marks.
     */
    static Distrib mmap(const Matrix<T>& D0, const std::vector<Matrix<T>>& D1k) {
        if (D1k.empty()) throw InputError("MMAP: no marked arrival block was given");
        Matrix<T> agg(D0.rows(), D0.cols(), num_traits<T>::from_int(0));
        for (const Matrix<T>& Dk : D1k) {
            if (Dk.rows() != D0.rows() || Dk.cols() != D0.cols())
                throw InputError("MMAP: every marked block must have the order of D0");
            for (std::size_t i = 0; i < agg.rows(); ++i)
                for (std::size_t j = 0; j < agg.cols(); ++j) agg(i, j) += Dk(i, j);
        }
        Distrib d = map_dist(D0, agg, ProcessType::MMAP);
        d.Dmark = D1k;
        return d;
    }

    /**
     * BMAP: the batch-size blocks D0, D1, ..., Dk, where Dj carries an arrival
     * of batch size j. The wire form is the whole list including D0, so the
     * head is split off here.
     */
    static Distrib bmap(const std::vector<Matrix<T>>& D) {
        if (D.size() < 2) throw InputError("BMAP: the block list must carry D0 and at least one batch block");
        const std::vector<Matrix<T>> batches(D.begin() + 1, D.end());
        Distrib d = mmap(D[0], batches);
        d.type = ProcessType::BMAP;
        return d;
    }

    bool is_immediate() const { return type == ProcessType::IMMEDIATE; }

    /** True when the type carries a (D0,D1) pair of its own. */
    bool has_map() const { return D0.rows() > 0; }

    /**
     * Phase rates, MATLAB's getMu: the total outgoing rate of each phase.
     *
     * Empty when the type carries no representation, which is what
     * refreshProcessPhases writes as NaN for a Fork or a Join.
     */
    std::vector<T> mu_vec() const {
        std::vector<T> v;
        if (!has_map()) {
            if (disabled) return v;
            v.push_back(rate());  // one phase at the mean rate, as MATLAB does
            return v;
        }
        for (std::size_t i = 0; i < D0.rows(); ++i) v.push_back(T(-D0(i, i)));
        return v;
    }

    /** Completion probabilities, MATLAB's getPhi: (D1 e) ./ (-diag(D0)). */
    std::vector<T> phi_vec() const {
        std::vector<T> v;
        if (!has_map()) {
            if (disabled) return v;
            v.push_back(num_traits<T>::from_int(1));
            return v;
        }
        const T zero = num_traits<T>::from_int(0);
        for (std::size_t i = 0; i < D1.rows(); ++i) {
            T s = zero;
            for (std::size_t j = 0; j < D1.cols(); ++j) s += D1(i, j);
            const T out = T(-D0(i, i));
            v.push_back(out == zero ? num_traits<T>::from_int(1) : T(s / out));
        }
        return v;
    }

    /** The order of the representation, MATLAB's sn.phases. */
    std::size_t phases() const {
        if (disabled) return 0;
        return has_map() ? D0.rows() : 1;
    }

    /**
     * The k-th raw moment of a phase-type (alpha, A): k! alpha (-A)^-k e.
     *
     * The inverse is never formed: the powers are accumulated by repeated
     * solves of (-A) x = b, which is exact in rational arithmetic and stable
     * in floating point.
     */
    static T ph_moment(const std::vector<T>& alpha, const Matrix<T>& A, unsigned k) {
        const std::size_t n = alpha.size();
        std::vector<T> x(n, num_traits<T>::from_int(1));
        for (unsigned i = 0; i < k; ++i) x = solve_neg(A, x);
        T acc = num_traits<T>::from_int(0);
        for (std::size_t i = 0; i < n; ++i) acc += alpha[i] * x[i];
        T fact = num_traits<T>::from_int(1);
        for (unsigned i = 2; i <= k; ++i) fact *= num_traits<T>::from_int(static_cast<long>(i));
        return T(fact * acc);
    }

    /** The same, for a representation entered in a single phase (1-based). */
    static T ph_moment_from(const Matrix<T>& A, std::size_t start, unsigned k) {
        std::vector<T> alpha(A.rows(), num_traits<T>::from_int(0));
        alpha[start - 1] = num_traits<T>::from_int(1);
        return ph_moment(alpha, A, k);
    }

  private:
    /** Solve (-A) x = b by Gaussian elimination with partial pivoting. */
    static std::vector<T> solve_neg(const Matrix<T>& A, const std::vector<T>& b) {
        const T zero = num_traits<T>::from_int(0);
        const std::size_t n = A.rows();
        if (A.cols() != n || b.size() != n)
            throw InputError("phase-type moment: the subgenerator is not square");
        Matrix<T> M(n, n, zero);
        std::vector<T> x = b;
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = 0; j < n; ++j) M(i, j) = T(-A(i, j));
        for (std::size_t col = 0; col < n; ++col) {
            std::size_t best = col;
            double bv = std::fabs(num_traits<T>::to_double(M(col, col)));
            for (std::size_t r = col + 1; r < n; ++r) {
                const double v = std::fabs(num_traits<T>::to_double(M(r, col)));
                if (v > bv) {
                    bv = v;
                    best = r;
                }
            }
            if (best != col) {
                for (std::size_t j = 0; j < n; ++j) std::swap(M(col, j), M(best, j));
                std::swap(x[col], x[best]);
            }
            if (M(col, col) == zero)
                throw NumericError("phase-type moment: the subgenerator is singular");
            for (std::size_t r = 0; r < n; ++r) {
                if (r == col) continue;
                const T f = T(M(r, col) / M(col, col));
                if (f == zero) continue;
                for (std::size_t j = 0; j < n; ++j) M(r, j) = T(M(r, j) - f * M(col, j));
                x[r] = T(x[r] - f * x[col]);
            }
        }
        for (std::size_t i = 0; i < n; ++i) x[i] = T(x[i] / M(i, i));
        return x;
    }

  public:


    /**
     * The rate MATLAB's refreshRates would store: 1/mean, with the Immediate
     * singleton short-circuited to its declared rate so that the reciprocal of
     * 1e-8 is exactly 1e8 in every arithmetic rather than 1e8 plus rounding.
     */
    T rate() const {
        if (type == ProcessType::IMMEDIATE) return num_traits<T>::from_double(GlobalConstants::Immediate);
        const T zero = num_traits<T>::from_int(0);
        if (mean <= zero) return num_traits<T>::from_double(GlobalConstants::Immediate);
        return T(num_traits<T>::from_int(1) / mean);
    }
};

/**
 * What a `Prior` carries, in either of its two forms.
 *
 * DISCRETE: an explicit set of alternative distributions and their prior
 * weights, which must sum to one. This is the form the model.json wire carries
 * (`{"type":"Prior","distributions":[...],"probabilities":[...]}`), because a
 * factory cannot cross JSON.
 *
 * CONTINUOUS: a density over a scalar parameter theta plus a map theta ->
 * Distribution, the form the epistemic propagation of Trivedi and Bobbio
 * (2017), Sec. 3.4 needs. It is reduced to the discrete form by
 * `prior_discretize` (lang/prior.h) before anything downstream sees it, so both
 * forms are consumed identically. It can only be BUILT programmatically.
 *
 * IT IS NOT A MIXTURE. Each alternative is a separate model realization whose
 * weight is a prior probability over models, not a branching probability inside
 * one model. The mixture moments are still computed (see `Distrib::prior`)
 * because MATLAB's `Prior.getMean`/`getSCV` do, but they are a summary of the
 * epistemic uncertainty and not the law any station serves.
 */
template <class T>
struct PriorSpec {
    /** True for the parameter-density form, false for the alternative-set form. */
    bool continuous = false;
    /** The alternatives and their weights; the discrete form only. */
    std::vector<Distrib<T>> alternatives;
    std::vector<T> probabilities;
    /** The law of the scalar parameter; the continuous form only. */
    Distrib<T> param_dist;
    /**
     * theta -> Distribution; the continuous form only.
     *
     * A `std::function` and not a serializable description, exactly as MATLAB's
     * `distFactory` is a function handle: the map is arbitrary code (a rate
     * becomes an Exp, a scale becomes an Erlang of fixed order), and no wire
     * format in this codebase encodes it. That is why the JSON reader builds
     * the discrete form only.
     */
    std::function<Distrib<T>(const T&)> factory;
};

}  // namespace lang
}  // namespace line

#endif  // LINE_LANG_LANG_TYPES_H
