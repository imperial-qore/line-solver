/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * The Environment model.json reader, on the bytes the reference writers emit.
 *
 * The JSON below is `save_model(env, ...)`'s output for
 * `matlab/examples/advanced/randomEnv/renv_basic.m` (native Python writer;
 * `environment2json` in `linemodel_save.m` writes the same keys), pasted
 * verbatim rather than regenerated, so a change to the reader that stopped
 * understanding the wire format fails here instead of in a parity sweep.
 *
 * What it asserts, in order: the environment the reader builds is the one the
 * file describes (stage names, stage networks, arc rates, and the zero-based
 * from/to convention BOTH writers use); solving it reproduces the MATLAB
 * ENV(FLD) means already pinned in `test_env.cpp`, which is what says the
 * reader hands the solver the same environment the programmatic API does; and
 * every construct the reader cannot honour is REFUSED BY NAME.
 */

#include <string>

#include "doctest.h"
#include "line/io/environment_reader.h"
#include "line/solvers/env/solver_env.h"

using namespace line;

namespace {

const char* kRenvBasic = R"JSON(
{
 "format": "line-model",
 "version": "1.0",
 "model": {
  "type": "Environment",
  "name": "ServerModes",
  "numStages": 2,
  "stages": [
   {"name": "Fast",
    "model": {"type": "Network", "name": "Fast",
     "nodes": [
      {"name": "ThinkTime", "type": "Delay", "scheduling": "INF",
       "service": {"Jobs": {"type": "Exp", "params": {"lambda": 1.0}}}},
      {"name": "Fast/Slow Server", "type": "Queue", "scheduling": "FCFS",
       "service": {"Jobs": {"type": "Exp", "params": {"lambda": 4.0}}}}],
     "classes": [{"name": "Jobs", "type": "Closed", "population": 5, "refNode": "ThinkTime"}],
     "routing": {"type": "matrix", "matrix": {"Jobs,Jobs": {
        "ThinkTime": {"Fast/Slow Server": 1.0},
        "Fast/Slow Server": {"ThinkTime": 1.0}}}}}},
   {"name": "Slow",
    "model": {"type": "Network", "name": "Slow",
     "nodes": [
      {"name": "ThinkTime", "type": "Delay", "scheduling": "INF",
       "service": {"Jobs": {"type": "Exp", "params": {"lambda": 1.0}}}},
      {"name": "Fast/Slow Server", "type": "Queue", "scheduling": "FCFS",
       "service": {"Jobs": {"type": "Exp", "params": {"lambda": 1.0}}}}],
     "classes": [{"name": "Jobs", "type": "Closed", "population": 5, "refNode": "ThinkTime"}],
     "routing": {"type": "matrix", "matrix": {"Jobs,Jobs": {
        "ThinkTime": {"Fast/Slow Server": 1.0},
        "Fast/Slow Server": {"ThinkTime": 1.0}}}}}}],
  "transitions": [
   {"from": 0, "to": 1, "distribution": {"type": "Exp", "params": {"lambda": 0.5}}},
   {"from": 1, "to": 0, "distribution": {"type": "Exp", "params": {"lambda": 1.0}}}]
 }
}
)JSON";

env::Environment<double> read_str(const std::string& text) {
    io::detail::json root = io::detail::json::parse(text);
    return io::build_environment_from_json<double>(root);
}

/** The renv_basic bytes with one key replaced, for the refusal cases. */
std::string with(const std::string& needle, const std::string& repl) {
    std::string s(kRenvBasic);
    const std::string::size_type at = s.find(needle);
    REQUIRE(at != std::string::npos);
    return s.replace(at, needle.size(), repl);
}

}  // namespace

TEST_CASE("environment_reader builds renv_basic from the writers' model.json") {
    env::Environment<double> e = read_str(kRenvBasic);

    REQUIRE(e.nstages() == 2);
    CHECK(e.name() == "ServerModes");
    CHECK(e.stage(0).name == "Fast");
    CHECK(e.stage(1).name == "Slow");

    // The stage networks: same structure, different server rate. Station 0 is
    // the Delay and station 1 the Queue, in declaration order.
    const qn::NetworkStruct<double>& fast = e.stage(0).model;
    const qn::NetworkStruct<double>& slow = e.stage(1).model;
    REQUIRE(fast.nstations == 2);
    REQUIRE(fast.nclasses == 1);
    CHECK(fast.rates(1, 0) == doctest::Approx(4.0));
    CHECK(slow.rates(1, 0) == doctest::Approx(1.0));
    CHECK(fast.rates(0, 0) == doctest::Approx(1.0));
    CHECK(fast.njobs()[0] == doctest::Approx(5.0));

    // The arcs, and the ZERO-BASED wire index: Fast -> Slow at 0.5 and
    // Slow -> Fast at 1.0. Reading them as 1-based would enable the two SELF
    // loops instead and leave the environment with no way to switch.
    CHECK(e.arc(0, 1).enabled);
    CHECK(e.arc(1, 0).enabled);
    CHECK_FALSE(e.arc(0, 0).enabled);
    CHECK_FALSE(e.arc(1, 1).enabled);
    CHECK(e.arc(0, 1).dist.mean == doctest::Approx(2.0));
    CHECK(e.arc(1, 0).dist.mean == doctest::Approx(1.0));

    // `init()` is the caller's, not the reader's; once run it reports the
    // holding times and stage probabilities MATLAB's env.init() does.
    e.init();
    CHECK(e.prob_env[0] == doctest::Approx(0.666666666667).epsilon(1e-9));
    CHECK(e.prob_env[1] == doctest::Approx(0.333333333333).epsilon(1e-9));
}

TEST_CASE("a read environment solves to the MATLAB ENV(FLD) means") {
    env::Environment<double> e = read_str(kRenvBasic);
    env::EnvOptions o;
    o.iter_max = 50;
    o.iter_tol = 0.01;
    o.timespan_end = 100.0;
    o.tran_points = 2001;
    env::SolverEnv<double> s(e, o);
    const env::EnvSolution sol = s.solve();

    // Same reference values, and same quadrature tolerance, as test_env.cpp's
    // programmatically built environment: MATLAB ENV(FLD).getAvg reports
    // QN = [2.99611792893, 2.00388207107].
    REQUIRE(sol.converged);
    CHECK(sol.QN(0, 0) == doctest::Approx(2.99611792893).epsilon(1e-2));
    CHECK(sol.QN(1, 0) == doctest::Approx(2.00388207107).epsilon(1e-2));
    CHECK(sol.QN(0, 0) + sol.QN(1, 0) == doctest::Approx(5.0).epsilon(1e-9));
}

TEST_CASE("environment_reader refuses what it cannot honour, by name") {
    // A Network envelope is not an Environment: -s env would otherwise build a
    // one-stage environment out of a model that never described one.
    CHECK_THROWS_AS(read_str(R"({"model": {"type": "Network", "name": "n"}})"), UnsupportedError);

    // A nodeFailures entry with no distributions describes nothing: the block
    // is the record of a breakdown, and a breakdown is its two times.
    CHECK_THROWS_AS(read_str(with("\"transitions\": [",
                                  "\"nodeFailures\": [{\"node\": \"Fast/Slow Server\"}],\n"
                                  "  \"transitions\": [")),
                    InputError);

    // An unknown key is a constraint this reader would otherwise drop.
    CHECK_THROWS_AS(read_str(with("\"numStages\": 2", "\"numStages\": 2, \"envState\": 1")),
                    UnsupportedError);

    // numStages disagreeing with the stage array means the transitions are
    // indexed against a count the file does not carry.
    CHECK_THROWS_AS(read_str(with("\"numStages\": 2", "\"numStages\": 3")), InputError);

    // A transition naming a stage that does not exist.
    CHECK_THROWS_AS(read_str(with("{\"from\": 1, \"to\": 0", "{\"from\": 2, \"to\": 0")),
                    InputError);

    // A stage with no network: ENV solves a network per stage, and an absent
    // one would surface as an empty drift rather than as a missing input.
    CHECK_THROWS_AS(read_str(with("{\"name\": \"Slow\",\n    \"model\":",
                                  "{\"name\": \"Slow\", \"unusedModel\":")),
                    UnsupportedError);
}

/**
 * A THREE-stage environment with a non-exponential arc, written by MATLAB.
 *
 * These are `linemodel_save`'s own bytes for the model `renv3_ref.m` builds: a
 * ring Fast -> Medium -> Slow -> Fast whose middle arc is an Erlang(2), so the
 * stage exit metrics integrate against a CDF with no memoryless shortcut, and
 * whose stage networks differ only in the server rate (4, 2, 1). MATLAB's own
 * writer is used here, against the Python writer's bytes above, because the two
 * are the only producers of this format and a reader that understood one but
 * not the other would be found by nothing else.
 */
static const char* kRenv3 = R"JSON(
{
 "format": "line-model",
 "version": "1.0",
 "model": {
  "name": "ThreeModes",
  "numStages": 3,
  "stages": [
   {
    "model": {
     "classes": [
      {
       "name": "Jobs",
       "population": 5,
       "refNode": "ThinkTime",
       "type": "Closed"
      }
     ],
     "name": "Fast",
     "nodes": [
      {
       "name": "ThinkTime",
       "scheduling": "INF",
       "service": {
        "Jobs": {
         "params": {
          "lambda": 1
         },
         "type": "Exp"
        }
       },
       "type": "Delay"
      },
      {
       "name": "Server",
       "scheduling": "FCFS",
       "service": {
        "Jobs": {
         "params": {
          "lambda": 4
         },
         "type": "Exp"
        }
       },
       "type": "Queue"
      }
     ],
     "routing": {
      "matrix": {
       "Jobs,Jobs": {
        "Server": {
         "ThinkTime": 1
        },
        "ThinkTime": {
         "Server": 1
        }
       }
      },
      "type": "matrix"
     },
     "type": "Network"
    },
    "name": "Fast"
   },
   {
    "model": {
     "classes": [
      {
       "name": "Jobs",
       "population": 5,
       "refNode": "ThinkTime",
       "type": "Closed"
      }
     ],
     "name": "Medium",
     "nodes": [
      {
       "name": "ThinkTime",
       "scheduling": "INF",
       "service": {
        "Jobs": {
         "params": {
          "lambda": 1
         },
         "type": "Exp"
        }
       },
       "type": "Delay"
      },
      {
       "name": "Server",
       "scheduling": "FCFS",
       "service": {
        "Jobs": {
         "params": {
          "lambda": 2
         },
         "type": "Exp"
        }
       },
       "type": "Queue"
      }
     ],
     "routing": {
      "matrix": {
       "Jobs,Jobs": {
        "Server": {
         "ThinkTime": 1
        },
        "ThinkTime": {
         "Server": 1
        }
       }
      },
      "type": "matrix"
     },
     "type": "Network"
    },
    "name": "Medium"
   },
   {
    "model": {
     "classes": [
      {
       "name": "Jobs",
       "population": 5,
       "refNode": "ThinkTime",
       "type": "Closed"
      }
     ],
     "name": "Slow",
     "nodes": [
      {
       "name": "ThinkTime",
       "scheduling": "INF",
       "service": {
        "Jobs": {
         "params": {
          "lambda": 1
         },
         "type": "Exp"
        }
       },
       "type": "Delay"
      },
      {
       "name": "Server",
       "scheduling": "FCFS",
       "service": {
        "Jobs": {
         "params": {
          "lambda": 1
         },
         "type": "Exp"
        }
       },
       "type": "Queue"
      }
     ],
     "routing": {
      "matrix": {
       "Jobs,Jobs": {
        "Server": {
         "ThinkTime": 1
        },
        "ThinkTime": {
         "Server": 1
        }
       }
      },
      "type": "matrix"
     },
     "type": "Network"
    },
    "name": "Slow"
   }
  ],
  "transitions": [
   {
    "distribution": {
     "params": {
      "lambda": 0.5
     },
     "type": "Exp"
    },
    "from": 0,
    "to": 1
   },
   {
    "distribution": {
     "params": {
      "k": 2,
      "lambda": 2
     },
     "type": "Erlang"
    },
    "from": 1,
    "to": 2
   },
   {
    "distribution": {
     "params": {
      "lambda": 1
     },
     "type": "Exp"
    },
    "from": 2,
    "to": 0
   }
  ],
  "type": "Environment"
 }
}
)JSON";

TEST_CASE("environment_reader: three stages, an Erlang arc, against MATLAB ENV(FLD)") {
    env::Environment<double> e = read_str(kRenv3);
    REQUIRE(e.nstages() == 3);
    CHECK(e.stage(0).name == "Fast");
    CHECK(e.stage(1).name == "Medium");
    CHECK(e.stage(2).name == "Slow");
    // The ring, and the Erlang(2) with mean 1 on its middle arc.
    CHECK(e.arc(0, 1).enabled);
    CHECK(e.arc(1, 2).enabled);
    CHECK(e.arc(2, 0).enabled);
    CHECK(e.arc(1, 2).dist.mean == doctest::Approx(1.0));
    CHECK(e.stage(0).model.rates(1, 0) == doctest::Approx(4.0));
    CHECK(e.stage(1).model.rates(1, 0) == doctest::Approx(2.0));
    CHECK(e.stage(2).model.rates(1, 0) == doctest::Approx(1.0));

    env::EnvOptions o;
    o.iter_max = 50;
    o.iter_tol = 0.01;
    o.timespan_end = 100.0;
    o.tran_points = 2001;
    env::SolverEnv<double> s(e, o);
    const env::EnvSolution sol = s.solve();
    REQUIRE(sol.converged);

    // Exact against MATLAB, and independent of the quadrature: the holding
    // means are 2, 1, 1 around a deterministic ring, so the environment sits in
    // Fast half the time; and every stage saturates its server, so the blended
    // throughput is 0.5*4 + 0.25*2 + 0.25*1 exactly.
    CHECK(e.prob_env[0] == doctest::Approx(0.5).epsilon(1e-9));
    CHECK(e.prob_env[1] == doctest::Approx(0.25).epsilon(1e-9));
    CHECK(e.prob_env[2] == doctest::Approx(0.25).epsilon(1e-9));
    CHECK(sol.UN(1, 0) == doctest::Approx(1.0).epsilon(1e-6));
    CHECK(sol.TN(1, 0) == doctest::Approx(2.75).epsilon(1e-6));
    CHECK(sol.QN(0, 0) + sol.QN(1, 0) == doctest::Approx(5.0).epsilon(1e-9));

    // The quadrature-sensitive half. MATLAB ENV(FLD).getAvg on these bytes
    // reports QN = [2.704598317158, 2.295401682842]; this port at 2001 points
    // reports 2.711508. Refining the grid moves it 2.705057 -> 2.711508 ->
    // 2.713747 -> 2.714195 for 501 -> 2001 -> 10001 -> 40001, i.e. to a limit
    // just ABOVE the reference, which sums the same rule over LSODA's coarser
    // adaptive grid. The 4e-3 below is that difference and nothing else. Note
    // that BOTH codebases report Tput 2.70 at the Delay against 2.75 at the
    // Queue, a flow-balance violation the mean-field blend has in the REFERENCE
    // (see git show 449847e7b:_kb/log.md, the exact ENV oracle), so neither
    // number is exact here.
    CHECK(sol.QN(0, 0) == doctest::Approx(2.704598317158).epsilon(4e-3));
    CHECK(sol.QN(1, 0) == doctest::Approx(2.295401682842).epsilon(4e-3));
}

/**
 * The `nodeFailures` block, in the MACRO form a hand-written file takes.
 *
 * These are the bytes of `gallery_renv_breakdown`: an M/M/1 whose Server breaks
 * down at rate 0.1, is repaired at rate 1.0, and serves at 0.5 while it is
 * down. One stage is declared -- the base model -- and the block is what turns
 * it into the two stages the reference's `addNodeFailureRepair` builds.
 */
static const char* kRenvBreakdownMacro = R"JSON(
{
 "format": "line-model",
 "version": "1.0",
 "model": {
  "type": "Environment",
  "name": "ServerEnv",
  "numStages": 1,
  "stages": [
   {"name": "Base",
    "model": {"type": "Network", "name": "ServerWithFailures",
     "nodes": [
      {"name": "Arrivals", "type": "Source",
       "service": {"Jobs": {"type": "Exp", "params": {"lambda": 0.8}}}},
      {"name": "Server", "type": "Queue", "scheduling": "FCFS", "servers": 1,
       "service": {"Jobs": {"type": "Exp", "params": {"lambda": 2.0}}}},
      {"name": "Departures", "type": "Sink"}],
     "classes": [{"name": "Jobs", "type": "Open"}],
     "routing": {"type": "matrix", "matrix": {"Jobs,Jobs": {
        "Arrivals": {"Server": 1.0},
        "Server": {"Departures": 1.0}}}}}}],
  "nodeFailures": [
   {"node": "Server",
    "breakdownRate": {"type": "Exp", "params": {"lambda": 0.1}},
    "repairRate": {"type": "Exp", "params": {"lambda": 1.0}},
    "downService": {"type": "Exp", "params": {"lambda": 0.5}},
    "breakdownResetPolicy": "keep",
    "repairResetPolicy": "keep"}]
 }
}
)JSON";

TEST_CASE("environment_reader expands a nodeFailures macro into UP and DOWN stages") {
    env::Environment<double> e = read_str(kRenvBreakdownMacro);

    // The stage names are the reference's and not the declared one: the base
    // stage becomes `UP` whatever the file called it, because `addNodeRepair`
    // finds the stage to repair into by that name.
    REQUIRE(e.nstages() == 2);
    CHECK(e.stage(0).name == "UP");
    CHECK(e.stage(0).type == "operational");
    CHECK(e.stage(1).name == "DOWN_Server");
    CHECK(e.stage(1).type == "failed");

    // Station 0 is the Source and station 1 the Queue; the Sink holds no jobs
    // and is not a station. Only the Server's service differs between stages.
    const qn::NetworkStruct<double>& up = e.stage(0).model;
    const qn::NetworkStruct<double>& down = e.stage(1).model;
    REQUIRE(up.nstations == 2);
    REQUIRE(down.nstations == 2);
    CHECK(up.rates(1, 0) == doctest::Approx(2.0));
    CHECK(down.rates(1, 0) == doctest::Approx(0.5));
    CHECK(up.rates(0, 0) == doctest::Approx(0.8));
    CHECK(down.rates(0, 0) == doctest::Approx(0.8));

    // The two implied arcs, and no others.
    CHECK(e.arc(0, 1).enabled);
    CHECK(e.arc(1, 0).enabled);
    CHECK_FALSE(e.arc(0, 0).enabled);
    CHECK_FALSE(e.arc(1, 1).enabled);
    CHECK(e.arc(0, 1).dist.mean == doctest::Approx(10.0));
    CHECK(e.arc(1, 0).dist.mean == doctest::Approx(1.0));

    // `keep` is the identity, which this port spells as no reset at all.
    CHECK_FALSE(static_cast<bool>(e.arc(0, 1).reset));
    CHECK_FALSE(static_cast<bool>(e.arc(1, 0).reset));

    // The descriptor survives the read, which is what lets the environment be
    // written back as the model it came from.
    REQUIRE(e.node_failures().size() == 1);
    CHECK(e.node_failures()[0].node == "Server");
    CHECK(e.node_failures()[0].has_repair);
    CHECK(e.node_failures()[0].breakdown_reset == "keep");
    CHECK(e.node_failures()[0].down_service.mean == doctest::Approx(2.0));

    // Availability: MTTF 10 against MTTR 1, so the server is up 10/11 of the
    // time. This is the reliability reading of probEnv and is exact.
    e.init();
    CHECK(e.prob_env[0] == doctest::Approx(10.0 / 11.0).epsilon(1e-9));
    CHECK(e.prob_env[1] == doctest::Approx(1.0 / 11.0).epsilon(1e-9));

    // `getReliabilityTable` on the same environment, which is what
    // renv_node_breakdown prints: MTTF 10, MTTR 1, MTBF 11, availability 10/11.
    const env::Environment<double>::Reliability rel = e.reliability();
    CHECK(rel.mttf == doctest::Approx(10.0).epsilon(1e-9));
    CHECK(rel.mttr == doctest::Approx(1.0).epsilon(1e-9));
    CHECK(rel.mtbf == doctest::Approx(11.0).epsilon(1e-9));
    CHECK(rel.availability == doctest::Approx(10.0 / 11.0).epsilon(1e-9));
}

TEST_CASE("environment_reader: a clear policy empties the queues on the switch") {
    // Same model, but the breakdown flushes the buffer. `clear` is one of the
    // two policies the wire can carry, and it is the one that changes numbers.
    std::string s(kRenvBreakdownMacro);
    const std::string::size_type at = s.find("\"breakdownResetPolicy\": \"keep\"");
    REQUIRE(at != std::string::npos);
    s = s.replace(at, std::string("\"breakdownResetPolicy\": \"keep\"").size(),
                  "\"breakdownResetPolicy\": \"clear\"");
    env::Environment<double> e = read_str(s);

    REQUIRE(static_cast<bool>(e.arc(0, 1).reset));
    Matrix<double> q(2, 1, 3.0);
    const Matrix<double> after = e.arc(0, 1).reset(q);
    CHECK(after.rows() == 2);
    CHECK(after(0, 0) == doctest::Approx(0.0));
    CHECK(after(1, 0) == doctest::Approx(0.0));
    // The repair arc kept its own policy, which is the identity.
    CHECK_FALSE(static_cast<bool>(e.arc(1, 0).reset));
    CHECK(e.node_failures()[0].breakdown_reset == "clear");
}

/**
 * The EXPANDED form, which is what both writers actually emit: the UP and
 * DOWN_<node> stages and both arcs are declared, and the block adds only what
 * they cannot express. Nothing is expanded a second time.
 */
static const char* kRenvBreakdownExpanded = R"JSON(
{
 "format": "line-model",
 "version": "1.0",
 "model": {
  "type": "Environment",
  "name": "ServerEnv",
  "numStages": 2,
  "stages": [
   {"name": "UP",
    "model": {"type": "Network", "name": "ServerWithFailures",
     "nodes": [
      {"name": "Arrivals", "type": "Source",
       "service": {"Jobs": {"type": "Exp", "params": {"lambda": 0.8}}}},
      {"name": "Server", "type": "Queue", "scheduling": "FCFS", "servers": 1,
       "service": {"Jobs": {"type": "Exp", "params": {"lambda": 2.0}}}},
      {"name": "Departures", "type": "Sink"}],
     "classes": [{"name": "Jobs", "type": "Open"}],
     "routing": {"type": "matrix", "matrix": {"Jobs,Jobs": {
        "Arrivals": {"Server": 1.0},
        "Server": {"Departures": 1.0}}}}}},
   {"name": "DOWN_Server",
    "model": {"type": "Network", "name": "ServerWithFailures",
     "nodes": [
      {"name": "Arrivals", "type": "Source",
       "service": {"Jobs": {"type": "Exp", "params": {"lambda": 0.8}}}},
      {"name": "Server", "type": "Queue", "scheduling": "FCFS", "servers": 1,
       "service": {"Jobs": {"type": "Exp", "params": {"lambda": 0.5}}}},
      {"name": "Departures", "type": "Sink"}],
     "classes": [{"name": "Jobs", "type": "Open"}],
     "routing": {"type": "matrix", "matrix": {"Jobs,Jobs": {
        "Arrivals": {"Server": 1.0},
        "Server": {"Departures": 1.0}}}}}}],
  "transitions": [
   {"from": 0, "to": 1, "distribution": {"type": "Exp", "params": {"lambda": 0.1}}},
   {"from": 1, "to": 0, "distribution": {"type": "Exp", "params": {"lambda": 1.0}}}],
  "nodeFailures": [
   {"node": "Server",
    "breakdownRate": {"type": "Exp", "params": {"lambda": 0.1}},
    "repairRate": {"type": "Exp", "params": {"lambda": 1.0}},
    "downService": {"type": "Exp", "params": {"lambda": 0.5}},
    "breakdownResetPolicy": "clear",
    "repairResetPolicy": "keep"}]
 }
}
)JSON";

TEST_CASE("environment_reader: nodeFailures on expanded stages carries only the policies") {
    env::Environment<double> e = read_str(kRenvBreakdownExpanded);
    REQUIRE(e.nstages() == 2);
    CHECK(e.stage(0).name == "UP");
    CHECK(e.stage(1).name == "DOWN_Server");
    CHECK(e.stage(0).model.rates(1, 0) == doctest::Approx(2.0));
    CHECK(e.stage(1).model.rates(1, 0) == doctest::Approx(0.5));

    // The policies, which are the only thing the block contributed here: the
    // breakdown flushes, the repair carries the jobs over.
    REQUIRE(static_cast<bool>(e.arc(0, 1).reset));
    CHECK_FALSE(static_cast<bool>(e.arc(1, 0).reset));
    REQUIRE(e.node_failures().size() == 1);
    CHECK(e.node_failures()[0].node == "Server");
    CHECK(e.node_failures()[0].breakdown_reset == "clear");
    CHECK(e.node_failures()[0].repair_reset == "keep");
    CHECK(e.node_failures()[0].has_repair);

    // The two stage networks are the file's own and were not rebuilt from the
    // block: an expansion here would have overwritten the DOWN stage with one
    // derived from the UP model, which on this file is the same network but on
    // any file whose DOWN stage differs by more than the service rate is not.
    CHECK(e.stage(1).model.nstations == 2);
    CHECK(e.stage(1).model.rates(0, 0) == doctest::Approx(0.8));
}

TEST_CASE("environment_reader: the DOWN stage NAME is what selects the block's role") {
    // The expanded form with the DOWN stage renamed. No `DOWN_Server` stage is
    // declared any more, so the block reads as a MACRO -- and the macro form
    // expands ONE base stage, which this file is not. The name is load bearing
    // in all three readers, and this is what happens when it does not match.
    std::string s(kRenvBreakdownExpanded);
    const std::string::size_type at = s.find("\"name\": \"DOWN_Server\"");
    REQUIRE(at != std::string::npos);
    s = s.replace(at, std::string("\"name\": \"DOWN_Server\"").size(), "\"name\": \"Degraded\"");
    CHECK_THROWS_AS(read_str(s), InputError);
}

TEST_CASE("environment_reader refuses a nodeFailures block it cannot expand") {
    auto macro_with = [](const std::string& needle, const std::string& repl) {
        std::string s(kRenvBreakdownMacro);
        const std::string::size_type at = s.find(needle);
        REQUIRE(at != std::string::npos);
        return s.replace(at, needle.size(), repl);
    };

    // The macro form IMPLIES the transitions, so a file that also declares them
    // is describing the environment twice and the two could disagree.
    CHECK_THROWS_AS(read_str(macro_with("\"nodeFailures\": [",
                                        "\"transitions\": [{\"from\": 0, \"to\": 0, "
                                        "\"distribution\": {\"type\": \"Exp\", \"params\": "
                                        "{\"lambda\": 1.0}}}],\n  \"nodeFailures\": [")),
                    InputError);

    // A macro block expands ONE base stage; a second declared stage that is no
    // DOWN stage is a stage the expansion has no place for.
    CHECK_THROWS_AS(
        read_str(macro_with(
            "{\"name\": \"Base\",",
            "{\"name\": \"Other\",\n"
            "    \"model\": {\"type\": \"Network\", \"name\": \"Other\",\n"
            "     \"nodes\": [{\"name\": \"D\", \"type\": \"Delay\", \"scheduling\": \"INF\",\n"
            "       \"service\": {\"J\": {\"type\": \"Exp\", \"params\": {\"lambda\": 1.0}}}}],\n"
            "     \"classes\": [{\"name\": \"J\", \"type\": \"Closed\", \"population\": 1, "
            "\"refNode\": \"D\"}],\n"
            "     \"routing\": {\"type\": \"matrix\", \"matrix\": {\"J,J\": {\"D\": {\"D\": "
            "1.0}}}}}},\n"
            "   {\"name\": \"Base\",")),
        InputError);

    // An unknown key at envelope level, alongside the block: the same refusal
    // as on the plain path, and it must not be softened by the expansion.
    CHECK_THROWS_AS(read_str(macro_with("\"nodeFailures\": [",
                                        "\"envState\": 0,\n  \"nodeFailures\": [")),
                    UnsupportedError);

    // `custom` is what both writers print when the policy is a function handle,
    // and they omit the key rather than emit it; a file carrying it is claiming
    // a policy that no file can carry.
    CHECK_THROWS_AS(read_str(macro_with("\"breakdownResetPolicy\": \"keep\"",
                                        "\"breakdownResetPolicy\": \"custom\"")),
                    InputError);

    // A breakdown of a node the base model does not have.
    CHECK_THROWS_AS(read_str(macro_with("\"node\": \"Server\"", "\"node\": \"Serv3r\"")),
                    InputError);

    // A breakdown of the Sink, which is not a station and has no service.
    CHECK_THROWS_AS(read_str(macro_with("\"node\": \"Server\"", "\"node\": \"Departures\"")),
                    InputError);

    // An unknown key inside an entry: the same rule as everywhere else on this
    // path, because a dropped one would change the model silently.
    CHECK_THROWS_AS(read_str(macro_with("\"node\": \"Server\"",
                                        "\"node\": \"Server\", \"downCapacity\": 2")),
                    UnsupportedError);

    // The same node twice: a node has one DOWN stage, and two entries would
    // silently keep whichever was read last.
    CHECK_THROWS_AS(
        read_str(macro_with("\"nodeFailures\": [",
                            "\"nodeFailures\": [{\"node\": \"Server\", \"breakdownRate\": "
                            "{\"type\": \"Exp\", \"params\": {\"lambda\": 0.2}}, \"downService\": "
                            "{\"type\": \"Exp\", \"params\": {\"lambda\": 0.3}}},")),
        InputError);
}

TEST_CASE("a node-breakdown environment read from JSON solves") {
    env::Environment<double> e = read_str(kRenvBreakdownMacro);
    env::EnvOptions o;
    o.iter_max = 50;
    o.iter_tol = 0.01;
    o.timespan_end = 100.0;
    o.tran_points = 2001;
    env::SolverEnv<double> s(e, o);
    const env::EnvSolution sol = s.solve();
    REQUIRE(sol.converged);

    // The Source is offered 0.8 jobs per time unit in both stages, so its
    // throughput is that whatever the environment does.
    CHECK(sol.TN(0, 0) == doctest::Approx(0.8).epsilon(1e-6));

    // The UP stage alone is a stable open queue whose fluid fixed point is
    // 0.8 / 2 = 0.4 jobs; the DOWN stage is offered 0.8 against a service rate
    // of 0.5 and builds up. The blend must therefore sit ABOVE the UP stage's
    // own steady state -- which is the whole reason to model the breakdown
    // rather than to solve the UP model and be done.
    REQUIRE(sol.QExit.size() == 2);
    CHECK(sol.QN(1, 0) > 0.4);
    CHECK(sol.QExit[1](1, 0) > 0.0);
}
