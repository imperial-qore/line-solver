/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * `ldes`: the native LDES engine behind the interface of `jline.cli.LdesCLI`.
 *
 * WHAT THIS BINARY IS FOR. `common/ldes` has been a GraalVM image of the JAVA
 * engine; this program is its replacement, built from
 * `cpp/include/line/solvers/ldes/`. Every client of that binary -- MATLAB's
 * `@SolverLDES/solveCli.m`, native Python's `wrappers/solver_ldes`, and the C++
 * client `wrappers/ldes/solver_ldes.h` -- speaks to it through exactly two
 * artefacts: the ARGUMENT VECTOR and the `ldes-result` DOCUMENT. Neither may
 * drift. A client is not recompiled when the engine changes hands, so an
 * argument this program rejects and the Java one accepted, or a key it spells
 * differently, is a silent breakage in three codebases at once.
 *
 * THE INTERFACE IS THEREFORE COPIED, not designed:
 *
 *   ldes solve <model.json> -o <result.json> [flags]
 *
 * with the flag set of `LdesCLI.handleSolveCommand` and the document shape of
 * `LDESResultIO.write`. Flags this engine cannot honour are REFUSED BY NAME
 * rather than ignored, which is the one deliberate difference from the AOT
 * image it replaces: that image ignored flags it predated, and an ignored
 * `--slotlength` silently ran a continuous-time simulation under a slotted
 * model's name. A refusal is visible; a silent default is not.
 *
 * WHAT IS NOT HERE. `--rest`/server mode belongs to `line_cli` and is refused
 * with the reason rather than accepted and dropped.
 */

#include <cmath>
#include <cstdio>
#include <ctime>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "line/io/network_reader.h"
#include "line/solvers/ldes/ldes_engine.h"
#include "line/solvers/wrappers/ldes/ldes_options.h"
#include "line/util/error.h"

namespace {

using line::ldes::LdesOptions;
using line::ldes::LdesResult;

void print_help() {
    std::cout
        << "LDES: discrete-event simulation engine for LINE models.\n\n"
        << "Usage: ldes solve <model.json> [-o <result.json>] [options]\n\n"
        << "Options:\n"
        << "  -s, --samples N          service completions to simulate (default 200000)\n"
        << "  -e, --maxevents N        event budget; overrides --samples\n"
        << "      --maxtime SECONDS    wall-clock budget\n"
        << "      --seed N             run seed (default 23000; -1 draws one)\n"
        << "      --method NAME        engine method (default 'default')\n"
        << "      --cnvgon             stop on the convergence check\n"
        << "      --cnvgtol X          relative tolerance of that check\n"
        << "      --cnvgbatch N        batches before the first check\n"
        << "      --cnvgchk N          events between checks\n"
        << "      --tranfilter NAME    mser5 | fixed | none\n"
        << "      --mserbatch N        MSER batch size\n"
        << "      --warmupfrac X       warmup fraction, for --tranfilter fixed\n"
        << "      --cimethod NAME      obm | bm | spectral | none\n"
        << "      --obmoverlap X       overlap of the batch means\n"
        << "      --ciminbatch N       batches below which no CI is reported\n"
        << "      --ciminobs N         observations below which no CI is reported\n"
        << "      --spectrallowfreqfrac X   low-frequency fraction of the spectral CI\n"
        << "      --slotted            run on the slot lattice\n"
        << "      --slotlength X       slot length (implies --slotted)\n"
        << "      --replications N     independent replications\n"
        << "      --numthreads N       threads across replications\n"
        << "      --timespan T0,T1     transient run over [T0,T1]\n"
        << "      --busyperiod K       busy-period orders 1..K\n"
        << "      --busyperiod-subnet i,j,...   an extra busy-period target\n"
        << "      --initsol V          initial state, station-major, comma separated\n"
        << "      --export-histogram   fill the state histogram block\n"
        << "      --trajectory         fill the state trajectory block\n"
        << "      --respt-samples      record every per-visit response time\n"
        << "  -h, --help               this text\n";
}

[[noreturn]] void fail(const std::string& msg) {
    std::cerr << "Error: " << msg << "\n";
    std::exit(1);
}

/** A flag the reference accepts and this engine cannot honour: say so. */
[[noreturn]] void refuse(const std::string& flag, const std::string& why) {
    std::cerr << "Error: " << flag << " is accepted by the Java engine but not by this one: "
              << why << ".\nRun the model through common/ldes.jar for it, rather than "
                        "taking a silently different simulation.\n";
    std::exit(1);
}

std::string need_value(int argc, char** argv, int& i, const std::string& flag) {
    if (i + 1 >= argc) fail(flag + " requires a value.");
    return std::string(argv[++i]);
}

std::vector<double> parse_doubles(const std::string& csv) {
    std::vector<double> out;
    std::stringstream ss(csv);
    std::string tok;
    while (std::getline(ss, tok, ',')) {
        if (!tok.empty()) out.push_back(std::atof(tok.c_str()));
    }
    return out;
}

std::vector<std::size_t> parse_indices(const std::string& csv) {
    std::vector<std::size_t> out;
    std::stringstream ss(csv);
    std::string tok;
    while (std::getline(ss, tok, ',')) {
        if (!tok.empty()) out.push_back(static_cast<std::size_t>(std::atol(tok.c_str())));
    }
    return out;
}

/** A JSON number the way the reference writes it: full round-trip precision. */
std::string num(double v) {
    if (std::isnan(v)) return "null";
    if (std::isinf(v)) return v > 0 ? "1e999" : "-1e999";
    char buf[40];
    std::snprintf(buf, sizeof(buf), "%.17g", v);
    return std::string(buf);
}

/** A (rows x cols) matrix as the nested array `LDESResultIO` emits. */
std::string mat_json(const line::Matrix<double>& m, std::size_t rows, std::size_t cols) {
    std::string s = "[";
    for (std::size_t i = 0; i < rows; ++i) {
        if (i) s += ", ";
        s += "[";
        for (std::size_t j = 0; j < cols; ++j) {
            if (j) s += ", ";
            const bool inside = (i < m.rows() && j < m.cols());
            s += num(inside ? m(i, j) : 0.0);
        }
        s += "]";
    }
    s += "]";
    return s;
}

/** The per-visit response times as the nested `[station][class][sample]` array. */
std::string respt_samples_json(const LdesResult& r) {
    if (r.respTimeSamples.empty()) return std::string();
    std::string s = "[";
    for (std::size_t i = 0; i < r.respTimeSamples.size(); ++i) {
        if (i) s += ", ";
        s += "[";
        for (std::size_t k = 0; k < r.respTimeSamples[i].size(); ++k) {
            if (k) s += ", ";
            s += "[";
            for (std::size_t j = 0; j < r.respTimeSamples[i][k].size(); ++j) {
                if (j) s += ", ";
                s += num(r.respTimeSamples[i][k][j]);
            }
            s += "]";
        }
        s += "]";
    }
    s += "]";
    return s;
}

std::string str_json(const std::string& s) {
    std::string out = "\"";
    for (char c : s) {
        if (c == '"' || c == '\\') {
            out += '\\';
            out += c;
        } else if (c == '\n') {
            out += "\\n";
        } else {
            out += c;
        }
    }
    out += "\"";
    return out;
}

/**
 * The `ldes-result` document, key for key with `LDESResultIO.write`.
 *
 * The order of the keys is the reference's, and the blocks a run did not fill
 * are still emitted with their zero matrices: every client indexes this
 * document by key and by station row, and a missing block reads as a parse
 * failure rather than as an absent measurement.
 */
std::string result_json(const LdesResult& r, const LdesOptions& o, double runtime,
                        std::size_t events) {
    const std::size_t S = r.nstations, C = r.nclasses;
    std::string s;
    s += "{\n";
    s += "  \"format\": \"ldes-result\",\n";
    s += "  \"version\": \"1.0\",\n";
    s += "  \"solver\": \"SolverLDES\",\n";
    s += "  \"method\": " + str_json(o.method) + ",\n";
    s += "  \"runtime\": " + num(runtime) + ",\n";
    s += "  \"converged\": false,\n";
    s += "  \"stoppingReason\": \"max_events\",\n";
    s += "  \"convergenceBatches\": 0,\n";
    s += "  \"totalSimulatedEvents\": " + std::to_string(events) + ",\n";

    s += "  \"dimensions\": {\"nstations\": " + std::to_string(S) +
         ", \"nclasses\": " + std::to_string(C) + ", \"nchains\": " + std::to_string(r.nchains) +
         ", \"stationNames\": [";
    for (std::size_t i = 0; i < r.station_names.size(); ++i) {
        if (i) s += ", ";
        s += str_json(r.station_names[i]);
    }
    s += "], \"classNames\": [";
    for (std::size_t i = 0; i < r.class_names.size(); ++i) {
        if (i) s += ", ";
        s += str_json(r.class_names[i]);
    }
    s += "]},\n";

    s += "  \"metrics\": {";
    s += "\"QN\": " + mat_json(r.QN, S, C);
    s += ", \"UN\": " + mat_json(r.UN, S, C);
    s += ", \"RN\": " + mat_json(r.RN, S, C);
    s += ", \"TN\": " + mat_json(r.TN, S, C);
    s += ", \"AN\": " + mat_json(r.AN, S, C);
    s += ", \"WN\": " + mat_json(r.WN, S, C);
    s += ", \"CN\": " + mat_json(r.CN, 1, C);
    s += ", \"XN\": " + mat_json(r.XN, 1, C);
    // The sibling-drop rate at a Join, where LossRate = ArvR - Tput does NOT
    // hold: ArvR counts the SIBLINGS offered and Tput the PARENT jobs released.
    // `LDESResultIO.write` puts it inside `metrics`, so it goes here and not in
    // an optional block of its own.
    if (!r.DropRateJoin.empty()) s += ", \"DropRateJoin\": " + mat_json(r.DropRateJoin, S, C);
    s += "},\n";

    // THE OPTIONAL BLOCKS, each written ONLY when the run measured it, exactly
    // as `LDESResultIO.write` gates them. A block emitted empty is not the same
    // as an absent one: `ldes_prob_from_histogram` and `getAvgBusyPeriod` both
    // treat presence as the claim that the measurement was taken.
    // GATED ON THE FLAG, not on the measurement. A transient run computes the
    // histogram whether or not anyone asked for it, and `LDESResultIO` writes
    // the block only when the flag was given; emitting it anyway would make
    // this engine's document differ from the reference's on every transient run.
    if (o.export_histogram && r.histogram_space.rows() > 0 && r.histogram_time.rows() > 0) {
        s += "  \"stateHistogram\": {";
        s += "\"space\": " + mat_json(r.histogram_space, r.histogram_space.rows(),
                                      r.histogram_space.cols());
        s += ", \"time\": " + mat_json(r.histogram_time, r.histogram_time.rows(), 1);
        if (r.traj_space.rows() > 0 && r.traj_time.rows() > 0) {
            s += ", \"trajSpace\": " +
                 mat_json(r.traj_space, r.traj_space.rows(), r.traj_space.cols());
            s += ", \"trajTime\": " + mat_json(r.traj_time, r.traj_time.rows(), 1);
        }
        s += "},\n";
    }

    if (!r.busy_periods.empty()) {
        s += "  \"busyPeriods\": {\"orders\": " +
             std::to_string(r.busy_periods[0].mean.size()) + ", \"targets\": [";
        for (std::size_t ti = 0; ti < r.busy_periods.size(); ++ti) {
            const LdesResult::BusyPeriodTarget& t = r.busy_periods[ti];
            if (ti) s += ", ";
            s += "{\"name\": " + str_json(t.name) + ", \"stations\": [";
            for (std::size_t k = 0; k < t.stations.size(); ++k) {
                if (k) s += ", ";
                s += std::to_string(t.stations[k]);
            }
            s += "], \"class\": " + std::to_string(t.job_class) + ", \"mean\": [";
            for (std::size_t k = 0; k < t.mean.size(); ++k) {
                if (k) s += ", ";
                s += num(t.mean[k]);
            }
            s += "], \"count\": [";
            for (std::size_t k = 0; k < t.count.size(); ++k) {
                if (k) s += ", ";
                s += num(t.count[k]);
            }
            s += "]}";
        }
        s += "]},\n";
    }

    if (!r.cache_metrics.empty()) {
        s += "  \"cacheMetrics\": {";
        bool first = true;
        for (std::map<std::string, line::ldes::LdesCacheMetrics>::const_iterator it =
                 r.cache_metrics.begin();
             it != r.cache_metrics.end(); ++it) {
            if (!first) s += ", ";
            first = false;
            const line::ldes::LdesCacheMetrics& cm = it->second;
            s += str_json(it->first) + ": {";
            s += "\"hit\": " + mat_json(cm.hit, cm.hit.rows(), cm.hit.cols());
            s += ", \"delayed\": " + mat_json(cm.delayed, cm.delayed.rows(), cm.delayed.cols());
            s += ", \"miss\": " + mat_json(cm.miss, cm.miss.rows(), cm.miss.cols());
            s += ", \"latency\": " + mat_json(cm.latency, cm.latency.rows(), cm.latency.cols());
            s += ", \"hitList\": " + mat_json(cm.hitList, cm.hitList.rows(), cm.hitList.cols());
            s += ", \"itemProb\": " + mat_json(cm.itemProb, cm.itemProb.rows(), cm.itemProb.cols());
            s += ", \"listCost\": " + mat_json(cm.listCost, cm.listCost.rows(), cm.listCost.cols());
            s += "}";
        }
        s += "},\n";
    }

    if (r.QNfcr.rows() > 0) {
        const std::size_t G = r.QNfcr.rows();
        s += "  \"fcr\": {\"nregions\": " + std::to_string(G);
        s += ", \"QNfcr\": " + mat_json(r.QNfcr, G, C);
        s += ", \"UNfcr\": " + mat_json(r.UNfcr, G, C);
        s += ", \"RNfcr\": " + mat_json(r.RNfcr, G, C);
        s += ", \"TNfcr\": " + mat_json(r.TNfcr, G, C);
        s += ", \"ANfcr\": " + mat_json(r.ANfcr, G, C);
        s += ", \"WNfcr\": " + mat_json(r.WNfcr, G, C);
        s += ", \"WeightNfcr\": " + mat_json(r.WeightNfcr, G, C);
        s += ", \"MemOccNfcr\": " + mat_json(r.MemOccNfcr, G, C);
        s += ", \"DropRateNfcr\": " + mat_json(r.DropRateNfcr, G, C);
        s += "},\n";
    }

    // The reference gates the whole transient block on `--trajectory`, not on
    // the run being transient.
    // The samples appear at TOP LEVEL under --respt-samples and INSIDE the
    // transient block under --trajectory. Both spellings are the reference's,
    // and both are read: the C++ client takes whichever is present, MATLAB's
    // `sample()` and `sampleSys()` read the nested one and returned nothing
    // while only the top-level spelling existed.
    const std::string respt_json = respt_samples_json(r);
    if (o.export_trajectory && !r.t.empty()) {
        s += "  \"transient\": {\"t\": [";
        for (std::size_t k = 0; k < r.t.size(); ++k) {
            if (k) s += ", ";
            s += num(r.t[k]);
        }
        s += "]";
        const char* keys[3] = {"QNt", "UNt", "TNt"};
        const std::vector<std::vector<line::Matrix<double>>>* src[3] = {&r.QNt, &r.UNt, &r.TNt};
        for (int q = 0; q < 3; ++q) {
            s += std::string(", \"") + keys[q] + "\": [";
            for (std::size_t i = 0; i < src[q]->size(); ++i) {
                if (i) s += ", ";
                s += "[";
                for (std::size_t k = 0; k < (*src[q])[i].size(); ++k) {
                    if (k) s += ", ";
                    const line::Matrix<double>& mm = (*src[q])[i][k];
                    s += mat_json(mm, mm.rows(), mm.cols());
                }
                s += "]";
            }
            s += "]";
        }
        if (!respt_json.empty()) s += ", \"respTimeSamples\": " + respt_json;
        s += "},\n";
    }

    if (o.export_respt && !respt_json.empty())
        s += "  \"respTimeSamples\": " + respt_json + ",\n";
    s += "  \"sampleCounts\": {";
    s += "\"QNSamples\": " + mat_json(r.QNSamples, S, C);
    s += ", \"UNSamples\": " + mat_json(r.UNSamples, S, C);
    s += ", \"RNSamples\": " + mat_json(r.RNSamples, S, C);
    s += ", \"TNSamples\": " + mat_json(r.TNSamples, S, C);
    s += "},\n";

    s += "  \"confidenceIntervals\": {";
    s += "\"QNCI\": " + mat_json(r.QNCI, S, C);
    s += ", \"UNCI\": " + mat_json(r.UNCI, S, C);
    s += ", \"RNCI\": " + mat_json(r.RNCI, S, C);
    s += ", \"TNCI\": " + mat_json(r.TNCI, S, C);
    s += ", \"ANCI\": " + mat_json(r.ANCI, S, C);
    s += ", \"WNCI\": " + mat_json(r.WNCI, S, C);
    s += "},\n";

    s += "  \"relativePrecision\": {";
    s += "\"QNRelPrec\": " + mat_json(r.QNRelPrec, S, C);
    s += ", \"UNRelPrec\": " + mat_json(r.UNRelPrec, S, C);
    s += ", \"RNRelPrec\": " + mat_json(r.RNRelPrec, S, C);
    s += ", \"TNRelPrec\": " + mat_json(r.TNRelPrec, S, C);
    s += "},\n";

    s += "  \"impatience\": {";
    s += "\"renegedCustomers\": " + mat_json(r.renegedCustomers, S, C);
    s += ", \"avgRenegingWaitTime\": " + mat_json(r.avgRenegingWaitTime, S, C);
    s += ", \"renegingRate\": " + mat_json(r.renegingRate, S, C);
    s += ", \"balkedCustomers\": " + mat_json(r.balkedCustomers, S, C);
    s += ", \"balkingProbability\": " + mat_json(r.balkingProbability, S, C);
    s += ", \"retriedCustomers\": " + mat_json(r.retriedCustomers, S, C);
    s += ", \"retrialDropped\": " + mat_json(r.retrialDropped, S, C);
    s += ", \"avgOrbitSize\": " + mat_json(r.avgOrbitSize, S, C);
    s += "}\n";
    s += "}\n";
    return s;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        print_help();
        return 1;
    }
    const std::string cmd = argv[1];
    if (cmd == "-h" || cmd == "--help") {
        print_help();
        return 0;
    }
    if (cmd != "solve") {
        std::cerr << "Error: unknown command '" << cmd << "'. The only command is 'solve'.\n";
        return 1;
    }
    if (argc < 3) fail("solve requires a model file.");

    const std::string model_path = argv[2];
    std::string output_path;
    LdesOptions o;
    bool have_samples = false, have_events = false;

    for (int i = 3; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "-o") {
            output_path = need_value(argc, argv, i, a);
        } else if (a == "-s" || a == "--samples") {
            o.samples = static_cast<std::size_t>(std::atol(need_value(argc, argv, i, a).c_str()));
            have_samples = true;
        } else if (a == "-e" || a == "--maxevents") {
            o.events = static_cast<std::size_t>(std::atol(need_value(argc, argv, i, a).c_str()));
            have_events = true;
        } else if (a == "--maxtime") {
            o.timeout = std::atof(need_value(argc, argv, i, a).c_str());
        } else if (a == "--seed") {
            o.seed = std::atol(need_value(argc, argv, i, a).c_str());
        } else if (a == "--method") {
            o.method = need_value(argc, argv, i, a);
        } else if (a == "--cnvgon") {
            o.cnvgon = true;
        } else if (a == "--cnvgtol") {
            o.cnvgtol = std::atof(need_value(argc, argv, i, a).c_str());
        } else if (a == "--cnvgbatch") {
            o.cnvgbatch = std::atoi(need_value(argc, argv, i, a).c_str());
        } else if (a == "--cnvgchk") {
            o.cnvgchk = std::atoi(need_value(argc, argv, i, a).c_str());
        } else if (a == "--tranfilter") {
            o.tranfilter = need_value(argc, argv, i, a);
        } else if (a == "--mserbatch") {
            o.mserbatch = std::atoi(need_value(argc, argv, i, a).c_str());
        } else if (a == "--warmupfrac") {
            o.warmupfrac = std::atof(need_value(argc, argv, i, a).c_str());
        } else if (a == "--cimethod") {
            o.cimethod = need_value(argc, argv, i, a);
        } else if (a == "--obmoverlap") {
            o.obmoverlap = std::atof(need_value(argc, argv, i, a).c_str());
        } else if (a == "--ciminbatch") {
            o.ciminbatch = std::atoi(need_value(argc, argv, i, a).c_str());
        } else if (a == "--ciminobs") {
            o.ciminobs = std::atoi(need_value(argc, argv, i, a).c_str());
        } else if (a == "--spectrallowfreqfrac") {
            o.spectral_low_freq_frac = std::atof(need_value(argc, argv, i, a).c_str());
        } else if (a == "--slotted") {
            o.slotted = true;
        } else if (a == "--slotlength") {
            o.slot_length = std::atof(need_value(argc, argv, i, a).c_str());
            o.slotted = true;  // the reference lets --slotlength imply --slotted
        } else if (a == "--replications") {
            o.replications = std::atoi(need_value(argc, argv, i, a).c_str());
        } else if (a == "--numthreads") {
            o.numthreads = std::atoi(need_value(argc, argv, i, a).c_str());
        } else if (a == "--timespan") {
            const std::vector<double> t = parse_doubles(need_value(argc, argv, i, a));
            if (t.size() != 2) fail("--timespan takes T0,T1.");
            o.has_timespan = true;
            o.t0 = t[0];
            o.t1 = t[1];
        } else if (a == "--busyperiod") {
            o.busy_period_orders = std::atoi(need_value(argc, argv, i, a).c_str());
        } else if (a == "--busyperiod-subnet") {
            o.busy_period_subnets.push_back(parse_indices(need_value(argc, argv, i, a)));
        } else if (a == "--initsol") {
            o.init_sol = parse_doubles(need_value(argc, argv, i, a));
        } else if (a == "--export-histogram") {
            o.export_histogram = true;
        } else if (a == "--trajectory") {
            o.export_trajectory = true;
        } else if (a == "--respt-samples") {
            o.export_respt = true;
        } else if (a == "--rest") {
            refuse(a, "server mode belongs to line-cli, not to this binary");
        } else if (a == "-h" || a == "--help") {
            print_help();
            return 0;
        } else {
            fail("unknown option '" + a + "'.");
        }
    }
    if (have_events && !have_samples) o.samples = o.events;

    try {
        line::qn::Network<double> model = line::io::read_network_json<double>(model_path);
        const line::qn::NetworkStruct<double>& sn = model.get_struct();

        const std::clock_t t_start = std::clock();
        const LdesResult r = line::ldes::ldes_engine_solve(sn, o);
        const double runtime = static_cast<double>(std::clock() - t_start) / CLOCKS_PER_SEC;

        const std::string doc = result_json(r, o, runtime, o.samples);
        if (output_path.empty()) {
            std::cout << doc;
        } else {
            std::ofstream out(output_path.c_str());
            if (!out) fail("cannot write '" + output_path + "'.");
            out << doc;
        }
        std::cerr << "LDES analysis [method: " << o.method
                  << "; type: approximate, randomized; lang: cpp] completed in " << runtime
                  << "s. Iterations: " << o.samples << ".\n";
        return 0;
    } catch (const line::UnsupportedError& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 2;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}
