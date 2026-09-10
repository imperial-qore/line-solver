/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * mp_pfqn oracle for the normalizing-constant family.
 *
 * mp_pfqn's bin/ca, bin/comom, bin/mom, bin/gld and bin/comomld print the
 * EXACT normalizing constant with -e, numerator on one line and denominator on
 * the next. Those rationals are the oracle: they come from an independent C
 * implementation over GMP, so agreement on the exact value is not something two
 * wrong implementations fall into.
 *
 * The comparison is on the RATIONALS, never on the printed doubles: mp_pfqn's
 * double output goes through a truncating mpf_get_d, so a last-ulp difference
 * there is an artifact of the printing and says nothing about the algorithms.
 *
 * Two halves:
 *  - a dynamic sweep over the model library, asserting that every ported exact
 *    algorithm returns the identical rational as pfqn_ca on the same model;
 *  - a static table of oracle values transcribed from the mp_pfqn binaries, so
 *    the cross-repository check needs no mp_pfqn build, only its models.
 *
 * THE MODELS ARE IN THE TREE, under cpp/tests/models/. They used to be read
 * from an absolute path in a personal checkout that no worker host could see,
 * which made all three cases pass having compared NOTHING -- 23 "file absent"
 * rows and a green doctest. Set LINE_MP_PFQN_MODELS to read an upstream
 * checkout instead.
 */
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>

#include "doctest.h"
#include "line/api/pfqn/pfqn_ca.h"
#include "line/api/pfqn/pfqn_comomrm.h"
#include "line/api/pfqn/pfqn_comomrm_ld.h"
#include "line/api/pfqn/pfqn_conv.h"
#include "line/api/pfqn/pfqn_gld.h"
#include "line/api/pfqn/pfqn_nc.h"
#include "line/api/pfqn/pfqn_ncld.h"
#include "line/io/qn_reader.h"

using line::Matrix;
using line::Rational;
using line::io::QnModel;
using line::io::read_qn;

namespace {

/** Directory holding the `.qn` models, with a trailing separator. */
std::string modelDir() {
    const char* env = std::getenv("LINE_MP_PFQN_MODELS");
    if (env != nullptr && *env != '\0') {
        std::string dir(env);
        if (dir.back() != '/') dir += '/';
        return dir;
    }
    return std::string(LINE_MP_REPO_ROOT) + "/cpp/tests/models/";
}

bool exists(const std::string& p) {
    std::ifstream f(p);
    return static_cast<bool>(f);
}

/** Build a rational from decimal numerator and denominator strings. */
Rational rat(const char* num, const char* den) {
    return Rational(line::BigInt(num), line::BigInt(den));
}

struct OracleRow {
    const char* model;
    const char* num;
    const char* den;
};

// Load-independent constants, from `bin/ca -e` (bin/comom and bin/mom agree).
const OracleRow kLiOracle[] = {
    {"01_single", "199902343750000000000", "1"},
    {"13_gld_small", "14560000", "1"},
    {"lcfs_2class", "18977", "1"},
    {"lcfs_3class", "67167920", "1"},
    {"test_singular7", "904931635000", "1"},
    {"02_bottleneck",
     "8334180091660472250391753671717605624531615977553389484678828681823214344627553248887580348"
     "37578651178821597701740220287547344923950731754302978515625000",
     "1"},
};

// Load-dependent constants, from `bin/gld -e` (bin/comomld agrees where it runs).
const OracleRow kLdOracle[] = {
    {"14_ld_multi", "3440000", "3"},
    {"15_repairman", "13541728", "27"},
};

}  // namespace

TEST_CASE("mp_pfqn oracle: load-independent constants match exactly") {
    const std::string dir = modelDir();
    int checked = 0, skipped = 0, missing = 0;
    for (const OracleRow& row : kLiOracle) {
        const std::string path = dir + row.model + ".qn";
        if (!exists(path)) {
            ++missing;
            continue;
        }
        const QnModel<Rational> m = read_qn<Rational>(path);
        REQUIRE(m.isClosedLoadIndependent());
        bool multiserver = false;
        for (int k : m.mi)
            if (k != 1) multiserver = true;
        if (multiserver) {
            // mp_pfqn's ca expands the multiplicity into replicas; that
            // convention is exercised in test_qn_models, not here.
            ++skipped;
            continue;
        }
        const Rational oracle = rat(row.num, row.den);
        INFO("model ", row.model);
        CHECK(line::pfqn::pfqn_ca(m.L, m.N, m.Z).G == oracle);
        // Every ported exact route must land on the same rational.
        CHECK(line::pfqn::pfqn_conv(m.L, m.N, m.Z).G == oracle);
        CHECK(line::pfqn::pfqn_nc(m.L, m.N, m.Z, line::pfqn::NcMethod::Ca).G == oracle);
        CHECK(line::pfqn::pfqn_nc(m.L, m.N, m.Z, line::pfqn::NcMethod::Exact).G == oracle);
        ++checked;
    }
    // A case that compared nothing must not report success. The models ship in
    // the tree, so an absent one is a broken checkout, not an absent machine.
    REQUIRE_MESSAGE(missing == 0, "load-independent oracle: " << missing
                                  << " model file(s) absent under " << dir);
    CHECK_MESSAGE(checked > 0, "load-independent oracle compared nothing (skipped "
                               << skipped << ")");
}

TEST_CASE("mp_pfqn oracle: load-dependent constants match exactly") {
    const std::string dir = modelDir();
    int checked = 0, missing = 0;
    for (const OracleRow& row : kLdOracle) {
        const std::string path = dir + row.model + ".qn";
        if (!exists(path)) {
            ++missing;
            continue;
        }
        const QnModel<Rational> m = read_qn<Rational>(path);
        REQUIRE(m.isLD);
        const Rational oracle = rat(row.num, row.den);
        INFO("model ", row.model);

        // mp_pfqn's gld takes no think time: a delay is a row whose rate
        // lattice is k. Build that model explicitly and check pfqn_gld, then
        // check that the pfqn_ncld dispatcher folds the delay the same way.
        const std::size_t D = m.Z.cols() == 0 ? 0 : 1;
        bool anyZ = false;
        for (int r = 0; r < m.R; ++r)
            if (m.Z(0, static_cast<std::size_t>(r)) != Rational(0)) anyZ = true;
        const std::size_t Drows = anyZ ? D : 0;
        Matrix<Rational> Lg(static_cast<std::size_t>(m.M) + Drows, static_cast<std::size_t>(m.R));
        Matrix<Rational> mug(Lg.rows(), static_cast<std::size_t>(m.Nt));
        for (std::size_t i = 0; i < static_cast<std::size_t>(m.M); ++i) {
            for (std::size_t r = 0; r < static_cast<std::size_t>(m.R); ++r) Lg(i, r) = m.L(i, r);
            for (std::size_t k = 0; k < static_cast<std::size_t>(m.Nt); ++k) mug(i, k) = m.mu(i, k);
        }
        for (std::size_t d = 0; d < Drows; ++d) {
            for (std::size_t r = 0; r < static_cast<std::size_t>(m.R); ++r)
                Lg(static_cast<std::size_t>(m.M) + d, r) = m.Z(d, r);
            for (std::size_t k = 0; k < static_cast<std::size_t>(m.Nt); ++k)
                mug(static_cast<std::size_t>(m.M) + d, k) = Rational(static_cast<long>(k) + 1);
        }
        CHECK(line::pfqn::pfqn_gld(Lg, m.N, mug).G == oracle);
        CHECK(line::pfqn::pfqn_ncld(m.L, m.N, m.Z, m.mu).G == oracle);
        ++checked;
    }
    REQUIRE_MESSAGE(missing == 0, "load-dependent oracle: " << missing
                                  << " model file(s) absent under " << dir);
    CHECK_MESSAGE(checked > 0, "load-dependent oracle compared nothing");
}

TEST_CASE("mp_pfqn model library: every ported exact route agrees with pfqn_ca") {
    // Dynamic sweep. Broader than the static table above: it asserts the ported
    // routes agree with each other on every model, without an oracle constant.
    const char* names[] = {"01_single",  "02_bottleneck",  "03_think",       "05_sparse",
                           "07_mixed",   "08_multiclass",  "09_asymmetric",  "10_diverse",
                           "12_expanded", "13_gld_small",  "lcfs_1class",    "lcfs_2class",
                           "lcfs_3class", "test_singular7", "test_singular8"};
    const std::string dir = modelDir();
    int checked = 0, skipped = 0, missing = 0;
    for (const char* name : names) {
        const std::string path = dir + name + ".qn";
        if (!exists(path)) {
            ++missing;
            continue;
        }
        const QnModel<Rational> m = read_qn<Rational>(path);
        if (!m.isClosedLoadIndependent()) {
            ++skipped;
            continue;
        }
        bool multiserver = false;
        for (int k : m.mi)
            if (k != 1) multiserver = true;
        if (multiserver) {
            ++skipped;
            continue;
        }
        const Rational gca = line::pfqn::pfqn_ca(m.L, m.N, m.Z).G;
        INFO("model ", name);
        CHECK(line::pfqn::pfqn_conv(m.L, m.N, m.Z).G == gca);
        CHECK(line::pfqn::pfqn_nc(m.L, m.N, m.Z, line::pfqn::NcMethod::Ca).G == gca);
        // The repairman routines only accept a single queueing station.
        if (m.M == 1) {
            CHECK(line::pfqn::pfqn_comomrm(m.L, m.N, m.Z).G == gca);
            CHECK(line::pfqn::pfqn_nc(m.L, m.N, m.Z, line::pfqn::NcMethod::Comom).G == gca);
        }
        // Unit rates make the load-dependent route the load-independent one.
        Matrix<Rational> mu1(static_cast<std::size_t>(m.M), static_cast<std::size_t>(m.Nt));
        for (std::size_t i = 0; i < mu1.rows(); ++i)
            for (std::size_t k = 0; k < mu1.cols(); ++k) mu1(i, k) = Rational(1);
        CHECK(line::pfqn::pfqn_ncld(m.L, m.N, m.Z, mu1).G == gca);
        ++checked;
    }
    REQUIRE_MESSAGE(missing == 0, "model-library sweep: " << missing
                                  << " model file(s) absent under " << dir);
    CHECK_MESSAGE(checked > 0, "model-library sweep compared nothing (skipped "
                               << skipped << ")");
}
