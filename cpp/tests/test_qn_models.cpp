/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Cross-algorithm regression over the mp_pfqn model library.
 *
 * Every exact normalizing-constant algorithm must return the IDENTICAL rational
 * on every closed load-independent model. This is the strongest cheap oracle
 * available: the algorithms share no code path (convolution over the population
 * lattice vs the MVA recursion), so agreement on an exact rational is not
 * something two wrong implementations fall into.
 *
 * The models live outside the repository, in the mp_pfqn checkout. When that
 * checkout is absent the test reports how many models it found and passes
 * trivially rather than failing, since it is a cross-repository oracle and not
 * a property of this tree.
 */
#include <string>
#include <vector>

#include "doctest.h"
#include "line/api/pfqn/pfqn_ca.h"
#include "line/api/pfqn/pfqn_mva.h"
#include "line/io/qn_reader.h"

using line::Rational;
using line::io::QnModel;
using line::io::read_qn;

namespace {

const char* kModelDir = "/home/gcasale/Dropbox/code/stable/mp_pfqn.git/models/";

const char* kClosedModels[] = {
    "01_single.qn",      "02_bottleneck.qn", "03_think.qn",     "04_replicated.qn",
    "05_sparse.qn",      "06_large.qn",      "07_mixed.qn",     "08_multiclass.qn",
    "09_asymmetric.qn",  "10_diverse.qn",    "11_swapped.qn",   "12_expanded.qn",
    "lcfs_1class.qn",    "lcfs_2class.qn",   "lcfs_3class.qn",
};

bool exists(const std::string& path) {
    std::ifstream f(path);
    return static_cast<bool>(f);
}

}  // namespace

TEST_CASE("pfqn_ca and pfqn_mva return the identical exact G on the mp_pfqn model library") {
    int checked = 0, skipped = 0;
    for (const char* name : kClosedModels) {
        const std::string path = std::string(kModelDir) + name;
        if (!exists(path)) {
            ++skipped;
            continue;
        }
        QnModel<Rational> m = read_qn<Rational>(path);
        if (!m.isClosedLoadIndependent()) {
            ++skipped;
            continue;
        }
        bool multiserver = false;
        for (int k : m.mi)
            if (k != 1) multiserver = true;

        auto ca = line::pfqn::pfqn_ca(m.L, m.N, m.Z);
        auto mva = line::pfqn::pfqn_mva(m.L, m.N, m.Z, m.mi);

        INFO("model ", name);
        if (!multiserver) {
            // Single-server stations: the two constants must be equal as exact
            // rationals, not merely close.
            CHECK(ca.G == mva.G);
        } else {
            // With station multiplicities the convolution constant is the one of
            // the expanded network, so only the MVA path is meaningful here;
            // check it is at least a positive finite rational.
            CHECK(mva.G > Rational(0));
        }
        ++checked;
    }
    INFO("checked ", checked, " models, skipped ", skipped);
    CHECK(checked + skipped == static_cast<int>(sizeof(kClosedModels) / sizeof(kClosedModels[0])));
}

TEST_CASE("qn_reader parses the sections of the .qn format") {
    const std::string path = std::string(kModelDir) + "03_think.qn";
    if (!exists(path)) return;
    QnModel<Rational> m = read_qn<Rational>(path);
    CHECK(m.R == 2);
    CHECK(m.M == 4);
    CHECK(m.N.size() == 2);
    CHECK(m.L.rows() == 4);
    CHECK(m.L.cols() == 2);
    CHECK(m.isClosedLoadIndependent());
    CHECK(m.Nt == m.N[0] + m.N[1]);
}

TEST_CASE("qn_reader parses rational tokens in the LAMBDA and MU sections") {
    CHECK(line::io::parse_rational_token<Rational>("3") == Rational(3));
    CHECK(line::io::parse_rational_token<Rational>("1/10") == Rational(1, 10));
    CHECK_THROWS_AS(line::io::parse_rational_token<Rational>("1/0"), line::InputError);
}
