/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * The citation registry (`line/io/line_citations.h`).
 *
 * Attribution in LINE is pull-based: nothing is printed during a solve, and a
 * method that reaches a user without an entry is a silent loss of attribution.
 * C++ carried NO registry at all until this one, so the `.citations()` contract
 * the other three codebases honour had no counterpart here.
 *
 * The table is generated from `matlab/src/io/line_citations.m`, the reference
 * implementation, so these assertions are about the LOOKUP CONTRACT rather than
 * about particular wording: normalisation, the family-qualified fallback,
 * de-duplication by bibliography key, and the `default/<actual>` compound a
 * dispatching solver reports.
 *
 * The MATLAB twin is line_citations.m, the JAR twin LineCitations.java and the
 * python twin api/io/citations.py.
 */
#include "doctest.h"
#include "line/io/line_citations.h"

using line::io::Citation;
using line::io::citation_registry;
using line::io::line_citations;
using line::io::line_citations_for_method;

TEST_CASE("the registry carries the reference implementation's entries") {
    // Generated from line_citations.m, which is the superset: python and the JAR
    // each carry 498 of these, lacking only method names for MATLAB-only features
    // (the dae.petri family, pfqn_divdiff_ld, asy).
    CHECK(citation_registry().size() == 505);
}

TEST_CASE("a known method name resolves to its bibliography key") {
    CHECK(line_citations({"bs"}).at(0).key == "Sch79");
    CHECK(line_citations({"ca"}).at(0).key == "ReiK75");
    CHECK(line_citations({"lc"}).at(0).key == "BirK92");
    CHECK(line_citations({"chains"}).at(0).key == "ReiL80");
}

TEST_CASE("an unknown method name is ignored rather than an error") {
    // A caller passes whatever it knows about a run, so an unrecognised
    // method name must not be fatal.
    CHECK(line_citations({"bogus"}).empty());
    CHECK(line_citations({}).empty());
    CHECK(line_citations({"", "   "}).empty());
}

TEST_CASE("a method name is normalised before lookup") {
    CHECK(line_citations({"  BS  "}).at(0).key == "Sch79");
}

TEST_CASE("a family-qualified method name falls back to the bare one") {
    // `bs` has an entry; `zz.bs` has none, so the lookup falls back to the bare
    // method name after the first dot. This is what lets a caller pass a
    // solver-qualified name it read off a result without knowing whether the
    // registry happens to carry the qualified form.
    REQUIRE(citation_registry().count("bs") == 1);
    REQUIRE(citation_registry().count("zz.bs") == 0);
    CHECK(line_citations({"zz.bs"}).at(0).key == "Sch79");
}

TEST_CASE("one paper reached through several method names is reported once") {
    const std::vector<Citation> v = line_citations({"ca", "nc.ca", "ca"});
    CHECK(v.size() == 1);
    CHECK(v.at(0).key == "ReiK75");
}

TEST_CASE("the reported compound method name resolves both halves") {
    // A dispatching solver reports `default/<actual>`, and a transformed solve
    // reports `<method>/<method name>`; both must reach the underlying paper.
    const std::vector<Citation> v = line_citations_for_method("default/lc");
    REQUIRE(v.size() == 1);
    CHECK(v.at(0).key == "BirK92");
    CHECK(line_citations_for_method("lc").at(0).key == "BirK92");
}

TEST_CASE("every entry is fully populated") {
    // An entry with an empty key or reference is worse than none: it reports
    // attribution that cannot be looked up.
    for (const auto& kv : citation_registry()) {
        CHECK_FALSE(kv.first.empty());
        CHECK_FALSE(kv.second.key.empty());
        CHECK_FALSE(kv.second.ref.empty());
        CHECK_FALSE(kv.second.covers.empty());
    }
}
