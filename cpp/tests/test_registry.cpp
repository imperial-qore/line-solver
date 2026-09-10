/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * Well-formedness of the coverage registry.
 *
 * The registry is a single brace-initialized list that several concurrent
 * contributors append to, and its failure mode is not obvious: an entry whose
 * description line is orphaned by an insertion in the middle of it still
 * PARSES, because the following entry's opening brace is accepted as the
 * description's place, and the whole 450-element initializer then fails to
 * convert with one diagnostic that dumps every element. That reads as "the
 * registry is broken" without saying where. It happened twice in one day, to
 * "pfqn_mva" and to "cache_miss".
 *
 * These cases turn that into a located failure. They also enforce the
 * invariants the CLI relies on: find_api() returns the first match, so a
 * duplicated name silently shadows one of the two implementations, and
 * --list-api prints the domain and reference verbatim, so neither may be empty.
 */

#include <algorithm>
#include <cctype>
#include <filesystem>
#include <fstream>

#include "doctest.h"
#include "line/reg/registry.h"

#include <map>
#include <set>
#include <string>

using namespace line;

TEST_CASE("registry: every entry is fully populated") {
    for (const ApiEntry& e : api_registry()) {
        CHECK_MESSAGE(!e.name.empty(), "entry with an empty name");
        CHECK_MESSAGE(!e.domain.empty(), "empty domain for " << e.name);
        CHECK_MESSAGE(!e.reference.empty(), "empty reference for " << e.name);
        CHECK_MESSAGE(!e.arith.empty(), "no arithmetic mode for " << e.name);
    }
}

TEST_CASE("registry: names are unique") {
    std::set<std::string> seen;
    for (const ApiEntry& e : api_registry()) {
        CHECK_MESSAGE(seen.insert(e.name).second, "duplicate registry name: " << e.name);
    }
}

/**
 * A reference must name a source path, not another entry's leftovers. The
 * orphaning bug produces a reference that is a function name rather than a
 * path, which this catches even when the braces happen to balance.
 */
TEST_CASE("registry: references are source paths") {
    for (const ApiEntry& e : api_registry()) {
        CHECK_MESSAGE(e.reference.find('/') != std::string::npos,
                      "reference for " << e.name << " is not a path: " << e.reference);
    }
}

/** The name prefix and the declared domain must agree, so --list-api groups correctly. */
TEST_CASE("registry: domain is consistent within a name prefix") {
    std::map<std::string, std::string> domainOfPrefix;
    for (const ApiEntry& e : api_registry()) {
        std::string::size_type us = e.name.find('_');
        if (us == std::string::npos) continue;
        std::string prefix = e.name.substr(0, us);
        std::map<std::string, std::string>::iterator it = domainOfPrefix.find(prefix);
        if (it == domainOfPrefix.end()) {
            domainOfPrefix[prefix] = e.domain;
        } else {
            CHECK_MESSAGE(it->second == e.domain,
                          e.name << " is in domain " << e.domain << " but prefix " << prefix
                                 << " was registered under " << it->second);
        }
    }
}

TEST_CASE("registry: lookup and arithmetic queries agree with the table") {
    const ApiEntry* ca = find_api("pfqn_ca");
    REQUIRE(ca != nullptr);
    CHECK(ca->domain == "pfqn");
    CHECK(api_supports(*ca, Arith::Exact));

    CHECK(find_api("no_such_function") == nullptr);
}

/**
 * Every `.m` path a registry entry claims to be ported from must exist.
 *
 * The reference is not decoration: it is the only cross-reference between a
 * ported function and the MATLAB it must agree with, and it is what a reader
 * follows to check the port. A fabricated path therefore fails silently in the
 * worst way, by looking authoritative in --list-api. Two were found by hand on
 * 2026-07-21: ljd_linearize claimed matlab/src/api/ljd/ (it lives under pfqn/),
 * and pfqn_perm claimed a matlab/src/api/pfqn/pfqn_perm.m that has never
 * existed, the MATLAB function being `perm` in matlab/util.
 *
 * Only `.m` tokens are checked. A reference may list several sources separated
 * by ';' (e.g. a MATLAB file and an mp_pfqn module), and a few legitimately
 * name a directory rather than a file.
 */
TEST_CASE("registry: every referenced MATLAB source exists") {
    const std::string root = LINE_MP_REPO_ROOT;
    for (const ApiEntry& e : api_registry()) {
        std::string::size_type start = 0;
        while (start <= e.reference.size()) {
            std::string::size_type semi = e.reference.find(';', start);
            std::string tok = e.reference.substr(
                start, semi == std::string::npos ? std::string::npos : semi - start);
            // trim
            std::string::size_type b = tok.find_first_not_of(" \t");
            std::string::size_type f = tok.find_last_not_of(" \t");
            if (b != std::string::npos) tok = tok.substr(b, f - b + 1);
            if (tok.size() > 2 && tok.compare(tok.size() - 2, 2, ".m") == 0) {
                std::ifstream in((root + "/" + tok).c_str());
                CHECK_MESSAGE(in.good(), e.name << " references a file that does not exist: "
                                                << tok);
            }
            if (semi == std::string::npos) break;
            start = semi + 1;
        }
    }
}

/**
 * A direct same-name header is the strongest mechanical evidence that an
 * active MATLAB API entry point has a C++ implementation. Keep that surface
 * represented in the coverage registry; MDD is the sole same-name class, not
 * a function, and therefore belongs to an explicit exclusion rather than a
 * misleading registry row.
 */
TEST_CASE("registry: every direct active MATLAB API port is represented") {
    namespace fs = std::filesystem;
    const fs::path root = LINE_MP_REPO_ROOT;
    std::set<std::string> matlab_names;
    std::set<std::string> cpp_names;
    const auto lower_stem = [](const fs::path& path) {
        std::string name = path.stem().string();
        std::transform(name.begin(), name.end(), name.begin(), [](unsigned char c) {
            return static_cast<char>(std::tolower(c));
        });
        return name;
    };

    for (const fs::directory_entry& entry :
         fs::recursive_directory_iterator(root / "matlab/src/api")) {
        if (entry.is_regular_file() && entry.path().extension() == ".m")
            matlab_names.insert(lower_stem(entry.path()));
    }
    for (const fs::directory_entry& entry :
         fs::recursive_directory_iterator(root / "cpp/include/line/api")) {
        if (entry.is_regular_file() && entry.path().extension() == ".h")
            cpp_names.insert(lower_stem(entry.path()));
    }

    std::size_t direct_functions = 0;
    bool saw_mdd_class = false;
    for (const std::string& name : cpp_names) {
        if (matlab_names.count(name) == 0) continue;
        if (name == "mdd") {
            saw_mdd_class = true;
            continue;
        }
        ++direct_functions;
        CHECK_MESSAGE(find_api(name) != nullptr,
                      "same-name MATLAB/C++ API port is absent from registry: " << name);
    }
    CHECK(saw_mdd_class);
    CHECK(direct_functions >= 82);
}
