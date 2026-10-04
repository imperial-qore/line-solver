/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * `line-examples`: run one example of the gallery by name, or all of them.
 *
 *   line-examples --list              every example, grouped by its directory
 *   line-examples gallery_mm1         one, by the name its .py/.m file carries
 *   line-examples --group basic/openQN   one directory
 *   line-examples --all               all of them, in registration order
 *
 * The exit status is the number of examples that threw, so a sweep is usable
 * from a shell without parsing the output.
 */

#include <algorithm>
#include <cstdio>
#include <cstring>
#include <exception>
#include <string>
#include <vector>

#include "examples_common.h"
#include "parity_recorder.h"

namespace {

using line::examples::Example;
using line::examples::registry;

std::vector<Example> sorted_registry() {
    std::vector<Example> all = registry();
    std::sort(all.begin(), all.end(), [](const Example& a, const Example& b) {
        if (a.group != b.group) return a.group < b.group;
        return a.name < b.name;
    });
    return all;
}

void list_all() {
    const std::vector<Example> all = sorted_registry();
    std::string group;
    for (const Example& e : all) {
        if (e.group != group) {
            group = e.group;
            std::printf("\n%s\n", group.c_str());
        }
        std::printf("  %s\n", e.name.c_str());
    }
    std::printf("\n%zu examples\n", all.size());
}

/** Run one example, reporting a throw as the example's own failure line. */
int run_one(const Example& e) {
    std::printf("=== %s (%s)\n", e.name.c_str(), e.group.c_str());
    try {
        e.run();
    } catch (const std::exception& ex) {
        std::fflush(stdout);
        std::fprintf(stderr, "FAILED %s: %s\n", e.name.c_str(), ex.what());
        return 1;
    }
    return 0;
}

int usage() {
    std::printf(
        "usage: line-examples [--list] [--all] [--group <dir>] [<name> ...]\n"
        "  --list           print every example grouped by reference directory\n"
        "  --all            run every example\n"
        "  --group <dir>    run one reference directory, e.g. basic/openQN\n"
        "  --record <file>  also write what the run computed to <file> as JSON,\n"
        "                   keyed by the solver that computed it (see\n"
        "                   parity_recorder.h). The printed output is unchanged.\n"
        "  <name>           run the example of that name, e.g. gallery_mm1\n");
    return 2;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 2) return usage();

    const std::vector<Example> all = sorted_registry();
    std::vector<const Example*> todo;
    int failures = 0;

    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--list") {
            list_all();
            return 0;
        }
        if (arg == "-h" || arg == "--help") return usage();
        if (arg == "--record") {
            if (++i >= argc) return usage();
            line::examples::parity::enable(argv[i]);
            continue;
        }
        if (arg == "--all") {
            for (const Example& e : all) todo.push_back(&e);
            continue;
        }
        if (arg == "--group") {
            if (++i >= argc) return usage();
            const std::string g = argv[i];
            std::size_t n = 0;
            for (const Example& e : all)
                if (e.group == g) {
                    todo.push_back(&e);
                    ++n;
                }
            if (n == 0) {
                std::fprintf(stderr, "no such group: %s (try --list)\n", g.c_str());
                return 2;
            }
            continue;
        }
        const Example* found = nullptr;
        for (const Example& e : all)
            if (e.name == arg) found = &e;
        if (!found) {
            std::fprintf(stderr, "no such example: %s (try --list)\n", arg.c_str());
            return 2;
        }
        todo.push_back(found);
    }

    if (todo.empty()) return usage();
    for (const Example* e : todo) failures += run_one(*e);
    if (todo.size() > 1) std::printf("\n%zu run, %d failed\n", todo.size(), failures);
    // THE DUMP IS WRITTEN EVEN WHEN A TWIN THREW, because a run that produced
    // some tables and then died is exactly the case a consumer has to tell from
    // a run that produced none: the first has numbers to check and a defect
    // after them, the second has nothing. A dump that cannot be written is a
    // failure of the run, not a silent omission.
    if (!line::examples::parity::dump()) return failures + 1;
    return failures;
}
