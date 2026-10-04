#!/usr/bin/env bash
# Copyright (c) 2012-2026, QORE Lab, Imperial College London
# All rights reserved.
#
# Build of the C++ port, installing line-cli and ldes into ../common/.
#
# Builds into a directory OUTSIDE the source tree so a stale glob can never
# shadow a new test file: cpp/CMakeLists collects tests/*.cpp at CONFIGURE time,
# so reusing a directory configured before a file existed compiles neither the
# file nor a warning about it. Pass -f to force a fresh configure when tests/
# has gained or lost a file.
#
# -O0 IS THE RULE IN A DEVELOPMENT TREE, and the ONLY exception in the whole
# repository is release/pre-release.sh, which pins the optimized build for the
# archive it stages. cpp/run-tests.sh CAN ask for the -O2 mode below (`-2`) and
# briefly did by default on 2026-09-07, but that build does not link, so the
# suite is back on -O0 too. src/cli/line_cli.cpp is one translation unit that includes
# every solver header, and at -O3 the optimizer, not the parser, is the cost:
# measured 88s at -O0 against ~900s at -O3 -march=native on the same file, a
# ~10x turnaround, with cc1plus peaking above 12 GB resident -- and it holds the
# shared cmake flock for all of it, stalling every other session.
#
# THE NUMBERS ARE NOT LESS TRUSTWORTHY FOR IT, which is the point. -O0 drops
# -march=native, so the binary uses base x86-64 instead of this box's widest
# vector unit, and cpp/CMakeLists already forces -ffp-contract=off in every
# configuration -- so the arithmetic is if anything CLOSER to the reference's
# than the optimized build. It is a slower binary, not a different one, so it is
# what common/line-cli and common/ldes are installed from. Time nothing here.
#
# -O2 IS THE MIDDLE MODE AND IT EXISTS FOR THE SUITE, which is RUN far more than
# it is compiled: at -O0 the run is the whole cost of the C++ test phase
# (measured 2026-09-07 on picard04: 4% of 3490 cases in 12 minutes, eta 2h49m,
# against a ~35 minute build), so paying the optimizer once is the cheaper half
# of the trade. It deliberately does NOT add -march=native: a suite binary is
# built on one roster host and read by others, and the portability rule that
# relocates a too-new build off the driver applies to the instruction set as much
# as to glibc. Its own build directory keeps it from invalidating the other two.
# IT LINKS, CONFIRMED 2026-09-11, and the reason nothing takes it by default is
# now the COMPILE rather than the link. The recorded blocker -- an undefined
# reference to the anonymous-namespace roscor_ of third_party/rodas.hpp -- was
# fixed on 2026-09-10 by namespace-scope forward declarations, and the U_fp /
# S_fp diagnosis that went with it was wrong anyway: both are `int (*)(...)`.
# What -O2 costs on this tree is 17.9 GB resident in cc1plus for line_cli.cpp
# and ~27 minutes of CPU, so a 16 GB host cannot build it at all. See the note
# at the top of cpp/run-tests.sh and _kb/08-build-and-test.md.
#
# Every compiler invocation goes through task-spooler (tsp) WHEN ONE IS INSTALLED.
# tsp ALONE DOES NOT SERIALIZE: it runs `tsp -S` jobs at once (4 here), so
# concurrent sessions were driving cmake into the same build directory
# simultaneously and racing on the object files and the linked binary. Each cmake
# invocation therefore also takes an flock on a shared lock file, which is what
# actually serializes them; other sessions in this repository use the same lock path.
#
# THE FLOCK IS THE MUTUAL EXCLUSION AND tsp IS ONLY THE QUEUE, so a box without
# task-spooler BUILDS DIRECTLY rather than refusing to build: the same lock is
# taken, in the same order, and the only thing lost is the machine-wide bound on
# how many compiles run at once. flock is the hard requirement here; tsp is not.

set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
build="${LINE_MP_BUILD_DIR:-/tmp/line-mp-build}"
common="$(cd "${here}/.." && pwd)/common"
jobs="${LINE_MP_JOBS:-4}"

force=0
tests=OFF
# `-t` means BUILD the doctest binary; whether to RUN it is a separate question,
# which `-n` answers. The suite driver splits the C++ work into a BUILD phase and
# a TEST phase so that the phases which merely consume common/line-cli do not
# wait out a multi-hour doctest run, and the build phase needs the test binary
# compiled -- that compile is most of what proves the tree builds -- without
# spending that run here.
run_tests=1

# UNOPTIMIZED IN A DEVELOPMENT TREE, agent session or human shell alike. Nearly
# every build here exists to read the NUMBERS -- a parity row, a refusal
# message, a table -- and not to time the binary, so paying ~900s of optimizer
# for a result that is looked at once is the wrong default for both. `-O` exists
# only so that pre-release.sh's pinned tree has something to pin; do not pass it
# by hand in a development tree.
#
# A DISTRIBUTED TREE PINS THIS TO 0. release/pre-release.sh rewrites the line
# below when it stages cpp/ into an archive, and verifies the rewrite took, so
# a released tree builds -O3 no matter who runs the build or what their
# environment says. A user is not iterating; they are building the binary that
# ships, and a silently unoptimized release binary is the one failure this
# default could otherwise cause. The line is its own statement with no leading
# whitespace so the release-time substitution is exact.
quick_default=1
quick="$quick_default"
# The third mode is never a default, in a released tree either, and since the
# -O2 link failure of 2026-09-07 nothing asks for it by name either.
o2=0
# THE FOURTH MODE IS WHAT A PUBLISHED BINARY IS BUILT WITH: -O3 exactly as -O,
# but at the BASE x86-64 instruction set rather than this box's widest. -march=native
# bakes in whatever the build host happens to carry (AVX-512, BMI2, ADX), and the
# first client with an older CPU dies of SIGILL before main -- a failure the
# glibc/GLIBCXX floor upload-line-cli.sh reports cannot see, because it is the
# instruction set and not a symbol version, and one the roster itself would hit:
# picard04 and picard05 have no BMI2/ADX while the newer hosts do. Every caller
# that hands a binary to someone else asks for it by name: upload-line-cli.sh,
# upload-ldes.sh, and the binaries release/pre-release.sh stages into an archive.
portable=0
# STATIC DEPENDENCIES: OFF unless asked for, and asked for only when the binary
# is going to another machine. See LINE_MP_STATIC_DEPS in CMakeLists.txt -- it
# carries HiGHS, LAPACK/BLAS and the C++/Fortran runtimes INTO the binary, which
# is what lets a downloader who has never heard of HiGHS run it.
static=0

# INSTALLING INTO common/ IS THE DEFAULT AND STAYS THE DEFAULT; -I is for the one
# caller that must not, the C++ TEST phase. See the install block below.
install_common=1
# An explicit flag always beats the environment, in both directions.
while [ $# -gt 0 ]; do
    case "$1" in
        -f|--fresh)  force=1 ;;
        -t|--tests)  tests=ON ;;
        -n|--no-run) run_tests=0 ;;
        -I|--no-install) install_common=0 ;;
        -q|--quick)  quick=1 ;;
        -O|--optimized) quick=0; o2=0 ;;
        -2|--o2)     quick=0; o2=1 ;;
        -P|--portable) quick=0; o2=0; portable=1 ;;
        -S|--static) static=1 ;;
        -h|--help)
            echo "usage: $0 [-f|--fresh] [-t|--tests] [-n|--no-run] [-I|--no-install] [-q|--quick] [-2|--o2] [-O|--optimized] [-P|--portable] [-S|--static]"
            echo "  -f  reconfigure from scratch (needed when tests/ gains a file)"
            echo "  -t  build and run the doctest suite as well"
            echo "  -n  with -t, BUILD the doctest binary but do not run it. This is"
            echo "      what the suite driver's C++ BUILD phase uses: it compiles"
            echo "      everything, installs into common/, and leaves the run to the"
            echo "      separate C++ TEST phase, which can then go to another machine."
            echo "  -I  build, but do NOT install into common/. For the C++ TEST"
            echo "      phase only: the BUILD phase installed those binaries"
            echo "      already, and the other phases are executing them by then."
            echo "  -q  unoptimized build (-O0, no -march=native): ~10x faster to"
            echo "      compile, slower to run, same numbers. This is the rule in a"
            echo "      development tree, and it installs over common/line-cli."
            echo "  -2  middle build (-O2, no -march=native): what cpp/run-tests.sh"
            echo "      builds the doctest suite with. Dearer to compile than -q,"
            echo "      much faster to RUN, and portable across the roster because"
            echo "      it takes no -march. Its own build dir, so it invalidates"
            echo "      neither of the other two."
            echo "  -O  optimized build (-O3 -march=native), ~900s against ~88s and"
            echo "      the shared cmake lock held for all of it. Not for iterating;"
            echo "      by hand use make-O3.sh, which adds a per-tree build dir."
            echo "  -P  portable optimized build (-O3, -march=x86-64 -mtune=generic):"
            echo "      what every PUBLISHED binary is built with. Same -O3 cost, but"
            echo "      it runs on any x86-64 CPU instead of only those as new as this"
            echo "      build host. Its own build dir."
            echo "  -S  link HiGHS, LAPACK/BLAS and the C++/Fortran runtimes STATICALLY,"
            echo "      leaving only libc and libm. The other half of publishing: a"
            echo "      downloader has no libhighs.so.1 and cannot get one from a"
            echo "      package. Combines with any mode and adds -static to the build"
            echo "      dir; the publisher uses -P -S."
            echo "env: LINE_MP_BUILD_DIR (default /tmp/line-mp-build), LINE_MP_JOBS (default 4)"
            # The help must describe THIS tree, not the one it was written in:
            # a released tree is pinned to the optimized build, and telling its
            # user that -q is the default would be false.
            if [ "$quick_default" = 1 ]; then
                echo "     default: -q (unoptimized), and it installs into common/"
            else
                echo "     default: optimized (-O3); this tree is pinned by pre-release.sh"
            fi
            exit 0 ;;
        *) echo "$0: unknown argument '$1'" >&2; exit 2 ;;
    esac
    shift
done

# `-t` IN A DISTRIBUTED TREE HAS NOTHING TO BUILD, and must say so itself.
# release/pre-release.sh strips cpp/tests (and cpp/bench) from every archive, as
# it already stripped the MATLAB, Java and Python suites, so cpp/CMakeLists
# skips the target when the glob comes back empty. Asking cmake to build a
# target that was never defined reports "No rule to make target", which reads as
# a broken build system rather than as an absent suite -- and the suite is
# absent BY DESIGN here: it needs cpp/run-tests.sh and the shared goldens/ tree,
# neither of which an archive carries.
#
# THE SUITE MOVED TO line-test.git (2026-09-12), so BOTH directories count: the
# in-tree one is empty in a development checkout now, and testing it alone
# refused a suite of 365 translation units that was sitting in the sibling.
# cpp/CMakeLists.txt globs the same pair, and cpp/run-tests.sh checksums it.
if [ -z "${LINE_TEST_DIR:-}" ]; then
    _mk_common="$(git -C "${here}" rev-parse --path-format=absolute --git-common-dir 2>/dev/null)"
    if [ -n "$_mk_common" ]; then
        LINE_TEST_DIR="$(dirname "$(dirname "$_mk_common")")/line-test.git"
    fi
    [ -d "${LINE_TEST_DIR:-}" ] || LINE_TEST_DIR="${here}/../../line-test.git"
fi
export LINE_TEST_DIR
if [ "$tests" = ON ] \
   && ! ls "${here}"/tests/*.cpp >/dev/null 2>&1 \
   && ! ls "${LINE_TEST_DIR}"/cpp/tests/*.cpp >/dev/null 2>&1; then
    echo "$0: no cpp/tests in this tree or in ${LINE_TEST_DIR}, so there is no" >&2
    echo "  doctest suite to build. The suite lives in line-test.git since" >&2
    echo "  2026-09-12; set LINE_TEST_DIR if that checkout is elsewhere." >&2
    echo "  Release archives ship cpp/ as source WITHOUT its tests, and for one" >&2
    echo "  of those this message is correct: drop -t to build line-cli, ldes" >&2
    echo "  and line-examples." >&2
    exit 2
fi

# A SEPARATE DIRECTORY PER MODE, not a flag flipped in one. Changing
# CMAKE_CXX_FLAGS_RELEASE invalidates every object in the directory, so sharing
# one would make each alternation between the three modes a full rebuild --
# exactly the cost -q exists to avoid.
if [ "$quick" = 1 ]; then
    build="${build}-quick"
    cxxflags="-O0 -DNDEBUG"
elif [ "$o2" = 1 ]; then
    build="${build}-o2"
    cxxflags="-O2 -DNDEBUG"
elif [ "$portable" = 1 ]; then
    build="${build}-portable"
    # -march=x86-64 is the BASELINE (SSE2 only), not a compromise level: a
    # published binary has no way to know what its reader runs, and -mtune keeps
    # the scheduling sensible without constraining the ISA. -ffp-contract=off is
    # forced for every configuration in cpp/CMakeLists.txt, so dropping the wide
    # vector units moves the arithmetic towards the MATLAB reference, not away.
    cxxflags="-O3 -DNDEBUG -march=x86-64 -mtune=generic"
else
    cxxflags="-O3 -DNDEBUG -march=native"
fi

# The static link is orthogonal to the optimization mode, so it appends rather
# than replaces -- `-P -S` builds in <dir>-portable-static. It still needs a
# directory of its own: the archives are chosen at CONFIGURE time, and sharing a
# directory with the dynamic build would silently serve whichever configure ran
# last (the cmake options below are only passed on a fresh configure).
cmake_static=()
if [ "$static" = 1 ]; then
    build="${build}-static"
    cmake_static=(-DLINE_MP_STATIC_DEPS=ON)
fi

# ccache turns a rebuild of UNCHANGED translation units into a copy, which is
# what a `-f` reconfigure or a branch switch otherwise pays full price for. It
# is a launcher, so it applies to whatever compiler CMake picked.
launcher=()
if command -v ccache >/dev/null 2>&1; then
    launcher=(-DCMAKE_CXX_COMPILER_LAUNCHER=ccache)
fi

command -v flock >/dev/null || { echo "$0: flock not found" >&2; exit 1; }

# Resolve a WORKING cmake rather than trusting the first name on PATH. A pip
# `cmake` console script outlives the interpreter it was installed for: the
# wrapper sits in ~/.local/bin, shadows /usr/bin/cmake, and imports a module
# that the current python3 no longer has -- so configure died on a Python
# traceback, left no cache and no object file, and was reported as a C++ test
# failure on a box whose own cmake was fine the whole time. Probe candidates in
# PATH order and take the first that actually answers --version.
cmake_candidates="$(type -aP cmake 2>/dev/null || true)"
cmake_bin=""
for cand in ${LINE_MP_CMAKE:+"$LINE_MP_CMAKE"} $cmake_candidates /usr/bin/cmake /usr/local/bin/cmake; do
    [ -x "$cand" ] || continue
    if "$cand" --version >/dev/null 2>&1; then
        cmake_bin="$cand"
        break
    fi
    echo "$0: ignoring broken cmake at ${cand}" >&2
done
if [ -z "$cmake_bin" ]; then
    echo "$0: no working cmake on PATH, /usr/bin or /usr/local/bin." >&2
    echo "  Set LINE_MP_CMAKE to one that runs." >&2
    exit 1
fi

# Resolve the queue once. Already being inside a tsp job means the slot is
# already held, so queueing again would deadlock against the same `-S` bound.
tsp_bin=""
if [ -n "${TS_SOCKET:-}" ]; then
    echo "$0: already inside a tsp job; running cmake directly"
elif command -v tsp >/dev/null 2>&1; then
    tsp_candidate="$(command -v tsp)"
    if "$tsp_candidate" -l >/dev/null 2>&1; then
        tsp_bin="$tsp_candidate"
    else
        echo "$0: tsp is installed but its server is unavailable; running cmake directly (unqueued, still flocked)" >&2
    fi
elif [ -x "${HOME:-}/.local/bin/tsp" ]; then
    tsp_candidate="${HOME}/.local/bin/tsp"
    if "$tsp_candidate" -l >/dev/null 2>&1; then
        tsp_bin="$tsp_candidate"
    else
        echo "$0: tsp is installed but its server is unavailable; running cmake directly (unqueued, still flocked)" >&2
    fi
else
    echo "$0: tsp (task-spooler) not found; running cmake directly (unqueued, still flocked)" >&2
fi

# Every cmake invocation goes through queue_cmake, and the caller passes the
# flock either way, so dropping tsp costs the queue slot and nothing else.
queue_cmake() {
    if [ -n "$tsp_bin" ]; then
        "$tsp_bin" -fn -L cmake "$@"
    else
        "$@"
    fi
}

# The test binary is not a compile, so it queues unlabelled and takes no lock.
queue_run() {
    if [ -n "$tsp_bin" ]; then
        "$tsp_bin" -fn "$@"
    else
        "$@"
    fi
}

lockdir=/tmp/line-tsp-locks
mkdir -p "$lockdir"
cmakelock="${lockdir}/cmake.lock"

if [ "$force" = 1 ]; then
    rm -rf "$build"
fi

# A failed configure can leave CMakeCache.txt behind without generating a build
# system. TargetDirectories.txt is written only once generation completes (for
# both Makefiles and Ninja), so require it as well or every retry degrades into
# the unrelated "No rule to make target 'Makefile'" error.
if [ ! -f "${build}/CMakeCache.txt" ] || [ ! -f "${build}/CMakeFiles/TargetDirectories.txt" ]; then
    queue_cmake flock "$cmakelock" "$cmake_bin" -S "$here" -B "$build" \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_CXX_FLAGS_RELEASE="$cxxflags" \
        "${launcher[@]}" \
        "${cmake_static[@]}" \
        -DLINE_MP_BUILD_TESTS="$tests" \
        -DLINE_MP_BUILD_BENCH=OFF
fi

queue_cmake flock "$cmakelock" "$cmake_bin" --build "$build" -j "$jobs" --target line-cli
# `ldes` is the simulation engine binary that common/ldes ships. It used to be a
# GraalVM image of the JAVA engine, fetched or built out of ldes/; it is now
# built here from cpp/src/cli/ldes_cli.cpp, so the same make.sh that produces
# line-cli produces it and the two cannot drift apart.
queue_cmake flock "$cmakelock" "$cmake_bin" --build "$build" -j "$jobs" --target ldes
# `line-examples` is the C++ twin of the reference example scripts, run by the
# name the script carries. The parity CPP row falls back to it for the examples
# whose golden holds a quantity the example DERIVES -- a state probability, a
# reward, a CDF-derived mean -- because those are printed by the example and by
# no `-a` view, so line-cli alone leaves the row nothing to compare and it went
# out as a SKIP. Installed here for the same reason line-cli is: a binary the
# harness resolves to must be built by the build everyone runs.
queue_cmake flock "$cmakelock" "$cmake_bin" --build "$build" -j "$jobs" --target line-examples

# EVERY build installs, unoptimized included. common/line-cli is what MATLAB
# lang='cpp', the Python bridge and the parity harness all resolve to, and
# leaving it to a mode nobody runs is how it went stale by weeks and
# manufactured failures that read as defects in the code under test. A slow
# binary that matches the source beats a fast one that does not. common/ldes is
# on the same footing, for MATLAB's and Python's SolverLDES.
#
# The install is STAGED AND RENAMED, never written onto the live path: `install`
# truncates the destination inode, so a parity row or a MATLAB lang='cpp' solve
# that execs common/line-cli during the write gets a half-written binary (or
# ETXTBSY), and a build then poisons any sweep running beside it. The rename is
# atomic within the filesystem, so a concurrent reader sees the old binary or
# the new one and never a partial one.
#
# THE ONE CALLER THAT MUST NOT INSTALL IS THE C++ TEST PHASE (-I), and the reason
# is the 2026-08-27 build/test split rather than anything about this build. Once
# BUILDING became its own phase, phases 3 to 8 START while holding nothing but
# what phases 1 and 2 left in common/, so a TEST phase that installs replaces the
# binaries its sibling phases are in the middle of measuring. Observed on
# 2026-09-07: the C++ BUILD phase installed line-cli, ldes and line-examples on
# picard01 at 08:10:28, wave 2 began at 08:14:43, and the C++ TEST phase on
# picard04 installed its own build of the same three at 08:15:30 -- 47s into the
# python and both parity phases, which resolve to exactly those paths. The rename
# is atomic, so nobody reads a half-written binary; what they read is a DIFFERENT
# binary from the one the rows before them read, built by another host, and a run
# that changes its engine halfway measures neither.
if [ "$install_common" = 1 ]; then
    mkdir -p "$common"
    install -m 0755 "${build}/line-cli" "${common}/.line-cli.new"
    mv -f "${common}/.line-cli.new" "${common}/line-cli"
    install -m 0755 "${build}/ldes" "${common}/.ldes.new"
    mv -f "${common}/.ldes.new" "${common}/ldes"
    install -m 0755 "${build}/line-examples" "${common}/.line-examples.new"
    mv -f "${common}/.line-examples.new" "${common}/line-examples"
fi

if [ "$tests" = ON ]; then
    queue_cmake flock "$cmakelock" "$cmake_bin" --build "$build" -j "$jobs" --target line_mp_tests
    if [ "$run_tests" = 1 ]; then
        queue_run "${build}/line_mp_tests"
    else
        echo "line_mp_tests -> ${build}/line_mp_tests   (built, not run: -n)"
    fi
fi

if [ "$install_common" != 1 ]; then
    # The log line a later reader uses to tell WHICH build a consumer executed,
    # so it must not name common/ when nothing was written there.
    echo "line-cli -> ${build}/line-cli   (-I: common/ left as the BUILD phase installed it)"
elif [ "$quick" = 1 ]; then
    echo "line-cli -> ${common}/line-cli   (unoptimized -O0, the rule in this tree)"
    echo "  build dir: ${build}   (export LINE_CLI=${build}/line-cli to pin it)"
elif [ "$o2" = 1 ]; then
    echo "line-cli -> ${common}/line-cli   (-O2, no -march=native: the suite build)"
    echo "  build dir: ${build}   (export LINE_CLI=${build}/line-cli to pin it)"
else
    echo "line-cli -> ${common}/line-cli   (optimized -O3: release build)"
fi
