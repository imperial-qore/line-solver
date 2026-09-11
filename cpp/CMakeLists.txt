# Copyright (c) 2012-2026, QORE Lab, Imperial College London
# All rights reserved.
#
# Multiprecision C++ port of the LINE API layer.
#
# Default backends are header-only Boost.Multiprecision (BSL-1.0), so the
# default build carries no LGPL obligation. -DLINE_MP_USE_GMP=ON switches the
# exact and high-precision backends to GMP/MPFR, which are faster but LGPL.

cmake_minimum_required(VERSION 3.16)

# CMake 3.30 removed the bundled FindBoost module; ask for BoostConfig.cmake
# directly where the policy exists, and keep module mode on older CMake.
if(POLICY CMP0167)
  cmake_policy(SET CMP0167 NEW)
endif()

project(line_mp VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# FLOATING POINT IS NOT REASSOCIATED HERE, and that is a correctness setting
# rather than a preference. GCC and Clang default to -ffp-contract=fast at -O2
# and above, which fuses a*b+c into an FMA ACROSS statements; two expression
# trees computing the same quantity then differ in the last ulp depending on
# which one the optimizer happened to contract. Measured on this port with
# -march=native: `qsys_gg1`'s general branch stopped comparing equal to
# `qsys_gig1_approx_allencunneen` (both 4.75), and `pfqn_clw` drifted from its
# documented 1e-11 aliasing floor to 9.2e-10 against the MATLAB constant -- an
# apparent regression in code nobody had touched. The tolerances in this port
# are quoted against MATLAB, which does not contract, so the port must not
# either. This costs a little speed and buys reproducibility against the
# reference; it does NOT disable vectorization or -march=native.
add_compile_options($<$<COMPILE_LANGUAGE:CXX>:-ffp-contract=off>)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)   # keeps a pybind11 module linkable
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

option(LINE_MP_USE_GMP "Use GMP/MPFR backends instead of Boost (LGPL)" OFF)
option(LINE_MP_USE_LAPACK "Link LAPACK for eigenvalues and singular values (BSD-3)" ON)
option(LINE_MP_USE_HIGHS "Link HiGHS for the sparse simplex and convex QP (MIT)" ON)
option(LINE_MP_USE_SYMENGINE "Link SymEngine for exact SYMBOLIC scalars (MIT)" OFF)
option(LINE_MP_BUILD_TESTS "Build the doctest suite" ON)
option(LINE_MP_BUILD_BENCH "Build the arithmetic-backend benchmarks" ON)
option(LINE_MP_BUILD_EXAMPLES "Build the example gallery binary" ON)

find_package(Boost REQUIRED)
# glibc 2.34 (Ubuntu 22.04) folded libpthread INTO libc, so a host that new links
# std::thread with no flag at all and hides this. glibc 2.31 (20.04) still ships
# a separate libpthread and its ld refuses the implicit dependency:
#   undefined reference to symbol 'pthread_create@@GLIBC_2.2.5'
#   /lib/x86_64-linux-gnu/libpthread.so.0: DSO missing from command line
# Carried on the INTERFACE library so every consumer inherits it.
find_package(Threads REQUIRED)

add_library(line_mp INTERFACE)
target_include_directories(line_mp INTERFACE
  ${CMAKE_CURRENT_SOURCE_DIR}/include
  ${CMAKE_CURRENT_SOURCE_DIR}/third_party
  ${Boost_INCLUDE_DIRS})
target_compile_features(line_mp INTERFACE cxx_std_17)
target_link_libraries(line_mp INTERFACE Threads::Threads)

if(LINE_MP_USE_GMP)
  find_library(GMP_LIB gmp REQUIRED)
  find_library(MPFR_LIB mpfr REQUIRED)
  target_compile_definitions(line_mp INTERFACE LINE_MP_USE_GMP)
  target_link_libraries(line_mp INTERFACE ${GMP_LIB} ${MPFR_LIB})
  message(STATUS "line_mp: GMP/MPFR backends (LGPL) enabled")
else()
  message(STATUS "line_mp: Boost.Multiprecision backends (BSL-1.0)")
endif()

if(LINE_MP_USE_LAPACK)
  find_library(LAPACK_LIB lapack)
  if(LAPACK_LIB)
    target_compile_definitions(line_mp INTERFACE LINE_MP_HAVE_LAPACK)
    target_link_libraries(line_mp INTERFACE ${LAPACK_LIB})
    message(STATUS "line_mp: LAPACK enabled (BSD-3) at ${LAPACK_LIB}")
  else()
    message(STATUS "line_mp: LAPACK not found; eig_values, svd_values, svd_full, "
                   "schur_decomposition, schur_reorder and sylvester_schur refuse by "
                   "name, which takes the fluid second-moment methods and the doctest "
                   "suite with them. Install liblapack-dev (apt), lapack-devel (dnf) "
                   "or lapack (brew, keg-only: export CMAKE_PREFIX_PATH=$(brew --prefix "
                   "lapack)).")
  endif()
endif()

# HiGHS (MIT) supplies a SPARSE simplex and a convex QP. Both are needed by the
# mapqn family: the QRF bound models run to tens of thousands of columns, which
# the built-in dense tableau in util/simplex.h cannot carry, and the QRF NLP
# entry points need a polytope projection for their prox step. Optional, and
# absent it the affected entry points refuse by name rather than run slowly.
# DO NOT use find_package(HIGHS CONFIG) here. Debian's libhighs-dev 1.10.0+ds-1
# ships a broken export: highs-targets.cmake sets the imported location of
# highs::highs to "/debian/tmp/usr/lib/x86_64-linux-gnu/libhighs.so.1.10.0", the
# packaging build root, which does not exist on an installed system. The failure
# is a FATAL_ERROR raised from inside the config file, so merely CALLING
# find_package -- even QUIET -- aborts configuration of the whole project rather
# than falling through to a fallback. Probing for the header and the library
# directly is both sufficient and immune to that bug.
if(LINE_MP_USE_HIGHS)
  find_library(HIGHS_LIB highs)
  find_path(HIGHS_INC Highs.h PATH_SUFFIXES highs)
  if(HIGHS_LIB AND HIGHS_INC)
    target_compile_definitions(line_mp INTERFACE LINE_MP_HAVE_HIGHS)
    target_include_directories(line_mp INTERFACE ${HIGHS_INC})
    target_link_libraries(line_mp INTERFACE ${HIGHS_LIB})
    message(STATUS "line_mp: HiGHS enabled (MIT) at ${HIGHS_LIB}")
  else()
    message(STATUS "line_mp: HiGHS not found; the dense simplex caps the QRF "
                   "bounds at a few hundred columns and the QRF NLP entry "
                   "points refuse by name, and the doctest suite cannot pass. "
                   "Ubuntu has no package: build github.com/ERGO-Code/HiGHS "
                   "with -DBUILD_SHARED_LIBS=ON, or install libhighs-dev "
                   "where it is packaged.")
  endif()
endif()

# SYMENGINE, the exact symbolic scalar behind num_traits<Sym>.
#
# OFF BY DEFAULT, unlike HiGHS, and the difference is deliberate: HiGHS decides
# how well an existing method solves, so a host without it computes a worse
# answer to the same question and the suite has to say so. SymEngine ADDS a
# numeric type and nothing that exists today asks for it, so a host without it
# builds exactly the binary it always did. Turn it on where symbolic input is
# wanted; nothing else changes when it is off.
#
# BUILD IT WITH -DINTEGER_CLASS=boostmp. That drops the GMP dependency entirely
# and leaves Boost -- which this project already requires for cpp_rational -- as
# the only one, keeping every licence permissive (MIT, BSD-3, Boost) and matching
# what the roster already carries. Ubuntu packages no libsymengine-dev, so this
# is a source build into /usr/local, exactly as HiGHS 1.15.1 is.
if(LINE_MP_USE_SYMENGINE)
  find_library(SYMENGINE_LIB symengine)
  find_path(SYMENGINE_INC symengine/basic.h)
  if(SYMENGINE_LIB AND SYMENGINE_INC)
    target_compile_definitions(line_mp INTERFACE LINE_MP_USE_SYMENGINE)
    target_include_directories(line_mp INTERFACE ${SYMENGINE_INC})
    target_link_libraries(line_mp INTERFACE ${SYMENGINE_LIB})
    message(STATUS "line_mp: SymEngine enabled (MIT) at ${SYMENGINE_LIB}")
  else()
    message(FATAL_ERROR "line_mp: LINE_MP_USE_SYMENGINE was asked for and no "
                        "SymEngine was found. Ubuntu packages none: build "
                        "github.com/symengine/symengine with "
                        "-DINTEGER_CLASS=boostmp -DBUILD_SHARED_LIBS=OFF, which "
                        "needs only Boost. Or leave the option OFF, which is "
                        "the default and costs nothing.")
  endif()
endif()

# The --api dispatch table, compiled once: three arithmetics times a dozen
# templated algorithms is too much instantiation to repeat in the CLI, the test
# binary and a future pybind11 module.
add_library(line_mp_api STATIC src/reg/api_dispatch.cpp src/api/mam/iltcme_table.cpp
                              src/io/line_citations.cpp
                              src/solvers/solver_facade.cpp)
target_link_libraries(line_mp_api PUBLIC line_mp)

# ONE CLI for both model kinds, as jline.cli.LineCLI is: a Network model
# (-i json) and a layered one (-i lqnx|xml) enter through the same binary. The
# layered path was a separate line-ln executable until it was folded in.
add_executable(line-cli src/cli/line_cli.cpp)
target_link_libraries(line-cli PRIVATE line_mp_api)

# `ldes`, the simulation engine as its own binary. It is SEPARATE from line-cli
# on purpose: this is the artefact `common/ldes` ships, and MATLAB, native
# Python and the C++ client all run it as a subprocess with the argument vector
# and result document of jline.cli.LdesCLI. Folding it into line-cli would make
# every client carry the whole analytical stack to run a simulation, and would
# tie the engine's interface to line-cli's, which answers to nobody but itself.
add_executable(ldes src/cli/ldes_cli.cpp)
target_link_libraries(ldes PRIVATE line_mp)

# A DISTRIBUTED TREE CARRIES NO tests/, and must still configure. The suite is
# maintained here and stripped from every release archive by
# release/pre-release.sh, as the MATLAB, Java and Python suites already were.
# The glob is therefore evaluated BEFORE the option is consulted: an empty one
# means this is such a tree, and add_executable with no sources is a CMake ERROR
# ("No SOURCES given to target"), so an unguarded block would make a released
# tree refuse to configure for anyone running cmake by hand. cpp/make.sh passes
# LINE_MP_BUILD_TESTS=OFF, so it never saw this; a user's IDE would have.
file(GLOB LINE_MP_TEST_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/tests/*.cpp)
if(LINE_MP_BUILD_TESTS AND LINE_MP_TEST_SOURCES)
  enable_testing()
  add_executable(line_mp_tests ${LINE_MP_TEST_SOURCES})
  target_link_libraries(line_mp_tests PRIVATE line_mp_api)
  target_compile_definitions(line_mp_tests PRIVATE
    LINE_MP_CLI_BINARY="$<TARGET_FILE:line-cli>"
    LINE_MP_REPO_ROOT="${CMAKE_CURRENT_SOURCE_DIR}/..")
  add_dependencies(line_mp_tests line-cli)
  add_test(NAME line_mp_tests COMMAND line_mp_tests)
endif()

# bench/ is stripped from the release archives for the same reason tests/ is, so
# the same guard applies: the option alone would name a source file that a
# distributed tree does not have.
if(LINE_MP_BUILD_BENCH AND EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/bench)
  add_executable(bench_backends bench/bench_backends.cpp)
  target_link_libraries(bench_backends PRIVATE line_mp)
  # Wall-clock of the native LDES engine, for comparison against the Java one
  # on the same model and the same completion budget.
  add_executable(bench_ldes_engine bench/bench_ldes_engine.cpp)
  target_link_libraries(bench_ldes_engine PRIVATE line_mp)
endif()

# The example gallery: `matlab/examples` and `python/examples` as one binary
# that runs an example by the name its script carries. One translation unit per
# reference DIRECTORY, since each instantiates the solvers it calls.
if(LINE_MP_BUILD_EXAMPLES)
  file(GLOB_RECURSE LINE_MP_EXAMPLE_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/examples/*.cpp)
  add_executable(line-examples ${LINE_MP_EXAMPLE_SOURCES})
  target_include_directories(line-examples PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/examples)
  target_compile_definitions(line-examples PRIVATE
    LINE_EXAMPLES_REPO_ROOT="${CMAKE_CURRENT_SOURCE_DIR}/..")
  target_link_libraries(line-examples PRIVATE line_mp_api)
endif()
