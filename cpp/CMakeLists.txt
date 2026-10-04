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
# STATIC DEPENDENCIES, FOR A BINARY THAT LEAVES THIS MACHINE. Default OFF, which
# is right for every build that runs where it was compiled. A PUBLISHED binary is
# not one of those: `line-cli` links HiGHS, LAPACK and BLAS dynamically, and
# HiGHS is not a distro package on most systems -- this roster builds it from
# source into /usr/local/lib -- so a downloader's fetch succeeds and the exec
# then dies with `error while loading shared libraries: libhighs.so.1`. Turning
# this ON links those, the Fortran runtime and libstdc++/libgcc INTO the binary,
# leaving libc and libm, which every Linux has. See upload-line-cli.sh, which
# builds `make.sh -P -S` and refuses to publish a binary still asking for
# anything else.
option(LINE_MP_STATIC_DEPS "Link HiGHS, LAPACK/BLAS and the C++/Fortran runtimes statically (published binaries)" OFF)
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

# WITH LINE_MP_STATIC_DEPS THE SAME find_library CALLS RESOLVE ARCHIVES, because
# the suffix list is what decides, not the call. Restored below, so nothing after
# this point is affected.
if(LINE_MP_STATIC_DEPS)
  set(_line_mp_saved_suffixes "${CMAKE_FIND_LIBRARY_SUFFIXES}")
  set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
endif()

if(LINE_MP_USE_LAPACK)
  # OpenBLAS FIRST, and deliberately: on this roster `liblapack.so.3` resolves
  # through the alternatives to OpenBLAS, so the archive of the same name is the
  # implementation the dynamic build was already using and the published binary
  # computes what the tested one computes. It is also self-contained -- BLAS and
  # LAPACK in one archive -- whereas netlib's `liblapack.a` carries no BLAS at
  # all: a shared liblapack.so.3 names libblas.so.3 in its DT_NEEDED and an
  # archive names nothing, so the static netlib path has to add it by hand.
  # Ubuntu builds OpenBLAS with DYNAMIC_ARCH, and the archive carries every
  # kernel from ATOM to COOPERLAKE, so linking it statically keeps the RUNTIME
  # CPU dispatch that makes a -P binary safe on an old client.
  if(LINE_MP_STATIC_DEPS)
    find_library(LAPACK_LIB NAMES openblas lapack)
  else()
    find_library(LAPACK_LIB lapack)
  endif()
  if(LAPACK_LIB)
    target_compile_definitions(line_mp INTERFACE LINE_MP_HAVE_LAPACK)
    target_link_libraries(line_mp INTERFACE ${LAPACK_LIB})
    message(STATUS "line_mp: LAPACK enabled (BSD-3) at ${LAPACK_LIB}")
    if(LINE_MP_STATIC_DEPS AND NOT LAPACK_LIB MATCHES "openblas")
      # A netlib archive resolves the LAPACK routines and leaves every BLAS call
      # undefined. Refuse rather than fail at link time with 200 lines of
      # `undefined reference to dgemm_`.
      find_library(BLAS_LIB NAMES blas)
      if(NOT BLAS_LIB)
        message(FATAL_ERROR "line_mp: LINE_MP_STATIC_DEPS is ON and the static LAPACK at "
                            "${LAPACK_LIB} is not OpenBLAS, so a static libblas.a is needed "
                            "beside it and none was found. Install libblas-dev, or install "
                            "libopenblas-dev whose archive carries both.")
      endif()
      target_link_libraries(line_mp INTERFACE ${BLAS_LIB})
      message(STATUS "line_mp: static BLAS at ${BLAS_LIB}")
    endif()
  elseif(LINE_MP_STATIC_DEPS)
    # A STATIC REQUEST THAT CANNOT BE MET IS FATAL, where the ordinary one is a
    # degradation. LINE_MP_STATIC_DEPS is only ever set to publish a binary, and
    # falling through here would publish one that refuses eig_values and the
    # fluid second moments on every client -- a worse outcome than a build that
    # stops and says which package is missing.
    message(FATAL_ERROR "line_mp: LINE_MP_STATIC_DEPS is ON but no LAPACK/BLAS ARCHIVE was "
                        "found (a .so does not count: this build must carry it). Install "
                        "libopenblas-dev, or liblapack-dev plus libblas-dev.")
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
  elseif(LINE_MP_STATIC_DEPS)
    # As for LAPACK above: publishing is the only reason this option is ON, and
    # `libhighs.a` is exactly what a published binary cannot do without, since
    # HiGHS is the dependency almost no client has installed. The roster builds
    # HiGHS from source; add a second configure with -DBUILD_SHARED_LIBS=OFF to
    # get the archive beside the .so (install-line-deps.yml does this).
    message(FATAL_ERROR "line_mp: LINE_MP_STATIC_DEPS is ON but no libhighs.a was found "
                        "(a .so does not count: this build must carry it). Build "
                        "github.com/ERGO-Code/HiGHS with -DBUILD_SHARED_LIBS=OFF and "
                        "install it, or pass -DLINE_MP_USE_HIGHS=OFF to publish without "
                        "the sparse simplex -- which changes what wide QRF bounds cost.")
  else()
    # WHAT THIS SAYS IS WHAT HAPPENS, checked against the code on 2026-09-11.
    # It used to promise that "the QRF NLP entry points refuse by name" and
    # that the dense simplex "caps the QRF bounds at a few hundred columns".
    # NEITHER IS IMPLEMENTED: `line::lp::highs_available()`
    # (util/lp_highs.h:51) has no caller outside its own test, and the primary
    # `lp_solve` template (:175-179) DISCARDS `dense_max_cols` and sends every
    # model to the dense tableau at any width. So a wide QRF bound is SLOW, not
    # refused and not capped, and nothing labels it. A promise of a refusal that
    # does not happen is worse than no promise: it invites shipping a
    # HiGHS-less build believing bad inputs would be caught.
    #
    # Implementing the refusals would change the numbers seven existing
    # consumers already return (api/mapqn/*, api/npfqn/npfqn_bnd_*,
    # api/spn/spn_lpbnd.h), which belongs in a change about those methods and
    # not in a build message. The message is corrected here; the refusals
    # remain an open item.
    message(STATUS "line_mp: HiGHS not found. The dense simplex then handles "
                   "EVERY width rather than handing wide models to the sparse "
                   "solver, so a wide QRF bound is slow rather than refused or "
                   "capped, and nothing labels it as such. The doctest suite "
                   "cannot pass without HiGHS (tests/test_lp_highs.cpp asserts "
                   "the backend is compiled in). Ubuntu has no package: build "
                   "github.com/ERGO-Code/HiGHS with -DBUILD_SHARED_LIBS=ON, or "
                   "install libhighs-dev where it is packaged.")
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

# WHAT A STATIC BUILD STILL HAS TO BE TOLD, once the archives above are chosen.
# An ARCHIVE NAMES NOTHING: a shared library records its own dependencies in
# DT_NEEDED and the loader follows them, while `libopenblas.a` and `libhighs.a`
# record nothing at all, so everything they call has to appear on the link line
# here -- the Fortran runtime for LAPACK, zlib for HiGHS's compressed readers.
# Order matters and this is the order: the archives above, then what they call.
#
# The C++ and GCC runtimes go in too, and that is what removes the GLIBCXX floor
# entirely rather than lowering it: a binary built here demanded GLIBCXX_3.4.26
# and now demands no GLIBCXX at all, so the remaining floor is glibc's and the
# only libraries left are libc and libm. Measured on the full combination
# (static HiGHS + static OpenBLAS + this block): `ldd` reports libm, libc and
# the loader, and nothing else.
if(LINE_MP_STATIC_DEPS)
  # ASK THE COMPILER WHERE ITS OWN RUNTIME IS. libgfortran.a and libquadmath.a
  # live in GCC's private directory (/usr/lib/gcc/<triple>/<ver>/), which is NOT
  # on CMake's default library search path, so find_library misses both on a
  # stock Ubuntu and the static build fails for a package that IS installed.
  # `-print-file-name` is the supported way to ask; it echoes the bare name back
  # when the compiler has no such file, which is what IS_ABSOLUTE tests for.
  set(LINE_MP_FORTRAN_RUNTIME "")
  foreach(_rtlib gfortran quadmath)
    execute_process(COMMAND "${CMAKE_CXX_COMPILER}" "-print-file-name=lib${_rtlib}.a"
                    OUTPUT_VARIABLE _rtpath OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
    if(NOT IS_ABSOLUTE "${_rtpath}" OR NOT EXISTS "${_rtpath}")
      find_library(_rtfound ${_rtlib})
      set(_rtpath "${_rtfound}")
      unset(_rtfound CACHE)
    endif()
    if(NOT _rtpath OR NOT EXISTS "${_rtpath}")
      # NAME THE VERSIONED PACKAGE, because the unversioned one is a trap: on
      # Ubuntu 20.04 `gfortran` is gfortran-9 and ships only .../gcc/<triple>/9/
      # libgfortran.a, while a g++-10 link looks in .../10/ and nowhere else. The
      # whole roster had `gfortran` installed and still failed here (2026-09-18).
      string(REGEX MATCH "^[0-9]+" _gxx_major "${CMAKE_CXX_COMPILER_VERSION}")
      message(FATAL_ERROR "line_mp: LINE_MP_STATIC_DEPS is ON and lib${_rtlib}.a was not "
                          "found for this compiler (GCC ${CMAKE_CXX_COMPILER_VERSION}). "
                          "LAPACK is Fortran and an archive names none of its own "
                          "dependencies, so the Fortran runtime has to be carried too. "
                          "Install the archive MATCHING THE COMPILER -- on Debian/Ubuntu "
                          "that is libgfortran-${_gxx_major}-dev (libquadmath.a comes with "
                          "libgcc-${_gxx_major}-dev). A plain `gfortran` is usually a "
                          "DIFFERENT major version and does not satisfy this.")
    endif()
    list(APPEND LINE_MP_FORTRAN_RUNTIME "${_rtpath}")
  endforeach()
  target_link_libraries(line_mp INTERFACE ${LINE_MP_FORTRAN_RUNTIME})
  message(STATUS "line_mp: static Fortran runtime at ${LINE_MP_FORTRAN_RUNTIME}")

  # zlib is the one dependency here that IS on every Linux, so a shared libz is
  # not a publishing defect -- it is taken statically when the archive is there
  # and left alone when it is not, rather than failing a release over it.
  find_library(Z_STATIC z)
  if(Z_STATIC)
    target_link_libraries(line_mp INTERFACE ${Z_STATIC})
    message(STATUS "line_mp: static zlib at ${Z_STATIC}")
  else()
    message(STATUS "line_mp: no libz.a; the published binary will name libz.so.1, which "
                   "every Linux has. Install zlib1g-dev to carry it as well.")
  endif()

  # THE DYNAMIC LOADER, FOR HiGHS. HighsDynamicLibrary.cpp calls dlopen, dlsym,
  # dlclose and dlerror, and the archive carrying it names no dependency of its
  # own, so -ldl has to be spelled out AFTER it. Invisible on a 22.04 driver,
  # whose glibc 2.34 folded libdl into libc; the 20.04 release host still has it
  # separate and the first dispatched -P -S build died exactly here, four
  # undefined references deep into libhighs.a. Same shape as the thread library
  # below, and CMAKE_DL_LIBS is empty on a platform that needs nothing.
  if(CMAKE_DL_LIBS)
    target_link_libraries(line_mp INTERFACE ${CMAKE_DL_LIBS})
  endif()

  # AND THE THREAD LIBRARY LAST. Threads::Threads is already on line_mp, but it
  # is linked BEFORE these archives and CMake will not reorder it. On a host
  # where the thread library is still separate -- glibc < 2.34, which is exactly
  # the 20.04 release host -- OpenBLAS's own pthread_create would then have to be
  # resolved by a library standing to its left, and Ubuntu links --as-needed. A
  # second, explicitly last -lpthread costs nothing where it is a stub (22.04
  # folded libpthread into libc) and removes the question where it is not.
  if(CMAKE_USE_PTHREADS_INIT)
    target_link_libraries(line_mp INTERFACE -lpthread)
  endif()

  target_link_options(line_mp INTERFACE -static-libstdc++ -static-libgcc)
endif()

# The suffix list is restored whatever the option did, so anything added after
# this point resolves the way it always has.
if(LINE_MP_STATIC_DEPS)
  set(CMAKE_FIND_LIBRARY_SUFFIXES "${_line_mp_saved_suffixes}")
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
#
# THE DISPATCH IS A LIBRARY AND THE BINARY IS TEN LINES. `line_cli.cpp` is the
# only place that serves the whole vocabulary -- every solver, every -a
# analysis, every input format, --api -- and a host that wants any of it had no
# route but to spawn the process. `line_cli_core` makes that dispatch callable
# in-process through include/line/cli/cli_run.h, which is what the C ABI in
# `r/capi/` links.
#
# THE SPLIT IS TWO FILES AND NOT AN `#ifdef`, because an object file is
# indivisible: a static library carrying `main` would drag that symbol into the
# shared object a host loads beside its own `main`. A guard macro would instead
# have meant compiling this translation unit TWICE, and it is the most
# expensive one in the project by a wide margin (measured 2026-09-11 at -O2:
# 17.9 GB peak in cc1plus and a 93 MB object). It is still compiled exactly
# once.
add_library(line_cli_core STATIC src/cli/line_cli.cpp)
target_link_libraries(line_cli_core PUBLIC line_mp_api)

add_executable(line-cli src/cli/line_cli_main.cpp)
target_link_libraries(line-cli PRIVATE line_cli_core)

# `ldes`, the simulation engine as its own binary. It is SEPARATE from line-cli
# on purpose: this is the artefact `common/ldes` ships, and MATLAB, native
# Python and the C++ client all run it as a subprocess with the argument vector
# and result document of jline.cli.LdesCLI. Folding it into line-cli would make
# every client carry the whole analytical stack to run a simulation, and would
# tie the engine's interface to line-cli's, which answers to nobody but itself.
add_executable(ldes src/cli/ldes_cli.cpp)
target_link_libraries(ldes PRIVATE line_mp)

# The example gallery: `matlab/examples` and `python/examples` as one binary
# that runs an example by the name its script carries. One translation unit per
# reference DIRECTORY, since each instantiates the solvers it calls.
#
# SPLIT INTO A LIBRARY AND A main() so the doctest suite can drive the same
# example bodies the binary runs, instead of only observing them through a
# subprocess as parity-static does. It is declared HERE, above the test section,
# because line_mp_tests links it.
#
# OBJECT, NOT STATIC, and the distinction is load-bearing: an example registers
# itself into line::examples::registry() from a static initializer and nothing
# ever names it, so a static archive would let the linker drop every object file
# -- leaving an empty registry and an example suite that passes having run
# nothing. An object library carries all of its objects into each consumer.
# CMake 3.24's $<LINK_LIBRARY:WHOLE_ARCHIVE> expresses the same intent, but the
# floor here is 3.16.
if(LINE_MP_BUILD_EXAMPLES)
  file(GLOB_RECURSE LINE_MP_EXAMPLE_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/examples/*.cpp)
  list(REMOVE_ITEM LINE_MP_EXAMPLE_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/examples/main.cpp)
  add_library(line_examples_lib OBJECT ${LINE_MP_EXAMPLE_SOURCES})
  # PUBLIC on both: a consumer's own translation units include examples_common.h
  # and must resolve example_data_file() against the SAME repo root the examples
  # were compiled with, or a test would read a different data file than the
  # example it is checking.
  target_include_directories(line_examples_lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/examples)
  target_compile_definitions(line_examples_lib PUBLIC
    LINE_EXAMPLES_REPO_ROOT="${CMAKE_CURRENT_SOURCE_DIR}/..")
  target_link_libraries(line_examples_lib PUBLIC line_mp_api)

  # Unchanged behaviour for the shipped binary: main.cpp over the same objects,
  # so common/line-examples and `--record` keep working exactly as before.
  add_executable(line-examples examples/main.cpp)
  target_link_libraries(line-examples PRIVATE line_examples_lib)
endif()

# A DISTRIBUTED TREE CARRIES NO tests/, and must still configure. The suite is
# maintained here and stripped from every release archive by
# release/pre-release.sh, as the MATLAB, Java and Python suites already were.
# The glob is therefore evaluated BEFORE the option is consulted: an empty one
# means this is such a tree, and add_executable with no sources is a CMake ERROR
# ("No SOURCES given to target"), so an unguarded block would make a released
# tree refuse to configure for anyone running cmake by hand. cpp/make.sh passes
# LINE_MP_BUILD_TESTS=OFF, so it never saw this; a user's IDE would have.
#
# THE SUITE MOVED TO line-test.git. cpp/tests/ is empty in a development tree
# now, so the glob below is joined by a second one over the sibling checkout.
# LINE_TEST_DIR overrides the location; a checkout without the sibling simply
# builds no suite, exactly as a released tree does.
#
# The in-tree glob is KEPT rather than deleted: a released archive has neither
# directory, and keeping both makes "no sources" mean the same thing in both
# cases. It also leaves a home for a test that genuinely belongs to cpp/ alone.
if(NOT DEFINED LINE_TEST_DIR)
  if(DEFINED ENV{LINE_TEST_DIR})
    set(LINE_TEST_DIR "$ENV{LINE_TEST_DIR}")
  else()
    set(LINE_TEST_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../line-test.git")
  endif()
endif()
file(GLOB LINE_MP_TEST_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/tests/*.cpp)
# The in-tree example tests link line_examples_lib and include examples_common.h,
# neither of which exists when the gallery is switched off. Drop them by name
# rather than guarding inside the TU: without the target there is no include
# directory either, so an #ifdef would not save the compile.
if(NOT LINE_MP_BUILD_EXAMPLES)
  list(REMOVE_ITEM LINE_MP_TEST_SOURCES
       ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_examples_basic.cpp
       ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_examples_advanced.cpp)
endif()
list(LENGTH LINE_MP_TEST_SOURCES LINE_MP_TEST_LOCAL_COUNT)
file(GLOB LINE_MP_TEST_SOURCES_EXT ${LINE_TEST_DIR}/cpp/tests/*.cpp)
list(APPEND LINE_MP_TEST_SOURCES ${LINE_MP_TEST_SOURCES_EXT})
# Report BOTH halves. This message named only the external count while
# cpp/tests/ was empty, so it read as the whole story; now that the example
# suite lives here it would silently under-report the denominator, which is the
# one way a suite shrinks without anyone noticing.
if(LINE_MP_TEST_SOURCES)
  list(LENGTH LINE_MP_TEST_SOURCES_EXT LINE_MP_TEST_EXT_COUNT)
  list(LENGTH LINE_MP_TEST_SOURCES LINE_MP_TEST_ALL_COUNT)
  message(STATUS "line_mp_tests: ${LINE_MP_TEST_ALL_COUNT} sources "
                 "(${LINE_MP_TEST_EXT_COUNT} from ${LINE_TEST_DIR}/cpp/tests, "
                 "${LINE_MP_TEST_LOCAL_COUNT} in cpp/tests)")
endif()
if(LINE_MP_BUILD_TESTS AND LINE_MP_TEST_SOURCES)
  enable_testing()
  add_executable(line_mp_tests ${LINE_MP_TEST_SOURCES})
  # line_cli_core as well as line_mp_api: tests/test_cli_run.cpp exercises the
  # in-process CLI seam, and it is the suite's job to prove that the binary and
  # the in-process run produce the same bytes.
  target_link_libraries(line_mp_tests PRIVATE line_cli_core line_mp_api)
  # The example bodies, so tests/test_examples_{basic,advanced}.cpp can walk the
  # registry in-process. Guarded: with the gallery off the two TUs above were
  # dropped from the glob and nothing here needs the objects.
  if(TARGET line_examples_lib)
    target_link_libraries(line_mp_tests PRIVATE line_examples_lib)
  endif()
  # LINE_MP_REPO_ROOT still points at line-dev.git: common/line-examples is
  # built here and does NOT move. LINE_MP_TEST_ROOT is where the suite and the
  # shared goldens now live. Both stay overridable at runtime by LINE_EXAMPLES
  # and LINE_GOLDENS respectively.
  target_compile_definitions(line_mp_tests PRIVATE
    LINE_MP_CLI_BINARY="$<TARGET_FILE:line-cli>"
    LINE_MP_REPO_ROOT="${CMAKE_CURRENT_SOURCE_DIR}/.."
    LINE_MP_TEST_ROOT="${LINE_TEST_DIR}")
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
