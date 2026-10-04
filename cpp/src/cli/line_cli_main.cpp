/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * `line-cli`'s entry point, and nothing else.
 *
 * ITS OWN TRANSLATION UNIT AND NOT AN `#ifdef` IN line_cli.cpp, because an
 * object file is indivisible. `line_cli.cpp` is compiled once into the static
 * library `line_cli_core`, which the C ABI links so that a host can reach the
 * whole CLI vocabulary in-process; a library carrying `main` would drag that
 * symbol into the shared object R loads, beside the R interpreter's own `main`.
 * A guard macro would have meant compiling that 10,700-line, ~15-minute,
 * multi-gigabyte translation unit TWICE, once per definition of the macro.
 *
 * So the split costs one file and zero rebuilds, and the file is this short on
 * purpose: everything it used to do now lives in `line::cli::main_body`, where
 * the in-process caller runs it too and the two cannot drift.
 */

#include "line/cli/cli_run.h"

int main(int argc, char** argv) { return line::cli::main_body(argc, argv); }
