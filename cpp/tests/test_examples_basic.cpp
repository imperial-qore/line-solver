/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Regression over the SHIPPED BASIC examples (cpp/examples/basic/), the C++
 * twin of jline.examples.ExampleRunnerBasicTest and of the MATLAB
 * testsExamples suite. Split from the advanced half for the reason Surefire
 * gives each Java runner class its own fork: neither half should dominate the
 * wall-clock of phase 4, and a failure localises to a package.
 *
 * See example_suite.h for what is asserted and why this lives in line-dev.git.
 */
#include "example_suite.h"

/**
 * 129 registrations across 11 source files as of 2026-09-12. The floor is what
 * guards against a link that drops the static registrations; it is not a
 * census, so adding an example needs no edit here.
 */
TEST_CASE("examples/basic: every shipped example runs and produces finite numbers") {
    line_tests::run_examples_under("basic/", 130);
}
