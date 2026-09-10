/**
 * SolverCTMC's memory pre-gate: the estimator and the refusal policy.
 *
 * The policy under test is REFUSAL, not warning. A model whose explicit state
 * space certainly cannot fit must be stopped before any generation, because the
 * alternative observed on 2026-08-01 was the OOM killer taking the whole
 * process and the terminal hosting it, with no diagnostic at all.
 */

#include "doctest.h"

#include <cmath>

#include "line/api/mc/ctmc_memory_gate.h"
#include "line/api/mc/ctmc_state_space_logsize.h"

using namespace line;

TEST_CASE("ctmc_memory_gate refuses what cannot fit and admits what can") {
    // exp(200) states is the `intractableCTMC` fixture's scale (8 PS stations,
    // N=400, Erlang-5). No host has this memory, so the verdict cannot depend
    // on the machine the test runs on.
    const mc::CtmcGateResult big = mc::ctmc_memory_gate(200.0, false);
    CHECK(big.ok == false);
    CHECK(big.message.find("exceeds the safe budget") != std::string::npos);

    // A handful of states fits anywhere, including the 1 GB fallback budget the
    // probe falls back to when it cannot read the host.
    const mc::CtmcGateResult small = mc::ctmc_memory_gate(std::log(50.0), false);
    CHECK(small.ok == true);
    CHECK(small.message.empty());
}

TEST_CASE("force downgrades the refusal to a warning but keeps the message") {
    const mc::CtmcGateResult forced = mc::ctmc_memory_gate(200.0, true);
    CHECK(forced.ok == true);
    // The message still has to be produced: forcing means the caller accepted a
    // known risk, not that the risk went unreported.
    CHECK(forced.message.find("exceeds the safe budget") != std::string::npos);
}

TEST_CASE("the refusal names a route that actually exists") {
    // 'mdd' is exempt from this gate by construction, so it is the one method
    // the message can honestly recommend for an oversized explicit space.
    const mc::CtmcGateResult big = mc::ctmc_memory_gate(200.0, false);
    CHECK(big.message.find("'mdd'") != std::string::npos);
    CHECK(big.message.find("force=true") != std::string::npos);
}

TEST_CASE("the gate compares in log space, so an overflowing prediction still ranks") {
    // exp(1e4) overflows a double in linear space. A gate that compared byte
    // counts would see inf on both sides of every such model and lose the
    // margin it exists to report; these two must still be ordered.
    const mc::CtmcGateResult a = mc::ctmc_memory_gate(1.0e4, false);
    const mc::CtmcGateResult b = mc::ctmc_memory_gate(1.0e5, false);
    CHECK(a.ok == false);
    CHECK(b.ok == false);
    CHECK(std::isfinite(a.budget_gb));
    CHECK(std::isfinite(b.budget_gb));
}

TEST_CASE("a safety fraction of zero refuses everything with a positive size") {
    const mc::CtmcGateResult r = mc::ctmc_memory_gate(std::log(2.0), false, 0.0);
    CHECK(r.ok == false);
}
