#pragma once
/**
 * @file fluid_nonhyperbolic.h
 * @brief The one exception the fluid fallback ladder catches.
 *
 * Split out of `fluid_moments.h` on 2026-08-28 so that `solver_fluid.h` can
 * raise it too. `fluid_moments.h` includes `solver_fluid.h`, so the definition
 * could not stay there and be visible in the window loop: the conservation
 * guard raises it from inside the integration, where a divergence is first
 * detectable, and `fluid_runner.h` catches it in exactly the same place as
 * before.
 */

#include <string>

#include "line/util/error.h"

namespace line {
namespace fluid {

/**
 * Raised when the moment closure cannot serve this model: the linearization at
 * the fixed point is not hyperbolic, or the drift has left the simplex. A
 * non-hyperbolic fixed point cannot be detected before the mean is solved, so
 * `fluid_minnormal_applicable` cannot decline it in advance; the runner catches
 * THIS exception, and only this one, to fall back to a first-order method when
 * the moment closure was RESOLVED from `default` rather than asked for by name.
 * Catching a plain NumericError there would also swallow a singular Jacobian or
 * a failed integration, which are defects and not model properties.
 */
class FluidNonHyperbolicError : public NumericError {
public:
    explicit FluidNonHyperbolicError(const std::string& what) : NumericError(what) {}
};

}  // namespace fluid
}  // namespace line
