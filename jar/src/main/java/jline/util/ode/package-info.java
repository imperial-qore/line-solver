/**
 * Ordinary-differential and differential-algebraic integrators.
 *
 * <p>{@link jline.util.ode.Rodas} is a transliteration of Hairer and Wanner's
 * RODAS, the index-1 DAE integrator the fluid {@code dae} method runs its
 * transient on. It is here rather than under {@code jline.solvers.fluid}
 * because nothing about it is specific to a queueing model, and because the
 * MATLAB, C++ and native Python ports of the same Fortran sit in the same
 * general position in their own trees.
 *
 * <p>The plain-ODE workhorse of the fluid solver remains LSODA, reached through
 * {@link jline.solvers.fluid.LSODAExt}; it solves {@code y' = f} and cannot
 * carry a mass matrix, which is why RODAS exists alongside it rather than
 * replacing it.
 */
package jline.util.ode;
