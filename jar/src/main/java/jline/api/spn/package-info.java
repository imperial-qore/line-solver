/**
 * Stochastic Petri net analysis.
 *
 * <p>Algorithms specific to the SPN formalism, as opposed to the
 * formalism-agnostic decision-diagram machinery in {@link jline.api.mdd}.
 * {@link jline.api.spn.Spn_mdd} turns a Petri net into the reachable set and
 * Kronecker rate descriptor that {@link jline.api.mdd.Mdd_mcd} consumes.
 * {@link jline.api.spn.Spn_sinvariants}, {@link jline.api.spn.Spn_conv},
 * {@link jline.api.spn.Spn_rec_enabled} and {@link jline.api.spn.Spn_metrics}
 * carry the product-form side: the invariant basis, the convolution over its
 * load vector, and the stationary measures read off the MDD-rec masses.</p>
 *
 * <p>References:</p>
 * <ul>
 *   <li>A.S. Miner, G. Ciardo, "Efficient Reachability Set Generation and
 *       Storage Using Decision Diagrams", ICATPN 1999, LNCS 1639, pp.6-25.</li>
 *   <li>S. Balsamo, A. Marin, I. Stojic, "Computation of the normalising
 *       constant for product-form models of distributed systems with
 *       synchronisation", Future Generation Computer Systems 111 (2020)
 *       475-490.</li>
 *   <li>J.L. Coleman, W. Henderson, P.G. Taylor, "Product form equilibrium
 *       distributions and a convolution algorithm for stochastic Petri nets",
 *       Performance Evaluation 26(3), 1996, 159-180.</li>
 * </ul>
 *
 * <p>MATLAB twin: {@code matlab/src/api/spn/}. Python twin:
 * {@code python/line_solver/api/spn/}.</p>
 */
package jline.api.spn;
