/**
 * Decision-diagram state-space storage and Miner-Ciardo-Donatelli aggregation.
 *
 * <p>{@link jline.api.mdd.MDD} stores a CTMC reachable set in a quasi-reduced
 * ordered multi-valued decision diagram instead of an explicit state list, and
 * {@link jline.api.mdd.Mdd_mcd} solves the chain by level aggregation without
 * ever forming the |S|-state generator. {@link jline.api.mdd.Mdd_rec} is the
 * exact counterpart for a model that already has a product form: one memoised
 * walk of the same diagram returns its normalising constant, and the same walk
 * under a per-level mask returns the marginals.</p>
 *
 * <p>References:</p>
 * <ul>
 *   <li>A.S. Miner, G. Ciardo, "Efficient Reachability Set Generation and
 *       Storage Using Decision Diagrams", ICATPN 1999, LNCS 1639, pp.6-25.</li>
 *   <li>A.S. Miner, G. Ciardo, S. Donatelli, "Using the exact state space of a
 *       Markov model to compute approximate stationary measures", ACM
 *       SIGMETRICS 2000, pp.207-216.</li>
 *   <li>S. Balsamo, A. Marin, I. Stojic, "Computation of the normalising
 *       constant for product-form models of distributed systems with
 *       synchronisation", Future Generation Computer Systems 111 (2020)
 *       475-490.</li>
 * </ul>
 *
 * <p>MATLAB twin: {@code matlab/src/api/mdd/}. Python twin:
 * {@code python/line_solver/api/mdd/}.</p>
 */
package jline.api.mdd;
