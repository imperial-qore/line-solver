#ifndef LINE_SOLVERS_NC_SOLVER_NC_SPN_H
#define LINE_SOLVERS_NC_SOLVER_NC_SPN_H

/**
 * @file solver_nc_spn.h
 * @brief Stationary analysis of a PRODUCT-FORM stochastic Petri net by MDD-rec.
 *
 * The normalising constant is obtained from one memoised walk of the decision
 * diagram holding the reachable set, and every reported measure is a masked walk
 * of the same diagram.
 *
 * This is the `rec` method of SolverNC, and the first analytical route LINE
 * offers for a Petri net -- CTMC solves the explicit generator, SSA and LDES
 * simulate, FLD fluidises. Three functions do the work and each is the subject
 * of its own reference:
 *
 *   `spn_pf`      decides the product form and derives the per-place factors
 *                 g_l (Coleman-Henderson-Taylor complex balance)
 *   `mdd_rec`     G = sum_S prod_l g_l(s_l) in O(sum_l nodes_l * |S_l|) rather
 *                 than O(|S|)                        (Balsamo-Marin-Stojic)
 *   `spn_metrics` mean tokens, place and mode utilisation, and throughputs, all
 *                 from masked walks of the same diagram
 *
 * WHAT THIS REACHES THAT THE EXPLICIT GENERATOR DOES NOT. The diagram stores the
 * reachable set, never the generator, so the cost is set by the number of
 * diagram nodes and not by |S|. It also does not need the marking to be a
 * conserved job population: a mode may consume two tokens and produce one, or
 * consume one and produce two, which is the fork-join and batch case that the
 * MDD-rec paper exists to serve.
 *
 * UN FOLLOWS LINE, NOT THE PAPER. A Place is an INF station, and LINE reports
 * U = Q at an infinite server, which is what SolverCTMC returns for the same
 * net. The paper's place utilisation u(P_j) = 1 - P(m_j = 0) is a different
 * quantity and rides on the returned metrics block instead.
 *
 * ARITHMETIC. `spn_pf` needs a transcendental field and says so; everything
 * downstream of it -- `mdd_rec`, `spn_metrics` -- is rational.
 *
 * @see spn_pf, mdd_rec, spn_metrics, spn_mdd
 */

#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

#include "line/api/spn/spn_metrics.h"
#include "line/api/spn/spn_pf.h"
#include "line/lang/qn/network_struct.h"
#include "line/solvers/nc/nc_types.h"
#include "line/util/matrix.h"

namespace line {
namespace nc {

/** The product-form solve, with the certificate that produced it. */
template <class T>
struct NcSpnSolution {
    NcSolution<T> sol;
    spn::SpnPfResult<T> pf;
    spn::SpnMetrics<T> metrics;
};

/**
 * Analyse a product-form stochastic Petri net.
 *
 * @param sn  the struct of a net holding Places and Transitions
 * @param opt solver controls; `tol` sets the product-form checks' tolerance
 */
template <class T>
NcSpnSolution<T> solver_nc_spn_analyzer(const qn::NetworkStruct<T>& sn, const NcSolverOptions& opt) {
    const T zero = num_traits<T>::from_int(0);
    const std::size_t M = sn.nstations, R = sn.nclasses;

    spn::SpnPfOptions pfopt;
    if (opt.tol > 0) pfopt.tol = std::max(opt.tol, 1e-12);
    spn::SpnPfResult<T> pf = spn::spn_pf<T>(sn, pfopt);
    const spn::SpnMetrics<T> met =
        spn::spn_metrics<T>(pf.spn.mdds, pf.g, pf.spn.info);

    NcSpnSolution<T> out;
    out.pf = pf;
    out.metrics = met;
    out.sol.actualmethod = "rec";
    out.sol.sol.method = "rec";
    out.sol.sol.iter = 1;
    out.sol.sol.Q = Matrix<T>(M, R, zero);
    out.sol.sol.U = Matrix<T>(M, R, zero);
    out.sol.sol.R = Matrix<T>(M, R, zero);
    out.sol.sol.Tp = Matrix<T>(M, R, zero);
    out.sol.sol.X.assign(R, zero);
    out.sol.sol.C.assign(R, zero);

    for (std::size_t pp = 0; pp < pf.spn.info.places.size(); ++pp) {
        const std::size_t ist = sn.nodes[pf.spn.info.places[pp] - 1].station;
        if (ist < 1) continue;
        out.sol.sol.Q(ist - 1, 0) = met.tokens[pp];
        // INF station: LINE charges one server per resident token, so U = Q. The
        // paper's 1 - P(m = 0) is met.place_util, on the returned metrics block.
        out.sol.sol.U(ist - 1, 0) = met.tokens[pp];
        out.sol.sol.Tp(ist - 1, 0) = met.place_tput[pp];
        if (met.place_tput[pp] > zero)
            out.sol.sol.R(ist - 1, 0) = T(met.tokens[pp] / met.place_tput[pp]);
    }

    // System throughput at the reference station, and the response time Little's
    // law then fixes. A net whose class population is not conserved has no
    // meaningful N/X, so C stays zero there rather than reporting a ratio
    // against a moving population.
    const std::size_t ref = sn.classes.empty() ? 0 : sn.classes[0].refstat;
    if (ref >= 1 && ref <= M) out.sol.sol.X[0] = out.sol.sol.Tp(ref - 1, 0);
    T Nk = zero;
    for (std::size_t i = 0; i < M; ++i) Nk += out.sol.sol.Q(i, 0);
    if (out.sol.sol.X[0] > zero && Nk > zero) out.sol.sol.C[0] = T(Nk / out.sol.sol.X[0]);

    out.sol.sol.lG = std::log(num_traits<T>::to_double(met.G));
    return out;
}

}  // namespace nc
}  // namespace line

#endif  // LINE_SOLVERS_NC_SOLVER_NC_SPN_H
