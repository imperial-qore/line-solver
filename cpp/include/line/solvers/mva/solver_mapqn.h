#pragma once
/**
 * @file solver_mapqn.h
 * @brief SolverMVA method 'amva.mapqn': the horizontal-cut mean value analysis (mapqn_amva) of
 * a closed multiclass model with one exponential delay station and one FCFS single-server
 * queue whose class-r service is a MAP.
 *
 * `mva_mapqn_reason` is the one structural predicate behind the method: list_valid_methods
 * drops the name when it is nonempty and `solver_mapqn` raises it, so the offered and run
 * answers cannot drift apart. Port of matlab/src/solvers/MVA/mva_mapqn_reason.m and
 * solver_mva_mapqn_analyzer.m.
 */
#include <cmath>
#include <cstddef>
#include <limits>
#include <string>
#include <vector>

#include "line/api/mapqn/mapqn_amva.h"
#include "line/lang/distribution.h"
#include "line/lang/qn/network_struct.h"
#include "line/lang/qn/solver_feature_sets.h"
#include "line/solvers/mva/mva_types.h"
#include "line/util/error.h"

namespace line {
namespace mva {

namespace mapqn_detail {
template <class T>
bool is_delay(const qn::Station<T>& st) {
    return st.sched == lang::SchedStrategy::INF || std::isinf(st.nservers);
}
inline bool is_markovian(lang::ProcessType t) {
    using lang::ProcessType;
    return t == ProcessType::EXP || t == ProcessType::ERLANG || t == ProcessType::HYPEREXP ||
           t == ProcessType::PH || t == ProcessType::APH || t == ProcessType::COXIAN ||
           t == ProcessType::COX2 || t == ProcessType::MAP || t == ProcessType::MMPP2;
}
}  // namespace mapqn_detail

/** The reason 'amva.mapqn' cannot solve L, or "" when it can (and "" for any other method). */
template <class T>
std::string mva_mapqn_reason(const qn::NetworkStruct<T>& L, const std::string& method) {
    using mapqn_detail::is_delay;
    using mapqn_detail::is_markovian;
    if (qn::mva_base_method(method) != "mapqn") return "";
    const std::size_t M = L.nstations, R = L.nclasses;
    if (M != 2)
        return "solver_mapqn: Method 'amva.mapqn' requires exactly two stations: one delay (infinite server) "
               "and one FCFS queue";
    std::size_t id = 0, ndelay = 0;
    for (std::size_t i = 0; i < M; ++i)
        if (is_delay(L.stations[i])) { id = i; ++ndelay; }
    if (ndelay != 1)
        return "solver_mapqn: Method 'amva.mapqn' requires exactly one delay (infinite-server) station and one queue";
    const std::size_t iq = 1 - id;
    if (L.stations[iq].sched != lang::SchedStrategy::FCFS)
        return "solver_mapqn: Method 'amva.mapqn' requires FCFS scheduling at the queue; station " +
               std::to_string(iq + 1) + " is not FCFS";
    if (L.stations[iq].nservers != 1.0)
        return "solver_mapqn: Method 'amva.mapqn' supports a single-server queue only";
    const std::vector<double> njobs = L.njobs();
    for (std::size_t r = 0; r < R; ++r)
        if (std::isinf(njobs[r])) return "solver_mapqn: Method 'amva.mapqn' supports closed models only";
    if (L.nclosedjobs() <= 0) return "solver_mapqn: Method 'amva.mapqn' supports closed models only";
    const std::size_t nd = L.station_to_node[id], nq = L.station_to_node[iq];
    const T one = num_traits<T>::from_int(1);
    for (std::size_t r = 0; r < R; ++r) {
        if (njobs[r] <= 0) continue;
        if (L.procid(id + 1, r + 1) != lang::ProcessType::EXP)
            return "solver_mapqn: Method 'amva.mapqn' requires exponential think times; class " +
                   std::to_string(r + 1) + " has a non-exponential think time";
        if (!is_markovian(L.procid(iq + 1, r + 1)) || L.service[iq][r].disabled)
            return "solver_mapqn: Method 'amva.mapqn' requires a Markovian (MAP-representable) service process "
                   "at the queue; class " + std::to_string(r + 1) + " is not";
        if (std::fabs(num_traits<T>::to_double(L.route_eff(r + 1, r + 1, nd, nq) - one)) > 1e-12 ||
            std::fabs(num_traits<T>::to_double(L.route_eff(r + 1, r + 1, nq, nd) - one)) > 1e-12)
            return "solver_mapqn: Method 'amva.mapqn' requires every class to cycle delay -> queue -> delay "
                   "without class switching; class " + std::to_string(r + 1) + " does not";
    }
    return "";
}

template <class T>
MvaSolution<T> solver_mapqn(const qn::NetworkStruct<T>& L, const MvaOptions& opt) {
    using mapqn_detail::is_delay;
    const std::string reason = mva_mapqn_reason(L, "mapqn");
    if (!reason.empty()) throw UnsupportedError(reason);
    (void)opt;
    const std::size_t M = L.nstations, R = L.nclasses;
    std::size_t id = 0;
    for (std::size_t i = 0; i < M; ++i)
        if (is_delay(L.stations[i])) id = i;
    const std::size_t iq = 1 - id;
    const std::vector<double> pop = L.njobs();
    std::vector<int> N(R, 0);
    std::vector<T> mu(R, num_traits<T>::from_int(1));
    std::vector<Matrix<T>> D0s(R), D1s(R);
    for (std::size_t r = 0; r < R; ++r) {
        N[r] = static_cast<int>(std::lround(pop[r]));
        if (N[r] > 0) {
            mu[r] = L.rates(id, r);
            const mam::Map<T> m = lang::dist_to_map(L.service[iq][r]);
            D0s[r] = m.D0;
            D1s[r] = m.D1;
        } else {
            D0s[r] = Matrix<T>(1, 1, num_traits<T>::from_int(-1));   // absent class: inert single phase
            D1s[r] = Matrix<T>(1, 1, num_traits<T>::from_int(1));
        }
    }
    const mapqn::MapqnAmvaResult<T> res = mapqn::mapqn_amva<T>(mu, D0s, D1s, N);
    const T zero = num_traits<T>::from_int(0);
    MvaSolution<T> sol;
    sol.Q = Matrix<T>(M, R, zero);
    sol.U = Matrix<T>(M, R, zero);
    sol.R = Matrix<T>(M, R, zero);
    sol.Tp = Matrix<T>(M, R, zero);
    sol.C.assign(R, zero);
    sol.X.assign(R, zero);
    for (std::size_t r = 0; r < R; ++r) {
        if (N[r] <= 0 || num_traits<T>::to_double(res.X[r]) <= 0.0) continue;
        const T X = res.X[r];
        sol.X[r] = X;
        sol.Q(iq, r) = res.Qq[r];
        sol.U(iq, r) = res.U[r];
        sol.Tp(iq, r) = X;
        sol.R(iq, r) = res.Qq[r] / X;
        sol.Q(id, r) = X / mu[r];
        sol.U(id, r) = X / mu[r];
        sol.Tp(id, r) = X;
        sol.R(id, r) = num_traits<T>::from_int(1) / mu[r];
        sol.C[r] = num_traits<T>::from_int(N[r]) / X;
    }
    sol.method = "amva.mapqn";
    int iter = 1;
    for (std::size_t r = 0; r < R; ++r) iter *= (N[r] + 1);
    sol.iter = iter;
    sol.lG = std::numeric_limits<double>::quiet_NaN();
    return sol;
}

}  // namespace mva
}  // namespace line
