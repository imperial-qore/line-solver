/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_MAM_SOLVER_MAM_BASIC_H
#define LINE_SOLVERS_MAM_SOLVER_MAM_BASIC_H

/**
 * Port of `solver_mam_basic.m`, the `dec.source` analyzer and the default
 * algorithm of SolverMAM.
 *
 * THE METHOD. Each queueing station is solved IN ISOLATION as a matrix-analytic
 * queue, and the stations are coupled only through a fixed point on the
 * per-chain throughput. The arrival stream a station sees is not derived from
 * traffic equations: it is the chain's SOURCE process (or, for a closed chain,
 * a Poisson surrogate at the current throughput iterate) rescaled to the visit
 * rate of that station. That is what "dec.source" names, and it is why the
 * method is cheap and why it is an approximation.
 *
 * THE FIXED POINT. `lambda(c)` is the chain arrival rate. Open chains have it
 * fixed by their source. Closed chains start at the no-contention lower bound
 * N/sum(D) and are then driven by an ITERATION-AVERAGED regula falsi towards
 * QN = N; the averaging weight walks from the raw Newton step to a no-op as the
 * iteration count approaches iter_max, which is what damps the oscillation a
 * bare N/QN step shows on a saturated chain. A purely open or purely closed
 * network backs off uniformly by 1/Umax when the busiest station saturates; a
 * MIXED network instead backs off the CLOSED chains alone onto the capacity the
 * open traffic leaves free, because scaling the open chains too would report a
 * throughput below their own source rate.
 *
 * WHAT EACH STATION GETS. In branch order, as the reference writes them:
 *   INF                delay: U = Q = T S, R = S
 *   PS                 U = T S / c and the M/M/1-PS-style Q = U/(1-Utot)
 *   FCFS / HOL (the reference also lists FCFSPRPRIO, routed through BUTools'
 *              MMAPPH1PRPR; that analyzer is not ported, so an FCFSPRPRIO
 *              station reaches the station-ladder throw near the end of this
 *              file instead of a branch here -- the C++ SchedStrategy DOES
 *              have an FCFSPRPRIO enumerator, lang_types.h:144, only no MAM
 *              analyzer serves it)
 *     open:  PH/M/c (renewal arrival, exp service, exact), D/M/c, MAP/D/c,
 *            finite buffer (exact M/M/c/K, exact MMAP/G/1/K at one server,
 *            truncate-and-renormalize otherwise), RAP/RAP/1, MAP/MAP/1 for a
 *            correlated single-class service, else MMAP[K]/PH[K]/1 FCFS
 *     closed: the MMAP[K]/PH[K]/1 queue-length DISTRIBUTION truncated at the
 *            class population, which is what makes a closed chain's queue
 *            length respect its own population bound
 *
 * THE SURROGATE DELAY. A c-server station is solved as a single server of c
 * times the speed and the missing T S (c-1)/c jobs are added back afterwards.
 * The branches that are exact for c > 1 (PH/M/c, D/M/c, MAP/D/c, M/M/c/K) skip
 * that correction, which is what `mapdcStations` records.
 *
 * THE SETUP / DELAY-OFF BRANCHES ARE WRITTEN, both of them, keyed on
 * NetworkStruct::setupparam (which IS the reference's sn.hassetup). The open
 * one collapses the station to one M/G/1-with-setup and splits the QBD's queue
 * back by load share; the closed one is not a QBD at all but a cold-start race,
 * charging the setup with the probability that the delay-off timer expired
 * before the job's own think time did. SolverLN routes a setup-task layer
 * here through `dec.poisson`.
 *
 * WHAT THIS PORT DOES NOT REACH, and why it is not a gap. The reference also
 * carries a self-looping-class clamp and an ME/RAP branch.
 * `lang/lang_types.h` has no SelfLoopingClass (JobClassType is OPEN or CLOSED)
 * and no ME or RAP ProcessType, so no model the C++ layer can express reaches
 * either. They are noted rather than written, because writing an unreachable
 * branch is writing an untested one.
 */

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

#include "line/api/mam/aph_fit_moments.h"
#include "line/api/mam/map_moment.h"
#include "line/api/mam/map_transform.h"
#include "line/api/mam/mmap_assemble.h"
#include "line/api/mam/mmap_lambda.h"
#include "line/api/mam/mmapph1fcfs.h"
#include "line/api/mam/qbd_mapmap1.h"
#include "line/api/mam/qbd_setupdelayoff.h"
#include "line/api/mc/ctmc_solve.h"
#include "line/api/sn/sn_get_buffer_size.h"
#include "line/api/qsys/qsys_dmc.h"
#include "line/api/qsys/qsys_mapdc.h"
#include "line/api/qsys/qsys_mapmc.h"
#include "line/api/qsys/qsys_mapphc.h"
#include "line/api/qsys/qsys_mmapgk1.h"
#include "line/api/qsys/qsys_mmapg1k.h"
#include "line/api/qsys/qsys_mmck.h"
#include "line/api/qsys/qsys_phmc.h"
#include "line/lang/distribution.h"
#include "line/lang/qn/network_struct.h"
#include "line/solvers/mam/mam_types.h"
#include "line/solvers/mva/sn_chain.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace mam {

using lang::GlobalConstants;
using lang::ProcessType;
using lang::SchedStrategy;

namespace basic_detail {

/**
 * True for a process the MAM analyzer can read as a (D0,D1) pair as declared.
 *
 * ME AND RAP COUNT. The matrix-geometric algebra never needs D0 to be a
 * generator -- R solves a quadratic in the blocks and the stationary vector
 * follows from it -- so a rational arrival process is as solvable as a
 * Markovian one, which is why `SolverMAM.getFeatureSet` lists both. They are
 * also exactly what `sn_nonmarkov_toph`'s default two-moment fit produces, so
 * excluding them would refuse the conversion the analyzer had just performed.
 * The paths that DO need a generator -- the retrial solver -- test
 * `is_markovian_map` on the pair itself and refuse there.
 */
inline bool is_markovian_type(ProcessType p) {
    switch (p) {
        case ProcessType::EXP:
        case ProcessType::ERLANG:
        case ProcessType::HYPEREXP:
        case ProcessType::PH:
        case ProcessType::APH:
        case ProcessType::ME:
        case ProcessType::RAP:
        case ProcessType::MAP:
        case ProcessType::MMAP:
        case ProcessType::COXIAN:
        case ProcessType::COX2:
        case ProcessType::MMPP2:
        case ProcessType::IMMEDIATE:
        case ProcessType::DISABLED:
            return true;
        default:
            return false;
    }
}

/**
 * True when station `i0` should be answered by MMAP[K]/G[K]/1.
 *
 * The generic mmapph1fcfs path reads the service law out of the PHASE-TYPE FIT:
 * for a Uniform, Gamma, Pareto, Weibull, Lognormal or Det that fit matches the
 * mean and, once the SCV exceeds one, nothing else. He (2001) needs only the
 * TRANSFORM of the original law, which `lang::dist_lst` evaluates directly off
 * `L.service`, so wherever a class declares one of those laws the exact analysis
 * is available. A matrix-exponential service qualifies too, its transform being
 * rational; a RAP does NOT, He's analysis assuming INDEPENDENT service times, so
 * reading a correlated service through its marginal transform would discard
 * exactly the autocorrelation the RAP was declared to carry. Mirrors MATLAB
 * `mam_gk1_applicable`.
 */
template <class T>
const lang::Distrib<T>& mam_declared_law(const lang::Distrib<T>& d) {
    // sn_nonmarkov_toph has already replaced a Gamma or a Lognormal by its
    // phase-type fit and retagged it PH/ME, keeping the original under
    // `declared`. That original is the law whose transform this path reads.
    return d.declared ? *d.declared : d;
}

template <class T>
bool mam_gk1_applicable(const qn::NetworkStruct<T>& L, std::size_t i0, std::size_t K) {
    for (std::size_t r = 0; r < K; ++r) {
        const lang::ProcessType t = mam_declared_law(L.service[i0][r]).type;
        if (t == lang::ProcessType::DET || t == lang::ProcessType::UNIFORM ||
            t == lang::ProcessType::GAMMA || t == lang::ProcessType::PARETO ||
            t == lang::ProcessType::WEIBULL || t == lang::ProcessType::LOGNORMAL ||
            t == lang::ProcessType::ME) {
            return true;
        }
    }
    return false;
}

/**
 * Port of `mam_is_renewal_map`: D1 = (-D0 e) sigma to within tolerance, i.e.
 * the phase after an arrival does not depend on the phase before it.
 */
template <class T>
bool is_renewal_map(const Map<T>& m) {
    const std::size_t n = m.D0.rows();
    if (n == 1) return true;
    double scale = 1.0;
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j)
            scale = std::max(scale, std::fabs(num_traits<T>::to_double(m.D1(i, j))));
    const double tol = 1e-9 * scale;
    std::vector<double> exitrate(n, 0.0), sigma(n, 0.0);
    double tot = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) exitrate[i] += num_traits<T>::to_double(m.D1(i, j));
        tot += exitrate[i];
    }
    if (tot <= 0.0) return true;
    // sigma from the first row with a positive exit rate.
    std::size_t ref = n;
    for (std::size_t i = 0; i < n; ++i)
        if (exitrate[i] > tol) {
            ref = i;
            break;
        }
    if (ref == n) return true;
    for (std::size_t j = 0; j < n; ++j)
        sigma[j] = num_traits<T>::to_double(m.D1(ref, j)) / exitrate[ref];
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j)
            if (std::fabs(num_traits<T>::to_double(m.D1(i, j)) - exitrate[i] * sigma[j]) > tol)
                return false;
    return true;
}

/** True when (D0,D1) is a genuine Markovian pair, as `is_markovian_map` tests. */
template <class T>
bool is_markovian_map(const Map<T>& m) {
    const std::size_t n = m.D0.rows();
    double scale = 1.0;
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) {
            scale = std::max(scale, std::fabs(num_traits<T>::to_double(m.D0(i, j))));
            scale = std::max(scale, std::fabs(num_traits<T>::to_double(m.D1(i, j))));
        }
    const double tol = 1e-9 * scale;
    for (std::size_t i = 0; i < n; ++i) {
        double rowsum = 0.0;
        for (std::size_t j = 0; j < n; ++j) {
            const double a = num_traits<T>::to_double(m.D0(i, j));
            const double b = num_traits<T>::to_double(m.D1(i, j));
            if (i != j && a < -tol) return false;
            if (b < -tol) return false;
            rowsum += a + b;
        }
        if (std::fabs(rowsum) > tol) return false;
    }
    return true;
}

/** The PH mixture `mam_svc_mixture` builds: alpha = [w_k sigma_k], T = blkdiag(S_k). */
template <class T>
qsys::ServiceLaw<T> svc_mixture(const Mmap<T>& arv, const std::vector<PhService<T>>& svc) {
    const T zero = num_traits<T>::from_int(0);
    const std::size_t K = svc.size();
    Matrix<T> Q = arv.D0;
    for (std::size_t k = 0; k < arv.classes(); ++k)
        for (std::size_t i = 0; i < Q.rows(); ++i)
            for (std::size_t j = 0; j < Q.cols(); ++j) Q(i, j) += arv.Dc[k](i, j);
    const std::vector<T> theta = mc::ctmc_solve(Q);
    std::vector<T> lam(K, zero);
    T sumL = zero;
    for (std::size_t k = 0; k < K && k < arv.classes(); ++k) {
        const std::vector<T> t = vecmul(theta, arv.Dc[k]);
        for (const T& v : t) lam[k] += v;
        sumL += lam[k];
    }
    std::vector<T> w(K, T(num_traits<T>::from_int(1) / num_traits<T>::from_int((int)K)));
    if (sumL > zero)
        for (std::size_t k = 0; k < K; ++k) w[k] = T(lam[k] / sumL);

    std::size_t ntot = 0;
    for (const PhService<T>& s : svc) ntot += s.S.rows();
    std::vector<T> alpha(ntot, zero);
    Matrix<T> Tm(ntot, ntot, zero);
    std::size_t off = 0;
    for (std::size_t k = 0; k < K; ++k) {
        for (std::size_t i = 0; i < svc[k].S.rows(); ++i) {
            alpha[off + i] = T(w[k] * svc[k].sigma[i]);
            for (std::size_t j = 0; j < svc[k].S.cols(); ++j) Tm(off + i, off + j) = svc[k].S(i, j);
        }
        off += svc[k].S.rows();
    }
    return qsys::ServiceLaw<T>::phase_type(alpha, Tm);
}

/**
 * `cellsum(sn.visits)`: V(i,k), the visits summed over the chains, at STATION
 * level.
 *
 * `sn.visits` is stateful-indexed, so the station index has to be mapped
 * through `stateful_of_station` before the sum -- the two agree only when every
 * stateful node is a station, which is not the case in a fork-join model. Every
 * MAM analyzer needs it, so it lives here rather than being written out again in
 * each; `mna_detail` and `solver_mam_basic_mmap` both call this one.
 */
template <class T>
Matrix<T> station_visits(const qn::NetworkStruct<T>& L) {
    Matrix<T> V(L.nstations, L.nclasses, num_traits<T>::from_int(0));
    for (std::size_t c = 0; c < L.nchains; ++c)
        for (std::size_t i = 0; i < L.nstations; ++i) {
            const std::size_t sf = L.stateful_of_station(i + 1) - 1;
            for (std::size_t k = 0; k < L.nclasses; ++k) V(i, k) += L.visits[c](sf, k);
        }
    return V;
}

/** The reference's terminal `X(isnan(X))=0` sweep over the metric matrices. */
template <class T>
void zero_nans(Matrix<T>& A) {
    for (std::size_t i = 0; i < A.rows(); ++i)
        for (std::size_t j = 0; j < A.cols(); ++j)
            if (std::isnan(num_traits<T>::to_double(A(i, j)))) A(i, j) = num_traits<T>::from_int(0);
}

/** Output of `mam_truncate_renorm`. */
template <class T>
struct TruncRenorm {
    T meanQ;
    T lossProb;
    std::vector<T> p;
};

/**
 * Port of `mam_truncate_renorm`: solve the infinite-buffer MMAP[K]/PH[K]/1
 * FCFS queue, truncate the AGGREGATE marginal at the buffer capacity and
 * renormalize.
 *
 * Multi-class input is first aggregated to a single-class MMAP/PH/1, because
 * `ncDistr` returns the PER-CLASS marginal P(N_k = n) while the truncation
 * needs the joint P(N_total = n).
 */
template <class T>
TruncRenorm<T> truncate_renorm(const Mmap<T>& arv, const std::vector<PhService<T>>& svc,
                               std::size_t capK) {
    const T zero = num_traits<T>::from_int(0);
    Mmap<T> call = arv;
    std::vector<PhService<T>> scall = svc;
    if (svc.size() > 1) {
        const qsys::ServiceLaw<T> mix = svc_mixture(arv, svc);
        Matrix<T> Dsum(arv.order(), arv.order(), zero);
        for (std::size_t k = 0; k < arv.classes(); ++k)
            for (std::size_t i = 0; i < Dsum.rows(); ++i)
                for (std::size_t j = 0; j < Dsum.cols(); ++j) Dsum(i, j) += arv.Dc[k](i, j);
        call.D0 = arv.D0;
        call.D1 = Dsum;
        call.Dc.assign(1, Dsum);
        scall.assign(1, PhService<T>{mix.ph_alpha, mix.ph_T});
    }
    const std::vector<std::vector<T>> d = mmapph1fcfs_ncdistr(call, scall, capK + 1);
    TruncRenorm<T> out;
    out.p.assign(capK + 1, zero);
    T mass = zero;
    for (std::size_t n = 0; n <= capK; ++n) {
        out.p[n] = num_abs(d[0][n]);
        mass += out.p[n];
    }
    if (!(mass > zero)) {
        out.p.assign(capK + 1, zero);
        out.p[0] = num_traits<T>::from_int(1);
    } else {
        for (std::size_t n = 0; n <= capK; ++n) out.p[n] /= mass;
    }
    T m = zero;
    for (std::size_t n = 0; n <= capK; ++n) m += num_traits<T>::from_int((int)n) * out.p[n];
    if (m < zero) m = zero;
    if (num_traits<T>::to_double(m) > static_cast<double>(capK))
        m = num_traits<T>::from_int((int)capK);
    out.meanQ = m;
    out.lossProb = out.p[capK];
    return out;
}

/**
 * One FCFS / HOL station of the decomposition, the inner body of
 * the reference's `case {SchedStrategy.FCFS, HOL, FCFSPRPRIO}` branch. The C++
 * SchedStrategy enumerator for FCFSPRPRIO exists (lang_types.h:144); what is
 * missing is a MAM analyzer for it (BUTools' MMAPPH1PRPR is not ported), so
 * only FCFS and HOL are dispatched here and an FCFSPRPRIO station falls to the
 * station-ladder throw instead.
 *
 * The arrival stream reaching the station is assembled here, chain by chain, in
 * the reference's order: mark the chain stream by the per-class visit shares,
 * retarget the per-class mean inter-arrival times, then superpose across
 * chains. Chain 1 is COLLAPSED to a single mark on the way in -- the reference
 * writes `{aggr{1} aggr{2} aggr{2}}` after superposing with a zero-rate
 * exponential -- so the number of marks reaching the queue solver is
 * 1 + sum_{c>=2} |chain c|, not K. That equals K exactly when chain 1 holds one
 * class; when it does not, MATLAB itself fails (the utilization line divides a
 * 1 x R vector by a 1 x K one, and the queue solve requests K outputs from an
 * R-class analyzer), so the mismatch is refused BY NAME here rather than
 * answered with a silently different model.
 */
/**
 * The setup / delay-off pair of a station, or false when it has none.
 *
 * Presence in `setupparam` IS the reference's `sn.hassetup(ist)`, and the
 * reference reads ONE pair per station (`{end}`, the last class that declares
 * one) rather than a pair per class.
 */
template <class T>
bool station_setup_pair(const qn::NetworkStruct<T>& L, std::size_t ist, lang::Distrib<T>& setup,
                        lang::Distrib<T>& delayoff) {
    typename std::map<std::size_t, qn::SetupDelayOffParam<T>>::const_iterator it =
        L.setupparam.find(ist);
    if (it == L.setupparam.end()) return false;
    return it->second.last(setup, delayoff);
}

/**
 * The open-class setup/delay-off queue, solver_mam_basic.m:507-547.
 *
 * `mu_k` is the per-class service RATE the phase representation already carries
 * scaled by S/nservers, so rho_k is per-class utilization. The aggregate rate
 * is defined by Lambda / sum(rho_k), which makes the surrogate's utilization
 * agree with the station's by construction; the QBD then answers one scalar
 * mean queue length, split back by load share.
 */
template <class T>
void mam_setup_qbd(const qn::NetworkStruct<T>& L, std::size_t ist,
                   const std::vector<T>& aggrLambda,
                   const std::vector<std::vector<PhService<T>>>& svc, const Matrix<T>& S,
                   const std::vector<std::vector<T>>& rates, std::size_t R, double ns,
                   const lang::Distrib<T>& setup, const lang::Distrib<T>& delayoff,
                   std::vector<T>& Qret) {
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const std::size_t i0 = ist - 1, K = L.nclasses;
    for (std::size_t r = 0; r < K; ++r) Qret[r] = zero;
    if constexpr (!num_traits<T>::has_transcendental) {
        throw UnsupportedError(
            "solver_mam_basic: a setup/delay-off station is solved by a QBD whose G matrix needs "
            "a logarithm; rerun with --arith double or real");
    } else {
        const Map<T> su = lang::dist_to_map(setup), doff = lang::dist_to_map(delayoff);
        const T alpharate = map_lambda(su), alphascv = map_scv(su);
        const T betarate = map_lambda(doff), betascv = map_scv(doff);

        std::vector<T> rho(K, zero);
        T rho_total = zero;
        for (std::size_t r = 0; r < K; ++r) {
            if (L.disabled[i0][r]) continue;
            // mean of the already-scaled service phase: sigma (-S)^-1 e
            const Matrix<T>& Sr = svc[i0][r].S;
            if (Sr.rows() == 0) continue;
            Matrix<T> negS = Sr;
            for (std::size_t a = 0; a < negS.rows(); ++a)
                for (std::size_t b = 0; b < negS.cols(); ++b) negS(a, b) = -negS(a, b);
            const Matrix<T> negSinv = inverse(negS);
            T mean = zero;
            for (std::size_t a = 0; a < negSinv.rows() && a < svc[i0][r].sigma.size(); ++a)
                for (std::size_t b = 0; b < negSinv.cols(); ++b)
                    mean += svc[i0][r].sigma[a] * negSinv(a, b);
            if (!(mean > zero)) continue;
            std::size_t c = L.nchains;
            for (std::size_t cc = 0; cc < L.nchains; ++cc)
                if (L.chains[cc][r]) c = cc;
            if (c == L.nchains) continue;
            // rho = lambda / mu, and mu is 1/mean of the already-scaled phase
            rho[r] = T(rates[c][r] * mean);
            rho_total += rho[r];
        }
        if (!(rho_total > zero)) return;

        T lam_total = zero;
        for (std::size_t r = 0; r < K && r < aggrLambda.size(); ++r) lam_total += aggrLambda[r];
        if (R == 1 && !aggrLambda.empty()) lam_total = aggrLambda[0];
        if (!(lam_total > zero)) return;
        const T aggrRate = T(lam_total / rho_total);
        const T qtot =
            mam::qbd_setupdelayoff(lam_total, aggrRate, alpharate, alphascv, betarate, betascv);
        for (std::size_t r = 0; r < K; ++r)
            if (rho[r] > zero) Qret[r] = T(qtot * rho[r] / rho_total);
        (void)S;
        (void)ns;
        (void)one;
    }
}

/**
 * The closed-class setup/delay-off queue, solver_mam_basic.m.
 *
 * THE CLOSED VACATION QUEUE, SOLVED. What stood here was the per-instance
 * cold-start race `R = p_cold*E[setup] + S`: it raced the delay-off against the
 * per-instance idle time and carried NO queueing term, so it described a
 * serverless instance pool rather than a single-server vacation queue and
 * reported the SAME response time across a tenfold change in the setup mean
 * (BUG-78). `qbd_setupdelayoff_closed` solves the finite level-dependent chain
 * the simulator walks, and the aggregate response it returns is split back over
 * the classes by R_k = W + S_k -- the decomposition the finite-capacity branch
 * already uses, and an identity when the chain holds one class.
 */
template <class T>
void mam_setup_closed(const qn::NetworkStruct<T>& L, std::size_t ist,
                      const std::vector<std::vector<PhService<T>>>& svc, const Matrix<T>& S,
                      const std::vector<std::vector<T>>& rates, const std::vector<T>& ztchain,
                      const Matrix<T>& V, double ns, const lang::Distrib<T>& setup,
                      const lang::Distrib<T>& delayoff, const Matrix<T>& QNprev,
                      std::vector<T>& Qret) {
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const std::size_t i0 = ist - 1, K = L.nclasses;
    for (std::size_t r = 0; r < K; ++r) Qret[r] = zero;
    if constexpr (!num_traits<T>::has_transcendental) {
        throw UnsupportedError(
            "solver_mam_basic: the closed setup/delay-off chain fits a phase-type setup and "
            "delay-off, which needs a square root; rerun with --arith double or real");
    } else {
        const Map<T> su = lang::dist_to_map(setup), doff = lang::dist_to_map(delayoff);
        const T alpharate = map_lambda(su), alphascv = map_scv(su);
        const T betarate = map_lambda(doff), betascv = map_scv(doff);
        const T ft = num_traits<T>::from_double(GlobalConstants::FineTol);

        for (std::size_t c = 0; c < L.nchains; ++c) {
            T Nc = zero;
            bool finiteNc = true, anyActive = false;
            for (std::size_t r = 0; r < K; ++r) {
                if (!L.chains[c][r]) continue;
                if (!std::isfinite(L.classes[r].population)) { finiteNc = false; break; }
                Nc += num_traits<T>::from_double(L.classes[r].population);
                if (rates[c][r] > zero && std::isfinite(num_traits<T>::to_double(S(i0, r))))
                    anyActive = true;
            }
            if (!finiteNc || !anyActive || !(Nc > zero)) continue;

            T vtot = zero;
            for (std::size_t j = 0; j < K; ++j)
                if (L.chains[c][j]) vtot += V(i0, j);
            const T ZT = T(ztchain[c] / (vtot > ft ? vtot : ft));

            // THE COMPLEMENTARY DELAY, not the think demand alone.
            // lambda(n) = (Nc-n)/Z is exact only when everything away from this
            // station is a pure delay; with other queues in the network the think
            // demand OVERSTATES the arrival rate and saturates the station. Z is
            // the mean time a customer currently spends away,
            // (Nc - QN_here)/lambda_here at this iterate, floored at ZT so it can
            // never be shorter than the think time it contains. On a Delay+Queue
            // the two coincide at convergence.
            T lamHere = zero, qnHere = zero, tnS = zero, tnTot = zero, svcAny = zero;
            std::size_t svcAnyCount = 0;
            for (std::size_t r = 0; r < K; ++r) {
                if (!L.chains[c][r]) continue;
                T lamK = rates[c][r];
                if (!std::isfinite(num_traits<T>::to_double(lamK))) lamK = zero;
                T sk = S(i0, r);
                if (!std::isfinite(num_traits<T>::to_double(sk))) sk = zero;
                lamHere += lamK;
                if (i0 < QNprev.rows() && r < QNprev.cols() &&
                    std::isfinite(num_traits<T>::to_double(QNprev(i0, r))))
                    qnHere += QNprev(i0, r);
                tnTot += lamK;
                tnS += T(lamK * sk);
                if (sk > zero) { svcAny += sk; ++svcAnyCount; }
            }
            T Zc = ZT;
            if (lamHere > ft && T(Nc - qnHere) > zero) {
                const T zalt = T(T(Nc - qnHere) / lamHere);
                Zc = zalt > ZT ? zalt : ZT;
            }
            // One server serves the whole chain, so the vacation cycle is a
            // property of the STATION: the chain is solved on the aggregate.
            T Sbar = tnTot > ft ? T(tnS / tnTot)
                                : (svcAnyCount > 0
                                       ? T(svcAny / num_traits<T>::from_int(
                                                        static_cast<long>(svcAnyCount)))
                                       : zero);
            if (!(Sbar > ft)) continue;
            const mam::SetupDelayoffClosed<T> cr = mam::qbd_setupdelayoff_closed(
                Nc, Zc, T(one / Sbar), alpharate, alphascv, betarate, betascv);
            if (!std::isfinite(num_traits<T>::to_double(cr.QN)) || !(cr.XN > zero)) continue;
            T Wq = T(T(cr.QN / cr.XN) - Sbar);
            if (!(Wq > zero)) Wq = zero;
            for (std::size_t r = 0; r < K; ++r) {
                if (!L.chains[c][r]) continue;
                if (L.disabled[i0][r]) continue;
                // S/c, NOT S: the caller adds the surrogate-delay jobs
                // TN*S*(c-1)/c back, so a full S here counts the service term
                // (2c-1)/c times. The two together make S. At c=1 the division
                // is the identity. Wq still comes from a SINGLE-SERVER chain, so
                // a closed multiserver setup station is approximated, not solved.
                const T cserv = (std::isfinite(ns) && ns >= 1.0)
                                    ? num_traits<T>::from_double(ns)
                                    : one;
                if (rates[c][r] > zero && std::isfinite(num_traits<T>::to_double(S(i0, r))))
                    Qret[r] = T(rates[c][r] * T(Wq + T(S(i0, r) / cserv)));
            }
        }
        (void)svc;
    }
}

template <class T>
void solve_fcfs_station(const qn::NetworkStruct<T>& L, const MamOptions& opt, std::size_t ist,
                        const std::vector<T>& lambda, const Matrix<T>& V, const Matrix<T>& S,
                        const std::vector<std::vector<bool>>& Sknown,
                        const std::vector<std::vector<Map<T>>>& PH,
                        const std::vector<std::vector<bool>>& PHset,
                        const std::vector<std::vector<PhService<T>>>& svc,
                        const std::vector<Mmap<T>>& chainSysArrivals,
                        const std::vector<T>& ztchain, Matrix<T>& QN, Matrix<T>& UN,
                        Matrix<T>& RN, Matrix<T>& TN, std::vector<bool>& exact_station) {
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const std::size_t M = L.nstations, K = L.nclasses, C = L.nchains;
    const double ns = L.stations[ist - 1].nservers;
    const std::size_t i0 = ist - 1;

    // The reference also routes FCFSPRPRIO here, through BUTools' MMAPPH1PRPR,
    // which is not ported; this function is only ever called for FCFS or HOL
    // (see the dispatch below), so that arm is unreachable, not because the
    // C++ SchedStrategy lacks the enumerator (it does not, lang_types.h:144).
    if (L.stations[i0].sched == SchedStrategy::HOL && K > 1) {
        bool distinct = false;
        for (std::size_t r = 1; r < K; ++r)
            if (L.classes[r].prio != L.classes[0].prio) distinct = true;
        if (distinct)
            throw UnsupportedError(
                std::string("solver_mam_basic: station '") + L.stations[i0].name +
                "' has non-identical class priorities under " +
                lang::sched_to_text(L.stations[i0].sched) +
                ", which the reference solves with BUTools' MMAPPH1NPPR; the priority analyzer "
                "is not ported to C++");
    }

    // rates(ist,:) per chain: the visit rate of every class at this station.
    std::vector<std::vector<T>> rates(C, std::vector<T>(K, zero));
    for (std::size_t c = 0; c < C; ++c)
        for (std::size_t r = 0; r < K; ++r) rates[c][r] = T(V(i0, r) * lambda[c]);

    Mmap<T> aggr;
    for (std::size_t c = 0; c < C; ++c) {
        const std::vector<std::size_t>& ic = L.inchain[c];
        T tot = zero;
        for (std::size_t r : ic) tot += rates[c][r - 1];
        Matrix<T> markProb(1, ic.size(), zero);
        for (std::size_t j = 0; j < ic.size(); ++j)
            markProb(0, j) = (tot > zero) ? T(rates[c][ic[j] - 1] / tot) : zero;
        Mmap<T> cur = mmap_mark_probs(chainSysArrivals[c], markProb);
        // Every arrival process the C++ lang layer can express is Markovian
        // (there is no ME or RAP ProcessType), so mam_chain_arrival_is_markovian
        // is unconditionally true and the normalization always applies.
        cur = mmap_normalize(cur);
        bool anypos = false;
        for (std::size_t j = 0; j < ic.size(); ++j)
            if (rates[c][ic[j] - 1] > zero) anypos = true;
        if (anypos) {
            std::vector<T> tgt(ic.size(), zero);
            for (std::size_t j = 0; j < ic.size(); ++j)
                tgt[j] = (rates[c][ic[j] - 1] > zero)
                             ? T(one / rates[c][ic[j] - 1])
                             : num_traits<T>::from_double(1.0 / GlobalConstants::Zero);
            cur = mmap_scale_perclass(cur, tgt);
        }
        if (c == 0) {
            // Chain 1 keeps its own per-class marks; lumping them into one gave
            // the aggregate the wrong mark count. See _kb/06-solver-catalog.md.
            aggr = cur;
        } else {
            aggr = mmap_super_safe(std::vector<Mmap<T>>{aggr, cur}, opt.space_max);
        }
    }

    // The marks come out chain by chain and every reader below indexes them by
    // CLASS, so permute them into class order.
    {
        std::vector<std::size_t> markorder;
        for (std::size_t c = 0; c < C; ++c)
            for (std::size_t r : L.inchain[c]) markorder.push_back(r);
        if (markorder.size() == aggr.Dc.size() &&
            !std::is_sorted(markorder.begin(), markorder.end())) {
            std::vector<std::size_t> perm(markorder.size());
            for (std::size_t j = 0; j < perm.size(); ++j) perm[j] = j;
            std::stable_sort(perm.begin(), perm.end(),
                             [&markorder](std::size_t a, std::size_t b) {
                                 return markorder[a] < markorder[b];
                             });
            std::vector<Matrix<T>> reordered(perm.size());
            for (std::size_t j = 0; j < perm.size(); ++j) reordered[j] = aggr.Dc[perm[j]];
            aggr.Dc.swap(reordered);
        }
    }

    const std::size_t R = aggr.classes();
    lang::Distrib<T> setup_dist, delayoff_dist;
    const bool has_setup = station_setup_pair(L, ist, setup_dist, delayoff_dist);

    /**
     * WHO ACTUALLY NEEDS ONE MARK PER CLASS, which is not everyone.
     *
     * The loop above superposes every chain, each contributing one mark per
     * class, and permutes the result into class order, so R is K on any model
     * whose chains partition the classes. It used to rewrite chain 1 as
     * `{aggr{1} aggr{2} aggr{2}}` -- ONE aggregate mark, whatever that chain
     * held -- which made R equal 1 + sum_{c>1}|inchain_c|, and the reference
     * carried the same defect until 2026-08-16.
     *
     * The setup/delay-off branches do not read the marking at all: the closed
     * one (solver_mam_basic.m:596-650) works from `rates`, `S` and the think
     * demand, and the open one takes the SCALAR aggregate rate. That is the
     * shape SolverLN sends here for a setup layer -- a single chain over the
     * client, entry and call classes -- and MATLAB solves it, so refusing it
     * was a WRONG REFUSAL and is what made lqn_setup unsolvable in this port.
     *
     * MMAPPH1FCFS and the finite-capacity analyzers do read it: each needs the
     * r-th mark to be the arrival stream of the r-th service law. There the
     * requirement stands, and it is enforced at those branches instead of here,
     * so a shape one branch cannot use no longer refuses the ones that can.
     */
    const std::string mark_refusal =
        std::string("solver_mam_basic: the assembled arrival stream at station '") +
        L.stations[i0].name + "' carries " + std::to_string(R) +
        " marked classes against the model's " + std::to_string(K) +
        "; the reference collapses chain 1 to a single mark, so this analyzer has no arrival "
        "stream to pair with each service law. A class-switching chain at an FCFS station is "
        "refused rather than solved with a mismatched marking";

    const std::vector<T> aggrLambda = mmap_lambda(aggr);
    double aggrUtil = 0.0;
    for (std::size_t r = 0; r < K; ++r) {
        // A class this station never serves has a NaN rate in the reference and
        // its term is dropped by 'omitnan'; the flag replaces the sentinel here.
        if (L.disabled[i0][r]) continue;
        const double lam = num_traits<T>::to_double(aggrLambda[R == 1 ? 0 : r]);
        const double mu = num_traits<T>::to_double(L.rates(i0, r));
        const double den = GlobalConstants::FineTol + mu * ns;
        if (den > 0.0 && std::isfinite(den)) aggrUtil += lam / den;
    }

    std::vector<T> Qret(K, zero);
    bool exact_here = false;   // the reference's mapdcUsed
    T exactRespT = zero;
    bool finiteCapUsed = false;
    T finiteCapMeanQ = zero, finiteCapLossProb = zero;
    std::vector<T> finiteCapLossPerClass;

    bool anyopen = false;
    for (std::size_t r = 0; r < K; ++r)
        if (std::isinf(L.classes[r].population)) anyopen = true;

    if (!(aggrUtil < 1.0 - GlobalConstants::FineTol)) {
        for (std::size_t r = 0; r < K; ++r)
            Qret[r] = num_traits<T>::from_double(L.classes[r].population);
    } else if (anyopen) {
        bool isMapDc = (K == 1) && (L.service[i0][0].type == ProcessType::DET);
        bool isDMc = false;
        std::size_t dmcSource = 0;
        if (K == 1 && !isMapDc && L.service[i0][0].type == ProcessType::EXP) {
            for (std::size_t j = 1; j <= M; ++j)
                if (j != ist && L.service[j - 1][0].type == ProcessType::DET) {
                    isDMc = true;
                    dmcSource = j;
                    break;
                }
        }
        bool isPhM1 = false;
        std::size_t phSource = 0;
        if (K == 1 && !isMapDc && !isDMc && L.service[i0][0].type == ProcessType::EXP &&
            std::isfinite(ns) && ns >= 1.0) {
            for (std::size_t j = 1; j <= M; ++j) {
                if (j == ist) continue;
                const ProcessType sp = L.service[j - 1][0].type;
                if (sp != ProcessType::EXP && sp != ProcessType::DET &&
                    sp != ProcessType::IMMEDIATE && sp != ProcessType::DISABLED) {
                    // qsys_phmc reads the arrival through its (pie, D0) marginal
                    // alone, so the gate is RENEWAL and not merely
                    // non-exponential: a correlated MAP would be answered as its
                    // renewal marginal.
                    if (!is_renewal_map(lang::dist_to_map(L.service[j - 1][0]))) continue;
                    isPhM1 = true;
                    phSource = j;
                    break;
                } else if (sp == ProcessType::EXP && ns > 1.0 && M == 2 &&
                           L.nodes[L.node_of_station(j) - 1].nodetype == qn::NodeType::Source) {
                    // M/M/c: the single-fast-server surrogate is inexact for
                    // c > 1, so take the exact Erlang-C route instead.
                    isPhM1 = true;
                    phSource = j;
                    break;
                }
            }
        }
        // An infinite-buffer closed form would silently discard sn.cap and
        // report a lossless queue, so a finite buffer takes precedence. The
        // gate is a buffer that can BIND, not a finite sn.cap: refreshCapacity
        // derives one from the chain population for every closed model.
        const bool isFiniteCap = std::isfinite(sn::sn_get_buffer_size(L, i0 + 1));
        if (isFiniteCap) {
            isPhM1 = false;
            isDMc = false;
        }
        // EACH CLOSED FORM IS A SURROGATE THAT MAY NOT APPLY, and the reference
        // says so by WRAPPING EACH CALL IN `try ... catch <flag> = false`
        // (solver_mam_basic.m:417-445): the gates above are structural (which
        // laws are present), while the surrogate's own stability is only known
        // once it is evaluated. The PH/M/c arm takes ANOTHER station's service
        // law as the arrival process, so its rho is the ratio of two service
        // rates and can exceed one on a perfectly stable network -- `oqn_basic`
        // (Source Exp(0.1) -> Delay HyperExp -> Queue1 Exp(1)) is exactly that,
        // and this port propagated `qsys_phmc: load rho must be strictly less
        // than 1` out of SolverMAM instead of answering with the generic
        // MMAPPH1FCFS path the reference falls back to. The catch is the GATE,
        // not a mask: a refused surrogate leaves `exact_here` false and the
        // model is still solved, by the analyzer the reference would use.
        if (isPhM1) {
            try {
                const T muQ = T(one / S(i0, 0));
                const Map<T> src = lang::dist_to_map(L.service[phSource - 1][0]);
                const qsys::PhMcResult<T> r = qsys::qsys_phmc(
                    map_pie(src), src.D0, muQ, static_cast<unsigned>(std::llround(ns)));
                Qret[0] = r.meanQueueLength;
                exactRespT = r.meanSojournTime;
                exact_here = true;
            } catch (const std::exception&) {
                isPhM1 = false;
            }
        }
        if (!exact_here && isDMc) {
            try {
                const T muQ = T(one / S(i0, 0));
                const qsys::DmcResult<T> r = qsys::qsys_dmc(
                    L.rates(dmcSource - 1, 0), muQ, static_cast<unsigned>(std::llround(ns)));
                Qret[0] = r.meanQueueLength;
                exactRespT = r.meanSojournTime;
                exact_here = true;
            } catch (const std::exception&) {
                isDMc = false;
            }
        }
        if (!exact_here && !isFiniteCap && isMapDc) {
            try {
                Map<T> arv;
                arv.D0 = aggr.D0;
                arv.D1 = aggr.D1;
                const qsys::MapDcResult<T> r =
                    qsys::qsys_mapdc(arv, S(i0, 0), static_cast<unsigned>(std::llround(ns)));
                Qret[0] = r.meanQueueLength;
                exactRespT = r.meanSojournTime;
                exact_here = true;
            } catch (const std::exception&) {
                isMapDc = false;
            }
        }
        // MAP/M/c: the PH/M/c gate above refused this station because its
        // aggregate arrival stream is NOT renewal. The single-fast-server
        // surrogate the generic path would use ignores the arrival
        // correlation; the level-dependent QBD is exact. c = 1 already goes to
        // the exact MAP/MAP/1 path, so only c > 1 is claimed here.
        bool isMapMc = !exact_here && !isFiniteCap && (K == 1) && !isMapDc && !isDMc &&
                       !isPhM1 && L.service[i0][0].type == ProcessType::EXP &&
                       std::isfinite(ns) && ns > 1.0;
        if (isMapMc) {
            try {
                Map<T> arv;
                arv.D0 = aggr.D0;
                arv.D1 = aggr.D1;
                const T muQ = T(one / S(i0, 0));
                const qsys::MapMcResult<T> r =
                    qsys::qsys_mapmc(arv, muQ, static_cast<unsigned>(std::llround(ns)));
                Qret[0] = r.meanQueueLength;
                exactRespT = r.meanSojournTime;
                exact_here = true;
            } catch (const std::exception&) {
                isMapMc = false;
            }
        }
        // Exact MAP/PH/c. The arms above cover c > 1 only for EXPONENTIAL
        // service; with a phase-type service law the generic path scales the
        // service by nservers and adds a surrogate delay, an approximation that
        // discards the service SHAPE. The service must be RENEWAL, since each
        // freed server restarts at alpha. DET goes to MAP/D/c above, and ME/RAP
        // have no phase-type configuration space at all.
        bool isMapPhc = false;
        if (!exact_here && !isFiniteCap && K == 1 && std::isfinite(ns) && ns > 1.0) {
            const ProcessType st = L.service[i0][0].type;
            if (st != ProcessType::EXP && st != ProcessType::DET && st != ProcessType::ME &&
                st != ProcessType::RAP && st != ProcessType::IMMEDIATE &&
                st != ProcessType::DISABLED) {
                isMapPhc = is_renewal_map(lang::dist_to_map(L.service[i0][0]));
            }
        }
        if (isMapPhc) {
            try {
                Map<T> arv;
                arv.D0 = aggr.D0;
                arv.D1 = aggr.D1;
                const Map<T> svc0 = lang::dist_to_map(L.service[i0][0]);
                const qsys::MapPhcResult<T> r = qsys::qsys_mapphc(
                    arv, map_pie(svc0), svc0.D0, static_cast<unsigned>(std::llround(ns)),
                    static_cast<std::size_t>(500), static_cast<std::size_t>(1), std::vector<T>());
                Qret[0] = r.meanQueueLength;
                exactRespT = r.meanSojournTime;
                exact_here = true;
            } catch (const std::exception&) {
                isMapPhc = false;
            }
        }
        if (exact_here) {
            // one of the closed forms above answered the station
        } else if (isFiniteCap) {
            // The buffer analyzers split the loss by class, so each mark has to
            // be one class's stream; M/M/c/K is exempt, it uses the total rate.
            if (R != K && !(aggr.order() == 1)) throw UnsupportedError(mark_refusal);
            const std::size_t capK = static_cast<std::size_t>(std::llround(L.cap[i0]));
            // mam_detect_mmck: single-phase arrivals and one shared exponential
            // service rate across every active class.
            bool isMmck = (aggr.order() == 1);
            T muMmck = zero;
            if (isMmck) {
                bool any = false;
                double lo = 0.0, hi = 0.0;
                for (std::size_t r = 0; r < K; ++r) {
                    if (L.disabled[i0][r]) continue;
                    if (L.service[i0][r].type != ProcessType::EXP) {
                        isMmck = false;
                        break;
                    }
                    const double v = num_traits<T>::to_double(L.rates(i0, r));
                    if (!(v > 0.0)) continue;
                    if (!any) {
                        lo = hi = v;
                        any = true;
                        muMmck = L.rates(i0, r);
                    } else {
                        lo = std::min(lo, v);
                        hi = std::max(hi, v);
                    }
                }
                if (!any) isMmck = false;
                if (isMmck && hi - lo > 1e-9 * std::max(1.0, hi)) isMmck = false;
            }
            if (isMmck) {
                T lamTot = zero;
                for (const T& v : aggrLambda) lamTot += v;
                const qsys::MmckResult<T> r =
                    qsys::qsys_mmck(lamTot, muMmck, static_cast<unsigned>(std::llround(ns)),
                                    static_cast<unsigned>(capK));
                finiteCapMeanQ = r.meanQueueLength;
                finiteCapLossProb = r.lossProbability;
            } else if (ns == 1.0) {
                // Exact MMAP[K]/G/1/K: the embedded chain at departure epochs
                // resolves the buffer level jointly with the arrival phase, so
                // each class gets its own loss ratio.
                std::vector<PhService<T>> sl;
                for (std::size_t r = 0; r < K; ++r) sl.push_back(svc[i0][r]);
                const qsys::ServiceLaw<T> mix = svc_mixture(aggr, sl);
                const qsys::MmapG1kResult<T> r = qsys::qsys_mmapg1k(aggr.D0, aggr.Dc, mix, capK);
                finiteCapMeanQ = r.meanQueueLength;
                finiteCapLossProb = r.lossAggregate;
                finiteCapLossPerClass = r.lossRatio;
            } else {
                std::vector<PhService<T>> sl;
                for (std::size_t r = 0; r < K; ++r) sl.push_back(svc[i0][r]);
                const TruncRenorm<T> r = truncate_renorm(aggr, sl, capK);
                finiteCapMeanQ = r.meanQ;
                finiteCapLossProb = r.lossProb;
            }
            finiteCapUsed = true;
            exact_station[i0] = true;
        } else if (has_setup) {
            // OPEN SETUP / DELAY-OFF, solver_mam_basic.m:507-547. The station is
            // collapsed to ONE M/G/1 with setup: the aggregate arrival rate, and
            // an aggregate service rate chosen so the aggregate utilization is
            // exactly the sum of the per-class ones. The queue that QBD returns
            // is then split back over the classes by their share of the load.
            mam_setup_qbd(L, ist, aggrLambda, svc, S, rates, R, ns, setup_dist, delayoff_dist,
                          Qret);
        } else {
            // MMAPPH1FCFS treats service as a renewal phase type, discarding
            // service autocorrelation. A single-class single-server station with
            // a genuinely correlated service takes the exact MAP/MAP/1 QBD,
            // which carries the service phase across departures.
            const bool corr =
                (K == 1) && (ns == 1.0) && PHset[i0][0] &&
                std::fabs(num_traits<T>::to_double(map_acf(PH[i0][0], std::vector<unsigned>{1})[0])) >
                    GlobalConstants::CoarseTol;
            if (corr) {
                Map<T> arv;
                arv.D0 = aggr.D0;
                arv.D1 = aggr.Dc[0];
                const QbdMapMap1Result<T> r = qbd_mapmap1(arv, PH[i0][0]);
                Qret[0] = r.QN;
            } else {
                // One arrival mark per service law, which is what the analyzer
                // pairs up; a collapsed marking has no such pairing.
                if (R != K) throw UnsupportedError(mark_refusal);
                // MMAP[K]/G[K]/1 whenever a class carries a service law that is
                // NOT phase type. mmapph1fcfs_ncmean below reads the PH FIT,
                // which matches the mean and, above SCV 1, nothing else; He's
                // transform analysis takes the ORIGINAL law, which L.service
                // still holds. The result is a /1, hence the single-server gate.
                bool gk_done = false;
                if (ns == 1.0 && mam_gk1_applicable(L, i0, K)) {
                    try {
                        std::vector<Matrix<T>> MM;
                        MM.push_back(aggr.D0);
                        Matrix<T> D1sum(aggr.D0.rows(), aggr.D0.cols(),
                                        num_traits<T>::from_int(0));
                        for (std::size_t r = 0; r < K; ++r)
                            for (std::size_t a = 0; a < D1sum.rows(); ++a)
                                for (std::size_t b = 0; b < D1sum.cols(); ++b)
                                    D1sum(a, b) += aggr.Dc[r](a, b);
                        MM.push_back(D1sum);
                        for (std::size_t r = 0; r < K; ++r) MM.push_back(aggr.Dc[r]);
                        std::vector<lang::Distrib<T>> laws;
                        for (std::size_t r = 0; r < K; ++r)
                            laws.push_back(mam_declared_law(L.service[i0][r]));
                        const qsys::MmapGk1Result<T> gk = qsys::qsys_mmapgk1(
                            MM, laws, std::vector<T>(), static_cast<std::size_t>(1), 1e-12,
                            static_cast<std::size_t>(10000));
                        for (std::size_t r = 0; r < K; ++r)
                            Qret[r] = gk.lambdas[r] * gk.meanSojournTime[r];
                        // NOT exact_here: that flag carries ONE scalar response
                        // time for the whole station, which the single-class
                        // arms above own. Here the answer is per class, and with
                        // ns == 1 the generic tail computes RN = QN/TN and adds
                        // no surrogate delay, which is exactly right.
                        gk_done = true;
                    } catch (const std::exception&) {
                        gk_done = false;
                    }
                }
                if (!gk_done) {
                    std::vector<PhService<T>> sl;
                    for (std::size_t r = 0; r < R; ++r) sl.push_back(svc[i0][r]);
                    const std::vector<T> m = mmapph1fcfs_ncmean(aggr, sl);
                    for (std::size_t r = 0; r < K; ++r) Qret[r] = m[r];
                }
            }
        }
    } else {
        // Every class closed. The queue-length DISTRIBUTION is truncated at each
        // class's own population, which is what stops the decomposition from
        // reporting more jobs than the chain owns.
        std::size_t maxLevel = 1;
        for (std::size_t r = 0; r < K; ++r)
            if (std::isfinite(L.classes[r].population))
                maxLevel += static_cast<std::size_t>(std::llround(L.classes[r].population));
        // "the station receives no arrivals" is a property of the AGGREGATE
        // stream. The reference passes {D0, Dc[0], ...} to map_lambda, which
        // reads only its second argument, so it tested CLASS 1 alone and a
        // station whose class 1 is disabled fell here however busy the rest was.
        const Map<T> probe = aggr.map();
        if (num_traits<T>::to_double(map_lambda(probe)) < GlobalConstants::FineTol) {
            for (std::size_t r = 0; r < K; ++r)
                Qret[r] = (L.rates(i0, 0) > zero)
                              ? T(num_traits<T>::from_double(GlobalConstants::FineTol) /
                                  L.rates(i0, 0))
                              : zero;
        } else if (has_setup) {
            // CLOSED SETUP / DELAY-OFF, solver_mam_basic.m:596-650. No QBD here:
            // a closed job returns after its own chain's think time, so what
            // matters is the race between that gap and the delay-off timer.
            // p_cold is the delay-off LST at 1/ZT, the probability the server
            // has already shut down when the job comes back, and the class holds
            // (expected setup paid + its own service) jobs by Little.
            mam_setup_closed(L, ist, svc, S, rates, ztchain, V, ns, setup_dist, delayoff_dist,
                             QN, Qret);
        } else {
            // Same pairing requirement as the open branch above.
            if (R != K) throw UnsupportedError(mark_refusal);
            std::vector<PhService<T>> sl;
            for (std::size_t r = 0; r < R; ++r) sl.push_back(svc[i0][r]);
            const std::vector<std::vector<T>> pd = mmapph1fcfs_ncdistr(aggr, sl, maxLevel);
            for (std::size_t r = 0; r < K; ++r) {
                const std::size_t Nk =
                    static_cast<std::size_t>(std::llround(L.classes[r].population));
                std::vector<T> p(Nk + 1, zero);
                const std::vector<T>& src = pd[r];
                T acc = zero;
                for (std::size_t n = 0; n < Nk; ++n) {
                    p[n] = num_abs(src[n]);
                    acc += p[n];
                }
                // Truncating at N(k) folds the tail into that level, so the
                // complement must be taken over the TRUNCATED vector.
                p[Nk] = num_abs(T(one - acc));
                T mass = zero;
                for (std::size_t n = 0; n <= Nk; ++n) mass += p[n];
                T m = zero;
                if (mass > zero)
                    for (std::size_t n = 0; n <= Nk; ++n)
                        m += num_traits<T>::from_int((int)n) * T(p[n] / mass);
                if (m < zero) m = zero;
                if (num_traits<T>::to_double(m) > static_cast<double>(Nk))
                    m = num_traits<T>::from_int((int)Nk);
                Qret[r] = m;
            }
        }
    }

    // ---- write the station's metrics back --------------------------------
    if (finiteCapUsed) {
        // Under FCFS the wait in queue is common to every class, so per class
        // R_k = Wq + S_k with Wq recovered from the aggregate queue length.
        std::vector<T> inflow(K, zero), eff(K, zero);
        for (std::size_t r = 0; r < K; ++r) {
            std::size_t c = 0;
            for (std::size_t cc = 0; cc < C; ++cc)
                if (L.chains[cc][r]) {
                    c = cc;
                    break;
                }
            inflow[r] = rates[c][r];
            const T loss = finiteCapLossPerClass.empty() ? finiteCapLossProb
                                                         : finiteCapLossPerClass[r];
            eff[r] = T(inflow[r] * T(one - loss));
        }
        T sumTN = zero;
        for (std::size_t r = 0; r < K; ++r) sumTN += eff[r];
        T Wq = zero;
        if (sumTN > zero) {
            T sw = zero;
            for (std::size_t r = 0; r < K; ++r)
                if (Sknown[i0][r]) sw += eff[r] * S(i0, r);
            const T w = T(T(finiteCapMeanQ / sumTN) - T(sw / sumTN));
            Wq = (w > zero) ? w : zero;
        }
        for (std::size_t r = 0; r < K; ++r) {
            TN(i0, r) = eff[r];
            UN(i0, r) = T(TN(i0, r) * S(i0, r) / num_traits<T>::from_double(ns));
            if (TN(i0, r) > zero) {
                RN(i0, r) = T(Wq + S(i0, r));
                QN(i0, r) = T(TN(i0, r) * RN(i0, r));
            } else {
                RN(i0, r) = zero;
                QN(i0, r) = zero;
            }
        }
    } else {
        for (std::size_t r = 0; r < K; ++r) {
            std::size_t c = 0;
            for (std::size_t cc = 0; cc < C; ++cc)
                if (L.chains[cc][r]) {
                    c = cc;
                    break;
                }
            TN(i0, r) = rates[c][r];
            UN(i0, r) = T(TN(i0, r) * S(i0, r) / num_traits<T>::from_double(ns));
            QN(i0, r) = Qret[r];
            if (exact_here) {
                RN(i0, r) = exactRespT;
                exact_station[i0] = true;
            } else {
                // Add back the jobs at the surrogate delay server the c-fold
                // service speedup removed.
                if (Sknown[i0][r] && std::isfinite(ns))
                    QN(i0, r) = T(QN(i0, r) + TN(i0, r) * S(i0, r) *
                                                  num_traits<T>::from_double((ns - 1.0) / ns));
                RN(i0, r) = (TN(i0, r) > zero) ? T(QN(i0, r) / TN(i0, r)) : zero;
            }
        }
    }
}

}  // namespace basic_detail

/**
 * Port of `solver_mam_basic.m`.
 *
 * @param L   the refreshed struct, with non-Markovian processes already gated
 *            out by the runner
 * @param opt the MAM options; `space_max` is the arrival superposition budget
 */
template <class T>
mva::MvaSolution<T> solver_mam_basic(const qn::NetworkStruct<T>& L, const MamOptions& opt) {
    if constexpr (!num_traits<T>::has_transcendental) {
        throw UnsupportedError(
            "solver_mam_basic: the matrix-analytic station solves run tolerance-terminated "
            "iterations (the Riccati doubling behind MMAP[K]/PH[K]/1, the QBD cyclic reduction) "
            "and need transcendental arithmetic; rerun this model with --arith double or "
            "--arith real");
    } else {
    using namespace basic_detail;
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const std::size_t I = L.nof_nodes(), M = L.nstations, K = L.nclasses, C = L.nchains;
    const double tol = opt.tol;

    // ---- service times, visits and chain demands -------------------------
    Matrix<T> S(M, K, zero);
    std::vector<std::vector<bool>> Sknown(M, std::vector<bool>(K, false));
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t r = 0; r < K; ++r)
            if (!L.disabled[i][r] && L.rates(i, r) > zero) {
                S(i, r) = T(one / L.rates(i, r));
                Sknown[i][r] = true;
            }
    const Matrix<T> V = station_visits(L);
    const mva::ChainDemands<T> dem = mva::sn_get_demands_chain(L);
    // Think demand of each chain, the sum of its demand at the INFINITE-server
    // stations. The closed setup branch divides it by the chain's visits to get
    // the gap between two visits of one job, which is what the cold-start race
    // is against (solver_mam_basic.m:630).
    std::vector<T> ztchain(C, zero);
    for (std::size_t c = 0; c < C; ++c)
        for (std::size_t i = 0; i < M; ++i)
            if (std::isinf(L.stations[i].nservers)) ztchain[c] += dem.Lchain(i, c);

    Matrix<T> QN(M, K, zero), UN(M, K, zero), RN(M, K, zero), TN(M, K, zero);
    std::vector<T> CN(K, zero), XN(K, zero);
    std::vector<bool> exact_station(M, false);  // the reference's mapdcStations

    // ---- per-station service phase-type pairs ----------------------------
    // Det is left alone (preserveDet), so the exact MAP/D/c branch can claim it.
    std::vector<std::vector<Map<T>>> PH(M, std::vector<Map<T>>(K));
    std::vector<std::vector<bool>> PHset(M, std::vector<bool>(K, false));
    std::vector<std::vector<PhService<T>>> svc(M, std::vector<PhService<T>>(K));
    for (std::size_t i = 0; i < M; ++i) {
        const SchedStrategy sc = L.stations[i].sched;
        if (!(sc == SchedStrategy::FCFS || sc == SchedStrategy::HOL ||
              sc == SchedStrategy::PS))
            continue;
        for (std::size_t r = 0; r < K; ++r) {
            if (L.service[i][r].type == ProcessType::DET && opt.preserve_det) continue;
            const double ns = L.stations[i].nservers;
            const T target = std::isfinite(ns) ? T(S(i, r) / num_traits<T>::from_double(ns))
                                               : S(i, r);
            if (L.disabled[i][r] || !Sknown[i][r] || !(target > zero)) {
                // The reference's `any(isnan(D0))` guard: a class this station
                // never serves gets an Immediate service rather than a NaN pair,
                // so it contributes nothing and cannot poison the phase-type
                // checks downstream.
                const T imm = num_traits<T>::from_double(GlobalConstants::Immediate);
                // MATLAB map_exponential takes a MEAN; the rate spelling gives mean 1e-8.
                PH[i][r] = map_exponential_mean(imm);
                PHset[i][r] = true;
                svc[i][r].sigma.assign(1, one);
                svc[i][r].S = Matrix<T>(1, 1, T(-imm));
                continue;
            }
            // ME and RAP are rescaled by the rate ALONE (solver_mam_basic.m:82-88):
            // map_scale finishes with map_normalize, whose clamp on the negative
            // entries is a repair for a MAP and destruction for a matrix
            // exponential, whose negative off-diagonals ARE the representation.
            const ProcessType pt = L.service[i][r].type;
            PH[i][r] = (pt == ProcessType::ME || pt == ProcessType::RAP)
                           ? map_scale_rate(lang::dist_to_map(L.service[i][r]), target)
                           : map_scale(lang::dist_to_map(L.service[i][r]), target);
            PHset[i][r] = true;
            svc[i][r].sigma = map_pie(PH[i][r]);
            svc[i][r].S = PH[i][r].D0;
        }
    }

    // ---- chain classification and the open-chain arrival MMAPs -----------
    std::vector<bool> isopenchain(C, false), isclosedchain(C, false);
    std::vector<T> lambda(C, zero);
    std::vector<Mmap<T>> chainSysArrivals(C);
    for (std::size_t c = 0; c < C; ++c) {
        double tot = 0.0;
        bool inf = false;
        for (std::size_t r : L.inchain[c]) {
            const double n = L.classes[r - 1].population;
            if (std::isinf(n)) inf = true;
            else tot += n;
        }
        (void)tot;
        isopenchain[c] = inf;
        isclosedchain[c] = !inf;
        const std::size_t ist = L.classes[L.inchain[c][0] - 1].refstat;
        if (!inf) continue;
        // Open chain: the arrival stream is the superposition of the per-class
        // source processes, each marked as its own class.
        std::vector<Mmap<T>> parts;
        for (std::size_t r : L.inchain[c]) {
            Mmap<T> m;
            if (L.disabled[ist - 1][r - 1] || !(L.rates(ist - 1, r - 1) > zero)) {
                // No arrivals from this class: the reference substitutes
                // map_exponential(Inf), a zero-rate stream.
                m.D0 = Matrix<T>(1, 1, zero);
                m.D1 = Matrix<T>(1, 1, zero);
                m.Dc.assign(1, Matrix<T>(1, 1, zero));
            } else {
                const Map<T> src = lang::dist_to_map(L.service[ist - 1][r - 1]);
                m.D0 = src.D0;
                m.D1 = src.D1;
                m.Dc.assign(1, src.D1);
                lambda[c] += L.rates(ist - 1, r - 1);
            }
            parts.push_back(m);
        }
        chainSysArrivals[c] = parts[0];
        for (std::size_t p = 1; p < parts.size(); ++p)
            chainSysArrivals[c] = mmap_super_safe(
                std::vector<Mmap<T>>{chainSysArrivals[c], parts[p]}, opt.space_max);
        for (std::size_t r : L.inchain[c])
            TN(ist - 1, r - 1) = L.disabled[ist - 1][r - 1] ? zero : L.rates(ist - 1, r - 1);
    }

    std::vector<bool> finite_srv(M, false);
    for (std::size_t i = 0; i < M; ++i) finite_srv[i] = std::isfinite(L.stations[i].nservers);

    bool ismixed = false;
    {
        bool anyc = false, anyo = false;
        for (std::size_t c = 0; c < C; ++c) {
            anyc = anyc || isclosedchain[c];
            anyo = anyo || isopenchain[c];
        }
        ismixed = anyc && anyo;
    }
    const double Ulim = 1.0 - GlobalConstants::CoarseTol;

    // Floor R at one full service time and restate Q = R*T; see _kb/06-solver-catalog.md (MAM closed-chain population)
    std::vector<std::size_t> all_classes(K);
    for (std::size_t r = 0; r < K; ++r) all_classes[r] = r + 1;
    auto resptime_floor = [&](const std::vector<std::size_t>& rlist) {
        for (std::size_t i = 0; i < M; ++i) {
            if (exact_station[i]) continue;
            for (std::size_t r : rlist) {
                const std::size_t k = r - 1;
                if (V(i, k) > zero) {
                    if (!finite_srv[i]) {
                        RN(i, k) = S(i, k);
                    } else {
                        const T byLittle = (TN(i, k) > zero) ? T(QN(i, k) / TN(i, k)) : zero;
                        RN(i, k) = (S(i, k) > byLittle) ? S(i, k) : byLittle;
                    }
                } else {
                    RN(i, k) = zero;
                }
                QN(i, k) = T(RN(i, k) * TN(i, k));
            }
        }
    };

    // ---- the throughput fixed point --------------------------------------
    Matrix<T> TN_1(M, K, zero);
    double delta = std::numeric_limits<double>::infinity();
    int it = 0;
    while (delta > tol && it <= opt.iter_max) {
        ++it;
        TN_1 = TN;
        double Umax = 0.0;
        for (std::size_t i = 0; i < M; ++i) {
            if (!finite_srv[i]) continue;
            double s = 0.0;
            for (std::size_t r = 0; r < K; ++r) s += num_traits<T>::to_double(UN(i, r));
            if (s > Umax) Umax = s;
        }
        if (ismixed || Umax < 1.0) {
            for (std::size_t c = 0; c < C; ++c) {
                if (!isclosedchain[c]) continue;
                double Nc = 0.0;
                for (std::size_t r : L.inchain[c]) Nc += L.classes[r - 1].population;
                double QNc = 0.0;
                for (std::size_t i = 0; i < M; ++i)
                    for (std::size_t r : L.inchain[c]) QNc += num_traits<T>::to_double(QN(i, r - 1));
                QNc = std::max(tol, QNc);
                T Dsum = zero;
                for (std::size_t i = 0; i < M; ++i) Dsum += dem.Lchain(i, c);
                if (it == 1) {
                    lambda[c] = (Dsum > zero) ? T(num_traits<T>::from_double(Nc) / Dsum) : zero;
                } else {
                    // Iteration-averaged regula falsi: the raw N/QN Newton step
                    // is blended with a no-op, the no-op's weight growing to one
                    // at iter_max.
                    const double w = static_cast<double>(it) / opt.iter_max;
                    const T step = T(num_traits<T>::from_double(Nc / QNc) * lambda[c]);
                    lambda[c] = T(lambda[c] * num_traits<T>::from_double(w) +
                                  step * num_traits<T>::from_double(1.0 - w));
                }
            }
        }
        if (ismixed) {
            double theta = 1.0;
            bool binding = false;
            for (std::size_t i = 0; i < M; ++i) {
                if (!finite_srv[i]) continue;
                double Uopen = 0.0, Uclosed = 0.0;
                for (std::size_t c = 0; c < C; ++c) {
                    const double u =
                        num_traits<T>::to_double(dem.Lchain(i, c)) * num_traits<T>::to_double(lambda[c]);
                    if (isclosedchain[c]) Uclosed += u;
                    else Uopen += u;
                }
                if (Uclosed > tol) {
                    binding = true;
                    theta = std::min(theta, (Ulim - Uopen) / Uclosed);
                }
            }
            if (binding && theta < 1.0)
                for (std::size_t c = 0; c < C; ++c)
                    if (isclosedchain[c])
                        lambda[c] = T(lambda[c] * num_traits<T>::from_double(std::max(0.0, theta)));
        } else if (Umax >= 1.0) {
            for (std::size_t c = 0; c < C; ++c)
                lambda[c] = T(lambda[c] / num_traits<T>::from_double(Umax));
        }

        for (std::size_t c = 0; c < C; ++c) {
            if (isclosedchain[c]) {
                // Poisson surrogate at the current throughput iterate. An open
                // chain keeps its source MMAP, whose rate is re-imposed per
                // station by the per-class scaling below.
                chainSysArrivals[c] =
                    mmap_exponential_vec(std::vector<T>(L.inchain[c].size(), lambda[c]), 1);
            }
            for (std::size_t i = 0; i < M; ++i)
                for (std::size_t r : L.inchain[c]) TN(i, r - 1) = T(V(i, r - 1) * lambda[c]);
        }

        for (std::size_t ind = 1; ind <= I; ++ind) {
            const qn::NodeDef& nd = L.nodes[ind - 1];
            if (nd.station == 0) continue;
            const std::size_t ist = nd.station;
            if (nd.nodetype == qn::NodeType::Join) {
                for (std::size_t c = 0; c < C; ++c)
                    for (std::size_t r : L.inchain[c]) {
                        std::size_t fanin = 0;
                        if (L.rtnodes.rows() > 0)
                            for (std::size_t row = 0; row < L.rtnodes.rows(); ++row)
                                if (L.rtnodes(row, (ind - 1) * K + (r - 1)) != zero) ++fanin;
                        if (fanin == 0) fanin = 1;
                        TN(ist - 1, r - 1) =
                            T(lambda[c] * V(ist - 1, r - 1) / num_traits<T>::from_int((int)fanin));
                        UN(ist - 1, r - 1) = zero;
                        QN(ist - 1, r - 1) = zero;
                        RN(ist - 1, r - 1) = zero;
                    }
                continue;
            }
            const SchedStrategy sched = L.stations[ist - 1].sched;
            const double ns = L.stations[ist - 1].nservers;
            if (sched == SchedStrategy::INF) {
                for (std::size_t c = 0; c < C; ++c)
                    for (std::size_t r : L.inchain[c]) {
                        const std::size_t k = r - 1;
                        if (!(V(ist - 1, k) > zero)) {
                            TN(ist - 1, k) = UN(ist - 1, k) = QN(ist - 1, k) = RN(ist - 1, k) = zero;
                            continue;
                        }
                        TN(ist - 1, k) = T(lambda[c] * V(ist - 1, k));
                        // An infinite server reports U = QLen = T S: dividing by
                        // Inf would annihilate it.
                        UN(ist - 1, k) = T(S(ist - 1, k) * TN(ist - 1, k));
                        QN(ist - 1, k) = T(TN(ist - 1, k) * S(ist - 1, k));
                        RN(ist - 1, k) = (TN(ist - 1, k) > zero)
                                             ? T(QN(ist - 1, k) / TN(ist - 1, k))
                                             : zero;
                    }
            } else if (sched == SchedStrategy::PS) {
                for (std::size_t c = 0; c < C; ++c) {
                    for (std::size_t r : L.inchain[c]) {
                        const std::size_t k = r - 1;
                        if (!(V(ist - 1, k) > zero)) {
                            TN(ist - 1, k) = UN(ist - 1, k) = zero;
                            continue;
                        }
                        TN(ist - 1, k) = T(lambda[c] * V(ist - 1, k));
                        UN(ist - 1, k) = T(S(ist - 1, k) * TN(ist - 1, k) /
                                           num_traits<T>::from_double(ns));
                    }
                    double Uden = 0.0;
                    for (std::size_t k = 0; k < K; ++k) Uden += num_traits<T>::to_double(UN(ist - 1, k));
                    Uden = std::min(1.0 - GlobalConstants::FineTol, Uden);
                    for (std::size_t r : L.inchain[c]) {
                        const std::size_t k = r - 1;
                        if (!(V(ist - 1, k) > zero)) {
                            QN(ist - 1, k) = RN(ist - 1, k) = zero;
                            continue;
                        }
                        QN(ist - 1, k) =
                            T(UN(ist - 1, k) / num_traits<T>::from_double(1.0 - Uden));
                        RN(ist - 1, k) = (TN(ist - 1, k) > zero)
                                             ? T(QN(ist - 1, k) / TN(ist - 1, k))
                                             : zero;
                    }
                }
            } else if (sched == SchedStrategy::FCFS || sched == SchedStrategy::HOL) {
                solve_fcfs_station(L, opt, ist, lambda, V, S, Sknown, PH, PHset, svc,
                                   chainSysArrivals, ztchain, QN, UN, RN, TN, exact_station);
            } else if (sched != SchedStrategy::EXT) {
                // The Source is skipped: it has no queue of its own, and the
                // dispatch overwrites its throughput with sn.rates afterwards.
                // The reference's switch has no default arm, so it skips every
                // other discipline SILENTLY; the featset gate in the runner
                // rejects them first, so reaching here is a defect and is named.
                throw UnsupportedError(
                    std::string("solver_mam_basic: station '") + L.stations[ist - 1].name +
                    "' uses the " + lang::sched_to_text(sched) +
                    " discipline, which the dec.source decomposition does not model (it solves "
                    "INF, PS, FCFS and HOL stations only)");
            }
        }

        // Calibrate the fixed point on the REPORTED QN; see _kb/06-solver-catalog.md (MAM closed-chain population)
        resptime_floor(all_classes);

        delta = 0.0;
        for (std::size_t i = 0; i < M; ++i)
            for (std::size_t r = 0; r < K; ++r)
                delta = std::max(delta, std::fabs(num_traits<T>::to_double(TN(i, r)) -
                                                  num_traits<T>::to_double(TN_1(i, r))));
    }

    // ---- the two rescaling passes ----------------------------------------
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t r = 0; r < K; ++r) QN(i, r) = num_abs(QN(i, r));
    for (int pass = 0; pass < 2; ++pass) {
        for (std::size_t c = 0; c < C; ++c) {
            double Nc = 0.0;
            bool closed = true;
            for (std::size_t r : L.inchain[c]) {
                if (std::isinf(L.classes[r - 1].population)) closed = false;
                else Nc += L.classes[r - 1].population;
            }
            if (closed) {
                T QNc = zero;
                for (std::size_t i = 0; i < M; ++i)
                    for (std::size_t r : L.inchain[c]) QNc += QN(i, r - 1);
                if (QNc > zero) {
                    const T f = T(num_traits<T>::from_double(Nc) / QNc);
                    for (std::size_t i = 0; i < M; ++i)
                        for (std::size_t r : L.inchain[c]) QN(i, r - 1) *= f;
                }
            }
            resptime_floor(L.inchain[c]);
            if (closed && Nc == 0.0) {
                for (std::size_t i = 0; i < M; ++i)
                    for (std::size_t r : L.inchain[c]) {
                        QN(i, r - 1) = UN(i, r - 1) = RN(i, r - 1) = TN(i, r - 1) = zero;
                    }
                for (std::size_t r : L.inchain[c]) CN[r - 1] = XN[r - 1] = zero;
            }
        }
    }

    for (std::size_t r = 0; r < K; ++r) {
        T s = zero;
        for (std::size_t i = 0; i < M; ++i) s += RN(i, r);
        CN[r] = s;
    }
    for (std::size_t c = 0; c < C; ++c)
        for (std::size_t r : L.inchain[c]) XN[r - 1] = TN(L.classes[r - 1].refstat - 1, r - 1);

    mva::MvaSolution<T> out;
    out.Q = QN;
    out.U = UN;
    out.R = RN;
    out.Tp = TN;
    out.C = CN;
    out.X = XN;
    out.method = opt.method;
    out.iter = it + 2;
    return out;
    }  // if constexpr has_transcendental
}

}  // namespace mam
}  // namespace line

#endif  // LINE_SOLVERS_MAM_SOLVER_MAM_BASIC_H
