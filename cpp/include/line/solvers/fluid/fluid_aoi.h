/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_FLUID_FLUID_AOI_H
#define LINE_SOLVERS_FLUID_FLUID_AOI_H

/**
 * Age of Information by Markovian fluid queues: a port of `solver_mfq_aoi.m`
 * (identical to `solver_fluid_aoi.m`) with the gate `aoi_is_aoi.m`, the
 * parameter map `aoi_extract_params.m`, and the two aoi-fluid algorithms
 * `solveBufferless.m` and `solveSingleBuffer.m`.
 *
 * WHAT AGE OF INFORMATION IS. Not a delay of a job but the freshness of the
 * information a monitor holds: at time t, the age is t minus the generation
 * time of the most recent update DELIVERED so far. It grows at unit rate and
 * drops on every delivery, so its mean depends on the whole delivery process
 * and not only on the response time -- a queue that delivers late but often can
 * beat one that delivers fast and rarely. Peak AoI (PAoI) is the value reached
 * just before a drop.
 *
 * WHY A FLUID QUEUE COMPUTES IT. The age is a sawtooth: it climbs with slope +1
 * and resets. That is exactly the sample path of a Markov-modulated fluid queue
 * with drift +1 in every state but one, so the stationary age distribution is
 * the stationary distribution of an MFQ level, and it comes out as a
 * MATRIX-EXPONENTIAL triple (g, A, h): P(AoI > t) = g exp(A t) h. Both
 * algorithms build the modulating chain of the age process, then solve one
 * linear system for g and read the moments off A. The reference for the pair is
 * the aoi-fluid toolbox of Dogan, Akar and Atay (BSD 2-Clause, 2020).
 *
 * THE TWO SYSTEMS COVERED, and no others: a BUFFERLESS PH/PH/1/1 where an
 * arrival meeting a busy server is discarded (p = 0) or preempts it (p = 1),
 * and a SINGLE-BUFFER M/PH/1/2 where a waiting update is kept (r = 0) or
 * REPLACED by a fresher one (r = 1). Anything else -- more capacity, more
 * servers, more classes, a second queue -- is refused by name.
 *
 * WHERE p AND r COME FROM. `aoi_extract_params.m` reads the scheduling policy:
 * FCFS gives no preemption and no replacement, LCFS-PR preempts, and LCFS
 * replaces in the buffered system while behaving non-preemptively in the
 * bufferless one. `FluidOptions::aoi_preemption` overrides both, as
 * `options.config.aoi_preemption` does in the reference.
 *
 * A CAVEAT THE REFERENCE CARRIES AND THIS PORT KEEPS. `aoi_dist2ph` builds the
 * PH pair from the (D0, D1) MAP with alpha proportional to theta .* (D1 e),
 * which is the phase distribution AT A COMPLETION rather than at a start. For a
 * multi-phase Erlang the two differ, so the mean service time the AoI branch
 * uses is not the distribution's mean (Erlang(2) with mean 1 is read as mean
 * 0.5). The standard metrics reported alongside AoI inherit that reading. This
 * is reproduced exactly, because MATLAB is the reference and the AoI numbers
 * would otherwise not match; it is called out here so it is never mistaken for
 * a defect of this port.
 *
 * THE STANDARD METRICS ARE AN M/M/1 APPROXIMATION, as the reference states in
 * so many words: QN = rho/(1 - rho) and RN = 1/(mu - lambda) with mu the
 * reciprocal of the mean service time above. They are NOT the metrics of the
 * finite-capacity system being analyzed, which by construction holds at most
 * one or two jobs; the AoI numbers are the output that means something here.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <string>
#include <vector>

#include "line/api/aoi/aoi_dist2ph.h"
#include "line/api/mam/mmap_lambda.h"
#include "line/lang/qn/network_struct.h"
#include "line/util/eig.h"
#include "line/util/error.h"
#include "line/util/expm.h"
#include "line/util/linalg.h"
#include "line/util/lstsq.h"
#include "line/util/lu.h"
#include "line/util/matrix.h"

namespace line {
namespace fluid {

/** What the AoI gate found, when it matches. */
struct AoiTopology {
    bool ok = false;
    std::string error;
    std::size_t source = 0;  ///< 0-based station index
    std::size_t queue = 0;
    std::size_t cls = 0;      ///< the single open class
    double capacity = 0.0;    ///< 1 = bufferless, 2 = single buffer
    lang::SchedStrategy sched = lang::SchedStrategy::FCFS;
};

/**
 * Port of `aoi_is_aoi.m`.
 *
 * Returns the reason rather than throwing, because the reference uses this as
 * a DISPATCH test inside the `mfq` method: a model that fails it is not an
 * error, it is a model for the ordinary single-queue fluid branch.
 */
template <class T>
AoiTopology aoi_is_aoi(const qn::NetworkStruct<T>& sn) {
    AoiTopology t;
    const std::size_t M = sn.nstations, K = sn.nclasses;
    std::vector<std::size_t> open_cls;
    for (std::size_t r = 0; r < K; ++r)
        if (!std::isfinite(sn.classes[r].population)) open_cls.push_back(r);
    if (open_cls.empty()) {
        t.error = "Not an open model - all classes are closed";
        return t;
    }
    if (open_cls.size() > 1) {
        t.error = "Multiple open classes found - AoI analysis requires single class";
        return t;
    }
    t.cls = open_cls[0];

    std::size_t nsrc = 0, nq = 0, nsink = 0;
    for (const qn::NodeDef& nd : sn.nodes) {
        if (nd.nodetype == qn::NodeType::Source) ++nsrc;
        if (nd.nodetype == qn::NodeType::Sink) ++nsink;
        if (nd.nodetype == qn::NodeType::Queue) ++nq;
    }
    if (nsrc != 1) {
        t.error = nsrc == 0 ? "No source node found" : "Multiple source nodes found";
        return t;
    }
    if (nsink != 1) {
        t.error = nsink == 0 ? "No sink node found" : "Multiple sink nodes found";
        return t;
    }
    if (nq != 1) {
        t.error = nq == 0 ? "No queue node found"
                          : "Multiple queue nodes found - AoI analysis supports single queue only";
        return t;
    }
    bool has_src = false, has_q = false;
    for (std::size_t i = 0; i < M; ++i) {
        if (sn.stations[i].nodetype == qn::NodeType::Source) {
            t.source = i;
            has_src = true;
        } else if (sn.stations[i].nodetype == qn::NodeType::Queue) {
            t.queue = i;
            has_q = true;
        } else if (sn.stations[i].nodetype != qn::NodeType::Sink) {
            t.error = "AoI analysis supports Source -> Queue -> Sink only";
            return t;
        }
    }
    if (!has_src || !has_q) {
        t.error = "The Source and the Queue must both be stations";
        return t;
    }

    if (sn.stations[t.queue].nservers != 1.0) {
        t.error = "Queue is not single-server - AoI analysis requires c = 1";
        return t;
    }
    const double cap = (t.queue < sn.cap.size()) ? sn.cap[t.queue]
                                                 : std::numeric_limits<double>::infinity();
    if (!std::isfinite(cap) || cap > 2.0 || cap < 1.0) {
        t.error =
            "AoI analysis requires queue capacity 1 (bufferless) or 2 (single-buffer)";
        return t;
    }
    t.capacity = cap;

    t.sched = sn.stations[t.queue].sched;
    if (t.sched != lang::SchedStrategy::FCFS && t.sched != lang::SchedStrategy::LCFS &&
        t.sched != lang::SchedStrategy::LCFSPR) {
        t.error = "AoI analysis supports FCFS, LCFS or LCFSPR scheduling only";
        return t;
    }

    // A single buffer is analyzed under Poisson arrivals only.
    if (cap == 2.0 && sn.service[t.source][t.cls].D0.rows() > 1) {
        t.error = "Single-buffer (capacity 2) requires exponential arrivals";
        return t;
    }

    const std::size_t rt_idx = t.queue * K + t.cls;
    if (rt_idx < sn.rt.rows() && num_traits<T>::to_double(sn.rt(rt_idx, rt_idx)) > 0.0) {
        t.error = "Self-loop detected at queue - violates AoI model assumptions";
        return t;
    }
    t.ok = true;
    return t;
}

/** The (tau, T) / (sigma, S) pairs and the preemption probability. */
struct AoiParams {
    Matrix<double> Tarr, Ssvc;
    std::vector<double> tau, sigma;
    double p = 0.0;       ///< preemption (bufferless) or replacement (single buffer)
    double lambda = 0.0;  ///< arrival rate, the single-buffer input
};

/**
 * Port of `aoi_extract_params.m`.
 *
 * @param preempt_override `options.config.aoi_preemption`; negative selects the
 *                         policy-driven default
 * @param sn the refreshed network struct
 * @param top the detected age-of-information topology
 */
template <class T>
AoiParams aoi_extract_params(const qn::NetworkStruct<T>& sn, const AoiTopology& top,
                             double preempt_override) {
    AoiParams par;
    const lang::Distrib<T>& svc = sn.service[top.queue][top.cls];
    const std::size_t ls = svc.D0.rows();
    if (ls == 0) throw InputError("aoi_extract_params: the queue has no service process");
    {
        Matrix<double> D0(ls, ls, 0.0), D1(ls, ls, 0.0);
        for (std::size_t a = 0; a < ls; ++a)
            for (std::size_t b = 0; b < ls; ++b) {
                D0(a, b) = num_traits<T>::to_double(svc.D0(a, b));
                D1(a, b) = num_traits<T>::to_double(svc.D1(a, b));
            }
        const aoi::AoiPh<double> ph = aoi::aoi_dist2ph(D0, D1);
        par.sigma = ph.alpha;
        par.Ssvc = ph.Tmat;
    }

    par.lambda = num_traits<T>::to_double(sn.rates(top.source, top.cls));
    if (top.capacity == 1.0) {
        const lang::Distrib<T>& arr = sn.service[top.source][top.cls];
        const std::size_t la = arr.D0.rows();
        if (la == 0) throw InputError("aoi_extract_params: the source has no arrival process");
        Matrix<double> D0(la, la, 0.0), D1(la, la, 0.0);
        for (std::size_t a = 0; a < la; ++a)
            for (std::size_t b = 0; b < la; ++b) {
                D0(a, b) = num_traits<T>::to_double(arr.D0(a, b));
                D1(a, b) = num_traits<T>::to_double(arr.D1(a, b));
            }
        const aoi::AoiPh<double> ph = aoi::aoi_dist2ph(D0, D1);
        par.tau = ph.alpha;
        par.Tarr = ph.Tmat;
    }

    if (preempt_override >= 0.0) {
        par.p = preempt_override;
    } else if (top.capacity == 1.0) {
        // Bufferless: only LCFS-PR preempts the update in service.
        par.p = (top.sched == lang::SchedStrategy::LCFSPR) ? 1.0 : 0.0;
    } else {
        // Single buffer: both LCFS variants replace the waiting update.
        par.p = (top.sched == lang::SchedStrategy::FCFS) ? 0.0 : 1.0;
    }
    return par;
}

/** A matrix-exponential age law: P(age > t) = g exp(A t) h. */
struct AoiMe {
    std::vector<double> g;
    Matrix<double> A;
    std::vector<double> h;
    double mean = std::numeric_limits<double>::quiet_NaN();
    double var = std::numeric_limits<double>::quiet_NaN();
};

/** Both age laws of one system, with the policy parameter that produced them. */
struct AoiSolution {
    AoiMe aoi, paoi;
    std::string system_type;  ///< "bufferless" or "singlebuffer"
    double preemption = std::numeric_limits<double>::quiet_NaN();
};

namespace aoi_detail {

/** Write B into A at (r0, c0). */
inline void put(Matrix<double>& A, std::size_t r0, std::size_t c0, const Matrix<double>& B) {
    for (std::size_t i = 0; i < B.rows(); ++i)
        for (std::size_t j = 0; j < B.cols(); ++j) A(r0 + i, c0 + j) = B(i, j);
}

/** A row vector as a 1 x n matrix, and a column vector as n x 1. */
inline Matrix<double> row(const std::vector<double>& v) {
    Matrix<double> M(1, v.size(), 0.0);
    for (std::size_t j = 0; j < v.size(); ++j) M(0, j) = v[j];
    return M;
}
inline Matrix<double> col(const std::vector<double>& v) {
    Matrix<double> M(v.size(), 1, 0.0);
    for (std::size_t i = 0; i < v.size(); ++i) M(i, 0) = v[i];
    return M;
}

/** MATLAB `g / A` for a row vector g: the x solving x A = g. */
inline std::vector<double> rdivide(const std::vector<double>& g, const Matrix<double>& A) {
    const std::size_t n = A.rows();
    if (g.size() != n) throw InputError("aoi: row division size mismatch");
    Matrix<double> At(n, n, 0.0);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) At(i, j) = A(j, i);
    return solve(At, g);
}

/** g A^{-k} h, the moment form the reference writes as g/(A^k)*h. */
inline double moment_form(const std::vector<double>& g, const Matrix<double>& A,
                          const std::vector<double>& h, unsigned k) {
    std::vector<double> x = g;
    for (unsigned i = 0; i < k; ++i) x = rdivide(x, A);
    double s = 0.0;
    for (std::size_t i = 0; i < x.size() && i < h.size(); ++i) s += x[i] * h[i];
    return s;
}

/**
 * Line 2 of Algorithm 1: the Householder reflector that maps the drift signs
 * onto the first coordinate, used when the level-zero split is known a priori.
 */
inline Matrix<double> householder_p(std::size_t z) {
    std::vector<double> u1(z, 1.0);
    u1[z - 1] = -1.0;
    double nrm = 0.0;
    for (std::size_t i = 0; i < z; ++i) nrm += u1[i] * u1[i];
    nrm = std::sqrt(nrm);
    std::vector<double> u = u1;
    u[0] -= nrm;
    double uu = 0.0;
    for (std::size_t i = 0; i < z; ++i) uu += u[i] * u[i];
    Matrix<double> P = line::eye<double>(z);
    if (uu > 0.0)
        for (std::size_t i = 0; i < z; ++i)
            for (std::size_t j = 0; j < z; ++j) P(i, j) -= 2.0 * u[i] * u[j] / uu;
    return P;
}

/**
 * Lines 3 and 4 of Algorithm 1, shared by both systems: split QR by P, solve
 * the boundary system for (g, d), and return g.
 *
 * @param QR   Q R^{-1}
 * @param P    the orthogonal splitter
 * @param Rd   the drift matrix R
 * @param Qtil the modified generator Qtilde
 * @param a    the number of negative-drift states
 */
struct AlgOneSplit {
    Matrix<double> A, H;
    std::vector<double> g, d;
};
inline AlgOneSplit alg_one(const Matrix<double>& QR, const Matrix<double>& P,
                           const Matrix<double>& Rd, const Matrix<double>& Qtil, std::size_t a) {
    const std::size_t z = QR.rows();
    const std::size_t b = z - a;
    AlgOneSplit s;
    // Ae = P' QR P, then A is its trailing (z-a) block and H the trailing
    // columns of P, transposed.
    Matrix<double> Pt(z, z, 0.0);
    for (std::size_t i = 0; i < z; ++i)
        for (std::size_t j = 0; j < z; ++j) Pt(i, j) = P(j, i);
    const Matrix<double> Ae = matmul(matmul(Pt, QR), P);
    s.A = Matrix<double>(b, b, 0.0);
    for (std::size_t i = 0; i < b; ++i)
        for (std::size_t j = 0; j < b; ++j) s.A(i, j) = Ae(a + i, a + j);
    s.H = Matrix<double>(b, z, 0.0);
    for (std::size_t i = 0; i < b; ++i)
        for (std::size_t j = 0; j < z; ++j) s.H(i, j) = P(j, a + i);

    // EqnMatrix = [H R, -A^{-1} H e; -Qtilde(b+1:z, :), e], solved transposed.
    const Matrix<double> HR = matmul(s.H, Rd);
    std::vector<double> He(b, 0.0);
    for (std::size_t i = 0; i < b; ++i)
        for (std::size_t j = 0; j < z; ++j) He[i] += s.H(i, j);
    const Matrix<double> Ainv = inverse(s.A);
    std::vector<double> AinvHe(b, 0.0);
    for (std::size_t i = 0; i < b; ++i)
        for (std::size_t j = 0; j < b; ++j) AinvHe[i] += Ainv(i, j) * He[j];

    Matrix<double> Eq(z, z + 1, 0.0);
    for (std::size_t i = 0; i < b; ++i) {
        for (std::size_t j = 0; j < z; ++j) Eq(i, j) = HR(i, j);
        Eq(i, z) = -AinvHe[i];
    }
    for (std::size_t i = 0; i < a; ++i) {
        for (std::size_t j = 0; j < z; ++j) Eq(b + i, j) = -Qtil(b + i, j);
        Eq(b + i, z) = 1.0;
    }
    // The system is EqnMatrix' x = [0 ... 0 1]', overdetermined by one row.
    Matrix<double> EqT(z + 1, z, 0.0);
    for (std::size_t i = 0; i < z; ++i)
        for (std::size_t j = 0; j < z + 1; ++j) EqT(j, i) = Eq(i, j);
    std::vector<double> rhs(z + 1, 0.0);
    rhs[z] = 1.0;
    const LstsqResult<double> sol = lstsq(EqT, rhs);
    s.g.assign(sol.x.begin(), sol.x.begin() + static_cast<std::ptrdiff_t>(b));
    s.d.assign(sol.x.begin() + static_cast<std::ptrdiff_t>(b), sol.x.end());
    return s;
}

/** Normalize g by -g A^{-1} h and read the first two moments. */
inline AoiMe finish_me(const std::vector<double>& g, const Matrix<double>& A,
                       const std::vector<double>& h) {
    AoiMe me;
    me.A = A;
    me.h = h;
    const double nc = -moment_form(g, A, h, 1);
    me.g.assign(g.size(), 0.0);
    for (std::size_t i = 0; i < g.size(); ++i) me.g[i] = g[i] / nc;
    me.mean = moment_form(me.g, A, h, 2);
    me.var = -2.0 * moment_form(me.g, A, h, 3) - me.mean * me.mean;
    return me;
}

}  // namespace aoi_detail

/**
 * Port of `solveBufferless.m`: the AoI and PAoI laws of a PH/PH/1/1 system in
 * which an arrival meeting a busy server preempts it with probability p.
 */
inline AoiSolution aoi_solve_bufferless(const std::vector<double>& tau, const Matrix<double>& Tm,
                                        const std::vector<double>& sigma, const Matrix<double>& Sm,
                                        double p) {
    using namespace aoi_detail;
    const std::size_t k = Tm.cols(), l = Sm.cols();
    if (k == 0 || l == 0) throw InputError("aoi_solve_bufferless: empty representation");
    if (tau.size() != k || sigma.size() != l)
        throw InputError("aoi_solve_bufferless: the initial vectors do not match their generators");
    const std::size_t z = 2 * k * l + k + 1, a = 1, b = z - 1;

    std::vector<double> kappa(k, 0.0), nu(l, 0.0);
    for (std::size_t i = 0; i < k; ++i)
        for (std::size_t j = 0; j < k; ++j) kappa[i] -= Tm(i, j);
    for (std::size_t i = 0; i < l; ++i)
        for (std::size_t j = 0; j < l; ++j) nu[i] -= Sm(i, j);
    const Matrix<double> kap = col(kappa), nuc = col(nu);
    const Matrix<double> taur = row(tau), sigr = row(sigma);
    const Matrix<double> Ik = line::eye<double>(k), Il = line::eye<double>(l);
    const Matrix<double> ones_l(l, 1, 1.0), ones_k(k, 1, 1.0);

    // Q11 and Q33 of Eqn (13): service and arrival running together, with the
    // non-preempting and the preempting restart respectively.
    Matrix<double> Q11 = mam::kron(Ik, Sm);
    {
        const Matrix<double> t2 = mam::kron(Tm, Il);
        const Matrix<double> t3 = mam::kron(mam::kron(kap, taur), Il);
        for (std::size_t i = 0; i < k * l; ++i)
            for (std::size_t j = 0; j < k * l; ++j) Q11(i, j) += t2(i, j) + (1.0 - p) * t3(i, j);
    }
    Matrix<double> Q33 = Q11;
    {
        const Matrix<double> t4 = mam::kron(kap, mam::kron(ones_l, mam::kron(taur, sigr)));
        for (std::size_t i = 0; i < k * l; ++i)
            for (std::size_t j = 0; j < k * l; ++j) Q33(i, j) += p * t4(i, j);
    }

    Matrix<double> Q(z, z, 0.0);
    put(Q, 0, 0, Q11);
    put(Q, 0, k * l, mam::kron(Ik, nuc));
    {
        Matrix<double> lastcol = mam::kron(kap, ones_l);
        for (std::size_t i = 0; i < k * l; ++i) Q(i, z - 1) = p * lastcol(i, 0);
    }
    put(Q, k * l, k * l, Tm);
    put(Q, k * l, k * l + k, mam::kron(kap, mam::kron(taur, sigr)));
    put(Q, k * l + k, k * l + k, Q33);
    {
        const Matrix<double> lastcol = mam::kron(ones_k, nuc);
        for (std::size_t i = 0; i < k * l; ++i) Q(k * l + k + i, z - 1) = lastcol(i, 0);
    }

    Matrix<double> Rd = line::eye<double>(z);
    Rd(z - 1, z - 1) = -1.0;
    Matrix<double> Qtil = Q;
    {
        const Matrix<double> ts = mam::kron(taur, sigr);
        for (std::size_t j = 0; j < k * l; ++j) Qtil(z - 1, j) = ts(0, j);
        Qtil(z - 1, z - 1) = -1.0;
    }
    // Q / R with R diagonal +-1 is Q with its last column negated.
    Matrix<double> QR = Q;
    for (std::size_t i = 0; i < z; ++i) QR(i, z - 1) = -QR(i, z - 1);

    const AlgOneSplit s = alg_one(QR, householder_p(z), Rd, Qtil, a);

    AoiSolution out;
    out.system_type = "bufferless";
    out.preemption = p;
    // AoI lives in phases 2 and 3; PAoI in phase 3, weighted by the exit rates.
    std::vector<double> sel(z, 0.0);
    for (std::size_t i = k * l; i < k * l + k * l + k; ++i) sel[i] = 1.0;
    std::vector<double> h(b, 0.0);
    for (std::size_t i = 0; i < b; ++i)
        for (std::size_t j = 0; j < z; ++j) h[i] += s.H(i, j) * sel[j];
    out.aoi = finish_me(s.g, s.A, h);

    std::vector<double> selp(z, 0.0);
    {
        const Matrix<double> kn = mam::kron(ones_k, nuc);
        for (std::size_t i = 0; i < k * l; ++i) selp[k * l + k + i] = kn(i, 0);
    }
    std::vector<double> hp(b, 0.0);
    for (std::size_t i = 0; i < b; ++i)
        for (std::size_t j = 0; j < z; ++j) hp[i] += s.H(i, j) * selp[j];
    out.paoi = finish_me(s.g, s.A, hp);
    return out;
}

/**
 * Port of `solveSingleBuffer.m`: the AoI and PAoI laws of an M/PH/1/2 system in
 * which a waiting update is replaced by a fresher arrival with probability r.
 *
 * Two fluid queues are solved in sequence. The FIRST gives the waiting-time law
 * of an update that finds the server busy, and needs the SCHUR splitter rather
 * than the Householder one: its drift has two negative states and their
 * invariant subspace is not known a priori. Lemma 1 of the reference then
 * turns that law into a PH pair (beta, B), which drives the SECOND queue --
 * the age process proper -- where the split is again explicit.
 */
inline AoiSolution aoi_solve_singlebuffer(double lambda, const std::vector<double>& sigma,
                                          const Matrix<double>& Sm, double r) {
    using namespace aoi_detail;
    const std::size_t l = Sm.cols();
    if (l == 0) throw InputError("aoi_solve_singlebuffer: empty representation");
    if (sigma.size() != l)
        throw InputError("aoi_solve_singlebuffer: sigma does not match its generator");
    if (!(lambda > 0.0)) throw InputError("aoi_solve_singlebuffer: the arrival rate must be positive");

    std::vector<double> nu(l, 0.0);
    for (std::size_t i = 0; i < l; ++i)
        for (std::size_t j = 0; j < l; ++j) nu[i] -= Sm(i, j);
    const Matrix<double> nuc = col(nu), sigr = row(sigma);

    // ---- the waiting-time fluid queue of Eqn (16) --------------------------
    std::size_t z = l + 2;
    const std::size_t a1 = 2, b1 = l;
    Matrix<double> Q(z, z, 0.0);
    put(Q, 0, 0, Sm);
    for (std::size_t i = 0; i < l; ++i) Q(i, l) = nu[i];
    Q(l, l) = -lambda;
    Q(l, l + 1) = lambda;

    Matrix<double> Rd = line::eye<double>(z);
    Rd(z - 2, z - 2) = -1.0;
    Rd(z - 1, z - 1) = -1.0;

    Matrix<double> Qtil(z, z, 0.0);
    for (std::size_t j = 0; j < l; ++j) {
        Qtil(l, j) = lambda * sigma[j];
        Qtil(l + 1, j) = sigma[j];
    }
    Qtil(l, l) = -lambda;
    Qtil(l + 1, l + 1) = -1.0;

    Matrix<double> QR = Q;
    for (std::size_t i = 0; i < z; ++i) {
        QR(i, z - 2) = -QR(i, z - 2);
        QR(i, z - 1) = -QR(i, z - 1);
    }

    // pik solves pik (Q + e e') = e', the reference's rank-one regularization.
    std::vector<double> pik;
    {
        Matrix<double> Aug(z, z, 0.0);
        for (std::size_t i = 0; i < z; ++i)
            for (std::size_t j = 0; j < z; ++j) Aug(i, j) = Q(i, j) + 1.0;
        pik = rdivide(std::vector<double>(z, 1.0), Aug);
    }
    // A1 shifts the singular generator so its stable subspace is separated.
    Matrix<double> A1 = QR;
    {
        std::vector<double> xR(z, 0.0);
        for (std::size_t i = 0; i < z; ++i)
            for (std::size_t j = 0; j < z; ++j) xR[i] += Rd(i, j);
        double den = 0.0;
        for (std::size_t i = 0; i < z; ++i) den += pik[i] * xR[i];
        for (std::size_t i = 0; i < z; ++i)
            for (std::size_t j = 0; j < z; ++j) A1(i, j) += xR[i] * pik[j] / den;
    }
    // ordschur(..., 'rhp'): the right-half-plane eigenvalues come first.
    Matrix<double> P;
    {
        const RealSchur sc = schur_decomposition(A1);
        std::vector<double> key(z, 0.0);
        for (std::size_t i = 0; i < z; ++i) key[i] = (sc.T(i, i) >= 0.0) ? 1.0 : 0.0;
        P = schur_reorder(sc, key).Z;
    }

    const AlgOneSplit s1 = alg_one(QR, P, Rd, Qtil, a1);
    const double c_0 = s1.d.empty() ? 0.0 : s1.d[0];

    // Lemma 1: the waiting law as a PH pair (beta, B).
    Matrix<double> wait_A = s1.A;
    for (std::size_t i = 0; i < b1; ++i) wait_A(i, i) -= r * lambda;
    std::vector<double> selw(z, 0.0);
    selw[l] = 1.0;
    selw[l + 1] = r;
    std::vector<double> wait_H(b1, 0.0);
    for (std::size_t i = 0; i < b1; ++i)
        for (std::size_t j = 0; j < z; ++j) wait_H[i] += s1.H(i, j) * selw[j];
    std::vector<double> wait_g = s1.g;
    {
        const double n1 = 1.0 / (-moment_form(wait_g, wait_A, wait_H, 1) + c_0);
        for (std::size_t i = 0; i < b1; ++i) wait_g[i] *= n1;
    }
    const std::vector<double> Mdiag = [&]() {
        std::vector<double> y = solve(wait_A, wait_H);
        for (std::size_t i = 0; i < y.size(); ++i) y[i] = -y[i];
        return y;
    }();
    Matrix<double> B(l, l, 0.0);
    for (std::size_t i = 0; i < l; ++i)
        for (std::size_t j = 0; j < l; ++j) B(i, j) = wait_A(i, j) * Mdiag[j] / Mdiag[i];
    std::vector<double> beta(l, 0.0);
    for (std::size_t i = 0; i < l; ++i) beta[i] = wait_g[i] * Mdiag[i];
    double beta_0 = 1.0;
    for (std::size_t i = 0; i < l; ++i) beta_0 -= beta[i];
    std::vector<double> psi(l, 0.0);
    for (std::size_t i = 0; i < l; ++i)
        for (std::size_t j = 0; j < l; ++j) psi[i] -= B(i, j);

    // ---- the age fluid queue of Eqn (19) ----------------------------------
    z = 4 * l + 2;
    const std::size_t a2 = 1, b2 = z - 1;
    Matrix<double> Q2(z, z, 0.0);
    put(Q2, 0, 0, B);
    put(Q2, 0, l, mam::kron(col(psi), sigr));
    for (std::size_t i = 0; i < l; ++i) {
        for (std::size_t j = 0; j < l; ++j) {
            Q2(l + i, l + j) = Sm(i, j) - (i == j ? lambda : 0.0);
            Q2(l + i, 2 * l + j) = (i == j) ? lambda : 0.0;
            Q2(2 * l + i, 2 * l + j) = Sm(i, j);
            Q2(3 * l + 1 + i, 3 * l + 1 + j) = Sm(i, j);
        }
        Q2(l + i, 3 * l) = nu[i];
        Q2(3 * l + 1 + i, z - 1) = nu[i];
    }
    put(Q2, 2 * l, 3 * l + 1, mam::kron(nuc, sigr));
    Q2(3 * l, 3 * l) = -lambda;
    for (std::size_t j = 0; j < l; ++j) Q2(3 * l, 3 * l + 1 + j) = lambda * sigma[j];

    Matrix<double> Rd2 = line::eye<double>(z);
    Rd2(z - 1, z - 1) = -1.0;
    Matrix<double> Qtil2 = Q2;
    for (std::size_t j = 0; j < l; ++j) {
        Qtil2(z - 1, j) = beta[j];
        Qtil2(z - 1, l + j) = beta_0 * sigma[j];
    }
    Qtil2(z - 1, z - 1) = -1.0;
    Matrix<double> QR2 = Q2;
    for (std::size_t i = 0; i < z; ++i) QR2(i, z - 1) = -QR2(i, z - 1);

    const AlgOneSplit s2 = alg_one(QR2, householder_p(z), Rd2, Qtil2, a2);

    AoiSolution out;
    out.system_type = "singlebuffer";
    out.preemption = r;
    std::vector<double> sel(z, 0.0);
    for (std::size_t i = 3 * l; i < 4 * l + 1; ++i) sel[i] = 1.0;
    std::vector<double> h(b2, 0.0);
    for (std::size_t i = 0; i < b2; ++i)
        for (std::size_t j = 0; j < z; ++j) h[i] += s2.H(i, j) * sel[j];
    out.aoi = finish_me(s2.g, s2.A, h);

    std::vector<double> selp(z, 0.0);
    for (std::size_t i = 0; i < l; ++i) selp[3 * l + 1 + i] = nu[i];
    std::vector<double> hp(b2, 0.0);
    for (std::size_t i = 0; i < b2; ++i)
        for (std::size_t j = 0; j < z; ++j) hp[i] += s2.H(i, j) * selp[j];
    out.paoi = finish_me(s2.g, s2.A, hp);
    return out;
}

/**
 * `getCdfAoI`: F(t) = 1 - S(t) with S the survival function of the age law.
 *
 * (g, A, h) IS A DENSITY TRIPLE, not a survival one: `finish_me` normalizes g by
 * `-g A^-1 h` so that the density f(t) = g exp(A t) h integrates to 1 and the
 * mean is `g A^-2 h`. The survival function of such a law is
 * S(t) = -g exp(A t) A^-1 h, so F(t) = 1 + g exp(A t) A^-1 h, which is 0 at
 * t = 0 (`g A^-1 h = -1`) and rises to 1. Subtracting the DENSITY instead --
 * the form all four codebases carried until 2026-07-31 -- gives a curve that
 * falls before it rises and is not a distribution function at all.
 */
inline double aoi_cdf(const AoiMe& me, double t) {
    if (t <= 0.0) return 0.0;
    const Matrix<double> E = expm(me.A, t);
    // gE = g exp(A t), then (gE) A^-1 by the same row solve the moments use.
    std::vector<double> gE(me.g.size(), 0.0);
    for (std::size_t j = 0; j < gE.size(); ++j)
        for (std::size_t i = 0; i < me.g.size(); ++i) gE[j] += me.g[i] * E(i, j);
    const std::vector<double> y = aoi_detail::rdivide(gE, me.A);
    double f = 1.0;
    for (std::size_t i = 0; i < y.size() && i < me.h.size(); ++i) f += y[i] * me.h[i];
    return std::max(0.0, std::min(1.0, f));
}

/** What the AoI branch of `mfq` returns: the age laws and the metrics beside them. */
struct FluidAoiResult {
    AoiSolution age;
    std::vector<double> QN, UN, RN, TN;  ///< per class, at the queue
    double lambda = 0.0, mu = 0.0;
};

/** Port of `solver_mfq_aoi.m`. */
template <class T>
FluidAoiResult fluid_aoi(const qn::NetworkStruct<T>& sn, const AoiTopology& top,
                         double preempt_override) {
    if (!top.ok)
        throw UnsupportedError("fluid aoi: " +
                               (top.error.empty() ? std::string("not an AoI topology") : top.error));
    const AoiParams par = aoi_extract_params(sn, top, preempt_override);

    FluidAoiResult out;
    if (top.capacity == 1.0)
        out.age = aoi_solve_bufferless(par.tau, par.Tarr, par.sigma, par.Ssvc, par.p);
    else
        out.age = aoi_solve_singlebuffer(par.lambda, par.sigma, par.Ssvc, par.p);

    const std::size_t K = sn.nclasses;
    out.QN.assign(K, 0.0);
    out.UN.assign(K, 0.0);
    out.RN.assign(K, 0.0);
    out.TN.assign(K, 0.0);

    // The mean service time of the PH pair the AoI algorithm was given.
    const std::size_t l = par.Ssvc.rows();
    const std::vector<double> y = aoi_detail::rdivide(par.sigma, par.Ssvc);
    double mean_svc = 0.0;
    for (std::size_t i = 0; i < l; ++i) mean_svc -= y[i];
    const double lambda = par.lambda;
    const double mu = 1.0 / mean_svc;
    out.lambda = lambda;
    out.mu = mu;
    const double rho = lambda / mu;
    const std::size_t c = top.cls;
    out.UN[c] = std::min(1.0, rho);
    out.TN[c] = (rho < 1.0) ? lambda : mu;
    if (rho < 1.0) {
        out.QN[c] = rho / (1.0 - rho);
        out.RN[c] = 1.0 / (mu - lambda);
    } else {
        out.QN[c] = std::numeric_limits<double>::infinity();
        out.RN[c] = std::numeric_limits<double>::infinity();
    }
    return out;
}

}  // namespace fluid
}  // namespace line

#endif  // LINE_SOLVERS_FLUID_FLUID_AOI_H
