/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_ENV_ENV_GENERATOR_H
#define LINE_SOLVERS_ENV_ENV_GENERATOR_H

/**
 * @file
 * @ingroup line_solvers
 * Port of `@@SolverENV/getGenerator.m`: the infinitesimal generator of the
 * joint (environment, network) chain, assembled from CTMC stage generators.
 *
 * BLOCK LAYOUT, exactly as the reference builds it. Block (e,e) starts as stage
 * e's generator and is Kronecker-summed with the D0 of every outgoing arc
 * e -> h in increasing h. Block (e,h) starts as a reset matrix (identity on the
 * first min(n_e, n_h) states) and, for each h' != e in increasing order, is
 * Kronecker-multiplied by `D1_(e,h) 1 pie_(h,e)` when h' = h and by
 * `1 pie_(h',h)` otherwise. The flattened matrix is closed with
 * `ctmc_makeinfgen`. The column factors of an off-diagonal block follow the
 * reference's loop order, not the target block's; that is the reference's
 * construction and is reproduced as is.
 *
 * A MISSING ARC is `Exp(0)`, as the SolverENV constructor substitutes for
 * `Disabled`: one phase, D0 = D1 = [0], and a NaN `map_pie` read as ones. The
 * stage solvers must be CTMC, which here is not a choice: the stage generator
 * is always `ctmc_get_generator` on the stage struct.
 */

#include <cmath>
#include <cstddef>
#include <vector>

#include "line/api/mam/map_moment.h"
#include "line/lang/qn/environment.h"
#include "line/solvers/ctmc/solver_ctmc.h"
#include "line/solvers/ctmc/solver_ctmc_getters.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace env {

/** One environment transition event, the reference's `Event(STAGE, node, NaN, NaN, [e,h])`. */
struct EnvStageEvent {
    std::size_t node;  ///< 1-based node of stage `from`'s model
    std::size_t from;  ///< 0-based source stage
    std::size_t to;    ///< 0-based target stage
};

/** The outputs of `getGenerator`, in the reference's order. */
template <class T>
struct EnvGenerator {
    Matrix<T> Q;                                    ///< renvInfGen, flattened and closed
    std::vector<Matrix<T>> stage_Q;                 ///< stageInfGen
    std::vector<std::vector<Matrix<T>>> filt;       ///< renvEventFilt[e][h]
    std::vector<ctmc::CtmcGenerator<T>> stage;      ///< stageEventFilt / stageEvents
    std::vector<EnvStageEvent> events;              ///< renvEvents
};

namespace gen_detail {

template <class T>
Matrix<T> kron(const Matrix<T>& A, const Matrix<T>& B) {
    Matrix<T> K(A.rows() * B.rows(), A.cols() * B.cols(), num_traits<T>::from_int(0));
    for (std::size_t i = 0; i < A.rows(); ++i)
        for (std::size_t j = 0; j < A.cols(); ++j) {
            const T a = A(i, j);
            if (a == num_traits<T>::from_int(0)) continue;
            for (std::size_t k = 0; k < B.rows(); ++k)
                for (std::size_t l = 0; l < B.cols(); ++l)
                    K(i * B.rows() + k, j * B.cols() + l) = T(a * B(k, l));
        }
    return K;
}

template <class T>
Matrix<T> eye(std::size_t n) {
    Matrix<T> I(n, n, num_traits<T>::from_int(0));
    for (std::size_t i = 0; i < n; ++i) I(i, i) = num_traits<T>::from_int(1);
    return I;
}

/** `krons(A,B) = kron(A, I) + kron(I, B)`. */
template <class T>
Matrix<T> krons(const Matrix<T>& A, const Matrix<T>& B) {
    Matrix<T> S = kron(A, eye<T>(B.rows()));
    const Matrix<T> R = kron(eye<T>(A.rows()), B);
    for (std::size_t i = 0; i < S.rows(); ++i)
        for (std::size_t j = 0; j < S.cols(); ++j) S(i, j) = T(S(i, j) + R(i, j));
    return S;
}

/** The arc's (D0, D1), or `Exp(0)` when the arc is not declared. */
template <class T>
void arc_process(const Environment<T>& env, std::size_t e, std::size_t h, Matrix<T>& D0,
                 Matrix<T>& D1) {
    const EnvArc<T>& a = env.arc(e, h);
    if (!a.enabled) {
        D0 = Matrix<T>(1, 1, num_traits<T>::from_int(0));
        D1 = D0;
        return;
    }
    if (!a.dist.has_map())
        throw UnsupportedError("SolverENV.getGenerator: the transition " + env.stage(e).name +
                               " -> " + env.stage(h).name +
                               " has no Markovian (D0,D1) representation");
    D0 = a.dist.D0;
    D1 = a.dist.D1;
}

/** `map_pie` of the arc as a row, or ones(1, n) when it is NaN (an `Exp(0)` arc). */
template <class T>
Matrix<T> arc_pie(const Environment<T>& env, std::size_t f, std::size_t h, std::size_t n) {
    Matrix<T> D0, D1;
    arc_process(env, f, h, D0, D1);
    T rate = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < D1.rows(); ++i)
        for (std::size_t j = 0; j < D1.cols(); ++j) rate += D1(i, j);
    Matrix<T> p(1, n, num_traits<T>::from_int(1));
    if (rate == num_traits<T>::from_int(0)) return p;  // pi D1 e = 0: map_pie is NaN
    mam::Map<T> m;
    m.D0 = D0;
    m.D1 = D1;
    const std::vector<T> v = mam::map_pie(m);
    p = Matrix<T>(1, v.size(), num_traits<T>::from_int(0));
    for (std::size_t j = 0; j < v.size(); ++j) p(0, j) = v[j];
    return p;
}

}  // namespace gen_detail

/**
 * Port of `@@SolverENV/getGenerator.m`.
 *
 * @param env the random environment; every stage must be a queueing network
 * @param opt the CTMC options of the stage solvers (`CTMC(model, 'exact', soptions)`)
 */
template <class T>
EnvGenerator<T> env_get_generator(const Environment<T>& env, const ctmc::CtmcOptions& opt) {
    using namespace gen_detail;
    env.reject_lqn_stages("SolverENV.getGenerator",
                          "the joint generator is assembled from CTMC stage generators");
    const std::size_t E = env.nstages();
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    EnvGenerator<T> g;
    std::vector<std::size_t> nstates(E);
    for (std::size_t e = 0; e < E; ++e) {
        g.stage.push_back(ctmc::ctmc_get_generator(env.stage(e).model, opt));
        g.stage_Q.push_back(g.stage.back().Q);
        nstates[e] = g.stage_Q[e].rows();
    }
    // nphases(i,j) for i != j; the diagonal is never read.
    std::vector<std::vector<std::size_t>> nph(E, std::vector<std::size_t>(E, 1));
    for (std::size_t i = 0; i < E; ++i)
        for (std::size_t j = 0; j < E; ++j)
            if (i != j && env.arc(i, j).enabled) nph[i][j] = env.arc(i, j).dist.phases();

    std::vector<std::vector<Matrix<T>>> B(E, std::vector<Matrix<T>>(E));
    for (std::size_t e = 0; e < E; ++e)
        for (std::size_t h = 0; h < E; ++h) {
            if (h == e) {
                B[e][e] = g.stage_Q[e];
                continue;
            }
            B[e][h] = Matrix<T>(nstates[e], nstates[h], zero);
            for (std::size_t i = 0; i < std::min(nstates[e], nstates[h]); ++i) B[e][h](i, i) = one;
        }

    for (std::size_t e = 0; e < E; ++e)
        for (std::size_t h = 0; h < E; ++h) {
            if (h == e) continue;
            Matrix<T> D0, D1;
            arc_process(env, e, h, D0, D1);
            B[e][e] = krons(B[e][e], D0);
            const Matrix<T> pie = arc_pie(env, h, e, nph[h][e]);
            // D1 * ones(nph(e,h),1) * pie
            Matrix<T> arg(D1.rows(), pie.cols(), zero);
            for (std::size_t i = 0; i < D1.rows(); ++i) {
                T s = zero;
                for (std::size_t j = 0; j < D1.cols(); ++j) s += D1(i, j);
                for (std::size_t j = 0; j < pie.cols(); ++j) arg(i, j) = T(s * pie(0, j));
            }
            B[e][h] = kron(B[e][h], arg);
            const std::size_t nodes = env.stage(e).model.nodes.size();
            for (std::size_t i = 1; i <= nodes; ++i) g.events.push_back(EnvStageEvent{i, e, h});
            for (std::size_t f = 0; f < E; ++f) {
                if (f == h || f == e) continue;
                const Matrix<T> pfh = arc_pie(env, f, h, nph[f][h]);
                Matrix<T> ones_pie(nph[e][h], pfh.cols(), zero);
                for (std::size_t i = 0; i < ones_pie.rows(); ++i)
                    for (std::size_t j = 0; j < pfh.cols(); ++j) ones_pie(i, j) = pfh(0, j);
                B[e][f] = kron(B[e][f], ones_pie);
            }
        }

    // cell2mat: every block row must agree on its height, every column on its width.
    std::vector<std::size_t> roff(E + 1, 0), coff(E + 1, 0);
    for (std::size_t e = 0; e < E; ++e) {
        roff[e + 1] = roff[e] + B[e][e].rows();
        coff[e + 1] = coff[e] + B[e][e].cols();
    }
    for (std::size_t e = 0; e < E; ++e)
        for (std::size_t h = 0; h < E; ++h)
            if (B[e][h].rows() != B[e][e].rows() || B[e][h].cols() != B[h][h].cols())
                throw InputError("SolverENV.getGenerator: block (" + std::to_string(e + 1) + "," +
                                 std::to_string(h + 1) +
                                 ") does not conform to the diagonal blocks (cell2mat)");
    auto flatten = [&](const std::vector<std::vector<Matrix<T>>>& C) {
        Matrix<T> M(roff[E], coff[E], zero);
        for (std::size_t e = 0; e < E; ++e)
            for (std::size_t h = 0; h < E; ++h)
                for (std::size_t i = 0; i < C[e][h].rows(); ++i)
                    for (std::size_t j = 0; j < C[e][h].cols(); ++j)
                        M(roff[e] + i, coff[h] + j) = C[e][h](i, j);
        return M;
    };

    // renvEventFilt{e,h}: only block (e,h) survives, and nothing when e == h.
    g.filt.assign(E, std::vector<Matrix<T>>(E));
    for (std::size_t e = 0; e < E; ++e)
        for (std::size_t h = 0; h < E; ++h) {
            std::vector<std::vector<Matrix<T>>> C(E, std::vector<Matrix<T>>(E));
            for (std::size_t a = 0; a < E; ++a)
                for (std::size_t b = 0; b < E; ++b)
                    C[a][b] = (a != b && a == e && b == h)
                                  ? B[a][b]
                                  : Matrix<T>(B[a][b].rows(), B[a][b].cols(), zero);
            g.filt[e][h] = flatten(C);
        }

    g.Q = flatten(B);
    ctmc::make_infgen(g.Q);
    return g;
}

}  // namespace env
}  // namespace line

#endif  // LINE_SOLVERS_ENV_ENV_GENERATOR_H
