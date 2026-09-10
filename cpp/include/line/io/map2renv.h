/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_IO_MAP2RENV_H
#define LINE_IO_MAP2RENV_H

/**
 * Port of matlab/src/io/map2renv.m and matlab/src/io/MAPQN2RENV.m (python twin
 * in api/io/converters.py): the Markov-modulated image of a network with
 * MAP/MMPP/MMAP arrival or service processes as a queueing network in a
 * random environment.
 *
 * Every non-renewal process is a point process modulated by the CTMC with
 * generator Q = D0 + D1, whose conditional intensity in phase k is
 * lambda(k) = sum_j D1(k,j). The transformation freezes each phase into an
 * environment stage in which the process is the Poisson process of that
 * intensity, and lets the environment switch stages at the rates of Q. With P
 * modulated processes the stage set is the Cartesian product of their phase
 * spaces and the environment generator is the Kronecker sum of the individual
 * Q's, so only one process changes phase at a time. The image is exact in
 * structure for an MMPP (diagonal D1); for a general MAP the phase jumps AT an
 * event epoch are aggregated into Q and their correlation with the event
 * stream is lost. Populations are carried across stage switches unchanged
 * (identity reset), as a phase switch moves no job.
 *
 * ARITHMETIC: field, through `sn_map_modulation` and the struct refresh.
 */

#include <cstddef>
#include <string>
#include <vector>

#include "line/api/sn/sn_map_modulation.h"
#include "line/lang/lang_types.h"
#include "line/lang/qn/environment.h"
#include "line/lang/qn/network_struct.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace io {

/** The INFO output of map2renv: how the stage set was assembled. */
template <class T>
struct Map2RenvInfo {
    std::size_t nstages = 0;
    /** Phase order of each modulated process. */
    std::vector<std::size_t> orders;
    /** True when every process was an MMPP, so the image is exact in structure. */
    bool is_mmpp = true;
    /** The modulation records themselves. */
    std::vector<api::SnMapModulation<T>> mods;
};

/**
 * @param base the network, as its NetworkStruct
 * @param info optional INFO output
 * @param max_stages cap on the environment stage count, the reference's
 *        options.config.map_env_maxstages (default 64)
 */
template <class T>
env::Environment<T> map2renv(const qn::NetworkStruct<T>& base, Map2RenvInfo<T>* info = nullptr,
                             std::size_t max_stages = 64) {
    const double zero_tol = lang::GlobalConstants::Zero;
    const std::vector<api::SnMapModulation<T>> mods = api::sn_map_modulation(base);
    if (mods.empty())
        throw InputError(
            "map2renv: the model declares no MAP, MMPP2 or MMAP process, so it has no "
            "random-environment image");

    const std::size_t P = mods.size();
    std::vector<std::size_t> orders(P);
    std::size_t nstages = 1;
    for (std::size_t p = 0; p < P; ++p) {
        orders[p] = mods[p].order;
        nstages *= orders[p];
    }
    if (nstages > max_stages)
        throw InputError("map2renv: the random-environment image of this model has " +
                         std::to_string(nstages) +
                         " stages, above the max_stages cap of " + std::to_string(max_stages) +
                         "; reduce the order of the modulating processes or raise the cap");

    // stage s enumerates the phase tuples in column-major order, 0-based
    std::vector<std::vector<std::size_t>> phase_of(nstages, std::vector<std::size_t>(P, 0));
    for (std::size_t s = 0; s < nstages; ++s) {
        std::size_t rem = s;
        for (std::size_t p = 0; p < P; ++p) {
            phase_of[s][p] = rem % orders[p];
            rem /= orders[p];
        }
    }

    env::Environment<T> envModel(base.name + "_renv", nstages);
    for (std::size_t s = 0; s < nstages; ++s) {
        std::string nm = "Phase";
        for (std::size_t p = 0; p < P; ++p) nm += "_" + std::to_string(phase_of[s][p] + 1);

        // copy of the base model in which every modulated process is the
        // exponential process of its phase-conditional intensity
        qn::NetworkStruct<T> stage = base;
        stage.name = base.name + "_" + nm;
        for (std::size_t p = 0; p < P; ++p) {
            for (std::size_t c = 0; c < mods[p].classes.size(); ++c) {
                const std::size_t r = mods[p].classes[c];
                const Matrix<T>& D1 = mods[p].D1[c];
                T rate = num_traits<T>::from_int(0);
                for (std::size_t j = 0; j < D1.cols(); ++j) rate += D1(phase_of[s][p], j);
                if (mods[p].arrival) {
                    // a silent phase (zero intensity) stays an ON/OFF source:
                    // the class must exist in every stage so the rate-averaged
                    // limit averages a zero instead of skipping the station
                    if (!(num_traits<T>::to_double(rate) > zero_tol))
                        rate = num_traits<T>::from_double(zero_tol);
                    stage.set_service(mods[p].ist, r, lang::Distrib<T>::exp_rate(rate));
                } else {
                    if (!(num_traits<T>::to_double(rate) > zero_tol))
                        throw InputError(
                            "map2renv: phase " + std::to_string(phase_of[s][p] + 1) +
                            " of the service process of class " + std::to_string(r) +
                            " at station " + std::to_string(mods[p].ist) +
                            " has zero completion rate: the station never empties while the "
                            "environment sits in that stage, so the stage has no steady state "
                            "and the random-environment image is not defined; model the "
                            "stalled server as a breakdown stage instead");
                    stage.set_service(mods[p].ist, r, lang::Distrib<T>::exp_rate(rate));
                }
            }
        }
        stage.refresh_struct();
        envModel.set_stage(s, nm, "item", stage);
    }

    // Kronecker sum of the phase generators: a transition changes the phase of
    // one process only, at the rate that process assigns to it
    for (std::size_t s = 0; s < nstages; ++s) {
        std::size_t stride = 1;
        for (std::size_t p = 0; p < P; ++p) {
            Matrix<T> Qp = mods[p].D0;
            for (const Matrix<T>& D1 : mods[p].D1)
                for (std::size_t a = 0; a < Qp.rows(); ++a)
                    for (std::size_t b = 0; b < Qp.cols(); ++b) Qp(a, b) += D1(a, b);
            const std::size_t k = phase_of[s][p];
            for (std::size_t l = 0; l < orders[p]; ++l) {
                if (l == k) continue;
                const double q = num_traits<T>::to_double(Qp(k, l));
                if (q > zero_tol) {
                    const std::size_t t = static_cast<std::size_t>(
                        static_cast<std::ptrdiff_t>(s) +
                        (static_cast<std::ptrdiff_t>(l) - static_cast<std::ptrdiff_t>(k)) *
                            static_cast<std::ptrdiff_t>(stride));
                    envModel.add_transition(s, t, lang::Distrib<T>::exp_rate(Qp(k, l)));
                }
            }
            stride *= orders[p];
        }
    }

    envModel.init();

    if (info != nullptr) {
        info->nstages = nstages;
        info->orders = orders;
        info->is_mmpp = true;
        for (const api::SnMapModulation<T>& m : mods)
            if (!m.is_mmpp) info->is_mmpp = false;
        info->mods = mods;
    }
    return envModel;
}

/**
 * Retained name for the transformation now implemented by `map2renv`, which
 * generalizes it from a single MMPP2 service process to any number of MAP,
 * MMPP2 or MMAP processes of arbitrary phase order. New code should call
 * `map2renv` directly.
 */
template <class T>
env::Environment<T> mapqn2renv(const qn::NetworkStruct<T>& base,
                               Map2RenvInfo<T>* info = nullptr, std::size_t max_stages = 64) {
    return map2renv(base, info, max_stages);
}

}  // namespace io
}  // namespace line

#endif  // LINE_IO_MAP2RENV_H
