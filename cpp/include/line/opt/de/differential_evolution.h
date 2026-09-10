/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_OPT_DE_DIFFERENTIAL_EVOLUTION_H
#define LINE_OPT_DE_DIFFERENTIAL_EVOLUTION_H

/**
 * Port of `matlab/src/opt/+opt/+de/DifferentialEvolution.m`.
 *
 * This is the self-contained line-opt subset of SciPy Differential Evolution:
 * binomial mutation strategies, immediate updating, dithered mutation,
 * Latin-hypercube initialization, no polish, and objective-side constraint
 * penalties.  Together with NumpyRandomState it reproduces the reference
 * trajectory for a fixed seed, including the data-dependent draw count of
 * out-of-bounds repair.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

#include "line/opt/de/numpy_random_state.h"
#include "line/util/error.h"

namespace line {
namespace opt {
namespace de {

struct DifferentialEvolutionResult {
    std::vector<double> x;
    double fun = std::numeric_limits<double>::infinity();
    std::size_t nit = 0;
    std::size_t nfev = 0;
    bool success = false;
    std::vector<std::vector<double>> best_per_generation;
};

class DifferentialEvolution {
public:
    using Objective = std::function<double(const std::vector<double>&)>;
    using Callback = std::function<bool(const std::vector<double>&, std::size_t)>;

    DifferentialEvolution(Objective objective, std::vector<double> low,
                          std::vector<double> high, std::string strategy = "best1bin",
                          std::size_t popsize_multiplier = 15, std::size_t max_iterations = 100,
                          double mutation_low = 0.5, double mutation_high = 1.0,
                          double recombination = 0.7, double tolerance = 0.01,
                          std::uint64_t seed = 0)
        : objective_(std::move(objective)),
          low_(std::move(low)),
          high_(std::move(high)),
          strategy_(std::move(strategy)),
          popsize_multiplier_(popsize_multiplier),
          max_iterations_(max_iterations),
          dither_low_(std::min(mutation_low, mutation_high)),
          dither_high_(std::max(mutation_low, mutation_high)),
          recombination_(recombination),
          tolerance_(tolerance),
          scale_(mutation_low),
          rng_(seed) {
        validate();
        initialize_latin_hypercube();
    }

    void set_callback(Callback callback) { callback_ = std::move(callback); }
    NumpyRandomState& random_state() { return rng_; }
    const NumpyRandomState& random_state() const { return rng_; }

    DifferentialEvolutionResult solve() {
        DifferentialEvolutionResult result;
        bool stopped = false;
        if (any_infinite()) {
            calculate_initial_energies();
            promote_lowest_energy();
        }

        for (std::size_t nit = 1; nit <= max_iterations_; ++nit) {
            next_generation();
            result.nit = nit;
            if (callback_ && callback_(best(), nit)) stopped = true;
            result.best_per_generation.push_back(best());
            if (stopped || converged()) break;
        }

        if (result.nit == max_iterations_ && !converged()) stopped = true;
        result.x = best();
        result.fun = energies_[0];
        result.nfev = evaluations_;
        result.success = !stopped;
        return result;
    }

private:
    void validate() const {
        if (!objective_) throw InputError("DifferentialEvolution: the objective is empty");
        if (low_.empty()) throw InputError("DifferentialEvolution: no parameters were supplied");
        if (low_.size() != high_.size())
            throw InputError("DifferentialEvolution: lower and upper bounds differ in size");
        for (std::size_t j = 0; j < low_.size(); ++j)
            if (!std::isfinite(low_[j]) || !std::isfinite(high_[j]) || low_[j] > high_[j])
                throw InputError("DifferentialEvolution: every bound must be finite and low <= high");
        if (popsize_multiplier_ == 0)
            throw InputError("DifferentialEvolution: popsize multiplier must be positive");
        if (!(recombination_ >= 0.0 && recombination_ <= 1.0))
            throw InputError("DifferentialEvolution: recombination must lie in [0,1]");
        if (tolerance_ < 0.0)
            throw InputError("DifferentialEvolution: tolerance must be non-negative");
    }

    std::vector<double> scale_parameters(const std::vector<double>& trial) const {
        std::vector<double> out(low_.size());
        for (std::size_t j = 0; j < low_.size(); ++j) {
            const double midpoint = 0.5 * (low_[j] + high_[j]);
            out[j] = midpoint + (trial[j] - 0.5) * std::fabs(low_[j] - high_[j]);
        }
        return out;
    }

    void initialize_latin_hypercube() {
        std::size_t equal = 0;
        for (std::size_t j = 0; j < low_.size(); ++j)
            if (low_[j] == high_[j]) ++equal;
        members_ = std::max<std::size_t>(5, popsize_multiplier_ *
                                               std::max<std::size_t>(1, low_.size() - equal));
        const double segment = 1.0 / static_cast<double>(members_);
        const std::vector<double> flat = rng_.uniform(0.0, 1.0, members_ * low_.size());
        std::vector<std::vector<double>> samples(members_, std::vector<double>(low_.size()));
        std::size_t p = 0;
        for (std::size_t i = 0; i < members_; ++i) {
            const double offset = static_cast<double>(i) / static_cast<double>(members_);
            for (std::size_t j = 0; j < low_.size(); ++j)
                samples[i][j] = segment * flat[p++] + offset;
        }
        population_.assign(members_, std::vector<double>(low_.size()));
        for (std::size_t j = 0; j < low_.size(); ++j) {
            const std::vector<std::size_t> order = rng_.permutation(members_);
            for (std::size_t i = 0; i < members_; ++i) population_[i][j] = samples[order[i]][j];
        }
        energies_.assign(members_, std::numeric_limits<double>::infinity());
        random_index_.resize(members_);
        std::iota(random_index_.begin(), random_index_.end(), std::size_t(0));
        evaluations_ = 0;
    }

    std::vector<double> best() const { return scale_parameters(population_[0]); }

    void calculate_initial_energies() {
        for (std::size_t i = 0; i < members_; ++i) {
            energies_[i] = objective_(scale_parameters(population_[i]));
            ++evaluations_;
        }
    }

    void promote_lowest_energy() {
        const auto pos = std::min_element(energies_.begin(), energies_.end());
        const std::size_t best_i = static_cast<std::size_t>(pos - energies_.begin());
        if (best_i != 0) {
            std::swap(energies_[0], energies_[best_i]);
            std::swap(population_[0], population_[best_i]);
        }
    }

    bool any_infinite() const {
        return std::any_of(energies_.begin(), energies_.end(),
                           [](double x) { return std::isinf(x); });
    }

    bool converged() const {
        if (any_infinite()) return false;
        const double mean =
            std::accumulate(energies_.begin(), energies_.end(), 0.0) / energies_.size();
        double variance = 0.0;
        for (double value : energies_) variance += (value - mean) * (value - mean);
        const double stddev = std::sqrt(variance / energies_.size());
        return stddev <= tolerance_ * std::fabs(mean);  // atol is zero in MATLAB
    }

    std::vector<std::size_t> select_samples(std::size_t candidate, std::size_t count) {
        rng_.shuffle(random_index_);
        std::vector<std::size_t> out;
        out.reserve(count);
        for (std::size_t i = 0; i < count + 1 && out.size() < count; ++i)
            if (random_index_[i] != candidate) out.push_back(random_index_[i]);
        return out;
    }

    std::vector<double> donor(std::size_t candidate, const std::vector<std::size_t>& s) const {
        std::vector<double> b(low_.size());
        for (std::size_t j = 0; j < low_.size(); ++j) {
            if (strategy_ == "rand1bin" || strategy_ == "rand1exp")
                b[j] = population_[s[0]][j] + scale_ *
                       (population_[s[1]][j] - population_[s[2]][j]);
            else if (strategy_ == "randtobest1bin" || strategy_ == "randtobest1exp") {
                b[j] = population_[s[0]][j];
                b[j] += scale_ * (population_[0][j] - b[j]);
                b[j] += scale_ * (population_[s[1]][j] - population_[s[2]][j]);
            } else if (strategy_ == "currenttobest1bin" || strategy_ == "currenttobest1exp")
                b[j] = population_[candidate][j] + scale_ *
                       (population_[0][j] - population_[candidate][j] +
                        population_[s[0]][j] - population_[s[1]][j]);
            else if (strategy_ == "best2bin" || strategy_ == "best2exp")
                b[j] = population_[0][j] + scale_ *
                       (population_[s[0]][j] + population_[s[1]][j] -
                        population_[s[2]][j] - population_[s[3]][j]);
            else if (strategy_ == "rand2bin" || strategy_ == "rand2exp")
                b[j] = population_[s[0]][j] + scale_ *
                       (population_[s[1]][j] + population_[s[2]][j] -
                        population_[s[3]][j] - population_[s[4]][j]);
            else
                b[j] = population_[0][j] +
                       scale_ * (population_[s[0]][j] - population_[s[1]][j]);
        }
        return b;
    }

    std::vector<double> mutate(std::size_t candidate) {
        const std::size_t fill = static_cast<std::size_t>(rng_.randint(low_.size()));
        const std::vector<std::size_t> samples = select_samples(candidate, 5);
        const std::vector<double> b = donor(candidate, samples);
        std::vector<double> trial = population_[candidate];
        const std::vector<double> cross = rng_.uniform(0.0, 1.0, low_.size());
        for (std::size_t j = 0; j < low_.size(); ++j)
            if (cross[j] < recombination_ || j == fill) trial[j] = b[j];
        return trial;
    }

    void repair_bounds(std::vector<double>& trial) {
        std::size_t outside = 0;
        for (double value : trial)
            if (value < 0.0 || value > 1.0) ++outside;
        if (outside == 0) return;
        const std::vector<double> replacements = rng_.uniform(0.0, 1.0, outside);
        std::size_t k = 0;
        for (double& value : trial)
            if (value < 0.0 || value > 1.0) value = replacements[k++];
    }

    void next_generation() {
        scale_ = rng_.uniform(dither_low_, dither_high_);
        for (std::size_t candidate = 0; candidate < members_; ++candidate) {
            std::vector<double> trial = mutate(candidate);
            repair_bounds(trial);
            const double energy = objective_(scale_parameters(trial));
            ++evaluations_;
            if (energy <= energies_[candidate]) {
                population_[candidate] = std::move(trial);
                energies_[candidate] = energy;
                if (energy <= energies_[0]) promote_lowest_energy();
            }
        }
    }

    Objective objective_;
    std::vector<double> low_, high_;
    std::string strategy_;
    std::size_t popsize_multiplier_, max_iterations_;
    double dither_low_, dither_high_, recombination_, tolerance_, scale_;
    Callback callback_;
    NumpyRandomState rng_;
    std::size_t members_ = 0, evaluations_ = 0;
    std::vector<std::vector<double>> population_;
    std::vector<double> energies_;
    std::vector<std::size_t> random_index_;
};

}  // namespace de
}  // namespace opt
}  // namespace line

#endif
