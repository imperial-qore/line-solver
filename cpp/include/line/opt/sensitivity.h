/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_OPT_SENSITIVITY_H
#define LINE_OPT_SENSITIVITY_H

#include <cmath>
#include <map>
#include <optional>
#include <string>
#include <vector>

#include "line/api/pfqn/pfqn_sens.h"
#include "line/lang/qn/network_builder.h"
#include "line/opt/results.h"
#include "line/solvers/ctmc/solver_ctmc_sens.h"
#include "line/solvers/nc/sn_pf_params.h"

namespace line {
namespace opt {

inline std::optional<SensitivityData> open_sensitivities(qn::Network<double>& model) {
    const qn::NetworkStruct<double>& sn = model.get_struct();
    if (sn.classes.empty()) return std::nullopt;
    for (const qn::JobClass& jobclass : sn.classes)
        if (!std::isinf(jobclass.population)) return std::nullopt;

    const nc::PfParams<double> pf = nc::sn_get_product_form_params(sn);
    SensitivityData sensitivity;
    for (std::size_t i = 0; i < sn.nstations; ++i) {
        const qn::Station<double>& station = sn.stations[i];
        if (station.sched == lang::SchedStrategy::EXT) continue;
        const bool delay = station.sched == lang::SchedStrategy::INF;
        if (!delay && station.nservers > 1.0) return std::nullopt;

        std::vector<double> demand(sn.nclasses, 0.0), rho(sn.nclasses, 0.0);
        for (std::size_t r = 0; r < sn.nclasses; ++r) {
            const double rate = sn.rates(i, r);
            const double visits = pf.V(i, r);
            if (std::isfinite(rate) && rate > 0.0 && visits > 0.0) {
                demand[r] = visits / rate;
                rho[r] = pf.lambda[r] * demand[r];
            }
        }
        double utilization = 0.0;
        if (!delay)
            for (double value : rho) utilization += value;
        const double denominator = 1.0 - utilization;
        if (!delay && denominator <= 0.0) return std::nullopt;

        const std::string& station_name = station.name;
        for (std::size_t s = 0; s < sn.nclasses; ++s) {
            const double rate = sn.rates(i, s);
            if (std::isfinite(rate) && rate > 0.0 && pf.V(i, s) > 0.0)
                sensitivity.add("Util", station_name,
                    SensitivityData::parameter_key(station_name, sn.classes[s].name),
                    -rho[s] / rate);
        }

        for (std::size_t r = 0; r < sn.nclasses; ++r) {
            const double rate_r = sn.rates(i, r);
            if (pf.V(i, r) <= 0.0 || !std::isfinite(rate_r) || rate_r <= 0.0) continue;
            const std::string metric = SensitivityData::metric_key(
                station_name, sn.classes[r].name);
            for (std::size_t s = 0; s < sn.nclasses; ++s) {
                const double rate_s = sn.rates(i, s);
                if (!std::isfinite(rate_s) || rate_s <= 0.0 || pf.V(i, s) <= 0.0)
                    continue;
                const std::string parameter = SensitivityData::parameter_key(
                    station_name, sn.classes[s].name);
                const double d_utilization = delay ? 0.0 : -rho[s] / rate_s;
                const double d_demand = s == r ? -demand[r] / rate_s : 0.0;
                const double d_rho = s == r ? -rho[r] / rate_s : 0.0;
                double d_response, d_queue;
                if (delay) {
                    d_response = d_demand;
                    d_queue = d_rho;
                } else {
                    d_response = (d_demand * denominator +
                                  demand[r] * d_utilization) /
                                 (denominator * denominator);
                    d_queue = (d_rho * denominator + rho[r] * d_utilization) /
                              (denominator * denominator);
                }
                sensitivity.add("RespT", metric, parameter, d_response);
                sensitivity.add("QLen", metric, parameter, d_queue);
            }
            sensitivity.add("Tput", metric,
                SensitivityData::parameter_key(station_name, sn.classes[r].name), 0.0);
        }
    }
    return sensitivity;
}

inline std::optional<SensitivityData> closed_sensitivities(qn::Network<double>& model) {
    const qn::NetworkStruct<double>& sn = model.get_struct();
    if (sn.classes.empty()) return std::nullopt;
    std::vector<int> population(sn.nclasses, 0);
    for (std::size_t r = 0; r < sn.nclasses; ++r) {
        if (!std::isfinite(sn.classes[r].population)) return std::nullopt;
        population[r] = static_cast<int>(std::llround(sn.classes[r].population));
    }

    const nc::PfParams<double> pf = nc::sn_get_product_form_params(sn);
    for (double servers : pf.S)
        if (std::isfinite(servers) && servers > 1.0) return std::nullopt;
    std::vector<double> think_time(sn.nclasses, 0.0);
    for (std::size_t z = 0; z < pf.Z.rows(); ++z)
        for (std::size_t r = 0; r < sn.nclasses; ++r)
            think_time[r] += pf.Z(z, r);
    const pfqn::SensResult<double> result = pfqn::pfqn_sens(
        pf.D, population, think_time);

    struct Parameter {
        std::size_t index;
        double chain;
        std::string key;
    };
    std::vector<Parameter> parameters;
    for (std::size_t j = 0; j < pf.queue_stations.size(); ++j) {
        const std::size_t station = pf.queue_stations[j] - 1;
        for (std::size_t s = 0; s < sn.nclasses; ++s) {
            const double rate = sn.rates(station, s);
            if (!std::isfinite(rate) || rate <= 0.0 || pf.D(j, s) <= 0.0) continue;
            parameters.push_back({j * sn.nclasses + s, -pf.D(j, s) / rate,
                SensitivityData::parameter_key(sn.stations[station].name,
                                               sn.classes[s].name)});
        }
    }

    SensitivityData sensitivity;
    for (std::size_t i = 0; i < pf.queue_stations.size(); ++i) {
        const std::size_t station = pf.queue_stations[i] - 1;
        const std::string& station_name = sn.stations[station].name;
        for (std::size_t r = 0; r < sn.nclasses; ++r) {
            if (pf.D(i, r) <= 0.0) continue;
            const std::string metric = SensitivityData::metric_key(
                station_name, sn.classes[r].name);
            for (const Parameter& parameter : parameters) {
                sensitivity.add("RespT", metric, parameter.key,
                                result.dR[parameter.index](i, r) * parameter.chain);
                sensitivity.add("QLen", metric, parameter.key,
                                result.dQ[parameter.index](i, r) * parameter.chain);
                sensitivity.add("Tput", metric, parameter.key,
                                result.dX(r, parameter.index) * parameter.chain);
                sensitivity.add("Util", station_name, parameter.key,
                                result.dU[parameter.index](i, r) * parameter.chain);
            }
        }
    }
    return sensitivity;
}

inline std::optional<SensitivityData> ctmc_sensitivities(qn::Network<double>& model,
                                                         ctmc::CtmcOptions options = {}) {
    const qn::NetworkStruct<double>& sn = model.get_struct();
    ctmc::CtmcSolution<double> base;
    Matrix<double> state_space;
    try {
        base = ctmc::solver_ctmc_analyzer(sn, options);
        state_space = ctmc::ctmc_state_space_aggr(sn, base.chain.space);
    } catch (const std::exception&) {
        return std::nullopt;
    }
    if (state_space.empty()) return std::nullopt;

    SensitivityData sensitivity;
    for (std::size_t pstation = 0; pstation < sn.nstations; ++pstation) {
        for (std::size_t pclass = 0; pclass < sn.nclasses; ++pclass) {
            const double rate = sn.rates(pstation, pclass);
            if (!std::isfinite(rate) || rate <= 0.0 ||
                !sn.has_service_law(pstation, pclass) ||
                sn.service[pstation][pclass].type != lang::ProcessType::EXP)
                continue;
            const std::string parameter = SensitivityData::parameter_key(
                sn.stations[pstation].name, sn.classes[pclass].name);
            ctmc::CtmcSensParam<double> specification;
            specification.name = parameter;
            specification.value = rate;
            specification.set = [pstation, pclass](qn::NetworkStruct<double>& changed,
                                                   double value) {
                changed.service[pstation][pclass] = lang::Distrib<double>::exp_rate(value);
                changed.refresh_rates();
            };
            ctmc::CtmcSens<double> derivative;
            try {
                derivative = ctmc::solver_ctmc_sensitivity(sn, options, specification);
            } catch (const std::exception&) {
                continue;
            }
            for (std::size_t i = 0; i < sn.nstations; ++i) {
                const std::string& station = sn.stations[i].name;
                for (std::size_t r = 0; r < sn.nclasses; ++r) {
                    const std::size_t column = i * sn.nclasses + r;
                    if (column >= state_space.cols()) continue;
                    double value = 0.0;
                    bool finite = true;
                    for (std::size_t state = 0; state < state_space.rows(); ++state) {
                        if (!std::isfinite(state_space(state, column))) {
                            finite = false;
                            break;
                        }
                        value += derivative.dpi[state] * state_space(state, column);
                    }
                    if (finite)
                        sensitivity.add("QLen", SensitivityData::metric_key(
                            station, sn.classes[r].name), parameter, value);
                }
            }
        }
    }
    return sensitivity;
}

inline std::optional<SensitivityData> compute_model_sensitivities(
    qn::Network<double> model, bool use_ctmc = false) {
    try {
        std::optional<SensitivityData> result = open_sensitivities(model);
        if (result) return result;
        result = closed_sensitivities(model);
        if (result) return result;
    } catch (const std::exception&) {
    }
    return use_ctmc ? ctmc_sensitivities(model) : std::nullopt;
}

}  // namespace opt
}  // namespace line

#endif
