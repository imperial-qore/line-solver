/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_INFERENCE_PARAM_ESTIMATOR_H
#define LINE_INFERENCE_PARAM_ESTIMATOR_H

/**
 * Sample collection and estimator orchestration for queueing-network parameters.
 *
 * This is the workflow contract of MATLAB `ParamEstimator`: sampled metrics are
 * stored by 1-based node and class, aggregate metrics have their own node row,
 * interpolation aligns every series on the union of timestamps, `auto_method`
 * applies the reference precedence, and `estimate_at` writes inferred means back
 * into the model. Estimators are typed callables; all eleven MATLAB methods are
 * installed by default.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "line/api/infer/infer_compute_ql_at_arrival.h"
#include "line/api/infer/infer_fmlps.h"
#include "line/api/infer/infer_gibbs.h"
#include "line/api/infer/infer_mlps.h"
#include "line/api/infer/infer_qmle.h"
#include "line/api/infer/infer_variational.h"
#include "line/api/sn/sn_rt_stations.h"
#include "line/inference/sampled_metric.h"
#include "line/lang/dist_scale_rate.h"
#include "line/lang/qn/network_builder.h"
#include "line/solvers/mva/solver_mva_runner.h"
#include "line/util/error.h"
#include "line/util/lu.h"
#include "line/util/matrix.h"
#include "line/util/neldermead.h"

namespace line {
namespace infer {

struct EstimatorOptions {
    int verbose = 1;
    std::string method = "ubr";
    std::string variant = "default";
    std::size_t iter_max = 1000;
    double tol = 1e-3;
    std::size_t open_population = 100;
    std::uint64_t random_seed = 5489u;
    std::vector<double> x0;
    GibbsOptions gibbs;

    /** Options of the variational estimator ("vi"). */
    VariationalOptions<double> variational;
    /** Probability that a queue-length reading is faulty ("vi"). */
    double epsilon = 0.05;
    /** Shape of the Gamma prior placed on each estimated rate ("vi"). */
    double prior_shape = 1.0;
    /** Posterior Gamma parameters left behind by the variational estimator. */
    std::vector<double> posterior_alpha, posterior_beta;
    /** Evidence lower bound per iteration, left behind by the same. */
    std::vector<double> bound;
};

namespace detail {

inline std::vector<double> spline_not_a_knot(const std::vector<double>& x,
                                              const std::vector<double>& y,
                                              const std::vector<double>& xq) {
    const std::size_t n = x.size();
    if (n != y.size())
        throw InputError("ParamEstimator.interpolate: sample times and data disagree in length");
    if (n < 2)
        throw InputError("ParamEstimator.interpolate: spline interpolation needs at least two samples");
    for (std::size_t i = 1; i < n; ++i)
        if (!(x[i] > x[i - 1]))
            throw InputError(
                "ParamEstimator.interpolate: sample timestamps must be strictly increasing");

    std::vector<double> out(xq.size(), 0.0);
    if (n == 2) {
        const double slope = (y[1] - y[0]) / (x[1] - x[0]);
        for (std::size_t q = 0; q < xq.size(); ++q) out[q] = y[0] + slope * (xq[q] - x[0]);
        return out;
    }
    if (n == 3) {
        const double d01 = (y[1] - y[0]) / (x[1] - x[0]);
        const double d12 = (y[2] - y[1]) / (x[2] - x[1]);
        const double d012 = (d12 - d01) / (x[2] - x[0]);
        for (std::size_t q = 0; q < xq.size(); ++q) {
            const double z = xq[q];
            out[q] = y[0] + d01 * (z - x[0]) + d012 * (z - x[0]) * (z - x[1]);
        }
        return out;
    }

    std::vector<double> h(n - 1, 0.0), delta(n - 1, 0.0);
    for (std::size_t i = 0; i + 1 < n; ++i) {
        h[i] = x[i + 1] - x[i];
        delta[i] = (y[i + 1] - y[i]) / h[i];
    }

    Matrix<double> A(n, n, 0.0);
    std::vector<double> rhs(n, 0.0);
    A(0, 0) = h[1];
    A(0, 1) = -(h[0] + h[1]);
    A(0, 2) = h[0];
    for (std::size_t i = 1; i + 1 < n; ++i) {
        A(i, i - 1) = h[i - 1];
        A(i, i) = 2.0 * (h[i - 1] + h[i]);
        A(i, i + 1) = h[i];
        rhs[i] = 6.0 * (delta[i] - delta[i - 1]);
    }
    A(n - 1, n - 3) = h[n - 2];
    A(n - 1, n - 2) = -(h[n - 3] + h[n - 2]);
    A(n - 1, n - 1) = h[n - 3];
    const std::vector<double> second = solve(A, rhs);

    for (std::size_t q = 0; q < xq.size(); ++q) {
        const double z = xq[q];
        std::size_t i = 0;
        if (z >= x[n - 1]) {
            i = n - 2;
        } else if (z > x[0]) {
            i = static_cast<std::size_t>(std::upper_bound(x.begin(), x.end(), z) - x.begin() - 1);
        }
        const double s = z - x[i];
        const double b = delta[i] - h[i] * (2.0 * second[i] + second[i + 1]) / 6.0;
        const double c = second[i] / 2.0;
        const double d = (second[i + 1] - second[i]) / (6.0 * h[i]);
        out[q] = y[i] + s * (b + s * (c + s * d));
    }
    return out;
}

inline std::vector<double> least_squares_passive(const Matrix<double>& A,
                                                  const std::vector<double>& b,
                                                  std::vector<bool>& passive) {
    std::vector<std::size_t> cols;
    for (std::size_t j = 0; j < passive.size(); ++j)
        if (passive[j]) cols.push_back(j);
    std::vector<double> z(passive.size(), 0.0);
    if (cols.empty()) return z;

    // Modified Gram-Schmidt avoids squaring the condition number through the
    // normal equations. Rank-deficient columns leave the passive set; this is
    // the Lawson-Hanson active-set convention for a singular subproblem.
    Matrix<double> r(cols.size(), cols.size(), 0.0);
    std::vector<std::vector<double>> q;
    std::vector<std::size_t> independent;
    for (const std::size_t col : cols) {
        std::vector<double> v(A.rows(), 0.0);
        double column_norm2 = 0.0;
        for (std::size_t i = 0; i < A.rows(); ++i) {
            v[i] = A(i, col);
            column_norm2 += v[i] * v[i];
        }

        std::vector<double> projection(q.size(), 0.0);
        for (int pass = 0; pass < 2; ++pass) {
            for (std::size_t k = 0; k < q.size(); ++k) {
                double coefficient = 0.0;
                for (std::size_t i = 0; i < A.rows(); ++i) coefficient += q[k][i] * v[i];
                projection[k] += coefficient;
                for (std::size_t i = 0; i < A.rows(); ++i) v[i] -= coefficient * q[k][i];
            }
        }

        double residual_norm2 = 0.0;
        for (const double value : v) residual_norm2 += value * value;
        const double residual_norm = std::sqrt(residual_norm2);
        const double rank_tol = 1e-12 * std::max(1.0, std::sqrt(column_norm2));
        if (residual_norm <= rank_tol) {
            passive[col] = false;
            continue;
        }

        const std::size_t a = independent.size();
        for (std::size_t k = 0; k < q.size(); ++k) r(k, a) = projection[k];
        r(a, a) = residual_norm;
        for (double& value : v) value /= residual_norm;
        q.push_back(std::move(v));
        independent.push_back(col);
    }

    std::vector<double> zp(independent.size(), 0.0);
    for (std::size_t k = 0; k < q.size(); ++k)
        for (std::size_t i = 0; i < A.rows(); ++i) zp[k] += q[k][i] * b[i];
    for (std::size_t offset = 0; offset < independent.size(); ++offset) {
        const std::size_t row = independent.size() - 1 - offset;
        for (std::size_t col = row + 1; col < independent.size(); ++col)
            zp[row] -= r(row, col) * zp[col];
        zp[row] /= r(row, row);
    }
    for (std::size_t a = 0; a < independent.size(); ++a) z[independent[a]] = zp[a];
    return z;
}

inline std::vector<double> nnls(const Matrix<double>& A, const std::vector<double>& b) {
    if (A.rows() != b.size()) throw InputError("ParamEstimator.ubr: NNLS dimensions disagree");
    const std::size_t n = A.cols();
    std::vector<double> x(n, 0.0), w(n, 0.0);
    std::vector<bool> passive(n, false);
    const double eps = 1e-12;
    const std::size_t limit = 30 * (n + 1) * (n + 1);

    const auto gradient = [&]() {
        std::vector<double> residual(b);
        for (std::size_t i = 0; i < A.rows(); ++i)
            for (std::size_t j = 0; j < n; ++j) residual[i] -= A(i, j) * x[j];
        for (std::size_t j = 0; j < n; ++j) {
            w[j] = 0.0;
            for (std::size_t i = 0; i < A.rows(); ++i) w[j] += A(i, j) * residual[i];
        }
    };

    gradient();
    for (std::size_t outer = 0; outer < limit; ++outer) {
        std::size_t enter = n;
        double best = eps;
        for (std::size_t j = 0; j < n; ++j)
            if (!passive[j] && w[j] > best) {
                best = w[j];
                enter = j;
            }
        if (enter == n) return x;
        passive[enter] = true;

        for (std::size_t inner = 0; inner < limit; ++inner) {
            const std::vector<bool> passive_before = passive;
            const std::vector<double> z = least_squares_passive(A, b, passive);
            for (std::size_t j = 0; j < n; ++j)
                if (passive_before[j] && !passive[j]) x[j] = 0.0;
            bool positive = true;
            for (std::size_t j = 0; j < n; ++j)
                if (passive[j] && z[j] <= eps) positive = false;
            if (positive) {
                x = z;
                break;
            }

            double alpha = std::numeric_limits<double>::infinity();
            for (std::size_t j = 0; j < n; ++j)
                if (passive[j] && z[j] <= eps)
                    alpha = std::min(alpha, x[j] / (x[j] - z[j]));
            for (std::size_t j = 0; j < n; ++j) x[j] += alpha * (z[j] - x[j]);
            for (std::size_t j = 0; j < n; ++j)
                if (passive[j] && x[j] <= eps) {
                    passive[j] = false;
                    x[j] = 0.0;
                }
        }
        gradient();
    }
    throw NumericError("ParamEstimator.ubr: NNLS did not converge");
}

}  // namespace detail

class ParamEstimator {
  public:
    using Samples = std::vector<std::vector<std::vector<SampledMetric>>>;
    using AggregateSamples = std::vector<std::vector<SampledMetric>>;
    using Estimator =
        std::function<Matrix<double>(ParamEstimator&, const std::vector<std::size_t>&)>;

    explicit ParamEstimator(qn::Network<double>& model,
                            const EstimatorOptions& options = EstimatorOptions())
        : options(options), model_(model) {
        const qn::NetworkStruct<double>& sn = model_.raw_struct();
        samples.resize(sn.nodes.size(),
                       std::vector<std::vector<SampledMetric>>(sn.classes.size()));
        samples_aggr.resize(sn.nodes.size());
        register_estimator("ubr", [](ParamEstimator& self, const std::vector<std::size_t>& nodes) {
            return self.estimator_ubr(nodes);
        });
        register_estimator("qmle", [](ParamEstimator& self, const std::vector<std::size_t>& nodes) {
            return self.estimator_qmle(nodes);
        });
        register_estimator("gibbs", [](ParamEstimator& self, const std::vector<std::size_t>& nodes) {
            return self.estimator_gibbs(nodes);
        });
        register_estimator("mlps", [](ParamEstimator& self, const std::vector<std::size_t>& nodes) {
            return self.estimator_mlps(nodes);
        });
        register_estimator("fmlps", [](ParamEstimator& self, const std::vector<std::size_t>& nodes) {
            return self.estimator_fmlps(nodes);
        });
        register_estimator("ubo", [](ParamEstimator& self, const std::vector<std::size_t>& nodes) {
            return self.estimator_ubo(nodes);
        });
        register_estimator("erps", [](ParamEstimator& self, const std::vector<std::size_t>& nodes) {
            return self.estimator_erps(nodes);
        });
        register_estimator("ekf", [](ParamEstimator& self, const std::vector<std::size_t>& nodes) {
            return self.estimator_ekf(nodes);
        });
        register_estimator("mcmc", [](ParamEstimator& self, const std::vector<std::size_t>& nodes) {
            return self.estimator_mcmc(nodes);
        });
        register_estimator("mle", [](ParamEstimator& self, const std::vector<std::size_t>& nodes) {
            return self.estimator_mle(nodes);
        });
        register_estimator("vi", [](ParamEstimator& self, const std::vector<std::size_t>& nodes) {
            return self.estimator_variational(nodes);
        });
    }

    EstimatorOptions options;
    Samples samples;
    AggregateSamples samples_aggr;

    void add_samples(const SampledMetric& sample) {
        require_node(sample.node, "add_samples");
        if (sample.is_aggregate()) {
            samples_aggr[sample.node - 1].push_back(sample);
        } else {
            require_class(sample.jobclass, "add_samples");
            samples[sample.node - 1][sample.jobclass - 1].push_back(sample);
        }
    }

    Samples& get_data() { return samples; }
    const Samples& get_data() const { return samples; }
    AggregateSamples& get_data_aggr() { return samples_aggr; }
    const AggregateSamples& get_data_aggr() const { return samples_aggr; }

    SampledMetric* get_arvr(std::size_t node, std::size_t jobclass) {
        return first_metric(node, jobclass, lang::MetricType::ArvR);
    }
    SampledMetric* get_util(std::size_t node, std::size_t jobclass) {
        return first_metric(node, jobclass, lang::MetricType::Util);
    }
    SampledMetric* get_respt(std::size_t node, std::size_t jobclass) {
        return first_metric(node, jobclass, lang::MetricType::RespT);
    }
    SampledMetric* get_tput(std::size_t node, std::size_t jobclass) {
        return first_metric(node, jobclass, lang::MetricType::Tput);
    }
    SampledMetric* get_aggr_util(std::size_t node) {
        return first_aggr_metric(node, lang::MetricType::Util, nullptr);
    }
    SampledMetric* get_aggr_qlen(std::size_t node, const ConditionEvent* event = nullptr) {
        return first_aggr_metric(node, lang::MetricType::QLen, event);
    }

    std::vector<SampledMetric*> get_qlen(std::size_t node, std::size_t jobclass,
                                         const ConditionEvent* event = nullptr) {
        require_node(node, "get_qlen");
        require_class(jobclass, "get_qlen");
        std::vector<SampledMetric*> out;
        std::vector<SampledMetric>& cell = samples[node - 1][jobclass - 1];
        for (SampledMetric& sample : cell) {
            if (sample.type != lang::MetricType::QLen) continue;
            if (event != nullptr && (!sample.cond.has_value() || *sample.cond != *event)) continue;
            out.push_back(&sample);
            if (event != nullptr) break;
        }
        return out;
    }

    std::string auto_method() {
        bool has_arvr = false, has_respt = false, has_util = false;
        bool has_qlen = false, has_tput = false, has_trace = false;
        bool has_aggr_util = false, has_aggr_qlen = false;
        for (const auto& node : samples)
            for (const auto& cls : node)
                for (const SampledMetric& sample : cls) {
                    if (sample.type == lang::MetricType::ArvR) has_arvr = true;
                    if (sample.type == lang::MetricType::RespT) has_respt = true;
                    if (sample.type == lang::MetricType::Util) has_util = true;
                    if (sample.type == lang::MetricType::QLen) has_qlen = true;
                    if (sample.type == lang::MetricType::Tput) has_tput = true;
                    if (sample.is_trace()) has_trace = true;
                }
        for (const auto& node : samples_aggr)
            for (const SampledMetric& sample : node) {
                if (sample.type == lang::MetricType::Util) has_aggr_util = true;
                if (sample.type == lang::MetricType::QLen) has_aggr_qlen = true;
            }
        (void)has_tput;

        if (has_trace && has_respt && has_aggr_qlen)
            options.method = "erps";
        else if (has_trace && has_respt && has_arvr)
            options.method = "mlps";
        else if (has_arvr && has_respt && (has_util || has_aggr_util))
            options.method = "ubo";
        else if (has_arvr && (has_util || has_aggr_util))
            options.method = "ubr";
        else if (has_qlen)
            options.method = "qmle";
        else
            throw InputError(
                "ParamEstimator.auto_method: insufficient data to select an estimation method");
        return options.method;
    }

    void interpolate() {
        std::set<double> union_set;
        for_each_sample([&](SampledMetric& sample) {
            union_set.insert(sample.t.begin(), sample.t.end());
        });
        if (union_set.empty()) return;
        const std::vector<double> tunion(union_set.begin(), union_set.end());
        for_each_sample([&](SampledMetric& sample) {
            sample.data = detail::spline_not_a_knot(sample.t, sample.data, tunion);
            sample.t = tunion;
        });
    }

    void register_estimator(const std::string& method, const Estimator& estimator) {
        estimators_[method] = estimator;
    }

    bool has_estimator(const std::string& method) const {
        return estimators_.find(method) != estimators_.end();
    }

    Matrix<double> estimate_at(const std::vector<std::size_t>& nodes) {
        const auto it = estimators_.find(options.method);
        if (it == estimators_.end())
            throw UnsupportedError("ParamEstimator: estimator '" + options.method +
                                   "' has no C++ adapter");
        Matrix<double> values = it->second(*this, nodes);
        const std::size_t classes = model_.raw_struct().classes.size();
        if (values.rows() != nodes.size() || values.cols() != classes)
            throw InputError(
                "ParamEstimator.estimate_at: estimator result must be nodes by classes");

        for (std::size_t n = 0; n < nodes.size(); ++n) {
            const std::size_t node = nodes[n];
            require_node(node, "estimate_at");
            qn::NetworkStruct<double>& sn = model_.raw_struct();
            const std::size_t station = sn.nodes[node - 1].station;
            if (station == 0)
                throw InputError("ParamEstimator.estimate_at: node '" +
                                 sn.nodes[node - 1].name + "' is not a station");
            if (sn.stations[station - 1].nodetype == lang::NodeType::Source) continue;
            for (std::size_t r = 0; r < classes; ++r) {
                const double target = values(n, r);
                if (!std::isfinite(target) || !(target > 0.0)) continue;
                const lang::Distrib<double> current = sn.service[station - 1][r];
                if (current.disabled)
                    throw InputError("ParamEstimator.estimate_at: node '" +
                                     sn.nodes[node - 1].name + "' has no class " +
                                     std::to_string(r + 1) + " service process");
                if (current.type == lang::ProcessType::IMMEDIATE ||
                    !std::isfinite(current.mean)) {
                    model_.set_service(node, r + 1, lang::Distrib<double>::exp_mean(target));
                } else {
                    switch (current.type) {
                        case lang::ProcessType::EXP:
                        case lang::ProcessType::ERLANG:
                        case lang::ProcessType::HYPEREXP:
                        case lang::ProcessType::COXIAN:
                        case lang::ProcessType::COX2:
                        case lang::ProcessType::APH:
                        case lang::ProcessType::PH:
                        case lang::ProcessType::MAP:
                        case lang::ProcessType::MMPP2:
                            break;
                        default:
                            throw UnsupportedError(
                                "ParamEstimator.estimate_at: " +
                                std::string(lang::process_to_text(current.type)) +
                                " has no MATLAB setMean contract");
                    }
                    const double factor = current.mean / target;
                    model_.set_service(node, r + 1, lang::dist_scale_rate(current, factor));
                }
            }
        }
        return values;
    }

    static EstimatorOptions default_options() { return EstimatorOptions(); }

    static std::string get_required_metrics(const std::string& method) {
        static const std::map<std::string, std::string> descriptions = {
            {"ubr", "ArvR (per-class) + Util (per-class or aggregate)"},
            {"ubo", "ArvR (per-class) + RespT (per-class) + Util (aggregate)"},
            {"erps",
             "RespT (per-class) + QLen (aggregate, conditional on class arrivals). PS "
             "stations only."},
            {"ekf", "RespT (per-class) + Util (aggregate). Sequential/recursive estimation."},
            {"mcmc", "QLen (aggregate). Gibbs sampling with MCMC. Open/mixed via closed "
                     "equivalence."},
            {"mle", "ArvR (per-class) + RespT (per-class) + Util (aggregate)"},
            {"vi", "QLen (per-class, timeseries) at every station. Variational "
                   "inference over transition counts; noisy readings, Gamma posteriors."},
            {"mlps", "ArvR (per-class, trace) + RespT (per-class, trace). PS stations only. "
                     "Open/mixed via closed equivalence."},
            {"fmlps", "ArvR (per-class, trace) + RespT (per-class, trace). PS stations only. "
                      "Open/mixed via closed equivalence."},
            {"qmle", "QLen (per-class). Open/mixed via closed equivalence (Z_r = N_r / "
                     "lambda_r)."},
            {"gibbs", "ArvR (per-class, trace) + RespT (per-class, trace) + Tput "
                      "(per-class). Gibbs sampling."}};
        const auto it = descriptions.find(method);
        return it == descriptions.end() ? "Unknown method: " + method : it->second;
    }

  private:
    qn::Network<double>& model_;
    std::map<std::string, Estimator> estimators_;

    void require_node(std::size_t node, const char* method) const {
        if (node == 0 || node > samples.size())
            throw InputError(std::string("ParamEstimator.") + method + ": node index is out of range");
    }

    void require_class(std::size_t jobclass, const char* method) const {
        const std::size_t classes = samples.empty() ? model_.raw_struct().classes.size()
                                                   : samples[0].size();
        if (jobclass == 0 || jobclass > classes)
            throw InputError(std::string("ParamEstimator.") + method +
                             ": class index is out of range");
    }

    SampledMetric* first_metric(std::size_t node, std::size_t jobclass,
                                lang::MetricType type) {
        require_node(node, "metric lookup");
        require_class(jobclass, "metric lookup");
        for (SampledMetric& sample : samples[node - 1][jobclass - 1])
            if (sample.type == type) return &sample;
        return nullptr;
    }

    SampledMetric* first_aggr_metric(std::size_t node, lang::MetricType type,
                                     const ConditionEvent* event) {
        require_node(node, "aggregate metric lookup");
        for (SampledMetric& sample : samples_aggr[node - 1]) {
            if (sample.type != type) continue;
            if (event != nullptr && (!sample.cond.has_value() || *sample.cond != *event)) continue;
            return &sample;
        }
        return nullptr;
    }

    template <class F>
    void for_each_sample(F fn) {
        for (auto& node : samples)
            for (auto& cls : node)
                for (SampledMetric& sample : cls) fn(sample);
        for (auto& node : samples_aggr)
            for (SampledMetric& sample : node) fn(sample);
    }

    void require_single_station(const std::vector<std::size_t>& nodes, const char* method,
                                std::size_t& node, std::size_t& station) const {
        if (nodes.size() != 1)
            throw InputError(std::string("ParamEstimator.") + method +
                             ": the estimator accepts exactly one station");
        node = nodes[0];
        require_node(node, method);
        const qn::NetworkStruct<double>& sn = model_.raw_struct();
        station = sn.nodes[node - 1].station;
        if (station == 0)
            throw InputError(std::string("ParamEstimator.") + method +
                             ": the target node is not a station");
    }

    void effective_population_and_think(std::vector<double>& population,
                                        std::vector<double>& think) {
        const qn::NetworkStruct<double>& sn = model_.get_struct();
        const std::size_t classes = sn.classes.size();
        population.assign(classes, 0.0);
        think.assign(classes, 0.0);
        for (std::size_t r = 0; r < classes; ++r)
            population[r] = std::isfinite(sn.classes[r].population)
                                ? sn.classes[r].population
                                : static_cast<double>(options.open_population);

        for (const qn::NodeDef& node : sn.nodes) {
            if (node.station == 0) continue;
            if (node.nodetype == lang::NodeType::Delay) {
                for (std::size_t r = 0; r < classes; ++r)
                    if (std::isfinite(sn.classes[r].population) &&
                        !sn.service[node.station - 1][r].disabled)
                        think[r] += sn.service[node.station - 1][r].mean;
            } else if (node.nodetype == lang::NodeType::Source) {
                for (std::size_t r = 0; r < classes; ++r)
                    if (!std::isfinite(sn.classes[r].population)) {
                        const lang::Distrib<double>& arrival = sn.service[node.station - 1][r];
                        if (arrival.disabled || !(arrival.mean > 0.0) ||
                            !std::isfinite(arrival.mean))
                            throw InputError("ParamEstimator: an open class has no finite source "
                                             "interarrival mean");
                        think[r] = population[r] * arrival.mean;
                    }
            }
        }
        for (std::size_t r = 0; r < classes; ++r)
            if (!(think[r] > 0.0) || !std::isfinite(think[r]))
                throw InputError("ParamEstimator: class " + std::to_string(r + 1) +
                                 " has no positive finite think time or open equivalent");
    }

    qn::Network<double> build_closed_equivalent(std::size_t node,
                                                 std::size_t& equivalent_queue) {
        std::vector<double> population, think;
        effective_population_and_think(population, think);
        const qn::NetworkStruct<double>& sn = model_.get_struct();
        const std::size_t station = sn.nodes[node - 1].station;
        qn::Network<double> equivalent("closed_equiv");
        const std::size_t delay = equivalent.add_delay("Think");
        equivalent_queue = equivalent.add_queue("Queue1", lang::SchedStrategy::PS);
        equivalent.set_number_of_servers(equivalent_queue, sn.stations[station - 1].nservers);
        std::vector<std::size_t> classes(population.size(), 0);
        for (std::size_t r = 0; r < population.size(); ++r) {
            classes[r] = equivalent.add_closed_class("Class" + std::to_string(r + 1),
                                                     population[r], delay);
            equivalent.set_service(delay, classes[r],
                                   lang::Distrib<double>::exp_rate(1.0 / think[r]));
            equivalent.set_service(equivalent_queue, classes[r], sn.service[station - 1][r]);
        }
        qn::RoutingMatrix<double> routing;
        for (std::size_t r = 0; r < classes.size(); ++r) {
            routing.set(classes[r], classes[r], delay, equivalent_queue, 1.0);
            routing.set(classes[r], classes[r], equivalent_queue, delay, 1.0);
        }
        equivalent.link(routing);
        return equivalent;
    }

    std::vector<api::MlpsSample> mlps_samples(std::size_t node, const char* method) {
        const std::size_t classes = model_.raw_struct().classes.size();
        std::vector<double> arrivals, response;
        std::vector<std::size_t> labels;
        for (std::size_t r = 0; r < classes; ++r) {
            SampledMetric* arvr = get_arvr(node, r + 1);
            SampledMetric* respt = get_respt(node, r + 1);
            if (arvr == nullptr || respt == nullptr)
                throw InputError(std::string("ParamEstimator.") + method +
                                 ": arrival and response-time traces are required for class " +
                                 std::to_string(r + 1));
            if (!arvr->is_trace() || !respt->is_trace())
                throw InputError(std::string("ParamEstimator.") + method +
                                 ": arrival and response-time metrics must use trace format");
            if (arvr->data.size() != respt->data.size())
                throw InputError(std::string("ParamEstimator.") + method +
                                 ": arrival and response-time traces have different lengths");
            arrivals.insert(arrivals.end(), arvr->data.begin(), arvr->data.end());
            response.insert(response.end(), respt->data.begin(), respt->data.end());
            labels.insert(labels.end(), arvr->data.size(), r);
        }
        std::vector<long> ids(arrivals.size(), 0);
        for (std::size_t i = 0; i < ids.size(); ++i) ids[i] = static_cast<long>(i + 1);
        const Matrix<double> qlen =
            infer_compute_ql_at_arrival(arrivals, ids, response, ids, labels, classes);
        std::vector<std::size_t> order(arrivals.size(), 0);
        for (std::size_t i = 0; i < order.size(); ++i) order[i] = i;
        std::stable_sort(order.begin(), order.end(),
                         [&](std::size_t a, std::size_t b) { return arrivals[a] < arrivals[b]; });

        std::vector<api::MlpsSample> out;
        for (const std::size_t i : order) {
            if (!(response[i] > 0.0)) continue;
            api::MlpsSample sample;
            sample.rt = response[i];
            sample.cls = labels[i] + 1;
            sample.ql.assign(classes, 0.0);
            for (std::size_t r = 0; r < classes; ++r) sample.ql[r] = qlen(i, r);
            out.push_back(sample);
        }
        if (out.empty())
            throw InputError(std::string("ParamEstimator.") + method +
                             ": no positive response-time observations remain");
        return out;
    }

    /**
     * Variational inference for Markovian queueing networks (Perez-Casale, AAP
     * 53(3), 2021). The network is translated into the transition set
     * eta=(i,j,c) of the paper, with lambda_eta = mu_{i,c} p^c_{i,j}; routing
     * probabilities are taken as known from the model and only the station
     * rates of the requested nodes are estimated. The data are QLen
     * timeseries, one per (node, class), read as exact with probability
     * 1-epsilon and uniform over the remaining feasible values otherwise.
     */
    Matrix<double> estimator_variational(const std::vector<std::size_t>& nodes) {
        if (nodes.empty()) throw InputError("ParamEstimator.vi: no stations were requested");
        const qn::NetworkStruct<double>& sn = model_.get_struct();
        const std::size_t M = sn.nstations, R = sn.nclasses, MR = M * R;

        // station-to-station routing; the pseudo-closed sink-to-source
        // feedback is not a job transition
        const Matrix<double> rtst = api::sn_rt_stations(sn).rtst;
        std::vector<int> sched(M, 1);
        std::vector<bool> is_source(M, false);
        for (std::size_t i = 0; i < M; ++i) {
            const lang::SchedStrategy ss = sn.stations[i].sched;
            if (ss == lang::SchedStrategy::INF) {
                sched[i] = 0;
            } else if (ss == lang::SchedStrategy::EXT) {
                sched[i] = 2;
                is_source[i] = true;
            } else if (ss == lang::SchedStrategy::PS || ss == lang::SchedStrategy::FCFS ||
                       ss == lang::SchedStrategy::DPS || ss == lang::SchedStrategy::GPS ||
                       ss == lang::SchedStrategy::SIRO || ss == lang::SchedStrategy::LCFS) {
                sched[i] = 1;
            } else {
                throw InputError("ParamEstimator.vi: unsupported scheduling at station " +
                                 std::to_string(i + 1));
            }
        }

        infer::VariationalSpec<double> spec;
        std::vector<double> probs;
        for (std::size_t c = 0; c < R; ++c) {
            for (std::size_t i = 0; i < M; ++i) {
                for (std::size_t j = 0; j < M; ++j) {
                    if (i == j || is_source[j]) continue;
                    const double p = rtst(i * R + c, j * R + c);
                    if (!(p > 0)) continue;
                    spec.arcs.push_back({{i + 1, j + 1, c + 1}});
                    probs.push_back(p);
                }
            }
        }
        if (spec.arcs.empty())
            throw InputError("ParamEstimator.vi: the model has no job transitions to infer from");

        // which station-class rates are being estimated
        std::vector<std::size_t> estimated(MR, 0), node_station(nodes.size(), 0);
        std::size_t P = 0;
        for (std::size_t n = 0; n < nodes.size(); ++n) {
            require_node(nodes[n], "vi");
            std::size_t st = 0;
            for (std::size_t i = 1; i <= M; ++i)
                if (sn.node_of_station(i) == nodes[n]) st = i;
            if (st == 0) throw InputError("ParamEstimator.vi: node " + std::to_string(nodes[n]) +
                                          " is not a station");
            node_station[n] = st - 1;
            for (std::size_t r = 0; r < R; ++r) {
                const double rate = sn.rates(st - 1, r);
                if (rate > 0 && std::isfinite(rate)) estimated[r * M + st - 1] = ++P;
            }
        }
        if (P == 0)
            throw InputError("ParamEstimator.vi: no station-class pair with a positive rate");

        const std::size_t narcs = spec.arcs.size();
        spec.routeprob = probs;
        spec.arcparam.assign(narcs, 0);
        spec.arcrate.assign(narcs, 0.0);
        for (std::size_t e = 0; e < narcs; ++e) {
            const std::size_t i = spec.arcs[e][0] - 1, c = spec.arcs[e][2] - 1;
            const std::size_t p = estimated[c * M + i];
            if (p > 0) {
                spec.arcparam[e] = p;
            } else {
                spec.arcrate[e] = sn.rates(i, c);
                if (!(spec.arcrate[e] > 0) || !std::isfinite(spec.arcrate[e]))
                    throw InputError("ParamEstimator.vi: station " + std::to_string(i + 1) +
                                     " class " + std::to_string(c + 1) +
                                     " has no usable rate to hold fixed");
            }
        }
        spec.sched = sched;
        spec.nservers.assign(M, 1.0);
        for (std::size_t i = 0; i < M; ++i) {
            const double k = sn.stations[i].nservers;
            spec.nservers[i] = (std::isfinite(k) && k > 0) ? k : 1.0;
        }

        // observations: QLen timeseries, one column per (station, class) pair
        std::set<double> tset;
        std::map<std::size_t, std::pair<std::vector<double>, std::vector<double>>> series;
        for (std::size_t i = 1; i <= M; ++i) {
            const std::size_t nd = sn.node_of_station(i);
            if (nd == 0) continue;
            for (std::size_t r = 0; r < R; ++r) {
                const std::vector<SampledMetric*> data = get_qlen(nd, r + 1);
                if (data.empty() || data[0]->t.empty()) continue;
                series[r * M + i - 1] = std::make_pair(data[0]->t, data[0]->data);
                for (const double tv : data[0]->t) tset.insert(tv);
            }
        }
        if (series.empty())
            throw InputError("ParamEstimator.vi: queue-length timeseries data is missing");
        spec.obsTimes.assign(tset.begin(), tset.end());
        const std::size_t K = spec.obsTimes.size();
        const double unobs = infer::VariationalSpec<double>::unobserved();
        spec.obsData = Matrix<double>(K, MR, unobs);
        for (std::map<std::size_t, std::pair<std::vector<double>, std::vector<double>>>::const_iterator
                 it = series.begin(); it != series.end(); ++it) {
            for (std::size_t k = 0; k < it->second.first.size(); ++k) {
                const std::vector<double>::const_iterator pos =
                    std::lower_bound(spec.obsTimes.begin(), spec.obsTimes.end(),
                                     it->second.first[k]);
                if (pos != spec.obsTimes.end() && *pos == it->second.first[k]) {
                    spec.obsData(static_cast<std::size_t>(pos - spec.obsTimes.begin()),
                                 it->first) = std::round(it->second.second[k]);
                }
            }
        }

        // population per class bounds both the contamination support and the load
        spec.obsRange.assign(MR, 1.0);
        spec.capacity.assign(MR, std::numeric_limits<double>::infinity());
        spec.x0 = Matrix<double>(M, R, 0.0);
        for (std::size_t r = 0; r < R; ++r) {
            const double njobs = sn.classes[r].population;
            double pop;
            const bool closed = std::isfinite(njobs);
            if (closed) {
                pop = njobs;
            } else {
                double peak = 1.0;
                for (std::size_t k = 0; k < K; ++k)
                    for (std::size_t i = 0; i < M; ++i)
                        if (spec.obsData(k, r * M + i) != unobs)
                            peak = std::max(peak, spec.obsData(k, r * M + i));
                pop = std::max(1.0, 2.0 * peak);
            }
            for (std::size_t i = 0; i < M; ++i) {
                spec.obsRange[r * M + i] = pop;
                if (closed) spec.capacity[r * M + i] = pop;
            }
            if (closed && njobs > 0) {
                const std::size_t ref = sn.classes[r].refstat;
                spec.x0((ref >= 1 && ref <= M) ? ref - 1 : 0, r) = njobs;
            }
        }

        // Gamma priors centred on the model's current rates
        spec.alpha0.assign(P, 0.0);
        spec.beta0.assign(P, 0.0);
        for (std::size_t i = 0; i < M; ++i) {
            for (std::size_t r = 0; r < R; ++r) {
                const std::size_t p = estimated[r * M + i];
                if (p > 0) {
                    spec.alpha0[p - 1] = options.prior_shape;
                    spec.beta0[p - 1] = options.prior_shape / sn.rates(i, r);
                }
            }
        }
        spec.epsilon = options.epsilon;

        const infer::VariationalResult<double> out =
            infer::infer_variational(spec, options.variational);
        options.posterior_alpha = out.alpha;
        options.posterior_beta = out.beta;
        options.bound = out.bound;

        Matrix<double> est(nodes.size(), R, 0.0);
        for (std::size_t n = 0; n < nodes.size(); ++n) {
            const std::size_t i = node_station[n];
            for (std::size_t r = 0; r < R; ++r) {
                const std::size_t p = estimated[r * M + i];
                if (p > 0) {
                    est(n, r) = out.mean_service_time[p - 1];
                } else if (sn.rates(i, r) > 0) {
                    est(n, r) = 1.0 / sn.rates(i, r);
                }
            }
        }
        return est;
    }

    Matrix<double> estimator_qmle(const std::vector<std::size_t>& nodes) {
        if (nodes.empty()) throw InputError("ParamEstimator.qmle: no stations were requested");
        std::vector<double> population, think;
        effective_population_and_think(population, think);
        const std::size_t classes = population.size();
        Matrix<double> qlen(nodes.size(), classes, 0.0);
        for (std::size_t n = 0; n < nodes.size(); ++n) {
            require_node(nodes[n], "qmle");
            for (std::size_t r = 0; r < classes; ++r) {
                const std::vector<SampledMetric*> data = get_qlen(nodes[n], r + 1);
                if (data.empty() || data[0]->data.empty())
                    throw InputError("ParamEstimator.qmle: queue-length data is missing for node " +
                                     std::to_string(nodes[n]) + " class " +
                                     std::to_string(r + 1));
                double sum = 0.0;
                for (const double value : data[0]->data) sum += value;
                qlen(n, r) = sum / static_cast<double>(data[0]->data.size());
            }
        }
        return infer_qmle(qlen, population, think);
    }

    Matrix<double> estimator_gibbs(const std::vector<std::size_t>& nodes) {
        std::size_t node = 0, station = 0;
        require_single_station(nodes, "gibbs", node, station);
        const qn::NetworkStruct<double>& sn = model_.get_struct();
        std::vector<GibbsTrace<double>> traces(sn.classes.size());
        for (std::size_t r = 0; r < sn.classes.size(); ++r) {
            SampledMetric* arvr = get_arvr(node, r + 1);
            SampledMetric* respt = get_respt(node, r + 1);
            SampledMetric* tput = get_tput(node, r + 1);
            if (arvr == nullptr || respt == nullptr || tput == nullptr)
                throw InputError("ParamEstimator.gibbs: arrival, response-time, and throughput "
                                 "data are required for class " + std::to_string(r + 1));
            if (!arvr->is_trace() || !respt->is_trace())
                throw InputError("ParamEstimator.gibbs: arrival and response-time metrics must "
                                 "use trace format");
            traces[r].arrival_ms.resize(arvr->data.size());
            for (std::size_t i = 0; i < arvr->data.size(); ++i)
                traces[r].arrival_ms[i] = arvr->data[i] * 1000.0;
            traces[r].respt_s = respt->data;
            traces[r].think_obs = tput->data;
        }
        GibbsOptions gibbs = options.gibbs;
        gibbs.tol = options.tol;
        pfqn::McRng rng(options.random_seed);
        const std::vector<double> demand =
            infer_gibbs(traces, sn.stations[station - 1].nservers, gibbs, rng);
        Matrix<double> out(1, demand.size(), 0.0);
        for (std::size_t r = 0; r < demand.size(); ++r) out(0, r) = demand[r];
        return out;
    }

    Matrix<double> estimator_mlps(const std::vector<std::size_t>& nodes) {
        std::size_t node = 0, station = 0;
        require_single_station(nodes, "mlps", node, station);
        const qn::NetworkStruct<double>& sn = model_.get_struct();
        if (sn.stations[station - 1].sched != lang::SchedStrategy::PS)
            throw InputError("ParamEstimator.mlps: the target station must use PS scheduling");
        std::vector<double> population, think;
        effective_population_and_think(population, think);
        std::vector<double> rates(think.size(), 0.0);
        for (std::size_t r = 0; r < think.size(); ++r) rates[r] = 1.0 / think[r];
        const std::vector<double> demand =
            api::infer_mlps(rates, sn.stations[station - 1].nservers,
                            mlps_samples(node, "mlps"));
        Matrix<double> out(1, demand.size(), 0.0);
        for (std::size_t r = 0; r < demand.size(); ++r) out(0, r) = demand[r];
        return out;
    }

    Matrix<double> estimator_fmlps(const std::vector<std::size_t>& nodes) {
        std::size_t node = 0, station = 0;
        require_single_station(nodes, "fmlps", node, station);
        const qn::NetworkStruct<double>& sn = model_.get_struct();
        if (sn.stations[station - 1].sched != lang::SchedStrategy::PS)
            throw InputError("ParamEstimator.fmlps: the target station must use PS scheduling");
        std::size_t equivalent_queue = 0;
        qn::Network<double> equivalent = build_closed_equivalent(node, equivalent_queue);
        const qn::NetworkStruct<double>& eqsn = equivalent.get_struct();
        const std::size_t eqstation = eqsn.nodes[equivalent_queue - 1].station;
        const std::vector<double> demand =
            api::infer_fmlps(eqsn, eqstation, mlps_samples(node, "fmlps"));
        Matrix<double> out(1, demand.size(), 0.0);
        for (std::size_t r = 0; r < demand.size(); ++r) out(0, r) = demand[r];
        return out;
    }

    Matrix<double> estimator_ubo(const std::vector<std::size_t>& nodes) {
        if (nodes.empty()) throw InputError("ParamEstimator.ubo: no stations were requested");
        const qn::NetworkStruct<double>& sn = model_.get_struct();
        const std::size_t stations = nodes.size(), classes = sn.classes.size();
        std::size_t samples_count = 0;
        std::vector<const SampledMetric*> util(stations, nullptr);
        std::vector<std::vector<const SampledMetric*>> arrivals(
            stations, std::vector<const SampledMetric*>(classes, nullptr));
        std::vector<std::vector<const SampledMetric*>> response(
            stations, std::vector<const SampledMetric*>(classes, nullptr));

        for (std::size_t i = 0; i < stations; ++i) {
            require_node(nodes[i], "ubo");
            const std::size_t station = sn.nodes[nodes[i] - 1].station;
            if (station == 0 || !std::isfinite(sn.stations[station - 1].nservers))
                throw InputError("ParamEstimator.ubo: every target must be a finite-server station");
            util[i] = get_aggr_util(nodes[i]);
            if (util[i] == nullptr)
                throw InputError("ParamEstimator.ubo: aggregate utilization is missing for node " +
                                 std::to_string(nodes[i]));
            if (i == 0) samples_count = util[i]->data.size();
            if (util[i]->data.size() != samples_count)
                throw InputError("ParamEstimator.ubo: sampled metrics have different sample "
                                 "counts; call interpolate first");
            for (std::size_t r = 0; r < classes; ++r) {
                arrivals[i][r] = get_arvr(nodes[i], r + 1);
                response[i][r] = get_respt(nodes[i], r + 1);
                if (arrivals[i][r] == nullptr || response[i][r] == nullptr)
                    throw InputError("ParamEstimator.ubo: arrival-rate and response-time data are "
                                     "required for node " + std::to_string(nodes[i]) + " class " +
                                     std::to_string(r + 1));
                if (arrivals[i][r]->data.size() != samples_count ||
                    response[i][r]->data.size() != samples_count)
                    throw InputError("ParamEstimator.ubo: sampled metrics have different sample "
                                     "counts; call interpolate first");
            }
        }

        std::vector<std::size_t> valid;
        for (std::size_t n = 0; n < samples_count; ++n) {
            bool keep = true;
            double total_arrival = 0.0;
            for (std::size_t i = 0; i < stations; ++i) {
                if (!std::isfinite(util[i]->data[n])) keep = false;
                for (std::size_t r = 0; r < classes; ++r)
                    total_arrival += arrivals[i][r]->data[n];
            }
            if (keep && total_arrival > 0.0) valid.push_back(n);
        }
        if (valid.empty()) throw InputError("ParamEstimator.ubo: no usable experiments remain");

        const std::size_t variables = stations * classes;
        Matrix<double> design(valid.size() * (classes + stations), variables, 0.0);
        std::vector<double> target(design.rows(), 0.0);
        std::size_t row = 0;
        for (const std::size_t n : valid) {
            std::vector<double> rho(stations, 0.0), beta(stations, 0.0);
            std::vector<double> lambda_class(classes, 0.0), end_to_end(classes, 0.0);
            for (std::size_t i = 0; i < stations; ++i) {
                const std::size_t station = sn.nodes[nodes[i] - 1].station;
                rho[i] = util[i]->data[n] * sn.stations[station - 1].nservers;
                if (rho[i] == 1.0)
                    throw NumericError("ParamEstimator.ubo: utilization reaches one");
                beta[i] = 1.0 / (1.0 - rho[i]);
                for (std::size_t r = 0; r < classes; ++r) {
                    lambda_class[r] += arrivals[i][r]->data[n];
                    end_to_end[r] += response[i][r]->data[n];
                }
            }
            double total_lambda = 0.0;
            for (const double value : lambda_class) total_lambda += value;
            for (std::size_t r = 0; r < classes; ++r) {
                const double weight = std::sqrt(lambda_class[r] / total_lambda);
                for (std::size_t i = 0; i < stations; ++i)
                    design(row, r * stations + i) = weight * beta[i];
                target[row++] = weight * end_to_end[r];
            }
            for (std::size_t i = 0; i < stations; ++i) {
                for (std::size_t r = 0; r < classes; ++r)
                    design(row, r * stations + i) = arrivals[i][r]->data[n];
                target[row++] = rho[i];
            }
        }
        const std::vector<double> fit = detail::nnls(design, target);
        Matrix<double> out(stations, classes, 0.0);
        for (std::size_t r = 0; r < classes; ++r)
            for (std::size_t i = 0; i < stations; ++i)
                out(i, r) = fit[r * stations + i];
        return out;
    }

    Matrix<double> estimator_erps(const std::vector<std::size_t>& nodes) {
        std::size_t node = 0, station = 0;
        require_single_station(nodes, "erps", node, station);
        const qn::NetworkStruct<double>& sn = model_.get_struct();
        if (sn.stations[station - 1].sched != lang::SchedStrategy::PS)
            throw InputError("ParamEstimator.erps: the target station must use PS scheduling");
        const std::size_t classes = sn.classes.size();
        std::vector<const SampledMetric*> response(classes, nullptr), qlen(classes, nullptr);
        double busy_sum = 0.0;
        std::size_t busy_count = 0;
        for (std::size_t r = 0; r < classes; ++r) {
            response[r] = get_respt(node, r + 1);
            const ConditionEvent arrival(node, r + 1, lang::EventType::ARV);
            qlen[r] = get_aggr_qlen(node, &arrival);
            if (response[r] == nullptr || qlen[r] == nullptr)
                throw InputError("ParamEstimator.erps: response-time and arrival-conditional "
                                 "aggregate queue-length data are required for class " +
                                 std::to_string(r + 1));
            if (response[r]->data.size() != qlen[r]->data.size())
                throw InputError("ParamEstimator.erps: sampled metrics have different sample "
                                 "counts; call interpolate first");
            for (const double q : qlen[r]->data) {
                if (q < 1.0)
                    throw InputError("ParamEstimator.erps: an arrival queue length must include "
                                     "the arriving job");
                busy_sum += q;
                ++busy_count;
            }
        }
        if (busy_count == 0) throw InputError("ParamEstimator.erps: no observations");
        const double busy =
            std::min(busy_sum / static_cast<double>(busy_count),
                     sn.stations[station - 1].nservers);
        if (!(busy > 0.0)) throw NumericError("ParamEstimator.erps: zero average busy cores");
        Matrix<double> out(1, classes, 0.0);
        for (std::size_t r = 0; r < classes; ++r) {
            double aa = 0.0, ab = 0.0;
            for (std::size_t i = 0; i < qlen[r]->data.size(); ++i) {
                const double regressor = qlen[r]->data[i] / busy;
                aa += regressor * regressor;
                ab += regressor * response[r]->data[i];
            }
            if (!(aa > 0.0)) throw NumericError("ParamEstimator.erps: zero queue-length regressor");
            out(0, r) = std::max(0.0, ab / aa);
        }
        return out;
    }

    std::vector<double> solver_measurement(const qn::NetworkStruct<double>& base,
                                           std::size_t station,
                                           const std::vector<double>& demand) {
        qn::NetworkStruct<double> sn = base;
        for (std::size_t r = 0; r < demand.size(); ++r) {
            if (!(demand[r] > 0.0) || !std::isfinite(demand[r]))
                throw NumericError("ParamEstimator: a solver-backed estimate reached a "
                                   "non-positive service demand");
            sn.set_service(station, r + 1, lang::Distrib<double>::exp_mean(demand[r]));
        }
        sn.refresh_rates();
        mva::MvaOptions solver_options;
        solver_options.method = "default";
        const mva::AvgResult<double> solved =
            mva::solver_mva_run_analyzer(sn, solver_options, Matrix<double>());
        std::vector<double> measurement(demand.size() + 1, 0.0);
        for (std::size_t r = 0; r < demand.size(); ++r) {
            measurement[r] = solved.RN(station - 1, r);
            measurement.back() += solved.UN(station - 1, r);
        }
        return measurement;
    }

    Matrix<double> estimator_ekf(const std::vector<std::size_t>& nodes) {
        std::size_t node = 0, station = 0;
        require_single_station(nodes, "ekf", node, station);
        const qn::NetworkStruct<double>& sn = model_.get_struct();
        const double servers = sn.stations[station - 1].nservers;
        if (!std::isfinite(servers))
            throw InputError("ParamEstimator.ekf: the target station must have finite servers");
        SampledMetric* aggregate = get_aggr_util(node);
        if (aggregate == nullptr)
            throw InputError("ParamEstimator.ekf: aggregate utilization data is missing");
        const std::size_t classes = sn.classes.size(), count = aggregate->data.size();
        std::vector<const SampledMetric*> arrivals(classes, nullptr), response(classes, nullptr);
        for (std::size_t r = 0; r < classes; ++r) {
            arrivals[r] = get_arvr(node, r + 1);
            response[r] = get_respt(node, r + 1);
            if (arrivals[r] == nullptr || response[r] == nullptr)
                throw InputError("ParamEstimator.ekf: arrival-rate and response-time data are "
                                 "required for class " + std::to_string(r + 1));
            if (arrivals[r]->data.size() != count || response[r]->data.size() != count)
                throw InputError("ParamEstimator.ekf: sampled metrics have different sample "
                                 "counts; call interpolate first");
        }

        std::vector<std::size_t> valid;
        for (std::size_t n = 0; n < count; ++n) {
            double throughput = 0.0;
            for (std::size_t r = 0; r < classes; ++r) throughput += arrivals[r]->data[n];
            if (std::isfinite(aggregate->data[n]) && throughput != 0.0) valid.push_back(n);
        }
        if (valid.empty()) throw InputError("ParamEstimator.ekf: no usable experiments remain");

        std::vector<double> x(classes, 0.0);
        if (!options.x0.empty()) {
            if (options.x0.size() != classes)
                throw InputError("ParamEstimator.ekf: x0 must contain one demand per class");
            x = options.x0;
        } else {
            pfqn::McRng rng(options.random_seed);
            for (std::size_t r = 0; r < classes; ++r) {
                double maximum = 0.0;
                for (const double value : response[r]->data) maximum = std::max(maximum, value);
                x[r] = pfqn::mc_uniform<double>(rng) * maximum;
                if (!(x[r] > 0.0)) x[r] = std::max(1e-9, maximum * 0.5);
            }
        }

        Matrix<double> covariance(classes, classes, 0.0);
        for (std::size_t r = 0; r < classes; ++r) covariance(r, r) = x[r] * x[r];
        const double delta = 1e-6;
        const std::size_t iterations = std::min(valid.size(), options.iter_max);
        for (std::size_t step = 0; step < iterations; ++step) {
            const std::size_t n = valid[step];
            Matrix<double> predicted_cov = covariance;
            for (std::size_t r = 0; r < classes; ++r) predicted_cov(r, r) += 0.001;
            const std::vector<double> predicted = solver_measurement(sn, station, x);

            Matrix<double> jacobian(classes + 1, classes, 0.0);
            for (std::size_t c = 0; c < classes; ++c) {
                std::vector<double> perturbed(x);
                perturbed[c] += delta;
                const std::vector<double> moved = solver_measurement(sn, station, perturbed);
                for (std::size_t k = 0; k < classes + 1; ++k)
                    jacobian(k, c) = (moved[k] - predicted[k]) / delta;
            }

            Matrix<double> innovation_cov(classes + 1, classes + 1, 0.0);
            for (std::size_t i = 0; i < classes + 1; ++i)
                for (std::size_t j = 0; j < classes + 1; ++j) {
                    double value = i == j ? 0.01 : 0.0;
                    for (std::size_t a = 0; a < classes; ++a)
                        for (std::size_t b = 0; b < classes; ++b)
                            value += jacobian(i, a) * predicted_cov(a, b) * jacobian(j, b);
                    innovation_cov(i, j) = value;
                }

            Matrix<double> gain(classes, classes + 1, 0.0);
            for (std::size_t r = 0; r < classes; ++r) {
                std::vector<double> rhs(classes + 1, 0.0);
                for (std::size_t k = 0; k < classes + 1; ++k)
                    for (std::size_t a = 0; a < classes; ++a)
                        rhs[k] += predicted_cov(r, a) * jacobian(k, a);
                const std::vector<double> solved = solve(innovation_cov, rhs);
                for (std::size_t k = 0; k < classes + 1; ++k) gain(r, k) = solved[k];
            }

            std::vector<double> residual(classes + 1, 0.0);
            for (std::size_t r = 0; r < classes; ++r)
                residual[r] = response[r]->data[n] - predicted[r];
            residual.back() = aggregate->data[n] * servers - predicted.back();
            for (std::size_t r = 0; r < classes; ++r) {
                double update = x[r];
                for (std::size_t k = 0; k < classes + 1; ++k)
                    update += gain(r, k) * residual[k];
                x[r] = std::max(0.4 * update, update);
            }
            double sum = 0.0;
            for (const double value : x) sum += value;
            if (sum < 0.0)
                for (double& value : x) value = -value;

            Matrix<double> next(classes, classes, 0.0);
            for (std::size_t i = 0; i < classes; ++i)
                for (std::size_t j = 0; j < classes; ++j) {
                    double value = 0.0;
                    for (std::size_t a = 0; a < classes; ++a) {
                        double ikh = i == a ? 1.0 : 0.0;
                        for (std::size_t k = 0; k < classes + 1; ++k)
                            ikh -= gain(i, k) * jacobian(k, a);
                        value += ikh * predicted_cov(a, j);
                    }
                    next(i, j) = value;
                }
            covariance = next;
        }
        Matrix<double> out(1, classes, 0.0);
        for (std::size_t r = 0; r < classes; ++r) out(0, r) = x[r];
        return out;
    }

    Matrix<double> estimator_mcmc(const std::vector<std::size_t>& nodes) {
        if (nodes.empty()) throw InputError("ParamEstimator.mcmc: no stations were requested");
        const qn::NetworkStruct<double>& sn = model_.get_struct();
        const std::size_t classes = sn.classes.size();
        std::vector<double> population, think;
        effective_population_and_think(population, think);
        (void)population;
        (void)think;

        Matrix<double> average(nodes.size(), classes, 0.0);
        std::size_t experiments = 0;
        for (std::size_t i = 0; i < nodes.size(); ++i) {
            require_node(nodes[i], "mcmc");
            SampledMetric* qlen = get_aggr_qlen(nodes[i]);
            if (qlen == nullptr || qlen->data.empty())
                throw InputError("ParamEstimator.mcmc: aggregate queue-length data is missing for "
                                 "node " + std::to_string(nodes[i]));
            if (i == 0)
                experiments = qlen->data.size();
            else if (qlen->data.size() != experiments)
                throw InputError("ParamEstimator.mcmc: sampled metrics have different sample "
                                 "counts; call interpolate first");
            double sum = 0.0;
            for (const double value : qlen->data) {
                if (!(value >= 0.0) || !std::isfinite(value))
                    throw InputError("ParamEstimator.mcmc: queue lengths must be finite and "
                                     "non-negative");
                sum += value;
            }
            for (std::size_t r = 0; r < classes; ++r)
                average(i, r) = sum / static_cast<double>(experiments);
        }

        double maximum = 0.0;
        for (std::size_t i = 0; i < average.rows(); ++i)
            for (std::size_t r = 0; r < average.cols(); ++r)
                maximum = std::max(maximum, average(i, r));
        if (!(maximum > 0.0))
            throw InputError("ParamEstimator.mcmc: at least one mean queue length must be positive");

        // In the MATLAB/JAR/Python reference, pfqn_mci is evaluated at the
        // unchanged sampleTheta for every grid point. Its normalizing constant
        // and the uniform prior therefore cancel from this coordinate draw.
        const std::size_t samples_count = 100, grid_count = 400;
        std::vector<Matrix<double>> theta(samples_count + 1,
                                          Matrix<double>(nodes.size(), classes, 0.0));
        pfqn::McRng rng(options.random_seed);
        for (std::size_t s = 0; s < samples_count; ++s) {
            Matrix<double> sample = theta[s];
            for (std::size_t i = 0; i < nodes.size(); ++i) {
                for (std::size_t r = 0; r < classes; ++r) {
                    const double exponent = static_cast<double>(experiments) * average(i, r);
                    std::vector<double> weights(grid_count + 1, 0.0);
                    double total = 0.0;
                    for (std::size_t k = 0; k <= grid_count; ++k) {
                        if (k == 0 && exponent > 0.0) continue;
                        const double ratio = static_cast<double>(k) /
                                             static_cast<double>(grid_count);
                        weights[k] = exponent == 0.0 ? 1.0 : std::pow(ratio, exponent);
                        total += weights[k];
                    }
                    double draw = pfqn::mc_uniform01(rng) * total;
                    std::size_t picked = grid_count;
                    for (std::size_t k = 0; k <= grid_count; ++k) {
                        draw -= weights[k];
                        if (draw < 0.0) {
                            picked = k;
                            break;
                        }
                    }
                    sample(i, r) = maximum * static_cast<double>(picked) /
                                   static_cast<double>(grid_count);
                }
            }
            theta[s + 1] = sample;
        }

        Matrix<double> out(nodes.size(), classes, 0.0);
        const std::size_t cutoff = samples_count / 2;
        for (std::size_t i = 0; i < nodes.size(); ++i) {
            for (std::size_t r = 0; r < classes; ++r) {
                double visit = 1.0;
                for (std::size_t c = 0; c < sn.chains.size(); ++c) {
                    if (r >= sn.chains[c].size() || !sn.chains[c][r] ||
                        c >= sn.nodevisits.size())
                        continue;
                    if (nodes[i] - 1 < sn.nodevisits[c].rows() &&
                        r < sn.nodevisits[c].cols()) {
                        const double candidate = sn.nodevisits[c](nodes[i] - 1, r);
                        if (candidate > 0.0) visit = candidate;
                    }
                    break;
                }
                for (std::size_t s = cutoff; s <= samples_count; ++s)
                    out(i, r) += theta[s](i, r) / visit;
                out(i, r) /= static_cast<double>(samples_count - cutoff + 1);
            }
        }
        return out;
    }

    Matrix<double> estimator_mle(const std::vector<std::size_t>& nodes) {
        std::size_t node = 0, station = 0;
        require_single_station(nodes, "mle", node, station);
        const qn::NetworkStruct<double>& sn = model_.get_struct();
        const double servers = sn.stations[station - 1].nservers;
        if (!std::isfinite(servers))
            throw InputError("ParamEstimator.mle: the target station must have finite servers");
        SampledMetric* aggregate = get_aggr_util(node);
        if (aggregate == nullptr)
            throw InputError("ParamEstimator.mle: aggregate utilization data is missing");
        const std::size_t classes = sn.classes.size(), count = aggregate->data.size();
        std::vector<const SampledMetric*> arrivals(classes, nullptr), response(classes, nullptr);
        std::vector<double> upper(classes, 0.0);
        for (std::size_t r = 0; r < classes; ++r) {
            arrivals[r] = get_arvr(node, r + 1);
            response[r] = get_respt(node, r + 1);
            if (arrivals[r] == nullptr || response[r] == nullptr)
                throw InputError("ParamEstimator.mle: arrival-rate and response-time data are "
                                 "required for class " + std::to_string(r + 1));
            if (arrivals[r]->data.size() != count || response[r]->data.size() != count)
                throw InputError("ParamEstimator.mle: sampled metrics have different sample "
                                 "counts; call interpolate first");
            for (const double value : response[r]->data)
                if (std::isfinite(value)) upper[r] = std::max(upper[r], value);
            if (!(upper[r] >= 1e-8))
                throw InputError("ParamEstimator.mle: response-time upper bounds must be positive");
        }

        std::vector<std::size_t> valid;
        for (std::size_t n = 0; n < count; ++n) {
            double throughput = 0.0;
            for (std::size_t r = 0; r < classes; ++r) throughput += arrivals[r]->data[n];
            if (std::isfinite(aggregate->data[n]) && throughput > 0.0) valid.push_back(n);
        }
        if (valid.empty()) throw InputError("ParamEstimator.mle: no usable experiments remain");

        std::vector<double> x0(classes, 0.0);
        if (!options.x0.empty()) {
            if (options.x0.size() != classes)
                throw InputError("ParamEstimator.mle: x0 must contain one demand per class");
            x0 = options.x0;
        } else {
            pfqn::McRng rng(options.random_seed);
            for (std::size_t r = 0; r < classes; ++r)
                x0[r] = std::max(1e-8, pfqn::mc_uniform01(rng) * upper[r]);
        }

        std::vector<Bound<double>> bounds(classes);
        for (std::size_t r = 0; r < classes; ++r) bounds[r] = bound_box(1e-8, upper[r]);
        const auto objective = [&](const std::vector<double>& x) {
            const std::vector<double> predicted = solver_measurement(sn, station, x);
            double value = 0.0;
            for (const std::size_t n : valid) {
                double total = 0.0;
                for (std::size_t r = 0; r < classes; ++r) total += arrivals[r]->data[n];
                for (std::size_t r = 0; r < classes; ++r) {
                    const double weight = arrivals[r]->data[n] / total;
                    const double residual = predicted[r] - response[r]->data[n];
                    value += weight * residual * residual;
                }
                const double residual = predicted.back() - aggregate->data[n] * servers;
                value += residual * residual;
            }
            return value;
        };
        NelderMeadOptions<double> nm = nelder_mead_defaults<double>();
        nm.max_iter = static_cast<unsigned>(std::min<std::size_t>(
            options.iter_max, std::numeric_limits<unsigned>::max()));
        nm.max_eval = std::max<unsigned>(1000u, nm.max_iter *
                                                   static_cast<unsigned>(10 * (classes + 1)));
        const NelderMeadResult<double> fit = nelder_mead_box(objective, x0, bounds, nm);
        Matrix<double> out(1, classes, 0.0);
        for (std::size_t r = 0; r < classes; ++r) out(0, r) = fit.x[r];
        return out;
    }

    Matrix<double> estimator_ubr(const std::vector<std::size_t>& nodes) {
        if (nodes.size() != 1)
            throw InputError("ParamEstimator.ubr: the estimator accepts exactly one station");
        const std::size_t node = nodes[0];
        require_node(node, "ubr");
        const qn::NetworkStruct<double>& sn = model_.raw_struct();
        const std::size_t station = sn.nodes[node - 1].station;
        if (station == 0)
            throw InputError("ParamEstimator.ubr: the target node is not a station");
        const double servers = sn.stations[station - 1].nservers;
        if (!std::isfinite(servers))
            throw InputError("ParamEstimator.ubr: the target station must have finite servers");
        const std::size_t classes = sn.classes.size();

        std::vector<const SampledMetric*> arrivals(classes, nullptr), utils(classes, nullptr);
        std::size_t count = 0;
        for (std::size_t r = 0; r < classes; ++r) {
            arrivals[r] = get_arvr(node, r + 1);
            if (arrivals[r] == nullptr)
                throw InputError("ParamEstimator.ubr: arrival-rate data is missing for class " +
                                 std::to_string(r + 1));
            if (r == 0)
                count = arrivals[r]->data.size();
            else if (arrivals[r]->data.size() != count)
                throw InputError(
                    "ParamEstimator.ubr: sampled metrics have different sample counts; call "
                    "interpolate first");
            utils[r] = get_util(node, r + 1);
            if (utils[r] != nullptr && utils[r]->data.size() != count)
                throw InputError(
                    "ParamEstimator.ubr: sampled metrics have different sample counts; call "
                    "interpolate first");
        }

        Matrix<double> estimates(1, classes, 0.0);
        std::vector<std::size_t> unknown;
        std::vector<double> known_busy(count, 0.0);
        for (std::size_t r = 0; r < classes; ++r) {
            if (utils[r] == nullptr) {
                unknown.push_back(r);
                continue;
            }
            double aa = 0.0, au = 0.0;
            for (std::size_t i = 0; i < count; ++i) {
                const double busy = utils[r]->data[i] * servers;
                aa += arrivals[r]->data[i] * arrivals[r]->data[i];
                au += arrivals[r]->data[i] * busy;
                known_busy[i] += busy;
            }
            if (!(aa > 0.0))
                throw NumericError("ParamEstimator.ubr: a class has zero arrival-rate regressor");
            estimates(0, r) = std::max(0.0, au / aa);
        }

        if (!unknown.empty()) {
            SampledMetric* aggregate = get_aggr_util(node);
            if (aggregate == nullptr)
                throw InputError(
                    "ParamEstimator.ubr: aggregate utilization is required for classes without "
                    "per-class utilization");
            if (aggregate->data.size() != count)
                throw InputError(
                    "ParamEstimator.ubr: sampled metrics have different sample counts; call "
                    "interpolate first");
            Matrix<double> A(count, unknown.size(), 0.0);
            std::vector<double> residual(count, 0.0);
            for (std::size_t i = 0; i < count; ++i) {
                residual[i] = aggregate->data[i] * servers - known_busy[i];
                for (std::size_t j = 0; j < unknown.size(); ++j)
                    A(i, j) = arrivals[unknown[j]]->data[i];
            }
            const std::vector<double> fit = detail::nnls(A, residual);
            for (std::size_t j = 0; j < unknown.size(); ++j) estimates(0, unknown[j]) = fit[j];
        }
        return estimates;
    }
};

}  // namespace infer
}  // namespace line

#endif  // LINE_INFERENCE_PARAM_ESTIMATOR_H
