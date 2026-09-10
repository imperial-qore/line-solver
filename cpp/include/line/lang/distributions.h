/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_LANG_DISTRIBUTIONS_H
#define LINE_LANG_DISTRIBUTIONS_H

/**
 * The distributions a user names, spelled as their Python twins.
 *
 *     queue1.set_service(closed, Exp(1.0));
 *     delay.set_service(open, Erlang(3.0, 2));
 *     delay.set_service(open, HyperExp(0.5, 3.0, 10.0));
 *
 * against Python's `Exp(1.0)`, `Erlang(3.0, 2)`, `HyperExp(0.5, 3.0, 10.0)`.
 *
 * EACH IS A THIN NAME OVER THE FACTORY THAT ALREADY EXISTS. Every constructor
 * assigns the `Distrib<double>` that `lang_types.h` builds, so the numeric
 * behaviour is the factory's and nothing new is computed here; deriving from
 * `Distrib<double>` is what lets one slice cleanly into `set_service`, which
 * takes the base by const reference.
 *
 * The fitters keep Python's static-method spelling
 * (`HyperExp::fit_mean_and_scv`), forwarding to `line/lang/dist_fitters.h`.
 *
 * `double`-only, as Python is. Reach for `lang::Distrib<T>`'s own factories
 * directly for the multiprecision paths.
 */

#include <cstddef>
#include <string>
#include <vector>

#include "line/lang/dist_fitters.h"
#include "line/lang/lang_types.h"

namespace line {

typedef lang::Distrib<double> Distribution;

/** Every named law below is a `Distrib<double>` under another name. */
#define LINE_DIST_CTOR(Cls) \
    struct Cls : lang::Distrib<double>

/** `Exp(rate)`: the exponential law of the given RATE, as Python's `Exp`. */
LINE_DIST_CTOR(Exp) {
    explicit Exp(double rate) { static_cast<lang::Distrib<double>&>(*this) = lang::Distrib<double>::exp_rate(rate); }
    /** `Exp.fitMean(mean)`. */
    static Exp fit_mean(double m) { return Exp(1.0 / m); }
    /** `Exp.fitRate(rate)`. */
    static Exp fit_rate(double r) { return Exp(r); }
};

/** `Erlang(phase_rate, nphases)`. */
LINE_DIST_CTOR(Erlang) {
    Erlang(double phase_rate, std::size_t nphases) {
        static_cast<lang::Distrib<double>&>(*this) = lang::Distrib<double>::erlang(phase_rate, nphases);
    }
    /** `Erlang.fitMeanAndSCV(mean, scv)`. */
    static lang::Distrib<double> fit_mean_and_scv(double m, double scv) {
        return lang::Distrib<double>::erlang_fit(m, scv);
    }
    /** `Erlang.fitMeanAndOrder(mean, k)`. */
    static lang::Distrib<double> fit_mean_and_order(double m, std::size_t k) {
        return lang::erlang_fit_mean_order<double>(m, k);
    }
};

/** `HyperExp(p, lambda1, lambda2)`. */
LINE_DIST_CTOR(HyperExp) {
    HyperExp(double p, double lambda1, double lambda2) {
        static_cast<lang::Distrib<double>&>(*this) = lang::Distrib<double>::hyperexp(p, lambda1, lambda2);
    }
    HyperExp(const std::vector<double>& p, const std::vector<double>& lambda) {
        static_cast<lang::Distrib<double>&>(*this) = lang::Distrib<double>::hyperexp_n(p, lambda);
    }
    /** `HyperExp.fitMeanAndSCV(mean, scv)`. */
    static lang::Distrib<double> fit_mean_and_scv(double m, double scv) {
        return lang::hyperexp_fit_mean_scv<double>(m, scv);
    }
    /** `HyperExp.fitMeanAndSCVBalanced(mean, scv)`. */
    static lang::Distrib<double> fit_mean_and_scv_balanced(double m, double scv) {
        return lang::hyperexp_fit_mean_scv_balanced<double>(m, scv);
    }
};

/** `Coxian(mu, phi)`. */
LINE_DIST_CTOR(Coxian) {
    Coxian(const std::vector<double>& mu, const std::vector<double>& phi) {
        static_cast<lang::Distrib<double>&>(*this) = lang::Distrib<double>::coxian(mu, phi);
    }
    /** `Coxian.fitMeanAndSCV(mean, scv)`. */
    static lang::Distrib<double> fit_mean_and_scv(double m, double scv) {
        return lang::coxian_fit_mean_scv<double>(m, scv);
    }
    /** `Coxian.fitCentral(mean, scv, skew)`. */
    static lang::Distrib<double> fit_central(double m, double scv, double skew) {
        return lang::coxian_fit_central<double>(m, scv, skew);
    }
};

/** `Cox2(mu1, mu2, phi1)`. */
LINE_DIST_CTOR(Cox2) {
    Cox2(double mu1, double mu2, double phi1) {
        static_cast<lang::Distrib<double>&>(*this) = lang::Distrib<double>::cox2(mu1, mu2, phi1);
    }
    /** `Cox2.fitCentral(mean, scv, skew)`. */
    static lang::Distrib<double> fit_central(double m, double scv, double skew) {
        return lang::cox2_fit_central<double>(m, scv, skew);
    }
};

/** `PH(alpha, A)`: a general phase-type law. */
LINE_DIST_CTOR(PH) {
    PH(const std::vector<double>& alpha, const Matrix<double>& A) {
        static_cast<lang::Distrib<double>&>(*this) = lang::Distrib<double>::phase_type(alpha, A, false);
    }
};

/** `APH(alpha, A)`: the acyclic phase-type law. */
LINE_DIST_CTOR(APH) {
    APH(const std::vector<double>& alpha, const Matrix<double>& A) {
        static_cast<lang::Distrib<double>&>(*this) = lang::Distrib<double>::phase_type(alpha, A, true);
    }
    /** `APH.fitMeanAndSCV(mean, scv)`. */
    static lang::Distrib<double> fit_mean_and_scv(double m, double scv) {
        return lang::aph_fit_mean_scv<double>(m, scv);
    }
    /** `APH.fitCentral(mean, scv, skew)`. */
    static lang::Distrib<double> fit_central(double m, double scv, double skew) {
        return lang::aph_fit_central<double>(m, scv, skew);
    }
};

/** `MAP(D0, D1)`. */
LINE_DIST_CTOR(MAP) {
    MAP(const Matrix<double>& D0, const Matrix<double>& D1) {
        static_cast<lang::Distrib<double>&>(*this) = lang::Distrib<double>::map_dist(D0, D1, lang::ProcessType::MAP);
    }
};

/** `Det(t)`: the deterministic law. */
LINE_DIST_CTOR(Det) {
    explicit Det(double t) { static_cast<lang::Distrib<double>&>(*this) = lang::Distrib<double>::det(t); }
};

/** `Immediate()`: a zero-time transition. */
LINE_DIST_CTOR(Immediate) {
    Immediate() { static_cast<lang::Distrib<double>&>(*this) = lang::Distrib<double>::immediate(); }
};

/** `Disabled()`: the class is not served here. */
LINE_DIST_CTOR(Disabled) {
    Disabled() { static_cast<lang::Distrib<double>&>(*this) = lang::Distrib<double>::disabled_dist(); }
};

/** `Uniform(a, b)`. */
LINE_DIST_CTOR(Uniform) {
    Uniform(double a, double b) { static_cast<lang::Distrib<double>&>(*this) = lang::Distrib<double>::uniform(a, b); }
};

/** `Pareto(shape, scale)`. */
LINE_DIST_CTOR(Pareto) {
    Pareto(double shape, double scale) { static_cast<lang::Distrib<double>&>(*this) = lang::Distrib<double>::pareto(shape, scale); }
    /** `Pareto.fitMeanAndSCV(mean, scv)`. */
    static lang::Distrib<double> fit_mean_and_scv(double m, double scv) {
        return lang::pareto_fit_mean_scv<double>(m, scv);
    }
};

/** `Gamma(shape, scale)`. */
LINE_DIST_CTOR(Gamma) {
    Gamma(double shape, double scale) { static_cast<lang::Distrib<double>&>(*this) = lang::Distrib<double>::gamma_dist(shape, scale); }
    /** `Gamma.fitMeanAndSCV(mean, scv)`. */
    static lang::Distrib<double> fit_mean_and_scv(double m, double scv) {
        return lang::gamma_fit_mean_scv<double>(m, scv);
    }
};

/** `Weibull(scale, shape)`. */
LINE_DIST_CTOR(Weibull) {
    Weibull(double scale, double shape) { static_cast<lang::Distrib<double>&>(*this) = lang::Distrib<double>::weibull(scale, shape); }
};

/** `Lognormal(logmean, logsigma)`. */
LINE_DIST_CTOR(Lognormal) {
    Lognormal(double logmean, double logsigma) {
        static_cast<lang::Distrib<double>&>(*this) = lang::Distrib<double>::lognormal(logmean, logsigma);
    }
};

/** `Normal(mu, sigma)`. */
LINE_DIST_CTOR(Normal) {
    Normal(double mu, double sigma) { static_cast<lang::Distrib<double>&>(*this) = lang::Distrib<double>::normal(mu, sigma); }
};

/** `Geometric(p)`. */
LINE_DIST_CTOR(Geometric) {
    explicit Geometric(double p) { static_cast<lang::Distrib<double>&>(*this) = lang::Distrib<double>::geometric(p); }
};

/** `Bernoulli(p)`. */
LINE_DIST_CTOR(Bernoulli) {
    explicit Bernoulli(double p) { static_cast<lang::Distrib<double>&>(*this) = lang::Distrib<double>::bernoulli(p); }
};

/** `Binomial(n, p)`. */
LINE_DIST_CTOR(Binomial) {
    Binomial(double n, double p) { static_cast<lang::Distrib<double>&>(*this) = lang::Distrib<double>::binomial(n, p); }
};

/** `Poisson(lambda)`. */
LINE_DIST_CTOR(Poisson) {
    explicit Poisson(double lambda) { static_cast<lang::Distrib<double>&>(*this) = lang::Distrib<double>::poisson(lambda); }
};

/** `DiscreteUniform(a, b)`. */
LINE_DIST_CTOR(DiscreteUniform) {
    DiscreteUniform(double a, double b) {
        static_cast<lang::Distrib<double>&>(*this) = lang::Distrib<double>::discrete_uniform(a, b);
    }
};

/** `Zipf(s, n)`. */
LINE_DIST_CTOR(Zipf) {
    Zipf(double s, std::size_t n) { static_cast<lang::Distrib<double>&>(*this) = lang::Distrib<double>::zipf(s, n); }
};

/** `DiscreteSampler(p, x)`. */
LINE_DIST_CTOR(DiscreteSampler) {
    explicit DiscreteSampler(const std::vector<double>& p) {
        static_cast<lang::Distrib<double>&>(*this) = lang::Distrib<double>::discrete_sampler(p, std::vector<double>());
    }
    DiscreteSampler(const std::vector<double>& p, const std::vector<double>& x) {
        static_cast<lang::Distrib<double>&>(*this) = lang::Distrib<double>::discrete_sampler(p, x);
    }
};

/** `Replayer(samples)` / `Replayer(samples, path)`: the trace-driven law. */
LINE_DIST_CTOR(Replayer) {
    explicit Replayer(const std::vector<double>& samples) {
        static_cast<lang::Distrib<double>&>(*this) = lang::Distrib<double>::replayer(samples);
    }
    Replayer(const std::vector<double>& samples, const std::string& path) {
        static_cast<lang::Distrib<double>&>(*this) = lang::Distrib<double>::replayer_from(samples, path);
    }
};

#undef LINE_DIST_CTOR

}  // namespace line

#endif  // LINE_LANG_DISTRIBUTIONS_H
