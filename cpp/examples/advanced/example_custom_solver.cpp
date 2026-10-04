/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * `matlab/examples/advanced/example_custom_solver/`,
 * `python/examples/advanced/example_custom_solver/`: the template for a
 * user-written LINE solver.
 *
 * NOTHING HERE COMPUTES ANYTHING, AND THAT IS THE EXAMPLE. The reference ships
 * a solver skeleton whose algorithm is deliberately empty: it reads the
 * NetworkStruct, returns zeros and names the file that has to be filled in.
 * What it teaches is the SHAPE a solver has -- where the model struct is read,
 * where the algorithm goes, and what the wrapper around it owes the caller.
 *
 * THE REFERENCE SPREADS THAT SHAPE OVER FOUR FILES, and this port keeps the
 * same three pieces as three functions rather than three registrations, because
 * only one of them is a thing to run:
 *
 *   `solver_custom`           the bare algorithm, the piece to fill in
 *   `solver_custom_analyzer`  what happens around the call: timing, cleanup
 *   `SolverCustom`            the solver object, which reads the struct and
 *                             publishes the results
 *
 * MATLAB's fourth file is the class-folder split of the third (`@SolverCustom/`
 * holds the constructor and `runAnalyzer` in separate files), an idiom neither
 * C++ nor Python has; Python already merges the two into `SolverCustom.py` for
 * the same reason. The registered example is what the Python twin's own entry
 * point runs: build a small closed model and drive it through the template.
 *
 * THE ZERO TABLE IS PRINTED RATHER THAN SUPPRESSED. `print_avg` drops a row
 * whose six metrics are all zero, which is right for a real solve and wrong
 * here, since the zeros ARE the result this template returns; the matrices go
 * out through `print_matrix` so the reader sees the shape the analyzer filled.
 */

#include <cmath>
#include <cstdio>
#include <ctime>
#include <string>
#include <vector>

#include "example_util.h"
#include "examples_common.h"

namespace line {
namespace examples {

namespace {

/**
 * The three template pieces keep the reference's own names, which is why they
 * sit in a namespace of their own: the runnable example must be registered as
 * `solver_custom`, and a file-scope `solver_custom(const Sn&)` beside it would
 * make that registration an ambiguous overload rather than a name.
 */
namespace custom {

/** The six average metrics `solver_custom` returns, in the reference's order. */
struct CustomResult {
    Matrix<double> QN, UN, RN, TN;
    std::vector<double> CN, XN;
    double runtime = 0.0;
};

/**
 * `solver_custom(sn, options)`: THE PIECE TO FILL IN.
 *
 * It receives the NetworkStruct and must return the six average-metric
 * matrices. The template reads the quantities an algorithm starts from --
 * station and class counts, populations, rates, visits -- returns zeros and
 * says so, exactly as the MATLAB and Python twins do.
 */
CustomResult solver_custom(const Sn& sn) {
    const std::size_t M = sn.nstations;          // number of stations
    const std::size_t K = sn.nclasses;           // number of classes
    const std::vector<double> N = sn.njobs();    // job populations
    const Matrix<double>& rates = sn.rates;      // arrival and service rates
    const std::vector<Matrix<double> >& V = sn.visits;  // visits, per chain
    (void)N;
    (void)rates;
    (void)V;

    CustomResult r;
    r.QN = Matrix<double>(M, K);
    r.UN = Matrix<double>(M, K);
    r.RN = Matrix<double>(M, K);
    r.TN = Matrix<double>(M, K);
    r.CN.assign(K, 0.0);
    r.XN.assign(K, 0.0);

    note("The solution algorithm needs to be implemented in solver_custom: "
         "returning with no result.");
    return r;
}

/**
 * `solver_custom_analyzer(sn, options)`: everything around the algorithm.
 *
 * Any activity prior to or after launching the solution algorithm belongs here.
 * The template times the call and scrubs the NaNs the algorithm may have left
 * behind, which is the one piece of cleanup every LINE analyzer owes its caller:
 * a metric a solver could not compute reaches the tables as zero, not as NaN.
 */
CustomResult solver_custom_analyzer(const Sn& sn) {
    const std::clock_t t0 = std::clock();
    note("Any activity prior or after launching the solution algorithm needs to be "
         "implemented in solver_custom_analyzer.");
    CustomResult r = solver_custom(sn);

    Matrix<double>* mats[4] = {&r.QN, &r.UN, &r.RN, &r.TN};
    for (int k = 0; k < 4; ++k)
        for (std::size_t i = 0; i < mats[k]->rows(); ++i)
            for (std::size_t j = 0; j < mats[k]->cols(); ++j)
                if (std::isnan((*mats[k])(i, j))) (*mats[k])(i, j) = 0.0;
    std::vector<double>* vecs[2] = {&r.CN, &r.XN};
    for (int k = 0; k < 2; ++k)
        for (std::size_t j = 0; j < vecs[k]->size(); ++j)
            if (std::isnan((*vecs[k])[j])) (*vecs[k])[j] = 0.0;

    r.runtime = double(std::clock() - t0) / double(CLOCKS_PER_SEC);
    return r;
}

/**
 * `SolverCustom`: the solver object.
 *
 * It reads the struct, calls the analyzer and publishes the results. To write a
 * real solver, fill in `solver_custom` and then narrow what the solver CLAIMS:
 * the reference's `getFeatureSet` is deliberately broad, and a featset wider
 * than the algorithm makes `supports` accept a model the algorithm then gets
 * wrong. A refusal by name is the correct answer for a feature the algorithm
 * cannot analyze; a table of zeros is not.
 */
class SolverCustom {
  public:
    explicit SolverCustom(Net& m) : model_(m) {}

    /** The data structure summarizing the model; needs no initial state. */
    const Sn& get_struct() const { return model_.get_struct(); }

    /** The method names this solver accepts, `listValidMethods`. */
    static std::vector<std::string> list_valid_methods() {
        return std::vector<std::string>(1, "default");
    }

    CustomResult run_analyzer() { return solver_custom_analyzer(get_struct()); }

  private:
    Net& model_;
};

}  // namespace custom

}  // namespace

/** The template driven over a small closed model, as the Python twin's `__main__` does. */
void solver_custom() {
    banner("customExample");

    Net m("customExample");
    Delay think(m, "Think");
    Queue queue(m, "Queue1", SchedStrategy::PS);
    ClosedClass jobclass(m, "Class1", 2.0, think);
    think.set_service(jobclass, Exp(1.0));
    queue.set_service(jobclass, Exp(2.0));
    link_serial(m, {think, queue});

    custom::SolverCustom solver(m);
    section("Custom");
    const custom::CustomResult r = solver.run_analyzer();

    print_matrix("QN", r.QN);
    print_matrix("UN", r.UN);
    print_matrix("RN", r.RN);
    print_matrix("TN", r.TN);
    print_vector("CN", r.CN);
    print_vector("XN", r.XN);
    std::printf("runtime: %.6f s\n", r.runtime);
}

LINE_EXAMPLE("advanced/example_custom_solver", solver_custom);

}  // namespace examples
}  // namespace line
