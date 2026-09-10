/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_SOLVERS_AG_AG_DISPATCH_H
#define LINE_SOLVERS_AG_AG_DISPATCH_H

/**
 * @file ag_dispatch.h
 * @brief The `-s ag` entry point: gates, fixed point, mean measures.
 *
 * The AG twin of `solver_mam_run_analyzer`. It filters the same metric kinds through the
 * same mask rules, because what a solver reports is a property of the model and
 * not of the algorithm; everything above the filter -- which agents exist, how
 * they couple, who evaluates them -- is SolverAG's own.
 */

#include <string>
#include <vector>

#include "line/lang/lang_types.h"
#include "line/lang/qn/network_struct.h"
#include "line/solvers/ag/ag_types.h"
#include "line/solvers/ag/solver_ag_runner.h"
#include "line/solvers/mva/solver_mva_runner.h"
#include "line/util/error.h"

namespace line {
namespace ag {

/** `SolverAG.runAnalyzer`: the converged agents as mean measures. */
template <class T>
mva::AvgResult<T> solver_ag_run_analyzer(const qn::NetworkStruct<T>& L, const AgOptions& opt) {
    const std::string origmethod = opt.method;
    const AgResult<T> d = solver_ag_solve(L, opt);
    const mva::MvaSolution<T>& s = d.sol;

    const std::size_t M = L.nstations, K = L.nclasses;
    std::vector<std::vector<bool>> mask(M, std::vector<bool>(K, false));
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t k = 0; k < K; ++k)
            mask[i][k] = num_traits<T>::to_double(s.R(i, k)) < 10.0 * lang::GlobalConstants::FineTol;
    std::vector<std::vector<bool>> srcmask(M, std::vector<bool>(K, false));
    for (std::size_t i = 0; i < M; ++i)
        if (L.stations[i].nodetype == qn::NodeType::Source)
            for (std::size_t k = 0; k < K; ++k) srcmask[i][k] = true;

    mva::AvgResult<T> out;
    out.QN = mva::filter_metric(L, s.Q, mva::MetricKind::QLen, &mask);
    out.UN = mva::filter_metric(L, s.U, mva::MetricKind::Util, &mask);
    out.RN = mva::filter_metric(L, s.R, mva::MetricKind::RespT, nullptr);
    out.TN = mva::filter_metric(L, s.Tp, mva::MetricKind::Tput, nullptr);
    out.WN = mva::filter_metric(L, mva::sn_get_residt_from_respt(L, out.RN),
                                mva::MetricKind::ResidT, nullptr);
    out.AN = mva::filter_metric(L, mva::sn_get_arvr_from_tput(L, out.TN), mva::MetricKind::ArvR,
                                &srcmask);
    out.CN = s.C;
    out.XN = s.X;
    // TWO FIELDS, NOT ONE: `method` is what the caller asked for and
    // `actualmethod` is the algorithm that produced the numbers, which is the
    // convention `solver_nc_run_analyzer` keeps and the one every banner reads. Writing
    // the resolved name into `method` and leaving `actualmethod` EMPTY is what
    // made the example twin print `AG (method=)` with no name at all, where the
    // JAR records "inap" for the same solve.
    //
    // 'default' resolves to inap and 'exact' falls back to it, so the resolved
    // name is reported rather than the asked-for one: 'exact' is classified
    // globally as an exact method, and leaving the name in place would banner an
    // iterative approximation as exact.
    out.method = origmethod;
    out.actualmethod = (origmethod == "default" && !d.actualmethod.empty() &&
                        d.actualmethod != "default")
                           ? "default/" + d.actualmethod
                           : (origmethod == "exact" ? d.actualmethod : origmethod);
    out.iter = s.iter;
    return out;
}

}  // namespace ag
}  // namespace line

#endif  // LINE_SOLVERS_AG_AG_DISPATCH_H
