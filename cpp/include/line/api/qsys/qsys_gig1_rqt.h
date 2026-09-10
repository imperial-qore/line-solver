/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_QSYS_QSYS_GIG1_RQT_H
#define LINE_API_QSYS_QSYS_GIG1_RQT_H

/**
 * Robust Queueing Theory (RQT) worst-case system time of a G/G/1 FCFS queue,
 * the single-server case of qsys_gigk_rqt (Theorem 2 and eq. 12).
 *
 * Templated port of matlab/src/api/qsys/qsys_gig1_rqt.m, cross-checked against
 * jar/src/main/java/jline/api/qsys/Qsys_gig1_rqt.java.
 *
 * Reference: C. Bandi, D. Bertsimas, N. Youssef (2015). Robust Queueing Theory.
 * Operations Research 63(3), 676-700.
 */

#include "line/api/qsys/qsys_gigk_rqt.h"

namespace line {
namespace qsys {

/**
 * @param lambda  arrival rate
 * @param mu      service rate
 * @param Gamma_a variability parameter of the arrival uncertainty set
 * @param Gamma_s variability parameter of the service uncertainty set
 * @param alpha_a arrival tail coefficient in (1,2]
 * @param alpha_s service tail coefficient in (1,2]
 */
template <class T>
GigkRqtResult<T> qsys_gig1_rqt(const T& lambda, const T& mu, const T& Gamma_a, const T& Gamma_s,
                               const T& alpha_a, const T& alpha_s) {
    return qsys_gigk_rqt(lambda, mu, Gamma_a, Gamma_s, static_cast<std::size_t>(1), alpha_a, alpha_s);
}

/** Finite-variance case, alpha_a = alpha_s = 2. */
template <class T>
GigkRqtResult<T> qsys_gig1_rqt(const T& lambda, const T& mu, const T& Gamma_a, const T& Gamma_s) {
    const T two = num_traits<T>::from_int(2);
    return qsys_gigk_rqt(lambda, mu, Gamma_a, Gamma_s, static_cast<std::size_t>(1), two, two);
}

}  // namespace qsys
}  // namespace line

#endif  // LINE_API_QSYS_QSYS_GIG1_RQT_H
