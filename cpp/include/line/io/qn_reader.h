/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_IO_QN_READER_H
#define LINE_IO_QN_READER_H

/**
 * Reader for the .qn closed-network model format of mp_pfqn.
 *
 * Format (mp_pfqn/util/readmodel.c):
 *   R                    number of classes
 *   N1 ... NR            population per class, -1 marks an open class
 *   Z1 ... ZR            think times, integers
 *   M                    number of queueing stations
 *   mi Li1 ... LiR       per station: multiplicity then demands, integers
 * followed by optional keyword sections:
 *   LAMBDA               arrival rates for the open classes, rationals allowed
 *   MU                   load-dependent rates, M x Nt, rationals allowed
 *
 * Demands and think times are integers on purpose: that is what keeps the
 * exact rational solvers cheap, and it is why this format is the interchange
 * used to diff against mp_pfqn's binaries. Note the trap documented in
 * _kb/14-cpp-multiprecision: a rational demand like 1/10 in the L section
 * parses as 1 and then desynchronizes the scan, in mp_pfqn as well as here.
 */

#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace io {

template <class T>
struct QnModel {
    int R = 0;                 ///< classes
    int M = 0;                 ///< queueing stations
    std::vector<int> N;        ///< population per class, -1 for an open class
    Matrix<T> Z;               ///< (1 x R) think times
    Matrix<T> L;               ///< (M x R) demands
    std::vector<int> mi;       ///< (M) multiplicities
    bool hasOpen = false;      ///< LAMBDA section present
    std::vector<T> lambda;     ///< (R) arrival rates, empty when absent
    bool isLD = false;         ///< MU section present
    Matrix<T> mu;              ///< (M x Nt) load-dependent rates
    int Nt = 0;                ///< total closed population

    /** True when the model is a plain closed network the exact solvers accept. */
    bool isClosedLoadIndependent() const { return !hasOpen && !isLD; }
};

/** Parse a rational token of the form "3" or "3/10". */
template <class T>
T parse_rational_token(const std::string& tok) {
    const std::size_t slash = tok.find('/');
    if (slash == std::string::npos) return num_traits<T>::from_int(std::stol(tok));
    const long num = std::stol(tok.substr(0, slash));
    const long den = std::stol(tok.substr(slash + 1));
    if (den == 0) throw InputError("qn_reader: zero denominator in '" + tok + "'");
    return num_traits<T>::from_rational(num, den);
}

template <class T>
QnModel<T> read_qn(const std::string& path) {
    std::ifstream f(path);
    if (!f) throw InputError("qn_reader: cannot open " + path);

    QnModel<T> m;
    if (!(f >> m.R) || m.R <= 0) throw InputError("qn_reader: bad class count in " + path);
    m.N.resize(m.R);
    for (int r = 0; r < m.R; ++r)
        if (!(f >> m.N[r])) throw InputError("qn_reader: truncated population vector");
    m.Z = Matrix<T>(1, m.R);
    for (int r = 0; r < m.R; ++r) {
        long z;
        if (!(f >> z)) throw InputError("qn_reader: truncated think-time vector");
        m.Z(0, r) = num_traits<T>::from_int(z);
    }
    if (!(f >> m.M) || m.M < 0) throw InputError("qn_reader: bad station count");
    m.L = Matrix<T>(m.M, m.R);
    m.mi.assign(m.M, 1);
    for (int i = 0; i < m.M; ++i) {
        if (!(f >> m.mi[i])) throw InputError("qn_reader: truncated station line");
        for (int r = 0; r < m.R; ++r) {
            long d;
            if (!(f >> d)) throw InputError("qn_reader: truncated demand row");
            m.L(i, r) = num_traits<T>::from_int(d);
        }
    }
    m.Nt = 0;
    for (int r = 0; r < m.R; ++r)
        if (m.N[r] > 0) m.Nt += m.N[r];

    std::string keyword;
    while (f >> keyword) {
        if (keyword == "LAMBDA") {
            m.hasOpen = true;
            m.lambda.resize(m.R);
            for (int r = 0; r < m.R; ++r) {
                std::string tok;
                if (!(f >> tok)) throw InputError("qn_reader: truncated LAMBDA section");
                m.lambda[r] = parse_rational_token<T>(tok);
            }
        } else if (keyword == "MU") {
            m.isLD = true;
            if (m.Nt <= 0) throw InputError("qn_reader: MU section with an empty closed population");
            m.mu = Matrix<T>(m.M, m.Nt);
            for (int i = 0; i < m.M; ++i)
                for (int k = 0; k < m.Nt; ++k) {
                    std::string tok;
                    if (!(f >> tok)) throw InputError("qn_reader: truncated MU section");
                    m.mu(i, k) = parse_rational_token<T>(tok);
                }
        }
        // Unknown keywords are skipped, as mp_pfqn's reader does.
    }
    return m;
}

}  // namespace io
}  // namespace line

#endif  // LINE_IO_QN_READER_H
