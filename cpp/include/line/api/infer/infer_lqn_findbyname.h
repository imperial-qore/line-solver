/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_INFER_INFER_LQN_FINDBYNAME_H
#define LINE_API_INFER_INFER_LQN_FINDBYNAME_H

/**
 * First element of a named LQN container matching a name.
 *
 * Port of matlab/src/api/infer/infer_lqn_findbyname.m. No JAR counterpart.
 *
 * MATLAB walks a cell array of LQN element handles and returns the first whose
 * .name equals the argument, or [] when none does. The handle itself belongs
 * to the LayeredNetwork model layer, which this phase of the port does not
 * carry, so what crosses is the lookup: the container is the name list
 * (LayeredNetworkStruct.names, the same list infer_lqn_getobs indexes) and the
 * result is the POSITION of the match, which the caller uses to index whatever
 * parallel array holds its elements.
 *
 * The MATLAB is a linear scan returning the FIRST match, and duplicates are
 * not diagnosed; the port keeps both properties, since an LQN with two
 * elements of one name resolves consistently between the two implementations
 * only if they agree on which one wins.
 *
 * ARITHMETIC: none. String comparison only, so the function is not templated
 * on the number type and is usable from every instantiation.
 */

#include <cstddef>
#include <string>
#include <vector>

namespace line {
namespace infer {

/** Returned when no element carries the requested name, MATLAB's []. */
static const std::size_t INFER_LQN_NOT_FOUND = static_cast<std::size_t>(-1);

/**
 * @param names container of element names, in element order
 * @param name  name to locate
 * @return      index of the first match, or INFER_LQN_NOT_FOUND
 */
inline std::size_t infer_lqn_findbyname(const std::vector<std::string>& names,
                                        const std::string& name) {
    for (std::size_t k = 0; k < names.size(); ++k)
        if (names[k] == name) return k;
    return INFER_LQN_NOT_FOUND;
}

}  // namespace infer
}  // namespace line

#endif  // LINE_API_INFER_INFER_LQN_FINDBYNAME_H
