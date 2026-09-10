/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_SN_SN_REGION_MEMBERS_H
#define LINE_API_SN_SN_REGION_MEMBERS_H

/**
 * Port of matlab/src/api/sn/sn_region_members.m.
 *
 * Which stations belong to finite-capacity region `f`, 0-based. Membership is
 * DECLARED (`sn.regionmembers`, here `Region::members`) and only inferred when
 * the declaration is absent: a station is then a member iff some entry of its
 * capacity row, or its memory cap, is not the -1 "unbounded" sentinel.
 *
 * The inference cannot see a member station left wholly unbounded, which is
 * precisely why the declaration exists and is preferred. A caller that infers
 * membership when the flags are present would silently shrink the region and
 * admit jobs the region forbids.
 *
 * ARITHMETIC: none. Structural.
 */

#include <cstddef>
#include <vector>

#include "line/lang/qn/network_struct.h"

namespace line {
namespace api {

/**
 * The membership mask of region `f` (0-based), one entry per station.
 *
 * Returns an empty mask when the model declares no such region, which is the
 * reference's behaviour of leaving `mask` empty rather than raising.
 */
template <class T>
std::vector<bool> sn_region_members(const qn::NetworkStruct<T>& sn, std::size_t f) {
    std::vector<bool> mask;
    if (f >= sn.regions.size()) return mask;
    const typename qn::NetworkStruct<T>::Region& rg = sn.regions[f];
    mask.assign(sn.nstations, false);
    bool declared = false;
    for (std::size_t i = 0; i < rg.members.size(); ++i)
        if (rg.members[i]) declared = true;
    for (std::size_t i = 0; i < sn.nstations; ++i) {
        if (declared) {
            mask[i] = i < rg.members.size() && rg.members[i];
            continue;
        }
        bool in = false;
        if (i < rg.cap.size())
            for (std::size_t r = 0; r < rg.cap[i].size(); ++r)
                if (rg.cap[i][r] != -1.0) in = true;
        if (i < rg.maxmem.size() && rg.maxmem[i] != -1.0) in = true;
        mask[i] = in;
    }
    return mask;
}

/** The same membership as 1-based station indices, the form most callers want. */
template <class T>
std::vector<std::size_t> sn_region_member_list(const qn::NetworkStruct<T>& sn, std::size_t f) {
    const std::vector<bool> mask = sn_region_members(sn, f);
    std::vector<std::size_t> out;
    for (std::size_t i = 0; i < mask.size(); ++i)
        if (mask[i]) out.push_back(i + 1);
    return out;
}

}  // namespace api
}  // namespace line

#endif  // LINE_API_SN_SN_REGION_MEMBERS_H
