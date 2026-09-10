/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MDD_MDD_H
#define LINE_API_MDD_MDD_H

/**
 * Quasi-reduced ordered Multi-valued Decision Diagram.
 *
 * Port of matlab/src/api/mdd/MDD.m, jline.api.mdd.MDD and
 * python/line_solver/api/mdd/mdd.py, after A.S. Miner, G. Ciardo, "Efficient
 * Reachability Set Generation and Storage Using Decision Diagrams", ICATPN
 * 1999, LNCS 1639, pp.6-25.
 *
 * A global state is a K-tuple of LOCAL state values, one per level/submodel,
 * state[k] in {0,...,domain[k]-1}. The set is stored as a directed acyclic
 * graph with K variable levels plus a terminal level: level 0 is the top
 * (root), a node at level k has domain[k] outgoing arcs to level k+1 nodes, and
 * a state belongs to the set iff its path of arcs reaches the TRUE terminal.
 * Canonicity is enforced by a per-level unique table (no duplicate nodes) and
 * by collapsing the all-FALSE node to the FALSE terminal. Storage is
 * O(#nodes), typically O(K * #local-states), instead of O(|S|).
 *
 * Node ids are 1-based in every codebase so that 0 can serve as TERM_FALSE, and
 * the unique table is keyed by the arc row itself: the encodings differ across
 * codebases (a byte-packed string in MATLAB, a tuple key in python) but the
 * canonical form does not.
 *
 * The diagram is pure combinatorics, so unlike the rest of the api layer it is
 * NOT templated on the numeric type: no rate ever enters it.
 */

#include <algorithm>
#include <cstddef>
#include <map>
#include <string>
#include <vector>

#include "line/util/error.h"

namespace line {
namespace mdd {

/** Terminal node "1": a completed path is accepted. */
const int TERM_TRUE = -1;
/** Terminal node "0": empty subgraph. */
const int TERM_FALSE = 0;

/**
 * Plain-array export of an MDD, the input contract of `mdd_mcd`.
 *
 * Mirrors MDD.toStruct in MATLAB/JAR and MDD.to_struct in python.
 */
struct MddStruct {
    /** Number of variable levels. */
    std::size_t K = 0;
    /** domain[k] is the number of local states at level k. */
    std::vector<int> domain;
    /** Id of the top (level-0) node; TERM_FALSE for the empty set. */
    int root = TERM_FALSE;
    /** nnodes[k] is the live node count at level k. */
    std::vector<int> nnodes;
    /**
     * node[k][p][v] is the child of arc v of level-k node id p+1: a level-(k+1)
     * node id when k < K-1, or a terminal when k == K-1.
     */
    std::vector<std::vector<std::vector<int>>> node;
};

/** Storage description of the set held in an MDD. */
struct MddStats {
    std::size_t levels = 0;
    /** Reachable node count per level. */
    std::vector<int> nodes_per_level;
    /** Reachable non-terminal nodes. */
    int num_nodes = 0;
    /** Nodes physically held in the tables, dead ones included. */
    int table_nodes = 0;
    /** |S|. */
    long long num_states = 0;
    /** Integers in the reachable arc arrays, the diagram footprint. */
    long long mdd_ints = 0;
    /** Integers an explicit state list would need, |S| * K. */
    long long explicit_ints = 0;

    /** Explicit footprint divided by the diagram footprint. */
    double compression() const {
        const long long den = mdd_ints > 1 ? mdd_ints : 1;
        return static_cast<double>(explicit_ints) / static_cast<double>(den);
    }
};

/** The diagram: insert / member / index / enumerate / cardinality. */
class MDD {
public:
    /** An empty set over the given per-level domains. */
    explicit MDD(const std::vector<int>& domain)
        : domain_(domain), K_(domain.size()), root_(TERM_FALSE), dirty_(true) {
        for (std::size_t k = 0; k < K_; ++k) {
            if (domain[k] <= 0) throw InputError("MDD: a level domain must be positive");
        }
        node_.assign(K_, std::vector<std::vector<int>>());
        uniq_.assign(K_, std::map<std::vector<int>, int>());
        cnt_.assign(K_, std::vector<long long>());
    }

    /** Build from a set of 0-based state tuples. */
    static MDD from_states(const std::vector<int>& domain,
                           const std::vector<std::vector<int>>& states) {
        MDD obj(domain);
        for (std::size_t i = 0; i < states.size(); ++i) obj.insert(states[i]);
        return obj;
    }

    std::size_t K() const { return K_; }
    const std::vector<int>& domain() const { return domain_; }
    /** Id of the root node, or TERM_FALSE for the empty set. */
    int root() const { return root_; }
    /** Number of live nodes at level k. */
    int node_count(std::size_t k) const { return static_cast<int>(node_[k].size()); }
    /** Arc row of level-k node id (1-based id). */
    const std::vector<int>& arcs(std::size_t k, int id) const { return node_[k][id - 1]; }

    /** Add a K-tuple of 0-based local values to the set. */
    void insert(const std::vector<int>& state) {
        if (state.size() != K_) throw InputError("MDD::insert: state has the wrong length");
        root_ = add_state(0, root_, state);
        dirty_ = true;
    }

    /** True iff state is in the set; O(K). */
    bool member(const std::vector<int>& state) const {
        int id = root_;
        for (std::size_t k = 0; k < K_; ++k) {
            if (id == TERM_FALSE) return false;
            id = node_[k][id - 1][state[k]];
        }
        return id == TERM_TRUE;
    }

    /** |S|, the number of stored states. */
    long long cardinality() const {
        ensure_counts();
        return child_count(0, root_);
    }

    /**
     * 0-based lexicographic rank of state among the stored set, level 0 most
     * significant, or -1 when the state is not stored.
     *
     * This is a bijection S -> {0,...,|S|-1}, so a generator matrix can be
     * assembled without an explicit state list.
     */
    long long index(const std::vector<int>& state) const {
        ensure_counts();
        long long idx = 0;
        int id = root_;
        for (std::size_t k = 0; k < K_; ++k) {
            if (id == TERM_FALSE) return -1;
            const std::vector<int>& a = node_[k][id - 1];
            const int v = state[k];
            for (int vv = 0; vv < v; ++vv) idx += child_count(k + 1, a[vv]);
            id = a[v];
        }
        return id == TERM_TRUE ? idx : -1;
    }

    /** All stored states as rows, in index() order. */
    std::vector<std::vector<int>> enumerate() const {
        std::vector<std::vector<int>> out;
        if (root_ == TERM_FALSE) return out;
        std::vector<int> prefix(K_, 0);
        enum_below(0, root_, prefix, out);
        return out;
    }

    /** Export the diagram as plain arrays for downstream algorithms. */
    MddStruct to_struct() const {
        MddStruct s;
        s.K = K_;
        s.domain = domain_;
        s.root = root_;
        s.nnodes.assign(K_, 0);
        s.node.assign(K_, std::vector<std::vector<int>>());
        for (std::size_t k = 0; k < K_; ++k) {
            s.nnodes[k] = static_cast<int>(node_[k].size());
            s.node[k] = node_[k];
        }
        return s;
    }

    /**
     * Reclaim dead nodes left by the append-only build.
     *
     * Membership, index and enumerate are unchanged. A production MDD would
     * reference-count instead and never accumulate dead nodes; this is the
     * basic sweep.
     */
    void compact() {
        const std::vector<std::vector<bool>> vis = reachable_ids();
        std::vector<std::vector<std::vector<int>>> newnode(K_);
        std::vector<std::vector<int>> remap(K_);
        for (std::size_t k = 0; k < K_; ++k) {
            remap[k].assign(node_[k].size() + 1, 0);
            for (std::size_t p = 0; p < node_[k].size(); ++p) {
                if (vis[k][p]) {
                    newnode[k].push_back(node_[k][p]);
                    remap[k][p + 1] = static_cast<int>(newnode[k].size());
                }
            }
        }
        for (std::size_t k = 0; k + 1 < K_; ++k) {
            for (std::size_t p = 0; p < newnode[k].size(); ++p) {
                for (int v = 0; v < domain_[k]; ++v) {
                    if (newnode[k][p][v] > 0) newnode[k][p][v] = remap[k + 1][newnode[k][p][v]];
                }
            }
        }
        node_ = newnode;
        if (root_ != TERM_FALSE) root_ = remap[0][root_];
        for (std::size_t k = 0; k < K_; ++k) {
            uniq_[k].clear();
            for (std::size_t p = 0; p < node_[k].size(); ++p)
                uniq_[k][node_[k][p]] = static_cast<int>(p + 1);
        }
        dirty_ = true;
    }

    /** Storage description of the current set; only reachable nodes are counted. */
    MddStats stats() const {
        const std::vector<std::vector<bool>> vis = reachable_ids();
        MddStats s;
        s.levels = K_;
        s.nodes_per_level.assign(K_, 0);
        for (std::size_t k = 0; k < K_; ++k) {
            int c = 0;
            for (std::size_t p = 0; p < vis[k].size(); ++p)
                if (vis[k][p]) ++c;
            s.nodes_per_level[k] = c;
            s.num_nodes += c;
            s.mdd_ints += static_cast<long long>(c) * domain_[k];
            s.table_nodes += static_cast<int>(node_[k].size());
        }
        s.num_states = cardinality();
        s.explicit_ints = s.num_states * static_cast<long long>(K_);
        return s;
    }

private:
    /** Canonical node creation through the per-level unique table. */
    int make_node(std::size_t k, const std::vector<int>& arc_row) {
        bool all_false = true;
        for (std::size_t v = 0; v < arc_row.size(); ++v)
            if (arc_row[v] != TERM_FALSE) {
                all_false = false;
                break;
            }
        if (all_false) return TERM_FALSE;  // collapse the empty node
        const std::map<std::vector<int>, int>::const_iterator it = uniq_[k].find(arc_row);
        if (it != uniq_[k].end()) return it->second;
        node_[k].push_back(arc_row);
        const int id = static_cast<int>(node_[k].size());
        uniq_[k][arc_row] = id;
        return id;
    }

    /**
     * Recursively add one state below node id at level k.
     *
     * Nodes are immutable and shared, so this rebuilds the path bottom-up
     * rather than mutating in place.
     */
    int add_state(std::size_t k, int id, const std::vector<int>& state) {
        if (k >= K_) return TERM_TRUE;
        std::vector<int> arc_row;
        if (id == TERM_FALSE)
            arc_row.assign(domain_[k], TERM_FALSE);
        else
            arc_row = node_[k][id - 1];
        const int v = state[k];
        if (v < 0 || v >= domain_[k])
            throw InputError("MDD::insert: a local value is outside its level domain");
        arc_row[v] = add_state(k + 1, arc_row[v], state);
        return make_node(k, arc_row);
    }

    long long child_count(std::size_t k, int child_id) const {
        if (k >= K_) return child_id == TERM_TRUE ? 1 : 0;
        if (child_id == TERM_FALSE) return 0;
        return count_node(k, child_id);
    }

    long long count_node(std::size_t k, int id) const {
        long long c = cnt_[k][id - 1];
        if (c >= 0) return c;
        const std::vector<int>& a = node_[k][id - 1];
        c = 0;
        for (int v = 0; v < domain_[k]; ++v) c += child_count(k + 1, a[v]);
        cnt_[k][id - 1] = c;
        return c;
    }

    void ensure_counts() const {
        bool sized = !dirty_;
        if (sized) {
            for (std::size_t k = 0; k < K_; ++k)
                if (cnt_[k].size() != node_[k].size()) {
                    sized = false;
                    break;
                }
        }
        if (sized) return;
        for (std::size_t k = 0; k < K_; ++k) cnt_[k].assign(node_[k].size(), -1);
        dirty_ = false;
    }

    /** Per-level masks of nodes reachable from the root. */
    std::vector<std::vector<bool>> reachable_ids() const {
        std::vector<std::vector<bool>> vis(K_);
        for (std::size_t k = 0; k < K_; ++k) vis[k].assign(node_[k].size(), false);
        if (root_ == TERM_FALSE) return vis;
        vis[0][root_ - 1] = true;
        std::vector<std::pair<std::size_t, int>> stack;
        stack.push_back(std::make_pair(static_cast<std::size_t>(0), root_));
        while (!stack.empty()) {
            const std::pair<std::size_t, int> top = stack.back();
            stack.pop_back();
            const std::size_t k = top.first;
            if (k + 1 == K_) continue;  // children are terminals
            const std::vector<int>& a = node_[k][top.second - 1];
            for (int v = 0; v < domain_[k]; ++v) {
                const int ch = a[v];
                if (ch > 0 && !vis[k + 1][ch - 1]) {
                    vis[k + 1][ch - 1] = true;
                    stack.push_back(std::make_pair(k + 1, ch));
                }
            }
        }
        return vis;
    }

    void enum_below(std::size_t k, int id, std::vector<int>& prefix,
                    std::vector<std::vector<int>>& out) const {
        const std::vector<int>& a = node_[k][id - 1];
        if (k + 1 == K_) {
            for (int v = 0; v < domain_[k]; ++v)
                if (a[v] == TERM_TRUE) {
                    prefix[k] = v;
                    out.push_back(prefix);
                }
            return;
        }
        for (int v = 0; v < domain_[k]; ++v) {
            if (a[v] != TERM_FALSE) {
                prefix[k] = v;
                enum_below(k + 1, a[v], prefix, out);
            }
        }
    }

    std::vector<int> domain_;
    std::size_t K_;
    /** node_[k][id-1] holds the arc row of level-k node id. */
    std::vector<std::vector<std::vector<int>>> node_;
    /** Per-level unique table, arc row -> node id. */
    std::vector<std::map<std::vector<int>, int>> uniq_;
    int root_;
    /** Memoized per-node state counts; -1 marks "not yet computed". */
    mutable std::vector<std::vector<long long>> cnt_;
    mutable bool dirty_;
};

/** Human-readable storage summary, the twin of MDD.toString. */
inline std::string mdd_to_string(const MDD& m) {
    const MddStats s = m.stats();
    std::string out = "  MDD  " + std::to_string(m.K()) + " levels, " +
                      std::to_string(s.num_states) + " states in " +
                      std::to_string(s.num_nodes) + " nodes, footprint " +
                      std::to_string(s.mdd_ints) + " ints vs " +
                      std::to_string(s.explicit_ints) + " explicit";
    if (s.table_nodes > s.num_nodes)
        out += " (" + std::to_string(s.table_nodes - s.num_nodes) +
               " dead nodes in tables; call compact() to reclaim)";
    return out;
}

}  // namespace mdd
}  // namespace line

#endif  // LINE_API_MDD_MDD_H
