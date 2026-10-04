classdef MDD < handle
    % MDD - Quasi-reduced ordered Multi-valued Decision Diagram.
    %
    % Compact symbolic store for a set of discrete states, after
    %   A.S. Miner, G. Ciardo, "Efficient Reachability Set Generation and
    %   Storage Using Decision Diagrams", ICATPN 1999, LNCS 1639, pp.6-25.
    %
    % A global state is a K-tuple of *local* state values (one per level/
    % submodel), state(k) in {0,...,domain(k)-1}. The set is stored as a
    % directed acyclic graph with K variable levels plus a terminal level:
    % level 1 is the top (root), a node at level k has domain(k) outgoing
    % arcs to level k+1 nodes, and a state belongs to the set iff its path
    % of arcs reaches the TRUE terminal. Canonicity is enforced by a
    % per-level unique table (no duplicate nodes) and by collapsing the
    % all-FALSE node to the FALSE terminal. Storage is O(#nodes), typically
    % O(K * #local-states), instead of O(|S|) as in an explicit state list.
    %
    % Two constant terminals encode the boolean value of a completed path:
    %   TERM_FALSE = 0  (empty subgraph / state not in set)
    %   TERM_TRUE  = -1 (state in set)
    % Arcs of a level-k node hold ids of level-(k+1) nodes when k<K, or a
    % terminal (TERM_TRUE/TERM_FALSE) when k==K.
    %
    % Example:
    %   m = MDD([3 3 3]);        % 3 levels, local values 0..2
    %   m.insert([0 1 2]);
    %   m.insert([2 0 0]);
    %   m.member([0 1 2])        % -> true
    %   m.cardinality()          % -> 2
    %   m.index([2 0 0])         % -> 0-based rank in the set
    %   S = m.enumerate();       % all states, one per row, in index order
    %
    % See also: mdd_reachset, ctmc_ssg.

    properties (Constant)
        TERM_TRUE  = -1;   % terminal node "1": completed path is accepted
        TERM_FALSE =  0;   % terminal node "0": empty subgraph
    end

    properties
        K         % number of variable levels
        domain    % 1 x K, domain(k) = number of local states at level k
        node      % 1 x K cell, node{k} is (nnodes_k x domain(k)) of child ids
        uniq      % 1 x K cell of dictionary(string->double), unique table per level
        root      % id of the top (level-1) node; TERM_FALSE for the empty set
    end

    properties (Access = private)
        nnodes    % 1 x K live row count per level (node{k} is preallocated with slack)
        cnt       % 1 x K cell of per-node state counts (-1 = not yet computed)
        dirty     % true when cnt must be rebuilt after an insert
    end

    methods
        function obj = MDD(domain)
            % MDD(domain) creates an empty set over the given per-level domains.
            % domain is a 1 x K vector; state values at level k are 0..domain(k)-1.
            obj.domain = double(domain(:)');
            obj.K = numel(obj.domain);
            obj.node = cell(1, obj.K);
            obj.uniq = cell(1, obj.K);
            for k = 1:obj.K
                obj.node{k} = zeros(0, obj.domain(k));
                obj.uniq{k} = configureDictionary('string', 'double');
            end
            obj.nnodes = zeros(1, obj.K);
            obj.root = obj.TERM_FALSE;
            obj.cnt = cell(1, obj.K);
            obj.dirty = true;
        end

        function insert(obj, state)
            % INSERT(state) adds a K-tuple (row vector, 0-based values) to the set.
            obj.root = obj.addState(1, obj.root, state);
            obj.dirty = true;
        end

        function tf = member(obj, state)
            % MEMBER(state) returns true iff state is in the set (O(K)).
            id = obj.root;
            for k = 1:obj.K
                if id == obj.TERM_FALSE
                    tf = false; return
                end
                id = obj.node{k}(id, state(k) + 1);
            end
            tf = (id == obj.TERM_TRUE);
        end

        function n = cardinality(obj)
            % CARDINALITY() returns |S|, the number of stored states.
            obj.ensureCounts();
            n = obj.childCount(1, obj.root);
        end

        function idx = index(obj, state)
            % INDEX(state) returns the 0-based lexicographic rank of state
            % among the stored set (level 1 most significant), or -1 if state
            % is not in the set. This is a bijection S <-> {0,...,|S|-1}, so a
            % generator matrix can be assembled without an explicit state list.
            obj.ensureCounts();
            idx = 0;
            id = obj.root;
            for k = 1:obj.K
                if id == obj.TERM_FALSE
                    idx = -1; return
                end
                arcs = obj.node{k}(id, :);
                v = state(k);
                for vv = 0:(v - 1)
                    idx = idx + obj.childCount(k + 1, arcs(vv + 1));
                end
                id = arcs(v + 1);
            end
            if id ~= obj.TERM_TRUE
                idx = -1;
            end
        end

        function S = enumerate(obj)
            % ENUMERATE() returns all stored states as rows, in index() order.
            if obj.root == obj.TERM_FALSE
                S = zeros(0, obj.K); return
            end
            S = obj.enumBelow(1, obj.root);
        end

        function s = toStruct(obj)
            % TOSTRUCT() exports the diagram as plain arrays for downstream
            % algorithms (e.g. mdd_mcd), trimmed to live rows. Fields:
            %   K, domain, root  - as the properties
            %   nnodes           - 1 x K live node count per level
            %   node             - 1 x K cell, node{k} is nnodes(k) x domain(k)
            %                      of child ids (level k+1 ids, or terminals at k=K)
            s.K = obj.K;
            s.domain = obj.domain;
            s.root = obj.root;
            s.nnodes = obj.nnodes;
            s.node = cell(1, obj.K);
            for k = 1:obj.K
                s.node{k} = obj.node{k}(1:obj.nnodes(k), :);
            end
        end

        function s = stats(obj)
            % STATS() returns a struct describing the storage of the current
            % set. Only nodes reachable from the root are counted (the build is
            % append-only, so superseded nodes may linger in the tables until
            % COMPACT is called; see COMPACT):
            %   numStates    - |S|
            %   numNodes     - reachable non-terminal nodes
            %   nodesPerLevel- 1 x K reachable node counts
            %   liveNodes    - reachable nodes (= numNodes)
            %   tableNodes   - nodes physically held in the tables (incl. dead)
            %   mddInts      - integers in the reachable arc arrays (footprint)
            %   explicitInts - integers an explicit state list needs (|S|*K)
            %   compression  - explicitInts / mddInts
            vis = obj.reachableIds();
            s.levels = obj.K;
            s.nodesPerLevel = cellfun(@(b) sum(b), vis);
            s.numNodes = sum(s.nodesPerLevel);
            s.liveNodes = s.numNodes;
            s.tableNodes = sum(obj.nnodes);
            s.numStates = obj.cardinality();
            s.mddInts = sum(arrayfun(@(k) s.nodesPerLevel(k) * obj.domain(k), 1:obj.K));
            s.explicitInts = s.numStates * obj.K;
            s.compression = s.explicitInts / max(s.mddInts, 1);
        end

        function compact(obj)
            % COMPACT() reclaims dead nodes left by the append-only build,
            % rebuilding the level tables and unique tables so that only nodes
            % reachable from the root remain. Membership/index/enumerate are
            % unchanged. (A production MDD would reference-count instead and
            % never accumulate dead nodes; this is the basic sweep.)
            vis = obj.reachableIds();
            newnode = cell(1, obj.K);
            remap = cell(1, obj.K);
            newcount = zeros(1, obj.K);
            for k = 1:obj.K
                ids = find(vis{k});
                remap{k} = zeros(obj.nnodes(k), 1);
                remap{k}(ids) = 1:numel(ids);
                newnode{k} = obj.node{k}(ids, :);
                newcount(k) = numel(ids);
            end
            for k = 1:(obj.K - 1)
                A = newnode{k};
                mask = A > 0;                 % positive entries are child ids
                A(mask) = remap{k + 1}(A(mask));
                newnode{k} = A;
            end
            obj.node = newnode;
            obj.nnodes = newcount;
            if obj.root ~= obj.TERM_FALSE
                obj.root = remap{1}(obj.root);
            end
            for k = 1:obj.K
                obj.uniq{k} = configureDictionary('string', 'double');
                for p = 1:obj.nnodes(k)
                    % key encoding MUST match MAKENODE
                    key = string(char(uint16(typecast(int32(obj.node{k}(p, :)), 'uint8')) + 256));
                    obj.uniq{k}(key) = p;
                end
            end
            obj.dirty = true;
        end

        function disp(obj)
            s = obj.stats();
            fprintf('  MDD  %d levels, domains [%s]\n', obj.K, strtrim(sprintf('%d ', obj.domain)));
            fprintf('       %d states stored in %d nodes (%s per level)\n', ...
                s.numStates, s.numNodes, strtrim(sprintf('%d ', s.nodesPerLevel)));
            fprintf('       footprint %d ints vs %d explicit (%.1fx compression)\n', ...
                s.mddInts, s.explicitInts, s.compression);
            if s.tableNodes > s.liveNodes
                fprintf('       (%d dead nodes in tables; call compact() to reclaim)\n', ...
                    s.tableNodes - s.liveNodes);
            end
        end
    end

    methods (Static)
        function obj = fromStates(domain, S)
            % MDD.fromStates(domain, S) builds an MDD from a matrix S whose
            % rows are 0-based state tuples over the given per-level domains.
            obj = MDD(domain);
            for i = 1:size(S, 1)
                obj.insert(S(i, :));
            end
        end
    end

    methods (Access = private)
        function id = makeNode(obj, k, arcs)
            % Canonical node creation with the per-level unique table.
            if all(arcs == obj.TERM_FALSE)
                id = obj.TERM_FALSE; return   % collapse empty node
            end
            % Exact string key of the arc row: each int32 is split into 4 bytes
            % offset into the 256..511 code-point range (null/surrogate-safe and
            % valid for any node count). This exact encoding MUST match the one
            % used when the unique table is rebuilt in COMPACT.
            key = string(char(uint16(typecast(int32(arcs), 'uint8')) + 256));
            if isKey(obj.uniq{k}, key)
                id = obj.uniq{k}(key); return
            end
            % append into a table preallocated by doubling to avoid O(n^2) grow
            id = obj.nnodes(k) + 1;
            if id > size(obj.node{k}, 1)
                obj.node{k}(max(16, 2 * size(obj.node{k}, 1)), 1) = 0;
            end
            obj.node{k}(id, :) = arcs;
            obj.nnodes(k) = id;
            obj.uniq{k}(key) = id;        % dictionary is a value type: write in place
        end

        function id = addState(obj, k, id, state)
            % Recursively add one state below node id at level k, returning the
            % (possibly new) canonical id. Nodes are immutable/shared, so this
            % rebuilds the path bottom-up rather than mutating in place.
            if k > obj.K
                id = obj.TERM_TRUE; return
            end
            if id == obj.TERM_FALSE
                arcs = zeros(1, obj.domain(k));   % all TERM_FALSE
            else
                arcs = obj.node{k}(id, :);
            end
            v = state(k);
            arcs(v + 1) = obj.addState(k + 1, arcs(v + 1), state);
            id = obj.makeNode(k, arcs);
        end

        function c = childCount(obj, k, childId)
            % Number of accepted states below a child reference at level k
            % (childId is a terminal when k > K).
            if k > obj.K
                c = double(childId == obj.TERM_TRUE); return
            end
            if childId == obj.TERM_FALSE
                c = 0; return
            end
            c = obj.countNode(k, childId);
        end

        function c = countNode(obj, k, id)
            c = obj.cnt{k}(id);
            if c >= 0, return; end
            arcs = obj.node{k}(id, :);
            c = 0;
            for v = 1:obj.domain(k)
                c = c + obj.childCount(k + 1, arcs(v));
            end
            obj.cnt{k}(id) = c;
        end

        function ensureCounts(obj)
            if ~obj.dirty && ~isempty(obj.cnt) && all(cellfun(@numel, obj.cnt) == obj.nnodes)
                return
            end
            for k = 1:obj.K
                obj.cnt{k} = -ones(obj.nnodes(k), 1);
            end
            obj.dirty = false;
        end

        function vis = reachableIds(obj)
            % Per-level logical masks of nodes reachable from the root.
            vis = cell(1, obj.K);
            for k = 1:obj.K
                vis{k} = false(obj.nnodes(k), 1);
            end
            if obj.root == obj.TERM_FALSE
                return
            end
            vis{1}(obj.root) = true;
            stack = {[1, obj.root]};
            while ~isempty(stack)
                top = stack{end}; stack(end) = [];
                k = top(1); id = top(2);
                if k == obj.K
                    continue                 % children are terminals
                end
                arcs = obj.node{k}(id, :);
                for v = 1:obj.domain(k)
                    ch = arcs(v);
                    if ch > 0 && ~vis{k + 1}(ch)
                        vis{k + 1}(ch) = true;
                        stack{end + 1} = [k + 1, ch]; %#ok<AGROW>
                    end
                end
            end
        end

        function S = enumBelow(obj, k, id)
            % All states over columns k..K reachable to TRUE below node id.
            arcs = obj.node{k}(id, :);
            if k == obj.K
                vals = find(arcs == obj.TERM_TRUE) - 1;   % 0-based accepted values
                S = vals(:);
                return
            end
            S = zeros(0, obj.K - k + 1);
            for v = 0:(obj.domain(k) - 1)
                child = arcs(v + 1);
                if child ~= obj.TERM_FALSE
                    sub = obj.enumBelow(k + 1, child);
                    m = size(sub, 1);
                    S = [S; [repmat(v, m, 1), sub]]; %#ok<AGROW>
                end
            end
        end
    end
end
