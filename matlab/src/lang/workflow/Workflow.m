classdef Workflow < Model
    % Workflow - A computational workflow that can be converted to a PH distribution.
    %
    % Workflow allows declaring computational workflows using the same
    % activity graph syntax as LayeredNetwork. The workflow can then be
    % converted to an equivalent phase-type (PH) distribution via explicit
    % CTMC construction.
    %
    % Supported precedence patterns:
    %   - Serial: Sequential execution
    %   - AndFork/AndJoin: Parallel execution with synchronization
    %   - OrFork/OrJoin: Probabilistic branching
    %   - Loop: Repeated execution
    %
    % Example:
    %   wf = Workflow('myWorkflow');
    %   A = wf.addActivity('A', Exp.fitMean(1.0));
    %   B = wf.addActivity('B', Exp.fitMean(2.0));
    %   wf.addPrecedence(Workflow.Serial(A, B));
    %   ph = wf.toPH();
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    properties
        activities = {};        % Cell array of WorkflowActivity objects
        precedences = [];       % Array of ActivityPrecedence objects
        activityMap;            % Map from activity name to index
    end

    properties (Hidden)
        cachedPH;               % Cached phase-type distribution result
    end

    methods
        function self = Workflow(name)
            % WORKFLOW Create a new Workflow instance
            %
            % SELF = WORKFLOW(NAME)
            %
            % Parameters:
            %   name - Name for the workflow

            self@Model(name);
            self.activities = {};
            self.precedences = [];
            self.activityMap = containers.Map();
            self.cachedPH = [];
        end

        function act = addActivity(self, name, hostDemand)
            % ADDACTIVITY Add an activity to the workflow
            %
            % ACT = ADDACTIVITY(SELF, NAME, HOSTDEMAND)
            %
            % Parameters:
            %   name       - Activity name (string)
            %   hostDemand - Service time distribution or numeric mean
            %
            % Returns:
            %   act - WorkflowActivity object

            if nargin < 3
                hostDemand = GlobalConstants.FineTol;
            end

            act = WorkflowActivity(self, name, hostDemand);
            self.activities{end+1} = act;
            act.index = length(self.activities);
            self.activityMap(name) = act.index;
            self.cachedPH = [];  % Invalidate cache
        end

        function self = addPrecedence(self, prec)
            % ADDPRECEDENCE Add precedence constraints to the workflow
            %
            % SELF = ADDPRECEDENCE(SELF, PREC)
            %
            % Parameters:
            %   prec - ActivityPrecedence object or cell array (from Serial)

            if iscell(prec)
                for m = 1:length(prec)
                    self.precedences = [self.precedences; prec{m}];
                end
            else
                self.precedences = [self.precedences; prec];
            end
            self.cachedPH = [];  % Invalidate cache
        end

        function act = getActivity(self, name)
            % GETACTIVITY Get activity by name
            %
            % ACT = GETACTIVITY(SELF, NAME)
            %
            % Parameters:
            %   name - Activity name (string)
            %
            % Returns:
            %   act - WorkflowActivity object

            if isa(name, 'WorkflowActivity')
                act = name;
                return;
            end
            if ~isKey(self.activityMap, name)
                line_error(mfilename, sprintf('Activity "%s" not found in workflow.', name));
            end
            idx = self.activityMap(name);
            act = self.activities{idx};
        end

        function idx = getActivityIndex(self, name)
            % GETACTIVITYINDEX Get activity index by name
            %
            % IDX = GETACTIVITYINDEX(SELF, NAME)

            if isa(name, 'WorkflowActivity')
                idx = name.index;
            elseif isKey(self.activityMap, name)
                idx = self.activityMap(name);
            else
                idx = -1;
            end
        end

        function [isValid, msg] = validate(self)
            % VALIDATE Validate workflow structure
            %
            % [ISVALID, MSG] = VALIDATE(SELF)
            %
            % Returns:
            %   isValid - true if workflow is valid
            %   msg     - Error message if invalid

            isValid = true;
            msg = '';

            % Check 1: At least one activity
            if isempty(self.activities)
                isValid = false;
                msg = 'Workflow must have at least one activity.';
                return;
            end

            % Check 2: All activities in precedences exist
            for p = 1:length(self.precedences)
                prec = self.precedences(p);
                for a = 1:length(prec.preActs)
                    if ~isKey(self.activityMap, prec.preActs{a})
                        isValid = false;
                        msg = sprintf('Activity "%s" referenced in precedence not found in workflow.', prec.preActs{a});
                        return;
                    end
                end
                for a = 1:length(prec.postActs)
                    if ~isKey(self.activityMap, prec.postActs{a})
                        isValid = false;
                        msg = sprintf('Activity "%s" referenced in precedence not found in workflow.', prec.postActs{a});
                        return;
                    end
                end
            end

            % Check 3: OR-fork probabilities sum to 1
            for p = 1:length(self.precedences)
                prec = self.precedences(p);
                if prec.postType == ActivityPrecedenceType.POST_OR
                    if isempty(prec.postParams)
                        isValid = false;
                        msg = 'OR-fork must have probabilities specified.';
                        return;
                    end
                    if abs(sum(prec.postParams) - 1.0) > GlobalConstants.FineTol
                        isValid = false;
                        msg = sprintf('OR-fork probabilities must sum to 1 (got %.4f).', sum(prec.postParams));
                        return;
                    end
                end
            end

            % Check 4: Loop counts must be positive
            for p = 1:length(self.precedences)
                prec = self.precedences(p);
                if prec.postType == ActivityPrecedenceType.POST_LOOP
                    if isempty(prec.postParams) || prec.postParams <= 0
                        isValid = false;
                        msg = 'Loop count must be a positive number.';
                        return;
                    end
                end
            end
        end

        function ph = toPH(self)
            % TOPH Convert workflow to phase-type distribution
            %
            % PH = TOPH(SELF)
            %
            % Returns:
            %   ph - APH distribution representing the workflow execution time

            if ~isempty(self.cachedPH)
                ph = self.cachedPH;
                return;
            end

            [isValid, msg] = self.validate();
            if ~isValid
                line_error(mfilename, msg);
            end

            [alpha, T] = self.buildCTMC();
            ph = APH(alpha, T);
            self.cachedPH = ph;
        end

        function [alpha, T] = buildCTMC(self)
            % BUILDCTMC Build the CTMC representation of the workflow
            %
            % [ALPHA, T] = BUILDCTMC(SELF)
            %
            % Returns:
            %   alpha - Initial probability vector
            %   T     - Subgenerator matrix

            % Handle trivial case: single activity
            if length(self.activities) == 1
                [alpha, T] = self.activities{1}.getPHRepresentation();
                return;
            end

            % Build adjacency structures from precedences
            [adjList, inDegree, outDegree, forkInfo, joinInfo, loopInfo] = self.analyzeStructure();

            % Find start activities (no predecessors in precedence graph)
            startActs = find(inDegree == 0);
            if isempty(startActs)
                % If all activities have predecessors, use first activity
                startActs = 1;
            end

            % Find end activities (no successors in precedence graph)
            endActs = find(outDegree == 0);
            if isempty(endActs)
                endActs = length(self.activities);
            end

            % Build the workflow CTMC using recursive block composition
            [alpha, T] = self.composeWorkflow(startActs, endActs, adjList, forkInfo, joinInfo, loopInfo);
        end
    end

    methods (Access = private)
        function [adjList, inDegree, outDegree, forkInfo, joinInfo, loopInfo] = analyzeStructure(self)
            % ANALYZESTRUCTURE Analyze workflow structure from precedences
            %
            % Build adjacency list and identify fork/join/loop structures

            n = length(self.activities);
            adjList = cell(n, 1);
            inDegree = zeros(n, 1);
            outDegree = zeros(n, 1);
            forkInfo = struct('type', {}, 'preAct', {}, 'postActs', {}, 'probs', {});
            joinInfo = struct('type', {}, 'preActs', {}, 'postAct', {}, 'quorum', {});
            loopInfo = struct('preAct', {}, 'loopActs', {}, 'endAct', {}, 'count', {});

            for p = 1:length(self.precedences)
                prec = self.precedences(p);

                % Get activity indices
                preInds = zeros(length(prec.preActs), 1);
                for a = 1:length(prec.preActs)
                    preInds(a) = self.getActivityIndex(prec.preActs{a});
                end

                postInds = zeros(length(prec.postActs), 1);
                for a = 1:length(prec.postActs)
                    postInds(a) = self.getActivityIndex(prec.postActs{a});
                end

                % Update adjacency and degrees
                for i = 1:length(preInds)
                    for j = 1:length(postInds)
                        adjList{preInds(i)} = [adjList{preInds(i)}, postInds(j)];
                        outDegree(preInds(i)) = outDegree(preInds(i)) + 1;
                        inDegree(postInds(j)) = inDegree(postInds(j)) + 1;
                    end
                end

                % Record fork/join/loop info
                switch prec.postType
                    case ActivityPrecedenceType.POST_AND
                        % AND-Fork
                        forkInfo(end+1) = struct('type', 'and', ...
                            'preAct', preInds(1), ...
                            'postActs', postInds(:)', ...
                            'probs', []);

                    case ActivityPrecedenceType.POST_OR
                        % OR-Fork
                        forkInfo(end+1) = struct('type', 'or', ...
                            'preAct', preInds(1), ...
                            'postActs', postInds(:)', ...
                            'probs', prec.postParams(:)');

                    case ActivityPrecedenceType.POST_LOOP
                        % Loop
                        count = prec.postParams;
                        if length(postInds) >= 2
                            loopInfo(end+1) = struct('preAct', preInds(1), ...
                                'loopActs', postInds(1:end-1)', ...
                                'endAct', postInds(end), ...
                                'count', count);
                        else
                            loopInfo(end+1) = struct('preAct', preInds(1), ...
                                'loopActs', postInds(1)', ...
                                'endAct', [], ...
                                'count', count);
                        end
                end

                switch prec.preType
                    case ActivityPrecedenceType.PRE_AND
                        % AND-Join
                        quorum = prec.preParams;
                        if isempty(quorum)
                            quorum = length(preInds);
                        end
                        joinInfo(end+1) = struct('type', 'and', ...
                            'preActs', preInds(:)', ...
                            'postAct', postInds(1), ...
                            'quorum', quorum);

                    case ActivityPrecedenceType.PRE_OR
                        % OR-Join
                        joinInfo(end+1) = struct('type', 'or', ...
                            'preActs', preInds(:)', ...
                            'postAct', postInds(1), ...
                            'quorum', []);
                end
            end
        end

        function [alpha, T] = composeWorkflow(self, startActs, endActs, adjList, forkInfo, joinInfo, loopInfo)
            % COMPOSEWORKFLOW Compose the workflow CTMC
            %
            % Uses a simplified approach: process precedences in order,
            % building up the CTMC incrementally.

            n = length(self.activities);

            % Handle special cases
            if n == 1
                [alpha, T] = self.activities{1}.getPHRepresentation();
                return;
            end

            % Determine workflow structure and compose accordingly
            % Check for simple serial workflow (no forks/joins/loops)
            if isempty(forkInfo) && isempty(joinInfo) && isempty(loopInfo)
                [alpha, T] = self.composeSerialWorkflow(adjList, startActs);
                return;
            end

            % For complex workflows, build block-by-block
            [alpha, T] = self.composeComplexWorkflow(adjList, forkInfo, joinInfo, loopInfo);
        end

        function [alpha, T] = composeSerialWorkflow(self, adjList, startActs)
            % COMPOSESERIALWORKFLOW Compose a purely serial workflow
            %
            % Activities are concatenated in topological order.

            n = length(self.activities);

            % Topological sort
            order = self.topologicalSort(adjList);

            % Compose in order
            [alpha, T] = self.activities{order(1)}.getPHRepresentation();

            for i = 2:length(order)
                [alpha_i, T_i] = self.activities{order(i)}.getPHRepresentation();
                [alpha, T] = Workflow.composeSerial(alpha, T, alpha_i, T_i);
            end
        end

        function order = topologicalSort(self, adjList)
            % TOPOLOGICALSORT Perform topological sort on activity graph

            n = length(self.activities);
            inDeg = zeros(n, 1);

            % Calculate in-degrees
            for i = 1:n
                for j = adjList{i}
                    inDeg(j) = inDeg(j) + 1;
                end
            end

            % BFS-based topological sort
            queue = find(inDeg == 0);
            if isempty(queue)
                queue = 1;  % Start with first activity if no clear start
            end
            order = [];

            while ~isempty(queue)
                curr = queue(1);
                queue(1) = [];
                order = [order, curr];

                for next = adjList{curr}
                    inDeg(next) = inDeg(next) - 1;
                    if inDeg(next) == 0
                        queue = [queue, next];
                    end
                end
            end

            % Add any remaining activities not in precedence graph
            remaining = setdiff(1:n, order);
            order = [order, remaining];
        end

        function [alpha, T] = composeComplexWorkflow(self, adjList, forkInfo, joinInfo, loopInfo)
            % COMPOSECOMPLEXWORKFLOW Compose workflow with forks/joins/loops
            %
            % Strategy: Build a stage-based representation where each stage
            % is either a single activity, a loop block, or a fork-join block.
            % Stages are then composed in topological order.
            %
            % The key insight is that fork-join creates parallel branches that
            % need Kronecker sum expansion for proper state space modeling.
            % The map_max structure is used for join synchronization:
            %   [krons(A,B), trans_to_onlyA, trans_to_onlyB]
            %   [0,          A,              0             ]
            %   [0,          0,              B             ]

            n = length(self.activities);

            % Build a representation where each activity or block has a PH
            blockAlpha = cell(n, 1);
            blockT = cell(n, 1);

            % blockRepresentative(i) = j means activity i's result is stored in block j
            % This handles the case where multiple activities are merged into one block
            blockRepresentative = (1:n)';

            % isBlockHead(i) = true means block i should be included in final composition
            isBlockHead = true(n, 1);

            % Initialize each activity as its own block
            for i = 1:n
                [blockAlpha{i}, blockT{i}] = self.activities{i}.getPHRepresentation();
            end

            % Process loops first (they modify the activity's PH)
            for k = 1:length(loopInfo)
                loop = loopInfo(k);
                preIdx = loop.preAct;
                loopActInds = loop.loopActs;
                endIdx = loop.endAct;
                count = loop.count;

                % Compose loop body
                if length(loopActInds) == 1
                    [alpha_loop, T_loop] = self.activities{loopActInds(1)}.getPHRepresentation();
                else
                    % Serial composition of loop activities
                    [alpha_loop, T_loop] = self.activities{loopActInds(1)}.getPHRepresentation();
                    for j = 2:length(loopActInds)
                        [a_j, T_j] = self.activities{loopActInds(j)}.getPHRepresentation();
                        [alpha_loop, T_loop] = Workflow.composeSerial(alpha_loop, T_loop, a_j, T_j);
                    end
                end

                % Create count-fold convolution
                [alpha_conv, T_conv] = Workflow.composeRepeat(alpha_loop, T_loop, count);

                % Compose with pre-activity
                [alpha_result, T_result] = Workflow.composeSerial(...
                    blockAlpha{preIdx}, blockT{preIdx}, alpha_conv, T_conv);

                % Compose with end activity if present
                if ~isempty(endIdx)
                    [alpha_end, T_end] = self.activities{endIdx}.getPHRepresentation();
                    [alpha_result, T_result] = Workflow.composeSerial(alpha_result, T_result, alpha_end, T_end);
                    % Merge end activity into preIdx block
                    blockRepresentative(endIdx) = preIdx;
                    isBlockHead(endIdx) = false;
                end

                blockAlpha{preIdx} = alpha_result;
                blockT{preIdx} = T_result;

                % Merge loop activities into preIdx block
                for j = 1:length(loopActInds)
                    blockRepresentative(loopActInds(j)) = preIdx;
                    isBlockHead(loopActInds(j)) = false;
                end
            end

            % Process AND-forks and their matching AND-joins
            % These create parallel blocks with Kronecker sum state expansion
            for k = 1:length(forkInfo)
                fork = forkInfo(k);
                if strcmp(fork.type, 'and')
                    % Find matching AND-join
                    joinIdx = self.findMatchingJoin(fork.postActs, joinInfo, 'and');

                    if ~isempty(joinIdx)
                        join = joinInfo(joinIdx);
                        preIdx = fork.preAct;
                        postIdx = join.postAct;
                        parallelInds = fork.postActs;

                        % Compose parallel activities using Kronecker sum (map_max structure)
                        [alpha_par, T_par] = self.composeAndForkBlock(parallelInds, blockAlpha, blockT);

                        % Store parallel block result - DO NOT compose with pre/post here
                        % The pre/post connections will be handled by topological composition

                        % Create a new "virtual" block for the parallel section
                        % We store it in the first parallel activity's slot
                        firstParallel = parallelInds(1);
                        blockAlpha{firstParallel} = alpha_par;
                        blockT{firstParallel} = T_par;

                        % Mark other parallel activities as merged into firstParallel
                        for j = 2:length(parallelInds)
                            blockRepresentative(parallelInds(j)) = firstParallel;
                            isBlockHead(parallelInds(j)) = false;
                        end

                        % Update adjacency: preIdx -> firstParallel -> postIdx
                        % The fork and join activities keep their blocks separate
                    end
                end
            end

            % Process OR-forks and their matching OR-joins
            for k = 1:length(forkInfo)
                fork = forkInfo(k);
                if strcmp(fork.type, 'or')
                    % Find matching OR-join
                    joinIdx = self.findMatchingJoin(fork.postActs, joinInfo, 'or');

                    branchInds = fork.postActs;
                    probs = fork.probs;

                    % Compose branches as probabilistic mixture
                    [alpha_or, T_or] = self.composeOrForkBlock(branchInds, probs, blockAlpha, blockT);

                    % Store in first branch's slot
                    firstBranch = branchInds(1);
                    blockAlpha{firstBranch} = alpha_or;
                    blockT{firstBranch} = T_or;

                    % Mark other branches as merged
                    for j = 2:length(branchInds)
                        blockRepresentative(branchInds(j)) = firstBranch;
                        isBlockHead(branchInds(j)) = false;
                    end
                end
            end

            % Build modified adjacency list that respects block merging
            modAdjList = cell(n, 1);
            for i = 1:n
                if isBlockHead(i)
                    % Find successors, mapping through representatives
                    successors = [];
                    for succ = adjList{i}
                        rep = blockRepresentative(succ);
                        if rep ~= i && isBlockHead(rep) && ~ismember(rep, successors)
                            successors = [successors, rep];
                        end
                    end
                    modAdjList{i} = successors;
                end
            end

            % Topological sort on block heads only
            order = self.topologicalSortBlocks(modAdjList, isBlockHead);

            % Compose blocks in topological order
            alpha = [];
            T = [];

            for i = 1:length(order)
                idx = order(i);
                if isempty(alpha)
                    alpha = blockAlpha{idx};
                    T = blockT{idx};
                else
                    [alpha, T] = Workflow.composeSerial(alpha, T, blockAlpha{idx}, blockT{idx});
                end
            end

            % If still empty (shouldn't happen), use first activity
            if isempty(alpha)
                [alpha, T] = self.activities{1}.getPHRepresentation();
            end
        end

        function order = topologicalSortBlocks(self, adjList, isBlockHead)
            % TOPOLOGICALSORTBLOCKS Topological sort considering only block heads

            n = length(self.activities);
            inDeg = zeros(n, 1);

            % Calculate in-degrees for block heads only
            for i = 1:n
                if isBlockHead(i)
                    for j = adjList{i}
                        if isBlockHead(j)
                            inDeg(j) = inDeg(j) + 1;
                        end
                    end
                end
            end

            % Find starting nodes (block heads with in-degree 0)
            queue = [];
            for i = 1:n
                if isBlockHead(i) && inDeg(i) == 0
                    queue = [queue, i];
                end
            end

            if isempty(queue)
                % Fallback: use first block head
                for i = 1:n
                    if isBlockHead(i)
                        queue = i;
                        break;
                    end
                end
            end

            order = [];
            while ~isempty(queue)
                curr = queue(1);
                queue(1) = [];
                order = [order, curr];

                for next = adjList{curr}
                    if isBlockHead(next)
                        inDeg(next) = inDeg(next) - 1;
                        if inDeg(next) == 0
                            queue = [queue, next];
                        end
                    end
                end
            end

            % Add any remaining block heads not reached
            for i = 1:n
                if isBlockHead(i) && ~ismember(i, order)
                    order = [order, i];
                end
            end
        end

        function joinIdx = findMatchingJoin(self, postActs, joinInfo, joinType)
            % FINDMATCHINGJOIN Find the join that matches a fork's post activities

            joinIdx = [];
            for j = 1:length(joinInfo)
                join = joinInfo(j);
                if strcmp(join.type, joinType)
                    % Check if this join's pre-activities match the fork's post-activities
                    if isequal(sort(join.preActs), sort(postActs))
                        joinIdx = j;
                        return;
                    end
                end
            end
        end

        function [alpha, T] = composeAndForkBlock(self, parallelInds, blockAlpha, blockT)
            % COMPOSEANDFORKBLOCK Compose parallel activities using Kronecker sum

            % Start with first parallel activity
            alpha = blockAlpha{parallelInds(1)};
            T = blockT{parallelInds(1)};

            % Combine with remaining activities using Kronecker sum
            for i = 2:length(parallelInds)
                idx = parallelInds(i);
                alpha_i = blockAlpha{idx};
                T_i = blockT{idx};

                [alpha, T] = Workflow.composeParallel(alpha, T, alpha_i, T_i);
            end
        end

        function [alpha, T] = composeOrForkBlock(self, branchInds, probs, blockAlpha, blockT)
            % COMPOSEORFORKBLOCK Compose branches as probabilistic mixture

            % Calculate total phases
            totalPhases = 0;
            for i = 1:length(branchInds)
                totalPhases = totalPhases + size(blockT{branchInds(i)}, 1);
            end

            T = zeros(totalPhases);
            alpha = zeros(1, totalPhases);

            offset = 0;
            for i = 1:length(branchInds)
                idx = branchInds(i);
                alpha_i = blockAlpha{idx};
                T_i = blockT{idx};
                n_i = size(T_i, 1);

                T(offset+1:offset+n_i, offset+1:offset+n_i) = T_i;
                alpha(offset+1:offset+n_i) = probs(i) * alpha_i;

                offset = offset + n_i;
            end
        end
    end

    methods (Static)
        function [alpha_out, T_out] = composeSerial(alpha1, T1, alpha2, T2)
            % COMPOSESERIAL Serial composition of two PH distributions
            %
            % The second PH starts when the first one absorbs.
            %
            % T_serial = [T1,  -T1*e*alpha2]
            %            [0,   T2         ]

            n1 = size(T1, 1);
            n2 = size(T2, 1);

            e1 = ones(n1, 1);
            absRate1 = -T1 * e1;  % Absorption rate from each state

            T_out = zeros(n1 + n2);
            T_out(1:n1, 1:n1) = T1;
            T_out(1:n1, n1+1:end) = absRate1 * alpha2;
            T_out(n1+1:end, n1+1:end) = T2;

            alpha_out = [alpha1, zeros(1, n2)];
        end

        function [alpha_out, T_out] = composeParallel(alpha1, T1, alpha2, T2)
            % COMPOSEPARALLEL Parallel fork-join composition
            %
            % Models two PH distributions running in parallel with synchronization
            % at completion (AND-join). The resulting PH represents the time until
            % BOTH activities complete (i.e., the maximum).
            %
            % State space:
            %   - States (i,j) where both are active: i=1..n1, j=1..n2
            %   - States where only activity 1 is active (2 completed): i=1..n1
            %   - States where only activity 2 is active (1 completed): j=1..n2
            %
            % Absorption occurs only when both have completed.

            n1 = size(T1, 1);
            n2 = size(T2, 1);
            e1 = ones(n1, 1);
            e2 = ones(n2, 1);

            % Absorption rates for each activity
            absRate1 = -T1 * e1;  % Rate of leaving each phase in activity 1
            absRate2 = -T2 * e2;  % Rate of leaving each phase in activity 2

            % Total state space: n1*n2 (both active) + n1 (only 1 active) + n2 (only 2 active)
            % States 1..n1*n2: both active, indexed as (i-1)*n2 + j for phase (i,j)
            % States n1*n2+1..n1*n2+n1: only activity 1 active (activity 2 completed)
            % States n1*n2+n1+1..n1*n2+n1+n2: only activity 2 active (activity 1 completed)

            nBoth = n1 * n2;
            nOnly1 = n1;
            nOnly2 = n2;
            nTotal = nBoth + nOnly1 + nOnly2;

            T_out = zeros(nTotal);

            % Transitions within "both active" states
            % Use Kronecker sum for simultaneous evolution
            T_both = krons(T1, T2);
            T_out(1:nBoth, 1:nBoth) = T_both;

            % Transitions from "both active" to "only 1 active" (activity 2 completes)
            % When in state (i,j), activity 2 can complete with rate absRate2(j)
            % This moves to state "only 1 active, phase i"
            for i = 1:n1
                for j = 1:n2
                    bothIdx = (i-1)*n2 + j;
                    only1Idx = nBoth + i;
                    T_out(bothIdx, only1Idx) = T_out(bothIdx, only1Idx) + absRate2(j);
                end
            end

            % Transitions from "both active" to "only 2 active" (activity 1 completes)
            % When in state (i,j), activity 1 can complete with rate absRate1(i)
            % This moves to state "only 2 active, phase j"
            for i = 1:n1
                for j = 1:n2
                    bothIdx = (i-1)*n2 + j;
                    only2Idx = nBoth + nOnly1 + j;
                    T_out(bothIdx, only2Idx) = T_out(bothIdx, only2Idx) + absRate1(i);
                end
            end

            % Transitions within "only 1 active" states
            T_out(nBoth+1:nBoth+nOnly1, nBoth+1:nBoth+nOnly1) = T1;

            % Transitions within "only 2 active" states
            T_out(nBoth+nOnly1+1:end, nBoth+nOnly1+1:end) = T2;

            % Absorption from "only 1 active" and "only 2 active" states happens
            % when the remaining activity completes - this is handled by the
            % negative row sums (implicit absorption)

            % Initial probability: start in "both active" states
            % with probability alpha1(i) * alpha2(j) for state (i,j)
            alpha_out = zeros(1, nTotal);
            for i = 1:n1
                for j = 1:n2
                    bothIdx = (i-1)*n2 + j;
                    alpha_out(bothIdx) = alpha1(i) * alpha2(j);
                end
            end
        end

        function [alpha_out, T_out] = composeRepeat(alpha, T, count)
            % COMPOSEREPEAT Repeat a PH distribution count times (convolution)
            %
            % Equivalent to serial composition of the same PH count times.

            if count <= 0
                alpha_out = 1;
                T_out = -1e10;  % Immediate
                return;
            end

            if count == 1
                alpha_out = alpha;
                T_out = T;
                return;
            end

            % Use map_sum for efficient convolution if available
            % Otherwise, compose serially
            try
                % Convert to MAP format
                n = size(T, 1);
                e = ones(n, 1);
                D0 = T;
                D1 = -T * e * alpha;
                MAP = {D0, D1};

                % Use map_sum for count-fold convolution
                MAP_conv = map_sum(MAP, count);

                T_out = MAP_conv{1};
                D1_conv = MAP_conv{2};

                % Extract alpha from D1
                n_out = size(T_out, 1);
                e_out = ones(n_out, 1);
                absRate = -T_out * e_out;
                idx = find(absRate > GlobalConstants.FineTol, 1);
                if ~isempty(idx)
                    alpha_out = D1_conv(idx, :) / absRate(idx);
                else
                    alpha_out = zeros(1, n_out);
                    alpha_out(1) = 1;
                end
            catch
                % Fallback: serial composition
                alpha_out = alpha;
                T_out = T;
                for i = 2:count
                    [alpha_out, T_out] = Workflow.composeSerial(alpha_out, T_out, alpha, T);
                end
            end
        end

        % Static precedence factory methods (delegate to ActivityPrecedence)
        function ap = Serial(varargin)
            % SERIAL Create serial precedence
            %
            % AP = WORKFLOW.SERIAL(A1, A2, ...)

            % Convert WorkflowActivity to names
            args = cell(1, nargin);
            for i = 1:nargin
                if isa(varargin{i}, 'WorkflowActivity')
                    args{i} = varargin{i}.name;
                else
                    args{i} = varargin{i};
                end
            end
            ap = ActivityPrecedence.Serial(args{:});
        end

        function ap = AndFork(preAct, postActs)
            % ANDFORK Create AND-fork precedence
            %
            % AP = WORKFLOW.ANDFORK(PREACT, {POSTACT1, POSTACT2, ...})

            if isa(preAct, 'WorkflowActivity')
                preAct = preAct.name;
            end
            for a = 1:length(postActs)
                if isa(postActs{a}, 'WorkflowActivity')
                    postActs{a} = postActs{a}.name;
                end
            end
            ap = ActivityPrecedence.AndFork(preAct, postActs);
        end

        function ap = AndJoin(preActs, postAct, quorum)
            % ANDJOIN Create AND-join precedence
            %
            % AP = WORKFLOW.ANDJOIN({PREACT1, PREACT2, ...}, POSTACT)
            % AP = WORKFLOW.ANDJOIN({PREACT1, PREACT2, ...}, POSTACT, QUORUM)

            for a = 1:length(preActs)
                if isa(preActs{a}, 'WorkflowActivity')
                    preActs{a} = preActs{a}.name;
                end
            end
            if isa(postAct, 'WorkflowActivity')
                postAct = postAct.name;
            end
            if nargin < 3
                quorum = [];
            end
            ap = ActivityPrecedence.AndJoin(preActs, postAct, quorum);
        end

        function ap = OrFork(preAct, postActs, probs)
            % ORFORK Create OR-fork (probabilistic) precedence
            %
            % AP = WORKFLOW.ORFORK(PREACT, {POSTACT1, POSTACT2, ...}, [P1, P2, ...])

            if isa(preAct, 'WorkflowActivity')
                preAct = preAct.name;
            end
            for a = 1:length(postActs)
                if isa(postActs{a}, 'WorkflowActivity')
                    postActs{a} = postActs{a}.name;
                end
            end
            ap = ActivityPrecedence.OrFork(preAct, postActs, probs);
        end

        function ap = Xor(preAct, postActs, probs)
            % XOR Create XOR precedence (alias for OrFork)
            %
            % AP = WORKFLOW.XOR(PREACT, {POSTACT1, POSTACT2, ...}, [P1, P2, ...])

            ap = Workflow.OrFork(preAct, postActs, probs);
        end

        function ap = OrJoin(preActs, postAct)
            % ORJOIN Create OR-join precedence
            %
            % AP = WORKFLOW.ORJOIN({PREACT1, PREACT2, ...}, POSTACT)

            for a = 1:length(preActs)
                if isa(preActs{a}, 'WorkflowActivity')
                    preActs{a} = preActs{a}.name;
                end
            end
            if isa(postAct, 'WorkflowActivity')
                postAct = postAct.name;
            end
            ap = ActivityPrecedence.OrJoin(preActs, postAct);
        end

        function ap = Loop(preAct, postActs, counts)
            % LOOP Create loop precedence
            %
            % AP = WORKFLOW.LOOP(PREACT, {LOOPACT, ENDACT}, COUNT)
            % AP = WORKFLOW.LOOP(PREACT, LOOPACT, ENDACT, COUNT)

            if isa(preAct, 'WorkflowActivity')
                preAct = preAct.name;
            end
            if iscell(postActs)
                for a = 1:length(postActs)
                    if isa(postActs{a}, 'WorkflowActivity')
                        postActs{a} = postActs{a}.name;
                    end
                end
            elseif isa(postActs, 'WorkflowActivity')
                postActs = postActs.name;
            end
            ap = ActivityPrecedence.Loop(preAct, postActs, counts);
        end

        function wf = fromWfCommons(jsonFile, options)
            % FROMWFCOMMONS Load a workflow from a WfCommons JSON file.
            %
            % WF = WORKFLOW.FROMWFCOMMONS(JSONFILE)
            % WF = WORKFLOW.FROMWFCOMMONS(JSONFILE, OPTIONS)
            %
            % Loads a workflow trace from the WfCommons format
            % (https://github.com/wfcommons/workflow-schema) into a
            % LINE Workflow object for queueing analysis.
            %
            % Parameters:
            %   jsonFile - Path to WfCommons JSON file
            %   options  - Optional struct with:
            %     .distributionType - 'exp' (default), 'det', 'aph', 'hyperexp'
            %     .defaultSCV       - Default SCV for APH/HyperExp (default: 1.0)
            %     .defaultRuntime   - Default runtime when missing (default: 1.0)
            %     .useExecutionData - Use execution data if available (default: true)
            %     .storeMetadata    - Store WfCommons metadata (default: true)
            %
            % Returns:
            %   wf - Workflow object
            %
            % Example:
            %   wf = Workflow.fromWfCommons('montage-workflow.json');
            %   ph = wf.toPH();
            %   fprintf('Mean execution time: %.2f\n', ph.getMean());

            if nargin < 2
                options = struct();
            end
            wf = WfCommonsLoader.load(jsonFile, options);
        end
    end
end
