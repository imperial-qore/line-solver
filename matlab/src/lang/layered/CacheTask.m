classdef CacheTask < Task
    % A software server in a LayeredNetwork.
    %
    % Supports multi-level cache hierarchies with configurable capacities
    % per level (e.g., L1, L2, L3 caches).
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    properties

        items;
        itemLevelCap;  % Scalar or array for multi-level cache capacities
        replacestrategy;
        retrieval = false;  % delayed-hit retrieval on the miss path

    end

    methods
        %public methods, including constructor

        %constructor
        function self = CacheTask(model, name, nitems, itemLevelCap, replStrat, multiplicity, scheduling)
            %self = CacheTask(model, name, nitems, itemLevelCap, replStrat, multiplicity, scheduling)
            %
            % itemLevelCap can be:
            %   - Scalar: Single-level cache with given capacity
            %   - Vector: Multi-level cache with capacity per level

            if ~exist('name','var')
                line_error(mfilename,'Constructor requires to specify at least a name.');
            end

            if nargin < 6
                multiplicity = 1;
            end

            if nargin < 7
                scheduling = SchedStrategy.FCFS;
            end
            self@Task(model, name, multiplicity, scheduling);

            self.items = nitems;

            % Validate and store itemLevelCap (scalar or array)
            if isscalar(itemLevelCap)
                self.itemLevelCap = itemLevelCap;
            elseif isvector(itemLevelCap)
                self.itemLevelCap = itemLevelCap(:)';  % Ensure row vector
            else
                line_error(mfilename, 'itemLevelCap must be scalar or vector');
            end

            self.replacestrategy = replStrat;
        end

        % Get item level capacity
        function cap = getItemLevelCap(self, level)
            % CAP = GETITEMLEVELCAP(SELF, LEVEL)
            %
            % Get capacity for a specific cache level or all levels.
            %
            % Args:
            %   level: (optional) Cache level index (1-based)
            %
            % Returns:
            %   cap: Capacity value(s) - scalar if level specified, array otherwise

            if nargin == 1
                % Return all levels
                cap = self.itemLevelCap;
            else
                % Return specific level
                if isscalar(self.itemLevelCap)
                    % Single-level cache
                    if level ~= 1
                        line_error(mfilename, 'Single-level cache only has level 1');
                    end
                    cap = self.itemLevelCap;
                else
                    % Multi-level cache
                    if level < 1 || level > length(self.itemLevelCap)
                        line_error(mfilename, sprintf('Level %d out of range [1, %d]', level, length(self.itemLevelCap)));
                    end
                    cap = self.itemLevelCap(level);
                end
            end
        end

        % Get total capacity across all levels
        function total = getTotalCapacity(self)
            % TOTAL = GETTOTALCAPACITY(SELF)
            %
            % Get the sum of capacities across all cache levels.
            %
            % Returns:
            %   total: Sum of all level capacities

            total = sum(self.itemLevelCap);
        end

        % Enable/disable delayed-hit retrieval on the miss path
        function self = setRetrieval(self, retrieval)
            % SETRETRIEVAL(SELF, RETRIEVAL)
            %
            % Enable a retrieval system with delayed-hit coalescing on the cache
            % miss path. When set, concurrent misses for the same item arriving
            % while a fetch (the miss-branch activity and its backend calls) is in
            % flight are parked and released together as delayed hits when the
            % fetch completes, instead of each triggering an independent fetch.
            if nargin < 2
                retrieval = true;
            end
            self.retrieval = logical(retrieval);
        end

        function tf = hasRetrieval(self)
            tf = self.retrieval;
        end

        % Set item level capacity
        function setItemLevelCap(self, caps)
            % SETITEMLEVELCAP(SELF, CAPS)
            %
            % Set cache level capacities.
            %
            % Args:
            %   caps: Scalar for single-level or vector for multi-level

            if isscalar(caps)
                self.itemLevelCap = caps;
            elseif isvector(caps)
                self.itemLevelCap = caps(:)';  % Ensure row vector
            else
                line_error(mfilename, 'caps must be scalar or vector');
            end
        end

        % Get number of cache levels
        function n = getNumberOfLevels(self)
            % N = GETNUMBEROFLEVELS(SELF)
            %
            % Get the number of cache levels.
            %
            % Returns:
            %   n: Number of levels (1 for single-level, >1 for multi-level)

            if isscalar(self.itemLevelCap)
                n = 1;
            else
                n = length(self.itemLevelCap);
            end
        end
    end
end