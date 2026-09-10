classdef CacheClassSwitcher < StatefulClassSwitcher
    % A class switcher section based on cache hits and misses
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.
    
    properties
        items;
        cap;
        levels;
        hitClass;
        missClass;
        retrievalClasses;          % [items x classes] -1 where no retrieval class is defined
        actualHitProb;
        actualMissProb;
        actualDelayedHitProb;      % delayed-hit fraction per arrival class (retrieval system; filled after solve)
        actualHitProbList;         % [classes x lists] per-list (per-level) hit fraction (filled after solve, where available)
        actualItemProb;            % [items x (lists+1)] per-item occupancy: col 1 = miss, cols 2..end = per-list (filled after solve, where available)
        actualListCost;            % [1 x lists] mean storage cost held by each list (filled after solve when item sizes are set)
        actualResidT;     % expected latency per arrival class (filled after solve)
        actualDelayedHitQLen;      % [1 x items] mean secondary requests waiting on the fetch of each item
        actualDelayedHitQLenFull;  % [1 x items] as above, including the request that triggered the fetch
    end

    methods
        function self = CacheClassSwitcher(classes, items, capacity, levels)
            % SELF = CACHECLASSSWITCHER(CLASSES, ITEMS, CAPACITY, LEVELS)

            if nargin<4 %~exist('levels','var')
                levels = 1;
            end
            self@StatefulClassSwitcher(classes, 'Cache');
            self.classes = classes;
            self.items = items;
            self.cap = capacity;
            self.levels = levels;
            self.csFun = @(r, s, state, statep) self.simpleHitMiss(r, s, state, statep); % do nothing by default
            self.hitClass = sparse([]);
            self.missClass = sparse([]);
            self.retrievalClasses = -ones(items, numel(classes));
            self.actualHitProb = sparse([]); % this field is filled after model solution
            self.actualMissProb = sparse([]); % this field is filled after model solution
            self.actualDelayedHitProb = sparse([]); % filled after model solution (retrieval system)
            self.actualHitProbList = sparse([]); % filled after model solution (per-list hit fractions)
            self.actualItemProb = sparse([]); % filled after model solution (per-item per-list occupancy)
            self.actualListCost = []; % filled after model solution (mean per-list storage cost)
            self.actualResidT = sparse([]); % filled after model solution
        end
    end
    
    methods
        function prob = simpleHitMiss(self, r, s, state, statep)
            % PROB = SIMPLEHITMISS(R, S, STATE, STATEP)
            
            if nargin <= 3
                state = []; %local server state
                statep = []; %local server state
            end
            if isempty(state) % get csMask (B matrix)
                % A retrieval class is an output class of the cache: it departs
                % unchanged (r==s) toward the retrieval system, so it must pass the
                % class-switch filter just like hitClass / missClass.
                isReceived = ~isempty(self.retrievalClasses) && any(self.retrievalClasses(:) == r);
                if (r==s  ... % hit and miss in the cache can depart in the same class
                        || ((r <= length(self.hitClass) && r <= length(self.missClass)) ... % since hitClass and missClass are sparse, check entry for r exists
                        && (s == self.hitClass(r) || s == self.missClass(r)))) ... % route out hit or miss classes
                        && (~isempty(find(r == self.hitClass)) || ~isempty(find(r == self.missClass)) || isReceived) % don't route out classes that are not hit/miss/retrieval-pending
                    prob = 1;
                else
                    prob = 0;
                end
            else
                % un-comment to restore class-switching in routing
                %                 if sum(state) == sum(statep)+1 % hit
                %                     if (r <= length(self.hitClass)) && s == self.hitClass(r)
                %                         prob = 1;
                %                     else
                %                         prob = 0;
                %                     end
                %                 else % miss
                %                     if (r <= length(self.missClass)) && s == self.missClass(r)
                %                         prob = 1;
                %                     else
                %                         prob = 0;
                %                     end
                %                 end
            end
        end
    end
    
end
