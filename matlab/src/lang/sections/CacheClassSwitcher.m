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
        actualHitProb;
        actualMissProb;
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
            self.actualHitProb = sparse([]); % this field is filled after model solution
            self.actualMissProb = sparse([]); % this field is filled after model solution
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
                if (r==s  ... % hit and miss in the cache can depart in the same class
                        || ((r <= length(self.hitClass) && r <= length(self.missClass)) ... % since hitClass and missClass are sparse, check entry for r exists
                        && (s == self.hitClass(r) || s == self.missClass(r)))) ... % route out hit or miss classes
                        && (~isempty(find(r == self.hitClass)) || ~isempty(find(r == self.missClass))) % don't route out classes that are not hit or miss
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
