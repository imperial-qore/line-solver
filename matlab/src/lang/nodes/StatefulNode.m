classdef StatefulNode < Node
    % An abstract class for nodes that under some 
    % parametrizations can be stateful
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.
    
    properties
        state;
        statePrior;
        space;
    end
    
    methods(Hidden)
        %Constructor
        function self = StatefulNode(name)
            % SELF = STATEFULNODE(NAME)
            self@Node(name);
            self.statePrior = [];
            self.state = [];
            self.space = {};
        end
        
        function prior = getStatePrior(self)
            % PRIOR = GETSTATEPRIOR()
            % prior(j) = probability that the initial state is state(j)
            prior = self.statePrior;
        end
        
        function self = setStatePrior(self, prior)
            % SELF = SETSTATEPRIOR(PRIOR)
            % the prior is marginalized on the station state and thus
            % assumed independent of the priors for other stations
            self.statePrior = prior(:); % we do not normalize to allow the user to manually run a model for each point of an external prior
            if size(self.statePrior,1) ~= size(self.space,1)
                line_error(mfilename,'The prior probability vector must have the same rows of the station state space.');
            end
        end
        
        function self = setState(self, state)
            % SELF = SETSTATE(STATE)
            % state can be stacked in a matrix of states, a state space
            self.state = state;
        end
        
        function state = getState(self)
            % STATE = GETSTATE()
            state = self.state;
        end
        
        function state = getStateSpace(self)
            % STATE = GETSTATESPACE()
            state = self.space;
        end
        
        function self = setStateSpace(self, space)
            % SELF = SETSTATESPACE(SPACE)
            self.space = space;
            %            mState.getHash = memoize(@(x) self.State.getHash(x));
        end
        
        function self = resetStateSpace(self)
            % SELF = RESETSTATESPACE()
            self.space = {};
        end                
        
    end
    
    
end
