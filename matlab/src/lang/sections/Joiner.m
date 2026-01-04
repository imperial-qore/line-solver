classdef Joiner < InputSection
    % An input section that joins siblings tasks
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.
    
    properties
        joinStrategy;
        joinRequired;
        joinJobClasses;
    end
    
    methods
        %Constructor
        function self = Joiner(customerClasses)
            % SELF = JOINER(CUSTOMERCLASSES)
            
            self@InputSection('Joiner');
            self.joinJobClasses = {};
            initJoinJobClasses(self, customerClasses);
        end
    end
    
    methods (Access = 'public')
        
        function setRequired(self, customerClass, nJobs)
            % SETREQUIRED(CUSTOMERCLASS, NJOBS)
            
            self.joinRequired{customerClass.index} = nJobs;
        end
        
        function setStrategy(self, customerClass, joinStrat)
            % SETSTRATEGY(CUSTOMERCLASS, JOINSTRAT)
            
            self.joinJobClasses{customerClass.index} = customerClass;
            self.joinStrategy{customerClass.index} = joinStrat;
        end
        
        function initJoinJobClasses(self, customerClasses)
            % INITJOINJOBCLASSES(CUSTOMERCLASSES)
            
            for i = 1 : length(customerClasses)
                self.joinJobClasses{i} = customerClasses{i};
                self.joinRequired{i} = -1;
                self.joinStrategy{i} = JoinStrategy.STD;
            end
        end
    end
    
    methods(Access = protected)
        % Override copyElement method:
        function clone = copyElement(self)
            % CLONE = COPYELEMENT()
            
            % Make a shallow copy of all properties
            clone = copyElement@Copyable(self);
            % Make a deep copy of each object
            for i=1:length(self.joinJobClasses)
                clone.joinJobClasses{i} = self.joinJobClasses{i}.copy;
            end
        end
    end
    
end

