classdef CacheTask < Task 
    % A software server in a LayeredNetwork.
    %
    % Copyright (c) 2012-2023, Imperial College London
    % All rights reserved.
    
    properties
        
        items;
        itemLevelCap;
        replacementPolicy           
        
    end

    methods
        %public methods, including constructor
        
        %constructor
        function self = CacheTask(model, name, nitems, itemLevelCap, replPolicy, multiplicity, scheduling)
            %self = CacheTask(model, name, nitems, itemLevelCap, replPolicy, multiplicity, scheduling)

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
            self.itemLevelCap = itemLevelCap; % item capacity            
            self.replacementPolicy = ReplacementStrategy.toId(replPolicy);
        end
    end
end