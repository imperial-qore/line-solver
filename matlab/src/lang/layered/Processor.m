classdef Processor  < Host
    % A hardware server in a LayeredNetwork.
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.
    
    
    methods
        %public methods, including constructor
        
        %constructor
        function self = Processor(model, name, multiplicity, scheduling, quantum, speedFactor)
            % OBJ = PROCESSOR(MODEL, NAME, MULTIPLICITY, SCHEDULING, QUANTUM, SPEEDFACTOR)
            name = char(name);

            if nargin<2 %~exist('name','var')
                line_error(mfilename,'Constructor requires to specify at least a name.');
            end            
            if nargin<3 %~exist('multiplicity','var')
                multiplicity = 1;
            end
            if nargin<4 %~exist('scheduling','var')
                scheduling = SchedStrategy.PS;
            end
            if nargin<5 %~exist('quantum','var')
                quantum = 0.001;
            end
            if nargin<6 %~exist('speedFactor','var')
                speedFactor = 1;
            end            
            self@Host(model, name, multiplicity, scheduling, quantum, speedFactor)
        end
        
    end
    
end
