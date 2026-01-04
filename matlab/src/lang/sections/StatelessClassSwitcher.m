classdef StatelessClassSwitcher < ClassSwitcher
    % A class switch node basic on class switch probabilities
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    properties 
        csMatrix;
    end
    
    methods
        %Constructor
        function self = StatelessClassSwitcher(jobclasses, csMatrix)
            % SELF = STATELESSCLASSSWITCHER(CLASSES, CSMATRIX)
            
            self@ClassSwitcher(jobclasses, 'StatelessClassSwitcher');
            self.csMatrix = csMatrix;
            % this is slower than indexing the matrix, but it is a small
            % matrix anyway
            self.updateClasses(jobclasses);
            self.csFun = @(r,s,state,statep) csMatrix(r,s); % state parameter if present is ignored
        end
        
        function self = updateClasses(self, jobclasses)
            self.classes = jobclasses; % state parameter if present is ignored
        end
        
        function self = updateClassSwitch(self, csMatrix)
            self.csMatrix = csMatrix;
            self.csFun = @(r,s,state,statep) csMatrix(r,s); % state parameter if present is ignored
        end
    end
    
end
