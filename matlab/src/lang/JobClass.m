classdef JobClass < NetworkElement
    % An abstract class for a collection of indistinguishable jobs
    %
    % Copyright (c) 2012-2024, Imperial College London
    % All rights reserved.
    
    properties
        priority;
        refstat; % reference station
        isrefclass; % is this a reference class within a chain?
        index; % node index
        type;
        completes; % true if passage through reference station is a completion        
    end
    
    methods (Hidden)
        %Constructor
        function self = JobClass(type, name)
            % SELF = JOBCLASS(TYPE, NAME)
            
            self@NetworkElement(name);
            self.priority = 0;
            self.refstat = Node('Unallocated');
            self.isrefclass = false;
            self.index = 1;
            self.type=type;
            self.completes = true;            
        end
        
        function self = setReferenceStation(self, source)
            % SELF = SETREFERENCESTATION(SOURCE)
            
            self.refstat = source;
        end
        
        function self = setReferenceClass(self, bool)
            % SELF = SETREFERENCECLASS(BOOL)            
            self.isrefclass = bool;
        end
        
        function boolIsa = isReferenceStation(self, node)
            % BOOLISA = ISREFERENCESTATION(NODE)
            
            boolIsa = strcmp(self.refstat.name,node.name);
        end
        
        function boolIsa = isReferenceClass(self)
            % BOOLISA = ISREFERENCECLASS()
            
            boolIsa = self.isrefclass;
        end
        
        %         function self = set.priority(self, priority)
        % SELF = SET.PRIORITY(PRIORITY)
        
        %             if ~(rem(priority,1) == 0 && priority >= 0)
        %                 line_error(mfilename,'Priority must be an integer.\n');
        %             end
        %             self.priority = priority;
        %         end
    end
    
    methods (Access=public)
        function ind = subsindex(self)
            % IND = SUBSINDEX()
            
            ind = double(self.index)-1; % 0 based
        end
        
        function summary(self)
            % SUMMARY()
            line_printf('Class (%s): <strong>%s</strong>',self.type,self.getName);
        end
    end
    
end
