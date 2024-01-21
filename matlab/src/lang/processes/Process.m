classdef Process < Copyable
    % An abstract class for stochastic processes
    %
    % Copyright (c) 2012-2024, Imperial College London
    % All rights reserved.
    
    properties
        name
        params
    end
    
    methods %(Abstract) % implemented with errors for Octave compatibility
        function X = sample(self)
            % X = SAMPLE()
            
            % Sample a value from the inter-arrival time distribution
            line_error(mfilename,'An abstract method was called. The function needs to be overridden by a subclass.');
        end        
    end
    
    methods (Hidden)
        %Constructor
        function self = Process(name, numParam)
            % SELF = POINTPROCESS(NAME, NUMPARAM)
            
            self.name = name;
            self.params = cell(1,numParam);
            for i=1:numParam
                self.params{i}=struct('paramName','','paramValue',-1);
            end
        end
        
        function nParam = getNumParams(self)
            % NPARAM = GETNUMPARAMS()
            
            nParam = length(self.params);
        end
        
        function setParam(self, id, name, value,typeClass)
            % SETPARAM(ID, NAME, VALUE,TYPECLASS)
            
            self.params{id}.paramName=name;
            self.params{id}.paramValue=value;
        end
        
        
        function param = getParam(self,id)
            % PARAM = GETPARAM(SELF,ID)
            
            param = self.params{id};
        end
    end
    
end
