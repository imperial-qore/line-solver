classdef PointProcess < Process
    % An abstract class for stochastic point processes
    %
    % Copyright (c) 2012-2024, Imperial College London
    % All rights reserved.
    
    properties
    end
    
    methods %(Abstract) % implemented with errors for Octave compatibility

        function ex = getMean(self)
            % EX = GETMEAN()
            
            % Returns the mean of the inter-arrival times
            line_error(mfilename,sprintf('An abstract method was called. The function needs to be overridden by a subclass.'));
        end
        
        function SCV = getSCV(self)
            % SCV = GETSCV()
            
            % Get squared coefficient of variation of the interarrival times (SCV = variance / mean^2)
            line_error(mfilename,sprintf('An abstract method was called. The function needs to be overridden by a subclass.'));
        end
        
        function SKEW = getSkewness(self)
            % SKEW = GETSKEWNESS()
            
            % Get skewness of the interarrival times
            line_error(mfilename,sprintf('An abstract method was called. The function needs to be overridden by a subclass.'));
        end
        
        function IDC = getIDC(self)
            % IDC = GETIDC()
            
            % Return the asymptotic index of dispersion for counts
            line_error(mfilename,sprintf('An abstract method was called. The function needs to be overridden by a subclass.'));
        end
        
        function IDI = getIDI(self)
            % IDI = GETIDI()
            
            % Return the asymptotic index of dispersion for intervals
            line_error(mfilename,sprintf('An abstract method was called. The function needs to be overridden by a subclass.'));
        end
        
        function lambda = getRate(self)
            % LAMBDA = GETRATE()
            % Return the inter-arrival rate
            line_error(mfilename,sprintf('An abstract method was called. The function needs to be overridden by a subclass.'));
        end
        
        function vart = evalVarT(self,t)
            % VART = EVALVART(SELF,T)
            
            % Evaluate the variance-time curve at t
            line_error(mfilename,sprintf('An abstract method was called. The function needs to be overridden by a subclass.'));
        end
    end
    
    methods (Hidden)
        %Constructor
        function self = PointProcess(name, numParam)
            % SELF = POINTPROCESS(NAME, NUMPARAM)
            self@Process(name, numParam);
        end
               
        function bool = isDisabled(self)
            % BOOL = ISDISABLED()
            
            bool = any(cellfun(@(c) ifthenelse(isstruct(c),false,isnan(c.paramValue)), self.params));
        end
        
        function bool = isImmediate(self)
            % BOOL = ISIMMEDIATE()
            
            bool = getMean(self) == 0;
        end
        
    end
    
end
