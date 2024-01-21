classdef RenewalProcess < SinkProcess
    % An abstract class for renewal point processes.
    %
    % Copyright (c) 2012-2024, Imperial College London
    % All rights reserved.
    
    properties
        distrib;
    end
    
    methods (Hidden)
        %Constructor
        function self = RenewalProcess(distrib)
            % SELF = RENEWALPROCESS(DISTRIB)
            
            self@SinkProcess('RenewalProcess',1);
            setParam(self, 1, 'distrib', distrib);
        end
    end
    
    methods
        function mt = evalMeanT(self, t, steps)
            % MT = EVALMEANT(T, STEPS)
            
            distrib =  self.getParam(1).paramValue;
            if nargin<3 %~exist('steps','var')
                steps = t/distrib.getMean()*20;
            end
            tmax = max(t);
            tmax = tmax + tmax/steps; % add one step for evalCDFInterval
            trange = unique([0:tmax/steps:tmax,t(:)']);
            
            mt = zeros(1,length(trange)-1);
            for ti=1:(length(trange)-1)
                mt(ti) = distrib.evalCDF(trange(ti));
                for si=1:(ti-1)
                    mt(ti) = mt(ti) +  mt(ti-si)*(distrib.evalCDFInterval(trange(si),trange(si+1)));
                end
            end
        end
        
        function vart = evalVarT(self,t)
            % VART = EVALVART(SELF,T)
            
            % Evaluate the variance-time curve at t
        end
        
        function mean = getMean(self)
            % MEAN = GETMEAN()
            
            distrib =  self.getParam(1).paramValue;
            mean = distrib.getMean;
        end
    end
end

