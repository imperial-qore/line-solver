classdef PH < MarkovianDistribution
    % An astract class for phase-type distributions
    %
    % Copyright (c) 2012-2023, Imperial College London
    % All rights reserved.
        
    methods
        function self = PH(alpha, T)
            % SELF = PH(ALPHA, T)
            
            % Abstract class constructor
            self@MarkovianDistribution('PH', 2);
            self.setParam(1, 'alpha', alpha);
            self.setParam(2, 'T', T);
        end
    end
    
    methods
        function alpha = getInitProb(self)
            % ALPHA = GETINITPROB()
            
            % Get vector of initial probabilities
            alpha = self.getParam(1).paramValue(:);
            alpha = reshape(alpha,1,length(alpha));
        end
        
        function T = getSubgenerator(self)
            % T = GETSUBGENERATOR()            
            
            % Get subgenerator
            T = self.getParam(2).paramValue;
        end               
        
        function X = sample(self, n)
            % X = SAMPLE(N)
            
            % Get n samples from the distribution
            if nargin<2 %~exist('n','var'), 
                n = 1; 
            end
            X = map_sample(self.getRepres,n);
        end
    end
    
    methods
        function update(self,varargin)
            % UPDATE(SELF,VARARGIN)
            
            % Update parameters to match the first n central moments
            % (n<=4)
            MEAN = varargin{1};
            SCV = varargin{2};
            SKEW = varargin{3};
            if length(varargin) > 3
                line_warning(mfilename,'Update in %s distributions can only handle 3 moments, ignoring higher-order moments.\n',class(self));
            end
            e1 = MEAN;
            e2 = (1+SCV)*e1^2;
            e3 = -(2*e1^3-3*e1*e2-SKEW*(e2-e1^2)^(3/2));
            
            if SCV<1
                [alpha,T] = APHFrom3Moments([e1,e2,e3]);
            else
                options = kpcfit_ph_options([e1,e2,e3],'MinNumStates',2,'MaxNumStates',2,'Verbose',false);
            end
            ph = kpcfit_ph_auto([e1,e2,e3],options);
            [alpha,T]=map2ph(ph{1,1});
            
            self.setParam(1, 'alpha', alpha);
            self.setParam(2, 'T', T);
        end
        
        function updateMean(self,MEAN)
            % UPDATEMEAN(SELF,MEAN)
            
            % Update parameters to match the given mean
            ph = self.getRepres;
            ph = map_scale(ph,MEAN);
            self.setParam(1, 'alpha', map_pie(ph));
            self.setParam(2, 'T', ph{1});
        end
        
        function updateMeanAndSCV(self, MEAN, SCV)
            % UPDATEMEANANDSCV(MEAN, SCV)
            
            % Fit phase-type distribution with given mean and squared coefficient of
            % variation (SCV=variance/mean^2)
            e1 = MEAN;
            e2 = (1+SCV)*e1^2;
            [alpha,T] = APHFrom2Moments([e1,e2]);
            %options = kpcfit_ph_options([e1,e2],'MinNumStates',2,'MaxNumStates',2,'Verbose',false);
            %ph = kpcfit_ph_auto([e1,e2],options);
            %[alpha,T]=map2ph(ph{1,1});            
            self.setParam(1, 'alpha', alpha);
            self.setParam(2, 'T', T);
        end
        
        function ph = getRepres(self)
            % PH = GETREPRESENTATION()
            
            % Return the renewal process associated to the distribution
            T = self.getSubgenerator;
            ph = {T,-T*ones(length(T),1)*self.getInitProb};
        end
        
    end
    
    methods (Static)
        function ex = fit(MEAN, SCV, SKEW)
            % EX = FIT(MEAN, SCV, SKEW)
            %
            % Fit the distribution from first three standard moments (mean,
            % SCV, skewness)
            ex = PH(1.0, [1]);
            ex.update(MEAN, SCV, SKEW);
        end
        
        function ex = fitCentral(MEAN, VAR, SKEW)
            % EX = FITCENTRAL(MEAN, VAR, SKEW)
            %
            % Fit the distribution from first three central moments (mean,
            % variance, skewness)
            ex = PH(1.0, [1]);
            ex.update(MEAN, VAR/MEAN^2, SKEW);
        end
        
        function ex = fitMeanAndSCV(MEAN, SCV)
            % EX = FITMEANANDSCV(MEAN, SCV)
            %
            % Fit the distribution from first three central moments (mean,
            % variance, skewness)
            ex = PH(1.0, [1]);
            ex.updateMeanAndSCV(MEAN, SCV);
        end
        
        function ex = fromRepresentation(dep)
            % EX = FROMREPRESENTATION(DEP)
            %
            % Build a phase-type distribution from a cell array {D0,D1}
            
            ex = PH(map_pie(dep),dep{1});
        end
    end
    
end
