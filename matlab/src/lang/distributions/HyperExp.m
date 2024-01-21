classdef HyperExp < MarkovianDistribution
    % The hyper-exponential statistical distribution
    %
    % Copyright (c) 2012-2024, Imperial College London
    % All rights reserved.

    methods
        function self = HyperExp(varargin)
            % SELF = HYPEREXP(VARARGIN)
            % Constructs a two-phase exponential distribution from
            % probability of selecting phase 1 and the two phase rates
            self@MarkovianDistribution('HyperExp',nargin);
            if length(varargin)==2
                p = varargin{1};
                lambda = varargin{2};
                setParam(self, 1, 'p', p);
                setParam(self, 2, 'lambda1', lambda);
                setParam(self, 3, 'lambda2', lambda);
                self.obj = jline.lang.distributions.HyperExp(p, lambda, lambda);
            elseif length(varargin)==3
                p1 = varargin{1};
                lambda1 = varargin{2};
                lambda2 = varargin{3};
                setParam(self, 1, 'p', p1);
                setParam(self, 2, 'lambda1', lambda1);
                setParam(self, 3, 'lambda2', lambda2);
                self.obj = jline.lang.distributions.HyperExp(p1, lambda1, lambda2);
            end
        end

        function ex = getMean(self)
            % EX = GETMEAN()
            % Get distribution mean
            if self.getNumberOfPhases == 2
                p = self.getParam(1).paramValue;
                mu1 = self.getParam(2).paramValue;
                mu2 = self.getParam(3).paramValue;
                ex = p/mu1 + (1-p)/mu2;
            else
                MarkovianDistribution.getMean(self);
            end
        end

        function SCV = getSCV(self)
            % SCV = GETSCV()
            % Get the squared coefficient of variation of the distribution (SCV = variance / mean^2)
            if self.getNumberOfPhases == 2
                p = self.getParam(1).paramValue;
                mu1 = self.getParam(2).paramValue;
                mu2 = self.getParam(3).paramValue;
                SCV = (2*(p/mu1^2 + (1-p)/mu2^2) - (p/mu1 + (1-p)/mu2)^2)/(p/mu1 + (1-p)/mu2)^2;
            else
                MarkovianDistribution.getSCV(self);
            end
        end

        function Ft = evalCDF(self,t)
            % FT = EVALCDF(SELF,T)
            % Evaluate the cumulative distribution function at t
            % AT T

            if self.getNumberOfPhases == 2
                p = self.getParam(1).paramValue;
                mu1 = self.getParam(2).paramValue;
                mu2 = self.getParam(3).paramValue;
                Ft = p*(1-exp(-mu1*t))+(1-p)*(1-exp(-mu2*t));
            else
                MarkovianDistribution.evalCDF(self,t);
            end
        end

        function PH = getPH(self)
            % PH = GETPH()
            % Return the renewal process associated to the distribution
            p = self.getParam(1).paramValue;
            n = length(p);
            if n == 1
                mu1 = self.getParam(2).paramValue;
                mu2 = self.getParam(3).paramValue;
                PH={[-mu1,0;0,-mu2],[mu1*p,mu1*(1-p);mu2*p,mu2*(1-p)]};
            else
                mu = self.getParam(2).paramValue;
                D0 = -diag(mu);
                D1 = -D0*p(:)*ones(1,n);
                PH = {D0,D1};
            end
        end
    end

    methods(Static)
        function he = fit(MEAN, SCV, SKEW)
            % HE = FIT(MEAN, SCV, SKEW)
            % Fit distribution from first three standard moments
            e1 = MEAN;
            e2 = (1+SCV)*e1^2;
            e3 = -(2*e1^3-3*e1*e2-SKEW*(e2-e1^2)^(3/2));
            for n=[2] % larger ph requires more moments
                PH = kpcfit_ph_prony([e1,e2,e3],n);
                if map_isfeasible(PH)
                    p = map_pie(PH);
                    lambda = -diag(PH{1});
                    he = HyperExp(p,lambda);
                    return
                end
            end
            he = HyperExp.fitMeanAndSCV(MEAN,SCV);
            he.immediate = MEAN < GlobalConstants.CoarseTol;
        end

        function he = fitRate(RATE)
            % HE = FITRATE(RATE)
            % Fit distribution with given rate
            he = HyperExp(p, RATE, RATE);
            he.immediate = 1/RATE < GlobalConstants.CoarseTol;
        end

        function he = fitMean(MEAN)
            % HE = FITMEAN(MEAN)
            % Fit distribution with given mean
            he = HyperExp(p, 1/MEAN, 1/MEAN);
            he.immediate = MEAN < GlobalConstants.CoarseTol;
        end

        function he = fitMeanAndSCV(MEAN, SCV)
            % HE = FITMEANANDSCV(MEAN, SCV)
            % Fit distribution with given mean and squared coefficient of variation (SCV=variance/mean^2)
            [~,mu1,mu2,p] = map_hyperexp(MEAN,SCV);
            he = HyperExp(p, mu1, mu2);
            he.immediate = MEAN < GlobalConstants.CoarseTol;
        end

        function he = fitMeanAndSCVBalanced(MEAN, SCV)
            % HE = FITMEANANDSCV(MEAN, SCV)
            % Fit distribution with given mean and squared coefficient of
            % variation (SCV=variance/mean^2) and balanced means, i.e.,
            % p/mu1 = (1-p)/mu2
            mu1 =  -(2*(((SCV - 1)/(SCV + 1))^(1/2)/2 - 1/2))/MEAN;
            p= 1/2 - ((SCV - 1)/(SCV + 1))^(1/2)/2;
            if mu1<0 || p<0 || p>1
                p = ((SCV - 1)/(SCV + 1))^(1/2)/2 + 1/2;
                mu1 = (2*(((SCV - 1)/(SCV + 1))^(1/2)/2 + 1/2))/MEAN;
            end
            mu2=(1-p)/p*mu1;
            p = real(p);
            mu1 = real(mu1);
            mu2 = real(mu2);
            he = HyperExp(p, mu1, mu2);
            he.immediate = MEAN < GlobalConstants.CoarseTol;
        end
    end

end

