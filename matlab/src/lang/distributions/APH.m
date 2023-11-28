classdef APH < MarkovianDistribution
    % An astract class for acyclic phase-type distributions
    %
    % Copyright (c) 2012-2023, Imperial College London
    % All rights reserved.

    methods
        function self = APH(alpha, T)
            % SELF = APH(ALPHA, T)

            % Abstract class constructor
            self@MarkovianDistribution('APH', 2);
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
            [alpha,T] = APHFrom3Moments([e1,e2,e3]);
            self.setParam(1, 'alpha', alpha);
            self.setParam(2, 'T', T);
            self.immediate = NaN;
        end

        function updateFromRawMoments(self,varargin)
            % UPDATE(SELF,VARARGIN)

            % Update parameters to match the first n moments
            % (n<=4)
            e1 = varargin{1};
            e2 = varargin{2};
            e3 = varargin{3};
            if length(varargin) > 3
                line_warning(mfilename,'Update in %s distributions can only handle 3 moments, ignoring higher-order moments.\n',class(self));
            end
            try
                [alpha,T] = APHFrom3Moments([e1,e2,e3]);
            catch
                m1 = e1;
                m2 = e2;
                m3 = e3;
                if m2<1.5*m1^2
                    m2 = 1.5*m1^2;
                else
                    m2 = m2;
                end

                scv = (m2/(m1^2))-1;

                if scv == 0.5
                    satisfy = 2;
                elseif 0.5<scv && scv<1
                    if 6*(m1^3)*scv<m3 && m3<3*(m1^3)*(3*scv-1+sqrt(2)*(1-scv)^(3/2))
                        satisfy = 1;
                        m2 = 3*(m1^3)*(3*scv-1+sqrt(2)*(1-scv)^(3/2));
                        m3 = 6*(m1^3)*scv;
                    elseif m3<min(3*(m1^3)*(3*scv-1+sqrt(2)*(1-scv)^(3/2)),6*(m1^3)*scv)
                        satisfy = 0;
                    elseif m3>max(6*(m1^3)*scv,3*(m1^3)*(3*scv-1+sqrt(2)*(1-scv)^(3/2)))
                        satisfy = 0;
                    else
                        satisfy = 1;
                        m3 = m3;
                    end
                elseif scv == 1
                    satisfy = 3;
                    % one dimensional exponential distribution with lambda = 1/m1
                elseif scv>1
                    satisfy = 1;
                    if m3<=(3/2)*m1^3*(1+scv)^2
                        m3 = (3/2)*m1^3*(1+scv)^2;
                    else
                        m3 = m3;
                    end
                else
                    satisfy = 1;
                    m3 = m3;
                end

                if satisfy == 2
                    c = 0.75*m1^4;
                    d = 0.5*m1^2;
                    b = 1.5*m1^3;
                    a = 0;
                    mu = (-b+6*m1*d+sqrt(a))/(b+sqrt(a));
                    lambda1 = (b-sqrt(a))/c;
                    lambda2 = (b+sqrt(a))/c;
                elseif satisfy == 1
                    c = 3*m2^2-2*m1*m3;
                    d = 2*m1^2-m2;
                    b = 3*m1*m2-m3;
                    a = b^2-6*c*d;

                    if c>0
                        mu = (-b+6*m1*d+sqrt(a))/(b+sqrt(a));
                        lambda1 = (b-sqrt(a))/c;
                        lambda2 = (b+sqrt(a))/c;
                    elseif c<0
                        mu = (b-6*m1*d+sqrt(a))/(-b+sqrt(a));
                        lambda1 = (b+sqrt(a))/c;
                        lambda2 = (b-sqrt(a))/c;
                    else
                        mu = 1/(2*scv);
                        lambda1 = 1/(scv*m1);
                        lambda2 = 2/m1;
                        % one dimensional exponential distribution with lambda = 1/m1
                    end
                else
                    mu = 1/(2*scv);
                    lambda1 = 1/(scv*m1);
                    lambda2 = 2/m1;
                end
                alpha = [mu,1-mu];
                T = [-lambda1,lambda1;0,-lambda2];
            end
            self.setParam(1, 'alpha', alpha);
            self.setParam(2, 'T', T);
            self.immediate = NaN;
        end

        function updateMean(self,MEAN)
            % UPDATEMEAN(SELF,MEAN)

            % Update parameters to match the given mean
            updateMean@MarkovianDistribution(self,MEAN);
        end

        function updateMeanAndSCV(self, MEAN, SCV)
            % UPDATEMEANANDSCV(MEAN, SCV)

            % Fit phase-type distribution with given mean and squared coefficient of
            % variation (SCV=variance/mean^2)
            e1 = MEAN;
            e2 = (1+SCV)*e1^2;
            [alpha,T] = APHFrom2Moments([e1,e2]);
            %e3=(3/2+1e-3)*e2^2/e1;
            %[alpha,T] = APHFrom3Moments([e1,e2,e3]);
            self.setParam(1, 'alpha', alpha);
            self.setParam(2, 'T', T);
            self.immediate = NaN;
        end

        function APH = getRepres(self)
            % APH = GETREPRESENTATION()

            % Return the renewal process associated to the distribution
            T = self.getSubgenerator;
            APH = {T,-T*ones(length(T),1)*self.getInitProb};
        end

        function Ft = evalCDF(self,varargin)
            alpha = self.getParam(1).paramValue;
            T = self.getParam(2).paramValue;
            e = ones(length(alpha),1);
            if isempty(varargin)
                mean = self.getMean;
                var = self.getVariance;
                sigma = sqrt(var);
                F = [];
                i = 1;
                if year(matlabRelease.Date)>2023 || strcmpi(matlabRelease.Release,"R2023b")
                    % use exmpmv
                    t = linspace(0,mean+10*sigma,500);
                    F(1:length(t)) = 1-sum(expmv(T',alpha',t)',2);
                else
                    for t = linspace(0,mean+10*sigma,500)
                        F(i) = 1-alpha*expm(T.*t)*e;
                        i = i+1;
                    end
                    t = linspace(0,mean+10*sigma,500);
                end
                Ft = [F',t'];
            else
                t = varargin{1};
                Ft = [];
                if year(matlabRelease.Date)>2023 || strcmpi(matlabRelease.Release,"R2023b")
                    % use exmpmv
                    F(1:length(t)) = 1-sum(expmv(T',alpha',t)',2);
                else
                    for i = 1:1:length(t)
                        Ft(i) = 1-alpha*expm(T.*t(i))*e;
                    end
                end
            end
        end
    end

    methods (Static)

        function aph = rand(order)
            map = aph_rand(order);
            aph = APH(map_pie(map),map{1});
        end

        function ex = fit(MEAN, SCV, SKEW)
            % EX = FIT(MEAN, SCV, SKEW)

            % Fit the distribution from first three standard moments (mean,
            % SCV, skewness)
            if MEAN <= GlobalConstants.FineTol
                ex = APH(1.0, [1]);
            else
                ex = APH(1.0, [1]);
                ex.update(MEAN, SCV, SKEW);
                ex.immediate = false;
            end
        end

        function ex = fitRawMoments(m1, m2, m3)

            % Fit the distribution from first three moments
            if m1 <= GlobalConstants.FineTol
                ex = Exp(Inf);
            else
                ex = APH(1.0, [1]);
                ex.updateFromRawMoments(m1, m2, m3);
                ex.immediate = false;
            end
        end

        function ex = fitCentral(MEAN, VAR, SKEW)
            % EX = FITCENTRAL(MEAN, VAR, SKEW)

            % Fit the distribution from first three central moments (mean,
            % variance, skewness)
            if MEAN <= GlobalConstants.FineTol
                ex = Exp(Inf);
            else
                ex = APH(1.0, [1]);
                ex.update(MEAN, VAR/MEAN^2, SKEW);
                ex.immediate = false;
            end
        end

        function ex = fitMeanAndSCV(MEAN, SCV)
            % EX = FITMEANANDSCV(MEAN, SCV)

            % Fit the distribution from first three central moments (mean,
            % variance, skewness)
            if MEAN <= GlobalConstants.FineTol
                ex = Exp(Inf);
            elseif SCV==1.0
                ex = Exp(1/MEAN);
            else
                ex = APH(1.0, [1]);
                ex.updateMeanAndSCV(MEAN, SCV);
                ex.immediate = false;
            end
        end
    end

end


