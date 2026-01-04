classdef MMPP2 < MarkovModulated
    % 2-phase Markov-Modulated Poisson Process - MMPP(2)
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    methods
        %Constructor
        function self = MMPP2(lambda0,lambda1,sigma0,sigma1)
            % SELF = MMPP2(LAMBDA0,LAMBDA1,SIGMA0,SIGMA1)

            self@MarkovModulated('MMPP2',4);
            setParam(self, 1, 'lambda0', lambda0);
            setParam(self, 2, 'lambda1', lambda1);
            setParam(self, 3, 'sigma0', sigma0);
            setParam(self, 4, 'sigma1', sigma1);
            lambda0 =  self.getParam(1).paramValue;
            lambda1 =  self.getParam(2).paramValue;
            sigma0 =  self.getParam(3).paramValue;
            sigma1 =  self.getParam(4).paramValue;
            D0 = [-sigma0-lambda0,sigma0;sigma1,-sigma1-lambda1];
            D1 = [lambda0,0;0,lambda1];
            self.process = {D0,D1};
        end

        function Di = D(self, i, wantSparse)
            % Di = D(i)
            if nargin<3
                wantSparse = false;
            end

            % Return representation matrix, e.g., D0=MAP.D(0)
            if wantSparse
                if nargin<2
                    Di=self.getProcess;
                else
                    Di=self.getProcess{i+1};
                end
            else
                if nargin<2
                    Di=self.getProcess;
                    for i=1:length(Di)
                        Di(i)=full(Di(i));
                    end
                else
                    Di=full(self.getProcess{i+1});
                end
            end
        end

        function meant = evalMeanT(self,t)
            % MEANT = EVALMEANT(SELF,T)

            lambda0 =  self.getParam(1).paramValue;
            lambda1 =  self.getParam(2).paramValue;
            sigma0 =  self.getParam(3).paramValue;
            sigma1 =  self.getParam(4).paramValue;
            lambda = (lambda0*sigma1 + lambda1*sigma0) / (sigma0+sigma1);
            meant = lambda * t;
        end

        function vart = evalVarT(self,t)
            % VART = EVALVART(SELF,T)

            % Evaluate the variance-time curve at t
            D0 = self.D(0);
            D1 = self.D(1);
            MAP = {D0,D1};
            vart = map_varcount(MAP, t);
        end

        function rate = getRate(self)
            rate = 1 / getMean(self);
        end

        function X = sample(self, n)
            % X = SAMPLE(N)
            if nargin<2 %~exist('n','var'),
                n = 1;
            end
            MAP = self.getProcess;
            if map_isfeasible(MAP)
                X = map_sample(MAP,n);
            else
                line_error(mfilename,'This process is infeasible (negative rates).');
            end
        end

        % inter-arrival time properties
        function mean = getMean(self)
            % MEAN = GETMEAN()

            lambda0 =  self.getParam(1).paramValue;
            lambda1 =  self.getParam(2).paramValue;
            sigma0 =  self.getParam(3).paramValue;
            sigma1 =  self.getParam(4).paramValue;
            lambda = (lambda0*sigma1 + lambda1*sigma0) / (sigma0+sigma1);
            mean = 1 / lambda;
        end

        function scv = getSCV(self)
            % SCV = GETSCV()

            lambda0 =  self.getParam(1).paramValue;
            lambda1 =  self.getParam(2).paramValue;
            sigma0 =  self.getParam(3).paramValue;
            sigma1 =  self.getParam(4).paramValue;
            scv = (2*lambda0^2*sigma0*sigma1 + lambda0*lambda1*sigma0^2 - 2*lambda0*lambda1*sigma0*sigma1 + lambda0*lambda1*sigma1^2 + lambda0*sigma0^2*sigma1 + 2*lambda0*sigma0*sigma1^2 + lambda0*sigma1^3 + 2*lambda1^2*sigma0*sigma1 + lambda1*sigma0^3 + 2*lambda1*sigma0^2*sigma1 + lambda1*sigma0*sigma1^2)/((sigma0 + sigma1)^2*(lambda0*lambda1 + lambda0*sigma1 + lambda1*sigma0));
        end

        function skew = getSkewness(self)
            skew = map_skew(self.getProcess());
        end

        function id = getIDC(self)
            % IDC = GETIDC() % INDEX OF DISPERSION FOR COUNTS

            lambda0 =  self.getParam(1).paramValue;
            lambda1 =  self.getParam(2).paramValue;
            sigma0 =  self.getParam(3).paramValue;
            sigma1 =  self.getParam(4).paramValue;
            id = 1 + 2*(lambda0-lambda1)^2*sigma0*sigma1/(sigma0+sigma1)^2/(lambda0*sigma1+lambda1*sigma0);
        end

        function id = getIDI(self)
            % IDI = GETIDI() % INDEX OF DISPERSION FOR INTERVALS

            id = self.getIDC(self); % only asymptotic for now
        end

        function n = getNumberOfPhases(~)
            % N = GETNUMBEROFPHASES()
            n = 2;
        end

        function acf = getACF(self, lags)
            % ACF = GETACF(self, lags)

            acf = map_acf(self.D,lags);
        end

        function [gamma2, gamma] = getACFDecay(self)
            % [gamma2, gamma] = GETACFDECAY(self)
            %
            % gamma2: asymptotic decay rate of acf
            % gamma: interpolated decay rate of acf

            gamma2 = map_gamma2(self.D);
            if nargout>1
                gamma = map_gamma(self.D);
            end
        end

        function bool = isImmediate(self)
            % BOOL = ISIMMEDIATE()

            bool = getMean(self) == 0;
        end


        function acf = evalACFT(self, lags, timescale)
            % ACF = EVALACFT(self, lags)
            %
            % Evaluate the autocorrelation in counts at timescale t

            acf = map_acfc(self.D, lags, timescale);
        end

        function mmdp = toMMDP(self)
            % TOMMDP Convert MMPP2 to MMDP2 (deterministic representation)
            %
            % Converts this Markov-Modulated Poisson Process to a
            % Markov-Modulated Deterministic Process suitable for fluid
            % queue analysis.
            %
            % @return mmdp MMDP2 representation of this MMPP2

            lambda0 = self.getParam(1).paramValue;
            lambda1 = self.getParam(2).paramValue;
            sigma0 = self.getParam(3).paramValue;
            sigma1 = self.getParam(4).paramValue;
            mmdp = MMDP2(lambda0, lambda1, sigma0, sigma1);
        end
    end

    methods (Static)

        function mmpp2 = rand()
            % MMPP2 = RAND()
            %
            % Generate random MAP using uniform random numbers

            D = mmpp_rand(2);
            mmpp2 = MMPP2(D{1}(1,2), D{1}(2,1), D{2}(1,1), D{2}(2,2));
        end

        function mmpp2 = fitCentralAndIDC(mean, var, skew, idc)
            if mean <= GlobalConstants.FineTol
                mmpp2 = Exp(Inf);
            else
                scv = var/mean^2;
                m = mmpp2_fit1(mean,scv,skew,idc);
                mmpp2 = MMPP2(m{2}(1,1),m{2}(2,2),m{1}(1,2),m{1}(2,1));
            end
        end

        function mmpp2 = fitCentralAndACFLag1(mean, var, skew, rho1)
            if mean <= GlobalConstants.FineTol
                mmpp2 = Exp(Inf);
            else
                scv = var/mean^2;
                m = mmpp2_fit4(mean,scv,skew,rho1);
                mmpp2 = MMPP2(m{2}(1,1),m{2}(2,2),m{1}(1,2),m{1}(2,1));
            end
        end

        function mmpp2 = fitCentralAndACFDecay(mean, var, skew, gamma2)
            if m1 <= GlobalConstants.FineTol
                mmpp2 = Exp(Inf);
            else
                scv = var/mean^2;
                m = mmpp2_fit2(mean,scv,skew,gamma2);
                mmpp2 = MMPP2(m{2}(1,1),m{2}(2,2),m{1}(1,2),m{1}(2,1));
            end
        end

        function mmpp2 = fitRawMomentsAndIDC(m1, m2, m3, idc)
            if m1 <= GlobalConstants.FineTol
                mmpp2 = Exp(Inf);
            else
                scv = (m2-m1^2)/m1^2;
                gamma2=-(scv-idc)/(-1+idc);
                m = mmpp2_fit3(m1,m2,m3,gamma2);
                mmpp2 = MMPP2(m{2}(1,1),m{2}(2,2),m{1}(1,2),m{1}(2,1));
            end
        end

        function mmpp2 = fitRawMomentsAndACFLag1(m1, m2, m3, rho1)
            if m1 <= GlobalConstants.FineTol
                mmpp2 = Exp(Inf);
            else
                scv = (m2-m1^2)/m1^2;
                rho0 = (1-1/scv)/2;
                gamma2 = rho1/rho0;
                m = mmpp2_fit3(m1,m2,m3,gamma2);
                mmpp2 = MMPP2(m{2}(1,1),m{2}(2,2),m{1}(1,2),m{1}(2,1));
            end
        end

        function mmpp2 = fitRawMomentsAndACFDecay(m1, m2, m3, gamma2)
            if m1 <= GlobalConstants.FineTol
                mmpp2 = Exp(Inf);
            else
                m = mmpp2_fit3(m1,m2,m3,gamma2);
                mmpp2 = MMPP2(m{2}(1,1),m{2}(2,2),m{1}(1,2),m{1}(2,1));
            end
        end

    end
end
