classdef MMPP2 < MarkovModulated
    % 2-phase Markov-Modulated Poisson Process - MMPP(2)
    %
    % Copyright (c) 2012-2024, Imperial College London
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
        end

        function Di = D(self, i, wantSparse)
            % Di = D(i)
            if nargin<3
                wantSparse = false;
            end

            % Return representation matrix, e.g., D0=MAP.D(0)
            if wantSparse
                if nargin<2
                    Di=self.getRepres;
                else
                    Di=self.getRepres{i+1};
                end
            else
                if nargin<2
                    Di=self.getRepres;
                    for i=1:length(Di)
                        Di(i)=full(Di(i));
                    end
                else
                    Di=full(self.getRepres{i+1});
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
            lambda0 =  self.getParam(1).paramValue;
            lambda1 =  self.getParam(2).paramValue;
            sigma0 =  self.getParam(3).paramValue;
            sigma1 =  self.getParam(4).paramValue;
            MAP = getRepres;
            D0 = MAP{1};
            D1 = MAP{2};
            e = [1;1];
            pie = map_pie(MAP);
            I = eye(2);
            lambda = (lambda0*sigma1 + lambda1*sigma0) / (sigma0+sigma1);
            vart = lambda*t;
            vart = vart + 2*t*(lambda^2-pie*D1*inv(D0+D1+pie*e)*D1*e);
            vart = vart + 2*pi*D1*(expm((D0+D1)*t)-I)*inv(D0+D1+pie*e)^2*D1*e;
        end

        function rate = getRate(self)
            rate = 1 / getMean(self);
        end

        function X = sample(self, n)
            % X = SAMPLE(N)
            if nargin<2 %~exist('n','var'),
                n = 1;
            end
            MAP = self.getRepres;
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

        function MAP = getRepres(self)
            MAP = getProcess(self);
        end

        function PH = getPH(self)
            % PH = GETPH()
            PH = getProcess(self);
        end

        function MAP = getProcess(self)
            % MAP = GETPROCESS()

            lambda0 =  self.getParam(1).paramValue;
            lambda1 =  self.getParam(2).paramValue;
            sigma0 =  self.getParam(3).paramValue;
            sigma1 =  self.getParam(4).paramValue;
            D0 = [-sigma0-lambda0,sigma0;sigma1,-sigma1-lambda1];
            D1 = [lambda0,0;0,lambda1];
            MAP = {D0,D1};
        end

        function n = getNumberOfPhases(self)
            % N = GETNUMBEROFPHASES()
            n = 2;
        end

        function mu = getMu(self)
            % MU = GETMU()
            D = getProcess(self);
            mu = sum(D{2},2); % sum D1 rows / diag -D0
        end

        function phi = getPhi(self)
            % PHI = GETPHI()
            % Return the exit vector of the underlying PH
            D = getProcess(self);
            phi = sum(D{2},2) ./ diag(-D{1}); % sum D1 rows / diag -D0
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

        function bool = isImmmediate(self)
            % BOOL = ISIMMMEDIATE()

            bool = getMean(self) == 0;
        end


        function acf = evalACFT(self, lags, timescale)
            % ACF = EVALACFT(self, lags)
            %
            % Evaluate the autocorrelation in counts at timescale t

            acf = map_acfc(self.D, lags, timescale);
        end

        function map = toMAP(self)
            map = MAP(self.D);
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
