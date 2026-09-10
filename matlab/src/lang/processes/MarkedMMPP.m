classdef MarkedMMPP < MarkovModulated
    % Marked Markov-Modulated Poisson Process
    % Also referred to as a M3PP model
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    methods
        %Constructor
        function self = MarkedMMPP(D,K)
            % SELF = MarkedMMPP(D,K)
            %
            % LINE uses the M3A representation format
            % D={D0, D1, D11, D12, D13, ..., D1K}
            % K is the number of marking types

            self@MarkovModulated('MarkedMMPP',K+2);

            if K == length(D)-1
                setParam(self, 1, 'D0', D{1});
                D1 = cellsum({D{2:end}});
                if any(D1-diag(diag(D1)))
                    line_error(mfilename,'Non-diagonal elements in the D1k matrices, not a Marked MMPP.');
                end
                setParam(self, 2, 'D1', D1);
                for k=1:K
                    setParam(self, 2+k, sprintf('D1%d',k), D{1+k});
                end
            elseif K == length(D)-2
                setParam(self, 1, 'D0', D{1});
                D1 = D{2};
                if any(D1-diag(diag(D1)))
                    line_error(mfilename,'Non-diagonal  elements in the D1k matrices, not a Marked MMPP.');
                end
                setParam(self, 2, 'D1', D1);
                for k=1:K
                    setParam(self, 2+k, sprintf('D1%d',k), D{2+k});
                end
            else
                line_error(mfilename,'Inconsistency between the number of classes and the number of supplied D matrices.')
            end
        end

        function Di = D(self, i, j, wantSparse)
            % Di = D(i)
            if nargin == 1
                Di = self.getProcess;
                return
            end
            if nargin<4
                wantSparse = false;
            end
            if nargin<3
                if i<=1
                    j = 0;
                else
                    line_error(mfilename,'Invalid D matrix indexes');
                end
            end

            % Return representation matrix, e.g., D0=MAP.D(0)
            if wantSparse
                if nargin<2
                    Di=self.getProcess;
                else
                    Di=self.getProcess{1+i+j};
                end
            % else
                if nargin<2
                    Di=self.getProcess;
                    for i=1:length(Di)
                        Di(i)=full(Di(i));
                    end
                else
                    Di=full(self.getProcess{1+i+j});
                end
            end
        end

        function meant = evalMeanT(self, t)
            % MEANT = EVALMEANT(SELF,T)

            meant = self.toMAP.evalMeanT(t);
        end

        function vart = evalVarT(self, t)
            % VART = EVALVART(SELF,T)

            % Evaluate the variance-time curve at timescale t
            vart =  self.toMAP.evalVarT(t);
        end

        function acf = evalACFT(self, lags, timescale)
            % ACF = EVALACFT(self, lags)
            %
            % Evaluate the autocorrelation in counts at timescale t

            acf = self.toMAP.evalACFT(lags, timescale);
        end

        % inter-arrival time properties
        function mean = getMean(self)
            % MEAN = GETMEAN()

            mean = self.toMAP.getMean;
        end

        function scv = getSCV(self)
            % SCV = GETSCV()

            scv = self.toMAP.scv;
        end

        % inter-arrival time properties for each marking type
        function mean = getMarkedMeans(self)
            % MEAN = GETMARKEDMEAN()

            mean = 1 ./ mmap_lambda(self.D);
        end

        function acf = getACF(self, lags)
            % ACF = GETACF(self, lags)

            acf = self.toMAP.getACF(lags);
        end

        function map = toMAP(self)
            map = MAP(self.D(0), self.D(1));
        end

        function maps = toMAPs(self, types)
            if nargin<2
                types = 1:self.getNumberOfTypes;
            end
            maps = {};
            for k=types
                Df = self.D(0);
                Df = Df + (self.D(1) -  self.D(1,k));
                maps{end+1} = MAP(Df, self.D(1,k));
            end
        end

        function [gamma2, gamma] = getACFDecay(self)
            % [gamma2, gamma] = GETACFDECAY(self)
            %
            % gamma2: asymptotic decay rate of acf
            % gamma: interpolated decay rate of acf

            [gamma2, gamma] = self.toMAP.getACFDecay();
        end

        function id = getIDC(self, t) % index of dispersion for counts
            % ID = GETIDC() % INDEX OF DISPERSION
            if nargin < 2
                id = self.toMAP.getIDC;
            else
                id = self.toMAP.getIDC(t);
            end
        end

        function id = getMarkedIDC(self, t) % asymptotic index of dispersion for counts for each type
            % ID = GETMARKEDIDC() % ASYMPTOTIC INDEX OF DISPERSION
            if nargin < 2
                id = mmpp_count_idc(self.getProcess, GlobalConstants.Immediate);
            else
                id = mmpp_count_idc(self.getProcess, t);
            end
        end

        function lam = getRate(self)
            % LAMBDA = GETRATE()

            lam = self.toMAP.getRate;
        end

        function MarkedMMPP = getProcess(self)
            % MarkedMMPP = GETPROCESS()
            K = self.getNumParams;
            MarkedMMPP = cell(1,K-1);
            for k=1:K
                MarkedMMPP{k} = self.getParam(k).paramValue;
            end
        end

        function n = getNumberOfPhases(self)
            % N = GETNUMBEROFMAPASES()
            D0 =  self.getParam(1).paramValue;
            n = length(D0);
        end

        function mu = getMu(self)
            % MU = GETMU()
            % Aggregate departure rate from each state
            mu = sum(self.D(1),2); % sum D1 rows / diag -D0
        end

        function phi = getPhi(self)
            % MAPI = GETMAPI()
            % Return the exit vector of the underlying MAP
            phi = sum(self.D(1),2) ./ diag(-self.D(0)); % sum D1 rows / diag -D0
        end

        function K = getNumberOfTypes(self)
            % K = GETNUMBEROFTYPES()
            % Number of marking types

            K = length(self.D)-2;
        end

        function MMAPr = toTimeReversed(self)
            MMAPr = MarkedMMPP(mmap_timereverse(self.D), self.getNumberOfTypes);
        end

        function [X,C] = sample(self, n)
            % [X,C] = SAMPLE(N)
            if nargin<2 %~exist('n','var'),
                n = 1;
            end
            MarkedMMPP = self.getProcess;
            if mmpp_isfeasible(MarkedMMPP)
                [X,C] = mmpp_sample(MarkedMMPP,n);
            else
                line_error(mfilename,'This process is infeasible (negative rates).');
            end
        end

    end

    methods (Static)

        function m3pp = fit(trace, markings, order)
            % M3PP = FIT(TRACE, ORDER)
            T = m3afit_init(trace,markings);
            m3pp = m3afit_auto(T,'NumStates',order,'Method',1);
        end

        function m3pp = rand(order, nclasses)
            % M3PP = RAND(ORDER,NCLASSES)
            %
            % Generate random MarkedMMPP using uniform random numbers
            if nargin < 1
                order = 2;
            end
            m3pp = MarkedMMPP(m3pp_rand(order,nclasses),nclasses);
        end
    end
end
