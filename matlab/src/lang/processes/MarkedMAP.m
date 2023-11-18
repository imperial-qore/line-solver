classdef MarkedMAP < MarkovModulated
    % Markov Modulated Arrival Process
    %
    % Copyright (c) 2012-2023, Imperial College London
    % All rights reserved.

    methods
        %Constructor
        function self = MarkedMAP(D,K)
            % SELF = MarkedMAP(D,K)
            %
            % LINE uses the M3A representation format
            % D={D0, D1, D11, D12, D13, ..., D1K}
            % K is the number of marking types

            self@MarkovModulated('MarkedMAP',K+2);

            if K == length(D)-1
                setParam(self, 1, 'D0', D{1});
                setParam(self, 2, 'D1', cellsum({D{2:end}}));
                for k=1:K
                    setParam(self, 2+k, sprintf('D1%d',k), D{1+k});
                end
            elseif K == length(D)-2
                setParam(self, 1, 'D0', D{1});
                setParam(self, 2, 'D1', D{2});
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
                Di = self.getRepres;
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
                    Di=self.getRepres;
                else
                    Di=self.getRepres{1+i+j};
                end
            else
                if nargin<2
                    Di=self.getRepres;
                    for i=1:length(Di)
                        Di(i)=full(Di(i));
                    end
                else
                    Di=full(self.getRepres{1+i+j});
                end
            end
        end

        function meant = evalMeanT(self, t)
            % MEANT = EVALMEANT(SELF,T)

            meant = self.getMAP.evalMeanT(t);
        end

        function vart = evalVarT(self, t)
            % VART = EVALVART(SELF,T)

            % Evaluate the variance-time curve at timescale t
            vart =  self.getMAP.evalVarT(t);
        end

        function acf = evalACFT(self, lags, timescale)
            % ACF = EVALACFT(self, lags)
            %
            % Evaluate the autocorrelation in counts at timescale t

            acf = self.getMAP.evalACFT(lags, timescale);
        end

        % inter-arrival time properties
        function mean = getMean(self)
            % MEAN = GETMEAN()

            mean = self.getMAP.getMean;
        end

        function scv = getSCV(self)
            % SCV = GETSCV()

            scv = self.getMAP.scv;
        end

        % inter-arrival time properties for each marking type
        function mean = getTypeMeans(self)
            % MEAN = GETTYPEMEAN()

            mean = 1 ./ map_lambda(self.D);
        end

        function acf = getACF(self, lags)
            % ACF = GETACF(self, lags)

            acf = self.getMAP.getACF(lags);
        end

        function map = getMAP(self)
            map = MAP(self.D(0), self.D(1));
        end

        function maps = getMAPs(self, types)
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

            [gamma2, gamma] = self.getMAP.getACFDecay();
        end

        function id = getIDC(self, t) % index of dispersion for counts
            % ID = GETIDC() % INDEX OF DISPERSION
            if nargin < 2
                id = self.getMAP.getIDC;
            else
                id = self.getMAP.getIDC(t);
            end
        end

        function id = getTypeIDC(self, t) % asymptotic index of dispersion for counts for each type
            % ID = GETTYPEIDC() % ASYMPTOTIC INDEX OF DISPERSION
            if nargin < 2
                id = mmap_count_idc(self.getRepres, GlobalConstants.Immediate);
            else
                id = mmap_count_idc(self.getRepres, t);
            end
        end

        function D = getRepresentation(self)
            % D = GETREPRESENTATION()

            % Retrocompatibility mapping
            D = self.getRepres();
        end

        function D = getRepres(self)
            % D = GETREPRESENTATION()
            D = getProcess(self);
        end

        function lam = getRate(self)
            % LAMBDA = GETRATE()

            lam = self.getMAP.getRate;
        end

        function MarkedMAP = getProcess(self)
            % MarkedMAP = GETPROCESS()
            K = self.getNumParams;
            MarkedMAP = cell(1,K-1);
            for k=1:K
                MarkedMAP{k} = self.getParam(k).paramValue;
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

        function bool = isImmmediate(self)
            % BOOL = ISIMMMEDIATE()

            bool = getMean(self) == 0;
        end

        function K = getNumberOfTypes(self)
            % K = GETNUMBEROFTYPES()
            % Number of marking types

            K = length(self.D)-2;
        end

        function MMAPr = toTimeReversed(self)
            MMAPr = MarkedMAP(mmap_timereverse(self.D), self.getNumberOfTypes);
        end

        function [X,C] = sample(self, n)
            % [X,C] = SAMPLE(N)
            if nargin<2 %~exist('n','var'),
                n = 1;
            end
            MarkedMAP = self.getRepres;
            if mmap_isfeasible(MarkedMAP)
                [X,C] = mmap_sample(MarkedMAP,n);
            else
                line_error(mfilename,'This process is infeasible (negative rates).');
            end
        end

    end

    methods (Static)

        function mmap = fit(trace, markings, order)
            % M3PP = FIT(TRACE, ORDER)
            T = m3afit_init(trace,markings);
            mmap = m3afit_auto(T,'NumStates',order);
        end


        function mmap = rand(order, nclasses)
            % MarkedMAP = RAND(ORDER,NCLASSES)
            %
            % Generate random MarkedMAP using uniform random numbers
            if nargin < 1
                order = 2;
            end
            mmap = MarkedMAP(mmap_rand(order,nclasses),nclasses);
        end
    end
end
