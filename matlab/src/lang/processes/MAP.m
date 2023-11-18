classdef MAP < MarkovModulated
    % Markovian Arrival Process
    %
    % Copyright (c) 2012-2023, Imperial College London
    % All rights reserved.

    methods
        %Constructor
        function self = MAP(D0,D1)
            % SELF = MAP(D0,D1)

            self@MarkovModulated('MAP',2);
            if nargin < 2 && iscell(D0)
                M = D0;
                D0 = M{1};
                D1 = M{2};
            end
            setParam(self, 1, 'D0', D0);
            setParam(self, 2, 'D1', D1);

            if ~map_isfeasible(self.D)
                line_warning(mfilename,'MAP is infeasible.\n');
            end
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

        function meant = evalMeanT(self, t)
            % MEANT = EVALMEANT(SELF,T)

            meant = map_count_mean(self.D, t);
        end

        function vart = evalVarT(self, t)
            % VART = EVALVART(SELF,T)

            % Evaluate the variance-time curve at timescale t
            vart =  map_count_var(self.D, t);
        end

        function acf = evalACFT(self, lags, timescale)
            % ACF = EVALACFT(self, lags)
            %
            % Evaluate the autocorrelation in counts at timescale t

            acf = map_acfc(self.D, lags, timescale);
        end

        function n = getNumberOfPhases(self)
            % N = GETNUMBEROFPHASES()
            n = 2;
        end

        % inter-arrival time properties
        function mean = getMean(self)
            % MEAN = GETMEAN()

            mean = map_mean(self.D);
        end

        function scv = getSCV(self)
            % SCV = GETSCV()

            scv = map_scv(self.D);
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

        function id = getIDC(self, t) % index of dispersion for counts
            % IDC = GETIDC() % ASYMPTOTIC INDEX OF DISPERSION

            if nargin < 2
                id = map_idc(self.D);
            else
                id = map_count_var(self.D,t) / map_count_mean(self.D,t);
            end
        end

        function MAP = getRepresentation(self)
            % MAP = GETREPRESENTATION()

            % Retrocompatibility mapping
            MAP = self.getRepres();
        end

        function MAP = getRepres(self)
            % MAP = GETREPRESENTATION()
            MAP = getProcess(self);
        end

        function lam = getRate(self)
            % MAP = GETRATE()

            lam = map_lambda(self.D);
        end

        function MAP = getProcess(self)
            % MAP = GETPROCESS()

            D0 =  self.getParam(1).paramValue;
            D1 =  self.getParam(2).paramValue;
            MAP = {D0,D1};
        end

        function bool = isImmmediate(self)
            % BOOL = ISIMMMEDIATE()

            bool = getMean(self) == 0;
        end

        function mapr = toTimeReversed(self)
            mapr = MAP(map_timereverse(self.D));
        end

        function ph = toPH(self)
            ph = PH(self.getEmbProb,self.D(0));
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

        function self = updateMean(self,MEAN)
            % UPDATEMEAN(SELF,MEAN)
            % Update parameters to match the given mean
            newMAP = map_scale(self.D,MEAN);
            self.params{1}.paramValue = newMAP{1};
            self.params{2}.paramValue = newMAP{2};
        end

        function bool = isImmediate(self)
            bool = self.getMean < GlobalConstants.FineTol;
        end

    end

    methods (Static)

        function map = rand(order)
            % MAP = RAND(ORDER)
            %
            % Generate random MAP using uniform random numbers
            if nargin < 1
                order = 2;
            end
            map = MAP(map_rand(order));
        end

        function map = randn(order, mu, sigma)
            % MAP = RANDN(ORDER, MU, SIGMA)
            %
            % Generate random MAP using specified Gaussian parameter and
            % taking the absolute value of the resulting values
            map = MAP(map_randn(order, mu, sigma));
        end

        function [map] = fit(trace, order)
            % MAP = FIT(TRACE, ORDER)
            if nargin<2
                order = 16;
            end
            map = MAP.kpcfit(trace, order);
        end

        function [map, fac, fbc, kpcMAPs] = kpcfit(trace, order, ac_runs_max, bc_runs_max)
            % MAP = KPCFIT(TRACE, ORDER, AC_RUNS_MAX, BC_RUNS_MAX)

            % ORDER must be a power of 2.
            % MAP.kpcfit(S,8,5,2) % fit a MAP(8) with 5 autocorrelation
            % trial runs and 2 bicorrelation trial runs

            if nargin<3
                ac_runs_max = 10;
            end
            if nargin<4
                bc_runs_max = 1;
            end

            T = kpcfit_init(trace);
            [D, fac, fbc, kpcMAPs] = kpcfit_auto(T, 'NumMAPs',ceil(log2(order)), 'MaxRunsAC', ac_runs_max, 'MaxRunsBC',bc_runs_max);
            map = MAP(D);
        end


        function [map] = anfit(trace, order, iter_max, iter_tol)
            if order<4
                error('error: requested %d states, but anfit can return only maps with a minimum of4 states', order);
            end
            if mod(order, 2) > 0
                error('error: requested %d states, but anfit can return only a number of states that is a power of 2', order);
            end
            if nargin<3
                iter_max=100;
            end
            if nargin<4
                iter_tol=GlobalConstants.FineTol;
            end
            T = kpcfit_init(trace);
            % ls - mean arrival rate
            ls = 1/ T.E(1);
            % rho - lag-1 autocorrelation of the counting process
            rho = T.AC(1);
            % H - Hurst coefficient (counting process)
            H = hurst_estimate(trace,'aggvar');
            % n - number of time scales to be modeled
            n = log2(order);
            % ds - number of IPP to be used
            ds = log2(order);
            % SA - set of autocorrelation coefficients of interarrival times
            SA = T.ACFull;
            % SAlags - lags used in SA
            SAlags = T.ACLags;
            line_printf('running Andersen-Nielsen fitting method\n');
            D = map_anfit(ls,rho,H,n,ds,SA,SAlags,iter_max,iter_tol);
            map = MAP(D);
        end

        function [map, logL] = emfit(trace, order, iter_max, iter_tol, verbose)
            % MAP = EMFIT(TRACE, ORDERS)

            % MAP = EMFIT(TRACE, ORDERS, ITER_MAX, ITER_TOL, VERBOSE)
            %
            % X = MAP.emfit(S, 4) % attempt all possible ER-CHMM structures in a MAP(4)
            % X = MAP.emfit(S, [1,3]) % attemp a ER-CHMM structure with an exponential and an Erlang-3

            if nargin< 3
                iter_max = 100;
            end

            if nargin< 4
                iter_tol = GlobalConstants.CoarseTol; % stop condition on the iterations
            end

            if nargin< 5
                verbose = true;
            end

            [D, logL] = erchmm_emfit(trace, order, iter_max, iter_tol, verbose);
            map = MAP(D);
        end
    end
end
