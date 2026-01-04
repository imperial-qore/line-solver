classdef MAP < MarkovModulated
    % Markovian Arrival Process for correlated arrival modeling
    %
    % Models arrival streams with correlation and burstiness via D0 and D1 matrices.
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    methods
        %Constructor
        function self = MAP(D0,D1)
            % MAP Create a Markovian Arrival Process instance
            %
            % @brief Creates a MAP with specified D0 and D1 matrices
            % @param D0 Generator matrix for transitions without arrivals
            % @param D1 Rate matrix for transitions with arrivals  
            % @return self MAP instance with specified matrices

            self@MarkovModulated('MAP',2);
            if nargin < 2 && iscell(D0)
                M = D0;
                D0 = M{1};
                D1 = M{2};
            end
            setParam(self, 1, 'D0', D0);
            setParam(self, 2, 'D1', D1);
            self.process = {D0,D1};
            if ~map_isfeasible(self.D)
                line_warning(mfilename,'MAP is infeasible.\n');
            end
        end

        function n = getNumberOfPhases(self)
            n = length(D(0));
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

        function lam = getRate(self)
            % MAP = GETRATE()

            lam = map_lambda(self.D);
        end        

        function mapr = toTimeReversed(self)
            mapr = MAP(map_timereverse(self.D));
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

        function self = setMean(self,MEAN)
            % UPDATEMEAN(SELF,MEAN)
            % Update parameters to match the given mean
            newMAP = map_scale(self.D,MEAN);
            self.params{1}.paramValue = newMAP{1};
            self.params{2}.paramValue = newMAP{2};
        end

        function bool = isImmediate(self)
            bool = self.getMean < GlobalConstants.FineTol;
        end

        function mmdp = toMMDP(self)
            % TOMMDP Convert MAP to MMDP (deterministic representation)
            %
            % Converts this Markovian Arrival Process to a Markov-Modulated
            % Deterministic Process suitable for fluid queue analysis.
            %
            % @return mmdp MMDP representation of this MAP

            mmdp = MMDP.fromMAP(self);
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

    end
end
