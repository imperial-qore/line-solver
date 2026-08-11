classdef DMAP < MarkovModulated
    % Discrete-time Markovian Arrival Process
    %
    % Models discrete-time arrival streams with D0 and D1 sub-stochastic matrices.
    % D0+D1 is a stochastic matrix (row sums = 1).
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    methods
        function self = DMAP(D0, D1)
            self@MarkovModulated('DMAP', 2);
            if nargin < 2 && iscell(D0)
                M = D0;
                D0 = M{1};
                D1 = M{2};
            end
            setParam(self, 1, 'D0', D0);
            setParam(self, 2, 'D1', D1);
            self.process = {D0, D1};
            if ~dmap_isfeasible(self.D)
                line_warning(mfilename, 'DMAP is infeasible.\n');
            end
        end

        function n = getNumberOfPhases(self)
            n = length(self.D(0));
        end

        function Di = D(self, i, wantSparse)
            if nargin < 3
                wantSparse = false;
            end
            if wantSparse
                if nargin < 2
                    Di = self.getProcess;
                else
                    Di = self.getProcess{i+1};
                end
            else
                if nargin < 2
                    Di = self.getProcess;
                    for j = 1:length(Di)
                        Di(j) = full(Di(j));
                    end
                else
                    Di = full(self.getProcess{i+1});
                end
            end
        end

        function MEAN = getMean(self)
            D0 = self.D(0);
            D1 = self.D(1);
            N = size(D0, 1);
            al = dmap_pie({D0, D1});
            MEAN = al * inv(eye(N) - D0) * ones(N, 1);
        end

        function meant = evalMeanT(self, t)
            meant = t / self.getMean();
        end

        function lam = getRate(self)
            lam = 1.0 / self.getMean();
        end

        function X = sample(self, n)
            if nargin < 2
                n = 1;
            end
            X = dmap_sample(self.getProcess, n);
        end

        function self = setMean(self, MEAN)
            D0 = self.D(0);
            D1 = self.D(1);
            N = size(D0, 1);
            I = eye(N);
            al = dmap_pie({D0, D1});
            currentMean = al * inv(I - D0) * ones(N, 1);
            scale = currentMean / MEAN;
            newImD0 = (I - D0) * scale;
            newD0 = I - newImD0;
            newD1 = D1 * scale;
            for i = 1:N
                rs = sum(newD0(i,:)) + sum(newD1(i,:));
                newD1(i,:) = newD1(i,:) / rs;
                newD0(i,:) = newD0(i,:) / rs;
            end
            self.params{1}.paramValue = newD0;
            self.params{2}.paramValue = newD1;
            self.process = {newD0, newD1};
        end

        function bool = isImmediate(self)
            bool = self.getMean < GlobalConstants.FineTol;
        end
    end

    methods (Static)
        function dmap = rand(order)
            if nargin < 1
                order = 2;
            end
            D0 = rand(order);
            for i = 1:order
                D0(i,:) = D0(i,:) / (sum(D0(i,:)) + 0.5 + rand);
            end
            remaining = ones(order,1) - sum(D0, 2);
            D1 = rand(order);
            for i = 1:order
                D1(i,:) = D1(i,:) / sum(D1(i,:)) * remaining(i);
            end
            dmap = DMAP(D0, D1);
        end
    end
end
