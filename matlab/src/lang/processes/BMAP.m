classdef BMAP < MarkedMAP
    % Batch Markovian Arrival Process (BMAP)
    %
    % BMAP is a point process where arrivals occur in batches.
    % Uses the standard BMAP representation:
    % - D0: infinitesimal generator for transitions without arrivals
    % - D1: rate matrix for transitions generating 1 arrival
    % - D2: rate matrix for transitions generating 2 arrivals
    % - ...
    % - Dk: rate matrix for transitions generating k arrivals
    %
    % BMAP extends MarkedMAP where each "mark" k represents a batch size k.
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    methods
        % Constructor
        function self = BMAP(D)
            % SELF = BMAP(D)
            %
            % D is a cell array {D0, D1, D2, ..., Dk} where:
            % - D0: transitions without arrivals
            % - Dk: transitions generating k arrivals (for k >= 1)
            %
            % The number of marking types K equals the maximum batch size
            % (i.e., K = length(D) - 2 when D also includes D1_total,
            %  or K = length(D) - 1 otherwise)

            if nargin == 0
                D = {};
            end

            % Determine K (number of marking types = max batch size)
            % If D = {D0, D1, D2, ..., Dk}, then K = k (max batch size)
            % MarkedMAP expects either:
            %   - {D0, D1_total, D11, D12, ..., D1K} with K marking types
            %   - {D0, D11, D12, ..., D1K} with K marking types

            if isempty(D)
                K = 0;
            else
                % Assume format {D0, D1, D2, ..., Dk} where Dk generates k arrivals
                % Need to convert to MarkedMAP format {D0, D1_total, D1, D2, ..., Dk}
                % where D1_total = sum(D1, D2, ..., Dk)
                K = length(D) - 1;  % max batch size

                % Compute total arrival matrix
                D1_total = D{2};
                for k = 3:length(D)
                    D1_total = D1_total + D{k};
                end

                % Create MarkedMAP cell array
                D_mmap = cell(1, K+2);
                D_mmap{1} = D{1};  % D0
                D_mmap{2} = D1_total;  % D1_total
                for k = 1:K
                    D_mmap{2+k} = D{1+k};  % Dk (batch size k)
                end

                D = D_mmap;
            end

            % Call parent MarkedMAP constructor
            self@MarkedMAP(D, K);

            % Fix self.process which is incorrectly set by MarkedMAP constructor
            % MarkedMAP constructor has a bug on line 33-37 where it creates
            % a cell array of size K-1 instead of K+2
            self.process = D;

            % Validate generator
            self.validateGenerator();
        end

        % Validate that D0 + sum(Dk) forms a proper infinitesimal generator
        function validateGenerator(self)
            D0 = self.D(0);
            nPhases = length(D0);

            % Compute D0 + sum(Dk)
            generator = D0;
            for k = 1:self.getMaxBatchSize()
                generator = generator + self.D(1, k);  % D(1, k) returns process{2+k} = Dk
            end

            % Check row sums (should be close to zero)
            rowSums = sum(generator, 2);
            maxRowSum = max(abs(rowSums));

            if maxRowSum > 1e-6
                warning('BMAP:InvalidGenerator', ...
                    'BMAP generator row sums are not zero (max: %.2e). Consider calling normalize().', maxRowSum);
            end
        end

        % Get maximum batch size
        function k = getMaxBatchSize(self)
            % K = GETMAXBATCHSIZE()
            %
            % Returns the maximum batch size k where Dk is defined
            k = self.getNumberOfTypes();
        end

        % Get mean batch size
        function mean_bs = getMeanBatchSize(self)
            % MEAN_BS = GETMEANBATCHSIZE()
            %
            % Computes the mean batch size as a weighted average:
            % E[batch size] = sum(k * rate_k) / sum(rate_k)

            % Get stationary distribution using D0 and D1_total
            D0 = self.D(0);
            D1_total = self.D(1);
            pi = map_pie({D0, D1_total});

            % Compute weighted average
            totalRate = 0;
            weightedSum = 0;
            maxBatch = self.getMaxBatchSize();

            for k = 1:maxBatch
                Dk = self.D(1, k);  % D(1, k) returns process{1+1+k} = process{2+k} = Dk
                rate_k = sum(pi * Dk * ones(size(Dk, 1), 1));
                totalRate = totalRate + rate_k;
                weightedSum = weightedSum + k * rate_k;
            end

            if totalRate > 0
                mean_bs = weightedSum / totalRate;
            else
                mean_bs = 0;
            end
        end

        % Get batch rates
        function rates = getBatchRates(self)
            % RATES = GETBATCHRATES()
            %
            % Returns array where rates(k) is the rate of batch size k arrivals

            D0 = self.D(0);
            D1_total = self.D(1);
            maxBatch = self.getMaxBatchSize();

            % Get stationary distribution
            pi = map_pie({D0, D1_total});

            % Compute rate for each batch size
            rates = zeros(1, maxBatch);
            for k = 1:maxBatch
                Dk = self.D(1, k);  % D(1, k) returns process{1+1+k} = process{2+k} = Dk
                rates(k) = sum(pi * Dk * ones(size(Dk, 1), 1));
            end
        end

        % Override toString-like display
        function display(self)
            % DISPLAY(SELF)
            fprintf('BMAP(phases=%d, maxBatchSize=%d)\n', ...
                self.getNumberOfPhases(), self.getMaxBatchSize());
        end

        % Get distribution name
        function name = getName(self)
            % NAME = GETNAME()
            %
            % Returns the name of this distribution type
            name = 'BMAP';
        end

        % Sample from BMAP (returns inter-arrival times and batch sizes)
        function [X, B] = sample(self, n)
            % [X, B] = SAMPLE(N)
            %
            % Sample n batches from the BMAP
            % Returns:
            %   X: inter-arrival times (between batches)
            %   B: batch sizes for each arrival

            if nargin < 2
                n = 1;
            end

            % Use parent MarkedMAP sampling
            % This returns inter-arrival times and class labels
            [X, C] = sample@MarkedMAP(self, n);

            % Convert class labels to batch sizes
            % Class k corresponds to batch size k
            B = C;
        end

        % Sample only inter-batch times
        function X = sampleInterBatch(self, n)
            % X = SAMPLEINTERBATCH(N)
            %
            % Sample n inter-batch arrival times
            % (ignoring batch sizes)

            if nargin < 2
                n = 1;
            end

            % Sample from aggregate MAP
            map = self.toMAP();
            X = map.sample(n);
        end

        % Get inter-batch MAP
        function map = getInterBatchMAP(self)
            % MAP = GETINTERBATCHMAP()
            %
            % Returns the underlying MAP for inter-batch arrivals
            map = self.toMAP();
        end
    end

    methods (Static)
        % Factory method: create BMAP from MAP + batch size distribution
        function bmap = fromMAPWithBatchPMF(D0, D1, batchSizes, pmf)
            % BMAP = FROMMAPWITHBATCHPMF(D0, D1, BATCHSIZES, PMF)
            %
            % Create BMAP from a base MAP and batch size distribution
            %
            % Inputs:
            %   D0: base MAP's D0 matrix (transitions without batch arrivals)
            %   D1: base MAP's D1 matrix (inter-batch arrival transitions)
            %   batchSizes: array of batch sizes (e.g., [1, 2, 4, 8])
            %   pmf: probability mass function for batch sizes (must sum to 1)
            %
            % Output:
            %   bmap: BMAP constructed from the base MAP and batch distribution

            if length(batchSizes) ~= length(pmf)
                line_error(mfilename, 'Batch sizes and PMF must have the same length');
            end

            % Normalize PMF
            pmf = pmf / sum(pmf);

            % Find maximum batch size
            maxBatch = max(batchSizes);

            % Create D matrices: Dk = D1 * pmf(batchSize==k)
            D = cell(1, maxBatch + 1);
            D{1} = D0;  % D0

            % Initialize Dk matrices to zero
            for k = 1:maxBatch
                D{1+k} = zeros(size(D0));
            end

            % Set Dk based on batch sizes and PMF
            for i = 1:length(batchSizes)
                k = batchSizes(i);
                if k < 1 || k > maxBatch
                    line_error(mfilename, sprintf('Invalid batch size: %d', k));
                end
                D{1+k} = D{1+k} + D1 * pmf(i);
            end

            % Create BMAP
            bmap = BMAP(D);
        end

        % Generate random BMAP
        function bmap = rand(order, maxBatchSize)
            % BMAP = RAND(ORDER, MAXBATCHSIZE)
            %
            % Generate random BMAP using uniform random numbers
            %
            % Inputs:
            %   order: number of phases (default: 2)
            %   maxBatchSize: maximum batch size (default: 3)

            if nargin < 1
                order = 2;
            end
            if nargin < 2
                maxBatchSize = 3;
            end

            % Generate random MarkedMAP and interpret as BMAP
            mmap = MarkedMAP.rand(order, maxBatchSize);

            % Convert to BMAP format
            D = cell(1, maxBatchSize + 1);
            D{1} = mmap.D(0);  % D0
            for k = 1:maxBatchSize
                D{1+k} = mmap.D(1, k-1);  % Dk
            end

            bmap = BMAP(D);
        end
    end
end
