classdef MPH < MarkedMAP
    % MPH Marked phase-type distribution.
    %
    % A PH whose absorption is LABELLED: the chain runs on h transient phases
    % with entry law alpha and sub-generator S, and absorption happens through
    % one of K exit vectors s_1, ..., s_K, so a completion carries both a
    % duration and the mark of the exit that produced it. The partition identity
    %
    %   sum_k s_k = -S*e
    %
    % is what makes the K exits account for exactly the absorption the
    % sub-generator leaves, no more and no less.
    %
    % AN MPH IS THE RENEWAL SPECIAL CASE OF AN MMAP. Lowering it by
    %
    %   D0 = S,   D1k = s_k*alpha,   D1 = sum_k D1k = (-S*e)*alpha
    %
    % gives a marked Markovian arrival process whose successive intervals are
    % independent and PH(alpha,S) distributed, because D1k factorises through
    % alpha and so the phase after an event does not depend on the phase before
    % it. That is exactly the M3A cell {D0, D1, D11, ..., D1K} this class hands
    % to MarkedMAP, so every consumer of that layout serves an MPH unchanged.
    %
    % It nevertheless carries a ProcessType.MPH of its own rather than aliasing
    % MMAP, so a solver that cannot honour a marked renewal process refuses it
    % by name instead of inheriting MMAP's support.
    %
    % At a Source the mark selects the class of the arriving job
    % (Source.setMarkedArrival, sn.markidx). As a SERVICE process it is sampled
    % for its duration and the mark is discarded, exactly as an MMAP is.
    %
    % Reference: Q.-M. He, "The versatility of MMAP[K] and the MMAP[K]/G[K]/1
    % queue", Queueing Systems 38(4), 2001.
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    properties
        alpha;   % 1-by-h entry law
        subgen;  % h-by-h sub-generator S
        exit;    % 1-by-K cell of h-by-1 exit vectors
    end

    methods
        function self = MPH(alpha, S, exits)
            % SELF = MPH(ALPHA, S, EXITS)
            %
            % ALPHA is a 1-by-h row, S an h-by-h sub-generator and EXITS a cell
            % of K non-negative h-by-1 vectors summing to -S*e.
            tol = 1e-10;
            if nargin < 3 || isempty(alpha) || isempty(S) || isempty(exits)
                line_error(mfilename, 'MPH: alpha, S and a non-empty cell of exit vectors are all required');
            end
            if ~iscell(exits), exits = {exits}; end
            alpha = alpha(:).';
            h = numel(alpha);
            if ~isequal(size(S), [h h])
                line_error(mfilename, sprintf('MPH: S must be square of order %d to match alpha', h));
            end
            if any(alpha < 0)
                line_error(mfilename, 'MPH: alpha must be non-negative');
            end
            if abs(sum(alpha) - 1) > tol
                line_error(mfilename, sprintf(['MPH: alpha must sum to one; it sums to %g. A defective ' ...
                    'entry law would put mass on an instantaneous completion that carries no mark.'], sum(alpha)));
            end
            if any(diag(S) >= 0)
                line_error(mfilename, 'MPH: the diagonal of S must be negative');
            end
            off = S - diag(diag(S));
            if any(off(:) < 0)
                line_error(mfilename, 'MPH: off-diagonal entries of S must be non-negative');
            end
            K = numel(exits);
            marked = zeros(h, 1);
            for k = 1:K
                sk = exits{k}(:);
                if numel(sk) ~= h
                    line_error(mfilename, sprintf('MPH: exit vector %d must have %d entries', k, h));
                end
                if any(sk < 0)
                    line_error(mfilename, sprintf('MPH: exit vector %d must be non-negative', k));
                end
                exits{k} = sk;
                marked = marked + sk;
            end
            absorb = -sum(S, 2);
            if any(abs(marked - absorb) > tol)
                line_error(mfilename, ['MPH: the exit vectors must partition the absorption rate of S. ' ...
                    'Every way of leaving the transient phases must carry exactly one mark.']);
            end
            if all(marked == 0)
                line_error(mfilename, 'MPH: every exit vector is zero, so absorption can never occur');
            end

            % The lowering, in the M3A layout {D0, D1, D11, ..., D1K}.
            D = cell(1, K + 2);
            D{1} = S;
            agg = zeros(h, h);
            for k = 1:K
                D{2 + k} = exits{k} * alpha;
                agg = agg + D{2 + k};
            end
            D{2} = agg;
            self@MarkedMAP(D, K);
            self.name = 'MPH';
            self.alpha = alpha;
            self.subgen = S;
            self.exit = exits;
        end

        function a = getAlpha(self)
            % A = GETALPHA() The 1-by-h entry law.
            a = self.alpha;
        end

        function S = getSubgenerator(self)
            % S = GETSUBGENERATOR() The h-by-h sub-generator.
            S = self.subgen;
        end

        function s = getExitVector(self, k)
            % S = GETEXITVECTOR(K) The h-by-1 exit vector of 1-based mark K.
            if k < 1 || k > numel(self.exit)
                line_error(mfilename, sprintf('MPH: mark index out of range: %d', k));
            end
            s = self.exit{k};
        end

        function e = getExitVectors(self)
            % E = GETEXITVECTORS() The exit vectors, ordered by mark.
            e = self.exit;
        end

        function h = getNumberOfPhases(self)
            % H = GETNUMBEROFPHASES() Number of transient phases.
            h = numel(self.alpha);
        end

        function ph = toPH(self)
            % PH = TOPH()
            % The UNMARKED phase-type law of the duration, PH(alpha, S). Hiding
            % the marks leaves the interval law untouched, which is what
            % separates an MPH from a general MMAP: the aggregate of an MMAP is
            % a MAP, and only here is it a renewal process with a phase-type
            % marginal.
            ph = PH(self.alpha, self.subgen);
        end

        function p = getMarkProbabilities(self)
            % P = GETMARKPROBABILITIES()
            % The probability that a completion carries each mark,
            % alpha*(-S)^-1*s_k. Successive completions are independent, so this
            % is both the long-run fraction of mark-k completions and the
            % probability of any one of them.
            tau = self.alpha / (-self.subgen);
            p = zeros(1, numel(self.exit));
            for k = 1:numel(self.exit)
                p(k) = tau * self.exit{k};
            end
        end
    end
end
