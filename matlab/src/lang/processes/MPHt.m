classdef MPHt < ContinuousDistribution
    % MPHt Time-inhomogeneous marked phase-type distribution (MPH_t).
    %
    % An MPH whose three ingredients are functions of the wall clock: segment j
    % covers [breakpoints(j), breakpoints(j+1)) and carries an entry law
    % alpha{j}, a sub-generator S{j} and K exit vectors exits{1}{j}, ...,
    % exits{K}{j} satisfying the partition identity
    %
    %   sum_c exits{c}{j} = -S{j}*e
    %
    % in every segment. This is to PHt what MPH is to PH.
    %
    % IT IS STORED LOWERED TO MMAPt FORM, segment by segment, by
    %
    %   D0{j} = S{j},   D1k{c}{j} = exits{c}{j} * alpha{j}
    %
    % exactly as the C++ port stores a PHt as its equivalent MAPt pair. One
    % marked schedule shape therefore reaches sn.proc for both families, and
    % ProcessType.isMarkedSchedule covers both. The original alpha, S and exit
    % vectors are kept on the object for the getters and for serialization, so
    % a round trip returns an MPHt and not the MMAPt it lowers to.
    %
    % The lowering makes each segment's arrivals RENEWAL within that segment:
    % D1k{c}{j} factorises through alpha{j}, so the phase after a completion
    % does not depend on the phase before it. Across a breakpoint the process is
    % still time-varying, which is what separates an MPH_t from an MPH.
    %
    % At a Source the mark selects the class of the arriving job
    % (Source.setMarkedArrival, sn.markidx). As a SERVICE process it is sampled
    % for its duration and the mark is discarded, exactly as an MMAP is.
    %
    % References: Q.-M. He, "The versatility of MMAP[K] and the MMAP[K]/G[K]/1
    % queue", Queueing Systems 38(4), 2001; Y. M. Ko, J. Pender, "Diffusion
    % limits for the (MAP_t/Ph_t/inf)^N queueing network", Oper. Res. Lett.
    % 45(3), 2017.
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    properties
        breakpoints; % 1-by-(n+1) segment boundaries, strictly increasing
        alpha;       % 1-by-n cell of 1-by-h entry laws
        subgen;      % 1-by-n cell of h-by-h sub-generators
        exit;        % 1-by-K cell of 1-by-n cells of h-by-1 exit vectors
        cyclic;      % logical, whether the schedule repeats
        lowered;     % the equivalent MMAPt, which is what is sampled and stored
    end

    methods
        function self = MPHt(breakpoints, alpha, S, exits, cyclic)
            % SELF = MPHT(BREAKPOINTS, ALPHA, S, EXITS, CYCLIC)
            %
            % ALPHA and S are 1-by-n cells; EXITS is a 1-by-K cell whose entry c
            % is a 1-by-n cell of h-by-1 vectors, i.e. it is MARK-MAJOR.
            self@ContinuousDistribution('MPHt', 0, [0, Inf]);
            tol = 1e-10;
            if nargin < 5 || isempty(cyclic)
                cyclic = true;
            end
            if ~iscell(alpha), alpha = {alpha}; end
            if ~iscell(S), S = {S}; end
            if ~iscell(exits) || isempty(exits)
                line_error(mfilename, 'MPHt: exits must be a non-empty cell of per-mark segment cells');
            end
            if ~iscell(exits{1}), exits = {exits}; end
            n = numel(alpha);
            K = numel(exits);
            if n == 0 || numel(S) ~= n
                line_error(mfilename, 'MPHt: alpha and S must be non-empty cells of equal length');
            end
            for c = 1:K
                if ~iscell(exits{c}) || numel(exits{c}) ~= n
                    line_error(mfilename, sprintf(['MPHt: mark exit %d has %d segments against %d in alpha; ' ...
                        'every mark must be defined on the whole schedule'], c, numel(exits{c}), n));
                end
            end
            h = numel(alpha{1});
            D0 = cell(1, n);
            D1k = cell(1, K);
            for c = 1:K
                D1k{c} = cell(1, n);
            end
            for j = 1:n
                a = alpha{j}(:).';
                Sm = S{j};
                if numel(a) ~= h
                    line_error(mfilename, sprintf('MPHt: every alpha must have %d entries; segment %d differs', h, j));
                end
                if ~isequal(size(Sm), [h h])
                    line_error(mfilename, sprintf('MPHt: every S must be square of order %d; segment %d differs', h, j));
                end
                if any(a < 0)
                    line_error(mfilename, sprintf('MPHt: alpha must be non-negative in segment %d', j));
                end
                if abs(sum(a) - 1) > tol
                    line_error(mfilename, sprintf('MPHt: alpha must sum to one in every segment; segment %d sums to %g', j, sum(a)));
                end
                if any(diag(Sm) >= 0)
                    line_error(mfilename, sprintf('MPHt: the diagonal of S must be negative in segment %d', j));
                end
                off = Sm - diag(diag(Sm));
                if any(off(:) < 0)
                    line_error(mfilename, sprintf('MPHt: off-diagonal entries of S must be non-negative in segment %d', j));
                end
                marked = zeros(h, 1);
                for c = 1:K
                    sk = exits{c}{j}(:);
                    if numel(sk) ~= h
                        line_error(mfilename, sprintf('MPHt: exit vector of mark %d in segment %d must have %d entries', c, j, h));
                    end
                    if any(sk < 0)
                        line_error(mfilename, sprintf('MPHt: exit vectors must be non-negative; mark %d of segment %d is not', c, j));
                    end
                    exits{c}{j} = sk;
                    marked = marked + sk;
                    D1k{c}{j} = sk * a;
                end
                if any(abs(marked + sum(Sm, 2)) > tol)
                    line_error(mfilename, sprintf(['MPHt: the exit vectors must partition the absorption rate of S; ' ...
                        'segment %d does not'], j));
                end
                alpha{j} = a;
                D0{j} = Sm;
            end

            % The lowering. Its constructor re-checks the generator property and
            % the common support, so a schedule that passes here is a valid
            % MMAPt by construction.
            self.lowered = MMAPt(breakpoints, D0, D1k, cyclic);
            self.breakpoints = breakpoints(:).';
            self.alpha = alpha(:).';
            self.subgen = S(:).';
            self.exit = exits(:).';
            self.cyclic = logical(cyclic);
            self.mean = self.lowered.getMean();
            self.immediate = false;
        end

        function b = getBreakpoints(self)
            b = self.breakpoints;
        end

        function a = getAlphaSegments(self)
            a = self.alpha;
        end

        function s = getSSegments(self)
            s = self.subgen;
        end

        function e = getExitSegments(self, k)
            % E = GETEXITSEGMENTS() All exit vectors, mark-major.
            % E = GETEXITSEGMENTS(K) Per-segment exit vectors of 1-based mark K.
            if nargin < 2
                e = self.exit;
            else
                if k < 1 || k > numel(self.exit)
                    line_error(mfilename, sprintf('MPHt: mark index out of range: %d', k));
                end
                e = self.exit{k};
            end
        end

        function K = getNumberOfTypes(self)
            K = numel(self.exit);
        end

        function c = isCyclic(self)
            c = self.cyclic;
        end

        function n = getNumSegments(self)
            n = numel(self.alpha);
        end

        function h = getNumberOfPhases(self)
            h = numel(self.alpha{1});
        end

        function T = getPeriod(self)
            T = self.breakpoints(end) - self.breakpoints(1);
        end

        function idx = getSegmentIndexAt(self, t)
            idx = self.lowered.getSegmentIndexAt(t);
        end

        function m = toMMAPt(self)
            % M = TOMMAPT() The equivalent MMAPt, which is the form this process
            % is stored and sampled in.
            m = self.lowered;
        end

        function p = toPHt(self)
            % P = TOPHT() The UNMARKED schedule of the duration, PHt(alpha, S).
            p = PHt(self.breakpoints, self.alpha, self.subgen, self.cyclic);
        end

        function m = toMAPts(self, k)
            % M = TOMAPTS(K) The marginal schedule of mark K; see MMAPt.toMAPts.
            m = self.lowered.toMAPts(k);
        end

        function [D0bar, D1bar] = getTimeAverageProcess(self)
            % [D0BAR, D1BAR] = GETTIMEAVERAGEPROCESS()
            % Width-weighted average pair of the lowered process.
            [D0bar, D1bar] = self.lowered.getTimeAverageProcess();
        end

        function rate = getTimeAverageRate(self)
            rate = self.lowered.getTimeAverageRate();
        end

        function lam = getTimeAverageMarkRates(self)
            % LAM = GETTIMEAVERAGEMARKRATES() See MMAPt.
            lam = self.lowered.getTimeAverageMarkRates();
        end

        function mean = getMean(self)
            mean = self.lowered.getMean();
        end

        function scv = getSCV(self)
            % SCV = GETSCV() NaN: an MPH_t is time-varying, so no i.i.d.
            % interval law summarises it.
            scv = NaN;
        end

        function skew = getSkewness(self)
            % SKEW = GETSKEWNESS() NaN; see getSCV.
            skew = NaN;
        end

        function F = evalCDF(self, t)
            % F = EVALCDF(T) NaN; see getSCV.
            F = NaN(size(t));
        end

        function L = evalLST(self, s)
            % L = EVALLST(S) NaN; see getSCV.
            L = NaN(size(s));
        end

        function proc = getProcess(self)
            % PROC = GETPROCESS()
            % The MMAPt-shaped cell of the lowering, which is what reaches
            % sn.proc. The procid stays MPHT, so the featset gate and the wire
            % type still tell the two families apart; only the representation is
            % shared.
            proc = self.lowered.getProcess();
        end

        function resetSampleClock(self)
            % RESETSAMPLECLOCK() Restart the sample path at the schedule start.
            self.lowered.resetSampleClock();
        end

        function m = getLastMark(self)
            % M = GETLASTMARK() 1-based mark of the interval last sampled.
            m = self.lowered.getLastMark();
        end

        function X = sample(self, n)
            % X = SAMPLE(N) See MMAPt.sample.
            if nargin < 2
                n = 1;
            end
            X = self.lowered.sample(n);
        end

        function [x, phase, mark] = nextArrival(self, from, phase)
            % [X, PHASE, MARK] = NEXTARRIVAL(FROM, PHASE) See MMAPt.nextArrival.
            [x, phase, mark] = self.lowered.nextArrival(from, phase);
        end
    end
end
