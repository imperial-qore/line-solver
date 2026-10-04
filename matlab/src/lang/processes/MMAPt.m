classdef MMAPt < ContinuousDistribution
    % MMAPt Time-inhomogeneous marked Markovian arrival process (MMAP_t).
    %
    % The two axes of MAPt and MarkedMAP crossed: arrivals are labelled with one
    % of K marks, AND the matrices that generate them are functions of the wall
    % clock. Segment j covers [breakpoints(j), breakpoints(j+1)) and carries
    % D0{j} together with the K blocks D1k{1}{j}, ..., D1k{K}{j}; D0 holds
    % transition rates without an arrival, D1k{c} the rates that generate an
    % arrival of mark c, and
    %
    %   D0{j} + sum_c D1k{c}{j}
    %
    % is a generator in every segment. The aggregate sum_c D1k{c}{j} is the D1
    % of the underlying MAPt, so hiding the marks recovers exactly that process.
    %
    % Two horizon conventions, as for MAPt and NHPP:
    %   cyclic     : the schedule repeats with period
    %                breakpoints(end) - breakpoints(1).
    %   non-cyclic : outside the horizon the process is frozen and emits
    %                nothing, so a non-cyclic MMAP_t is a transient construct.
    %
    % REDUCTIONS. With K = 1 this is exactly the MAPt with the same matrices,
    % and it is held to the same constructor rules so the reduction is exact
    % rather than merely close. With one segment, or with every segment
    % identical, it is exactly the stationary MMAP.
    %
    % Like MAPt this is neither renewal nor time-homogeneous, so getSCV,
    % getSkewness, evalCDF and evalLST return NaN rather than a value that would
    % misreport the process as stationary. It deliberately does NOT extend
    % Markovian: code gated on isMarkovian reads getProcess as a single
    % stationary pair and would silently drop the schedule.
    %
    % CONSTANT SUPPORT. The constructor requires one sparsity pattern across
    % segments, per mark block and for the off-diagonal of D0, reusing
    % MAPt.checkCommonSupport. The rule keeps the K = 1 reduction to MAPt exact
    % and leaves the fluid path open.
    %
    % At a Source the mark selects the class of the arriving job
    % (Source.setMarkedArrival, sn.markidx). As a SERVICE process it is sampled
    % for its duration and the mark is discarded, exactly as an MMAP is.
    %
    % The process representation stores:
    %   process{1} : 1-by-(n+1) row vector of breakpoints
    %   process{2} : 1-by-n cell of D0 matrices
    %   process{3} : 1-by-n cell of aggregate D1 matrices
    %   process{4} : logical, true if cyclic
    %   process{5} : 1-by-K cell of 1-by-n cells, the per-mark blocks
    %
    % The first four entries are exactly a MAPt slot, so sn_schedule_nominal and
    % every other time-blind consumer read the aggregate schedule unchanged; the
    % marks are an APPENDED fifth entry rather than a different shape.
    %
    % References: Q.-M. He, "The versatility of MMAP[K] and the MMAP[K]/G[K]/1
    % queue", Queueing Systems 38(4), 2001, for the marked structure;
    % Y. M. Ko, J. Pender, "Diffusion limits for the (MAP_t/Ph_t/inf)^N queueing
    % network", Oper. Res. Lett. 45(3), 2017, for the time-inhomogeneous one.
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    properties
        breakpoints; % 1-by-(n+1) segment boundaries, strictly increasing
        D0;          % 1-by-n cell of no-arrival rate matrices
        D1k;         % 1-by-K cell of 1-by-n cells, the per-mark blocks
        D1;          % 1-by-n cell of aggregate arrival matrices
        cyclic;      % logical, whether the schedule repeats
        process;     % {breakpoints, D0, D1, cyclic, D1k}, for serialization
        sampleClock; % wall-clock position of the next sample (see sample)
        samplePhase; % phase of the modulating chain at sampleClock
        lastMark;    % 1-based mark of the interval last returned by sample
    end

    methods
        function self = MMAPt(breakpoints, D0, D1k, cyclic)
            % SELF = MMAPT(BREAKPOINTS, D0, D1K, CYCLIC)
            %
            % D1K is a 1-by-K cell whose entry c is a 1-by-n cell of the blocks
            % of mark c, i.e. it is MARK-MAJOR.
            self@ContinuousDistribution('MMAPt', 0, [0, Inf]);
            if nargin < 4 || isempty(cyclic)
                cyclic = true;
            end
            if ~iscell(D0), D0 = {D0}; end
            if ~iscell(D1k) || isempty(D1k)
                line_error(mfilename, 'MMAPt: D1k must be a non-empty cell of per-mark segment cells');
            end
            if ~iscell(D1k{1}), D1k = {D1k}; end
            n = numel(D0);
            K = numel(D1k);
            if n == 0
                line_error(mfilename, 'MMAPt: D0 must be a non-empty cell of segment matrices');
            end
            for c = 1:K
                if ~iscell(D1k{c}) || numel(D1k{c}) ~= n
                    line_error(mfilename, sprintf(['MMAPt: mark block %d has %d segments against %d in D0; ' ...
                        'every mark must be defined on the whole schedule'], c, numel(D1k{c}), n));
                end
            end
            breakpoints = breakpoints(:).';
            if numel(breakpoints) ~= n + 1
                line_error(mfilename, 'MMAPt: breakpoints must have one more entry than the number of segments');
            end
            if any(diff(breakpoints) <= 0)
                line_error(mfilename, 'MMAPt: breakpoints must be strictly increasing');
            end
            h = size(D0{1}, 1);
            D1 = cell(1, n);
            for j = 1:n
                if ~isequal(size(D0{j}), [h h])
                    line_error(mfilename, sprintf('MMAPt: every D0 must be square of order %d; segment %d differs', h, j));
                end
                agg = zeros(h, h);
                for c = 1:K
                    B = D1k{c}{j};
                    if ~isequal(size(B), [h h])
                        line_error(mfilename, sprintf(['MMAPt: every mark block must be square of order %d; ' ...
                            'mark %d of segment %d differs'], h, c, j));
                    end
                    if any(B(:) < 0)
                        line_error(mfilename, sprintf('MMAPt: mark block %d must be non-negative in segment %d', c, j));
                    end
                    agg = agg + B;
                end
                off = D0{j} - diag(diag(D0{j}));
                if any(off(:) < 0)
                    line_error(mfilename, sprintf('MMAPt: off-diagonal D0 entries must be non-negative in segment %d', j));
                end
                if any(abs(sum(D0{j} + agg, 2)) > 1e-10)
                    line_error(mfilename, sprintf(['MMAPt: D0 plus the mark blocks must have zero row sums ' ...
                        '(generator) in segment %d'], j));
                end
                D1{j} = agg;
            end
            MAPt.checkCommonSupport(D0, true, 'off-diagonal D0');
            for c = 1:K
                MAPt.checkCommonSupport(D1k{c}, false, sprintf('D1 of mark %d', c));
            end
            if all(cellfun(@(M) sum(M(:)) <= 0, D1))
                line_error(mfilename, 'MMAPt: every segment has zero arrival intensity, so no event can ever occur');
            end
            self.breakpoints = breakpoints;
            self.D0 = D0(:).';
            self.D1 = D1;
            self.D1k = D1k(:).';
            self.cyclic = logical(cyclic);
            self.process = {self.breakpoints, self.D0, self.D1, self.cyclic, self.D1k};
            self.sampleClock = breakpoints(1);
            self.samplePhase = 1;
            self.lastMark = 0;
            self.mean = 1.0 / self.getTimeAverageRate();
            self.immediate = false;
        end

        function b = getBreakpoints(self)
            b = self.breakpoints;
        end

        function d = getD0Segments(self)
            d = self.D0;
        end

        function d = getD1Segments(self, k)
            % D = GETD1SEGMENTS() Aggregate per-segment D1.
            % D = GETD1SEGMENTS(K) Per-segment blocks of 1-based mark K.
            if nargin < 2
                d = self.D1;
            else
                if k < 1 || k > numel(self.D1k)
                    line_error(mfilename, sprintf('MMAPt: mark index out of range: %d', k));
                end
                d = self.D1k{k};
            end
        end

        function d = getMarkSegments(self)
            % D = GETMARKSEGMENTS() All mark blocks, mark-major.
            d = self.D1k;
        end

        function K = getNumberOfTypes(self)
            K = numel(self.D1k);
        end

        function c = isCyclic(self)
            c = self.cyclic;
        end

        function n = getNumSegments(self)
            n = numel(self.D0);
        end

        function h = getNumberOfPhases(self)
            h = size(self.D0{1}, 1);
        end

        function T = getPeriod(self)
            % T = GETPERIOD() Horizon length, which is the period when cyclic.
            T = self.breakpoints(end) - self.breakpoints(1);
        end

        function idx = getSegmentIndexAt(self, t)
            % IDX = GETSEGMENTINDEXAT(T) Active segment, 0 past a non-cyclic horizon.
            T = self.getPeriod();
            offset = t - self.breakpoints(1);
            if self.cyclic
                offset = mod(offset, T);
            elseif offset < 0 || offset >= T
                idx = 0;
                return
            end
            pos = self.breakpoints(1) + offset;
            idx = find(pos < self.breakpoints(2:end), 1, 'first');
            if isempty(idx)
                idx = numel(self.D0);
            end
        end

        function m = toMAPt(self)
            % M = TOMAPT()
            % The UNMARKED schedule, i.e. the MAPt whose D1 is the per-segment
            % aggregate. Hiding the marks is exact: an arrival of the MMAPt is
            % an arrival of this process regardless of its label.
            m = MAPt(self.breakpoints, self.D0, self.D1, self.cyclic);
        end

        function m = toMAPts(self, k)
            % M = TOMAPTS(K)
            % The MARGINAL schedule of 1-based mark K: arrivals fire only on
            % that mark's blocks, while the other marks' transitions become
            % hidden phase changes. Segment by segment this is
            % MAPt(D0 + D1 - D1k, D1k), the time-varying analogue of
            % MarkedMAP.toMAPs.
            if k < 1 || k > numel(self.D1k)
                line_error(mfilename, sprintf('MMAPt: mark index out of range: %d', k));
            end
            n = numel(self.D0);
            hidden = cell(1, n);
            for j = 1:n
                hidden{j} = self.D0{j} + self.D1{j} - self.D1k{k}{j};
            end
            m = MAPt(self.breakpoints, hidden, self.D1k{k}, self.cyclic);
        end

        function [D0bar, D1bar] = getTimeAverageProcess(self)
            % [D0BAR, D1BAR] = GETTIMEAVERAGEPROCESS()
            % Width-weighted average pair over the horizon. This is the
            % stationary carrier of the phase structure where a solver needs a
            % time-homogeneous one; a convex combination of generators is a
            % generator, so it is itself a valid MAP.
            [D0bar, D1bar] = self.toMAPt().getTimeAverageProcess();
        end

        function B = getTimeAverageMark(self, k)
            % B = GETTIMEAVERAGEMARK(K)
            % Width-weighted average of mark K's blocks, i.e. the D1 of that
            % mark in the nominal marked process.
            if k < 1 || k > numel(self.D1k)
                line_error(mfilename, sprintf('MMAPt: mark index out of range: %d', k));
            end
            widths = diff(self.breakpoints);
            total = sum(widths);
            B = zeros(size(self.D0{1}));
            for j = 1:numel(self.D0)
                B = B + (widths(j) / total) * self.D1k{k}{j};
            end
        end

        function rate = getTimeAverageRate(self)
            % RATE = GETTIMEAVERAGERATE() Arrival rate of the time-averaged aggregate MAP.
            [D0bar, D1bar] = self.getTimeAverageProcess();
            rate = map_lambda({D0bar, D1bar});
        end

        function lam = getTimeAverageMarkRates(self)
            % LAM = GETTIMEAVERAGEMARKRATES()
            % Per-mark arrival rates of the time-averaged process. These sum to
            % getTimeAverageRate, which is the identity a marked stream has to
            % satisfy: labelling the arrivals cannot change how many there are.
            % Through the canonical mmap_lambda rather than a second implementation
            % of theta*D1k*e: the rate of a marked counting process is weighted by
            % the STATIONARY phase distribution map_prob, not by the embedded
            % departure law map_pie, and the two differ whenever the phases do not
            % all arrive at the same rate.
            [D0bar, D1bar] = self.getTimeAverageProcess();
            K = numel(self.D1k);
            cell_mmap = cell(1, K + 2);
            cell_mmap{1} = D0bar;
            cell_mmap{2} = D1bar;
            for c = 1:K
                cell_mmap{2 + c} = self.getTimeAverageMark(c);
            end
            lam = reshape(mmap_lambda(cell_mmap), 1, K);
        end

        function mean = getMean(self)
            % MEAN = GETMEAN() Palm mean interval of the time-averaged aggregate MAP.
            mean = 1.0 / self.getTimeAverageRate();
        end

        function scv = getSCV(self)
            % SCV = GETSCV()
            % NaN: an MMAP_t is neither renewal nor time-homogeneous, so there
            % is no i.i.d. interval distribution for an SCV to summarise.
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
            proc = self.process;
        end

        function resetSampleClock(self)
            % RESETSAMPLECLOCK() Restart the sample path at the schedule start.
            self.sampleClock = self.breakpoints(1);
            self.samplePhase = 1;
            self.lastMark = 0;
        end

        function m = getLastMark(self)
            % M = GETLASTMARK() 1-based mark of the interval last sampled.
            m = self.lastMark;
        end

        function X = sample(self, n)
            % X = SAMPLE(N)
            % Draws N successive interarrival times along ONE sample path. Both
            % the intensity and the phase depend on absolute time, so this
            % advances an internal clock and phase across calls. Use
            % resetSampleClock to restart. A non-cyclic schedule that runs out
            % returns 0 for every remaining sample.
            if nargin < 2
                n = 1;
            end
            X = zeros(n, 1);
            for i = 1:n
                [interval, phase, mark] = self.nextArrival(self.sampleClock, self.samplePhase);
                X(i) = interval;
                if interval <= 0
                    break % horizon exhausted: no further arrival can occur
                end
                self.sampleClock = self.sampleClock + interval;
                self.samplePhase = phase;
                self.lastMark = mark;
            end
        end

        function [x, phase, mark] = nextArrival(self, from, phase)
            % [X, PHASE, MARK] = NEXTARRIVAL(FROM, PHASE)
            % Time to the next arrival from wall clock FROM in the given phase,
            % the phase after it and the 1-based mark it carries. Exact: within
            % a segment the phase process is a homogeneous CTMC, and by the
            % memoryless property the residual holding time may be redrawn at a
            % breakpoint, so the boundary is crossed by advancing the clock and
            % resampling under the new matrices.
            %
            % The mark costs NO EXTRA DRAW: the competing transitions are
            % accumulated destination-major and mark-minor, so the running total
            % after all marks of a destination equals the aggregate total after
            % that destination. The winning destination is therefore the one the
            % unmarked walk would choose for the same uniform, which is what
            % makes a K = 1 MMAPt reproduce MAPt sample path for sample path.
            elapsed = 0.0;
            pos = from;
            mark = 0;
            h = self.getNumberOfPhases();
            K = numel(self.D1k);
            while true
                idx = self.getSegmentIndexAt(pos);
                if idx == 0
                    x = 0.0;
                    return
                end
                T = self.getPeriod();
                offset = pos - self.breakpoints(1);
                if self.cyclic
                    offset = mod(offset, T);
                end
                toBoundary = (self.breakpoints(idx + 1) - self.breakpoints(1)) - offset;
                Dz = self.D0{idx};
                total = -Dz(phase, phase);
                if total <= 0
                    % Absorbing phase in this segment: only a boundary frees it.
                    if ~self.cyclic && idx == numel(self.D0)
                        x = 0.0;
                        return
                    end
                    elapsed = elapsed + toBoundary;
                    pos = pos + toBoundary;
                    continue
                end
                holding = -log(1 - rand()) / total;
                if holding >= toBoundary
                    if ~self.cyclic && idx == numel(self.D0)
                        x = 0.0;
                        return
                    end
                    elapsed = elapsed + toBoundary;
                    pos = pos + toBoundary;
                    continue
                end
                elapsed = elapsed + holding;
                pos = pos + holding;
                u = rand() * total;
                cum = 0.0;
                for j = 1:h
                    for c = 1:K
                        cum = cum + self.D1k{c}{idx}(phase, j);
                        if u < cum
                            x = elapsed;
                            phase = j;
                            mark = c;
                            return
                        end
                    end
                end
                moved = false;
                for j = 1:h
                    if j == phase
                        continue
                    end
                    cum = cum + Dz(phase, j);
                    if u < cum
                        phase = j;
                        moved = true;
                        break
                    end
                end
                if ~moved
                    % Rounding left u at or past the total: fall back to the last
                    % mark with any mass out of this phase, as MAPt does.
                    for j = h:-1:1
                        for c = K:-1:1
                            if self.D1k{c}{idx}(phase, j) > 0
                                x = elapsed;
                                phase = j;
                                mark = c;
                                return
                            end
                        end
                    end
                end
            end
        end
    end
end
