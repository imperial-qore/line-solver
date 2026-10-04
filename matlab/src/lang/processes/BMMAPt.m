classdef BMMAPt < ContinuousDistribution
    % BMMAPt Batch marked time-inhomogeneous Markovian arrival process.
    %
    % The three axes LINE models for an arrival stream, crossed: arrivals are
    % labelled with one of K marks, they occur in BATCHES of up to B jobs, AND
    % the matrices that generate them are functions of the wall clock. Segment j
    % covers [breakpoints(j), breakpoints(j+1)) and carries D0{j} together with
    % the K*B blocks D1kb{c}{b}{j}, where D1kb{c}{b}{j} holds the rates that
    % release a batch of b jobs all of mark c, and
    %
    %   D0{j} + sum_c sum_b D1kb{c}{b}{j}
    %
    % is a generator in every segment. A BATCH IS HOMOGENEOUS IN ITS MARK: one
    % epoch releases b jobs that all carry mark c. This is the BMMAP[K] of the
    % queueing literature and it is what the simulation engines can release,
    % since a batch is dispatched under a single class.
    %
    % Two derived levels are stored beside the blocks and they are what keeps
    % every existing consumer working:
    %
    %   D1k{c}{j} = sum_b D1kb{c}{b}{j}   the per-mark schedule, batches hidden
    %   D1{j}     = sum_c D1k{c}{j}       the aggregate schedule, a MAPt
    %
    % So a consumer that ignores batches reads exactly the MMAPt of D1k, and one
    % that ignores marks as well reads exactly the MAPt of D1.
    %
    % Two horizon conventions, as for MAPt, MMAPt and NHPP:
    %   cyclic     : the schedule repeats with period
    %                breakpoints(end) - breakpoints(1).
    %   non-cyclic : outside the horizon the process is frozen and emits
    %                nothing, so a non-cyclic BMMAP_t is a transient construct.
    %
    % REDUCTIONS, and they are exact rather than merely close because the
    % constructor holds this class to the same rules as the families it reduces
    % to. With every batch size 1 (B = 1) this is the MMAPt with the same
    % blocks; with K = 1 it is the unmarked batch schedule; with both it is the
    % MAPt with the same matrices; with every segment identical it is the
    % stationary BMAP.
    %
    % Like MAPt and MMAPt this is neither renewal nor time-homogeneous, so
    % getSCV, getSkewness, evalCDF and evalLST return NaN rather than a value
    % that would misreport the process as stationary. It deliberately does NOT
    % extend Markovian, and it deliberately does NOT subclass MMAPt: a
    % consumer that tested isa(dist,'MMAPt') would then accept it and silently
    % DROP the batches, which is a wrong answer rather than a refusal.
    %
    % RATES. getTimeAverageRate is the EVENT (batch epoch) rate and getMean its
    % reciprocal, the Palm mean inter-batch interval, matching BMAP, whose
    % getInterBatchMAP is the aggregate pair. The JOB rate, which is what a
    % station's throughput must balance against, is getTimeAverageJobRate and
    % equals sum_c sum_b b*rate(c,b).
    %
    % CONSTANT SUPPORT. The constructor requires one sparsity pattern across
    % segments for the off-diagonal of D0, for each BATCH-AGGREGATED mark block
    % D1k{c}, and for the aggregate D1, reusing MAPt.checkCommonSupport. Those
    % three are precisely what toMAPt and toMMAPt pass to constructors that
    % enforce the rule themselves, so the reductions stay constructible.
    %
    % THE INDIVIDUAL (mark, batch) BLOCKS ARE DELIBERATELY EXEMPT. Requiring one
    % pattern there too would forbid the composition changing with the segment
    % -- pairs in the morning, singles at night -- which is the one thing this
    % family exists to express and which neither reduction needs. The batch axis
    % is also DENSE in b = 1..B, as BMAP's {D0, D1, ..., Dk} is, so an unused
    % batch size is declared as a zero block rather than omitted.
    %
    % The process representation stores:
    %   process{1} : 1-by-(n+1) row vector of breakpoints
    %   process{2} : 1-by-n cell of D0 matrices
    %   process{3} : 1-by-n cell of aggregate D1 matrices
    %   process{4} : logical, true if cyclic
    %   process{5} : 1-by-K cell of 1-by-n cells, the per-mark blocks
    %   process{6} : 1-by-K cell of 1-by-B cells of 1-by-n cells, the blocks
    %
    % The first five entries are exactly an MMAPt slot, and the first four
    % exactly a MAPt slot, so sn_schedule_nominal and every other time-blind or
    % batch-blind consumer reads the schedule it always did; the batch blocks
    % are an APPENDED sixth entry rather than a different shape.
    %
    % References: Q.-M. He, "The versatility of MMAP[K] and the MMAP[K]/G[K]/1
    % queue", Queueing Systems 38(4), 2001, for the marked structure;
    % D. M. Lucantoni, "New results on the single server queue with a batch
    % Markovian arrival process", Stochastic Models 7(1), 1991, for the batch
    % structure; Y. M. Ko, J. Pender, "Diffusion limits for the
    % (MAP_t/Ph_t/inf)^N queueing network", Oper. Res. Lett. 45(3), 2017, for
    % the time-inhomogeneous one.
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    properties
        breakpoints; % 1-by-(n+1) segment boundaries, strictly increasing
        D0;          % 1-by-n cell of no-arrival rate matrices
        D1kb;        % 1-by-K cell of 1-by-B cells of 1-by-n cells, the blocks
        D1k;         % 1-by-K cell of 1-by-n cells, batch-aggregated per mark
        D1;          % 1-by-n cell of aggregate arrival matrices
        cyclic;      % logical, whether the schedule repeats
        process;     % {breakpoints, D0, D1, cyclic, D1k, D1kb}, for serialization
        sampleClock; % wall-clock position of the next sample (see sample)
        samplePhase; % phase of the modulating chain at sampleClock
        lastMark;    % 1-based mark of the interval last returned by sample
        lastBatch;   % batch size of the interval last returned by sample
    end

    methods
        function self = BMMAPt(breakpoints, D0, D1kb, cyclic)
            % SELF = BMMAPT(BREAKPOINTS, D0, D1KB, CYCLIC)
            %
            % D1KB is a 1-by-K cell whose entry c is a 1-by-B cell whose entry b
            % is a 1-by-n cell of the segment blocks releasing a batch of b jobs
            % of mark c. It is MARK-MAJOR, then BATCH, then segment.
            self@ContinuousDistribution('BMMAPt', 0, [0, Inf]);
            if nargin < 4 || isempty(cyclic)
                cyclic = true;
            end
            if ~iscell(D0), D0 = {D0}; end
            if ~iscell(D1kb) || isempty(D1kb)
                line_error(mfilename, 'BMMAPt: D1kb must be a non-empty cell of per-mark batch cells');
            end
            % A single mark may be given as its 1-by-B cell of segment cells,
            % and a single batch size as its 1-by-n cell of segment matrices.
            if ~iscell(D1kb{1}), D1kb = {D1kb}; end
            for c = 1:numel(D1kb)
                if ~iscell(D1kb{c})
                    line_error(mfilename, sprintf('BMMAPt: mark %d must carry a cell of batch-size blocks', c));
                end
                if ~iscell(D1kb{c}{1}), D1kb{c} = {D1kb{c}}; end
            end
            n = numel(D0);
            K = numel(D1kb);
            if n == 0
                line_error(mfilename, 'BMMAPt: D0 must be a non-empty cell of segment matrices');
            end
            B = numel(D1kb{1});
            if B == 0
                line_error(mfilename, 'BMMAPt: mark 1 declares no batch size');
            end
            for c = 1:K
                if numel(D1kb{c}) ~= B
                    line_error(mfilename, sprintf(['BMMAPt: mark %d declares %d batch sizes against %d in mark 1; ' ...
                        'the batch axis is dense and every mark must span it (use a zero block for an ' ...
                        'unused batch size)'], c, numel(D1kb{c}), B));
                end
                for b = 1:B
                    if ~iscell(D1kb{c}{b}) || numel(D1kb{c}{b}) ~= n
                        line_error(mfilename, sprintf(['BMMAPt: block (mark %d, batch %d) has %d segments ' ...
                            'against %d in D0; every block must be defined on the whole schedule'], ...
                            c, b, numel(D1kb{c}{b}), n));
                    end
                end
            end
            breakpoints = breakpoints(:).';
            if numel(breakpoints) ~= n + 1
                line_error(mfilename, 'BMMAPt: breakpoints must have one more entry than the number of segments');
            end
            if any(diff(breakpoints) <= 0)
                line_error(mfilename, 'BMMAPt: breakpoints must be strictly increasing');
            end
            h = size(D0{1}, 1);
            D1 = cell(1, n);
            D1k = cell(1, K);
            for c = 1:K
                D1k{c} = repmat({zeros(h, h)}, 1, n);
            end
            for j = 1:n
                if ~isequal(size(D0{j}), [h h])
                    line_error(mfilename, sprintf('BMMAPt: every D0 must be square of order %d; segment %d differs', h, j));
                end
                agg = zeros(h, h);
                for c = 1:K
                    perMark = zeros(h, h);
                    for b = 1:B
                        M = D1kb{c}{b}{j};
                        if ~isequal(size(M), [h h])
                            line_error(mfilename, sprintf(['BMMAPt: every block must be square of order %d; ' ...
                                'block (mark %d, batch %d) of segment %d differs'], h, c, b, j));
                        end
                        if any(M(:) < 0)
                            line_error(mfilename, sprintf(['BMMAPt: block (mark %d, batch %d) must be ' ...
                                'non-negative in segment %d'], c, b, j));
                        end
                        perMark = perMark + M;
                    end
                    D1k{c}{j} = perMark;
                    agg = agg + perMark;
                end
                off = D0{j} - diag(diag(D0{j}));
                if any(off(:) < 0)
                    line_error(mfilename, sprintf('BMMAPt: off-diagonal D0 entries must be non-negative in segment %d', j));
                end
                if any(abs(sum(D0{j} + agg, 2)) > 1e-10)
                    line_error(mfilename, sprintf(['BMMAPt: D0 plus every batch block must have zero row sums ' ...
                        '(generator) in segment %d'], j));
                end
                D1{j} = agg;
            end
            MAPt.checkCommonSupport(D0, true, 'off-diagonal D0');
            % The support rule binds the PER-MARK AGGREGATES and the total, NOT
            % the individual (mark, batch) blocks. Those three are exactly what
            % toMAPt and toMMAPt hand to the MAPt and MMAPt constructors, which
            % enforce the rule themselves, so checking them here is what keeps
            % the reductions constructible; checking each block as well would
            % add a constraint neither reduction needs and would forbid the one
            % thing this family exists to express, a batch composition that
            % changes with the segment (pairs in the morning, singles at night).
            for c = 1:K
                MAPt.checkCommonSupport(D1k{c}, false, sprintf('batch-aggregated D1 of mark %d', c));
            end
            MAPt.checkCommonSupport(D1, false, 'aggregate D1');
            if all(cellfun(@(M) sum(M(:)) <= 0, D1))
                line_error(mfilename, 'BMMAPt: every segment has zero arrival intensity, so no event can ever occur');
            end
            self.breakpoints = breakpoints;
            self.D0 = D0(:).';
            self.D1 = D1;
            self.D1k = D1k(:).';
            self.D1kb = D1kb(:).';
            self.cyclic = logical(cyclic);
            self.process = {self.breakpoints, self.D0, self.D1, self.cyclic, self.D1k, self.D1kb};
            self.sampleClock = breakpoints(1);
            self.samplePhase = 1;
            self.lastMark = 0;
            self.lastBatch = 0;
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
            % D = GETD1SEGMENTS(K) Batch-aggregated per-segment blocks of mark K.
            if nargin < 2
                d = self.D1;
            else
                self.assertMark(k);
                d = self.D1k{k};
            end
        end

        function d = getBatchSegments(self, k, b)
            % D = GETBATCHSEGMENTS(K, B) Per-segment blocks of mark K, batch size B.
            self.assertMark(k);
            self.assertBatch(b);
            d = self.D1kb{k}{b};
        end

        function d = getMarkSegments(self)
            % D = GETMARKSEGMENTS() The batch-aggregated mark blocks, mark-major.
            d = self.D1k;
        end

        function d = getBatchBlocks(self)
            % D = GETBATCHBLOCKS() All blocks, mark-major then batch then segment.
            d = self.D1kb;
        end

        function K = getNumberOfTypes(self)
            K = numel(self.D1kb);
        end

        function B = getMaxBatchSize(self)
            % B = GETMAXBATCHSIZE() Largest batch size the schedule declares.
            B = numel(self.D1kb{1});
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
            % The UNMARKED, UNBATCHED schedule, i.e. the MAPt whose D1 is the
            % per-segment aggregate. Hiding both labels is exact: an epoch of
            % the BMMAPt is an epoch of this process whatever it released.
            m = MAPt(self.breakpoints, self.D0, self.D1, self.cyclic);
        end

        function m = toMMAPt(self)
            % M = TOMMAPT()
            % The BATCH-BLIND marked schedule, i.e. the MMAPt whose mark blocks
            % are the per-mark aggregates over batch size. Every epoch keeps its
            % mark and releases one job.
            m = MMAPt(self.breakpoints, self.D0, self.D1k, self.cyclic);
        end

        function m = toMAPts(self, k)
            % M = TOMAPTS(K)
            % The MARGINAL schedule of 1-based mark K: epochs fire only on that
            % mark's blocks, while the other marks' transitions become hidden
            % phase changes. Segment by segment this is
            % MAPt(D0 + D1 - D1k, D1k), the time-varying analogue of
            % MarkedMAP.toMAPs.
            self.assertMark(k);
            n = numel(self.D0);
            hidden = cell(1, n);
            for j = 1:n
                hidden{j} = self.D0{j} + self.D1{j} - self.D1k{k}{j};
            end
            m = MAPt(self.breakpoints, hidden, self.D1k{k}, self.cyclic);
        end

        function bm = toBMAP(self)
            % BM = TOBMAP()
            % The width-weighted time average as a stationary BMAP, i.e. the
            % batch structure with the schedule averaged out and the marks
            % hidden. BMAP takes {D0, D_1, ..., D_B} with D_b releasing b jobs.
            D0bar = self.getTimeAverageProcess();
            B = self.getMaxBatchSize();
            blocks = cell(1, B + 1);
            blocks{1} = D0bar;
            for b = 1:B
                acc = zeros(size(D0bar));
                for c = 1:self.getNumberOfTypes()
                    acc = acc + self.getTimeAverageBatch(c, b);
                end
                blocks{1 + b} = acc;
            end
            bm = BMAP(blocks);
        end

        function [D0bar, D1bar] = getTimeAverageProcess(self)
            % [D0BAR, D1BAR] = GETTIMEAVERAGEPROCESS()
            % Width-weighted average pair over the horizon. This is the
            % stationary carrier of the phase structure where a solver needs a
            % time-homogeneous one; a convex combination of generators is a
            % generator, so it is itself a valid MAP.
            [D0bar, D1bar] = self.toMAPt().getTimeAverageProcess();
        end

        function M = getTimeAverageMark(self, k)
            % M = GETTIMEAVERAGEMARK(K)
            % Width-weighted average of mark K's batch-aggregated blocks.
            self.assertMark(k);
            M = self.widthAverage(self.D1k{k});
        end

        function M = getTimeAverageBatch(self, k, b)
            % M = GETTIMEAVERAGEBATCH(K, B)
            % Width-weighted average of the (mark K, batch size B) block.
            self.assertMark(k);
            self.assertBatch(b);
            M = self.widthAverage(self.D1kb{k}{b});
        end

        function rate = getTimeAverageRate(self)
            % RATE = GETTIMEAVERAGERATE()
            % EVENT rate of the time-averaged aggregate MAP, i.e. batch epochs
            % per unit time. This is NOT the job rate when any batch exceeds
            % one; see getTimeAverageJobRate.
            [D0bar, D1bar] = self.getTimeAverageProcess();
            rate = map_lambda({D0bar, D1bar});
        end

        function rate = getTimeAverageJobRate(self)
            % RATE = GETTIMEAVERAGEJOBRATE()
            % JOBS per unit time, sum_c sum_b b*rate(c,b). This is the quantity
            % a station's throughput balances against, and it exceeds
            % getTimeAverageRate whenever a batch larger than one has mass.
            rate = sum(self.getTimeAverageBatchRates() .* (1:self.getMaxBatchSize()));
        end

        function lam = getTimeAverageMarkRates(self)
            % LAM = GETTIMEAVERAGEMARKRATES()
            % Per-mark EVENT rates of the time-averaged process. These sum to
            % getTimeAverageRate, which is the identity a marked stream has to
            % satisfy: labelling the epochs cannot change how many there are.
            % Through the canonical mmap_lambda rather than a second
            % implementation of theta*D1k*e: the rate of a marked counting
            % process is weighted by the STATIONARY phase distribution map_prob,
            % not by the embedded departure law map_pie, and the two differ
            % whenever the phases do not all arrive at the same rate.
            K = self.getNumberOfTypes();
            blocks = cell(1, K);
            for c = 1:K
                blocks{c} = self.getTimeAverageMark(c);
            end
            lam = reshape(self.markedRates(blocks), 1, K);
        end

        function lam = getTimeAverageMarkJobRates(self)
            % LAM = GETTIMEAVERAGEMARKJOBRATES()
            % Per-mark JOB rates, sum_b b*rate(c,b). These sum to
            % getTimeAverageJobRate.
            K = self.getNumberOfTypes();
            B = self.getMaxBatchSize();
            blocks = cell(1, K * B);
            for c = 1:K
                for b = 1:B
                    blocks{(c - 1) * B + b} = self.getTimeAverageBatch(c, b);
                end
            end
            flat = self.markedRates(blocks);
            lam = zeros(1, K);
            for c = 1:K
                for b = 1:B
                    lam(c) = lam(c) + b * flat((c - 1) * B + b);
                end
            end
        end

        function rates = getTimeAverageBatchRates(self)
            % RATES = GETTIMEAVERAGEBATCHRATES()
            % RATES(b) is the EVENT rate of batches of size b, marks hidden.
            % These sum to getTimeAverageRate.
            B = self.getMaxBatchSize();
            K = self.getNumberOfTypes();
            blocks = cell(1, B);
            for b = 1:B
                acc = zeros(self.getNumberOfPhases());
                for c = 1:K
                    acc = acc + self.getTimeAverageBatch(c, b);
                end
                blocks{b} = acc;
            end
            rates = reshape(self.markedRates(blocks), 1, B);
        end

        function m = getMeanBatchSize(self)
            % M = GETMEANBATCHSIZE() Mean jobs per epoch of the time average.
            rates = self.getTimeAverageBatchRates();
            total = sum(rates);
            if total > 0
                m = sum(rates .* (1:numel(rates))) / total;
            else
                m = 0;
            end
        end

        function mean = getMean(self)
            % MEAN = GETMEAN()
            % Palm mean interval BETWEEN EPOCHS of the time-averaged aggregate
            % MAP, matching BMAP, whose inter-batch process is the aggregate.
            mean = 1.0 / self.getTimeAverageRate();
        end

        function scv = getSCV(self)
            % SCV = GETSCV()
            % NaN: a BMMAP_t is neither renewal nor time-homogeneous, so there
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
            self.lastBatch = 0;
        end

        function m = getLastMark(self)
            % M = GETLASTMARK() 1-based mark of the interval last sampled.
            m = self.lastMark;
        end

        function b = getLastBatch(self)
            % B = GETLASTBATCH() Batch size of the interval last sampled.
            b = self.lastBatch;
        end

        function [X, C, Bsz] = sample(self, n)
            % [X, C, BSZ] = SAMPLE(N)
            % Draws N successive inter-epoch times along ONE sample path, with
            % their marks C and batch sizes BSZ. Both the intensity and the
            % phase depend on absolute time, so this advances an internal clock
            % and phase across calls. Use resetSampleClock to restart. A
            % non-cyclic schedule that runs out returns 0 for every remaining
            % sample.
            if nargin < 2
                n = 1;
            end
            X = zeros(n, 1);
            C = zeros(n, 1);
            Bsz = zeros(n, 1);
            for i = 1:n
                [interval, phase, mark, batch] = self.nextArrival(self.sampleClock, self.samplePhase);
                X(i) = interval;
                if interval <= 0
                    break % horizon exhausted: no further epoch can occur
                end
                C(i) = mark;
                Bsz(i) = batch;
                self.sampleClock = self.sampleClock + interval;
                self.samplePhase = phase;
                self.lastMark = mark;
                self.lastBatch = batch;
            end
        end

        function [x, phase, mark, batch] = nextArrival(self, from, phase)
            % [X, PHASE, MARK, BATCH] = NEXTARRIVAL(FROM, PHASE)
            % Time to the next epoch from wall clock FROM in the given phase,
            % the phase after it, the 1-based mark it carries and how many jobs
            % it releases. Exact: within a segment the phase process is a
            % homogeneous CTMC, and by the memoryless property the residual
            % holding time may be redrawn at a breakpoint, so the boundary is
            % crossed by advancing the clock and resampling under the new
            % matrices.
            %
            % NEITHER LABEL COSTS AN EXTRA DRAW: the competing transitions are
            % accumulated destination-major, then mark-minor, then batch-minor,
            % so the running total after all (mark, batch) pairs of a
            % destination equals the aggregate total after that destination. The
            % winning destination is therefore the one the unlabelled walk would
            % choose for the same uniform. Putting batch INSIDE mark is what
            % makes a B = 1 BMMAPt reproduce the MMAPt sample path for sample
            % path, exactly as mark-inside-destination makes a K = 1 MMAPt
            % reproduce the MAPt.
            elapsed = 0.0;
            pos = from;
            mark = 0;
            batch = 0;
            h = self.getNumberOfPhases();
            K = self.getNumberOfTypes();
            B = self.getMaxBatchSize();
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
                        for b = 1:B
                            cum = cum + self.D1kb{c}{b}{idx}(phase, j);
                            if u < cum
                                x = elapsed;
                                phase = j;
                                mark = c;
                                batch = b;
                                return
                            end
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
                    % Rounding left u at or past the total: fall back to the
                    % last block with any mass out of this phase, as MMAPt does.
                    for j = h:-1:1
                        for c = K:-1:1
                            for b = B:-1:1
                                if self.D1kb{c}{b}{idx}(phase, j) > 0
                                    x = elapsed;
                                    phase = j;
                                    mark = c;
                                    batch = b;
                                    return
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    methods (Access = private)
        function assertMark(self, k)
            if k < 1 || k > numel(self.D1kb)
                line_error(mfilename, sprintf('BMMAPt: mark index out of range: %d', k));
            end
        end

        function assertBatch(self, b)
            if b < 1 || b > self.getMaxBatchSize()
                line_error(mfilename, sprintf('BMMAPt: batch size out of range: %d', b));
            end
        end

        function M = widthAverage(self, segs)
            % M = WIDTHAVERAGE(SEGS) Width-weighted mean of a 1-by-n segment cell.
            widths = diff(self.breakpoints);
            total = sum(widths);
            M = zeros(size(self.D0{1}));
            for j = 1:numel(segs)
                M = M + (widths(j) / total) * segs{j};
            end
        end

        function lam = markedRates(self, blocks)
            % LAM = MARKEDRATES(BLOCKS)
            % Per-block rates of the time-averaged process, through the
            % canonical mmap_lambda so that the stationary phase weighting is
            % the one M3A uses rather than a second implementation of it. The
            % blocks must partition the aggregate D1bar.
            [D0bar, D1bar] = self.getTimeAverageProcess();
            cell_mmap = cell(1, numel(blocks) + 2);
            cell_mmap{1} = D0bar;
            cell_mmap{2} = D1bar;
            for c = 1:numel(blocks)
                cell_mmap{2 + c} = blocks{c};
            end
            lam = mmap_lambda(cell_mmap);
        end
    end

    methods (Static)
        function bm = fromMMAPtWithBatchPMF(mmapt, batchSizes, pmf)
            % BM = FROMMMAPTWITHBATCHPMF(MMAPT, BATCHSIZES, PMF)
            %
            % Split every mark block of an MMAPt across a batch-size law, the
            % time-varying analogue of BMAP.fromMAPWithBatchPMF. The batch size
            % is independent of the mark, of the phase and of the segment;
            % build the blocks directly when it is not.
            if ~isa(mmapt, 'MMAPt')
                line_error(mfilename, 'BMMAPt.fromMMAPtWithBatchPMF: the first argument must be an MMAPt');
            end
            if numel(batchSizes) ~= numel(pmf)
                line_error(mfilename, 'BMMAPt.fromMMAPtWithBatchPMF: batch sizes and PMF must have the same length');
            end
            if any(pmf < 0) || sum(pmf) <= 0
                line_error(mfilename, 'BMMAPt.fromMMAPtWithBatchPMF: the PMF must be non-negative with positive mass');
            end
            pmf = pmf(:).' / sum(pmf);
            batchSizes = batchSizes(:).';
            if any(batchSizes < 1) || any(batchSizes ~= round(batchSizes))
                line_error(mfilename, 'BMMAPt.fromMMAPtWithBatchPMF: batch sizes must be positive integers');
            end
            B = max(batchSizes);
            K = mmapt.getNumberOfTypes();
            n = mmapt.getNumSegments();
            h = mmapt.getNumberOfPhases();
            blocks = cell(1, K);
            for c = 1:K
                segs = mmapt.getD1Segments(c);
                blocks{c} = cell(1, B);
                for b = 1:B
                    blocks{c}{b} = repmat({zeros(h, h)}, 1, n);
                end
                for i = 1:numel(batchSizes)
                    b = batchSizes(i);
                    for j = 1:n
                        blocks{c}{b}{j} = blocks{c}{b}{j} + pmf(i) * segs{j};
                    end
                end
            end
            bm = BMMAPt(mmapt.getBreakpoints(), mmapt.getD0Segments(), blocks, mmapt.isCyclic());
        end
    end
end
