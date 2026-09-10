classdef PHt < ContinuousDistribution
    % PHt Time-inhomogeneous phase-type distribution (Ph_t).
    %
    % Following Ko and Pender (Oper. Res. Lett. 45, 2017), a Ph_t is an ordinary
    % phase-type distribution whose initial vector and sub-generator are
    % functions of the wall clock, alpha(t) and S(t), required only to be locally
    % integrable. This class realises that definition with a piecewise-constant
    % schedule: segment k covers [breakpoints(k), breakpoints(k+1)) and carries
    % the pair (alpha{k}, S{k}), so breakpoints has one more entry than the
    % lists. The exit vector is s(t) = -S(t)e.
    %
    % Because both the phase and the elapsed service depend on absolute time, a
    % Ph_t service time is a function of the epoch at which service starts:
    % sampleFrom(t0) is the operative sampler, and sample walks one path.
    %
    % Setting h = 1 with S = -mu_k recovers a time-varying exponential, whose
    % completion stream at a saturated server is the NHPP with rates mu_k.
    %
    % Like MAPt this does NOT extend Markovian, so that isMarkovian-gated code
    % cannot read it as a single stationary (alpha, S) pair; and the scalar
    % summaries getSCV, getSkewness, evalCDF, evalLST return NaN, the
    % distribution of a service time being different at every start epoch.
    %
    % Constant support. The fluid solver expresses a segment as a per-entry
    % multiplier on the time-averaged nominal, so the constructor requires the
    % sparsity pattern of alpha, of the off-diagonal S and of the exit vector to
    % be identical across segments.
    %
    % The process representation stores:
    %   process{1} : 1-by-(n+1) row vector of breakpoints
    %   process{2} : 1-by-n cell of alpha row vectors
    %   process{3} : 1-by-n cell of S matrices
    %   process{4} : logical, true if cyclic
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    properties
        breakpoints; % 1-by-(n+1) segment boundaries, strictly increasing
        alpha;       % 1-by-n cell of initial probability row vectors
        S;           % 1-by-n cell of sub-generators
        cyclic;      % logical, whether the schedule repeats
        process;     % {breakpoints, alpha, S, cyclic}, for serialization
        sampleClock; % wall-clock position of the next sample (see sample)
    end

    methods
        function self = PHt(breakpoints, alpha, S, cyclic)
            % SELF = PHT(BREAKPOINTS, ALPHA, S, CYCLIC)
            self@ContinuousDistribution('PHt', 0, [0, Inf]);
            if nargin < 4 || isempty(cyclic)
                cyclic = true;
            end
            if ~iscell(alpha), alpha = {alpha}; end
            if ~iscell(S), S = {S}; end
            n = numel(S);
            if n == 0 || numel(alpha) ~= n
                line_error(mfilename, 'PHt: alpha and S must be non-empty cells of equal length');
            end
            breakpoints = breakpoints(:).';
            if numel(breakpoints) ~= n + 1
                line_error(mfilename, 'PHt: breakpoints must have one more entry than the number of segments');
            end
            if any(diff(breakpoints) <= 0)
                line_error(mfilename, 'PHt: breakpoints must be strictly increasing');
            end
            h = size(S{1}, 1);
            for k = 1:n
                alpha{k} = alpha{k}(:).';
                if ~isequal(size(S{k}), [h h]) || numel(alpha{k}) ~= h
                    line_error(mfilename, sprintf('PHt: every S must be square of order %d with a matching alpha; segment %d differs', h, k));
                end
                if any(alpha{k} < 0) || abs(sum(alpha{k}) - 1) > 1e-10
                    line_error(mfilename, sprintf('PHt: alpha must be a probability vector in segment %d', k));
                end
                off = S{k} - diag(diag(S{k}));
                if any(off(:) < 0)
                    line_error(mfilename, sprintf('PHt: off-diagonal S entries must be non-negative in segment %d', k));
                end
                if any(-sum(S{k}, 2) < -1e-10)
                    line_error(mfilename, sprintf('PHt: S must have non-positive row sums in segment %d', k));
                end
            end
            MAPt.checkCommonSupport(S, true, 'off-diagonal S');
            MAPt.checkCommonSupport(cellfun(@(M) -sum(M, 2), S, 'UniformOutput', false), false, 'exit vector');
            MAPt.checkCommonSupport(alpha, false, 'alpha');
            if all(cellfun(@(M) sum(-sum(M, 2)) <= 0, S))
                line_error(mfilename, 'PHt: every segment has zero exit rate, so service can never complete');
            end
            self.breakpoints = breakpoints;
            self.alpha = alpha(:).';
            self.S = S(:).';
            self.cyclic = logical(cyclic);
            self.process = {self.breakpoints, self.alpha, self.S, self.cyclic};
            self.sampleClock = breakpoints(1);
            self.mean = 1.0 / self.getTimeAverageRate();
            self.immediate = false;
        end

        function b = getBreakpoints(self)
            b = self.breakpoints;
        end

        function a = getAlphaSegments(self)
            a = self.alpha;
        end

        function s = getSSegments(self)
            s = self.S;
        end

        function c = isCyclic(self)
            c = self.cyclic;
        end

        function n = getNumSegments(self)
            n = numel(self.S);
        end

        function h = getNumberOfPhases(self)
            h = size(self.S{1}, 1);
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
                idx = numel(self.S);
            end
        end

        function a = getAlphaAt(self, t)
            % A = GETALPHAAT(T) alpha in force at T; the last segment's vector
            % past a non-cyclic horizon, where getSAt is zero anyway.
            idx = self.getSegmentIndexAt(t);
            if idx == 0
                a = self.alpha{end};
            else
                a = self.alpha{idx};
            end
        end

        function M = getSAt(self, t)
            % M = GETSAT(T) S in force at T; zero past a non-cyclic horizon.
            idx = self.getSegmentIndexAt(t);
            if idx == 0
                M = zeros(size(self.S{1}));
            else
                M = self.S{idx};
            end
        end

        function [abar, Sbar] = getTimeAverageProcess(self)
            % [ABAR, SBAR] = GETTIMEAVERAGEPROCESS()
            % Width-weighted average over the horizon. A convex combination of
            % sub-generators is a sub-generator and of probability vectors a
            % probability vector, so the nominal is a valid phase-type.
            widths = diff(self.breakpoints);
            total = sum(widths);
            abar = zeros(size(self.alpha{1}));
            Sbar = zeros(size(self.S{1}));
            for k = 1:numel(self.S)
                abar = abar + (widths(k) / total) * self.alpha{k};
                Sbar = Sbar + (widths(k) / total) * self.S{k};
            end
        end

        function [D0bar, D1bar] = getTimeAverageProcessMAP(self)
            % [D0BAR, D1BAR] = GETTIMEAVERAGEPROCESSMAP()
            % The nominal as a (D0, D1) pair, D1 = s*alpha, for the fluid carrier.
            [abar, Sbar] = self.getTimeAverageProcess();
            D0bar = Sbar;
            D1bar = (-sum(Sbar, 2)) * abar;
        end

        function rate = getTimeAverageRate(self)
            % RATE = GETTIMEAVERAGERATE() Completion rate of the time-averaged PH.
            [abar, Sbar] = self.getTimeAverageProcess();
            rate = 1.0 / (-abar * (Sbar \ ones(size(Sbar, 1), 1)));
        end

        function r = getRateAt(self, t)
            % R = GETRATEAT(T) Completion rate of the PH in force at T.
            idx = self.getSegmentIndexAt(t);
            if idx == 0
                r = 0.0;
            else
                Sm = self.S{idx};
                r = 1.0 / (-self.alpha{idx} * (Sm \ ones(size(Sm, 1), 1)));
            end
        end

        function sched = getRateSchedule(self)
            % SCHED = GETRATESCHEDULE()
            % The parameterisation of the process; the scalar interval summaries
            % are not. Model compilation recognises a schedule-bearing process by
            % this method rather than by class name.
            sched = struct('breakpoints', self.breakpoints, 'alpha', {self.alpha}, ...
                'S', {self.S}, 'cyclic', self.cyclic);
        end

        function mean = getMean(self)
            % MEAN = GETMEAN() Mean of the time-averaged phase-type.
            mean = 1.0 / self.getTimeAverageRate();
        end

        function scv = getSCV(self)
            % SCV = GETSCV()
            % NaN: the service-time distribution differs at every start epoch, so
            % there is no single i.i.d. law for an SCV to summarise. Returning the
            % SCV of the time-averaged representation would report a time-varying
            % process as a stationary one to every consumer of sn.scv.
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
        end

        function X = sample(self, n)
            % X = SAMPLE(N)
            % Draws N successive service times along ONE sample path: sample i
            % starts where sample i-1 completed, not at a fixed epoch. Use
            % resetSampleClock to restart, or sampleFrom to draw a service time
            % starting at a chosen epoch.
            if nargin < 2
                n = 1;
            end
            X = zeros(n, 1);
            for i = 1:n
                interval = self.sampleFrom(self.sampleClock);
                X(i) = interval;
                if interval <= 0
                    break % horizon exhausted: service can never complete
                end
                self.sampleClock = self.sampleClock + interval;
            end
        end

        function x = sampleFrom(self, t0)
            % X = SAMPLEFROM(T0)
            % Service time for a job whose service starts at wall clock T0.
            % Exact: within a segment the phase process is a homogeneous
            % absorbing CTMC, and by the memoryless property the residual holding
            % time may be redrawn at a breakpoint. The initial phase is drawn from
            % alpha in force at T0. Returns 0 when a non-cyclic horizon is
            % exhausted before absorption.
            idx = self.getSegmentIndexAt(t0);
            if idx == 0
                x = 0.0;
                return
            end
            a = self.alpha{idx};
            phase = find(cumsum(a) > rand() * sum(a), 1, 'first');
            if isempty(phase)
                phase = numel(a);
            end
            elapsed = 0.0;
            pos = t0;
            h = self.getNumberOfPhases();
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
                Sm = self.S{idx};
                total = -Sm(phase, phase);
                if total <= 0
                    if ~self.cyclic && idx == numel(self.S)
                        x = 0.0;
                        return
                    end
                    elapsed = elapsed + toBoundary;
                    pos = pos + toBoundary;
                    continue
                end
                holding = -log(1 - rand()) / total;
                if holding >= toBoundary
                    if ~self.cyclic && idx == numel(self.S)
                        x = 0.0;
                        return
                    end
                    elapsed = elapsed + toBoundary;
                    pos = pos + toBoundary;
                    continue
                end
                elapsed = elapsed + holding;
                pos = pos + holding;
                % Competing transitions: absorption first, then phase changes.
                weights = [-sum(Sm(phase, :)), Sm(phase, :)];
                weights(1 + phase) = 0;
                u = rand() * total;
                j = find(cumsum(weights) > u, 1, 'first');
                if isempty(j)
                    j = find(weights > 0, 1, 'last');
                end
                if j == 1
                    x = elapsed;
                    return
                end
                phase = j - 1;
            end
        end
    end
end
