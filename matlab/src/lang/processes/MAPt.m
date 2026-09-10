classdef MAPt < ContinuousDistribution
    % MAPt Time-inhomogeneous Markovian arrival process (MAP_t).
    %
    % Following Ko and Pender (Oper. Res. Lett. 45, 2017), a MAP_t is an
    % ordinary MAP whose two matrices are functions of the wall clock, D0(t) and
    % D1(t), required only to be locally integrable. This class realises that
    % definition with a piecewise-constant schedule, which is dense in L1_loc and
    % is the form that serialises: segment k covers
    % [breakpoints(k), breakpoints(k+1)) and carries the pair (D0{k}, D1{k}), so
    % breakpoints has one more entry than the matrix lists. D0 holds transition
    % rates without an arrival, D1 the rates that generate one, and D0+D1 is a
    % generator in every segment.
    %
    % Two horizon conventions, as for NHPP:
    %   cyclic     : the schedule repeats with period
    %                T = breakpoints(end) - breakpoints(1).
    %   non-cyclic : outside [breakpoints(1), breakpoints(end)) the process is
    %                frozen in its last phase and emits nothing, so a non-cyclic
    %                MAP_t is a transient construct.
    %
    % Setting h = 1 with D0 = -lambda_k, D1 = lambda_k recovers exactly the NHPP
    % with the same breakpoints and rates.
    %
    % This is neither a renewal process nor a time-homogeneous one, so the scalar
    % summaries that presuppose an i.i.d. interval distribution -- getSCV,
    % getSkewness, evalCDF, evalLST -- are undefined and return NaN rather than a
    % representative value that would misreport the process as stationary. The
    % schedule is the parameterisation: read it with getRateSchedule.
    %
    % The class deliberately does NOT extend Markovian. Code gated on
    % isMarkovian reads getProcess as a single stationary (D0, D1) pair and would
    % silently drop the schedule; NHPP avoids the base class for the same reason.
    %
    % Constant support. The fluid solver expresses a segment as a per-entry
    % multiplier on the time-averaged nominal, so the constructor requires the
    % sparsity pattern of the matrices to be identical across segments. A
    % schedule that switches a transition on or off is refused outright rather
    % than silently losing it.
    %
    % The process representation stores:
    %   process{1} : 1-by-(n+1) row vector of breakpoints
    %   process{2} : 1-by-n cell of D0 matrices
    %   process{3} : 1-by-n cell of D1 matrices
    %   process{4} : logical, true if cyclic
    %
    % Solver support. SolverFLD honours the schedule in getTranAvg, and its 'kp'
    % method integrates the Ko-Pender fluid and diffusion limits. Every other
    % solver rejects a model using it via the standard unsupported-feature check.
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    properties
        breakpoints; % 1-by-(n+1) segment boundaries, strictly increasing
        D0;          % 1-by-n cell of no-arrival rate matrices
        D1;          % 1-by-n cell of arrival-generating rate matrices
        cyclic;      % logical, whether the schedule repeats
        process;     % {breakpoints, D0, D1, cyclic}, for serialization
        sampleClock; % wall-clock position of the next sample (see sample)
        samplePhase; % phase of the modulating chain at sampleClock
    end

    methods
        function self = MAPt(breakpoints, D0, D1, cyclic)
            % SELF = MAPT(BREAKPOINTS, D0, D1, CYCLIC)
            self@ContinuousDistribution('MAPt', 0, [0, Inf]);
            if nargin < 4 || isempty(cyclic)
                cyclic = true;
            end
            if ~iscell(D0), D0 = {D0}; end
            if ~iscell(D1), D1 = {D1}; end
            n = numel(D0);
            if n == 0 || numel(D1) ~= n
                line_error(mfilename, 'MAPt: D0 and D1 must be non-empty cells of equal length');
            end
            breakpoints = breakpoints(:).';
            if numel(breakpoints) ~= n + 1
                line_error(mfilename, 'MAPt: breakpoints must have one more entry than the number of segments');
            end
            if any(diff(breakpoints) <= 0)
                line_error(mfilename, 'MAPt: breakpoints must be strictly increasing');
            end
            h = size(D0{1}, 1);
            for k = 1:n
                if ~isequal(size(D0{k}), [h h]) || ~isequal(size(D1{k}), [h h])
                    line_error(mfilename, sprintf('MAPt: every D0 and D1 must be square of order %d; segment %d differs', h, k));
                end
                if any(D1{k}(:) < 0)
                    line_error(mfilename, sprintf('MAPt: D1 must be non-negative in segment %d', k));
                end
                off = D0{k} - diag(diag(D0{k}));
                if any(off(:) < 0)
                    line_error(mfilename, sprintf('MAPt: off-diagonal D0 entries must be non-negative in segment %d', k));
                end
                if any(abs(sum(D0{k} + D1{k}, 2)) > 1e-10)
                    line_error(mfilename, sprintf('MAPt: D0+D1 must have zero row sums (generator) in segment %d', k));
                end
            end
            MAPt.checkCommonSupport(D0, true, 'off-diagonal D0');
            MAPt.checkCommonSupport(D1, false, 'D1');
            if all(cellfun(@(M) sum(M(:)) <= 0, D1))
                line_error(mfilename, 'MAPt: every segment has zero arrival intensity, so no event can ever occur');
            end
            self.breakpoints = breakpoints;
            self.D0 = D0(:).';
            self.D1 = D1(:).';
            self.cyclic = logical(cyclic);
            self.process = {self.breakpoints, self.D0, self.D1, self.cyclic};
            self.sampleClock = breakpoints(1);
            self.samplePhase = 1;
            self.mean = 1.0 / self.getTimeAverageRate();
            self.immediate = false;
        end

        function b = getBreakpoints(self)
            b = self.breakpoints;
        end

        function d = getD0Segments(self)
            d = self.D0;
        end

        function d = getD1Segments(self)
            d = self.D1;
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

        function M = getD0At(self, t)
            % M = GETD0AT(T) D0 in force at T; zero past a non-cyclic horizon.
            idx = self.getSegmentIndexAt(t);
            if idx == 0
                M = zeros(size(self.D0{1}));
            else
                M = self.D0{idx};
            end
        end

        function M = getD1At(self, t)
            % M = GETD1AT(T) D1 in force at T; zero past a non-cyclic horizon.
            idx = self.getSegmentIndexAt(t);
            if idx == 0
                M = zeros(size(self.D1{1}));
            else
                M = self.D1{idx};
            end
        end

        function [D0bar, D1bar] = getTimeAverageProcess(self)
            % [D0BAR, D1BAR] = GETTIMEAVERAGEPROCESS()
            % Width-weighted average pair over the horizon. This is the
            % stationary carrier of the phase structure where a solver needs a
            % time-homogeneous one; a convex combination of generators is a
            % generator, so it is itself a valid MAP.
            widths = diff(self.breakpoints);
            total = sum(widths);
            D0bar = zeros(size(self.D0{1}));
            D1bar = zeros(size(self.D1{1}));
            for k = 1:numel(self.D0)
                D0bar = D0bar + (widths(k) / total) * self.D0{k};
                D1bar = D1bar + (widths(k) / total) * self.D1{k};
            end
        end

        function rate = getTimeAverageRate(self)
            % RATE = GETTIMEAVERAGERATE() Arrival rate of the time-averaged MAP.
            % For h = 1 this is exactly the NHPP width-weighted average intensity.
            [D0bar, D1bar] = self.getTimeAverageProcess();
            rate = map_lambda({D0bar, D1bar});
        end

        function r = getRateAt(self, t)
            % R = GETRATEAT(T) Arrival rate of the MAP in force at T, zero past a
            % non-cyclic horizon. This is the stationary rate of that segment, not
            % the instantaneous conditional intensity, which depends on the phase.
            idx = self.getSegmentIndexAt(t);
            if idx == 0
                r = 0.0;
            else
                r = map_lambda({self.D0{idx}, self.D1{idx}});
            end
        end

        function sched = getRateSchedule(self)
            % SCHED = GETRATESCHEDULE()
            % The parameterisation of the process; the scalar interval summaries
            % are not. Model compilation recognises a schedule-bearing process by
            % this method rather than by class name.
            sched = struct('breakpoints', self.breakpoints, 'D0', {self.D0}, ...
                'D1', {self.D1}, 'cyclic', self.cyclic);
        end

        function mean = getMean(self)
            % MEAN = GETMEAN() Palm mean interval of the time-averaged MAP.
            mean = 1.0 / self.getTimeAverageRate();
        end

        function scv = getSCV(self)
            % SCV = GETSCV()
            % NaN: a MAP_t is neither renewal nor time-homogeneous, so there is
            % no i.i.d. interval distribution for an SCV to summarise. Returning
            % the SCV of the time-averaged MAP would report a time-varying
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
            self.samplePhase = 1;
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
                [interval, phase] = self.nextArrival(self.sampleClock, self.samplePhase);
                X(i) = interval;
                if interval <= 0
                    break % horizon exhausted: no further arrival can occur
                end
                self.sampleClock = self.sampleClock + interval;
                self.samplePhase = phase;
            end
        end

        function [x, phase] = nextArrival(self, from, phase)
            % [X, PHASE] = NEXTARRIVAL(FROM, PHASE)
            % Time to the next arrival from wall clock FROM in the given phase,
            % and the phase after that arrival. Exact: within a segment the phase
            % process is a homogeneous CTMC, and by the memoryless property the
            % residual holding time may be redrawn at a breakpoint, so the
            % boundary is crossed by advancing the clock and resampling under the
            % new matrices. Returns 0 when a non-cyclic horizon is exhausted.
            elapsed = 0.0;
            pos = from;
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
                Dz = self.D0{idx};
                Da = self.D1{idx};
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
                % Competing transitions out of the current phase, arrivals first.
                weights = [Da(phase, :), Dz(phase, :)];
                weights(h + phase) = 0;
                u = rand() * total;
                j = find(cumsum(weights) > u, 1, 'first');
                if isempty(j)
                    j = find(weights > 0, 1, 'last');
                end
                if j <= h
                    x = elapsed;
                    phase = j;
                    return
                end
                phase = j - h;
            end
        end
    end

    methods (Static)
        function checkCommonSupport(mats, ignoreDiagonal, label)
            % CHECKCOMMONSUPPORT(MATS, IGNOREDIAGONAL, LABEL)
            % Reject a schedule whose matrices do not share one sparsity pattern.
            % The fluid solver expresses a segment as a per-entry multiplier on a
            % nominal matrix, and that multiplier is undefined where the nominal
            % entry is zero.
            ref = mats{1} ~= 0;
            if ignoreDiagonal
                ref(logical(eye(size(ref)))) = false;
            end
            for k = 2:numel(mats)
                pat = mats{k} ~= 0;
                if ignoreDiagonal
                    pat(logical(eye(size(pat)))) = false;
                end
                if ~isequal(pat, ref)
                    line_error(mfilename, sprintf(['MAPt/PHt: the %s sparsity pattern must be identical ' ...
                        'across segments; segment %d differs from segment 1. A schedule that switches a ' ...
                        'transition on or off cannot be expressed as a per-entry multiplier on the ' ...
                        'time-averaged process. Keep the entry present with a small positive rate instead.'], ...
                        label, k));
                end
            end
        end
    end
end
