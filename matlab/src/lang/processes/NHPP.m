classdef NHPP < ContinuousDistribution
    % NHPP Non-homogeneous Poisson process with a piecewise-constant intensity.
    %
    % The intensity is a step function of the wall clock: segment i covers
    % [breakpoints(i), breakpoints(i+1)) and carries rate rates(i), so
    % breakpoints has one more entry than rates. With cyclic=true the schedule
    % repeats, giving a cyclic Poisson process.
    %
    % Two horizon conventions:
    %   cyclic     : the schedule repeats with period
    %                T = breakpoints(end) - breakpoints(1); the active segment
    %                at time t follows from mod(t - breakpoints(1), T).
    %   non-cyclic : the intensity is zero outside
    %                [breakpoints(1), breakpoints(end)), so the process emits
    %                nothing once the schedule is exhausted. A non-cyclic NHPP
    %                is therefore a transient construct: run to steady state it
    %                converges to the empty system, so callers should use a time
    %                span within the horizon.
    %
    % This is NOT a renewal process. Successive intervals are dependent, because
    % the position within the schedule carries over from one event to the next.
    % Accordingly the scalar summaries that presuppose an i.i.d. interval
    % distribution -- getSCV, getSkewness, evalCDF, evalLST -- are undefined and
    % return NaN rather than a representative exponential value, which would
    % silently misreport the process as Poisson. The schedule is the
    % parameterisation: read it with getRateSchedule. getMean is well defined and
    % returns the arrival-stationary (Palm) mean interval 1/timeAverageRate.
    %
    % The process representation stores:
    %   process{1} : 1-by-(n+1) row vector of breakpoints
    %   process{2} : 1-by-n row vector of rates
    %   process{3} : logical, true if cyclic
    %
    % Solver support. The LDES simulation engine honours the exact schedule in
    % both steady state (cyclic only) and transient analysis. SolverFLD honours
    % it in getTranAvg, by injecting the intensity as a time-varying rate
    % multiplier on the closing ODE; SolverFLD.getAvg uses the time-average
    % rate, which is the steady state of a cyclic schedule. Every other solver
    % rejects a model using it via the standard unsupported-feature check.
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    properties
        breakpoints; % 1-by-(n+1) segment boundaries, strictly increasing
        rates;       % 1-by-n non-negative rate on each segment
        cyclic;      % logical, whether the schedule repeats
        process;     % {breakpoints, rates, cyclic}, for serialization
        sampleClock; % wall-clock position of the next sample (see sample)
    end

    methods
        function self = NHPP(breakpoints, rates, cyclic)
            % SELF = NHPP(BREAKPOINTS, RATES, CYCLIC)
            self@ContinuousDistribution('NHPP', 0, [0, Inf]);
            if nargin < 3 || isempty(cyclic)
                cyclic = true;
            end
            if nargin < 2 || isempty(breakpoints) || isempty(rates) ...
                    || numel(breakpoints) ~= numel(rates) + 1
                line_error(mfilename, 'NHPP: breakpoints must be non-empty with one more entry than rates');
            end
            breakpoints = breakpoints(:).';
            rates = rates(:).';
            if any(diff(breakpoints) <= 0)
                line_error(mfilename, 'NHPP: breakpoints must be strictly increasing');
            end
            if any(rates < 0) || any(isinf(rates))
                line_error(mfilename, 'NHPP: rates must be finite and non-negative');
            end
            if sum(rates .* diff(breakpoints)) <= 0
                line_error(mfilename, 'NHPP: the schedule has zero total intensity, so no event can ever occur');
            end
            self.breakpoints = breakpoints;
            self.rates = rates;
            self.cyclic = logical(cyclic);
            self.process = {self.breakpoints, self.rates, self.cyclic};
            self.sampleClock = breakpoints(1);
            self.mean = 1.0 / self.getTimeAverageRate();
            self.immediate = false;
        end

        function b = getBreakpoints(self)
            b = self.breakpoints;
        end

        function r = getRates(self)
            r = self.rates;
        end

        function c = isCyclic(self)
            c = self.cyclic;
        end

        function n = getNumSegments(self)
            n = numel(self.rates);
        end

        function T = getPeriod(self)
            % T = GETPERIOD() Horizon length, which is the period when cyclic.
            T = self.breakpoints(end) - self.breakpoints(1);
        end

        function rate = getTimeAverageRate(self)
            % RATE = GETTIMEAVERAGERATE()
            % sum(rates.*widths)/sum(widths). For a cyclic schedule this is the
            % long-run arrival rate; for a non-cyclic one it averages over the
            % active horizon only, the intensity being zero afterwards.
            rate = sum(self.rates .* diff(self.breakpoints)) / self.getPeriod();
        end

        function r = getRateAt(self, t)
            % R = GETRATEAT(T) Rate in force at T; zero past a non-cyclic horizon.
            T = self.getPeriod();
            offset = t - self.breakpoints(1);
            if self.cyclic
                offset = mod(offset, T);
            elseif offset < 0 || offset >= T
                r = 0.0;
                return
            end
            pos = self.breakpoints(1) + offset;
            idx = find(pos < self.breakpoints(2:end), 1, 'first');
            if isempty(idx)
                idx = numel(self.rates);
            end
            r = self.rates(idx);
        end

        function sched = getRateSchedule(self)
            % SCHED = GETRATESCHEDULE()
            % The parameterisation of the process; the scalar interval summaries
            % are not. Model compilation recognises a schedule-bearing process by
            % this method rather than by class name.
            sched = struct('breakpoints', self.breakpoints, ...
                'rates', self.rates, 'cyclic', self.cyclic);
        end

        function mean = getMean(self)
            % MEAN = GETMEAN() Arrival-stationary (Palm) mean interval.
            mean = 1.0 / self.getTimeAverageRate();
        end

        function scv = getSCV(self)
            % SCV = GETSCV()
            % NaN: an NHPP is not a renewal process, so there is no i.i.d.
            % interval distribution for an SCV to summarise. Returning a
            % representative value here would report a time-varying process as an
            % exponential one to every consumer of sn.scv.
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
            % Draws N successive interarrival times along ONE sample path. The
            % intensity depends on absolute time, so this advances an internal
            % clock across calls: consecutive samples form a realisation of the
            % process starting at breakpoints(1), not independent draws from a
            % marginal. Use resetSampleClock to restart. A non-cyclic schedule
            % that runs out returns 0 for every remaining sample, the intensity
            % there being zero.
            if nargin < 2
                n = 1;
            end
            X = zeros(n, 1);
            for i = 1:n
                residual = -log(1 - rand());
                interval = self.nextInterval(self.sampleClock, residual);
                X(i) = interval;
                if interval <= 0
                    break % horizon exhausted: no further event can occur
                end
                self.sampleClock = self.sampleClock + interval;
            end
        end

        function x = nextInterval(self, from, residual)
            % X = NEXTINTERVAL(FROM, RESIDUAL)
            % Solves int_{from}^{from+x} lambda(u) du = RESIDUAL for X, walking
            % the schedule forward and consuming the budget segment by segment.
            % Returns 0 when a non-cyclic horizon is exhausted first, which
            % callers read as "no further event".
            %
            % Exact for an NHPP: conditional on no event since the last one, the
            % residual is governed by the intensity from the current instant
            % onward, so a holding time drawn under a rate that has since changed
            % is not a sample from this process.
            T = self.getPeriod();
            offset = from - self.breakpoints(1);
            if self.cyclic
                offset = mod(offset, T);
            elseif offset >= T
                x = 0.0;
                return
            elseif offset < 0
                offset = 0.0;
            end
            pos = self.breakpoints(1) + offset;
            idx = 1;
            while idx < numel(self.rates) && pos >= self.breakpoints(idx+1)
                idx = idx + 1;
            end
            elapsed = 0.0;
            while true
                remainingInSegment = self.breakpoints(idx+1) - pos;
                massInSegment = self.rates(idx) * remainingInSegment;
                % The rate guard also keeps a zero-rate segment from dividing 0/0
                % on the measure-zero draw residual == 0.
                if self.rates(idx) > 0 && massInSegment >= residual
                    x = elapsed + residual / self.rates(idx);
                    return
                end
                residual = residual - massInSegment;
                elapsed = elapsed + remainingInSegment;
                idx = idx + 1;
                if idx > numel(self.rates)
                    if ~self.cyclic
                        x = 0.0;
                        return
                    end
                    idx = 1;
                end
                pos = self.breakpoints(idx);
            end
        end
    end
end
