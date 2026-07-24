function bool = convergedStoch(self, it)
% BOOL = CONVERGEDSTOCH(IT)
%
% Convergence controller for stochastic layer solvers (Robbins-Monro mode).
%
% When one or more layer solvers return noisy estimates (simulation, e.g.
% JMT/SSA/LDES, or Monte Carlo integration, e.g. NC with mci/imci/ls), the
% deterministic Picard iteration in converged.m cannot terminate: the
% successive-difference error is bounded below by the standard error of the
% layer estimates, and the layer-reset confirmation step merely resamples
% the noise. This routine implements a stochastic approximation iteration:
%
%  1. Burn-in: for the first stochiter_burnin iterations the plain Picard
%     iteration runs with the relaxation factor configured at init, to move
%     quickly toward the fixed point.
%  2. Robbins-Monro step: afterwards the relaxation factor applied by
%     updateMetrics to the fed-forward iterate (servt, residt, tput,
%     callservt) decays as omega_k = a0/k^alpha with alpha in (0.5,1].
%     Under the contraction assumption already made by the deterministic
%     iteration, and zero-mean noise with bounded variance, the iterate
%     converges almost surely to the true fixed point (Robbins and Monro,
%     1951). Layer seeds are rotated per iteration in pre() so successive
%     evaluations observe independent noise.
%  3. Polyak-Ruppert averaging: running averages of the layer results and
%     of the reported iterates are maintained and installed as the final
%     solution in finish(), giving the optimal O(1/sqrt(k)) rate and
%     robustness to the choice of a0 (Polyak and Juditsky, 1992).
%  4. Stopping: iteration stops when the drift of the averaged results
%     stays below iter_tol for stochiter_conseq consecutive iterations.
%     The drift of a running average decays like 1/k even under persistent
%     noise, so the test terminates, and it self-calibrates: larger noise
%     keeps the drift above tolerance longer, forcing more averaging.
%
% Note: the Robbins-Monro step acts through relax_omega, which is applied
% by updateMetricsDefault; the moment3 update path does not use relaxation,
% so this controller is primarily intended for method 'default'.

bool = false;
if it < 1
    return
end
E = self.nlayers;
burnin = self.options.config.stochiter_burnin;
a0 = self.options.config.stochiter_a0;
alpha = self.options.config.stochiter_alpha;

% Schedule the Robbins-Monro step used by updateMetrics at the next iteration
if it >= burnin
    self.relax_omega = min(1.0, a0 / max(1, it - burnin + 1)^alpha);
end

if it <= burnin
    % pure Picard burn-in; no averaging or convergence testing yet
    self.maxitererr(it) = Inf;
    if self.options.verbose
        line_printf(sprintf('Stochastic iteration burn-in %d/%d.', it, burnin));
    end
    return
end

if isempty(self.stochiter_start)
    self.stochiter_start = it;
    if self.options.verbose
        line_printf('\b Started Robbins-Monro averaging (stochastic layer solvers detected).');
    end
end

%% Polyak-Ruppert update of the layer result averages and drift metric
k = self.stoch_avg_count + 1;
err = 0;
fields = {'QN','UN','RN','TN','AN','WN'};
for e = 1:E
    raw = self.results{end,e};
    avg = struct();
    if k == 1
        for f = 1:length(fields)
            avg.(fields{f}) = raw.(fields{f});
        end
    else
        prev = self.stoch_avg{e};
        for f = 1:length(fields)
            avg.(fields{f}) = polyak(prev.(fields{f}), raw.(fields{f}), k);
        end
        % drift of the averaged queue lengths, normalized by population
        N = sum(self.ensemble{e}.getNumberOfJobs);
        if N > 0
            d = abs(avg.QN(:) - prev.QN(:));
            d(isnan(d)) = 0;
            err = err + max(d)/N;
        end
    end
    self.stoch_avg{e} = avg;
end
self.stoch_avg_count = k;

% Polyak-Ruppert averages of the fed-forward iterates used in reporting
if k == 1
    self.stoch_servt_avg = self.servt;
    self.stoch_residt_avg = self.residt;
else
    self.stoch_servt_avg = polyak(self.stoch_servt_avg, self.servt, k);
    self.stoch_residt_avg = polyak(self.stoch_residt_avg, self.residt, k);
end

self.maxitererr(it) = err;
if self.options.verbose
    line_printf(sprintf('RMIterErr=%.6e (tol=%.6e, omega=%.3f, k=%d)', ...
        err, self.options.iter_tol, self.relax_omega, k));
end

%% Stop when the averaged-iterate drift stays below tolerance
conseq = self.options.config.stochiter_conseq;
if k > conseq
    bool = all(self.maxitererr(it-conseq+1:it) < self.options.iter_tol);
    self.hasconverged = bool;
end
end

function m = polyak(prevv, raww, k)
% M = POLYAK(PREVV, RAWW, K)
% Running-mean update m = prev + (raw - prev)/k, robust to NaN entries
% in either operand (a NaN sample leaves the average untouched).
m = prevv + (raww - prevv)/k;
bad = isnan(m);
m(bad) = raww(bad);
bad = isnan(m);
m(bad) = prevv(bad);
end
