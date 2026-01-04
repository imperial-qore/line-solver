function bool = converged(self, it)
% BOOL = CONVERGED(IT)
%
% Apply convergence test to the SolverLN iterations. As the solver keeps
% iterating, this method maintains a moving avg of the recent results based
% on which it averages across the layer the maximum queue-length error.
% Convergence is tested by resetting all layers (to avoid caching) and
% doing an extra iteration. If the iteration keeps fulfilling the error
% requirements for convergence, then the solver completes.

bool = false;
iter_min = max([2*length(self.model.ensemble),ceil(self.options.iter_max/4)]);
E = self.nlayers;
results = self.results; %#ok<NASGU> % faster in matlab

%% Start moving average to help convergence

if false%it>self.averagingstart%<self.averagingstart+50
    % In the first 50 averaging iterations use Cesaro summation
    if ~isempty(self.averagingstart)
        if it>=iter_min % assume steady-state
            for e=1:E
                wnd_size_max = (it-self.averagingstart+1);
                sk_q = cell(1,wnd_size_max);
                sk_u = cell(1,wnd_size_max);
                sk_r = cell(1,wnd_size_max);
                sk_t = cell(1,wnd_size_max);
                sk_a = cell(1,wnd_size_max);
                sk_w = cell(1,wnd_size_max);
                % compute all partial sumbs of up to wnd_size_max elements
                for k= 1:wnd_size_max
                    if k==1
                        sk_q{k} = results{self.averagingstart,e}.QN;
                        sk_u{k} = results{self.averagingstart,e}.UN;
                        sk_r{k} = results{self.averagingstart,e}.RN;
                        sk_t{k} = results{self.averagingstart,e}.TN;
                        sk_a{k} = results{self.averagingstart,e}.AN;
                        sk_w{k} = results{self.averagingstart,e}.WN;
                    else
                        sk_q{k} = results{self.averagingstart+k-1,e}.QN/k + sk_q{k-1}*(k-1)/k;
                        sk_u{k} = results{self.averagingstart+k-1,e}.UN/k + sk_u{k-1}*(k-1)/k;
                        sk_r{k} = results{self.averagingstart+k-1,e}.RN/k + sk_r{k-1}*(k-1)/k;
                        sk_t{k} = results{self.averagingstart+k-1,e}.TN/k + sk_t{k-1}*(k-1)/k;
                        sk_a{k} = results{self.averagingstart+k-1,e}.AN/k + sk_a{k-1}*(k-1)/k;
                        sk_w{k} = results{self.averagingstart+k-1,e}.WN/k + sk_w{k-1}*(k-1)/k;
                    end
                end
                results{end,e}.QN = cellsum(sk_q)/wnd_size_max;
                results{end,e}.UN = cellsum(sk_u)/wnd_size_max;
                results{end,e}.RN = cellsum(sk_r)/wnd_size_max;
                results{end,e}.TN = cellsum(sk_t)/wnd_size_max;
                results{end,e}.AN = cellsum(sk_a)/wnd_size_max;
                results{end,e}.WN = cellsum(sk_w)/wnd_size_max;
            end
        end
    end
else
    wnd_size = max(5,ceil(iter_min/5)); % moving window size
    mov_avg_weight = 1/wnd_size;
    results = self.results; % faster in matlab
    if it>=iter_min % assume steady-state
        for e=1:E
            results{end,e}.QN = mov_avg_weight*results{end,e}.QN;
            results{end,e}.UN = mov_avg_weight*results{end,e}.UN;
            results{end,e}.RN = mov_avg_weight*results{end,e}.RN;
            results{end,e}.TN = mov_avg_weight*results{end,e}.TN;
            results{end,e}.AN = mov_avg_weight*results{end,e}.AN;
            results{end,e}.WN = mov_avg_weight*results{end,e}.WN;
            for k=1:(wnd_size-1)
                results{end,e}.QN = results{end,e}.QN + results{end-k,e}.QN * mov_avg_weight;
                results{end,e}.UN = results{end,e}.UN + results{end-k,e}.UN * mov_avg_weight;
                results{end,e}.RN = results{end,e}.RN + results{end-k,e}.RN * mov_avg_weight;
                results{end,e}.TN = results{end,e}.TN + results{end-k,e}.TN * mov_avg_weight;
                results{end,e}.AN = results{end,e}.AN + results{end-k,e}.AN * mov_avg_weight;
                results{end,e}.WN = results{end,e}.WN + results{end-k,e}.WN * mov_avg_weight;
            end
        end
    end
end
self.results = results;

%% Take as error metric the max qlen-error averaged across layers
if it>1
    self.maxitererr(it) = 0;
    for e = 1:E
        metric = results{end,e}.QN;
        metric_1 = results{end-1,e}.QN;
        N = sum(self.ensemble{e}.getNumberOfJobs);
        if N>0
            try
                IterErr = max(abs(metric(:) - metric_1(:)))/N;
            catch
                IterErr = 0;
            end
            self.maxitererr(it) = self.maxitererr(it) + IterErr;
        end
        % if self.options.verbose
        %     if self.solvers{e}.options.verbose
        %         line_printf(sprintf('QLen change: %f.\n',self.maxitererr(it)/E));
        %     elseif e==1
        %         line_printf('\n');
        %     end
        % end
    end
    if it==iter_min
        if self.options.verbose
            line_printf( '\b Started averaging to aid convergence.');
        end
        self.averagingstart = it;
    end

    %% Update relaxation factor for adaptive/auto modes
    relax_mode = self.options.config.relax;
    if strcmpi(relax_mode, 'adaptive') || strcmpi(relax_mode, 'auto')
        % Track error history
        self.relax_err_history = [self.relax_err_history, self.maxitererr(it)];
        wnd = self.options.config.relax_history;
        if length(self.relax_err_history) > wnd
            self.relax_err_history = self.relax_err_history(end-wnd+1:end);
        end

        if length(self.relax_err_history) >= 3
            % Detect oscillation by counting sign changes in error differences
            err = self.relax_err_history;
            diff_err = diff(err);
            sign_changes = sum(diff_err(1:end-1) .* diff_err(2:end) < 0);

            if strcmpi(relax_mode, 'auto') && self.relax_omega == 1.0
                % For 'auto' mode: enable relaxation when oscillation detected
                if sign_changes >= length(diff_err) * 0.5
                    self.relax_omega = self.options.config.relax_factor;
                    % Debug output removed
                    if self.options.verbose
                        % line_printf(sprintf(' [enabling relaxation, omega=%.2f]', self.relax_omega));
                    end
                end
            elseif strcmpi(relax_mode, 'adaptive')
                % For 'adaptive' mode: adjust omega based on error trajectory
                if sign_changes >= length(diff_err) * 0.5
                    % Oscillating - reduce omega
                    self.relax_omega = max(self.options.config.relax_min, self.relax_omega * 0.8);
                    % Debug output removed
                    if self.options.verbose
                        % line_printf(sprintf(' [omega=%.2f]', self.relax_omega));
                    end
                elseif sign_changes == 0 && self.maxitererr(it) < self.maxitererr(it-1)
                    % Monotonically decreasing - can increase omega slightly
                    self.relax_omega = min(1.0, self.relax_omega * 1.05);
                end
            end
        end
    end
end

%% Check convergence. Do not allow to converge in less than 2 iterations.
if it==0 && self.options.verbose
    % Debug output removed
elseif it>2 && self.maxitererr(it) < self.options.iter_tol && self.maxitererr(it-1) < self.options.iter_tol&& self.maxitererr(it-1) < self.options.iter_tol
    if ~self.hasconverged % if potential convergence has just been detected
        % do a hard reset of every layer to check that this is really the fixed point
        for e=1:E
            self.ensemble{e}.reset();
        end
        if self.options.verbose
           % Debug output removed
        end
        self.hasconverged = true; % if it passes the change again next time then complete
    else
        if self.options.verbose
            if self.solvers{end}.options.verbose
                % Debug output removed
            end
        end
        bool = true;
    end
else
    self.hasconverged = false;
end
end
