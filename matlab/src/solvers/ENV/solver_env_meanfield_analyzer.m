function varargout = solver_env_meanfield_analyzer(self, phase, it, e)
% SOLVER_ENV_MEANFIELD_ANALYZER  Mean-field (marginal mean queue-length) coupling for SolverENV.
%
% This analyzer implements the default (and historical) SolverENV coupling:
% every stage is solved transiently and reduced to its marginal mean queue
% lengths, which are carried across environment switches via initFromMarginal.
% This mean-field collapse discards the joint state distribution at the handoff;
% see solver_env_statevec_analyzer for the full state-vector alternative. It is
% the default analyzer (options.method other than 'statevec').
%
% Phase dispatch (called by the SolverENV EnsembleSolver hooks):
%   'pre'       pre_(self,it)            -> []
%   'analyze'   analyze_(self,it,e)      -> [results_e, runtime]
%   'post'      post_(self,it)           -> []
%   'finish'    finish_(self)            -> []
%   'converged' converged_(self,it)      -> bool
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

switch phase
    case 'pre'
        pre_(self, it);
    case 'analyze'
        [varargout{1}, varargout{2}] = analyze_(self, it, e);
    case 'post'
        post_(self, it);
    case 'finish'
        finish_(self);
    case 'converged'
        varargout{1} = converged_(self, it);
    otherwise
        line_error(mfilename, sprintf('Unknown meanfield-analyzer phase: %s', phase));
end
end

function bool = converged_(self, it) % convergence test at iteration it
            % BOOL = CONVERGED(IT) % CONVERGENCE TEST AT ITERATION IT
            % Computes max relative absolute difference of queue lengths between iterations
            % Aligned with JAR SolverEnv.converged() implementation

            bool = false;
            if it <= 1
                return
            end

            E = self.getNumberOfModels;
            M = self.ensemble{1}.getNumberOfStations;
            K = self.ensemble{1}.getNumberOfClasses;

            % Check convergence per class (aligned with JAR structure)
            for k = 1:K
                % Build QEntry and QExit matrices (M x E) for this class
                QEntry = zeros(M, E);
                QExit = zeros(M, E);
                for e = 1:E
                    % Skip stages where analysis failed
                    res_curr = self.results{it,e};
                    res_prev = self.results{it-1,e};
                    for i = 1:M
                        if ~isempty(res_curr) && isfield(res_curr, 'Tran') ...
                                && isstruct(res_curr.Tran.Avg) && isfield(res_curr.Tran.Avg, 'Q')
                            Qik_curr = res_curr.Tran.Avg.Q{i,k};
                            if isstruct(Qik_curr) && isfield(Qik_curr, 'metric') && ~isempty(Qik_curr.metric)
                                QExit(i,e) = Qik_curr.metric(1);
                            end
                        end
                        if ~isempty(res_prev) && isfield(res_prev, 'Tran') ...
                                && isstruct(res_prev.Tran.Avg) && isfield(res_prev.Tran.Avg, 'Q')
                            Qik_prev = res_prev.Tran.Avg.Q{i,k};
                            if isstruct(Qik_prev) && isfield(Qik_prev, 'metric') && ~isempty(Qik_prev.metric)
                                QEntry(i,e) = Qik_prev.metric(1);
                            end
                        end
                    end
                end

                % Compute max relative absolute difference using maxpe
                % maxpe computes max(abs(1 - approx./exact)) = max(abs((approx-exact)./exact))
                % This matches JAR's Matrix.maxAbsDiff() implementation
                maxDiff = maxpe(QExit(:), QEntry(:));
                if isempty(maxDiff)
                    maxDiff = 0;  % Handle case where all QEntry values are zero
                end
                if isnan(maxDiff) || isinf(maxDiff)
                    return  % Non-convergence on invalid values
                end
                if maxDiff >= self.options.iter_tol
                    return  % Not converged
                end
            end
            bool = true;
end

function pre_(self, it)
            % PRE(IT)

            if it==1
                for e=list(self)
                    try
                        if isinf(self.getSolver(e).options.timespan(2))
                            [QN,~,~,~] = self.getSolver(e).getAvg();
                        else
                            [QNt,~,~] = self.getSolver(e).getTranAvg();
                            % Handle case where getTranAvg returns NaN for disabled/empty results
                            QN = zeros(size(QNt));
                            for i = 1:size(QNt,1)
                                for k = 1:size(QNt,2)
                                    if isstruct(QNt{i,k}) && isfield(QNt{i,k}, 'metric')
                                        QN(i,k) = QNt{i,k}.metric(end);
                                    else
                                        QN(i,k) = 0; % Default to 0 for NaN/invalid entries
                                    end
                                end
                            end
                        end
                        if ~isa(self.solvers{e}, 'SolverFluid') && ~isa(self.ensemble{e}, 'LayeredNetwork')
                            QN = self.roundMarginalForDiscreteSolver(QN, self.sn{e});
                        end
                        self.ensemble{e}.initFromMarginal(QN);
                    catch
                        % Skip pre-initialization for stages where getAvg fails
                    end
                end
            end
end

function [results_e, runtime] = analyze_(self, it, e)
            % [RESULTS_E, RUNTIME] = ANALYZE(IT, E)
            results_e = struct();
            results_e.('Tran') = struct();
            results_e.Tran.('Avg') = [];
            T0 = tic;
            runtime = toc(T0);
            %% initialize
            try
                [Qt,Ut,Tt] = self.ensemble{e}.getTranHandles;
                self.solvers{e}.reset();
                [QNt,UNt,TNt] = self.solvers{e}.getTranAvg(Qt,Ut,Tt);
                results_e.Tran.Avg.Q = QNt;
                results_e.Tran.Avg.U = UNt;
                results_e.Tran.Avg.T = TNt;
            catch
                % Skip analysis for stages where transient analysis fails
            end
end

function post_(self, it)
            % POST(IT)

            E = self.getNumberOfModels;
            M = self.ensemble{1}.getNumberOfStations;
            K = self.ensemble{1}.getNumberOfClasses;
            [detSojourn, dvals] = sojournConfig_(self, E);
            for e=1:E
                if isempty(self.results{it,e}) || ~isfield(self.results{it,e}, 'Tran') ...
                        || ~isstruct(self.results{it,e}.Tran.Avg) || ~isfield(self.results{it,e}.Tran.Avg, 'Q')
                    for h = 1:E
                        Qexit{e,h} = zeros(M, K);
                        Uexit{e,h} = zeros(M, K);
                        Texit{e,h} = zeros(M, K);
                    end
                    continue
                end
                if detSojourn
                    % Deterministic sojourn: exit metrics are the transient at
                    % t=d_e. see _kb/06-solver-catalog.md for rationale
                    Qd = zeros(M, K); Ud = zeros(M, K); Td = zeros(M, K);
                    for i=1:size(self.results{it,e}.Tran.Avg.Q,1)
                        for r=1:size(self.results{it,e}.Tran.Avg.Q,2)
                            Qir = self.results{it,e}.Tran.Avg.Q{i,r};
                            if isstruct(Qir) && isfield(Qir,'t') && isfield(Qir,'metric') && ~isempty(Qir.t)
                                Qd(i,r) = detEval_(Qir.t, Qir.metric, dvals(e));
                                Uir = self.results{it,e}.Tran.Avg.U{i,r};
                                if isstruct(Uir) && isfield(Uir,'metric') && ~isempty(Uir.metric)
                                    Ud(i,r) = detEval_(Qir.t, Uir.metric, dvals(e));
                                end
                                Tir = self.results{it,e}.Tran.Avg.T{i,r};
                                if isstruct(Tir) && isfield(Tir,'metric') && ~isempty(Tir.metric)
                                    Td(i,r) = detEval_(Qir.t, Tir.metric, dvals(e));
                                end
                            end
                        end
                    end
                    for h = 1:E
                        Qexit{e,h} = Qd; Uexit{e,h} = Ud; Texit{e,h} = Td;
                    end
                    continue
                end
                for h = 1:E
                    Qexit{e,h} = zeros(size(self.results{it,e}.Tran.Avg.Q));
                    Uexit{e,h} = zeros(size(self.results{it,e}.Tran.Avg.U));
                    Texit{e,h} = zeros(size(self.results{it,e}.Tran.Avg.T));
                    for i=1:size(self.results{it,e}.Tran.Avg.Q,1)
                        for r=1:size(self.results{it,e}.Tran.Avg.Q,2)
                            Qir = self.results{it,e}.Tran.Avg.Q{i,r};
                            % Check if result is a valid struct with required fields
                            if isstruct(Qir) && isfield(Qir, 't') && isfield(Qir, 'metric') && ~isempty(Qir.t)
                                w{e,h} = [0, map_cdf(self.envObj.proc{e}{h}, Qir.t(2:end)) - map_cdf(self.envObj.proc{e}{h}, Qir.t(1:end-1))]';
                                if ~isnan(w{e,h})
                                    Qexit{e,h}(i,r) = Qir.metric'*w{e,h}/sum(w{e,h});
                                    Uir = self.results{it,e}.Tran.Avg.U{i,r};
                                    if isstruct(Uir) && isfield(Uir, 'metric') && ~isempty(Uir.metric)
                                        Uexit{e,h}(i,r) = Uir.metric'*w{e,h}/sum(w{e,h});
                                    end
                                    Tir = self.results{it,e}.Tran.Avg.T{i,r};
                                    if isstruct(Tir) && isfield(Tir, 'metric') && ~isempty(Tir.metric)
                                        Texit{e,h}(i,r) = Tir.metric'*w{e,h}/sum(w{e,h});
                                    end
                                else
                                    Qexit{e,h}(i,r) = 0;
                                    Uexit{e,h}(i,r) = 0;
                                    Texit{e,h}(i,r) = 0;
                                end
                            else
                                w{e,h} = 0;
                            end
                        end
                    end
                end
            end

            Qentry = cell(1,E); % average entry queue-length
            for e = 1:E
                % Skip stages where analysis failed (no valid Tran results)
                if isempty(self.results{it,e}) || ~isfield(self.results{it,e}, 'Tran') ...
                        || ~isstruct(self.results{it,e}.Tran.Avg) || ~isfield(self.results{it,e}.Tran.Avg, 'Q')
                    continue
                end
                Qentry{e} = zeros(size(Qexit{e}));
                for h=1:E
                    % probability of coming from h to e \times resetFun(Qexit from h to e
                    if self.envObj.probOrig(h,e) > 0
                        Qentry{e} = Qentry{e} + self.envObj.probOrig(h,e) * self.resetFromMarginal{h,e}(Qexit{h,e});
                    end
                end
                if ~isa(self.solvers{e}, 'SolverFluid') && ~isa(self.ensemble{e}, 'LayeredNetwork')
                    Qentry{e} = self.roundMarginalForDiscreteSolver(Qentry{e}, self.sn{e});
                end
                self.solvers{e}.reset();
                self.ensemble{e}.initFromMarginal(Qentry{e});
            end

            % Update transition rates between stages if state-dependent method
            if isfield(self.options, 'method') && strcmp(self.options.method, 'statedep')
                for e = 1:E
                    for h = 1:E
                        if ~isa(self.envObj.env{e,h}, 'Disabled') && ~isempty(self.resetEnvRates{e,h})
                            self.envObj.env{e,h} = self.resetEnvRates{e,h}(...
                                self.envObj.env{e,h}, Qexit{e,h}, Uexit{e,h}, Texit{e,h});
                        end
                    end
                end
                % Reinitialize environment after rate updates
                self.envObj.init();
            end
end

function finish_(self)
            % FINISH()

            it = size(self.results,1); % use last iteration
            E = self.getNumberOfModels;
            M = self.ensemble{1}.getNumberOfStations;
            K = self.ensemble{1}.getNumberOfClasses;
            [detSojourn, dvals] = sojournConfig_(self, E);
            for e=1:E
                QExit{e}=zeros(M, K);
                UExit{e}=zeros(M, K);
                TExit{e}=zeros(M, K);
                if it>0 && ~isempty(self.results{it,e}) && isfield(self.results{it,e}, 'Tran') ...
                        && isstruct(self.results{it,e}.Tran.Avg) && isfield(self.results{it,e}.Tran.Avg, 'Q')
                    for i=1:size(self.results{it,e}.Tran.Avg.Q,1)
                        for r=1:size(self.results{it,e}.Tran.Avg.Q,2)
                            Qir = self.results{it,e}.Tran.Avg.Q{i,r};
                            % Check if result is a valid struct with required fields
                            if isstruct(Qir) && isfield(Qir, 't') && isfield(Qir, 'metric') && ~isempty(Qir.t)
                                if detSojourn
                                    % Deterministic sojourn: evaluate at t=d_e.
                                    QExit{e}(i,r) = detEval_(Qir.t, Qir.metric, dvals(e));
                                    Uir = self.results{it,e}.Tran.Avg.U{i,r};
                                    if isstruct(Uir) && isfield(Uir, 'metric') && ~isempty(Uir.metric)
                                        UExit{e}(i,r) = detEval_(Qir.t, Uir.metric, dvals(e));
                                    else
                                        UExit{e}(i,r) = 0;
                                    end
                                    Tir = self.results{it,e}.Tran.Avg.T{i,r};
                                    if isstruct(Tir) && isfield(Tir, 'metric') && ~isempty(Tir.metric)
                                        TExit{e}(i,r) = detEval_(Qir.t, Tir.metric, dvals(e));
                                    else
                                        TExit{e}(i,r) = 0;
                                    end
                                    continue
                                end
                                w{e} = [0, map_cdf(self.envObj.holdTime{e}, Qir.t(2:end)) - map_cdf(self.envObj.holdTime{e}, Qir.t(1:end-1))]';
                                QExit{e}(i,r) = Qir.metric'*w{e}/sum(w{e});
                                Uir = self.results{it,e}.Tran.Avg.U{i,r};
                                if isstruct(Uir) && isfield(Uir, 'metric') && ~isempty(Uir.metric)
                                    UExit{e}(i,r) = Uir.metric'*w{e}/sum(w{e});
                                else
                                    UExit{e}(i,r) = 0;
                                end
                                Tir = self.results{it,e}.Tran.Avg.T{i,r};
                                if isstruct(Tir) && isfield(Tir, 'metric') && ~isempty(Tir.metric)
                                    TExit{e}(i,r) = Tir.metric'*w{e}/sum(w{e});
                                else
                                    TExit{e}(i,r) = 0;
                                end
                            else
                                QExit{e}(i,r) = 0;
                                UExit{e}(i,r) = 0;
                                TExit{e}(i,r) = 0;
                            end
                        end
                    end
                end
                %                 for h = 1:E
                %                     QE{e,h} = zeros(size(self.results{it,e}.Tran.Avg.Q));
                %                     for i=1:size(self.results{it,e}.Tran.Avg.Q,1)
                %                         for r=1:size(self.results{it,e}.Tran.Avg.Q,2)
                %                             w{e,h} = [0, map_cdf(self.envObj.proc{e}{h}, self.results{it,e}.Tran.Avg.Q{i,r}(2:end,2)) - map_cdf(self.envObj.proc{e}{h}, self.results{it,e}.Tran.Avg.Q{i,r}(1:end-1,2))]';
                %                             if ~isnan(w{e,h})
                %                                 QE{e,h}(i,r) = self.results{it,e}.Tran.Avg.Q{i,r}(:,1)'*w{e,h}/sum(w{e,h});
                %                             else
                %                                 QE{e,h}(i,r) = 0;
                %                             end
                %                         end
                %                     end
                %                 end
            end

            Qval=0*QExit{e};
            Uval=0*UExit{e};
            Tval=0*TExit{e};
            for e=1:E
                Qval = Qval + self.envObj.probEnv(e) * QExit{e}; % to check
                Uval = Uval + self.envObj.probEnv(e) * UExit{e}; % to check
                Tval = Tval + self.envObj.probEnv(e) * TExit{e}; % to check
            end
            self.result.Avg.Q = Qval;
            %    self.result.Avg.R = R;
            %    self.result.Avg.X = X;
            self.result.Avg.U = Uval;
            self.result.Avg.T = Tval;
            %    self.result.Avg.C = C;
            %self.result.runtime = runtime;
            %if self.options.verbose
            %    line_printf('\n');
            %end

            % Cache-hit aggregation across the environment (probEnv-weighted,
            % sojourn-averaged), written onto the reference model's cache nodes.
            % see _kb/06-solver-catalog.md for rationale
            aggregateCacheMeanfield_(self, it, E);
end

function aggregateCacheMeanfield_(self, it, E) %#ok<INUSD>
            % Cache-hit aggregation for fluid inner solvers. Runs a self-
            % contained mean-field fixed point that carries each cache's mean
            % occupancy across environment switches (the cache analog of the
            % queue-length initFromMarginal handoff): stage e is integrated over
            % its sojourn from an entry occupancy that mixes the exit occupancy
            % of its predecessors by probOrig. At convergence the per-class hit
            % throughput is the probEnv-weighted, sojourn-averaged arrival x
            % hit-prob, and the reported hit ratio is hit/(hit+miss).
            ref = self.ensemble{1};
            K = ref.getNumberOfClasses;

            % Only fluid stages expose the RMF cache transient used here.
            for e = 1:E
                if ~isa(self.solvers{e}, 'SolverFluid')
                    return
                end
            end
            sn1 = self.sn{1};
            if ~isfield(sn1, 'nodetype')
                return
            end
            cacheNodes = find(sn1.nodetype == NodeType.Cache);
            if isempty(cacheNodes)
                return
            end
            ncaches = numel(cacheNodes);

            % Finite integration window per stage (fall back to a few mean
            % holding times when the inner solver left the timespan open).
            tspan = cell(1, E);
            for e = 1:E
                ts = self.solvers{e}.options.timespan;
                if isempty(ts) || any(~isfinite(ts))
                    ts = [0, 20 * map_mean(self.envObj.holdTime{e})];
                end
                tspan{e} = ts;
            end

            entryOcc = repmat({cell(1, ncaches)}, 1, E); % entryOcc{e}{c}
            hp = cell(1, E); mp = cell(1, E); ar = cell(1, E);
            tc = cell(1, E); wmass = cell(1, E); exitOcc = cell(1, E);
            maxSweep = max(1, self.options.iter_max);
            tol = self.options.iter_tol;
            prevEntryFlat = [];
            for sweep = 1:maxSweep
                for e = 1:E
                    opts_e = self.solvers{e}.options;
                    opts_e.method = 'rmf';
                    opts_e.timespan = tspan{e};
                    [tc{e}, hp{e}, mp{e}, cnodes, ar{e}, xocc] = ...
                        solver_fld_cacheqn_tran(self.sn{e}, opts_e, entryOcc{e});
                    tt = tc{e}(:)';
                    % Sojourn-weighted exit occupancy per cache for the handoff
                    % (RMF drift mapped from per-request to real time).
                    % see _kb/06-solver-catalog.md for rationale
                    exitOcc{e} = cell(1, ncaches);
                    wmass{e} = cell(1, ncaches);
                    for cc = 1:ncaches
                        cidx = find(cnodes == cacheNodes(cc), 1);
                        if isempty(cidx)
                            continue
                        end
                        Lam = sum(ar{e}(cidx, :));
                        if ~(Lam > 0)
                            continue
                        end
                        treal = tt / Lam;
                        wm = [0, map_cdf(self.envObj.holdTime{e}, treal(2:end)) ...
                                 - map_cdf(self.envObj.holdTime{e}, treal(1:end-1))];
                        wmass{e}{cc} = wm(:);
                        sw = sum(wm);
                        if ~isempty(xocc{cidx}) && sw > 0
                            exitOcc{e}{cc} = xocc{cidx} * wmass{e}{cc} / sw;
                        end
                    end
                end
                % Update each stage's entry occupancy from its predecessors.
                newEntry = repmat({cell(1, ncaches)}, 1, E);
                for e = 1:E
                    for cc = 1:ncaches
                        acc = [];
                        for h = 1:E
                            po = self.envObj.probOrig(h, e);
                            if po > 0 && ~isempty(exitOcc{h}{cc})
                                if isempty(acc)
                                    acc = po * exitOcc{h}{cc};
                                else
                                    acc = acc + po * exitOcc{h}{cc};
                                end
                            end
                        end
                        newEntry{e}{cc} = acc;
                    end
                end
                entryFlat = [];
                for e = 1:E
                    for cc = 1:ncaches
                        entryFlat = [entryFlat; newEntry{e}{cc}(:)]; %#ok<AGROW>
                    end
                end
                entryOcc = newEntry;
                if ~isempty(prevEntryFlat) && numel(prevEntryFlat) == numel(entryFlat)
                    if norm(entryFlat - prevEntryFlat, inf) < tol
                        prevEntryFlat = entryFlat;
                        break
                    end
                end
                prevEntryFlat = entryFlat;
            end

            % Aggregate the converged hit/miss throughputs and write results.
            for cc = 1:ncaches
                hitT = zeros(1, K); missT = zeros(1, K);
                for e = 1:E
                    if numel(wmass{e}) < cc || isempty(wmass{e}{cc})
                        continue
                    end
                    w = wmass{e}{cc};
                    sw = sum(w);
                    if ~(sw > 0) || any(~isfinite(w))
                        continue
                    end
                    cidx = cc; % cnodes ordering matches cacheNodes across stages
                    pe = self.envObj.probEnv(e);
                    for k = 1:K
                        a = ar{e}(cidx, k);
                        if a <= 0
                            continue
                        end
                        hbar = reshape(hp{e}(cidx,k,:), 1, []) * w / sw;
                        mbar = reshape(mp{e}(cidx,k,:), 1, []) * w / sw;
                        hitT(k)  = hitT(k)  + pe * a * hbar;
                        missT(k) = missT(k) + pe * a * mbar;
                    end
                end
                hitprob = NaN(1, K); missprob = NaN(1, K);
                for k = 1:K
                    tot = hitT(k) + missT(k);
                    if tot > 0
                        hitprob(k)  = hitT(k) / tot;
                        missprob(k) = missT(k) / tot;
                    end
                end
                node = ref.getNodeByIndex(cacheNodes(cc));
                node.setResultHitProb(hitprob);
                node.setResultMissProb(missprob);
            end
end

function [detSojourn, dvals] = sojournConfig_(self, E)
% Sojourn option: 'stochastic' (default) averages each stage's transient over
% the random holding-time CDF; 'deterministic' (off by default) evaluates it at
% the fixed mean duration d_e = map_mean(holdTime{e}).
detSojourn = isfield(self.options,'sojourn') && strcmpi(self.options.sojourn,'deterministic');
dvals = zeros(1, E);
if detSojourn
    for e = 1:E
        dvals(e) = map_mean(self.envObj.holdTime{e});
    end
end
end

function v = detEval_(t, metric, d)
% Transient metric at the deterministic sojourn time d (clamped to the grid).
t = t(:); metric = metric(:);
d = max(t(1), min(d, t(end)));
v = interp1(t, metric, d, 'linear');
end
