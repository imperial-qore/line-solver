function varargout = getProbSysAggr(self,varargin)
% [PNIR,LOGPN] = GETPROBSYSSTATEAGGR()
% The result recorder captures the scalar this getter returned together
% with the solver that produced it -- see LineResultRecorder. Six of the
% statepr_* goldens hold exactly this number and nothing else, so recording
% it is what makes those goldens attributable instead of "the first bare
% number the example printed". The wrapper exists so that recording happens
% on EVERY exit path of the implementation below.
[scope, scopeGuard] = LineResultRecorder.enter(); %#ok<ASGLU>
[varargout{1:max(nargout,1)}] = getProbSysAggr_impl(self,varargin{:});
LineResultRecorder.captureScalar(scope, self, 'probSysAggr', varargout{1});
end

function [Pnir,logPn] = getProbSysAggr_impl(self)
% GETPROBSYSAGGR_IMPL Implementation of GETPROBSYSAGGR; see the wrapper above.


if isfield(self.options,'lang') && strcmp(self.options.lang,'cpp')
    p = CPPLINE.probAggr(self.name, self.model, self.options);
    if isempty(p.ProbSysAggr)
        CPPLINE.cppUnsupported(self.name, 'getProbSysAggr', ...
            'line-cli''s -a prob answer for this solver carries no ProbSysAggr value');
    end
    Pnir = p.ProbSysAggr;
    logPn = log(Pnir);
    return
end

if isempty(self.result)
    self.runAnalyzer;
end
Q = self.result.Avg.Q;
sn = self.getStruct;
N = sn.njobs;
if all(isfinite(N))
    switch self.options.method
        case 'exact'
            line_error(mfilename,'Exact joint state probabilities not available yet in SolverMVA.');
        otherwise
            state = sn.state;
            % Binomial approximation with mean fitted to queue-lengths.
            % Rainer Schmidt, "An approximate MVA ...", PEVA 29:245-254, 1997.
            logPn = sum(factln(N));
            for ist=1:sn.nstations
                [~, nir, ~, ~] = State.toMarginal(sn, ist, state{sn.stationToStateful(ist)});
                %                    logPn = logPn - log(sum(nir));
                for r=1:sn.nclasses
                    % A CLASS WITH NO POPULATION CONTRIBUTES NOTHING. Under
                    % class switching a class can carry N(r)=0 and still show a
                    % positive Q(ist,r), since the jobs in it arrived by
                    % switching; log(Q/0) is then +Inf and nir(r) is 0, so the
                    % term evaluated to 0*Inf = NaN and poisoned the whole
                    % product. Its binomial factor is C(0,0)=1, i.e. zero in
                    % logs, which is what skipping it records.
                    if N(r) == 0
                        continue
                    end
                    logPn = logPn - factln(nir(r));
                    if Q(ist,r)>0
                        logPn = logPn + nir(r)*log(Q(ist,r)/N(r));
                    end
                end
            end
            Pnir = real(exp(logPn));
    end
else
    % Mixed or open model: product of per-station probabilities
    U = self.result.Avg.U;
    openClasses = find(isinf(N));
    closedClasses = find(isfinite(N));
    state = sn.state;
    logPn = 0;

    % Add closed-class multinomial normalization
    if ~isempty(closedClasses)
        logPn = sum(factln(N(closedClasses)));
    end

    for ist=1:sn.nstations
        [~, nir, ~, ~] = State.toMarginal(sn, ist, state{sn.stationToStateful(ist)});

        % Open classes: product-form contribution
        if ~isempty(openClasses)
            if sn.sched(ist) == SchedStrategy.INF
                % Delay: independent Poisson per class
                for r = openClasses
                    if Q(ist,r) > 0
                        logPn = logPn + nir(r)*log(Q(ist,r)) - Q(ist,r) - gammaln(nir(r)+1);
                    elseif nir(r) > 0
                        logPn = -Inf;
                    end
                end
            elseif sn.sched(ist) ~= SchedStrategy.EXT
                % Queue: multinomial-geometric product form
                rho_total = sum(U(ist, openClasses));
                n_total = sum(nir(openClasses));
                if rho_total < 1
                    logPn = logPn + log(1 - rho_total) + gammaln(n_total + 1);
                    for r = openClasses
                        rho_r = U(ist, r);
                        if nir(r) > 0
                            if rho_r > 0
                                logPn = logPn + nir(r)*log(rho_r) - gammaln(nir(r)+1);
                            else
                                logPn = -Inf;
                            end
                        end
                    end
                else
                    logPn = -Inf;
                end
            end
        end

        % Closed classes: binomial approximation
        for r = closedClasses
            % N(r)=0 with Q(ist,r)>0 is reachable under class switching; see
            % the closed branch above for why 0*Inf must not be formed.
            if N(r) == 0
                continue
            end
            logPn = logPn - factln(nir(r));
            if Q(ist,r) > 0
                logPn = logPn + nir(r)*log(Q(ist,r)/N(r));
            end
        end
    end
    Pnir = real(exp(logPn));
end
end