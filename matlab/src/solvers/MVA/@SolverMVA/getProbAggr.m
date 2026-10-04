function varargout = getProbAggr(self,varargin)
% [PNIR,LOGPNIR] = GETPROBAGGR(IST)
%
% Probability of a SPECIFIC per-class job distribution at a station.
% Returns P(n1 jobs of class 1, n2 jobs of class 2, ...) for current state.
%
% Compare with getProbMarg: returns queue-length distribution for a
% single class, i.e., P(n jobs of class r) for n=0,1,...,N(r).
%
% Input:
%   ist - Station index
%
% Output:
%   Pnir    - Scalar probability in [0,1]
%   logPnir - Log probability for numerical stability
% The result recorder captures the scalar this getter returned together
% with the solver that produced it -- see LineResultRecorder. Six of the
% statepr_* goldens hold exactly this number and nothing else, so recording
% it is what makes those goldens attributable instead of "the first bare
% number the example printed". The wrapper exists so that recording happens
% on EVERY exit path of the implementation below.
[scope, scopeGuard] = LineResultRecorder.enter(); %#ok<ASGLU>
[varargout{1:max(nargout,1)}] = getProbAggr_impl(self,varargin{:});
LineResultRecorder.captureScalar(scope, self, 'probAggr', varargout{1});
end

function [Pnir,logPnir] = getProbAggr_impl(self, ist)
% GETPROBAGGR_IMPL Implementation of GETPROBAGGR; see the wrapper above.


if nargin<2 %~exist('ist','var')
    line_error(mfilename,'getProbAggr requires to pass a parameter the station of interest.');
end
if isfield(self.options,'lang') && strcmp(self.options.lang,'cpp')
    p = CPPLINE.probAggr(self.name, self.model, self.options);
    Pnir = CPPLINE.probEntry(p, 'ProbAggr', ist, self.name, 'getProbAggr');
    logPnir = log(Pnir);
    return
end
if isempty(self.result)
    self.runAnalyzer;
end
% Read the struct AFTER the analysis: sn.state is empty until the model is
% initialized, and State.toMarginal indexes the phase block off its width.
sn = self.getStruct;
if ist > sn.nstations
    line_error(mfilename,'Station number exceeds the number of stations in the model.');
end
Q = self.result.Avg.Q;
N = sn.njobs;
if all(isfinite(N))
    switch self.options.method
        case 'exact'
            line_error(mfilename,'Exact marginal state probabilities not available yet in SolverMVA.');
        otherwise
            state = sn.state{sn.stationToStateful(ist)};
            [~, nir, ~, ~] = State.toMarginal(sn, ist, state);
            % Binomial approximation with mean fitted to queue-lengths.
            % Rainer Schmidt, "An approximate MVA ...", PEVA 29:245-254, 1997.
            logPnir = 0;
            for r=1:size(nir,2)
                % A CLASS WITH NO POPULATION CONTRIBUTES NOTHING: its binomial
                % factor is C(0,0)=1. Under class switching N(r)=0 can coexist
                % with Q(ist,r)>0, and log(Q/0) then makes the term 0*Inf = NaN.
                if N(r) == 0
                    continue
                end
                logPnir = logPnir + nchoosekln(N(r),nir(r));
                logPnir = logPnir + nir(r)*log(Q(ist,r)/N(r));
                logPnir = logPnir + (N(r)-nir(r))*log(1-Q(ist,r)/N(r));
            end
            Pnir = real(exp(logPnir));
    end
else
    % Mixed or open model: use product-form distribution
    U = self.result.Avg.U;
    state = sn.state{sn.stationToStateful(ist)};
    [~, nir, ~, ~] = State.toMarginal(sn, ist, state);
    openClasses = find(isinf(N));
    closedClasses = find(isfinite(N));
    logPnir = 0;

    % Product-form probability for open classes
    if ~isempty(openClasses)
        if sn.sched(ist) == SchedStrategy.INF
            % Delay (infinite server): independent Poisson per class
            for r = openClasses
                if Q(ist,r) > 0
                    logPnir = logPnir + nir(r)*log(Q(ist,r)) - Q(ist,r) - gammaln(nir(r)+1);
                elseif nir(r) > 0
                    logPnir = -Inf;
                end
            end
        elseif sn.sched(ist) ~= SchedStrategy.EXT
            % Queue station: multinomial-geometric product form
            % P(n_1,...,n_R) = (1-rho) * n!/prod(n_r!) * prod(rho_r^n_r)
            rho_total = sum(U(ist, openClasses));
            n_total = sum(nir(openClasses));
            if rho_total < 1
                logPnir = logPnir + log(1 - rho_total) + gammaln(n_total + 1);
                for r = openClasses
                    rho_r = U(ist, r);
                    if nir(r) > 0
                        if rho_r > 0
                            logPnir = logPnir + nir(r)*log(rho_r) - gammaln(nir(r)+1);
                        else
                            logPnir = -Inf;
                        end
                    end
                end
            else
                logPnir = -Inf;
            end
        end
    end

    % Binomial approximation for closed classes
    % Rainer Schmidt, "An approximate MVA ...", PEVA 29:245-254, 1997.
    for r = closedClasses
        % See the closed branch above: N(r)=0 must not reach log(Q/N).
        if N(r) == 0
            continue
        end
        logPnir = logPnir + nchoosekln(N(r),nir(r));
        logPnir = logPnir + nir(r)*log(Q(ist,r)/N(r));
        logPnir = logPnir + (N(r)-nir(r))*log(1-Q(ist,r)/N(r));
    end

    Pnir = real(exp(logPnir));
end
end
