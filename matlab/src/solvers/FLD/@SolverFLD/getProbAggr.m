function [Pnir,logPnir] = getProbAggr(self, ist)
% [PNIR,LOGPNIR] = GETPROBAGGR(IST)
%
% Probability of a SPECIFIC per-class job distribution at a station.
% Returns P(n1 jobs of class 1, n2 jobs of class 2, ...) for current state.
%
% Compare with getProbMarg: returns total queue-length distribution,
% i.e., P(n total jobs) summed over all class combinations.
%
% Input:
%   ist - Station index
%
% Output:
%   Pnir    - Scalar probability in [0,1]
%   logPnir - Log probability for numerical stability

if nargin<2 %~exist('ist','var')
    line_error(mfilename,'getProbAggr requires to pass a parameter the station of interest.');
end
sn = self.getStruct;
if ist > sn.nstations
    line_error(mfilename,'Station number exceeds the number of stations in the model.');
end
if isempty(self.result)
    self.run;
end
Q = self.result.Avg.Q;
N = sn.njobs;

state = sn.state{sn.stationToStateful(ist)};
[~, nir, ~, ~] = State.toMarginal(sn, ist, state);

openClasses = find(isinf(N));
closedClasses = find(isfinite(N));
logPnir = 0;

% Product-form probability for open classes
if ~isempty(openClasses)
    U = self.result.Avg.U;
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
        % where rho_r = U(ist,r) is per-server utilization
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
    logPnir = logPnir + nchoosekln(N(r),nir(r));
    logPnir = logPnir + nir(r)*log(Q(ist,r)/N(r));
    logPnir = logPnir + (N(r)-nir(r))*log(1-Q(ist,r)/N(r));
end

Pnir = real(exp(logPnir));
end