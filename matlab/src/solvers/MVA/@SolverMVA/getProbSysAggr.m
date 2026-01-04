function [Pnir,logPn] = getProbSysAggr(self)
% [PNIR,LOGPN] = GETPROBSYSSTATEAGGR()

if isempty(self.result)
    self.run;
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
                [~, nir, ~, ~] = State.toMarginal(sn, ist, state{ist});
                %                    logPn = logPn - log(sum(nir));
                for r=1:sn.nclasses
                    logPn = logPn - factln(nir(r));
                    if Q(ist,r)>0
                        logPn = logPn + nir(r)*log(Q(ist,r)/N(r));
                    end
                end
            end
            Pnir = real(exp(logPn));
    end
else
    line_error(mfilename,'getProbAggr not yet implemented for models with open classes.');
end
end