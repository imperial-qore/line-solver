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
            logPn = logPn - factln(nir(r));
            if Q(ist,r) > 0
                logPn = logPn + nir(r)*log(Q(ist,r)/N(r));
            end
        end
    end
    Pnir = real(exp(logPn));
end
end