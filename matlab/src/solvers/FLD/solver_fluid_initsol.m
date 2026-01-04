function [init_sol, state] = solver_fluid_initsol(sn, options) %#ok<INUSD>
% [INIT_SOL, STATE] = SOLVER_FLUID_INITSOL(QN, OPTIONS) %#OK<INUSD>

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin<2 %~exist('options','var')
    options = Solver.defaultOptions; %#ok<NASGU>
end

init_sol = [];
for ind=1:sn.nnodes
    if sn.isstateful(ind)
        isf = sn.nodeToStateful(ind);
        ist = sn.nodeToStation(ind);
        % Skip non-station stateful nodes like Router
        if isnan(ist)
            continue
        end
        state_i = [];
        init_sol_i = [];  % compared to state_i, this does not track disabled classes and removes Inf entries in the Sources
        [~, nir, ~, kir_i] = State.toMarginal(sn, ind, sn.state{isf});
        switch sn.sched(ist)
            case {SchedStrategy.EXT}
                state_i(:,1) = Inf; % fluid does not model infinite buffer?
                for r=1:size(kir_i,2)
                    for k=1:length(sn.mu{ist}{r})
                        state_i(:,end+1) = kir_i(:,r,k);
                        if ~isnan(sn.rates(ist,r))
                            init_sol_i(:,end+1) = kir_i(:,r,k);
                        end
                    end
                end
            case {SchedStrategy.FCFS, SchedStrategy.SIRO, SchedStrategy.PS, SchedStrategy.INF, SchedStrategy.DPS, SchedStrategy.HOL}
                for r=1:size(kir_i,2)
                    for k=1:length(sn.mu{ist}{r})
                        if k==1
                            state_i(:,end+1) = nir(:,r) - sum(kir_i(:,r,2:end),3); % jobs in waiting buffer are re-started phase 1
                            if ~isnan(sn.rates(ist,r))
                                init_sol_i(:,end+1) = nir(:,r) - sum(kir_i(:,r,2:end),3); % jobs in waiting buffer are re-started phase 1
                            end
                        else
                            state_i(:,end+1) = kir_i(:,r,k);
                            if ~isnan(sn.rates(ist,r))
                                init_sol_i(:,end+1) = kir_i(:,r,k);
                            end
                        end
                    end
                end
            otherwise
                line_error(mfilename,sprintf('Unsupported scheduling policy at station %d',ist));
                return
        end
        init_sol = [init_sol, init_sol_i];
        state{isf} = state_i;
    end
end
end
