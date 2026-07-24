function [Pnir,runtime,fname] = solver_ctmc_joint(sn, options)
% [PNIR,RUNTIME,FNAME] = SOLVER_CTMC_JOINT(QN, OPTIONS)
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

M = sn.nstations;    %number of stations
K = sn.nclasses;    %number of classes
fname = '';
Tstart = tic;

[Q,SS,~,~,~,~,sn] = solver_ctmc(sn, options);
if options.keep
    fname = lineTempName;
    save([fname,'.mat'],'Q','SS')
    line_printf('\nCTMC generator and state space saved in: ');
    line_printf([fname, '.mat'])
end
pi = ctmc_solve_reducible(Q);
pi(pi<GlobalConstants.Zero)=0;

statevec = [];
state = sn.state;
for i=1:sn.nstations
    if sn.isstateful(i)
        isf = sn.nodeToStateful(i);
        state_i = [zeros(1,size(sn.space{isf},2)-length(state{isf})),state{isf}];
        statevec = [statevec, state_i];
    end
end
Pnir = pi(findrows(SS, statevec));

runtime = toc(Tstart);

%if options.verbose
%    line_printf('CTMC analysis completed. Runtime: %f seconds.\n',runtime);
%end
end
