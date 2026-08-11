function [Pnir,runtime,fname] = solver_ctmc_marg(sn, options)
% [PNIR,RUNTIME,FNAME] = SOLVER_CTMC_MARG(QN, OPTIONS)
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.


M = sn.nstations;    %number of stations
K = sn.nclasses;    %number of classes
state = sn.state;
fname = '';

Tstart = tic;

[Q,SS,SSq,~,~,~,sn] = solver_ctmc(sn, options);
if options.keep
    fname = lineTempName;
    save([fname,'.mat'],'Q','SSq')
    line_printf('\nCTMC generator and state space saved in: ');
    line_printf([fname, '.mat'])
end
pi = ctmc_solve_reducible(Q);
pi(pi<GlobalConstants.Zero)=0;
pi = pi/sum(pi);

statesz = [];
for ind=1:sn.nnodes
    if sn.isstateful(ind)
        isf = sn.nodeToStateful(ind);
        statesz(isf) = size(sn.space{isf},2);
    end
end
cstatesz = [0,cumsum(statesz)];
Pnir = zeros(1,sn.nstations);
for ind=1:sn.nnodes
    if sn.isstateful(ind)
        isf = sn.nodeToStateful(ind);
        ist = sn.nodeToStation(ind);
        state_i = [zeros(1,size(sn.space{isf},2)-length(state{isf})),state{isf}];
        Pnir(ist) = sum(pi(findrows(SS(:,(cstatesz(isf)+1):(cstatesz(isf)+length(state_i))), state_i)));
    end
end

runtime = toc(Tstart);

%if options.verbose
%    line_printf('\nCTMC analysis completed. Runtime: %f seconds.\n',runtime);
%end
end
