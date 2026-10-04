function [Pnir,runtime,fname] = solver_ctmc_jointaggr(sn, options)
% [PNIR,RUNTIME,FNAME] = SOLVER_CTMC_JOINTAGGR(QN, OPTIONS)
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

M = sn.nstations;    %number of stations
K = sn.nclasses;    %number of classes
fname = '';
Tstart = tic;

[Q,SS,SSq,~,~,~,sn] = solver_ctmc(sn, options);
% SSq is an aggregate state space
if options.keep
    fname = lineTempName;
    save([fname,'.mat'],'Q','SSq')
    line_printf('\nCTMC generator and aggregate state space saved in: ');
    line_printf([fname, '.mat'])
end
% Through the seeded entry point, as solver_ctmc_marg does: re-solving unseeded
% made getProb* disagree with the means from getAvg on a reducible chain, the
% two answers coming from different recurrent classes of the same generator.
pi = ctmc_stationary(Q, SS, sn, options);
pi(pi<GlobalConstants.Zero)=0;

state = sn.state;
nvec = [];
for i=1:sn.nstations
    if sn.isstateful(i)
        isf = sn.stationToStateful(i);
        [~,nir,~,~] = State.toMarginal(sn, isf, state{isf});
        nvec = [nvec, nir(:)'];
    end
end
Pnir = sum(pi(findrows(SSq,nvec)));

runtime = toc(Tstart);

%if options.verbose
%    line_printf('\nCTMC analysis completed. Runtime: %f seconds.\n',runtime);
%end
end
