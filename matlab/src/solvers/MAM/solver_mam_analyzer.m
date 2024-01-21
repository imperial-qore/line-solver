function [QN,UN,RN,TN,CN,XN,runtime] = solver_mam_analyzer(sn, options)
% [QN,UN,RN,TN,CN,XN,RUNTIME] = SOLVER_MAM_ANALYZER(QN, OPTIONS)

% Copyright (c) 2012-2024, Imperial College London
% All rights reserved.

Tstart = tic;

if nargin<2 || isempty(options.config) || ~isfield(options.config,'merge')
    options.config = struct();
    options.config.merge = 'super';
    options.config.compress = 'mixture.order1';
    %options.config.compress = 'none';
    options.config.space_max = 128;
end

switch options.method
    case 'dec.mmap'
        % service distribution per class scaled by utilization used as
        % departure process
        [QN,UN,RN,TN,CN,XN] = solver_mam(sn, options);
    case {'default', 'dec.source'}
        % arrival process per chain rescaled by visits at each node
        [QN,UN,RN,TN,CN,XN] = solver_mam_basic(sn, options);
    case 'dec.poisson'
        % analyze the network with Poisson streams
        options.config.space_max = 1;
        [QN,UN,RN,TN,CN,XN] = solver_mam_basic(sn, options);
    case {'qna'}
        if snIsOpenModel(sn)
            [QN,UN,RN,TN,CN,XN] = solver_qna_mam(sn, options);
        elseif snIsClosedModel(sn)
            [QN,UN,RN,TN,CN,XN] = solver_qna_mam_closed(sn, options);
        else
            line_error(mfilename,'The qna method in SolverMAM does not support mixed models.');
        end
    otherwise
        line_error(mfilename,'Unknown method.');
end

for i=1:sn.nstations
    switch sn.schedid(i)
        case SchedStrategy.ID_EXT
            TN(i,:) = sn.rates(i,:);
    end
end

QN(isnan(QN))=0;
CN(isnan(CN))=0;
RN(isnan(RN))=0;
UN(isnan(UN))=0;
XN(isnan(XN))=0;
TN(isnan(TN))=0;

runtime = toc(Tstart);
end
