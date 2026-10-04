function result = getMAMResult(self)
% RESULT = GETMAMRESULT()
%
% Intermediate quantities of the matrix-analytic analysis of a single-queue
% model, in addition to the mean performance measures returned by getAvg.
%
% Mean values alone hide the objects the method is actually built on, so a
% matrix-analytic result cannot be inspected, taught, or checked against a
% published derivation. This accessor returns them.
%
% For a BMAP (or MAP) arrival stream feeding an exponential single server the
% result is that of qsys_bmapm1 and carries the M/G/1-type quantities: the
% phase-process stationary vectors theta and alpha, the randomized blocks A0,
% A1, B0 and Bk, the matrix G, the drift, the measured decay rate and the level
% probabilities.
%
% For a retrial station the result is that of qsys_bmapphnn_retrial and carries
% the orbit-level stationary distribution together with the truncation level and
% its residual.
%
% See also qsys_bmapm1, qsys_bmapphnn_retrial, SolverMAM.getAvg
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

sn = self.model.getStruct();

% A retrial station carries its own engine, whose result object already exposes
% the orbit-level internals.
[isRetrial, retInfo] = qsys_is_retrial(sn);
if isRetrial
    options = self.getOptions();
    options.verbose = false;
    [~,~,~,~,~,~,~,result] = solver_mam_retrial(sn, options);
    return
end

% Otherwise: BMAP/MAP arrivals into a single exponential server.
sourceIdx = [];
queueIdx = [];
for ist = 1:sn.nstations
    nodeIdx = sn.stationToNode(ist);
    if sn.nodetype(nodeIdx) == NodeType.Source
        sourceIdx = ist;
    elseif sn.nodetype(nodeIdx) == NodeType.Queue
        if isempty(queueIdx)
            queueIdx = ist;
        else
            line_error(mfilename, 'getMAMResult exposes the matrix-analytic internals of a single-queue model only.');
        end
    end
end
if isempty(sourceIdx) || isempty(queueIdx)
    line_error(mfilename, 'getMAMResult requires an open model with one Source and one Queue.');
end
if sn.nclasses > 1
    line_error(mfilename, 'getMAMResult exposes the matrix-analytic internals of a single-class model only.');
end
if sn.nservers(queueIdx) ~= 1
    line_error(mfilename, 'getMAMResult requires a single-server queue.');
end

arrivalProc = sn.proc{sourceIdx}{1};
if isempty(arrivalProc) || ~iscell(arrivalProc) || numel(arrivalProc) < 2
    line_error(mfilename, 'The arrival process has no Markovian (D0,D1,...) representation.');
end

serviceProc = sn.proc{queueIdx}{1};
if isempty(serviceProc) || ~iscell(serviceProc) || size(serviceProc{1},1) ~= 1
    line_error(mfilename, ['getMAMResult exposes the M/G/1-type internals for exponential service only; ' ...
        'the queue has a multi-phase service process.']);
end
mu = -serviceProc{1}(1,1);

result = qsys_bmapm1(arrivalProc, mu);
end
