function varargout = solver_tr_fjtag_analyzer(self, phase, varargin)
% VARARGOUT = SOLVER_TR_FJTAG_ANALYZER(SELF, PHASE, ...)
%
% Fork-join TAG AUGMENTATION as a solver-agnostic model transformation.
%
% ModelAdapter.fjtag rewrites the model into an exact tag-augmented struct in
% which every (fork, class, branch, tag) carries its own transient closed class
% and a join fires only once all siblings of one tag are buffered. The struct it
% returns is solved by the CALLER'S OWN analyzer, unchanged, and the auxiliary
% columns are folded back afterwards by SN_FJ_FOLDBACK.
%
% Both halves used to be written out inside @SolverCTMC/runAnalyzer.m and
% @SolverSSA/runAnalyzer.m, including the sn_orig/Korig bookkeeping and the
% per-sibling join response time, so a fix to one never reached the other. They
% live here once.
%
% This transformation has NO callback seam: the transformed struct is solved
% INLINE by the rest of the caller's analyzer rather than by a nested solver, so
% the caller invokes the two phases itself instead of handing control to a
% driver. That is the difference between it and @NetworkSolver/fjFixedPoint.m.
%
% Phase dispatch, following the SolverENV analyzer convention:
%
%   [snAug, ctx] = solver_tr_fjtag_analyzer(self, 'expand', sn, options)
%       SN is the untransformed struct. SNAUG is the tag-augmented struct to
%       solve; CTX carries what the lift needs and is opaque to the caller.
%
%   [QN,UN,RN,TN,AN,CN,XN] = solver_tr_fjtag_analyzer(self, 'lift', ctx, ...
%                                QN,UN,RN,TN,CN,XN, T)
%       Fold the auxiliary classes back onto the original ones and derive the
%       arrival rates in the ORIGINAL station/class coordinates. T is the
%       throughput handle matrix from GETAVGTPUTHANDLES.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

switch phase
    case 'expand'
        [varargout{1}, varargout{2}] = expand_(self, varargin{:});
    case 'lift'
        [varargout{1}, varargout{2}, varargout{3}, varargout{4}, ...
         varargout{5}, varargout{6}, varargout{7}] = lift_(self, varargin{:});
    otherwise
        line_error(mfilename, sprintf('Unknown fjtag-transformation phase: %s', phase));
end
end

function [snAug, ctx] = expand_(self, sn, options)
% Tag-augment the model and retain the original coordinates for the lift.
ctx = struct('sn_orig', sn, 'Korig', sn.nclasses, 'fjclassmap', []);
[~, fjsn, fjclassmap] = ModelAdapter.fjtag(self.model);
ctx.fjclassmap = fjclassmap;
snAug = fjsn;
line_debug(options, '%s: fork-join tag augmentation, %d classes (%d auxiliary), %d fork firings', ...
    self.getName(), snAug.nclasses, snAug.nclasses - ctx.Korig, length(snAug.fjsync));
end

function [QN,UN,RN,TN,AN,CN,XN] = lift_(self, ctx, QN,UN,RN,TN,CN,XN, T)
% Fold the auxiliary-class columns back, then re-derive the rates in the
% ORIGINAL coordinates. sn_pn_avg_rates must see sn_orig: a Place counts tokens
% and not firings, and rescaling has to happen before the arrival rates are
% derived so that everything downstream sees one convention.
[QN,UN,RN,TN,CN,XN] = sn_fj_foldback(QN,UN,RN,TN,CN,XN,ctx.fjclassmap,ctx.Korig);
[TN,~,RN] = sn_pn_avg_rates(ctx.sn_orig, QN, TN, [], RN);
AN = sn_get_arvr_from_tput(ctx.sn_orig, TN, T);
% Join stations report the per-sibling waiting time (JMT convention): queue
% length over the sibling arrival rate rather than over the join firing rate.
for ist=1:ctx.sn_orig.nstations
    if ctx.sn_orig.nodetype(ctx.sn_orig.stationToNode(ist)) == NodeType.Join
        for r=1:ctx.Korig
            if AN(ist,r) > 0
                RN(ist,r) = QN(ist,r)/AN(ist,r);
            end
        end
    end
end
self.result.fjclassmap = ctx.fjclassmap;
end
