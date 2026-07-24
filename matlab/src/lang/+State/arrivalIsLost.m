function tf = arrivalIsLost(sn, ist, class)
% TF = ARRIVALISLOST(SN, IST, CLASS)
%
% True when an arrival of CLASS that finds no room at station IST is LOST, and
% false when it must BLOCK the upstream instead. This is the single predicate
% that decides the refusal semantics; every refusal path must branch on it.
%
% The rule is the CLASS TYPE, not the drop rule:
%
%   OPEN class   -> LOST. The external arrival stream is memoryless, so a job
%                   that finds the station full simply never enters. The caller
%                   must then leave the state UNCHANGED (a self-loop): the
%                   arrival event still fires, so the offered rate reaches the
%                   arrival-rate statistic and the loss shows up as
%                   ArvR - Tput. A self-loop cancels on the generator diagonal
%                   and therefore cannot perturb the stationary distribution,
%                   so QLen/Util/Tput are unaffected.
%   CLOSED class -> BLOCKED. A closed network's N jobs have nowhere to go;
%                   population conservation is a defining invariant, so a
%                   closed job can never be dropped. The caller must return an
%                   EMPTY outspace, which disables the upstream departure until
%                   room frees.
%
% An explicit blocking drop rule (BAS/BBS/RSRD) also asks for blocking, for any
% class: the user has said the job must wait rather than be lost.
%
% This is the same open/closed predicate as the CTMC analyzer's canDropClass
% and the BUG-12 utilization guard. Note the conventions are complementary, not
% contradictory: the arrival rate counts the OFFERED job (this predicate lets
% the event fire), while Util/QLen/Tput count only the CARRIED one.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if ~isempty(sn.droprule) && size(sn.droprule,1) >= ist && size(sn.droprule,2) >= class
    dr = sn.droprule(ist,class);
    if dr == DropStrategy.BAS || dr == DropStrategy.BBS || dr == DropStrategy.RSRD
        tf = false; % the user asked for blocking explicitly
        return
    end
end
% The rule above sees only THIS station's declaration. Under the upstream
% declaration form (cqn_bas_blocking.m) the BAS rule sits on the blocking
% station, not on the destination where the refusal happens, so the destination
% side is recorded separately by refreshLocalVars. Without this branch an open
% class refused here would be declared lost, the become-blocked edge would never
% fire, and the blocking station would behave as if its destination were
% unbounded. See BUG-83 and the R2 tandem-BAS fidelity fix.
if isfield(sn,'isbasdestination') && ~isempty(sn.isbasdestination) ...
        && size(sn.isbasdestination,1) >= ist && size(sn.isbasdestination,2) >= class ...
        && sn.isbasdestination(ist,class)
    tf = false; % refusing here must block the upstream BAS station
    return
end
% Open -> lost (self-loop), closed -> blocked. This decides only what happens
% once the arrival has already been refused (the state was not placed). The
% physical-vs-cutoff distinction is NOT made here: an open arrival always
% self-loops when refused, exactly as the pre-change producer did (a refused
% WAITQ arrival was left unchanged). What DID change with the physical cap is
% WHERE the refusal happens -- the hasRoom gate in afterEventStation stops
% placing the job at a physical cap so that it reaches this self-loop, whereas
% a cutoff-only class is still placed and truncated by the en_o filter. So the
% cutoff case never relies on this branch, and open self-looping here is the
% pre-change behaviour it must preserve.
tf = isinf(sn.njobs(class)); % open -> lost, closed -> blocked
end
