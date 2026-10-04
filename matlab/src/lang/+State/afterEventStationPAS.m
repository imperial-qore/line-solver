function [outspace, outrate, outprob, eventCache, outstart, outpreempt] = afterEventStationPAS(sn, ind, ist, inspace, event, class, isSimulation, eventCache, R, V, key) %#ok<INUSL>
% [OUTSPACE, OUTRATE, OUTPROB, EVENTCACHE, OUTSTART, OUTPREEMPT] = AFTEREVENTSTATIONPAS(...)
%
% Event handler for pass-and-swap (PAS) / order-independent (OI) stations.
%
% State layout (see State.fromMarginal): the local state is the full ordered
% list of class indices c=(c1,...,cn), c1 the oldest job, stored left-aligned
% in the first W = sn.cap(ist) columns and right zero-padded; the trailing V
% columns hold routing local variables (carried through unchanged). There is
% no server/buffer split: service is governed by the total rate function mu(c)
% (sn.nodeparam{ind}.svcRateFun) and the swapping graph G
% (sn.nodeparam{ind}.swapGraph), per Dorsman & Gardner (2024), Sect. 2.
%
% - ARV (passive): a class-`class` job joins at the back, c -> (c, class),
%   subject to capacity (lost when full).
% - DEP (active): for each position p, the service token of position p fires at
%   rate Delta_mu(c1..cp) = mu(c1..cp) - mu(c1..c_{p-1}); the pass-and-swap
%   mechanism (State.passAndSwap) then determines the departing class. A DEP of
%   `class` collects all positions whose pass-and-swap ejects a class-`class`
%   job.
% - PHASE: none (PAS service is exponential).
%
% START/PREEMPT tags. A PAS station has one clock for the whole station and no
% servers to hold, so "in service" means "at a position whose rate increment
% Delta_mu is positive". A job therefore starts service exactly when a position
% goes from a zero increment to a positive one, which is what sub_startedHere
% compares. The two canaries bracket the rule: with mu(c) = 1 only the head is
% served, so a job starts when it reaches the head (M/M/1); with mu(c) = |c|
% every position is served, so every arrival starts at once (M/M/inf). Under a
% swap the tag follows the POSITION rather than the job identity, since
% pass-and-swap redefines which job holds a position. Nothing is ever pushed
% back out of service here, so OUTPREEMPT is identically zero.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

outspace = [];
outrate = [];
outprob = [];
outstart = [];
outpreempt = [];

muFun = sn.nodeparam{ind}.svcRateFun;
G = sn.nodeparam{ind}.swapGraph;
if isempty(muFun)
    line_error(mfilename,'PAS station has no service rate function mu(c); set it via setService(@(c) ...).');
end

W = size(inspace,2) - V;          % width of the ordered-list region
list = inspace(:,1:W);
varcols = inspace(:,(W+1):end);   % routing local variables
nrows = size(inspace,1);
cap = sn.cap(ist);

switch event
    case EventType.ARV % passive: append the arriving job at the back of the list
        for row=1:nrows
            c = list(row, list(row,:)>0);   % current ordered list (contiguous, left-aligned)
            n = numel(c);
            if n >= cap
                continue;                   % buffer full: arrival is lost
            end
            newc = [c, class];
            outspace = [outspace; newc, zeros(1, W-numel(newc)), varcols(row,:)]; %#ok<AGROW>
            outrate = [outrate; -1];        % passive action, rate unspecified %#ok<AGROW>
            outprob = [outprob; 1];         %#ok<AGROW>
            outstart = [outstart; sub_startedHere(muFun, c, newc, R)]; %#ok<AGROW>
        end
    case EventType.DEP % active: a class-`class` job departs via pass-and-swap
        for row=1:nrows
            c = list(row, list(row,:)>0);
            n = numel(c);
            if n == 0
                continue;
            end
            muPrev = 0;                     % mu of the empty prefix = 0
            for p=1:n
                muCur = muFun(c(1:p));
                ratep = muCur - muPrev;     % Delta_mu(c1..cp): rate of position p
                muPrev = muCur;
                if ratep <= 0
                    continue;               % position p receives no service
                end
                [cnew, depClass] = State.passAndSwap(c, p, G);
                if depClass ~= class
                    continue;               % a job of a different class departs
                end
                outspace = [outspace; cnew, zeros(1, W-numel(cnew)), varcols(row,:)]; %#ok<AGROW>
                outrate = [outrate; ratep]; %#ok<AGROW>
                outprob = [outprob; 1];     %#ok<AGROW>
                outstart = [outstart; sub_startedHere(muFun, c, cnew, R)]; %#ok<AGROW>
            end
        end
    case EventType.PHASE
        % PAS service is exponential: no intra-service phase transitions.
end

outstart = State.tagPad(outstart, size(outspace,1), R);
outpreempt = State.tagPad(outpreempt, size(outspace,1), R);

if isSimulation
    if ~isnan(key) && isobject(eventCache)
        eventCache{key} = {outprob, outspace, outrate, outstart, outpreempt};
    end
    if size(outspace,1) > 1
        if event == EventType.DEP
            tot_rate = sum(outrate);
            cum_rate = cumsum(outrate) / tot_rate;
            firing_ctr = 1 + max([0, find(rand > cum_rate')]);
            outspace = outspace(firing_ctr,:);
            outrate = sum(outrate);
            outprob = outprob(firing_ctr,:);
        else % passive ARV
            cum_prob = cumsum(outprob) / sum(outprob);
            firing_ctr = 1 + max([0, find(rand > cum_prob')]);
            outspace = outspace(firing_ctr,:);
            outrate = -1;
            outprob = 1;
        end
        outstart = outstart(firing_ctr,:); % the tags of the sampled arc
        outpreempt = outpreempt(firing_ctr,:);
    end
end
end

function st = sub_startedHere(muFun, cold, cnew, R)
% Count, per class, the positions of CNEW that are served (Delta_mu > 0) and
% were not served in COLD. mu of the empty prefix is 0, so position 1 is served
% whenever mu(c1) > 0.
st = zeros(1,R);
incNew = sub_increments(muFun, cnew);
incOld = sub_increments(muFun, cold);
for p = 1:numel(incNew)
    if incNew(p) <= 0
        continue
    end
    if p <= numel(incOld) && incOld(p) > 0
        continue % the position was already being served
    end
    st(cnew(p)) = st(cnew(p)) + 1;
end
end

function inc = sub_increments(muFun, c)
% Per-position service rate increments Delta_mu(c1..cp) = mu(c1..cp) - mu(c1..c_{p-1}).
inc = zeros(1, numel(c));
muPrev = 0;
for p = 1:numel(c)
    muCur = muFun(c(1:p));
    inc(p) = muCur - muPrev;
    muPrev = muCur;
end
end
