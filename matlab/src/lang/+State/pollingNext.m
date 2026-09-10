function [q, mode, budget] = pollingNext(pinfo, pos, nbuf, R, arrived)
% [Q, MODE, BUDGET] = POLLINGNEXT(PINFO, POS, NBUF, R, ARRIVED)
%
% Resolve the tangible controller state a polling server reaches once it stops
% serving buffer POS. NBUF(r) is the number of class-r jobs waiting in the
% buffer (the service facility is empty at this point).
%
% ARRIVED (default false) selects where the cyclic walk starts. When false the
% server is LEAVING pos, so the walk starts at pos+1; when true the server has
% just ARRIVED at pos (a switchover into pos has completed, or the server is
% parked at pos) and pos itself is examined first, without charging its
% switchover a second time.
%
% MODE = 1  start a visit to buffer q: q has work and is reached in zero time.
%           BUDGET is the initial value of the ctr column for that visit.
% MODE = 2  enter the switchover into buffer q: a strictly positive timer, so
%           this is where the server dwells. BUDGET is 0 (the visit budget is
%           only set once the server arrives at q).
% MODE = 0  park at q: the walk completed a full lap without finding work and
%           without meeting a timed switchover, which can only happen when the
%           station is empty and every switchover is immediate. The server
%           would otherwise cycle in zero time forever, so it is held here
%           until the next arrival. q is canonical and unobservable.
%
% The walk is what makes an Immediate() switchover a zero-time leg rather than
% a state: the server passes straight through such a buffer when it has no
% work, and only stops at a buffer that either has work or costs time to reach.
% This is the classical cyclic-polling discipline, in which the server visits
% the buffers in strict cyclic order and pays the switchover of every buffer it
% moves to, whether or not that buffer turns out to have work.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 5 || isempty(arrived)
    arrived = false;
end

if arrived && pinfo.polled(pos) && nbuf(pos) > 0
    % the switchover into pos has already been paid, so a visit starts here
    q = pos; mode = 1; budget = State.pollingBudget(pinfo, nbuf(pos));
    return
end

p = pos;
for step = 1:R % a full lap, so that the last buffer examined is pos itself
    p = mod(p, R) + 1;
    if ~pinfo.polled(p)
        continue
    end
    if pinfo.hasSw(p)
        q = p; mode = 2; budget = 0;
        return
    end
    if nbuf(p) > 0
        q = p; mode = 1; budget = State.pollingBudget(pinfo, nbuf(p));
        return
    end
end

q = pos; mode = 0; budget = 0;
end
