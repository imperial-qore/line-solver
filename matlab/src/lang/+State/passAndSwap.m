function [cnew, depClass, chain] = passAndSwap(c, p, G)
% [CNEW, DEPCLASS, CHAIN] = PASSANDSWAP(C, P, G)
%
% Applies the pass-and-swap mechanism (Dorsman & Gardner 2024, Queueing
% Systems 107:205-256, Sect. 2.3) triggered by the service completion of the
% job in position P of an order-independent/pass-and-swap (PAS) queue.
%
% Inputs:
%   c - 1 x n ordered state vector of class indices; c(1) is the oldest job.
%   p - position (1..n) of the job whose service token completes.
%   G - (nclasses x nclasses) swapping graph adjacency (G(i,j)~=0 iff a class-i
%       job may take the place of a class-j job; undirected, self-loops allowed).
%
% Outputs:
%   cnew     - 1 x (n-1) ordered state vector after the pass-and-swap transition.
%   depClass - class index of the job that departs the system.
%   chain    - sequence of positions (in the original c) visited by the scan;
%              chain(1)=p is the vacated slot, c(chain(end)) is the departing job.
%
% Mechanism: starting at position p, the completing job scans backwards (towards
% newer/later positions) for the first job whose class is swappable with its own
% per G; it takes that job's place, ejecting it. The ejected job repeats the
% scan from its new position. The chain ends when an ejected job finds no
% swappable successor: that job departs. Classes shift one step along the chain
% and the head-of-chain slot is removed.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

n = numel(c);
if p < 1 || p > n
    line_error(mfilename, 'Completion position p is out of range for state c.');
end

chain = p;            % q0: position of the completing job (becomes vacated)
movingCls = c(p);     % class currently looking for a swappable successor
curpos = p;
while true
    q = 0;
    for j = curpos+1:n
        if G(movingCls, c(j))
            q = j;
            break;
        end
    end
    if q == 0
        break;        % movingCls finds no swappable successor -> it departs
    end
    chain(end+1) = q; %#ok<AGROW>
    movingCls = c(q);
    curpos = q;
end

depClass = c(chain(end));   % the last ejected job departs

% Shift classes one step along the chain: c(chain(i)) -> chain(i+1).
% The class previously at chain(end) is overwritten (it departed); the
% head-of-chain slot chain(1) is then removed.
cnew = c;
for i = 1:numel(chain)-1
    cnew(chain(i+1)) = c(chain(i));
end
cnew(chain(1)) = [];
end
