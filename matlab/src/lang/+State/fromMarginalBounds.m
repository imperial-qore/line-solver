function space = fromMarginalBounds(sn, ind, lb, ub, cap, options)
% SPACE = FROMMARGINALBOUNDS(QN, IND, LB, UB, CAP, OPTIONS)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin<6 %~exist('options','var')
    options = Solver.defaultOptions;
end

% ind: node index
ist = sn.nodeToStation(ind);
%isf = sn.nodeToStateful(ind);

% returns all states lb<= x<= ub, where ub/lb are either a scalar (total
% number of jobs) or a vector (per-class number of jobs)
space = [];
if isempty(lb), lb=0*ub; end
R = sn.nclasses;
if length(lb) == 1, isVectorLB =0; else, isVectorLB = 1; end
if length(ub) == 1, isVectorUB =0; else, isVectorUB = 1; end
if isVectorLB~=isVectorUB, line_error(mfilename,'Bounds must either be both vectors or both scalars'); end

if isVectorUB && isVectorLB
    nmax = State.fromMarginal(sn, ind, ub, options);
    if isempty(nmax)        
         nmax = State.fromMarginal(sn, ind, ub, options);
    end
    n = pprodcon(lb,ub);
    while n ~= -1
        state = State.fromMarginal(sn, ind, n, options);
        space(end+1:end+size(state,1),(size(nmax,2)-size(state,2)+1):size(nmax,2)) = state;
        n = pprodcon(n,lb,ub);
    end
else % both scalar
    if ub >= lb
        for bi=ub:-1:lb % reverse order so that largest vector determines size(space,2)
            nset = multichoose(R,bi);
            for j=1:size(nset,1)
                state = State.fromMarginal(sn, ind, nset(j,:));
                if bi==ub && j==1
                    space(end+1:end+size(state,1),:) = state;
                else
                    space(end+1:end+size(state,1),(end-size(state,2)+1):end) = state;
                end
            end
        end
    end
end
space = unique(space,'rows');
% now we remove states that are not reachable
if sn.isstateful(ind)
    keep = [];
    for s=1:size(space,1)
        [ni,nir] = State.toMarginal(sn,ind,space(s,:));
        if sn.isstation(ind)
            if all(nir <= sn.classcap(ist,:)) && ni <= cap
                keep = [keep; s];
            end
        else
            if ni <= cap
                keep = [keep; s];
            end
        end
    end
    space = space(keep,:);
end
space = space(end:-1:1,:); % so that states with jobs in phase 1 comes earlier
end
