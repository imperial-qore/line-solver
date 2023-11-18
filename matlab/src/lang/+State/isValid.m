function isValid = isValid(sn, n, s, options)
% ISVALID = ISVALID(QN, N, S, OPTIONS)

% Copyright (c) 2012-2023, Imperial College London
% All rights reserved.

% n(r): number of jobs at the station in class r
% s(r): jobs of class r that are running

isValid = true;

if isa(sn,'Network')
    sn = sn.getStruct();
end

if isempty(n) & ~isempty(s)
    isValid = false;
    return
end

if iscell(n) %then n is a cell array of states
    ncell = n;
    n = [];
    for isf=1:length(ncell)
        ist = sn.statefulToStation(isf);
        ind = sn.statefulToNode(isf);
        [~, n(ist,:), s(ist,:), ~] = State.toMarginal(sn, ind, ncell{isf});
    end
end

R = sn.nclasses;
K = zeros(1,R);
for ist=1:sn.nstations
    for r=1:R
        K(r) = sn.phases(ist,r);
        if sn.nodetype(sn.stationToNode(ist)) ~= NodeType.Place
            if ~isempty(sn.proc) && ~isempty(sn.proc{ist}{r}) && any(any(isnan(sn.proc{ist}{r}{1}))) && n(ist,r)>0 % if disabled
                isValid = false;
                %            line_error(mfilename,sprintf('Chain %d is initialized with an incorrect number of jobs: %f instead of %d.', nc, statejobs_chain, njobs_chain));
                return
            end
        end
    end
    if any(n(ist,:)>sn.classcap(ist,:))
        line_warning(mfilename,'Station %d is in a state with more jobs than its allowed capacity.\n');
        isValid = false;
        return
    end
end

if nargin > 2 && ~isempty(s)
    for ist=1:sn.nstations
        if sn.nservers(ist)>0
            % if more running jobs than servers
            if sum(s(ist,:)) > sn.nservers(ist)
                switch sn.schedid(ist) % don't flag invalid if ps
                    case {SchedStrategy.ID_FCFS,SchedStrategy.ID_SIRO,SchedStrategy.ID_LCFS,SchedStrategy.ID_HOL}
                        isValid = false;
                        return
                end
            end
            % if more running jobs than jobs at the node
            if any(n<s)
                isValid = false;
                return
            end
            % non-idling condition
            if sum(s(ist,:)) ~= min(sum(n(ist,:)), sn.nservers(ist))
                % commented because in ps queues s are in service as well
                % isValid = false;
                %return
            end
        end
    end
end

for nc=1:sn.nchains
    njobs_chain = sum(sn.njobs(find(sn.chains(nc,:))));
    if ~isinf(njobs_chain)
        statejobs_chain = sum(sum(n(:,find(sn.chains(nc,:))),2),1);
        %if ~options.force && abs(1-njobs_chain/statejobs_chain) > options.iter_tol
        if abs(1-njobs_chain/statejobs_chain) > GlobalConstants.CoarseTol
            isValid = false;
            line_error(mfilename,sprintf('Chain %d is initialized with an incorrect number of jobs: %f instead of %d.', nc, statejobs_chain, njobs_chain));
            return
        end
        %end
    end
end
end
