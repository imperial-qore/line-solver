function sn = PMIF2QN(filename,verbose)
% QN = PMIF2QN(FILENAME,VERBOSE)

% Copyright (c) 2012-2024, Imperial College London
% All rights reserved.

if nargin == 1; verbose = 0; end

parsed = PMIF.parseXML(filename,verbose);

if ~isempty(parsed)
    %% build CQN model
    % transformation of scheduling policies
    schedidTranslate = {  'IS',   SchedStrategy.ID_INF;
                        'FCFS', SchedStrategy.ID_FCFS;
                        'PS',   SchedStrategy.ID_PS};
    
    % extract basic information
    Ms = length(parsed.servers);
    Mw = length(parsed.workUnitServers);
    M = Ms + Mw;
    
    listServerIDs = cell(M,1);
    for i = 1:Ms
        listServerIDs{i} = parsed.servers(i).name;
    end
    for i = 1:Mw
        listServerIDs{Ms+i} = parsed.workUnitServers(i).name;
    end
    
    K = length(parsed.closedWorkloads);
    listClasses = cell(K,1);
    for i = 1:K
        listClasses{i} = parsed.closedWorkloads(i).name;
    end
    
    P = zeros(M*K,M*K);
    numservers = zeros(M,1);
    rates = zeros(M,K);
    schedid = nan(M,1);
    njobs = zeros(K,1);
    chains = eye(K);
    refstat = zeros(K,1);
    stationnames = listServerIDs;
    classnames = listClasses;
    
    % extract information from closedWorkloads
    for k = 1:K
        njobs(k) = parsed.closedWorkloads(k).numberJobs;
        myRefNode = findstring(stationnames, parsed.closedWorkloads(k).thinkDevice);
        refstat(k) =  myRefNode;
        rates(myRefNode, k) = 1/parsed.closedWorkloads(k).thinkTime;
        
        %routing
        for r = 1:size(parsed.closedWorkloads(k).transits,1)
            dest = findstring(stationnames,  parsed.closedWorkloads(k).transits{r,1} ) ;
            P((myRefNode-1)*K+k, (dest-1)*K+k) = parsed.closedWorkloads(k).transits{r,2};
        end
    end
    
    
    N = sum(njobs(isfinite(njobs)));
    
    % extract information from servers
    for i = 1:Ms
        schedid(i,1) = schedidTranslate{findstring(schedidTranslate(:,1),parsed.servers(i).scheduling), 2 };
        if schedid(i) == SchedStrategy.ID_INF
            numservers(i) = Inf;
        else
            numservers(i) = parsed.servers(i).quantity;
        end
    end
    
    for i = 1:Mw
        schedid(Ms+i,1) = schedidTranslate{findstring(schedidTranslate(:,1),parsed.workUnitServers(i).scheduling), 2 };
        if schedid(Ms+i) == SchedStrategy.ID_INF
            numservers(Ms+i) = Inf;
        else
            numservers(Ms+i) = parsed.workUnitServers(i).quantity;
        end
    end
    
    
    %extract information from demandServiceRequest
    for j = 1:length(parsed.demandServiceRequests)
        % service rate
        k = findstring(classnames, parsed.demandServiceRequests(j).workloadName );
        i = findstring(stationnames, parsed.demandServiceRequests(j).serverID );
        rates(i,k) = parsed.demandServiceRequests(j).serviceDemand / parsed.demandServiceRequests(j).numberVisits;
        
        %routing
        for r = 1:size(parsed.demandServiceRequests(j).transits,1)
            dest = findstring(stationnames,  parsed.demandServiceRequests(j).transits{r,1} ) ;
            P((i-1)*K+k, (dest-1)*K+k) = parsed.demandServiceRequests(j).transits{r,2};
        end
    end
    
    %extract information from workUnitServiceRequest
    for j = 1:length(parsed.workUnitServiceRequests)
        % service rate
        k = findstring(classnames, parsed.workUnitServiceRequests(j).workloadName );
        i = findstring(stationnames, parsed.workUnitServiceRequests(j).serverID );
        % work-unit-servers specify their own service time - indexed from Ms+1 to Ms + Mw in the list of servers
        rates(i,k) = 1 / ( parsed.workUnitServers(i-Ms).serviceTime * parsed.workUnitServiceRequests(j).numberVisits );
        
        %routing
        for r = 1:size(parsed.workUnitServiceRequests(j).transits,1)
            dest = findstring(stationnames,  parsed.workUnitServiceRequests(j).transits{r,1} ) ;
            P((i-1)*K+k, (dest-1)*K+k) = parsed.workUnitServiceRequests(j).transits{r,2};
        end
    end
    
    %extract information from timeServiceRequest
    for j = 1:length(parsed.timeServiceRequests)
        % service rate
        k = findstring(classnames, parsed.timeServiceRequests(j).workloadName );
        i = findstring(stationnames, parsed.timeServiceRequests(j).serverID );
        rates(i,k) = parsed.timeServiceRequests(j).serviceDemand / parsed.timeServiceRequests(j).numberVisits;
        
        %routing
        for r = 1:size(parsed.timeServiceRequests(j).transits,1)
            dest = findstring(stationnames,  parsed.timeServiceRequests(j).transits{r,1} ) ;
            P((i-1)*K+k, (dest-1)*K+k) = parsed.timeServiceRequests(j).transits{r,2};
        end
    end
    
    
    nodenames = stationnames;
    routing = RoutingStrategy.ID_PROB * ones(size(rates));
    for j=1:size(schedid,1)
        switch schedid(j)
            case SchedStrategy.ID_INF
                nodetype(j) = NodeType.Delay;
            otherwise
                nodetype(j) = NodeType.Queue;
        end
    end
            
    sn = NetworkStruct();
    sn.nnodes = numel(nodenames);
    sn.nclasses = length(classnames);
    sn.nclosedjobs = sum(njobs(isfinite(njobs)));
    sn.nservers = numservers;
    sn.nodetype = -1*ones(sn.nstations,1);
    sn.scv = ones(sn.nstations,sn.nclasses);
    %sn.forks = zeros(M,K);
    sn.njobs = njobs(:)';
    sn.refstat = refstat;
    sn.space = cell(sn.nstations,1);
    sn.dropid = -1* ones(sn.nstations,sn.nclasses);
    sn.routing = routing;
    sn.chains = [];
    sn.lst = {};
    sn.nodetype = nodetype;
    sn.isstation = NodeType.isStation(nodetype);
    sn.nstations = sum(sn.isstation);
    sn.isstateful = NodeType.isStateful(nodetype);
    sn.isstatedep = false(sn.nnodes,3); % col 1: buffer, col 2: srv, col 3: routing
    for ind=1:sn.nnodes
        switch sn.nodetype(ind)
            case NodeType.Cache
                sn.isstatedep(ind,2) = true; % state dependent service
                %                            sn.isstatedep(ind,3) = true; % state dependent routing
        end
        for r=1:sn.nclasses
            switch sn.routing(ind,r)
                case {RoutingStrategy.ID_RROBIN, RoutingStrategy.ID_WRROBIN, RoutingStrategy.ID_JSQ}
                    sn.isstatedep(ind,3) = true; % state dependent routing
            end
        end
    end
    sn.nstateful = sum(sn.isstateful);
    sn.state = cell(sn.nstations,1); for i=1:sn.nstateful sn.state{i} = []; end
    sn.nodenames = nodenames;
    sn.classnames = classnames;
    %sn.reindex();
    sn.schedid = schedid;
    sn.phi = ones(size(rates));
    sn.pie = ones(size(rates));
    sn.proc = cell(size(rates));
    sn.rt = P;
    for i=1:size(rates,1)
        sn.sched{i} = SchedStrategy.fromId(sn.schedid(i));
        for r=1:size(rates,2)
            if rates(i,r) == 0
                rates(i,r) = NaN;
                sn.proc{i}{r} = {[NaN],[NaN]};
            else
                sn.proc{i}{r} = map_exponential(1/rates(i,r));
            end
        end
    end
    sn.rates = rates;
else
    line_printf('Error: XML file empty or nonexistent');
    parsed = [];
end

