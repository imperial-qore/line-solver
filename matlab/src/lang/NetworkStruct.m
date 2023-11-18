function sn=NetworkStruct()
 % Data structure representation for a Network object
 %
 % Copyright (c) 2012-2023, Imperial College London
 % All rights reserved.
 
 sn=[]; %faster than sn=struct();
 sn.cap=[];     % total buffer size
 sn.cdscaling={}; % class-dependent scalings
 sn.chains=[];     % binary CxK matrix where 1 in entry (i,j) indicates that class j is in chain i.
 sn.classcap=[];    % buffer size for each class
 sn.classnames=string([]);  % name of each job class
 sn.classprio=[];       % scheduling priorities in each class (optional)
 sn.connmatrix=[]; % (i,j) entry if node i can route to node j
 sn.csmask=[]; % (r,s) entry if class r can switch into class s somewhere
 %forks;      % forks table from each station
 % (MKxMK matrix with integer entries), indexed first by
 % station, then by class
 sn.dropid=[]; % (i,r) gives the drop rule for class r at station i
 sn.fj=[]; % (i,j) is true if node j can join jobs forked by node i
 sn.inchain={}; % entry c is a vector with class ids in chain c
 sn.isstatedep=[]; % state dependent routing
 sn.isstation=[]; % element i is true if node i is a station
 sn.isstateful=[]; % element i is true if node i is stateful
 sn.isslc=[]; % element r is true if class r self-loops at its reference station
 sn.lldscaling={}; % limited load-dependent scalings
 sn.lst={}; % laplace-stieltjes transform
 sn.mu={};          % service rate in each service phase, for each job class in each station
 % (MxK cell with n_{i,k}x1 double entries)
 sn.nchains=[];           % number of chains (int)
 sn.nclasses=[];          % number of classes (int)
 sn.nclosedjobs=[];          % total population (int)
 sn.njobs=[];             % initial distribution of jobs in classes (Kx1 int)
 sn.nnodes=[]; % number of nodes (Mn int)
 sn.nservers=[];   % number of servers per station (Mx1 int)
 sn.nstations=[];  % number of stations (int)
 sn.nstateful=[];  % number of stations (int)
 sn.nvars=[]; % number of local variables
 sn.nodenames=string([]);   % name of each node
 sn.nodeparam={};     % parameters for local variables
 sn.nodetype=[]; % server type in each node
 sn.nodevisits={};  % visits placed by classes at the nodes
 sn.phases=[]; % number of phases in each service or arrival process
 sn.phasessz=[]; % number of phases
 sn.phaseshift=[]; % shift for phases
 sn.phi={};         % probability of service completion in each service phase,
 % for each job class in each station
 % (MxK cell with n_{i,k}x1 double entries)
 sn.pie={};        % probability of entry in each each service phase
 sn.proc={};     % cell matrix of service and arrival process representations
 sn.procid=[]; % service or arrival process type id
 sn.rates=[];       % service rate for each job class in each station
 sn.refstat=[];    % index of the reference node for each request class (Kx1 int)
 sn.routing=[];     % routing strategy type
 sn.rt=[];         % routing table with class switching
 % (M*K)x(M*K) matrix with double entries), indexed first by
 % station, then by class
 sn.rtorig={};         % linked routing table rtorig{r,s}(i,j)
 sn.rtnodes=[];         % routing table with class switching
 % (Mn*K)x(Mn*K) matrix with double entries), indexed first by
 % node, then by class
 sn.rtfun = @nan; % local routing functions
 % (Mn*K)x(Mn*K) matrix with double entries), indexed first by
 % station, then by class
 sn.sched=[];       % scheduling strategy in each station
 sn.schedid=[];       % scheduling strategy id in each station (optional)
 sn.schedparam=[];       % scheduling weights in each station and class (optional)
 sn.sync={};
 sn.space={};    % state space
 sn.state={};    % initial or current state
 sn.stateprior={};  % prior distribution of initial or current state
 sn.scv=[]; % squared coefficient of variation of service times (MxK)
 sn.visits={};           % visits placed by classes at the resources
 
 % hashing maps
 sn.nodeToStateful=[];
 sn.nodeToStation=[];
 sn.stationToNode=[];
 sn.stationToStateful=[];
 sn.statefulToNode=[];
end