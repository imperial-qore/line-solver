function sn=NetworkStruct()
 % Data structure representation for a Network object
 %
 % Copyright (c) 2012-2026, Imperial College London
 % All rights reserved.
 
 sn=[]; %faster than sn=struct();
 sn.cap=[];     % total buffer size
 sn.cdscaling={}; % class-dependent scalings
 sn.chains=[];     % binary CxK matrix where 1 in entry (i,j) indicates that class j is in chain i.
 sn.classcap=[];    % buffer size for each class
 sn.classnames=string([]);  % name of each job class
 sn.classprio=[];       % scheduling priorities in each class (optional)
 sn.classdeadline=[];   % deadline for each class (Inf = no deadline)
 sn.connmatrix=[]; % (i,j) entry if node i can route to node j
 sn.csmask=[]; % (r,s) entry if class r can switch into class s somewhere
 %forks;      % forks table from each station
 % (MKxMK matrix with integer entries), indexed first by
 % station, then by class
 sn.droprule=[]; % (i,r) gives the drop rule for class r at station i
 sn.fj=[]; % (i,j) is true if node j can join jobs forked by node i
 sn.gsync={};
 sn.inchain={}; % entry c is a vector with class ids in chain c
 sn.isstatedep=[]; % state dependent routing
 sn.isstation=[]; % element i is true if node i is a station
 sn.isstateful=[]; % element i is true if node i is stateful
 sn.isslc=[]; % element r is true if class r self-loops at its reference station
 sn.immfeed=[]; % (M x K) boolean matrix: immfeed(i,r) = true if class r uses immediate feedback at station i
 sn.issignal=[]; % element r is true if class r is a signal class (nclasses x 1)
 sn.signaltype={}; % signal type for each class, cell(nclasses,1) with NaN for non-signal classes
 sn.syncreply=[]; % (nclasses x 1) vector where entry r is reply signal class index for class r, -1 if no reply expected
 sn.signalRemovalDist={}; % cell(nclasses,1) with removal distribution for each signal class (empty for single removal)
 sn.signalRemovalPolicy=[]; % (nclasses x 1) with RemovalPolicy for each signal class
 sn.isCatastrophe=[]; % (nclasses x 1) where true indicates catastrophe signal
 sn.lldscaling={}; % limited load-dependent scalings
 sn.ljdscaling={}; % limited joint-dependent scalings (linearized per station)
 sn.ljdcutoffs=[]; % per-class cutoffs for joint dependence (M x K matrix)
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
 sn.schedparam=[];       % scheduling weights in each station and class (optional)
 sn.sync={};
 sn.space={};    % state space
 sn.state={};    % initial or current state
 sn.stateprior={};  % prior distribution of initial or current state
 sn.scv=[]; % squared coefficient of variation of service times (MxK)
 sn.visits={};           % visits placed by classes at the resources

 % finite capacity regions
 sn.nregions=[];         % number of finite capacity regions (F)
 sn.region={};           % cell array of size F; region{f} is Matrix(M, K+1) where entry (i,r) is max jobs of class r at station i in region f; (i,K+1) is global max at station i; -1 = infinite
 sn.regionrule=[];       % Matrix(F, K) where entry (f,r) is DropStrategy for class r in region f
 sn.regionweight=[];     % Matrix(F, K) where entry (f,r) is class weight for class r in region f (default 1.0)
 sn.regionsz=[];         % Matrix(F, K) where entry (f,r) is class size/memory for class r in region f (default 1)

 % hashing maps
 sn.nodeToStateful=[];
 sn.nodeToStation=[];
 sn.stationToNode=[];
 sn.stationToStateful=[];
 sn.statefulToStation=[];
 sn.statefulToNode=[];

 % reward definitions for CTMC reward computation
 sn.reward={};  % cell array of reward definitions
                % each entry is a struct with fields:
                %   .name - string identifier for the reward
                %   .fn   - function handle @(state, sn) -> double
                %   .type - reward type ('state' for state-dependent)

 % heterogeneous server fields
 sn.nservertypes=[];      % (M x 1) number of server types per station (0 = homogeneous)
 sn.servertypenames={};   % {station}{type} = name of server type
 sn.serverspertype={};    % {station}(type) = count of servers of that type
 sn.servercompat={};      % {station}(type, class) = 1 if type serves class, 0 otherwise
 sn.heterorates={};       % {station}{type, class} = service rate
 sn.heteroproc={};        % {station}{type, class} = {D0, D1} PH matrices
 sn.heteroprocid=[];      % (M x nTypes x K) or {station}(type, class) = ProcessType id
 sn.heteroschedpolicy=[]; % (M x 1) HeteroSchedPolicy per station
end