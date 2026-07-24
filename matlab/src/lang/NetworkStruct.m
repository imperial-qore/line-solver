function sn=NetworkStruct()
 % Data structure representation for a Network object
 %
 % Copyright (c) 2012-2026, Imperial College London
 % All rights reserved.
 
 sn=[]; %faster than sn=struct();
 sn.cap=[];     % total buffer size
 sn.cdscaling={}; % class-dependent (product-form) scalings beta_{i,r}(n): per-class output, argument the own-class marginal n_{i,r} or the total; see setClassDependence
 sn.cdscalingpeak=[]; % (nstations x nclasses) declared peak class-dependent rate scaling, for Util=T*S/peak normalization
 sn.jdscaling={}; % joint-dependent (non-product-form) scalings eta_i(n): scalar shared across classes or per-class, argument the joint vector (n_{i,1},...,n_{i,R}); see setJointDependence
 sn.jdscalingpeak=[]; % (nstations x nclasses) declared peak joint-dependent rate scaling, for Util=T*S/peak normalization
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
 sn.fjsync={}; % fork firing synchronizations (set only on FJ-augmented structs, see ModelAdapter.fjtag)
 sn.fjclassmap=[]; % (1,nclasses) original class of each FJ auxiliary class, 0 for originals (FJ-augmented structs only)
 sn.isfjaugmented=false; % true on FJ tag-augmented structs (Join/Fork carry count-vector states)
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
 sn.classspawn=[]; % (nclasses x 1) vector where entry r is the class injected at the same station on each completion of class r, -1 if none
 sn.signalremdist={}; % cell(nclasses,1) with removal distribution for each signal class (empty for single removal)
 sn.signalrempolicy=[]; % (nclasses x 1) with RemovalPolicy for each signal class
 sn.iscatastrophe=[]; % (nclasses x 1) where true indicates catastrophe signal
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
 sn.isbasblocking=[]; % (nnodes,1) 1 iff the node is the upstream/blocking side of a true-BAS relation (BUG-83)
 sn.isbasdestination=[]; % (nstations,nclasses) true iff refusing an arrival here must block an upstream BAS station
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
sn.regionmaxmem={};     % cell(F,1); regionmaxmem{f} is Matrix(M,1) with the region global memory budget replicated on member station rows, -1 = unbounded
 sn.regionlincon={};     % cell(F,2); regionlincon{f,1} is Matrix(C_f, K) linear constraint matrix and regionlincon{f,2} is Matrix(C_f, 1) capacity vector for region f
 sn.regionmembers={};    % cell(F,1); regionmembers{f} is logical(M,1), true where station i belongs to region f. Membership is NOT recoverable from region{f}: -1 there means unbounded, which is indistinguishable from not-a-member, so a region constrained only by regionlincon would read as empty

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

 % cache item state tracking (mirrors JAR varsparam)
 sn.varsparam=[];            % (nnodes x 1) item indices for cache state, -1 = none

 % marked (MMAP) source arrivals (populated by refreshStruct)
 sn.markidx=[];              % (nstations x nclasses) mark index (1-based) of class r
                             % at source station i, -1 = not a marked class

 % balking and retrial fields (populated by refreshStruct)
 sn.balkingStrategy=[];      % (M x K) BalkingStrategy id per (station,class), 0 = none
 sn.balkingThresholds={};    % (M x K) cell of {minJobs,maxJobs,prob} tuples per (station,class)
 sn.retrialType=[];          % (M x K) ProcessType id of retrial delay distribution, 0 = none
 sn.retrialMu=[];            % (M x K) retrial delay rate (1/mean)
 sn.retrialPhi=[];           % (M x K) retrial delay SCV
 sn.retrialProc={};          % (M x K) cell of {D0,D1} MAP/PH representation
 sn.retrialMaxAttempts=[];   % (M x K) max retrial attempts (-1 = unlimited)
 sn.retrialPolicy=[];        % (M x K) RetrialPolicy: LINEAR = per-customer rate n*nu, CONSTANT = orbit-wide rate nu
 sn.orbitMaxJobs=[];         % (M x K) orbit capacity (-1 = unbounded); a job finding the orbit full is lost
 sn.orbitImpatience={};      % (M x K) cell of {D0,D1} orbit abandonment representation
 % server breakdown / repair fields (populated by refreshStruct)
 sn.hasbreakdown=[];         % (nnodes,1) 1 iff the node's server is subject to breakdowns
 sn.breakdownMu=[];          % (M x 1) failure rate of an up server (0 = never fails)
 sn.repairMu=[];             % (M x 1) repair rate of a down server (0 = never repaired)
 sn.breakdownProc={};        % (M x 1) cell of {D0,D1} failure-time representation
 sn.repairProc={};           % (M x 1) cell of {D0,D1} repair-time representation
 sn.downServiceRates=[];     % (M x K) service rate while the server is down (0 = no service)

 % heterogeneous server fields are ragged, node-type-conditional parameters and
 % therefore live in the nodeparam container (indexed by node), not as flat root
 % fields. For a Queue node ind with server types, sn.nodeparam{ind} carries the
 % fields: nservertypes, servertypenames, serverspertype, servercompat,
 % heteroschedpolicy (see @MNetwork/refreshStruct.m).
end