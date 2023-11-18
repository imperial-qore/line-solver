function model = JSIM2LINE(filename,modelName)
% MODEL = JSIM2LINE(FILENAME,MODELNAME)

% Copyright (c) 2012-2023, Imperial College London
% All rights reserved.
T0=tic;
% import model
Pref.Str2Num = 'always';
xDoc = xml_read(filename,Pref);
try
    xDoc = xDoc.sim;
end

% create network
if nargin<2
    [~,modelName] = fileparts(xDoc.ATTRIBUTE.name);
end
model = Network(modelName);

% create stations
node_name = cellfun(@(x) x.name, {xDoc.node.ATTRIBUTE},'UniformOutput',false)';
orig_node_name = node_name;
for i=1:length(node_name)
    node_name{i}=strrep(node_name{i},'/','_');
    node_name{i}=strrep(node_name{i},'\','_');
end

xsection = {xDoc.node.section};
strategy = cell(1,length(node_name));
xsection_par = {};
xsection_i = {};
xsection_javaClass = {};
sink_idx = -1;
source_idx = -1;

% This is to create the cs elements last, unclear if it affects correctness
% isStation = ones(1,length(node_name));
% for i=1:length(node_name)
%     xsection_i{i} = {xsection{i}};
%     xsection_i{i} = xsection_i{i}{1}; %   input, service, and output sections of node i
%     xsection_class{i} = {xsection_i{i}.ATTRIBUTE};
%     switch xsection_class{i}{1}.className % input section
%         case {'Buffer'}
%             xsection_i_type{i} = {xsection{i}.ATTRIBUTE};
%             switch xsection_i_type{i}{2}.className
%                 case {'StatelessClassSwitcher'}
%                     isStation(i) = 0;
%             end
%     end
% end

%for i=[find(isStation==1), find(isStation==0)]
for i=1:length(node_name)
    xsection_i{i} = {xsection{i}};
    xsection_i{i} = xsection_i{i}{1}; %   input, service, and output sections of node i
    xsection_javaClass{i} = {xsection_i{i}.ATTRIBUTE};
    switch xsection_javaClass{i}{1}.className % input section
        case 'JobSink'
            node{i} = Sink(model, node_name{i});
            sink_idx = i;
        case 'RandomSource'
            node{i} = Source(model, node_name{i});
            source_idx = i;
            xrouting{i} = {xsection_i{i}(3).parameter.subParameter.ATTRIBUTE};
            nSources=1;
        case 'Join'
            forkMap=find(cellfun(@any,strfind(cellfun(@class,model.nodes,'UniformOutput',false),'Fork')));
            if length(forkMap)>1
                line_error(mfilename,'JSIM2LINE supports at most a single fork-join pair.');
            end
            node{i} = Join(model, node_name{i}, node{forkMap});
            xrouting{i} = {xsection_i{i}(3).parameter.subParameter.ATTRIBUTE};
        case 'Queue'
            switch xsection_javaClass{i}{3}.className
                case 'Fork'
                    node{i} = Fork(model, node_name{i});
                    node{i}.setTasksPerLink(xsection_i{i}(3).parameter(1).value); %jobsPerLink
                    xrouting{i} = {xsection_i{i}(3).parameter(4).subParameter.ATTRIBUTE};
                otherwise
                    switch xsection_javaClass{i}{2}.className
                        case 'ServiceTunnel'
                            node{i} = Router(model, node_name{i});
                            xrouting{i} = {xsection_i{i}(3).parameter.subParameter.ATTRIBUTE};
                        otherwise
                            xsection_par{i} = {xsection{i}.parameter};
                            xsection_i_par{i} = xsection_i{i}.parameter;
                            
                            xsection_i_value{i} = {xsection_i_par{i}.value};
                            xsection_i_par_attr{i} = {xsection_i_par{i}.ATTRIBUTE};
                            
                            xsection_i_subpar{i} = {xsection_i_par{i}.subParameter};
                            %if xsection_i_value{i}{1}==-1
                            %    node{i} = Router(model, node_name{i});
                            %else
                            
                            xsvc{i} = {xsection_i{i}(2).parameter.subParameter};
                            xrouting{i} = {xsection_i{i}(3).parameter.subParameter.ATTRIBUTE};
                            
                            %    xget_strategy{i} = {xsection_i_par{i}.ATTRIBUTE};
                            %     switch xget_strategy{i}{3}.name
                            %         case 'LCFSstrategy'
                            %             strategy{i} = SchedStrategy.LCFS;
                            %         case 'FCFSstrategy'
                            %             strategy{i} = SchedStrategy.FCFS;
                            %     end
                            
                            xput_strategy{i} = xsection_i_par{i};
                            switch xput_strategy{i}(3).ATTRIBUTE.name
                                case 'retrialDistributions'
                                    % new XML format from 1.2.0
                                    %xretrial_strategy{i}= {xput_strategy{i}(4)};
                                    xput_strategy{i}= {xput_strategy{i}(5).subParameter.ATTRIBUTE};
                                otherwise
                                    xput_strategy{i}= {xput_strategy{i}(4).subParameter.ATTRIBUTE};
                            end
                            switch xput_strategy{i}{1}.name
                                case 'TailStrategy'
                                    strategy{i} = SchedStrategy.FCFS;
                                case 'TailStrategyPriority'
                                    strategy{i} = SchedStrategy.HOL;
                                case 'HeadStrategy'
                                    strategy{i} = SchedStrategy.LCFS;
                                case 'RandStrategy'
                                    strategy{i} = SchedStrategy.SIRO;
                                case 'SJFStrategy'
                                    strategy{i} = SchedStrategy.SJF;
                                case 'SEPTStrategy'
                                    strategy{i} = SchedStrategy.SEPT;
                                case 'LJFStrategy'
                                    strategy{i} = SchedStrategy.LJF;
                                case 'LEPTStrategy'
                                    strategy{i} = SchedStrategy.LEPT;
                            end
                            
                            xsection_i_type{i} = {xsection{i}.ATTRIBUTE};
                            switch xsection_i_type{i}{2}.className
                                case 'Delay'
                                    node{i} = DelayStation(model, node_name{i});
                                    xcapacity = {xsection_i_par{i}.value};
                                    node{i}.setCapacity(xcapacity{1}); % buffer size
                                case 'Server'
                                    node{i} = Queue(model, node_name{i}, strategy{i});
                                    xcapacity = {xsection_i_par{i}.value};
                                    node{i}.setCapacity(xcapacity{1}); % buffer size
                                    xsection_par_val{i} = {xsection_par{end}{2}.value};
                                    node{i}.setNumServers(xsection_par_val{i}{1});
                                    switch SchedStrategy.toId(strategy{i})
                                        case SchedStrategy.ID_SEPT
                                            schedparams{i} = NaN;
                                    end
                                case 'PSServer' % requires JMT >= 1.0.2
                                    strategy_i_sub={xsection_par{i}{2}.subParameter};
                                    strategy_i_sub4=strategy_i_sub{4}; strategy_i_sub4={strategy_i_sub4.ATTRIBUTE};
                                    strategy_i_sub5=strategy_i_sub{5};
                                    schedparams{i} = cell2mat({strategy_i_sub5.value});
                                    r=1; % we assume the strategies are identical across classes
                                    switch strategy_i_sub4{r}.name
                                        case 'EPSStrategy'
                                            strategy{i} = SchedStrategy.PS;
                                        case 'DPSStrategy'
                                            strategy{i} = SchedStrategy.DPS;
                                        case 'GPSStrategy'
                                            strategy{i} = SchedStrategy.GPS;
                                    end
                                    node{i} = Queue(model, node_name{i}, strategy{i});
                                    xcapacity = {xsection_i_par{i}.value};
                                    node{i}.setCapacity(xcapacity{1}); % buffer size
                                    xsection_par_val{i} = {xsection_par{end}{2}.value};
                                    node{i}.setNumServers(xsection_par_val{i}{1});
                                case 'ClassSwitch'
                                    strategy_i_sub={xsection_par{i}{2}.subParameter};
                                    strategy_i_sub1=strategy_i_sub{1}; strategy_i_sub1={strategy_i_sub1.subParameter};
                                    csMatrix = zeros(length(strategy_i_sub1));
                                    for r=1:length(strategy_i_sub1)
                                        csMatrix(r,:) = cell2mat({strategy_i_sub1{r}.value});
                                    end
                                    node{i} = ClassSwitch(model, node_name{i}, csMatrix);
                            end
                    end
            end
        case 'Storage'
            node{i} = Place(model, node_name{i});
        case 'Enabling'
            node{i} = Transition(model, node_name{i});
    end
end

% create classes
classes = {xDoc.userClass.ATTRIBUTE};
for r=1:length(classes)
    ref = findstring(node_name,classes{r}.referenceSource);
    switch classes{r}.type
        case 'closed'
            jobclass{r} =  ClosedClass(model, classes{r}.name, classes{r}.customers, node{ref}, classes{r}.priority);
        case 'open'
            % sink and source have been created before
            jobclass{r} =  OpenClass(model, classes{r}.name, classes{r}.priority);
            if strcmpi(classes{r}.referenceSource,'StatelessClassSwitcher')
                sourceIdx = cellisa(node,'Source');
                node{sourceIdx}.setArrival(jobclass{r},Disabled.getInstance());
            end
    end
end


for i=1:length(node_name)
    xsection_i{i} = {xsection{i}};
    xsection_i{i} = xsection_i{i}{1}; %   input, service, and output sections of node i
    xsection_javaClass{i} = {xsection_i{i}.ATTRIBUTE};
    switch xsection_javaClass{i}{1}.className % input section
        
        case 'Storage'
            node{i}.init();
            if xsection_i{1,i}(1).parameter(1).value == -1
                node{i}.setCapacity(Inf);
            else
                node{i}.setCapacity(xsection_i{1,i}(1).parameter(1).value);
            end
            if isa(xsection_i{1, i}(1).parameter(2).refClass,'cell')
                nclasses = length(xsection_i{1, i}(1).parameter(2).refClass);
            else
                nclasses = 1;
            end
            for c=1:nclasses
                if xsection_i{1, i}(1).parameter(2).subParameter(c).value == -1
                    node{i}.setClassCapacity(c, Inf);
                else
                    node{i}.setClassCapacity(c, xsection_i{1, i}(1).parameter(2).subParameter(c).value);
                end
                switch xsection_i{1, i}(1).parameter(3).subParameter(c).value
                    case 'BAS blocking'
                        node{i}.setDropRule(c, DropStrategy.BlockingAfterService);
                    case 'drop'
                        node{i}.setDropRule(c, DropStrategy.Drop);
                    case 'waiting queue'
                        node{i}.setDropRule(c, DropStrategy.WaitingQueue);
                end
            end
            node{i}.setState(0);
            
        case 'Enabling'
            % Enabling Section
            nmodes = length(xsection_i{1, i}(1).parameter(1).subParameter);
            for m=1:nmodes
                node{i}.setModeNames(m, xsection_i{1, i}(2).parameter(1).subParameter(m).value);
            end
            node{i}.init();
            for m=1:nmodes
                ninputs = length(xsection_i{1, i}(1).parameter(1).subParameter(m).subParameter.subParameter);
                for j=1:ninputs
                    refClasses = xsection_i{1, i}(1).parameter(1).subParameter(m).subParameter.subParameter(j).subParameter(2).refClass;
                    if isa(refClasses,'cell')
                        nclasses = length(refClasses);
                    else
                        nclasses = 1;
                    end
                    nodeName = xsection_i{1, i}(1).parameter(1).subParameter(m).subParameter.subParameter(j).subParameter(1).value;
                    ind = model.getNodeIndex(nodeName);
                    for k=1:nclasses
                        enable = xsection_i{1, i}(1).parameter(1).subParameter(m).subParameter.subParameter(j).subParameter(2).subParameter(k).value;
                        if enable == -1
                            node{i}.setEnablingConditions(m,k,ind,Inf);
                        else
                            node{i}.setEnablingConditions(m,k,ind,enable);
                        end
                        inhibit = xsection_i{1, i}(1).parameter(2).subParameter(m).subParameter.subParameter(j).subParameter(2).subParameter(k).value;
                        if inhibit == -1
                            node{i}.setInhibitingConditions(m,k,ind,Inf);
                        else
                            node{i}.setInhibitingConditions(m,k,ind,inhibit);
                        end
                    end
                end
            end
            % Timing Section
            for m=1:nmodes
                numOfServers = xsection_i{1, i}(2).parameter(2).subParameter(m).value;
                if numOfServers == -1
                    node{i}.setNumberOfServers(m, Inf);
                else
                    node{i}.setNumberOfServers(m, numOfServers);
                end
                timingSt = xsection_i{1, i}(2).parameter(3).subParameter(m).ATTRIBUTE.classPath;
                if strcmp(timingSt,'jmt.engine.NetStrategies.ServiceStrategies.ZeroServiceTimeStrategy')
                    node{i}.setTimingStrategy(m,TimingStrategy.IMMEDIATE);
                else
                    node{i}.setTimingStrategy(m,TimingStrategy.TIMED);
                    distribution = xsection_i{1, i}(2).parameter(3).subParameter(m).subParameter(1).ATTRIBUTE.name;
                    lambda = xsection_i{1, i}(2).parameter(3).subParameter(m).subParameter(2).subParameter(1).value;
                    switch distribution
                        case 'Exponential'
                            node{i}.setDistribution(m,Exp(lambda));
                        case 'Erlang'
                            lambda1 = xsection_i{1, i}(2).parameter(3).subParameter(m).subParameter(2).subParameter(2).value;
                            node{i}.setDistribution(m,Erlang(lambda, lambda1));
                        case 'Hyperexponential'
                            lambda1 = xsection_i{1, i}(2).parameter(3).subParameter(m).subParameter(2).subParameter(2).value;
                            lambda2 = xsection_i{1, i}(2).parameter(3).subParameter(m).subParameter(2).subParameter(3).value;
                            node{i}.setDistribution(m,HyperExp(lambda, lambda1, lambda2));
                        case 'Coxian'
                            lambda1 = xsection_i{1, i}(2).parameter(3).subParameter(m).subParameter(2).subParameter(2).value;
                            lambda2 = xsection_i{1, i}(2).parameter(3).subParameter(m).subParameter(2).subParameter(3).value;
                            node{i}.setDistribution(m,Coxian(lambda, lambda1, lambda2));
                        case 'Deterministic'
                            node{i}.setDistribution(m,Det(lambda));
                        case 'Pareto'
                            lambda1 = xsection_i{1, i}(2).parameter(3).subParameter(m).subParameter(2).subParameter(2).value;
                            node{i}.setDistribution(m,Pareto(lambda, lambda1));
                        case 'Gamma'
                            lambda1 = xsection_i{1, i}(2).parameter(3).subParameter(m).subParameter(2).subParameter(2).value;
                            node{i}.setDistribution(m,Gamma(lambda, lambda1));
                        case 'Uniform'
                            lambda1 = xsection_i{1, i}(2).parameter(3).subParameter(m).subParameter(2).subParameter(2).value;
                            node{i}.setDistribution(m,Uniform(lambda, lambda1));
                        case 'Replayer'
                            node{i}.setDistribution(m,Replayer(lambda));
                        case 'Weibull'
                            lambda1 = xsection_i{1, i}(2).parameter(3).subParameter(m).subParameter(2).subParameter(2).value;
                            node{i}.setDistribution(m,Weibull(lambda1, lambda)); % scale and shape are inverted in the constructor
                        case 'Lognormal'
                            lambda1 = xsection_i{1, i}(2).parameter(3).subParameter(m).subParameter(2).subParameter(2).value;
                            node{i}.setDistribution(m,Lognormal(lambda, lambda1));
                        otherwise
                            error('The model includes an arrival distribution not supported by the model-to-model transformation from JMT.')
                            node{i}.setDistribution(m,Exp(1));
                    end
                end
                firingPriorities = xsection_i{1, i}(2).parameter(4).subParameter(m).value;
                node{i}.setFiringPriorities(m, firingPriorities);
                firingWeights = xsection_i{1, i}(2).parameter(5).subParameter(m).value;
                node{i}.setFiringWeights(m, firingWeights);
            end
            % Firing Section
            for m=1:nmodes
                if isfield(xsection_i{1, i}(3).parameter(1).subParameter(m).subParameter,'CONTENT')
                    noutputs = 0;
                else
                    noutputs = length(xsection_i{1, i}(3).parameter(1).subParameter(m).subParameter.subParameter);
                end
                for j=1:noutputs
                    refClasses = xsection_i{1, i}(3).parameter(1).subParameter(m).subParameter.subParameter(j).subParameter(2).refClass;
                    if isa(refClasses, 'cell')
                        nclasses = length(refClasses);
                    else
                        nclasses = 1;
                    end
                    nodeName = xsection_i{1, i}(3).parameter(1).subParameter(m).subParameter.subParameter(j).subParameter(1).value;
                    ind = model.getNodeIndex(nodeName);
                    for k=1:nclasses
                        outcome = xsection_i{1, i}(3).parameter(1).subParameter(m).subParameter.subParameter(j).subParameter(2).subParameter(k).value;
                        if outcome == -1
                            node{i}.setFiringOutcome(k,m,ind,Inf);
                        else
                            node{i}.setFiringOutcome(k,m,ind,outcome);
                        end
                    end
                end
            end
    end
end

schedparams = cell(1,length(node_name));
% set service distributions
for i=1:length(node_name)
    if isa(node{i},'Source')
        for r=1:length(classes)
            xsection_par{i} = {xsection{i}.parameter};
            xsection_i_par{i} = xsection_i{i}.parameter;
            xsection_i_subpar{i} = {xsection_i_par{i}.subParameter};
            xarv_statdistrib{i}{r}={xsection_i_subpar{i}{1}.subParameter};
            if isempty(xarv_statdistrib{i}{r}{r})
                node{i}.setArrival(jobclass{r}, Disabled.getInstance());
            else
                xarv_statdistrib{i}{r}={xarv_statdistrib{i}{r}{r}.ATTRIBUTE};
                xarv{i} = {xsection_i{i}(1).parameter.subParameter};
                xarv_sec{i} = {xarv{i}{1}.subParameter};
                switch xarv_statdistrib{i}{r}{1}.name
                    case 'Exponential'
                        par={xarv_sec{i}{r}.subParameter}; par=par{2};
                        node{i}.setArrival(jobclass{r}, Exp(par.value));
                    case 'Erlang'
                        par={xarv_sec{i}{r}.subParameter}; par=par{2};
                        node{i}.setArrival(jobclass{r}, Erlang(par(1).value,par(2).value));
                    case 'Hyperexponential'
                        par={xarv_sec{i}{r}.subParameter}; par=par{2};
                        node{i}.setArrival(jobclass{r}, HyperExp(par(1).value,par(2).value,par(3).value));
                    case 'Coxian'
                        par={xarv_sec{i}{r}.subParameter}; par=par{2};
                        node{i}.setArrival(jobclass{r}, Coxian(par(1).value,par(2).value,par(3).value));
                    case 'Deterministic'
                        par={xarv_sec{i}{r}.subParameter}; par=par{2};
                        node{i}.setArrival(jobclass{r}, Det(par.value));
                    case 'Pareto'
                        par={xarv_sec{i}{r}.subParameter}; par=par{2};
                        node{i}.setArrival(jobclass{r}, Pareto(par(1).value, par(2).value));
                    case 'Weibull'
                        par={xarv_sec{i}{r}.subParameter}; par=par{2};
                        node{i}.setArrival(jobclass{r}, Weibull(par(1).value, par(2).value));
                    case 'Lognormal'
                        par={xarv_sec{i}{r}.subParameter}; par=par{2};
                        node{i}.setArrival(jobclass{r}, Lognormal(par(1).value, par(2).value));
                    case 'Gamma'
                        par={xarv_sec{i}{r}.subParameter}; par=par{2};
                        node{i}.setArrival(jobclass{r}, Gamma(par(1).value, par(2).value));
                    case 'Uniform'
                        par={xarv_sec{i}{r}.subParameter}; par=par{2};
                        node{i}.setArrival(jobclass{r}, Uniform(par(1).value, par(2).value));
                    case 'Replayer'
                        par={xarv_sec{i}{r}.subParameter}; par=par{2};
                        node{i}.setArrival(jobclass{r}, Replayer(par.value));
                    otherwise
                        line_error(mfilename,'The model includes an arrival distribution not supported by the model-to-model transformation from JMT.')
                        xarv_statdistrib{i}{r}{1}.name
                        node{i}.setArrival(jobclass{r}, Exp(1)); %TODO
                end
            end
        end
    elseif isa(node{i},'Queue') || isa(node{i},'DelayStation')
        if isempty(schedparams{i})
            switch SchedStrategy.toId(strategy{i})
                case {SchedStrategy.ID_SEPT,SchedStrategy.ID_LEPT}
                    schedparams{i} = NaN*ones(1,length(classes));
                otherwise
                    schedparams{i} = ones(1,length(classes));
            end
        end
        for r=1:length(classes)
            switch xsection_i_type{i}{2}.className
                case 'StatelessClassSwitcher'
                    % do nothing
                    continue
                case 'Delay'
                    xsvc_sec{i} = {xsvc{i}{1}.subParameter};
                otherwise
                    xsvc_sec{i} = {xsvc{i}{3}.subParameter};
            end
            if isempty(xsvc_sec{i}{r})
                xsvc_statdistrib{i}{r}={struct('name','Disabled')};
            else
                xsvc_statdistrib{i}{r}={xsvc_sec{i}{r}.ATTRIBUTE};
            end
            para_ir = schedparams{i}(r);
            switch xsvc_statdistrib{i}{r}{1}.name
                case 'Disabled'
                    node{i}.setService(jobclass{r}, Disabled.getInstance());
                case 'Replayer'
                    par={xsvc_sec{i}{r}.subParameter}; par=par{2};
                    node{i}.setService(jobclass{r}, Replayer(par.value), para_ir);
                case 'Exponential'
                    par={xsvc_sec{i}{r}.subParameter}; par=par{2};
                    node{i}.setService(jobclass{r}, Exp(par.value), para_ir);
                case 'Erlang'
                    par={xsvc_sec{i}{r}.subParameter}; par=par{2};
                    node{i}.setService(jobclass{r}, Erlang(par(1).value,par(2).value), para_ir);
                case 'Hyperexponential'
                    par={xsvc_sec{i}{r}.subParameter}; par=par{2};
                    node{i}.setService(jobclass{r}, HyperExp(par(1).value,par(2).value,par(3).value), para_ir);
                case 'Coxian'
                    par={xsvc_sec{i}{r}.subParameter}; par=par{2};
                    node{i}.setService(jobclass{r}, Coxian(par(1).value,par(2).value,par(3).value), para_ir);
                case 'Deterministic'
                    par={xsvc_sec{i}{r}.subParameter}; par=par{2};
                    node{i}.setService(jobclass{r}, Det(par.value));
                case 'Pareto'
                    par={xsvc_sec{i}{r}.subParameter}; par=par{2};
                    node{i}.setService(jobclass{r}, Pareto(par(1).value, par(2).value));
                case 'Weibull'
                    par={xsvc_sec{i}{r}.subParameter}; par=par{2};
                    node{i}.setService(jobclass{r}, Weibull(par(1).value, par(2).value));
                case 'Lognormal'
                    par={xsvc_sec{i}{r}.subParameter}; par=par{2};
                    node{i}.setService(jobclass{r}, Lognormal(par(1).value, par(2).value));
                case 'Gamma'
                    par={xsvc_sec{i}{r}.subParameter}; par=par{2};
                    node{i}.setService(jobclass{r}, Gamma(par(1).value, par(2).value));
                case 'Phase-Type'
                    par={xsvc_sec{i}{r}.subParameter}; par=par{2};
                    alpha = [par(1).subParameter.subParameter.value];
                    pars = {par(2).subParameter.subParameter};
                    T = [];
                    for c=1:length(pars)
                        T = [T; pars{c}.value];
                    end
                    if any(any(tril(T,-1))>0) % not APH
                        line_warning(mfilename,'The input model uses a general PH distribution, which is not yet supported in LINE. Fitting the first three moments into an APH distribution.\n');
                        PH = {T,-T*ones(size(T,1),1)*alpha};
                        ax = APH.fitCentral(map_mean(PH), map_var(PH), map_skew(PH));
                    else % APH
                        ax = APH(alpha, T);
                    end
                    node{i}.setService(jobclass{r}, ax);
                case 'Uniform'
                    par={xsvc_sec{i}{r}.subParameter}; par=par{2};
                    node{i}.setService(jobclass{r}, Uniform(par(1).value, par(2).value));
                otherwise
                    xsvc_statdistrib{i}{r}{1}.name
                    node{i}.setService(jobclass{r}, Exp(1), para_ir); %TODO
            end
        end
        for c=1:length(xsection_i_par_attr{i})
            switch xsection_i_par_attr{i}{c}.name
                case 'size'
                    node{i}.input.setSize(xsection_i_value{i}{c}); % buffer size
                case 'dropStrategies'
                    for r=1:length(node{i}.input.inputJobClasses)
                        switch xsection_i_subpar{i}{c}(r).value
                            case 'drop'
                                node{i}.input.inputJobClasses{r}{3} = DropStrategy.Drop; % this can be stored by JMT even if buffer capacity is infinite
                            case 'waiting queue'
                                node{i}.input.inputJobClasses{r}{3} = DropStrategy.WaitingQueue;
                            case 'BAS blocking'
                                node{i}.input.inputJobClasses{r}{3} = DropStrategy.BlockingAfterService;
                        end
                    end
            end
        end
    end
end

% create links
C = zeros(length(node_name)); % connection matrix
links = {xDoc.connection.ATTRIBUTE};
for l=1:length(links)
    source = findstring(orig_node_name,links{l}.source);
    target = findstring(orig_node_name,links{l}.target);
    %    model.addLink(station{source},station{target});
    C(source,target) = 1;
end

% assign routing probabilities
P = zeros(length(node_name)*length(classes));
for from=1:length(node_name)
    for target=1:length(node_name)
        if C(from,target)
            model.addLink(node{from},node{target});
        end
    end
end

for from=1:length(node_name)
    switch  class(node{from})
        case 'Place'
        case 'Transition'
        case 'Sink'
        otherwise
            for r=1:length(classes)
                switch xrouting{from}{r}.name
                    case 'Random'
                        node{from}.setRouting(jobclass{r},RoutingStrategy.RAND);
                        %                     targets = find(C(from,:));
                        %                     if isa(jobclass{r},'Class')
                        %                         targets = setdiff(targets, [sink_idx, source_idx]);
                        %                     end
                        %                     for target = targets(:)'
                        % %                        node{from}.setProbRouting(jobclass{r}, node{target}, 1 / length(targets));
                        %                        P((from-1)*length(classes)+r, (target-1)*length(classes)+r) = 1 / length(targets);
                        %                     end
                    case 'Probabilities'
                        node{from}.setRouting(jobclass{r},RoutingStrategy.PROB);
                        xroutprobarray = {xsection_i{from}(3).parameter.subParameter.subParameter};
                        xroutprob = {xroutprobarray{r}.subParameter}; xroutprob = xroutprob{1};
                        xroutprobdest = {xroutprob.subParameter};
                        for j=1:length(xroutprobdest)
                            xprob={xroutprobdest{j}.value};
                            target = findstring(node_name,xprob{1});
                            prob = xprob{2};
                            node{from}.setProbRouting(jobclass{r}, node{target}, prob);
                            %                        P((from-1)*length(classes)+r, (target-1)*length(classes)+r) = prob;
                        end
                    case 'Power of k'
                        line_error(mfilename,'Power of k import not supported yet.')
                    case 'Round Robin'
                        node{from}.setRouting(jobclass{r},RoutingStrategy.RROBIN);
                    case 'Weighted Round Robin'
                        node{from}.setRouting(jobclass{r},RoutingStrategy.WRROBIN);
                        xroutprobarray = {xsection_i{from}(3).parameter.subParameter.subParameter};
                        xroutprob = {xroutprobarray{r}.subParameter}; xroutprob = xroutprob{1};
                        xroutprobdest = {xroutprob.subParameter};
                        for j=1:length(xroutprobdest)
                            xprob={xroutprobdest{j}.value};
                            target = findstring(node_name,xprob{1});
                            weight = xprob{2};
                            node{from}.setRouting(RoutingStrategy.WRROBIN, node{target}, jobclass{r}, weight);
                        end
                    case 'Join the Shortest Queue (JSQ)'
                        node{from}.setRouting(jobclass{r},RoutingStrategy.JSQ);
                    case 'Disabled'
                        node{from}.setRouting(jobclass{r},RoutingStrategy.DISABLED);
                end
            end
    end
end
%model.link(P);
%line_printf(['JMT2LINE parsing time: ',num2str(Ttot),' s\n']);

if length(model.getIndexSourceStation)>1
    txt = sprintf('LINE supports JMT models with at most a single source node. You can refactor your JMT model in several ways:\n - If you are mapping in JMT each class to a different source, this is not required. You can instead assign the same reference station to each class and configure class routing in the routing panel of the source node.\n - In more general cases, you may follow these three steps:\n    (1) give a different name to each class of arrival, assigning these classes to a single source as reference station.\n    (2) put a class-switch node after the source to switch the new classes into the original classes they were in the model with multiple sources.\n    (3) configure the routing section of this class-switch node to set the same routing for the classes as they were in the original model.\n');
    error(txt);    
end


% Preload
state = zeros(length(node),length(classes));
if isfield(xDoc,'preload') && ~isempty(xDoc.preload.stationPopulations)
    npreloadStates = length(xDoc.preload.stationPopulations);
    for st=1:npreloadStates
        nodeName = xDoc.preload.stationPopulations(st).ATTRIBUTE.stationName;
        ind = model.getNodeIndex(nodeName);
        for r=1:length(xDoc.preload.stationPopulations(st).classPopulation)
            c = model.getClassIndex(xDoc.preload.stationPopulations(st).classPopulation(r).ATTRIBUTE.refClass);
            state(ind,c) = xDoc.preload.stationPopulations(st).classPopulation(r).ATTRIBUTE.population;
        end
        if isa(node{ind},'Place')
            %node{ind}.setState(state(ind,:));
            if length(classes)>1
                line_error(mfilename,'Import failed: Colored Petri net models are not yet supported in LINE.\n');
            end
        end
    end
end
try
    model.initFromMarginal(state);
catch
    line_warning(mfilename,'Import failed to automatically initialize the model.\n');
end
Ttot=toc(T0);
end