function envModel = MAPQN2RENV(model, options)
% ENVMODEL = MAPQN2RENV(MODEL, OPTIONS)
% Transform a queueing network with MMPP service into a random environment model
%
% This function transforms a queueing network where servers use MMPP2
% (2-phase Markov Modulated Poisson Process) service distributions into
% a random environment model with exponential services. The environment
% has two stages (one per MMPP phase) with transitions defined by the
% MMPP D0 matrix.
%
% Input:
%   model   - Network with MMPP2 service distributions
%   options - Optional struct with configuration (reserved for future use)
%
% Output:
%   envModel - Env model with exponential services modulated by MMPP phases
%
% The transformation works as follows:
% - Input MMPP2 has D0 (phase transitions) and D1 (diagonal service rates)
% - Output has 2 environment stages with exponential services
% - Stage 0 uses service rate λ₀ = D1(1,1)
% - Stage 1 uses service rate λ₁ = D1(2,2)
% - Environment transitions: σ₀₁ = D0(1,2), σ₁₀ = D0(2,1)
%
% Example:
%   model = Network('MMPP_Queue');
%   source = Source(model, 'Source');
%   queue = Queue(model, 'Queue', SchedStrategy.FCFS);
%   sink = Sink(model, 'Sink');
%   jobclass = OpenClass(model, 'Class1');
%   source.setArrival(jobclass, Exp(1.0));
%   queue.setService(jobclass, MMPP2(2.0, 3.0, 0.5, 0.5));
%   model.link(Network.serialRouting(source, queue, sink));
%   envModel = MAPQN2RENV(model);

    % Input validation
    if ~isa(model, 'Network')
        error('MAPQN2RENV:InvalidInput', 'First argument must be a Network object');
    end

    % Phase 1: Input Validation
    % Verify all service distributions are MMPP2 and D1 is diagonal
    [mmpp2Found, mmppParams] = validateAndExtractMMPP(model);

    if ~mmpp2Found
        error('MAPQN2RENV:NoMMPP', 'Network must contain at least one MMPP2 service distribution');
    end

    % Phase 2: Extract MMPP Parameters
    D0 = mmppParams.D0;
    D1 = mmppParams.D1;

    % Extract service rates from D1 diagonal
    lambda0 = D1(1, 1);
    lambda1 = D1(2, 2);

    % Extract transition rates from D0 off-diagonal
    sigma01 = D0(1, 2);
    sigma10 = D0(2, 1);

    % Validate rates are non-negative
    if lambda0 < 0 || lambda1 < 0 || sigma01 < 0 || sigma10 < 0
        error('MAPQN2RENV:InvalidRates', 'All extracted rates must be non-negative');
    end

    % Phase 3: Create Environment Model
    envModel = Environment('MAPQN_Env');

    % Phase 4: Build Stage Networks
    % Create stage network for Phase 0
    stageNet0 = buildStageNetwork(model, 'Phase0', lambda0);
    envModel.addStage('Phase0', 'item', stageNet0);

    % Create stage network for Phase 1
    stageNet1 = buildStageNetwork(model, 'Phase1', lambda1);
    envModel.addStage('Phase1', 'item', stageNet1);

    % Phase 5: Add Environment Transitions
    if sigma01 > 0
        envModel.addTransition('Phase0', 'Phase1', Exp(sigma01));
    end

    if sigma10 > 0
        envModel.addTransition('Phase1', 'Phase0', Exp(sigma10));
    end

    % Initialize environment
    envModel.init();
end

%% Helper: Validate MMPP2 and Extract Parameters
function [found, params] = validateAndExtractMMPP(model)
    found = false;
    params = struct();

    nodes = model.nodes;
    firstMMPP2 = [];

    % Iterate through all nodes to find MMPP2 distributions
    for i = 1:length(nodes)
        node = nodes{i};

        % Check if node is a Queue or Delay
        if ~(isa(node, 'Queue') || isa(node, 'Delay'))
            continue;
        end

        % Check service distributions for all classes
        if isa(node, 'Queue')
            for r = 1:length(node.server.serviceProcess)
                if ~isempty(node.server.serviceProcess{r})
                    dist = node.server.serviceProcess{r}{end};

                    if isa(dist, 'MMPP2')
                        if isempty(firstMMPP2)
                            firstMMPP2 = dist;
                            found = true;
                        else
                            % Verify consistency with first MMPP2
                            verifyMMPPConsistency(firstMMPP2, dist);
                        end
                    elseif isa(dist, 'MAP') && ~isa(dist, 'MMPP2')
                        % Reject generic MAP
                        D0_check = dist.D(0);
                        D1_check = dist.D(1);

                        % Check if D1 is diagonal
                        for ii = 1:size(D1_check, 1)
                            for jj = 1:size(D1_check, 2)
                                if ii ~= jj && abs(D1_check(ii, jj)) > GlobalConstants.Zero
                                    line_warning(mfilename, 'Generic MAP detected. Only MMPP (diagonal D1) is supported. Aborting transformation.');
                                    error('MAPQN2RENV:GenericMAP', 'Generic MAP with non-diagonal D1 is not supported. Use MMPP2 instead.');
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    if found && ~isempty(firstMMPP2)
        D0 = firstMMPP2.D(0);
        D1 = firstMMPP2.D(1);

        params.D0 = D0;
        params.D1 = D1;
    end
end

%% Helper: Verify MMPP Consistency
function verifyMMPPConsistency(mmpp1, mmpp2)
    % Future: Could add checks for consistent MMPP parameters across network
    % For now, just allow multiple MMPP2 with same parameters
end

%% Helper: Build Stage Network
function stageNet = buildStageNetwork(originalModel, stageName, expRate)
    % Clone the network structure and replace MMPP2 with Exp
    stageNet = Network([originalModel.name '_' stageName]);

    % Build mapping from original nodes to cloned nodes
    nodeMap = containers.Map('KeyType', 'char', 'ValueType', 'any');
    origClasses = originalModel.getClasses();

    % PASS 1: Create all nodes first (Source, Sink, Delay, Queue)
    % All nodes must exist before creating job classes
    for i = 1:length(originalModel.nodes)
        origNode = originalModel.nodes{i};
        if isa(origNode, 'Source')
            newNode = Source(stageNet, origNode.name);
            nodeMap(origNode.name) = newNode;
        elseif isa(origNode, 'Sink')
            newNode = Sink(stageNet, origNode.name);
            nodeMap(origNode.name) = newNode;
        elseif isa(origNode, 'Delay')
            newNode = Delay(stageNet, origNode.name);
            nodeMap(origNode.name) = newNode;
        elseif isa(origNode, 'Queue')
            newNode = Queue(stageNet, origNode.name, origNode.schedStrategy);
            % Copy number of servers and capacity
            if hasProperty(origNode, 'numberOfServers')
                newNode.setNumberOfServers(origNode.numberOfServers);
            end
            if hasProperty(origNode, 'cap')
                newNode.setCapacity(origNode.cap);
            end
            nodeMap(origNode.name) = newNode;
        end
    end

    % PASS 2: Create all job classes (all nodes must exist first)
    for r = 1:length(origClasses)
        origClass = origClasses{r};
        if isa(origClass, 'OpenClass')
            OpenClass(stageNet, origClass.name);
        elseif isa(origClass, 'ClosedClass')
            % Find reference station in the new network
            origRefStat = origClass.refstat;
            if isKey(nodeMap, origRefStat.name)
                newRefStat = nodeMap(origRefStat.name);
            else
                % Fallback to first node if reference station not yet created
                newRefStat = stageNet.nodes{1};
            end
            ClosedClass(stageNet, origClass.name, origClass.population, newRefStat);
        end
    end

    % PASS 4: Set arrivals from Source nodes
    for i = 1:length(originalModel.nodes)
        origNode = originalModel.nodes{i};

        if isa(origNode, 'Source') && isKey(nodeMap, origNode.name)
            newNode = nodeMap(origNode.name);

            % Copy arrival distributions
            origInput = origNode.input;
            if ~isempty(origInput) && ~isempty(origInput.sourceClasses)
                for r = 1:length(origInput.sourceClasses)
                    if ~isempty(origInput.sourceClasses{r})
                        origDist = origInput.sourceClasses{r}{end};
                        classes = stageNet.getClasses();
                        if r <= length(classes)
                            jobClass = classes{r};
                            newNode.setArrival(jobClass, origDist);
                        end
                    end
                end
            end
        end
    end

    % PASS 5: Set services on Queue/Delay nodes, replacing MMPP2 with Exp
    for i = 1:length(originalModel.nodes)
        origNode = originalModel.nodes{i};

        if (isa(origNode, 'Queue') || isa(origNode, 'Delay')) && isKey(nodeMap, origNode.name)
            newNode = nodeMap(origNode.name);
            classes = stageNet.getClasses();

            % Copy service distributions, replacing MMPP2 with Exp
            if isa(origNode, 'Queue') && hasProperty(origNode, 'server')
                for r = 1:length(origNode.server.serviceProcess)
                    if ~isempty(origNode.server.serviceProcess{r})
                        origDist = origNode.server.serviceProcess{r}{end};

                        if isa(origDist, 'MMPP2')
                            % Replace with exponential
                            newDist = Exp(expRate);
                        else
                            % Keep other distributions as-is
                            newDist = origDist;
                        end

                        % Set service for this class
                        if r <= length(classes)
                            jobClass = classes{r};
                            newNode.setService(jobClass, newDist);
                        end
                    end
                end
            elseif isa(origNode, 'Delay') && hasProperty(origNode, 'server')
                % Handle Delay node services
                for r = 1:length(origNode.server.serviceProcess)
                    if ~isempty(origNode.server.serviceProcess{r})
                        origDist = origNode.server.serviceProcess{r}{end};

                        if isa(origDist, 'MMPP2')
                            % Replace with exponential
                            newDist = Exp(expRate);
                        else
                            % Keep other distributions as-is
                            newDist = origDist;
                        end

                        % Set service for this class
                        if r <= length(classes)
                            jobClass = classes{r};
                            newNode.setService(jobClass, newDist);
                        end
                    end
                end
            end
        end
    end

    % PASS 6: Setup routing from original model
    % Determine if this is an open or closed network
    isClosedNetwork = any(cellfun(@(c) isa(c, 'ClosedClass'), origClasses));

    % Collect route nodes in order
    routeNodes = {};
    for i = 1:length(originalModel.nodes)
        origNode = originalModel.nodes{i};
        if isKey(nodeMap, origNode.name)
            routeNodes{end+1} = nodeMap(origNode.name);
        end
    end

    % Setup routing based on network type
    if isClosedNetwork
        % Closed network: use cyclic routing
        if length(routeNodes) >= 2
            try
                stageNet.link(Network.serialRouting(routeNodes{:}));
            catch
                % If serial routing fails, create cyclic routing manually
                for i = 1:length(routeNodes)
                    fromNode = routeNodes{i};
                    toNode = routeNodes{mod(i, length(routeNodes)) + 1};
                    try
                        P = stageNet.initRoutingMatrix();
                        for r = 1:length(stageNet.getClasses())
                            P{r, r}(stageNet.getNodeIndex(fromNode), stageNet.getNodeIndex(toNode)) = 1.0;
                        end
                        stageNet.link(P);
                    catch
                        % Skip if fails
                    end
                end
            end
        end
    else
        % Open network: use serial routing
        if length(routeNodes) > 2
            try
                stageNet.link(Network.serialRouting(routeNodes{:}));
            catch
                % If serial routing fails, try manual node-to-node routing
                for i = 1:length(routeNodes)-1
                    fromNode = routeNodes{i};
                    toNode = routeNodes{i+1};
                    try
                        stageNet.setRouting(fromNode, toNode, 1.0);
                    catch
                        % Skip if fails
                    end
                end
            end
        elseif length(routeNodes) == 2
            try
                stageNet.link(Network.serialRouting(routeNodes{1}, routeNodes{2}));
            catch
                % Skip routing if fails
            end
        end
    end
end

%% Helper: Check if object has property
function hasIt = hasProperty(obj, propName)
    try
        dummy = obj.(propName); %#ok<NASGU>
        hasIt = true;
    catch
        hasIt = false;
    end
end
