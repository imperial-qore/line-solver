function [runtime, sruntime, results] = iterate(self, options)
% ITERATE Perform LN solver iteration with optional MOL hierarchical approach
%
% [RUNTIME, SRUNTIME, RESULTS] = ITERATE(OPTIONS)
%
% For 'mol' method: Uses LQNS-style Method of Layers with nested iteration:
%   - Outer loop: Iterate on host (processor) layers until convergence
%   - Inner loop: Iterate on task layers until convergence (host demands fixed)
%
% For other methods ('default', 'moment3'): Uses flat iteration over all layers

if nargin < 2
    options = self.options;
end

line_debug('LN solver iterate starting: method=%s, nlayers=%d', self.options.method, self.nlayers);

if ~strcmp(self.options.method, 'mol')
    % Use parent implementation for default/moment3
    line_debug('Using parent EnsembleSolver iterate method');
    [runtime, sruntime, results] = iterate@EnsembleSolver(self, options);
    return;
end

line_debug('Using MOL (Method of Layers) iteration');

%% MOL (Method of Layers) implementation
T0 = tic;
E = getNumberOfModels(self);
results = cell(1, E);
sruntime = zeros(1, E);

init(self);

% Get MOL configuration parameters
mol_task_inner_max = self.options.config.mol_task_inner_max;
mol_task_inner_tol = self.options.config.mol_task_inner_tol;
mol_host_outer_tol = self.options.config.mol_host_outer_tol;
mol_min_steps = self.options.config.mol_min_steps;

% Initialize MOL iteration state
self.mol_it_host_outer = 0;
delta_host = Inf;
it = 0; % Total iteration counter for results indexing

%% First iteration: analyze ALL layers to initialize results
it = it + 1;
if self.options.verbose
    line_printf('\nMOL Initialization: Analyzing all layers');
end
self.pre(it);
sruntime(it, 1:E) = 0;

% Analyze all layers (host + task) in first iteration
for e = 1:E
    [results{it, e}, solverTime] = self.analyze(it, e);
    sruntime(it, e) = solverTime;
end
self.results = results;

% Full post-update after initial pass
self.post(it);

if self.options.verbose
    line_printf(' Done.\n');
end

%% OUTER LOOP: Host (processor) layer iteration
while (self.mol_it_host_outer < mol_min_steps || delta_host > mol_host_outer_tol) ...
      && self.mol_it_host_outer < self.options.iter_max

    self.mol_it_task_inner = 0;
    delta_task = Inf;

    % INNER LOOP: Task layer iteration (host layer demands kept FIXED)
    while delta_task > mol_task_inner_tol && self.mol_it_task_inner < mol_task_inner_max
        self.mol_it_task_inner = self.mol_it_task_inner + 1;
        it = it + 1;

        if self.options.verbose
            line_printf('\nMOL Outer %d, Inner %d: ', self.mol_it_host_outer + 1, self.mol_it_task_inner);
        end

        self.pre(it);
        sruntime(it, 1:E) = 0;

        % Analyze TASK layers only (host demands remain fixed)
        for e = self.taskLayerIndices
            [results{it, e}, solverTime] = self.analyze(it, e);
            sruntime(it, e) = solverTime;
        end

        % Copy host layer results from previous iteration (keeping them fixed)
        for e = self.hostLayerIndices
            results{it, e} = results{it-1, e};
        end

        self.results = results;

        % Compute task delta
        delta_task = self.computeTaskDelta();

        % Update metrics and TASK layer demands only (not host layer demands)
        self.updateMetrics(it);
        self.updateThinkTimes(it);
        if self.options.config.interlocking
            self.updatePopulations(it);
        end
        % updateLayers updates service times - only refresh task layer solvers
        self.updateLayers(it);
        self.updateRoutingProbabilities(it);

        % Reset only TASK layer solvers (keep host layers unchanged)
        for e = self.taskLayerIndices
            self.ensemble{e}.refreshChains();
            switch self.solvers{e}.name
                case {'SolverMVA', 'SolverNC'}
                    self.ensemble{e}.refreshRates();
                otherwise
                    self.ensemble{e}.refreshProcesses();
            end
            self.solvers{e}.reset();
        end

        if self.options.verbose
            line_printf('Task delta=%.6f', delta_task);
        end
    end

    self.mol_it_host_outer = self.mol_it_host_outer + 1;
    it = it + 1;

    if self.options.verbose
        line_printf('\nMOL Outer %d: Analyzing host layers... ', self.mol_it_host_outer);
    end

    self.pre(it);
    sruntime(it, 1:E) = 0;

    % Analyze HOST layers
    for e = self.hostLayerIndices
        [results{it, e}, solverTime] = self.analyze(it, e);
        sruntime(it, e) = solverTime;
    end

    % Copy task layer results from previous iteration
    for e = self.taskLayerIndices
        results{it, e} = results{it-1, e};
    end

    self.results = results;

    % Compute host delta
    delta_host = self.computeHostDelta();

    % Update metrics and HOST layer demands
    self.updateMetrics(it);
    self.updateThinkTimes(it);
    if self.options.config.interlocking
        self.updatePopulations(it);
    end
    self.updateLayers(it);
    self.updateRoutingProbabilities(it);

    % Reset HOST layer solvers (update their demands for next outer iteration)
    for e = self.hostLayerIndices
        self.ensemble{e}.refreshChains();
        switch self.solvers{e}.name
            case {'SolverMVA', 'SolverNC'}
                self.ensemble{e}.refreshRates();
            otherwise
                self.ensemble{e}.refreshProcesses();
        end
        self.solvers{e}.reset();
    end

    if self.options.verbose
        line_printf('Host delta=%.6f (%d inner iters)\n', ...
            delta_host, self.mol_it_task_inner);
    end
end

finish(self);
runtime = toc(T0);

if self.options.verbose
    line_printf('\nMOL Summary: %d outer iterations, total runtime: %.3fs\n', ...
        self.mol_it_host_outer, runtime);
end
end
