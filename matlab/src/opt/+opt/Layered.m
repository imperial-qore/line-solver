classdef Layered
    % Layered  LayeredNetwork (LQN) support for line-opt. Static-method twin of
    % native-Python line_solver.opt.layered: model-type detection, element
    % resolution by name in a per-evaluation LQN model copy, the activity ->
    % processor mapping used to tag host-layer variables and key the per-layer
    % sensitivity table, and the SolverLN avg/sensitivity readers.
    %
    % IMPORTANT (see SolverLN.getSensitivityTable): the per-layer table holds
    % WITHIN-LAYER PARTIAL service-rate derivatives (fixed-point layer
    % parameters held constant); it omits cross-layer coupling and is thus a
    % biased estimate of the total derivative. lqnGradient='fd' finite-
    % differences the whole LayeredNetwork instead (correct total derivative);
    % 'partial_sens' uses this table directly; 'partial_plus_fd' corrects it
    % with a periodic full-model finite difference.

    methods (Static)
        function tf = isLayered(model)
            % True if MODEL is a LayeredNetwork (LQN), false for a flat Network.
            tf = ~isempty(model) && isa(model, 'LayeredNetwork');
        end

        function n = elemName(element)
            % Name of an LQN element (Processor/Task/Entry/Activity).
            if ischar(element)
                n = element;
            else
                n = element.getName();
            end
        end

        function m = distMean(value)
            % Mean of a think-time/demand that may be a distribution or scalar
            % ([] if unavailable). Used by LQN variables' currentValue.
            m = [];
            if isempty(value)
                return;
            end
            if isnumeric(value)
                m = double(value);
                return;
            end
            if isa(value, 'Distribution')
                try
                    m = value.getMean();
                catch
                    m = [];
                end
            end
        end

        function el = byName(elements, name)
            % First element of the cell/array ELEMENTS whose name is NAME ([]).
            el = [];
            for i = 1:numel(elements)
                if iscell(elements)
                    cand = elements{i};
                else
                    cand = elements(i);
                end
                if strcmp(opt.Layered.elemName(cand), name)
                    el = cand;
                    return;
                end
            end
        end

        function proc = resolveProcessor(model, name)
            % Resolve a Processor/Host by name inside a (copied) LQN model.
            proc = opt.Layered.byName(model.getHosts(), name);
        end

        function task = resolveTask(model, name)
            % Resolve a Task by name inside a (copied) LQN model.
            task = opt.Layered.byName(model.getTasks(), name);
        end

        function act = resolveActivity(model, name)
            % Resolve an Activity by name inside a (copied) LQN model.
            act = opt.Layered.byName(model.getActivities(), name);
        end

        function task = taskOfActivity(model, activityName)
            % The Task an activity belongs to, resolved in MODEL ([] if none).
            % Prefers the activity's own parent handle; falls back to the task
            % whose activity list or name matches.
            act = opt.Layered.resolveActivity(model, activityName);
            task = [];
            if isempty(act)
                return;
            end
            p = act.getParent();
            if ~isempty(p) && isa(p, 'Task')
                task = p;
                return;
            end
            % Fall back: activity.parentName is the owning task's name.
            pname = '';
            if isprop(act, 'parentName'), pname = act.parentName; end
            tasks = model.getTasks();
            for i = 1:numel(tasks)
                t = tasks{i};
                if ~isempty(pname) && strcmp(t.getName(), pname)
                    task = t; return;
                end
                acts = t.activities;
                for j = 1:numel(acts)
                    if strcmp(acts(j).getName(), activityName)
                        task = t; return;
                    end
                end
            end
        end

        function name = activityProcessorName(model, activityName)
            % Name of the processor an activity ultimately runs on ([] if the
            % Activity -> Task -> Processor chain is incomplete). Used both to
            % tag a HostDemand variable's host layer and to key its host-layer
            % sensitivity row (Layer=processor, Station=processor,
            % JobClass=activity).
            name = [];
            task = opt.Layered.taskOfActivity(model, activityName);
            if isempty(task)
                return;
            end
            proc = task.getParent();
            if isempty(proc)
                return;
            end
            name = proc.getName();
        end

        function name = taskProcessorName(model, taskName)
            % Name of the processor a task is deployed on ([] if undeployed).
            name = [];
            task = opt.Layered.resolveTask(model, taskName);
            if isempty(task)
                return;
            end
            proc = task.getParent();
            if ~isempty(proc)
                name = proc.getName();
            end
        end

        function tf = isRefTask(task)
            % True if TASK is a reference (workload-generating) task.
            s = task.getScheduling();
            if ~ischar(s)
                try
                    s = SchedStrategy.toText(s);
                catch
                    s = char(string(s));
                end
            end
            tf = strcmpi(s, 'ref') || strcmpi(s, 'reference');
        end

        function solver = makeSolver(model)
            % Construct a quiet SolverLN for a per-evaluation LQN model copy.
            solver = SolverLN(model, 'verbose', false);
        end

        function [solver, avgTable] = solveAvg(model)
            % Solve an LQN and return (solver, avgTable). The table has one row
            % per LQN node with columns Node, NodeType, QLen, Util, RespT,
            % ResidT, ArvR, Tput.
            solver = opt.Layered.makeSolver(model);
            avgTable = solver.getAvgTable();
        end

        function sens = computeSensitivities(solver)
            % Per-(Station,JobClass) within-layer service-rate partial
            % derivatives from SolverLN.getSensitivityTable, reshaped into a
            % containers.Map keyed 'Station||JobClass' -> struct with fields
            % Tput/RespT/QLen/Util (d(metric)/d(service rate)). [] on failure.
            sens = [];
            try
                SensTable = solver.getSensitivityTable();
            catch
                return;
            end
            if isempty(SensTable) || height(SensTable) == 0
                return;
            end
            sens = containers.Map('KeyType', 'char', 'ValueType', 'any');
            for r = 1:height(SensTable)
                st = SensTable.Station{r};
                cl = SensTable.JobClass{r};
                key = [st '||' cl];
                entry = struct('Tput', SensTable.dTput_dRate(r), ...
                    'RespT', SensTable.dRespT_dRate(r), ...
                    'QLen', SensTable.dQLen_dRate(r), ...
                    'Util', SensTable.dUtil_dRate(r));
                sens(key) = entry; %#ok<NASGU>
            end
        end
    end
end
