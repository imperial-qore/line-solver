function [nonfjmodel, fjclassmap, fjforkmap, fj_auxiliary_delays] = ht(model)
    % Transforms the queueing network containing a FJ subsystem into a queueing network without one.
    % Fork node replaced by routers
    % One artificial class is created for each parallel branch and for each class
    % Join node replaced by delay
    % We add another delay to model the sojourn time of the original classes.
    % ---
    % This approach is derived by PHILIP HEIDELBERGER and KISHOR S. TRIVEDI in
    % "Analytic Queueing Models for Programs with Internal Concurrency"
    sn = model.getStruct;
    fjclassmap = []; % s = fjclassmap(r) for auxiliary class r gives the index s of the original class
    fjforkmap = []; % f = fjforkmap(r) for auxiliary class r gives the associated fork node f
    nonfjmodel = model.copy();
    nonfjmodel.allowReplace = true;
    P = nonfjmodel.getLinkedRoutingMatrix;
    nonfjmodel.resetNetwork(true);
    nonfjmodel.resetStruct();
    Vnodes = cellsum(sn.nodevisits);
    forkedClasses = {};
    forkIndexes = find(sn.nodetype == NodeType.Fork)';
    fj_auxiliary_delays = {}; % d = fj_auxiliary_delays{j} for join j gives the additional delay d to mimic the sojourn time of the original classes

    %% Replace each fork with a router
    for f=forkIndexes
        nonfjmodel.nodes{f} = Router(nonfjmodel, nonfjmodel.nodes{f}.name);
        forkedClasses{f,1} = find(Vnodes(f,:)>0);
    end

    %% Replace each join with a delay
    for j=find(sn.nodetype == NodeType.Join)'
        nonfjmodel.nodes{j} = Delay(nonfjmodel, nonfjmodel.nodes{j}.name);
        nonfjmodel.stations{model.nodes{j}.stationIndex} = nonfjmodel.nodes{j};
        for c=1:length(nonfjmodel.classes)
            nonfjmodel.nodes{j}.setService(nonfjmodel.classes{c},Immediate());
        end

        %% Add another delay to mimic the sojourn time of the original classes for the artificial classes
        new_delay = Delay(nonfjmodel, ['Auxiliary Delay - ', nonfjmodel.nodes{j}.name]);
        fj_auxiliary_delays{j} = new_delay.index;
        for r=1:length(nonfjmodel.classes)
            new_delay.setService(nonfjmodel.classes{r}, Immediate());
        end
        for r=1:size(P, 1)
            for s=1:size(P, 2)
                P{r, s}(new_delay.index, :) = P{r,s}(j, :);
                P{r, s}(:, new_delay.index) = 0.0;
            end
            P{r, r}(j, :) = 0.0;
            P{r, r}(j, new_delay.index) = 1.0;
        end
    end

    % nonfjmodel.stations={nonfjmodel.stations{1:model.getNumberOfStations}}'; % remove automatically added station and put it where the join was
    nonfjmodel.connections = zeros(length(nonfjmodel.nodes));
    %% Create the auxiliary classes    
    for f=forkIndexes
        joinIdx = find(sn.fj(f,:));
        if length(joinIdx)>1
            line_error(mfilename,'SolverMVA supports at present only a single join station per fork node.');
        end
        forkedChains = find(sum(sn.chains(:,forkedClasses{f}),2));
        for fc=forkedChains'
            aux_classes = {}; % aux = aux_classes{r, par} gives the auxiliary class created for original class r and parallel branch par
            inchain = find(sn.chains(fc,:)); inchain = inchain(:)';
            for r=inchain
                if sn.nodevisits{fc}(f,r) == 0
                    continue;
                end
                parallel_branches = length(model.nodes{f}.output.outputStrategy{r}{3}); % Assumption: every class forks into exactly the same parallel branches
                for par=1:parallel_branches % One auxiliary class for each parallel branch
                    if model.nodes{f}.output.tasksPerLink > 1
                        line_error(mfilename, 'Multiple tasks per link are not supported in H-T.');
                    end
                    if isa(model.classes{r},'OpenClass')
                        line_error(mfilename, 'H-T method can be used only on closed models.');
                    end                    
                    aux_population = model.nodes{f}.output.tasksPerLink * model.classes{r}.population;
                    aux_classes{r, par} = ClosedClass(nonfjmodel, [nonfjmodel.classes{r}.name,'.',nonfjmodel.nodes{f}.name, '.B', int2str(par)], aux_population, nonfjmodel.nodes{fj_auxiliary_delays{joinIdx}}, 0);
                    fjclassmap(aux_classes{r, par}.index) = nonfjmodel.classes{r}.index;
                    fjforkmap(aux_classes{r, par}.index) = f;
                    % Set the service rates at the join node and at the stations
                    for i=1:sn.nnodes
                        if sn.isstation(i)
                            switch sn.nodetype(i)
                                case NodeType.Join
                                    nonfjmodel.nodes{i}.setService(aux_classes{r, par},Immediate());
                                case {NodeType.Source, NodeType.Fork}
                                    %no-op
                                otherwise
                                    nonfjmodel.nodes{i}.setService(aux_classes{r, par},model.nodes{i}.getService(model.classes{r}).copy());
                            end
                        end
                    end
                    nonfjmodel.nodes{fj_auxiliary_delays{joinIdx}}.setService(aux_classes{r, par}, Immediate());
                end
            end
            
            % Set the routing of the artificial classes
            for r=inchain
                if sn.nodevisits{fc}(f,r) == 0
                    continue;
                end
                parallel_branches = length(model.nodes{f}.output.outputStrategy{r}{3});
                for s=inchain
                    if sn.nodevisits{fc}(f,s) == 0
                        continue;
                    end
                    for par=1:parallel_branches
                        P{aux_classes{r, par}, aux_classes{s, par}} = P{r, s};
                        P{aux_classes{r, par}, aux_classes{s, par}}(f, :) = 0.0;
                        P{aux_classes{r, par}, aux_classes{s, par}}(f, model.nodes{f}.output.outputStrategy{r}{3}{par}{1}.index) = 1.0;
                        P{aux_classes{r, par}, aux_classes{s, par}}(joinIdx, fj_auxiliary_delays{joinIdx}) = 1.0;
                        P{aux_classes{r, par}, aux_classes{s, par}}(fj_auxiliary_delays{joinIdx}, :) = 0.0;
                        P{aux_classes{r, par}, aux_classes{s, par}}(fj_auxiliary_delays{joinIdx}, f) = 1.0;
                    end
                    % Route the original classes straight to the join to avoid the interference with the artificial classes
                    P{r,s}(f, :) = 0.0;
                    P{r,s}(f, joinIdx) = 1.0;
                end
            end
        end
    end
    nonfjmodel.relink(P);
    % snPrintRoutingMatrix(nonfjmodel.getStruct)
    % nonfjmodel.jsimgView
end