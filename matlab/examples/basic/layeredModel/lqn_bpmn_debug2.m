% More detailed debug to trace numerical explosion in lqn_bpmn
clear solver AvgTable
fprintf(1,'DEBUG: Tracing numerical explosion in lqn_bpmn model\n\n');

model = LayeredNetwork('myLayeredModel');

P{1} = Processor(model, 'R1_Processor', 100, SchedStrategy.FCFS);
P{2} = Processor(model, 'R2_Processor', Inf, SchedStrategy.INF);
P{3} = Processor(model, 'R3_Processor', 2, SchedStrategy.FCFS);
P{4} = Processor(model, 'R1A_Processor', 7, SchedStrategy.FCFS);
P{5} = Processor(model, 'R1B_Processor', 3, SchedStrategy.FCFS);
P{6} = Processor(model, 'R2A_Processor', 4, SchedStrategy.FCFS);
P{7} = Processor(model, 'R2B_Processor', 5, SchedStrategy.FCFS);

T{1} = Task(model, 'R1_Task', 100, SchedStrategy.REF).on(P{1}).setThinkTime(Exp.fitMean(20));
T{2} = Task(model, 'R2_Task', Inf, SchedStrategy.INF).on(P{2}).setThinkTime(Immediate());
T{3} = Task(model, 'R3_Task', 2, SchedStrategy.FCFS).on(P{3}).setThinkTime(Immediate());
T{4} = Task(model, 'R1A_Task', 7, SchedStrategy.FCFS).on(P{4}).setThinkTime(Immediate());
T{5} = Task(model, 'R1B_Task', 3, SchedStrategy.FCFS).on(P{5}).setThinkTime(Immediate());
T{6} = Task(model, 'R2A_Task', 4, SchedStrategy.FCFS).on(P{6}).setThinkTime(Immediate());
T{7} = Task(model, 'R2B_Task', 5, SchedStrategy.FCFS).on(P{7}).setThinkTime(Immediate());

E{1} = Entry(model, 'R1_Ref_Entry').on(T{1});
E{2} = Entry(model, 'R2_Synch_A2_Entry').on(T{2});
E{3} = Entry(model, 'R2_Synch_A5_Entry').on(T{2});
E{4} = Entry(model, 'R3_Synch_A9_Entry').on(T{3});
E{5} = Entry(model, 'R1A_Synch_A1_Entry').on(T{4});
E{6} = Entry(model, 'R1A_Synch_A2_Entry').on(T{4});
E{7} = Entry(model, 'R1A_Synch_A3_Entry').on(T{4});
E{8} = Entry(model, 'R1B_Synch_A4_Entry').on(T{5});
E{9} = Entry(model, 'R1B_Synch_A5_Entry').on(T{5});
E{10} = Entry(model, 'R1B_Synch_A6_Entry').on(T{5});
E{11} = Entry(model, 'R2A_Synch_A7_Entry').on(T{6});
E{12} = Entry(model, 'R2A_Synch_A8_Entry').on(T{6});
E{13} = Entry(model, 'R2A_Synch_A11_Entry').on(T{6});
E{14} = Entry(model, 'R2B_Synch_A9_Entry').on(T{7});
E{15} = Entry(model, 'R2B_Synch_A10_Entry').on(T{7});
E{16} = Entry(model, 'R2B_Synch_A12_Entry').on(T{7});

A{1} = Activity(model, 'A1_Empty', Immediate()).on(T{1}).boundTo(E{1}).synchCall(E{5},1);
A{2} = Activity(model, 'A2_Empty', Immediate()).on(T{1}).synchCall(E{6},1);
A{3} = Activity(model, 'A5_Empty', Immediate()).on(T{1}).synchCall(E{9},1);
A{4} = Activity(model, 'A6_Empty', Immediate()).on(T{1}).synchCall(E{10},1);
A{5} = Activity(model, 'A3_Empty', Immediate()).on(T{1}).synchCall(E{7},1);
A{6} = Activity(model, 'A4_Empty', Immediate()).on(T{1}).synchCall(E{8},1);
A{7} = Activity(model, 'E4_Empty', Immediate()).on(T{2}).boundTo(E{2});
A{8} = Activity(model, 'A7_Empty', Immediate()).on(T{2}).synchCall(E{11},1);
A{9} = Activity(model, 'A8_Empty', Immediate()).on(T{2}).synchCall(E{12},1);
A{10} = Activity(model, 'A9_Empty', Immediate()).on(T{2}).synchCall(E{14},1);
A{11} = Activity(model, 'A11_Empty', Immediate()).on(T{2}).synchCall(E{13},1).repliesTo(E{2});
A{12} = Activity(model, 'A12_Empty', Immediate()).on(T{2}).boundTo(E{3}).synchCall(E{16},1).repliesTo(E{3});
A{13} = Activity(model, 'A10_Empty', Immediate()).on(T{2}).synchCall(E{15},1);
A{14} = Activity(model, 'A13', Exp.fitMean(10)).on(T{3}).boundTo(E{4}).repliesTo(E{4});
A{15} = Activity(model, 'A1', Exp.fitMean(7)).on(T{4}).boundTo(E{5}).repliesTo(E{5});
A{16} = Activity(model, 'A2', Exp.fitMean(4)).on(T{4}).boundTo(E{6});
A{17} = Activity(model, 'A3', Exp.fitMean(5)).on(T{4}).boundTo(E{7}).repliesTo(E{7});
A{18} = Activity(model, 'A2_Res_Empty', Immediate()).on(T{4}).synchCall(E{2},1).repliesTo(E{6});
A{19} = Activity(model, 'A4', Exp.fitMean(8)).on(T{5}).boundTo(E{8}).repliesTo(E{8});
A{20} = Activity(model, 'A5', Exp.fitMean(4)).on(T{5}).boundTo(E{9});
A{21} = Activity(model, 'A6', Exp.fitMean(6)).on(T{5}).boundTo(E{10}).repliesTo(E{10});
A{22} = Activity(model, 'A5_Res_Empty', Immediate()).on(T{5}).synchCall(E{3},1).repliesTo(E{9});
A{23} = Activity(model, 'A7', Exp.fitMean(6)).on(T{6}).boundTo(E{11}).repliesTo(E{11});
A{24} = Activity(model, 'A8', Exp.fitMean(8)).on(T{6}).boundTo(E{12}).repliesTo(E{12});
A{25} = Activity(model, 'A11', Exp.fitMean(4)).on(T{6}).boundTo(E{13}).repliesTo(E{13});
A{26} = Activity(model, 'A9', Exp.fitMean(4)).on(T{7}).boundTo(E{14});
A{27} = Activity(model, 'A10', Exp.fitMean(6)).on(T{7}).boundTo(E{15}).repliesTo(E{15});
A{28} = Activity(model, 'A12', Exp.fitMean(8)).on(T{7}).boundTo(E{16}).repliesTo(E{16});
A{29} = Activity(model, 'A9_Res_Empty', Immediate()).on(T{7}).synchCall(E{4},1).repliesTo(E{14});

T{1}.addPrecedence(ActivityPrecedence.Serial(A{1}, A{2}));
T{1}.addPrecedence(ActivityPrecedence.Serial(A{3}, A{4}));
T{2}.addPrecedence(ActivityPrecedence.Serial(A{7}, A{8}));
T{2}.addPrecedence(ActivityPrecedence.Serial(A{10}, A{13}));
T{4}.addPrecedence(ActivityPrecedence.Serial(A{16}, A{18}));
T{5}.addPrecedence(ActivityPrecedence.Serial(A{20}, A{22}));
T{7}.addPrecedence(ActivityPrecedence.Serial(A{26}, A{29}));
T{1}.addPrecedence(ActivityPrecedence.OrFork(A{2},{A{5}, A{6}},[0.6,0.4]));
T{2}.addPrecedence(ActivityPrecedence.AndFork(A{8},{A{9}, A{10}}));
T{1}.addPrecedence(ActivityPrecedence.OrJoin({A{5}, A{6}}, A{3}));
T{2}.addPrecedence(ActivityPrecedence.AndJoin({A{9}, A{13}}, A{11}));

% Create solver
lnoptions = LN.defaultOptions;
lnoptions.verbose = 0;
lnoptions.iter_max = 10;  % Limit to see the explosion clearly

options = MVA.defaultOptions;
options.verbose = 0;

solver = LN(model, @(model) MVA(model, options), lnoptions);

% Manually iterate to trace the explosion
solver.init();

fprintf(1,'=== Starting manual iteration ===\n\n');

for it = 1:10
    fprintf(1,'\n========== ITERATION %d ==========\n', it);

    solver.pre(it);

    % Analyze each layer
    for e = 1:solver.nlayers
        [results{it,e}, ~] = solver.analyze(it, e);

        % Check for explosion in this layer
        if ~isempty(results{it,e})
            maxQN = max(results{it,e}.QN(:));
            maxTN = max(results{it,e}.TN(:));
            maxRN = max(results{it,e}.RN(:));
            if maxQN > 1e6 || maxRN > 1e6
                fprintf(1,'*** EXPLOSION DETECTED: Layer %d, Iter %d: max(QN)=%.6e, max(RN)=%.6e\n', ...
                    e, it, maxQN, maxRN);

                % Print layer details
                fprintf(1,'Layer %d structure:\n', e);
                fprintf(1,'  Number of jobs: %s\n', mat2str(solver.ensemble{e}.getNumberOfJobs));

                % Print service demands for each class
                for c = 1:solver.ensemble{e}.getNumberOfClasses
                    cls = solver.ensemble{e}.classes{c};
                    fprintf(1,'  Class %d (%s): ', c, cls.name);

                    % Get service time at server
                    serverIdx = solver.ensemble{e}.attribute.serverIdx;
                    node = solver.ensemble{e}.nodes{serverIdx};
                    if ~isempty(node.serviceProcess) && c <= length(node.serviceProcess)
                        proc = node.serviceProcess{c};
                        if ~isempty(proc) && ~isempty(proc{end})
                            fprintf(1,'servt=%.6e ', proc{end}.getMean());
                        end
                    end
                    fprintf(1,'\n');
                end
            end
        end
    end

    solver.results = results;
    solver.post(it);

    % Print key servt values after post
    fprintf(1,'\nAfter post(): key servt values:\n');
    lqn = solver.lqn;
    for eidx = (lqn.eshift+1):(lqn.eshift+min(5,lqn.nentries))
        fprintf(1,'  Entry %d servt=%.6e\n', eidx-lqn.eshift, solver.servt(eidx));
    end

    % Check for explosion in servt
    maxServt = max(solver.servt(:));
    if maxServt > 1e6
        fprintf(1,'*** servt EXPLOSION: max(servt)=%.6e\n', maxServt);
        [val, idx] = max(solver.servt(:));
        fprintf(1,'    at index %d\n', idx);
    end

    % Check convergence
    if solver.converged(it)
        fprintf(1,'Converged at iteration %d\n', it);
        break;
    end
end

fprintf(1,'\n=== Final servt values ===\n');
for eidx = (lqn.eshift+1):(lqn.eshift+lqn.nentries)
    fprintf(1,'Entry %d: servt=%.6e\n', eidx-lqn.eshift, solver.servt(eidx));
end
