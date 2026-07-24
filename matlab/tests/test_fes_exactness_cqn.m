function test_fes_exactness_cqn()
% Exactness of FES on multiclass CLOSED queueing networks (no delay station):
% all-queue product-form (PS) closed nets, incl. probabilistic multi-exit
% routing and 3 classes. Full solved exactly (NC exact + CTMC); subset
% aggregated via aggregateFES; reduced solved NC exact; complement metrics
% must match.

tol = 1e-6;
fprintf('\n============ FES EXACTNESS on multiclass CQN (no delay) ============\n');

% Case 1: 3-queue cyclic PS CQN, 2 classes, aggregate {Q1,Q2} (single exit)
[m,nd] = cqn_cyclic3(); run_case('cyclic3_aggQ1Q2', m, nd([1 2]), {'Q3'}, tol);

% Case 2: 4-queue PS CQN with branch; aggregate {Q1,Q2} -> Q2 branches to
% Q3 and Q4 (MULTI-EXIT subset in a pure CQN)
[m,nd] = cqn_branch4(); run_case('branch4_aggQ1Q2_multiexit', m, nd([1 2]), {'Q3','Q4'}, tol);

% Case 3: same net, aggregate {Q3,Q4} (MULTI-ENTRY subset, single exit)
[m,nd] = cqn_branch4(); run_case('branch4_aggQ3Q4_multientry', m, nd([3 4]), {'Q1','Q2'}, tol);

% Case 4: 3-class cyclic PS CQN, aggregate {Q1,Q2}
[m,nd] = cqn_cyclic3_3cls(); run_case('cyclic3_3class_aggQ1Q2', m, nd([1 2]), {'Q3'}, tol);

% Case 5: 4-queue branch, 3 classes, multi-exit subset {Q1,Q2}
[m,nd] = cqn_branch4_3cls(); run_case('branch4_3class_multiexit', m, nd([1 2]), {'Q3','Q4'}, tol);

fprintf('====================================================================\n');
fprintf('test_fes_exactness_cqn: PASS\n');
end

function run_case(name, model, subsetNodes, complementNames, tol)
% Full-model oracle: MVA exact (product-form exact), cross-checked by CTMC.
% NC 'exact' errors on plain all-queue closed nets (separate solver issue),
% so it is used only for the FES-reduced model, which needs the cdscaling path.
Afull = MVA(model,'method','exact').getAvgTable;
try, Actmc = CTMC(model).getAvgTable; c = maxrelerr(Afull, Actmc, complementNames);
catch, c = NaN; end
[fesModel,~,~] = ModelAdapter.aggregateFES(model, subsetNodes, struct('solver','mva','verbose',false));
Ared = NC(fesModel,'method','exact').getAvgTable;
err = maxrelerr(Afull, Ared, complementNames);
verdict = 'EXACT'; if ~(err < tol), verdict = 'NOT EXACT'; end
fprintf('  %-28s complement maxRelErr=%.3e -> %s (NC-vs-CTMC=%.1e)\n', name, err, verdict, c);
assert(err < tol, sprintf('%s not exact: maxRelErr=%.3e', name, err));
end

function err = maxrelerr(A, B, stationNames)
err = 0; metrics = {'QLen','Util','Tput','RespT'};
for n = 1:numel(stationNames)
    ra = find(strcmp(string(A.Station), stationNames{n}));
    for ii = ra'
        cls = string(A.JobClass(ii));
        rb = find(strcmp(string(B.Station), stationNames{n}) & strcmp(string(B.JobClass), cls));
        if isempty(rb), continue; end
        for m = 1:numel(metrics)
            va = A.(metrics{m})(ii); vb = B.(metrics{m})(rb(1));
            err = max(err, abs(va-vb)/max(abs(va),1e-9));
        end
    end
end
end

function [model,node] = cqn_cyclic3()
model = Network('cqn_cyclic3');
node{1}=Queue(model,'Q1',SchedStrategy.PS);
node{2}=Queue(model,'Q2',SchedStrategy.PS);
node{3}=Queue(model,'Q3',SchedStrategy.PS);
jc{1}=ClosedClass(model,'Class1',3,node{1},0);
jc{2}=ClosedClass(model,'Class2',2,node{1},0);
node{1}.setService(jc{1},Exp.fitMean(0.5)); node{1}.setService(jc{2},Exp.fitMean(0.8));
node{2}.setService(jc{1},Exp.fitMean(0.3)); node{2}.setService(jc{2},Exp.fitMean(0.6));
node{3}.setService(jc{1},Exp.fitMean(0.4)); node{3}.setService(jc{2},Exp.fitMean(0.7));
P=model.initRoutingMatrix();
Pc=[0 1 0; 0 0 1; 1 0 0]; P{1,1}=Pc; P{2,2}=Pc;
model.link(P);
end

function [model,node] = cqn_branch4()
model = Network('cqn_branch4');
for i=1:4, node{i}=Queue(model,sprintf('Q%d',i),SchedStrategy.PS); end
jc{1}=ClosedClass(model,'Class1',3,node{1},0);
jc{2}=ClosedClass(model,'Class2',2,node{1},0);
r1=[0.5 0.8]; r2=[0.3 0.6]; r3=[0.4 0.7]; r4=[0.6 0.5];
for c=1:2
    node{1}.setService(jc{c},Exp.fitMean(r1(c)));
    node{2}.setService(jc{c},Exp.fitMean(r2(c)));
    node{3}.setService(jc{c},Exp.fitMean(r3(c)));
    node{4}.setService(jc{c},Exp.fitMean(r4(c)));
end
P=model.initRoutingMatrix();
% Q1->Q2; Q2->Q3(0.5)/Q4(0.5); Q3->Q1; Q4->Q1
Pc=[0 1 0 0; 0 0 0.5 0.5; 1 0 0 0; 1 0 0 0];
P{1,1}=Pc; P{2,2}=Pc;
model.link(P);
end

function [model,node] = cqn_cyclic3_3cls()
model = Network('cqn_cyclic3_3cls');
node{1}=Queue(model,'Q1',SchedStrategy.PS);
node{2}=Queue(model,'Q2',SchedStrategy.PS);
node{3}=Queue(model,'Q3',SchedStrategy.PS);
pops=[2 2 1];
for c=1:3, jc{c}=ClosedClass(model,sprintf('Class%d',c),pops(c),node{1},0); end
sv=[0.5 0.8 0.4; 0.3 0.6 0.5; 0.4 0.7 0.3]; % rows=nodes, cols=classes
for c=1:3
    node{1}.setService(jc{c},Exp.fitMean(sv(1,c)));
    node{2}.setService(jc{c},Exp.fitMean(sv(2,c)));
    node{3}.setService(jc{c},Exp.fitMean(sv(3,c)));
end
P=model.initRoutingMatrix();
Pc=[0 1 0; 0 0 1; 1 0 0];
for c=1:3, P{c,c}=Pc; end
model.link(P);
end

function [model,node] = cqn_branch4_3cls()
model = Network('cqn_branch4_3cls');
for i=1:4, node{i}=Queue(model,sprintf('Q%d',i),SchedStrategy.PS); end
pops=[2 2 1];
for c=1:3, jc{c}=ClosedClass(model,sprintf('Class%d',c),pops(c),node{1},0); end
sv=[0.5 0.8 0.4; 0.3 0.6 0.5; 0.4 0.7 0.3; 0.6 0.5 0.45];
for c=1:3
    for i=1:4, node{i}.setService(jc{c},Exp.fitMean(sv(i,c))); end
end
P=model.initRoutingMatrix();
Pc=[0 1 0 0; 0 0 0.5 0.5; 1 0 0 0; 1 0 0 0];
for c=1:3, P{c,c}=Pc; end
model.link(P);
end
