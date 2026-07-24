function test_fes_exactness()
% TEST_FES_EXACTNESS  Exactness of FES aggregation on product-form multiclass nets.
%
% Reference oracle for jar/src/test/java/jline/api/fes/FESExactnessTest.java.
%
% Kritzinger, van Wyk and Krzesinski (1982), "A generalisation of Norton's
% theorem for multiclass queueing networks": replacing a subset of centres in
% a locally balanced network by the composite centre with state-dependent
% per-class throughput tau_r(n1,...,nR) leaves the marginal distribution of the
% complement centres unchanged. For a product-form model the per-class metrics
% at the complement stations must match the exact full-model solution.
%
% Result: EXACT (machine precision) for single-exit (tandem) subsets; NOT exact
% for multi-exit subsets because aggregateFES weights the FES exit routing
% uniformly instead of using the exact relative arrival rates xi_cr (Eq. 4.6).

tol = 1e-6;

% --- single-exit (tandem) cases: must be exact ---
[m,nd] = local_tandem();
check_case('tandem_aggQ1Q2',   m, nd([2 3]),   {'Delay','Q3'}, tol, true);
[m,nd] = local_tandem();
check_case('tandem_aggQ2Q3',   m, nd([3 4]),   {'Delay','Q1'}, tol, true);
[m,nd] = local_tandem();
check_case('tandem_aggQ1Q2Q3', m, nd([2 3 4]), {'Delay'},      tol, true);

% --- multi-exit case: also exact for a closed PF network (Kritzinger, arbitrary
%     entry/exit) once the FES exit routing is visit-ratio weighted ---
[m,nd] = local_prob();
check_case('multiExit_aggQ1Q2', m, nd([2 3]), {'Delay','Q3'}, tol, true);

fprintf('test_fes_exactness: PASS\n');
end

function check_case(name, model, subsetNodes, complementNames, tol, mustBeExact)
Afull = NC(model,'method','exact').getAvgTable;
[fesModel,~,~] = ModelAdapter.aggregateFES(model, subsetNodes, ...
    struct('solver','mva','verbose',false));
Ared = NC(fesModel,'method','exact').getAvgTable;
err = local_maxrelerr(Afull, Ared, complementNames);
if mustBeExact
    assert(err < tol, sprintf('%s: expected EXACT but maxRelErr=%.3e', name, err));
    fprintf('  %-16s EXACT (maxRelErr=%.2e)\n', name, err);
else
    % Known limitation: assert it is currently NOT exact so a future fix flips this.
    assert(err > tol, sprintf(['%s: expected inexact (uniform exit-routing) but ' ...
        'maxRelErr=%.3e is within tol -- has the aggregator been fixed?'], name, err));
    fprintf('  %-16s KNOWN-INEXACT (maxRelErr=%.2e) [multi-exit exit-routing]\n', name, err);
end
end

function err = local_maxrelerr(A, B, stationNames)
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

function [model,node] = local_tandem()
model = Network('PF_tandem');
node{1}=Delay(model,'Delay');
node{2}=Queue(model,'Q1',SchedStrategy.PS);
node{3}=Queue(model,'Q2',SchedStrategy.PS);
node{4}=Queue(model,'Q3',SchedStrategy.PS);
jc{1}=ClosedClass(model,'Class1',3,node{1},0);
jc{2}=ClosedClass(model,'Class2',2,node{1},0);
node{1}.setService(jc{1},Exp.fitMean(1.0)); node{1}.setService(jc{2},Exp.fitMean(1.5));
node{2}.setService(jc{1},Exp.fitMean(0.5)); node{2}.setService(jc{2},Exp.fitMean(0.8));
node{3}.setService(jc{1},Exp.fitMean(0.3)); node{3}.setService(jc{2},Exp.fitMean(0.6));
node{4}.setService(jc{1},Exp.fitMean(0.4)); node{4}.setService(jc{2},Exp.fitMean(0.7));
P=model.initRoutingMatrix();
P{1,1}=[0 1 0 0;0 0 1 0;0 0 0 1;1 0 0 0]; P{2,2}=P{1,1};
model.link(P);
end

function [model,node] = local_prob()
model = Network('PF_prob');
node{1}=Delay(model,'Delay');
node{2}=Queue(model,'Q1',SchedStrategy.PS);
node{3}=Queue(model,'Q2',SchedStrategy.PS);
node{4}=Queue(model,'Q3',SchedStrategy.PS);
jc{1}=ClosedClass(model,'Class1',2,node{1},0);
jc{2}=ClosedClass(model,'Class2',2,node{1},0);
node{1}.setService(jc{1},Exp.fitMean(1.0)); node{1}.setService(jc{2},Exp.fitMean(1.2));
node{2}.setService(jc{1},Exp.fitMean(0.5)); node{2}.setService(jc{2},Exp.fitMean(0.7));
node{3}.setService(jc{1},Exp.fitMean(0.4)); node{3}.setService(jc{2},Exp.fitMean(0.5));
node{4}.setService(jc{1},Exp.fitMean(0.6)); node{4}.setService(jc{2},Exp.fitMean(0.3));
P=model.initRoutingMatrix();
Pc=[0 1 0 0; 0 0 0.6 0.4; 1 0 0 0; 1 0 0 0];
P{1,1}=Pc; P{2,2}=Pc;
model.link(P);
end
