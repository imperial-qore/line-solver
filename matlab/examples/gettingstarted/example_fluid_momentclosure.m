function example_fluid_momentclosure()
% EXAMPLE_FLUID_MOMENTCLOSURE Second-order moment closures in SolverFLD.
%
% The default fluid methods close the moment hierarchy at first order: the
% drift of the mean uses min(E[X],c) in place of E[min(X,c)], so no second
% moment is ever computed and the mean is biased where min() bends. This
% example runs the three second-order methods against the exact CTMC on a
% closed two-station model swept through saturation, where that bias peaks.
%
%   'minnormal'  min-normal closure (Guenther, Stefanek, Bradley), mean and
%                covariance solved self-consistently
%   'refined'    O(1/N) refined mean field on the mean-field fixed point
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

Nvals = [2 4 6 8];
methods = {'closing','minnormal','refined'};

fprintf('Closed model: Delay(Z=1) -> Queue(PS, mu=1, c=2), sweeping population\n\n');
fprintf('%4s %10s %10s %10s %10s\n','N','CTMC Q2','closing','minnormal','refined');

for N = Nvals
    model = local_model(N);
    ctmc = SolverCTMC(model).getAvgTable;
    qexact = ctmc.QLen(2);

    row = zeros(1,numel(methods));
    for m = 1:numel(methods)
        solver = SolverFLD(model,'method',methods{m});
        avg = solver.getAvgTable;
        row(m) = avg.QLen(2);
    end
    fprintf('%4d %10.4f %10.4f %10.4f %10.4f\n', N, qexact, row);
end

% the covariance is only produced by the second-order methods; the reference
% is the stationary distribution of the birth-death chain this model reduces
% to, with birth rate (N-n) and death rate min(n,2)
fprintf('\nQueue-length standard deviation at the server (exact vs closures)\n');
fprintf('%4s %10s %10s %10s\n','N','exact','minnormal','refined');
for N = Nvals
    model = local_model(N);
    p = ones(1,N+1);
    for n = 1:N
        p(n+1) = p(n) * (N-n+1) / min(n,2);
    end
    p = p/sum(p);
    n = 0:N;
    stdExact = sqrt(sum(p.*n.^2) - sum(p.*n)^2);

    row = zeros(1,2);
    ms = {'minnormal','refined'};
    for m = 1:2
        solver = SolverFLD(model,'method',ms{m});
        solver.getAvg();
        mom = solver.getMoments();
        row(m) = mom.QStd(2,1);
    end
    fprintf('%4d %10.4f %10.4f %10.4f\n', N, stdExact, row);
end

% the same closure carries a limited load dependence: alpha(n) multiplies the
% scheduling share, so the closed term becomes E[min(X,c)*alpha(X)]. Only the
% closing family evaluates alpha; the other FLD methods refuse the model.
alpha = [1.0 1.7 2.2 2.5 2.6 2.65];
fprintf('\nLoad-dependent server, alpha = %s\n', mat2str(alpha));
fprintf('%4s %10s %10s %10s %10s\n','N','exact','closing','minnormal','refined');
for N = 2:6
    model = local_ld_model(N, alpha);
    qexact = SolverCTMC(model).getAvgTable.QLen(2);
    row = zeros(1,3); ms = {'closing','minnormal','refined'};
    for m = 1:3
        row(m) = SolverFLD(model,'method',ms{m}).getAvgTable.QLen(2);
    end
    fprintf('%4d %10.4f %10.4f %10.4f %10.4f\n', N, qexact, row);
end

% min() is not the only non-linear rate term. The capacity share of a DPS
% station, w_k*X_k/sum_j w_j*X_j, is a RATIO of populations, so evaluating it
% at the mean is a second closure: it biases the split towards the class with
% the larger weight while leaving the station total correct. The same
% covariance closes it, so the per-class utilization improves without the
% aggregate moving.
fprintf('\nDPS server, per-class utilization split (weights [1 w2])\n');
fprintf('%4s %21s %21s %21s\n','w2','exact','closing','minnormal');
for w2 = [1 2 4 8]
    model = local_dps_model(w2);
    ue = SolverCTMC(model).getAvgTable.Util(3:4)';
    uc = SolverFLD(model,'method','closing').getAvgTable.Util(3:4)';
    ug = SolverFLD(model,'method','minnormal').getAvgTable.Util(3:4)';
    fprintf('%4d %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n', w2, ue, uc, ug);
end

% GPS is the discipline where the second moment is not a correction but the
% ENTIRE mechanism. GPS divides the server by weight among the BACKLOGGED
% classes, so its share depends on the backlog INDICATOR, not on populations.
% A first-order closure cannot express it: with continuous x_k > 0 every class
% is always backlogged and the share collapses to the constant w_k/sum(w),
% which is the heavy-traffic limit and is wrong at any other load. The closure
% enumerates the 2^K backlog patterns weighted by P(N_k >= 1).
fprintf('\nGPS server, per-class utilization split (weights [1 w2])\n');
fprintf('%4s %21s %21s %21s\n','w2','exact','minnormal','first-order const');
for w2 = [1 2 4 8]
    model = local_gps_model(w2);
    ue = SolverCTMC(model).getAvgTable.Util(3:4)';
    ug = SolverFLD(model,'method','minnormal').getAvgTable.Util(3:4)';
    uc = [1 w2]/(1+w2);
    fprintf('%4d %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n', w2, ue, ug, uc);
end
end

function model = local_gps_model(w2)
model = Network('momentclosure_gps');
delay = Delay(model,'Think');
queue = Queue(model,'Server',SchedStrategy.GPS);
queue.setNumberOfServers(1);
c1 = ClosedClass(model,'Class1',2,delay);
c2 = ClosedClass(model,'Class2',2,delay);
delay.setService(c1,Exp(1.0)); delay.setService(c2,Exp(1.0));
queue.setService(c1,Exp(1.0)); queue.setService(c2,Exp(1.0));
queue.setStrategyParam(c1,1); queue.setStrategyParam(c2,w2);
P = model.initRoutingMatrix;
P{1} = Network.serialRouting(delay,queue);
P{2} = Network.serialRouting(delay,queue);
model.link(P);
end

function model = local_dps_model(w2)
model = Network('momentclosure_dps');
delay = Delay(model,'Think');
queue = Queue(model,'Server',SchedStrategy.DPS);
queue.setNumberOfServers(1);
c1 = ClosedClass(model,'Class1',2,delay);
c2 = ClosedClass(model,'Class2',2,delay);
delay.setService(c1,Exp(1.0)); delay.setService(c2,Exp(1.0));
queue.setService(c1,Exp(1.0)); queue.setService(c2,Exp(1.0));
queue.setStrategyParam(c1,1); queue.setStrategyParam(c2,w2);
P = model.initRoutingMatrix;
P{1} = Network.serialRouting(delay,queue);
P{2} = Network.serialRouting(delay,queue);
model.link(P);
end

function model = local_ld_model(N, alpha)
model = Network('momentclosure_ld');
delay = Delay(model,'Think');
queue = Queue(model,'Server',SchedStrategy.PS);
queue.setLoadDependence(alpha);
jobclass = ClosedClass(model,'Class1',N,delay);
delay.setService(jobclass,Exp(1.0));
queue.setService(jobclass,Exp(1.0));
model.link(Network.serialRouting(delay,queue));
end

function model = local_model(N)
model = Network('momentclosure');
delay = Delay(model,'Think');
queue = Queue(model,'Server',SchedStrategy.PS);
queue.setNumberOfServers(2);
jobclass = ClosedClass(model,'Class1',N,delay);
delay.setService(jobclass,Exp(1.0));
queue.setService(jobclass,Exp(1.0));
model.link(Network.serialRouting(delay,queue));
end
