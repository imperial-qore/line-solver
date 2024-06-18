function [XN,UN,QN,RN,TN,CN,tranSysState,tranSync]=solver_ssa_analyzer_taussa(network, options, tl_option, tau)
% [XN,UN,QN,RN,TN,CN]=SOLVER_SSA_ANALYZER_SERIAL(SN, OPTIONS)

jmodel = LINE2JLINE(network);
jsolver = JLINE.SolverTauSSA(jmodel);
import jline.solvers.ssa.*;
import jline.solvers.ssa.strategies.*;

M = jmodel.getNumberOfStatefulNodes; %number of stations
K = jmodel.getNumberOfClasses;    %number of classes

jsolver.options.samples = options.samples;
jsolver.options.seed = options.seed;
%jsolver.R5(7); % example_closedModel_1 gives wrong results with this enabled
%jsolver.MSER5(); % example_closedModel_1 gives wrong results with this enabled

if tl_option == 1    
    if tau <= 0
        tau = 1.5/jmodel.avgRate;
    end
    jsolver.tauLeapingType = TauLeapingType(TauLeapingVarType.Poisson, TauLeapingOrderStrategy.DirectedCycle, TauLeapingStateStrategy.Cutoff, tau);
    jsolver.useTauLeap = true;
elseif tl_option == 2    
    if tau <= 0
        tau = 1.5/jmodel.avgRate;
    end
    jsolver.tauLeapingType = TauLeapingType(TauLeapingVarType.Poisson, TauLeapingOrderStrategy.RandomEvent, TauLeapingStateStrategy.Cutoff, tau);
    jsolver.useTauLeap = true;
end

timeline = jsolver.solve();

tranSysState = [];
tranSync = [];

XN = NaN*zeros(1,K);
UN = NaN*zeros(M,K);
QN = NaN*zeros(M,K);
RN = NaN*zeros(M,K);
TN = NaN*zeros(M,K);
CN = NaN*zeros(1,K);

for n=1:M
    for m=1:K
        metrics = timeline.getMetrics(n-1, m-1);
        QN(n,m) = metrics.getMetricValueByName("Queue Length");
        UN(n,m) = metrics.getMetricValueByName("Utilization");
        TN(n,m) = metrics.getMetricValueByName("Throughput");
        RN(n,m) = metrics.getMetricValueByName("Response Time");
    end
end

QN(isnan(QN))=0;
CN(isnan(CN))=0;
RN(isnan(RN))=0;
UN(isnan(UN))=0;
XN(isnan(XN))=0;
TN(isnan(TN))=0;
end
