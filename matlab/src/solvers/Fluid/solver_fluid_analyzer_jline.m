function [QN, UN, RN, TN, CN, XN, t, QNt, UNt, TNt, xvec] = solver_fluid_analyzer_jline(network, options)
% [QN, UN, RN, TN, CN, XN, T, QNT, UNT, TNT, XVEC] = SOLVER_FLUID_ANALYZER_JLINE(QN, OPTIONS)

% Copyright (c) 2012-2023, Imperial College London
% All rights reserved.

%%%% Returning Result from JLINE %%%%

jmodel = LINE2JLINE(network);
jsolver = JLINE.SolverFluid(jmodel);
import jline.solvers.fluid.*;

jsolver.options.method = options.method;
jsolver.options.stiff = options.stiff;
result = jsolver.runMethodSpecificAnalyzerViaLINE();

%%%% Migrating Result from JLINE SolverResult to native MatLab data structures %%%%

M = jmodel.getNumberOfStatefulNodes; %number of stations
K = jmodel.getNumberOfClasses;       %number of classes

QN = NaN*zeros(M,K);
UN = NaN*zeros(M,K);
RN = NaN*zeros(M,K);
TN = NaN*zeros(M,K);
CN = NaN*zeros(1,K);
XN = NaN*zeros(1,K);

QNt = cell(M,K);
UNt = cell(M,K);
TNt = cell(M,K);

Tmax = result.t.length();
t = NaN*zeros(Tmax, 1); 

for i=1:M 
    for j=1:K
        QN(i,j) = result.QN.get(i-1, j-1);
        UN(i,j) = result.UN.get(i-1, j-1);
        RN(i,j) = result.RN.get(i-1, j-1);
        TN(i,j) = result.TN.get(i-1, j-1);
    end
end

for j=1:K
    CN(1,j) = result.CN.get(0, j-1);
    XN(1,j) = result.XN.get(0, j-1);
end

for i=1:M 
    for j=1:K
        for p=1:Tmax
            QNt{i,j}(p,1) = result.QNt(i,j).get(p-1, 0);
            UNt{i,j}(p,1) = result.UNt(i,j).get(p-1, 0);
            TNt{i,j}(p,1) = result.TNt(i,j).get(p-1, 0);
        end
    end
end

for p=1:Tmax
    t(p,1) = result.t.get(p-1, 0);
end

% NOTE: JLINE designed such that odeStateVec is not returned to LINE
xvec.odeStateVec = [];
xvec.sn = network;
end