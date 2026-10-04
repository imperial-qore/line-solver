%{ @file qsys_serves_method.m
 % @brief Whether solver_mva_qsys_analyzer has an arm for METHOD.
%}

function bool = qsys_serves_method(method)
% BOOL = QSYS_SERVES_METHOD(METHOD)
%
% ONE PREDICATE FOR THE INTERCEPTION AND THE RUN. The Source-Queue-Sink shape is
% claimed by solver_mva_qsys_analyzer, which answers this FIXED list of closed
% forms and refuses every other name -- so each general network method
% listValidMethods offers on an open model was advertised on the one open shape
% it could not run on, and raised 'Unsupported method for a model with 1 station
% and 1 class' the moment it was asked for: mva, amva, sum, esum, lin, gflin,
% egflin, qli, fli, qd and qdlin, eleven of them. They are NETWORK methods, and
% the general branch solves a one-queue network exactly as it solves a larger
% one, so mvaDispatch stands aside for them rather than claiming a model it
% cannot answer. 'qna' was the first name found this way and used to be excluded
% by hand in mvaDispatch; it needs no special case now, having no arm here
% either.
%
% Mirrors detail::qsys_serves_method in cpp/include/line/solvers/mva/mva_dispatch.h
% and its JAR and native python twins.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

served = {'default', 'exact', 'erlanga', 'mgisrgi', 'gigk.diffusion', ...
    'mm1', 'mmk', 'mg1', 'mgi1', 'gigk', 'gigk.kingman_approx', 'gigk.whitt', ...
    'gig1', 'gig1.allen', 'gig1.kingman', 'gig1.heyman', 'gig1.kobayashi', ...
    'gig1.klb', 'gig1.marchal', 'gig1.gelenbe', 'gig1.kimura', 'gig1.extremal', ...
    'qed', 'rqna', 'rqt', 'gm1', 'gim1'};
bool = any(strcmp(method, served));
end
