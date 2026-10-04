function mtrace = m3afit_init(S, C, varargin)
% mtrace = m3a_init(S, C)
%
% DESCRIPTION
% Prepare multiclass trace S for M3A fitting
% 
% INPUT
% S - inter-arrival times
% C - class number for each arrival
%
% OUTPUT
% mtrace.S - inter-arrival time trace
% mtrace.C - arrival class

%% options
NumClasses = length(unique(C));
mtrace = struct('S',S,'C',C,'NumClasses',NumClasses);
end
