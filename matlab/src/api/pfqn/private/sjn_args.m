%{
%{
 % @file sjn_args.m
 % @brief Argument checking and defaults shared by the SJN solvers.
%}
%}

function [M,R,N,Z,scv,sjnset,V,S,options] = sjn_args(caller,L,N,Z,scv,sjnset,V,options)
%{
%{
 % @brief Normalise the arguments of pfqn_mvasjn and pfqn_amvasjn, which take
 %        the same model description, and derive the per-visit service times
 %        S = L./V that the shortest-job-next equations are written in (the
 %        job size the discipline compares is one visit's service time, not
 %        the service demand accumulated over all visits).
 % @fn sjn_args(caller, L, N, Z, scv, sjnset, V, options)
 % @param caller Name of the calling function, used in error messages.
 % @param L Service demand matrix (M x R).
 % @param N Population vector (1 x R).
 % @param Z Think time vector (1 x R), may be empty.
 % @param scv Squared coefficients of variation (M x R), may be empty.
 % @param sjnset Indices of the SJN stations, may be empty.
 % @param V Visit ratios (M x R), may be empty.
 % @param options Struct of solver options, may be empty.
 % @return M Number of queueing stations.
 % @return R Number of classes.
 % @return N Population vector, rounded to integers and reshaped.
 % @return Z Think times, reshaped.
 % @return scv Squared coefficients of variation, defaulted to ones.
 % @return sjnset Row vector of SJN station indices.
 % @return V Visit ratios, defaulted to ones.
 % @return S Per-visit service times (M x R).
 % @return options Options with ns, Lfactor, prio, tol and iter_max filled in.
%}
%}
[M,R] = size(L);
N = round(N(:)');
if length(N) ~= R
    line_error(caller,'demand matrix and population vector have different number of classes');
end
if any(N < 0)
    line_error(caller,'negative class populations');
end
if isempty(Z)
    Z = zeros(1,R);
end
Z = Z(:)';
if isempty(scv)
    scv = ones(M,R);
end
if isempty(sjnset)
    sjnset = [];
end
sjnset = sjnset(:)';
if any(sjnset < 1 | sjnset > M)
    line_error(caller,'sjnset contains a station index outside 1..M');
end
if length(unique(sjnset)) ~= length(sjnset)
    line_error(caller,'sjnset repeats a station index');
end
if isempty(V)
    V = ones(M,R);
end
S = zeros(M,R);
nz = V > 0;
S(nz) = L(nz) ./ V(nz);
if ~isfield(options,'ns') || isempty(options.ns)
    options.ns = 32;
end
if mod(options.ns,2) ~= 0
    line_error(caller,'options.ns must be even, composite Simpson integrates over panels of two subdivisions');
end
if ~isfield(options,'Lfactor') || isempty(options.Lfactor)
    options.Lfactor = 8;
end
if ~isfield(options,'prio')
    options.prio = [];
end
if ~isfield(options,'tol') || isempty(options.tol)
    options.tol = 1e-8;
end
if ~isfield(options,'iter_max') || isempty(options.iter_max)
    options.iter_max = 1000;
end
if ~isfield(options,'umax') || isempty(options.umax)
    options.umax = 0.999;
end
if options.umax <= 0 || options.umax >= 1
    line_error(caller,'options.umax must lie strictly between zero and one');
end
if ~isempty(options.prio)
    options.prio = options.prio(:)';
    if length(options.prio) ~= R
        line_error(caller,'options.prio must have one priority level per class');
    end
    if length(unique(options.prio)) ~= R
        line_error(caller,'options.prio must assign distinct levels, ties across classes are not covered by the SJN priority equations');
    end
end
end
