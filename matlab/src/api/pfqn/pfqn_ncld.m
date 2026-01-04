%{
%{
 % @file pfqn_ncld.m
 % @brief Normalizing constant for load-dependent closed networks.
%}
%}

%{
%{
 % @brief Normalizing constant for load-dependent closed networks.
 % @fn pfqn_ncld(L, N, Z, mu, varargin)
 % @param L Service demand matrix.
 % @param N Population vector.
 % @param Z Think time vector.
 % @param mu Load-dependent rate matrix.
 % @param varargin Optional solver parameters.
 % @return lG Logarithm of normalizing constant.
 % @return G Normalizing constant.
 % @return method Method used for computation.
%}
%}
function [lG,G,method] = pfqn_ncld(L,N,Z,mu,varargin)
% [LGN,G,METHOD] = PFQN_NCLD(L,N,Z,VARARGIN)

options = Solver.parseOptions(varargin, SolverNC.defaultOptions);
lG = NaN;
G = NaN;
method = options.method;

% backup initial parameters

mu = mu(:,1:sum(N));
% first remove empty classes
nnzClasses = find(N);
L = L(:,nnzClasses);
N = N(:,nnzClasses);
Z = Z(:,nnzClasses);

% then scale demands in [0,1], importat that stays before the other
% simplications in case both D and Z are all very small or very large in a
% given class, in which case the may look to filter but not if all of them
% are at the same scale
R = length(N);
scalevec = ones(1,R);
for r=1:R
    scalevec(r) = max([L(:,r);Z(:,r)]);
end
L = L ./ repmat(scalevec,size(L,1),1);
Z = Z ./ scalevec;

% remove stations with no demand
Lsum = sum(L,2);
Lmax = max(L,[],2);
demStations = find((Lmax./Lsum)>GlobalConstants.FineTol);
L = L(demStations,:);
mu = mu(demStations,:);

% if there is a class with jobs but with L and Z all zero
if any(N((sum(L,1) + sum(Z,1)) == 0)>0)
    line_warning(mfilename,'The model has no positive demands in any class.\n');
    if isempty(Z) || sum(Z(:))<options.tol
        lG = 0;
    else
        lG = - sum(factln(N)) + sum(N.*log(sum(Z,1))) + N*log(scalevec)';
    end
    G = NaN;

    return
end

% update M and R
[M,R]=size(L);

% return immediately if the model is a degenerate case
if isempty(L) || sum(L(:))<options.tol % all demands are zero
    if isempty(Z) || sum(Z(:))<options.tol
        lG = 0;
    else
        lG = - sum(factln(N)) + sum(N.*log(sum(Z,1))) + N*log(scalevec)';
    end
    return
elseif M==1 && (isempty(Z) || sum(Z(:))<options.tol) % single node and no think time
    lG = factln(sum(N)) - sum(factln(N)) + sum(N.*log(sum(L,1))) + N*log(scalevec)' - sum(log(mu(1:sum(N))));
    return
end

% determine contribution from jobs that permanently loop at delay
zeroDemandClasses = find(sum(L,1)<options.tol); % all jobs in delay
nonzeroDemandClasses = setdiff(1:R, zeroDemandClasses);

if isempty(sum(Z,1)) || all(sum(Z(:,zeroDemandClasses),1)<options.tol)
    lGzdem = 0;
    Nz = 0;
else
    if isempty(zeroDemandClasses) % for old MATLAB release compatibility
        lGzdem = 0;
        Nz = 0;
    else
        Nz = N(zeroDemandClasses);
        lGzdem = - sum(factln(Nz)) + sum(Nz.*log(sum(Z(:,zeroDemandClasses),1))) + Nz*log(scalevec(zeroDemandClasses))';
    end
end
L = L(:,nonzeroDemandClasses);
N = N(nonzeroDemandClasses);
Z = Z(:,nonzeroDemandClasses);
scalevecz = scalevec(nonzeroDemandClasses);
% compute G for classes No with non-zero demand
if any(N<0)
    lGnnzdem = 0;
else
    [lGnnzdem,method] = compute_norm_const_ld(L, N, Z, mu, options);
end
% scale back to original demands
lG = lGnnzdem + lGzdem + N*log(scalevecz)';
G = exp(lG);
end

function [lG,method] = compute_norm_const_ld(L,N,Z,mu,options)
% LG = COMPUTE_NORM_CONST_LD(L,N,Z,OPTIONS)
[M,R] = size(L);
method = options.method;
switch options.method
    case {'default','exact'}
        if sum(Z(:))<GlobalConstants.FineTol
            Lz = L;
            muz = mu;
        else
            D = size(Z,1); % number of delays
            Lz = [L;Z];
            muz = [mu; repmat(1:size(mu,2),D,1)];
        end
        if R==1
            [lG] = pfqn_gldsingle(Lz, N, muz, options);
            method = 'exact/gld';
        elseif M==1 && any(Z>0)
            [~,lG]= pfqn_comomrm_ld(L, N, Z, mu, options);
            method = 'exact/comomld';
        elseif M==2 && max(Z(:)) < GlobalConstants.FineTol
            [~,lG] = pfqn_comomrm_ld(L, N, 0*N, mu, options);
            method = 'exact/comomld';
        else
            [~,lG] = pfqn_gld(Lz, N, muz, options);
            method = 'exact/gld';
        end        
    case 'rd'
        [lG] = pfqn_rd(L, N, Z, mu, options);
    case 'nrp'
        [lG] = pfqn_nrp(L, N, Z, mu, options);
    case 'nrl'
        [lG] = pfqn_nrl(L, N, Z, mu, options);
    case 'comomld'
        if M<=1 || sum(Z) < GlobalConstants.Zero
            [~,lG]= pfqn_comomrm_ld(L, N, Z, mu, options);
        else
            line_warning(mfilename,'Load-dependent CoMoM is available only in models with a delay and m identical stations, running the ''rd'' algorithm instead.\n');
            [lG] = pfqn_rd(L, N, Z, mu, options);
            method = 'rd';
        end
    otherwise
        line_error(mfilename,sprintf('Unrecognized method for solving load-dependent models: %s',options.method));
end
return
end
