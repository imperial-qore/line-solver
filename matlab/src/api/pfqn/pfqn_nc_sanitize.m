%{
%{
 % @file pfqn_nc_sanitize.m
 % @brief Sanitize and preprocess network parameters for NC solvers.
%}
%}

%{
%{
 % @brief Sanitize and preprocess network parameters for NC solvers.
 % @fn pfqn_nc_sanitize(lambda, L, N, Z, atol)
 % @param lambda Arrival rate vector.
 % @param L Service demand matrix.
 % @param N Population vector.
 % @param Z Think time vector.
 % @param atol Absolute tolerance.
 % @return lambda Sanitized arrival rates.
 % @return L Sanitized service demands (rescaled).
 % @return N Sanitized populations.
 % @return Z Sanitized think times (rescaled).
 % @return lGremaind Log normalization factor from removed classes.
%}
%}
function [lambda,L,N,Z,lGremaind] = pfqn_nc_sanitize(lambda,L,N,Z,atol)
% erase empty classes
L(isnan(L)) = 0;
Z(isnan(Z)) = 0;
nnzclasses=find(N);
L=L(:,nnzclasses);
N=N(:,nnzclasses);
Z=Z(:,nnzclasses);
lambda=lambda(:,nnzclasses);
% erase ill-defined classes
zeroclasses=find((sum(L,1)+sum(Z,1))<atol);
L(:,zeroclasses)=[];
N(:,zeroclasses)=[];
Z(:,zeroclasses)=[];
lambda(:,zeroclasses)=[];
%
lGremaind= 0;
% find zero demand classes
% see _kb/03-api-layer.md (pfqn/ family: scaling, log-domain switches, dispatch gates)
zerodemands=find(sum(L,1)<atol);
if ~isempty(zerodemands)
    lGremaind = lGremaind + N(zerodemands) * log(sum(Z(:,zerodemands),1))' - sum(factln(N(zerodemands)));
    L(:,zerodemands)=[];
    Z(:,zerodemands)=[];
    N(:,zerodemands)=[];
    lambda(:,zerodemands)=[];
end
% rescale demands
Lmax = max(L,[],1); % use L, which has been santized to always be ~=0
if isempty(Lmax)
    Lmax = ones(1,size(Z,2));
end
L = L./repmat(Lmax,size(L,1),1);
Z = Z./repmat(Lmax,size(Z,1),1);
lGremaind = lGremaind + N*log(Lmax)';
% sort from smallest to largest think time
% see _kb/03-api-layer.md (pfqn/ family: scaling, log-domain switches, dispatch gates)
if ~isempty(Z)
    [~,rsort] = sort(sum(Z,1),'ascend');
    if ~isempty(L)
        L=L(:,rsort);
    end
    Z=Z(:,rsort);
    N=N(:,rsort);
    if ~isempty(lambda)
        lambda=lambda(:,rsort);
    end
end
% ensure zero think time classes are anyway first. The test is on the COLUMN
% SUM: find() on the (K x R) matrix Z returns column-major LINEAR indices, and
% using them as column indices is only correct when Z has a single row.
zerothinktimes=find(sum(Z,1)<atol);
nonzerothinktimes = setdiff(1:size(L,2),zerothinktimes);
L=L(:,[zerothinktimes,nonzerothinktimes]);
N=N(:,[zerothinktimes,nonzerothinktimes]);
Z=Z(:,[zerothinktimes,nonzerothinktimes]);
if ~isempty(lambda)
    lambda=lambda(:,[zerothinktimes,nonzerothinktimes]);
end
end