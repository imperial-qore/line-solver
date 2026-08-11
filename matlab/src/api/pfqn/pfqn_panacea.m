%{
%{
 % @file pfqn_panacea.m
 % @brief PANACEA (PAth-based Normal Approximation for Closed networks Estimation Algorithm).
%}
%}

%{
%{
 % @brief PANACEA (PAth-based Normal Approximation for Closed networks Estimation Algorithm).
 % @fn pfqn_panacea(L, N, Z, terms)
 % @param L Service demand matrix.
 % @param N Population vector.
 % @param Z Think time vector.
 % @param terms Number of terms in the normal-usage asymptotic series
 %        (1, 2, or 3; default 3), as selectable in the original PANACEA
 %        package (Ramakrishnan-Mitra, BSTJ 61(10):2849-2872, 1982).
 % @return Gn Normalizing constant.
 % @return lGn Logarithm of normalizing constant.
%}
%}
function [Gn,lGn]=pfqn_panacea(L,N,Z,terms)
% [GN,LGN]=PFQN_PANACEA(L,N,Z,TERMS)

% K = population vector
[q,p]=size(L);
if nargin==2 || isempty(Z)
    Z=N*0+1e-8;
end
if nargin<4 || isempty(terms)
    terms=3;
end
if ~isscalar(terms) || terms<1 || terms>3 || terms~=round(terms)
    line_error(mfilename,'The terms parameter must be 1, 2, or 3 (higher-order coefficients are not implemented).');
end
if isempty(L) | sum(L,1)==zeros(1,p)
    lGn = - sum(factln(N)) + sum(N.*log(sum(Z,1)));
    Gn=exp(lGn);
    return
end
r = L./repmat(Z,q,1);
Nt = max(1./r(r>0));  % ignore structural zeros (classes not visiting a station)
beta = N/Nt;
gamma = r * Nt;
alpha = 1-N*r';
gammatilde = gamma ./ repmat(alpha',1,p);
if min(alpha)<0
    %    line_warning(mfilename,'Model is not in normal usage');
    Gn=NaN;
    lGn=NaN;
    return
end

A0 = 1;
A1 = 0;
if terms>=2
    for j=1:p
        m = zeros(1,p); m(j)=2;
        A1 = A1 -beta(j) * pfqn_ca(gammatilde,m);
    end
end

A2 = 0;
if terms>=3
    for j=1:p
        m = zeros(1,p); m(j)=3;
        A2 = A2 + 2 * beta(j) * pfqn_ca(gammatilde,m);
        m = zeros(1,p); m(j)=4;
        A2 = A2 + 3 * beta(j)^2 * pfqn_ca(gammatilde,m);
        for k=1:p
            if k~=j
                m = zeros(1,p); m(j)=2; m(k)=2;
                A2 = A2 + 0.5 * beta(j) * beta(k) * pfqn_ca(gammatilde,m);
            end
        end
    end
end

% if false
%     A3 = 0;
%     for j=1:p
%         m = zeros(1,p); m(j)=4;
%         A3 = A3 - 6 * beta(j) * pfqn_ca(gammatilde,m);
%         m = zeros(1,p); m(j)=5;
%         A3 = A3 - 20 * beta(j)^2 * pfqn_ca(gammatilde,m);
%         m = zeros(1,p); m(j)=6;
%         A3 = A3 - 15 * beta(j)^3 * pfqn_ca(gammatilde,m);
%         for k=setdiff(1:p,j)
%             m = zeros(1,p); m(j)=4; m(k)=2;
%             A3 = A3 - 2 * beta(j) * beta(k) * pfqn_ca(gammatilde,m);
%             m = zeros(1,p); m(j)=2; m(k)=3;
%             A3 = A3 - 3 * beta(j)^2 * beta(k) * pfqn_ca(gammatilde,m);
%             for l=setdiff(1:p,[j,k])
%                 m = zeros(1,p); m(j)=2; m(k)=2; m(l)=2;
%                 A3 = A3 - (1/6) * beta(j) * beta(k) * beta(l) * pfqn_ca(gammatilde,m);
%             end
%         end
%     end
% end
I = [A0, A1/Nt, A2/Nt^2];
I = I(1:terms);
%, A3/N^3*0];

lGn = - sum(factln(N)) + sum(N.*log(sum(Z,1))) + log(sum(I)) - sum(log(alpha));
Gn = exp(lGn);
if ~isfinite(lGn)
    Gn=NaN;
    lGn=NaN;
end
end
