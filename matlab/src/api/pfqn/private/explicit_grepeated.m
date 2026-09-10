%{
%{
 % @file explicit_grepeated.m
 % @brief Single-class fixed-rate normalizing constant at demands of arbitrary multiplicity.
%}
%}

%{
%{
 % @brief Single-class fixed-rate normalizing constant at demands of arbitrary multiplicity.
 %
 % Eq. (16) of Casale, "Accelerating Performance Inference over Closed Systems by
 % Asymptotic Methods", ACM SIGMETRICS 2017. With K' distinct demands th_j of
 % multiplicity m_j the general partial-fraction expansion is
 %
 %    g(Nt) = sum_j (-1)^(m_j-1) th_j^(Nt+K-m_j)
 %            * sum_{r>=0, |r|=m_j-1} (-1)^r_j nchoosek(Nt+r_j,r_j)
 %              prod_{k~=j} nchoosek(m_k+r_k-1,r_k) th_k^r_k / (th_j-th_k)^(m_k+r_k)
 %
 % which reduces to Eq. (15) when every m_j is one. Demands within tol of each
 % other, relatively to the largest one, are merged into one distinct value
 % carrying their count. Shared by pfqn_explicit and pfqn_explicit_ld.
 %
 % @fn explicit_grepeated(th, Nt, K, tol)
 % @param th Demands (Kx1), repetitions allowed.
 % @param Nt Total population.
 % @param K Number of queues.
 % @param tol Relative tolerance declaring two demands redundant.
 % @return lg Logarithm of |g|.
 % @return sg Sign of g.
 % @return lossDigits Decimal digits lost to cancellation.
%}
%}
function [lg,sg,lossDigits] = explicit_grepeated(th,Nt,K,tol)
ths = sort(th(:));
scale = ths(end);
if scale<=0
    scale = 1;
end
gid = cumsum([1; diff(ths) > tol*scale]);
Kp = gid(end);
thd = zeros(Kp,1);
m = zeros(Kp,1);
for j=1:Kp
    sel = (gid==j);
    thd(j) = mean(ths(sel)); % the centroid represents a cluster of near-ties
    m(j) = sum(sel);
end
lin = [];
sgv = [];
for j=1:Kp
    if thd(j)<=0
        % the exponent Nt+K-m_j is at least Nt>=1, so a zero cluster
        % contributes nothing
        continue
    end
    louter = (Nt+K-m(j))*log(thd(j));
    souter = (-1)^(m(j)-1);
    rs = multichoose(Kp,m(j)-1); % every K'-vector r>=0 with sum(r)=m_j-1
    for i=1:size(rs,1)
        r = rs(i,:);
        lval = louter + nchoosekln(Nt+r(j),r(j));
        sval = souter*(-1)^r(j);
        for k=[1:(j-1),(j+1):Kp]
            lval = lval + nchoosekln(m(k)+r(k)-1,r(k));
            if r(k)>0
                if thd(k)<=0
                    lval = -Inf; % theta_k^r_k vanishes, 0^0=1 is the r_k=0 case
                    break
                end
                lval = lval + r(k)*log(thd(k));
            end
            dd = thd(j)-thd(k);
            lval = lval - (m(k)+r(k))*log(abs(dd));
            sval = sval*sign(dd)^(m(k)+r(k));
        end
        lin(end+1,1) = lval; %#ok<AGROW>
        sgv(end+1,1) = sval; %#ok<AGROW>
    end
end
[lg,sg,lossDigits] = explicit_signedlogsumexp(lin,sgv);
end
