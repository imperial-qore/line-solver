%{ @file cache_erec.m
 %  @brief Recursive computation of the normalizing constant for cache models
 %
 %  @author LINE Development Team
%}

%{
 % @brief Computes the normalizing constant for cache models recursively
 %
 % @details
 % This function recursively computes the normalizing constant E for
 % cache models with given item popularity probabilities and cache capacity.
 % With item sizes and per-list storage cost caps it evaluates instead the
 % constrained normalizing constant E(m,k) of Casale-Gast (IEEE/ACM ToN,
 % 2021), Sec. IX, using the recursion
 %   E(m,k) = E_i(m,k) + sum_j m_j gamma_ij E_i(m-1_j, k-sigma_i 1_j)
 % with the extra boundary E(m,k)=0 whenever some k_j is negative.
 %
 % @par Syntax:
 % @code
 % E = cache_erec(gamma, m)
 % E = cache_erec(gamma, m, sigma, k)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>gamma<td>Item popularity probabilities
 % <tr><td>m<td>Cache capacity
 % <tr><td>sigma<td>(Optional) item storage costs (sizes), 1 x n, positive integers
 % <tr><td>k<td>(Optional) per-list storage cost caps, 1 x h, non-negative integers
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>E<td>Normalizing constant
 % </table>
%}
function E=cache_erec(gamma,m,sigma,k)
if nargin<3 || isempty(sigma) || isempty(k)
    E = sub_cache_erec(gamma,m,length(gamma));
else
    E = sub_cache_erec_cost(gamma,m,sigma(:).',k(:).');
end
end

function E=sub_cache_erec(gamma,m,k)
h=length(m);
if sum(m)==0
    E=1;
    return
end
if sum(m)>k || min(m)<0
    E=0;
    return
end

if k==1 && sum(m)==1
    j = find(m);
    E=gamma(1,j);
    return
end
E = sub_cache_erec(gamma,m,k-1);
for j=1:h
    if m(j)>0
        E = E + gamma(k,j)*m(j)*sub_cache_erec(gamma,oner(m,j),k-1);
    end
end
end

% Dynamic program over the (residual capacity, residual cost cap) lattice.
function E=sub_cache_erec_cost(gamma,m,sigma,k)
[n,h] = size(gamma);
m = m(:).';
if numel(sigma)~=n
    line_error(mfilename,'The item size vector must have one entry per item.');
end
if numel(k)~=h
    line_error(mfilename,'The cost cap vector must have one entry per cache list.');
end
if any(sigma<=0) || any(abs(sigma-round(sigma))>0)
    line_error(mfilename,'Item sizes must be positive integers.');
end
if any(k<0) || any(abs(k-round(k))>0)
    line_error(mfilename,'Storage cost caps must be non-negative integers.');
end
if sum(m)>n || min(m)<0
    E=0;
    return
end
if sum(m)==0
    E=1;
    return
end
dims = [m+1, k+1];
lattice = prod(dims);
if lattice > 1e7
    line_error(mfilename,'The cost-constrained normalizing constant lattice has %d states, which exceeds the exact method limit; use the sampling method.', lattice);
end
% strides for a column-major linear index over dims
stride = cumprod([1, dims(1:end-1)]);
% F(mm,kk) over items 1..t; at t=0 only the empty cache contributes
F = zeros(lattice,1);
mcount = zeros(lattice,1);  % sum(mm) for each lattice point
sub = zeros(1,2*h);
for idx = 1:lattice
    rem = idx-1;
    for d = 2*h:-1:1
        sub(d) = floor(rem/stride(d));
        rem = rem - sub(d)*stride(d);
    end
    mcount(idx) = sum(sub(1:h));
end
F(mcount==0) = 1;
for t = 1:n
    Fprev = F;
    for idx = 1:lattice
        if mcount(idx)>t
            F(idx) = 0;
            continue
        end
        val = Fprev(idx);
        rem = idx-1;
        for d = 2*h:-1:1
            sub(d) = floor(rem/stride(d));
            rem = rem - sub(d)*stride(d);
        end
        for j = 1:h
            mj = sub(j);
            kj = sub(h+j);
            if mj>0 && kj>=sigma(t) && gamma(t,j)~=0
                val = val + gamma(t,j)*mj*Fprev(idx - stride(j) - sigma(t)*stride(h+j));
            end
        end
        F(idx) = val;
    end
end
E = F(end);
end
