%{ @file dmap_thin.m
 %  @brief Bernoulli thinning of a discrete-time batch arrival stream
 %
 %  @author LINE Development Team
%}

%{
 % @brief Bernoulli thinning of a discrete-time batch arrival stream
 %
 % @details
 % Routes each event of a discrete batch MAP to the branch with probability p,
 % independently across events. Since a slot may carry a batch of n events, the
 % number routed is Binomial(n,p) and the thinned matrices are
 % B_k = sum_{n>=k} nchoosek(n,k) p^k (1-p)^(n-k) A_n. The phase process is
 % untouched, so the result is exact for Bernoulli (PROB/RAND) routing.
 %
 % @par Syntax:
 % @code
 % B = dmap_thin(A, p)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>A<td>Cell {A_0, A_1, ...} of batch matrices
 % <tr><td>p<td>Routing probability in [0,1]
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>B<td>Cell {B_0, B_1, ...} of thinned batch matrices
 % </table>
%}
function B = dmap_thin(A, p)

if p < 0 || p > 1
    line_error(mfilename, 'The routing probability must lie in [0,1], got %g.', p);
end

n = length(A) - 1;
B = cell(1, n + 1);
for k = 0:n
    Bk = zeros(size(A{1}));
    for j = k:n
        w = nchoosek(j, k) * p^k * (1-p)^(j-k);
        if w > 0
            Bk = Bk + w * A{j+1};
        end
    end
    B{k+1} = Bk;
end

% trailing zero batch levels carry no mass and only inflate the M/G/1 blocks
while length(B) > 2 && all(all(abs(B{end}) < GlobalConstants.Zero))
    B(end) = [];
end

end
