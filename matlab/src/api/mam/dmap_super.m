%{ @file dmap_super.m
 %  @brief Superposition of two discrete-time arrival streams
 %
 %  @author LINE Development Team
%}

%{
 % @brief Superposition of two discrete-time arrival streams
 %
 % @details
 % Superposing two DMAPs on a slotted time scale is NOT a DMAP: both streams
 % can fire in the same slot, so the merged stream carries batches. With
 % A = {A_0,...,A_p} and B = {B_0,...,B_q} the merged batch matrices are
 % E_k = sum_{i+j=k} kron(A_i, B_j), which is a discrete batch MAP of order
 % (p+q). Folding E_2 into E_1 would conserve neither the arrival rate nor the
 % slot in which the work appears, so the batch dimension is kept and the
 % downstream station is solved as an M/G/1-type chain instead of a QBD.
 %
 % @par Syntax:
 % @code
 % E = dmap_super(A, B)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>A<td>Cell {A_0, A_1, ...} of batch matrices (DMAP is {D0,D1})
 % <tr><td>B<td>Cell {B_0, B_1, ...} of batch matrices
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>E<td>Cell {E_0, E_1, ...} of merged batch matrices
 % </table>
%}
function E = dmap_super(A, B)

p = length(A) - 1;
q = length(B) - 1;
E = cell(1, p + q + 1);
for k = 0:(p+q)
    Ek = [];
    for i = max(0, k-q):min(p, k)
        term = kron(A{i+1}, B{k-i+1});
        if isempty(Ek)
            Ek = term;
        else
            Ek = Ek + term;
        end
    end
    E{k+1} = Ek;
end

end
