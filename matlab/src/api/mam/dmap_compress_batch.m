%{ @file dmap_compress_batch.m
 %  @brief Reduces the order of a discrete batch arrival stream
 %
 %  @author LINE Development Team
%}

%{
 % @brief Reduces the order of a discrete batch arrival stream
 %
 % @details
 % Superposing slotted streams multiplies phase dimensions, so a decomposition
 % over several stations has to bound them. The reduction keeps the two
 % features that the downstream M/G/1-type solve consumes: the law of the time
 % between NONEMPTY slots, matched to three moments by dmap_compress, and the
 % stationary batch-size distribution conditional on a nonempty slot, kept
 % exactly. The reduced stream is B_k = q_k * D1, B_0 = D0, whose event rate
 % equals the original one by construction.
 %
 % @par Syntax:
 % @code
 % Bc = dmap_compress_batch(B, maxOrder)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>B<td>Cell {B_0, B_1, ...} of batch matrices
 % <tr><td>maxOrder<td>Largest phase order left uncompressed
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>Bc<td>Cell {B_0, B_1, ...} of bounded phase order
 % </table>
%}
function Bc = dmap_compress_batch(B, maxOrder)

if size(B{1}, 1) <= maxOrder
    Bc = B;
    return;
end

nb = length(B) - 1;
m = size(B{1}, 1);

% stationary batch-size distribution, conditional on the slot being nonempty
Ptot = zeros(m, m);
for k = 1:length(B)
    Ptot = Ptot + B{k};
end
piPhase = dtmc_solve(Ptot);
e = ones(m, 1);
qraw = zeros(1, nb);
for k = 1:nb
    qraw(k) = piPhase * B{k+1} * e;
end
massNonEmpty = sum(qraw);
if massNonEmpty <= GlobalConstants.Zero
    line_error(mfilename, 'The batch stream carries no events, so it cannot be compressed.');
end
q = qraw / massNonEmpty;

% time between nonempty slots, compressed to a renewal DMAP
marked = {B{1}, Ptot - B{1}};
markedC = dmap_compress(marked, maxOrder);

Bc = cell(1, nb + 1);
Bc{1} = markedC{1};
for k = 1:nb
    Bc{k+1} = q(k) * markedC{2};
end

end
