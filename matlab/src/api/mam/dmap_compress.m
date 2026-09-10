%{ @file dmap_compress.m
 %  @brief Reduces the order of a DMAP by matching interevent moments
 %
 %  @author LINE Development Team
%}

%{
 % @brief Reduces the order of a DMAP by matching interevent moments
 %
 % @details
 % Leaves the process untouched while its order stays within maxOrder, and
 % otherwise replaces it by the two-phase discrete phase-type law with the same
 % first three interevent moments (BuTools DPH2From3Moments), read back as a
 % renewal DMAP. Correlation is NOT preserved, which mirrors what the
 % continuous-time 'mixture.order1' compression does, and it is the reason the
 % multi-station discrete-time path is an approximation. When the moment
 % triple is outside the DPH(2) region the fallback keeps the exact mean with a
 % Geometric, so the arrival rate of the decomposition is conserved in every
 % branch.
 %
 % @par Syntax:
 % @code
 % DMAPc = dmap_compress(DMAP, maxOrder)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>DMAP<td>Cell {D0, D1}
 % <tr><td>maxOrder<td>Largest order left uncompressed
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>DMAPc<td>Cell {D0, D1} of order at most max(maxOrder,2)
 % </table>
%}
function DMAPc = dmap_compress(DMAP, maxOrder)

if size(DMAP{1}, 1) <= maxOrder
    DMAPc = DMAP;
    return;
end

moms = dmap_moment(DMAP, [1 2 3]);
DMAPc = [];
try
    [alpha, A] = DPH2From3Moments(moms);
    if all(isfinite(alpha(:))) && all(isfinite(A(:)))
        cand = dph_to_dmap(alpha, A);
        if dmap_isfeasible(cand)
            DMAPc = cand;
        end
    end
catch
    DMAPc = [];
end

if isempty(DMAPc)
    % the moment triple is not DPH(2)-feasible: keep the rate, drop the shape
    p = 1 / moms(1);
    p = min(1, max(GlobalConstants.Zero, p));
    DMAPc = dph_to_dmap(1, 1 - p);
end

end
