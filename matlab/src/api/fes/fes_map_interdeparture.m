%{ @file fes_map_interdeparture.m
 %  @brief MAP of the inter-departure times of a two-resource subnetwork
 %
 %  @author LINE Development Team
%}

%{
 % @brief Builds the MAP (T0,T1) of the inter-departure times of a closed
 % subnetwork made of one MAP station and one MAP flow-equivalent server
 %
 % @details
 % Implements the block bidiagonal construction of Casale, Mi, Cherkasova
 % and Smirni, IEEE Trans. Soft. Eng. 37(5), 2011, Section 5.2.2. The
 % subnetwork holds n jobs that circulate between a station with MAP
 % service (D0,D1) and a flow-equivalent server whose MAP (F0^k,F1^k)
 % depends on the number k of jobs it holds. Level k of (T0,T1) is the
 % population of the flow-equivalent server, so the station holds n-k jobs
 % and both processes may be load dependent. Marked transitions are the
 % completions of the station, which are the departures fed to the rest of
 % the model.
 %
 %   T0 = [ D0 (x) I        0              0            0
 %          I (x) F1^1   D0 (+) F0^1       0            0
 %            ...           ...           ...          ...
 %            0        I (x) F1^{n-1} D0 (+) F0^{n-1}   0
 %            0            0        I (x) F1^n     I (x) F0^n ]
 %
 %   T1 = superdiagonal blocks D1 (x) I, last block row zero
 %
 % @par Syntax:
 % @code
 % [T0,T1] = fes_map_interdeparture(MAPs, FES, n)
 % [T0,T1] = fes_map_interdeparture(MAPs, FES, n, mi)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>MAPs<td>Service process of the station, either {D0,D1} or a cell
 %                 array with MAPs{j}={D0,D1} active when it holds j jobs
 % <tr><td>FES<td>Flow-equivalent server, either {F0,F1} or a cell array
 %                with FES{k}={F0,F1} active when it holds k jobs
 % <tr><td>n<td>Number of jobs circulating in the subnetwork
 % <tr><td>mi<td>(Optional) [mi_station mi_fes] servers used to scale a load
 %               independent process, Inf for a delay, default [1 1]
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>T0<td>Hidden transitions of the inter-departure MAP, sparse
 % <tr><td>T1<td>Marked transitions of the inter-departure MAP, sparse
 % </table>
 %
 % @see fes_map_moments, fes_map_aggregate, fes_map_levels
%}
function [T0,T1] = fes_map_interdeparture(MAPs, FES, n, mi)

if nargin < 3
    line_error(mfilename,'Three input arguments are required.');
end
if nargin < 4
    mi = [1 1];
end
if n < 1
    line_error(mfilename,'The subnetwork population n must be at least 1.');
end

MAPs = fes_map_levels(MAPs, n, mi(1));
FES = fes_map_levels(FES, n, mi(2));
ms = size(MAPs{1}{1},1);
mf = size(FES{1}{1},1);

blk = ms*mf;
dim = (n+1)*blk;
Ims = speye(ms);
Imf = speye(mf);

% triplet buffers, assembled once to avoid repeated sparse reallocation
i0 = []; j0 = []; v0 = [];
i1 = []; j1 = []; v1 = [];

for k = 0:n
    off = k*blk;
    j = n - k;
    if j > 0 && k > 0
        diagBlk = kron(MAPs{j}{1}, Imf) + kron(Ims, FES{k}{1});
    elseif j > 0
        diagBlk = kron(MAPs{j}{1}, Imf);
    else
        diagBlk = kron(Ims, FES{k}{1});
    end
    [bi,bj,bv] = find(diagBlk);
    i0 = [i0; off+bi]; j0 = [j0; off+bj]; v0 = [v0; bv];

    if k > 0
        [bi,bj,bv] = find(kron(Ims, FES{k}{2}));
        i0 = [i0; off+bi]; j0 = [j0; off-blk+bj]; v0 = [v0; bv];
    end

    if j > 0
        [bi,bj,bv] = find(kron(MAPs{j}{2}, Imf));
        i1 = [i1; off+bi]; j1 = [j1; off+blk+bj]; v1 = [v1; bv];
    end
end

T0 = sparse(i0, j0, v0, dim, dim);
T1 = sparse(i1, j1, v1, dim, dim);
end
