%{ @file fes_map_levels.m
 %  @brief Normalizes a flow-equivalent server descriptor to one MAP per level
 %
 %  @author LINE Development Team
%}

%{
 % @brief Expands a MAP into the per-level processes of a load-dependent server
 %
 % @details
 % A flow-equivalent server is described by one MAP (F0^k,F1^k) per
 % population level k=1..n. This function accepts either that cell array,
 % which it validates and returns unchanged, or a single MAP {F0,F1} which
 % it replicates over the levels. When the number of servers mi is given,
 % the replicated MAP is scaled by min(k,mi), which reproduces a queue with
 % mi servers and, for mi=Inf, a delay station serving at rate k*mu. The
 % scaling is exact for exponential service and is the load-dependent rate
 % approximation otherwise.
 %
 % @par Syntax:
 % @code
 % FESlev = fes_map_levels(FES, n)
 % FESlev = fes_map_levels(FES, n, mi)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>FES<td>Either {F0,F1} or a cell array with FES{k}={F0,F1}
 % <tr><td>n<td>Number of levels required
 % <tr><td>mi<td>(Optional) number of servers, Inf for a delay, default 1
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>FESlev<td>Cell array with FESlev{k}={F0,F1} for k=1..n
 % </table>
 %
 % @see fes_map_interdeparture
%}
function FESlev = fes_map_levels(FES, n, mi)

if nargin < 3
    mi = 1;
end

if numel(FES) == 2 && isnumeric(FES{1}) && isnumeric(FES{2})
    FESlev = cell(1,n);
    for k = 1:n
        s = min(k, mi);
        FESlev{k} = {s*FES{1}, s*FES{2}};
    end
    return
end

if numel(FES) < n
    line_error(mfilename, sprintf('The flow-equivalent server is defined for %d levels but %d are required.', numel(FES), n));
end

FESlev = cell(1,n);
mf = size(FES{1}{1},1);
for k = 1:n
    if size(FES{k}{1},1) ~= mf
        line_error(mfilename,'All levels of a flow-equivalent server must have the same number of phases.');
    end
    FESlev{k} = FES{k};
end
end
