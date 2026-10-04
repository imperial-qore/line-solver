%{ @file fes_map_aggregate.m
 %  @brief Recursive MAP flow-equivalent server for a station subset
 %
 %  @author LINE Development Team
%}

%{
 % @brief Aggregates a subnetwork into a load-dependent MAP flow-equivalent
 % server that reproduces mean, variability and burstiness of its output
 %
 % @details
 % Implements the recursion of Section 5.2.1 of Casale, Mi, Cherkasova and
 % Smirni, IEEE Trans. Soft. Eng. 37(5), 2011. The first station seeds the
 % flow equivalent server; every further station is folded against the
 % running server by building the inter-departure MAP of the resulting pair
 % at each population level and fitting a MAP(2) to its first three moments
 % and index of dispersion. The result is one MAP per level, which is the
 % service process of a single load-dependent station that replaces the
 % whole subnetwork. Unlike the classic flow-equivalent server, which keeps
 % only the mean throughput of the subnetwork, this one also carries the
 % burstiness of its departure stream, so a bottleneck switch across the
 % aggregated resources remains visible to the rest of the model.
 %
 % Levels are evaluated on a grid and the four descriptors are interpolated
 % between grid points, as MAPs fitted at neighbouring populations are
 % similar. The MAP is refitted at every level from the interpolated
 % descriptors, never interpolated entrywise.
 %
 % @par Syntax:
 % @code
 % [FES,info] = fes_map_aggregate(maps, servers, n)
 % [FES,info] = fes_map_aggregate(maps, servers, n, options)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>maps<td>Cell array of the station service processes {D0,D1},
 %                 already scaled by their visit ratios within the subset
 % <tr><td>servers<td>Number of servers of each station, Inf for a delay
 % <tr><td>n<td>Largest population the flow-equivalent server must serve
 % <tr><td>options<td>(Optional) struct with fields method ('ssolve' or
 %                    'euler'), grid (levels to evaluate) and verbose
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>FES<td>Cell array with FES{k}={F0,F1} for k=1..n
 % <tr><td>info<td>Struct with fields throughput (1/e1 per level), moments
 %                 (4 x n descriptors e1,e2,e3,idc), status (fit status per
 %                 level) and grid (levels actually evaluated)
 % </table>
 %
 % @see fes_map_interdeparture, fes_map_moments, map2_fit_idc
%}
function [FES,info] = fes_map_aggregate(maps, servers, n, options)

if nargin < 4
    options = struct();
end
if ~isfield(options,'verbose')
    options.verbose = false;
end
if ~isfield(options,'grid') || isempty(options.grid)
    options.grid = fes_map_grid(n);
end

M = numel(maps);
if M < 1
    line_error(mfilename,'At least one station is required.');
end
if numel(servers) ~= M
    line_error(mfilename,'One server count per station is required.');
end

grid = unique([options.grid(:)' n]);
grid = grid(grid >= 1 & grid <= n);

FES = fes_map_levels(maps{1}, n, servers(1));
moments = zeros(4,n);
status = zeros(1,n);
for k = 1:n
    moments(:,k) = [map_moment(FES{k},1); map_moment(FES{k},2); map_moment(FES{k},3); map_idc(FES{k})];
end

for i = 2:M
    gmom = zeros(4,numel(grid));
    for g = 1:numel(grid)
        k = grid(g);
        [T0,T1] = fes_map_interdeparture(maps{i}, FES, k, [servers(i) 1]);
        [e1,e2,e3,~,idc] = fes_map_moments(T0,T1,options);
        gmom(:,g) = [e1;e2;e3;idc];
    end

    if numel(grid) < n
        moments = fes_map_interp(grid, gmom.', 1:n).';
    else
        moments = gmom;
    end

    newFES = cell(1,n);
    for k = 1:n
        [newFES{k},status(k)] = map2_fit_idc(moments(1,k), moments(2,k), moments(3,k), moments(4,k));
    end
    FES = newFES;

    if options.verbose
        fprintf('FES fold %d/%d: %d levels on a grid of %d, %d exact fits, %d fallbacks\n', ...
            i, M, n, numel(grid), sum(status==0), sum(status>0));
    end
end

info = struct();
info.throughput = 1./moments(1,:);
info.moments = moments;
info.status = status;
info.grid = grid;
end
