%{ @file fes_map_grid.m
 %  @brief Population grid on which a MAP flow-equivalent server is fitted
 %
 %  @author LINE Development Team
%}

%{
 % @brief Returns the population levels at which the inter-departure MAP is
 % evaluated
 %
 % @details
 % Fitting one MAP per population is wasteful because the processes of
 % neighbouring populations are similar. Section 5.2.2 of Casale, Mi,
 % Cherkasova and Smirni, IEEE Trans. Soft. Eng. 37(5), 2011, evaluates the
 % first ten populations and ten further equispaced points, which is what
 % this function returns. Populations up to the grid size are returned in
 % full, so no approximation is introduced on small models.
 %
 % @par Syntax:
 % @code
 % grid = fes_map_grid(n)
 % grid = fes_map_grid(n, nhead, ntail)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>n<td>Largest population
 % <tr><td>nhead<td>(Optional) leading populations kept in full, default 10
 % <tr><td>ntail<td>(Optional) equispaced points above them, default 10
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>grid<td>Sorted vector of populations to evaluate
 % </table>
 %
 % @see fes_map_aggregate
%}
function grid = fes_map_grid(n, nhead, ntail)

if nargin < 2
    nhead = 10;
end
if nargin < 3
    ntail = 10;
end

if n <= nhead + ntail
    grid = 1:n;
    return
end

grid = unique([1:nhead, round(linspace(nhead+1, n, ntail))]);
end
