%{ @file weaklyconncomp.m
 %  @brief Finds weakly connected components in a directed graph
 %
 %  @author LINE Development Team
%}

%{
 % @brief Finds weakly connected components in a directed graph
 %
 % @details
 % Returns the weakly connected components in a graph G.
 %
 % @par Syntax:
 % @code
 % [S, C] = weaklyconncomp(G)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>G<td>Adjacency matrix of the graph
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>S<td>Number of connected components
 % <tr><td>C<td>Component assignment vector (C(i) is the component index of node i)
 % </table>
%}
function [S,C] = weaklyconncomp(G)
  G(isnan(G)) = 1;
  [p,~,r] = dmperm(G'+speye(size(G)));
  S = numel(r)-1;
  C = cumsum(full(sparse(1,r(1:end-1),1,1,size(G,1))));
  C(p) = C;
end