function [S,C] = weaklyconncomp(G)
% [S,C] = WEAKLYCONNCOMP(G)
% Return the weakly connected components in a graph from an undirected
% incidence matrix G
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
  [p,~,r] = dmperm(G'+speye(size(G))); %dmperm is a btf computation
  S = numel(r)-1;
  C = cumsum(full(sparse(1,r(1:end-1),1,1,size(G,1))));
  C(p) = C;
end