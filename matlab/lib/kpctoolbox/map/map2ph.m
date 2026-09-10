function [alpha,T,PHR]=map2ph(MAPIN)
% [MAPOUT,A,pi]=map2ph(MAPIN) - Returns a PH distribution
%
%  Input:
%  MAPIN: a MAP in the form of {D0,D1}
%
%  Output:
%  alpha: entry probability vector of the PH distribution
%  T: subgenerator of the PH distribution
%  PHR: PH-renewal process with distribution (pi,A) in (D0,D1) notation

T=MAPIN{1};
alpha=map_pie(MAPIN);
PHR{2}=MAPIN{2}*ones(length(MAPIN{2}),1)*alpha;
end

