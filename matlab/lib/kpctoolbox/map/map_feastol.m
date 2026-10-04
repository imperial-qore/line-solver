function FEASTOL=map_feastol()
% FEASTOL=map_feastol() - feasibility tolerance EXPONENT (default 8)
%
%  Output:
%  FEASTOL: the exponent k of the toolbox feasibility tolerance 10^-k, so
%  the tolerance itself is 10^-8. Callers use it as 10^(-map_feastol), or
%  compare a tolerance magnitude against it directly as in map_isfeasible.
%  It is NOT the tolerance: map_feastol()==8, not 1e-8.
%

FEASTOL=8;
end