function envModel = MAPQN2RENV(model, options)
% ENVMODEL = MAPQN2RENV(MODEL, OPTIONS)
% Random-environment image of a network with MAP/MMPP service or arrival
% processes.
%
% Retained name for the transformation now implemented by MAP2RENV, which
% generalizes it from a single MMPP2 service process to any number of MAP,
% MMPP2 or MMAP arrival and service processes of arbitrary phase order (the
% stage set is then the Cartesian product of the phase spaces). New code should
% call MAP2RENV directly.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 2
    options = [];
end
envModel = map2renv(model, options);
end
