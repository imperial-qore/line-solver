function retval = containsstr(varargin)
% R = CONTAINSSTR(STR)
% Determine if pattern is in string
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if length(varargin) == 2
    retval = builtin('contains',varargin{1},varargin{2});
else
    retval = builtin('contains',varargin{1},varargin{2},varargin{3},varargin{4});
end
end
