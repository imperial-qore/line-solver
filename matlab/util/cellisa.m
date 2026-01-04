function v=cellisa(c,className)
% s=CELLISA(c,class)
% Check if cell elements belong to class
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

v=cellfun(@(x) isa(x,className),c);
end