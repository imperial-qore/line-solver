function [M,nanM]=maxpe(approx, exact)
% M = MAXPE(approx, exact)
% Returns maximum absolute percentage error of approx with respect to exact
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
M = max(abs(1-approx(exact>0)./exact(exact>0)));
if nargout > 1
    nanM = nanmax(abs(1-approx(exact>0)./exact(exact>0)));
end
end