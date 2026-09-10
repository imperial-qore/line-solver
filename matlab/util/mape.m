function [M,nanM]=mape(approx, exact)
% M = MAPE(approx, exact)
% Return mean absolute percentage error of approx with respect to exact
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
M = mean(abs(1-approx(exact>0)./exact(exact>0)));
if nargout > 1
    nanM = nanmean(abs(1-approx(exact>0)./exact(exact>0)));
end
end