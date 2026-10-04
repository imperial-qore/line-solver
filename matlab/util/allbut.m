function y=allbut(y,ypos)
% y=ALLBUT(y,ypos)
% Returns all elements in y but the ones in ypos
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

y=y(setdiff(1:length(y),ypos));
end