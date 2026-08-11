function varargout = avg(self,varargin)
% VARARGOUT = AVG(SELF,VARARGIN)
% Alias for getAvg
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

[varargout{1:nargout}] = self.getAvg(varargin{:});
end
