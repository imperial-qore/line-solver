function reset(self, resetState)
% RESET(RESETSTATE)
%
% If RESETSTATE is true, the model requires re-initialization
% of its state
if nargin == 1
    resetModel(self);
else
    resetModel(self, resetState);
end
self.hasState = false;
end