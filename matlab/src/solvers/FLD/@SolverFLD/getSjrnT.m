function RD = getSjrnT(self, R)
% RD = GETSJRNT(R)
% Alias for getCdfRespT. Returns cumulative distribution functions of sojourn times.

if nargin < 2
    RD = self.getCdfRespT;
else
    RD = self.getCdfRespT(R);
end
end
