function RD = getSjrnT(self, R)
% RD = GETSJRNT(R)
% Alias for getCdfRespT. Returns cumulative distribution functions of sojourn times.


% lang='cpp' needs no arm here: this is an alias, and getCdfRespT serves it
% from line-cli under every lang.

if nargin < 2
    RD = self.getCdfRespT;
else
    RD = self.getCdfRespT(R);
end
end
