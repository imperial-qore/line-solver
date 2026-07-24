function RD = sjrnT(self, R)
% RD = SJRNT(R)
% Alias for getSjrnT. Lowercase Kotlin-style wrapper.

if nargin < 2
    RD = self.getSjrnT;
else
    RD = self.getSjrnT(R);
end
end
