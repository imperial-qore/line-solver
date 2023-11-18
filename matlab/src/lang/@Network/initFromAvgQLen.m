function initFromAvgQLen(self, AvgQLen)
% INITFROMAVGQLEN(AVGQLEN)

n = round(AvgQLen);
njobs = sum(n,1);
% we now address the problem that round([0.5,0.5]) = [1,1] so
% different from the total initial population
for r=1:size(AvgQLen,2)
    if njobs(r) > sum(AvgQLen,1) % error at most by 1
        i = maxpos(n(:,r));
        n(i,r) = n(i,r) - 1;
        njobs = sum(n,1)';
    end
end
try
    self.initFromMarginal(n);
catch
    self.initDefault;
end
end