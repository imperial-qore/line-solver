function [rates,scv,mu,phi,phases,lst,proctypes] = refreshProcesses(self, statSet, classSet)
% [RATES,SCV, MU,PHI,PHASES,LT,PROCTYPES] = REFRESHPROCESSES(STATSET,CLASSSET)
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin==1
    [rates, scv, hasRateChanged, hasSCVChanged] = refreshRates(self);
    if hasSCVChanged
        proctypes = refreshProcessTypes(self);
        [~,mu,phi,phases] = refreshProcessPhases(self);
        lst = refreshLST(self);
    end
elseif nargin==2
    [rates, scv, hasRateChanged, hasSCVChanged] = self.refreshRates(statSet);
    if hasSCVChanged
        proctypes = refreshProcessTypes(self);
        if any(scv(statSet,:)) % with immediate phases this block is not needed
            [~,mu,phi,phases] = self.refreshProcessPhases(statSet);
            lst = self.refreshLST(statSet);
        end
    end
elseif nargin==3
    [rates, scv, hasRateChanged, hasSCVChanged] = self.refreshRates(statSet, classSet);
    if hasSCVChanged
        proctypes = refreshProcessTypes(self);
        if any(scv(statSet,classSet))  % with immediate phases this block is not needed
            [~,mu,phi,phases] = self.refreshProcessPhases(statSet, classSet);
            lst = self.refreshLST(statSet, classSet);
        end
    end
end

if isempty(self.sn.sched) || (hasRateChanged && any(self.sn.sched == SchedStrategy.SEPT | self.sn.sched == SchedStrategy.LEPT))
    refreshScheduling(self); % SEPT and LEPT may be affected
end
end
