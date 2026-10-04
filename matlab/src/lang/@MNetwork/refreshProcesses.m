function [rates,scv,mu,phi,phases,lst,proctypes] = refreshProcesses(self, statSet, classSet)
% [RATES,SCV, MU,PHI,PHASES,LT,PROCTYPES] = REFRESHPROCESSES(STATSET,CLASSSET)
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin==1
    [rates, scv, hasRateChanged, hasSCVChanged] = refreshRates(self);
    % Always refresh process types since distribution type can change
    % even when rate/SCV stay the same (e.g., Replayer -> APH with same SCV)
    proctypes = refreshProcessTypes(self);
    if hasSCVChanged
        [~,mu,phi,phases] = refreshProcessPhases(self);
        lst = refreshLST(self);
    end
elseif nargin==2
    [rates, scv, hasRateChanged, hasSCVChanged] = self.refreshRates(statSet);
    % Always refresh process types since distribution type can change
    proctypes = refreshProcessTypes(self);
    if hasSCVChanged
        if any(scv(statSet,:)) % with immediate phases this block is not needed
            [~,mu,phi,phases] = self.refreshProcessPhases(statSet);
            lst = self.refreshLST(statSet);
        end
    end
elseif nargin==3
    [rates, scv, hasRateChanged, hasSCVChanged] = self.refreshRates(statSet, classSet);
    % Always refresh process types since distribution type can change
    proctypes = refreshProcessTypes(self);
    if hasSCVChanged
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
