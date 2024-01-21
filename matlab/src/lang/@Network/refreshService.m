function [rates,scv,mu,phi,phases,lst,proctypes] = refreshService(self, statSet, classSet)
% [RATES,SCV, MU,PHI,PHASES,LT,PROCTYPES] = REFRESHSERVICE(STATSET,CLASSSET)
% Copyright (c) 2012-2024, Imperial College London
% All rights reserved.

if nargin==1
    [rates, scv, hasRateChanged, hasSCVChanged] = refreshRates(self);
    if hasSCVChanged
        proctypes = refreshServiceTypes(self);
        [~,mu,phi,phases] = refreshServicePhases(self);
        lst = refreshLST(self);
    end
elseif nargin==2
    [rates, scv, hasRateChanged, hasSCVChanged] = self.refreshRates(statSet);
    if hasSCVChanged
        proctypes = refreshServiceTypes(self);
        if any(scv(statSet,:)) % with immediate phases this block is not needed
            [~,mu,phi,phases] = self.refreshServicePhases(statSet);
            lst = self.refreshLST(statSet);
        end
    end
elseif nargin==3
    [rates, scv, hasRateChanged, hasSCVChanged] = self.refreshRates(statSet, classSet);
    if hasSCVChanged
        proctypes = refreshServiceTypes(self);
        if any(scv(statSet,classSet))  % with immediate phases this block is not needed
            [~,mu,phi,phases] = self.refreshServicePhases(statSet, classSet);
            lst = self.refreshLST(statSet, classSet);
        end
    end
end

if isempty(self.sn.schedid) || (hasRateChanged && any(self.sn.schedid == SchedStrategy.ID_SEPT | self.sn.schedid == SchedStrategy.ID_LEPT))
    refreshScheduling(self); % SEPT and LEPT may be affected
end
end
