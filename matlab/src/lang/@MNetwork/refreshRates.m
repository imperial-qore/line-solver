function [rates, scv, hasRateChanged, hasSCVChanged] = refreshRates(self, statSet, classSet)
% [RATES, SCV] = REFRESHRATES()

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

hasRateChanged = false;
hasSCVChanged = false;
M = getNumberOfStations(self);
K = getNumberOfClasses(self);
if nargin<2
    statSet = 1:M;
    classSet = 1:K;
    rates = zeros(M,K);
    scv = nan(M,K);
    hasRateChanged = true;
    hasSCVChanged = true;
elseif nargin==2
    classSet = 1:K;
    rates = self.sn.rates;
    scv = self.sn.scv;
    rates_orig = self.sn.rates;
    scv_orig = self.sn.scv;
elseif nargin==3 % this is used only to update self.sn
    % we are only updating selected stations and classes so use the
    % existing ones for the others
    rates = self.sn.rates;
    scv = self.sn.scv;
    rates_orig = self.sn.rates;
    scv_orig = self.sn.scv;
end
hasOpenClasses = self.hasOpenClasses;

stations = self.stations;
% determine rates
for i=statSet
    station_i = stations{i};
    for r=classSet
        switch station_i.server.className
            case 'ServiceTunnel'
                % do nothing
                switch class(station_i)
                    case 'Source'
                        if isempty(station_i.input.sourceClasses{r})
                            rates(i,r) = NaN;
                            scv(i,r) = NaN;
                        else
                            distr = station_i.input.sourceClasses{r}{end};
                            rates(i,r) = distr.getRate();
                            scv(i,r) = distr.getSCV();
                        end
                    case 'Join'
                        rates(i,r) = Inf;
                        scv(i,r) = 0;
                    case 'Place'
                        rates(i,r) = NaN;
                        scv(i,r) = NaN;
                        %Draft SPN code:
                        %rates(i,r) = 1;
                        %scv(i,r) = 0;                        
                end
            otherwise
                if ~hasOpenClasses || i ~= self.getIndexSourceStation
                    if isempty(station_i.server.serviceProcess{r})
                        rates(i,r) = NaN;
                        scv(i,r) = NaN;
                    else
                        distr = station_i.server.serviceProcess{r}{end};
                        rates(i,r) = distr.getRate();
                        scv(i,r) = distr.getSCV();
                    end
                end
        end
    end
end

if ~hasRateChanged
    if any((abs(rates-rates_orig)>0))
        hasRateChanged = true;
    end
end

if ~hasSCVChanged
    if any((abs(scv-scv_orig)>0))
        hasSCVChanged = true;
    end
end

%if ~isempty(self.sn)
    if hasRateChanged
        self.sn.rates = rates;
    end
    if hasSCVChanged
        self.sn.scv = scv;
    end
%end

end
