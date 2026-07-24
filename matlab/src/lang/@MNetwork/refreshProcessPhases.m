function [ph, mu, phi, phases] = refreshProcessPhases(self, statSet, classSet)
% [PH, MU, PHI, PHASES] = REFRESHPROCESSPHASES()
% Obtain information about phases of service and arrival processes.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

M = getNumberOfStations(self);
K = getNumberOfClasses(self);
if nargin<2
    statSet = 1:M;
    classSet = 1:K;
    mu = cell(M,1);
    for ist=1:M
        mu{ist,1} = cell(1,K);
    end
    phi = cell(M,1);
    for ist=1:M
        phi{ist,1} = cell(1,K);
    end
    phases = zeros(M,K);
elseif nargin==2
    classSet = 1:K;
    mu = cell(M,1);
    for ist=1:M
        mu{ist,1} = cell(1,K);
    end
    phi = cell(M,1);
    for ist=1:M
        phi{ist,1} = cell(1,K);
    end
    phases = zeros(M,K);
elseif nargin==3 && isfield(self.sn, 'mu') && isfield(self.sn, 'phi') && isfield(self.sn, 'phases')
    % we are only updating selected stations and classes so use the
    % existing ones for the others
    mu = self.sn.mu;
    phi = self.sn.phi;
    phases = self.sn.phases;
else
    mu = cell(M,1);
    for ist=1:M
        mu{ist,1} = cell(1,K);
    end
    phi = cell(M,1);
    for ist=1:M
        phi{ist,1} = cell(1,K);
    end
    phases = zeros(M,K);
end

stations = self.stations;
for ist=statSet
    if ist == self.getIndexSourceStation
        [~,mu_i,phi_i] = stations{ist}.getSourceRates();
    else
        switch class(stations{ist})
            case 'Fork'
                mu_i = cell(1,K);
                phi_i = cell(1,K);
                for r=1:K
                    mu_i{r} = NaN;
                    phi_i{r} = NaN;
                end
            case 'Join'
                mu_i = cell(1,K);
                phi_i = cell(1,K);
                for r=1:K
                    mu_i{r} = NaN;
                    phi_i{r} = NaN;
                end
            otherwise
                [~,mu_i,phi_i] = stations{ist}.getServiceRates();
        end
    end
    for r=classSet
        if ~isempty(mu_i{r})
        mu{ist}{r} = mu_i{r};
        phi{ist}{r} = phi_i{r};
        if isnan(mu_i{r}) % disabled
            phases(ist,r) = 0;
        else
            phases(ist,r) = length(mu_i{r});
        end
    end
end

if ~isempty(self.sn) %&& isprop(self.sn,'mu')
    self.sn.mu = mu;
    self.sn.phi = phi;
    self.sn.phases = phases;
    self.sn.phasessz = max(self.sn.phases,ones(size(self.sn.phases)));
    self.sn.phaseshift = [zeros(size(phases,1),1),cumsum(self.sn.phasessz,2)];
end
[ph, phases] = refreshProcessRepresentations(self);
end