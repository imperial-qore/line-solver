function [ph, mu, phi, phases] = refreshServicePhases(self, statSet, classSet)
% [PH, MU, PHI, PHASES] = REFRESHSERVICEPHASES()
% Obtain information about phases of service and arrival processes.

% Copyright (c) 2012-2023, Imperial College London
% All rights reserved.

M = getNumberOfStations(self);
K = getNumberOfClasses(self);
if nargin<2
    statSet = 1:M;
    classSet = 1:K;
    mu = cell(M,1);
    for i=1:M
        mu{i,1} = cell(1,K);
    end
    phi = cell(M,1);
    for i=1:M
        phi{i,1} = cell(1,K);
    end
    phases = zeros(M,K);
elseif nargin==2
    classSet = 1:K;
    mu = cell(M,1);
    for i=1:M
        mu{i,1} = cell(1,K);
    end
    phi = cell(M,1);
    for i=1:M
        phi{i,1} = cell(1,K);
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
    for i=1:M
        mu{i,1} = cell(1,K);
    end
    phi = cell(M,1);
    for i=1:M
        phi{i,1} = cell(1,K);
    end
    phases = zeros(M,K);
end

stations = self.stations;
for i=statSet
    if i == self.getIndexSourceStation
        [~,mu_i,phi_i] = stations{i}.getSourceRates();
    else
        switch class(stations{i})
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
                [~,mu_i,phi_i] = stations{i}.getServiceRates();
        end
    end
    for r=classSet
        if ~isempty(mu_i{r})
        mu{i}{r} = mu_i{r};
        phi{i}{r} = phi_i{r};
        if isnan(mu_i{r}) % disabled
            phases(i,r) = 0;
        else
            phases(i,r) = length(mu_i{r});
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
[ph, phases] = refreshServiceRepresentations(self);
end