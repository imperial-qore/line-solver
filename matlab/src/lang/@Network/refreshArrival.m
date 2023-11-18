function [rates, mu, phi, phases] = refreshArrival(self) % LINE treats arrival distributions as service distributions of the Source object
% [RATES, MU, PHI, PHASES] = REFRESHARRIVAL() % LINE TREATS ARRIVAL DISTRIBUTIONS AS SERVICE DISTRIBUTIONS OF THE SOURCE OBJECT

[rates, mu, phi, phases, ~, ~] = refreshService(self);
end