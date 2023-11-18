function MMAP = mamap22_fit_gamma_bs(M1, M2, M3, GAMMA, P, B, S)
% Performs approximate fitting of a MMAP given the underlying MAP,
% the class probabilities (always fitted exactly), the backward moments,
% and the one-step class transition probabilities.
% Input
% - M1,M2,M3: moments of the inter-arrival times
% - GAMMA: auto-correlation decay rate of the inter-arrival times
% - P: class probabilities
% - B: first-order backward moments
% - S: one-step class transition probabilities
% Output
% - mmap: fitted MAMAP[2]

if size(B,1) == 1
    B = B';
end

% fit underlying AMAP(2)
[~,MAPS] = amap2_fit_gamma(M1, M2, M3, GAMMA);

% TODO: if it is a Poisson process then use second-order Poisson process
% with (h1 - h2 + h2 r1 = 0)

% fit marked Poisson process
if length(MAPS) == 1 && size(MAPS{1}{1},1) == 1
    MAP = MAPS{1};
    % number of classes
    m = length(P);
    % fit class probabilities
    MMAP = cell(1,2+m);
    MMAP{1} = MAP{1};
    MMAP{2} = MAP{2};
    for c = 1:m
        MMAP{2+c} = MMAP{2} .* P(c);
    end
    return;
end

% fit the MAMAP(2,m) using the underlying AMAP(2) form which produces the
% least error
MMAPS = cell(1,length(MAPS));
ERRORS = zeros(1,length(MAPS));
for j = 1:length(MAPS)
    [MMAPS{j},fB,fS] = mamap22_fit_bs_multiclass(MAPS{j}, P, B, S);
    ERRORS(j) = sum((fB(1)./B(1) - 1).^2) + sum((fS(1,1)./S(1,1) - 1).^2);
end
[~,BEST] = min(ERRORS);

MMAP = MMAPS{BEST};

end