function MMAP = mamap2m_fit_gamma_fb(M1, M2, M3, GAMMA, P, F, B)
% Computes the second-order MAMAP[m] fitting the given ordinary moments 
% (of order up to three), the autocorrelation decay rate,
% the class probabilities (always fitted exactly), the forward moments,
% and the backward moments.
%
% Input
% - M1,M2,M3: moments of the inter-arrival times
% - GAMMA: auto-correlation decay rate of the inter-arrival times
% - P: class probabilities
% - F: first-order forward moments
% - B: first-order backward moments
% Output
% - mmap: fitted MAMAP[m]

if size(F,1) == 1
    F = F';
end

if size(B,1) == 1
   B = B'; 
end

% fit underlying AMAP(2)
% up to four equivalent AMAP(2) representations may be found
[~,MAPS] = amap2_fit_gamma(M1, M2, M3, GAMMA);

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

% use the underlying AMAP(2) form which produces the least error
MMAPS = cell(1,length(MAPS));
ERRORS = zeros(1,length(MAPS));
for j = 1:length(MAPS)
    [MMAPS{j},fF,fB] = mamap2m_fit_fb_multiclass(MAPS{j}, P, F, B);
    ERRORS(j) = sum((fF./F - 1).^2) + sum((fB./B - 1).^2);
end
[~,BEST] = min(ERRORS);
%fprintf('Fitting MAMAP(2,m) F+B: minimum error is %e\n', ERRORS(BEST));

MMAP = MMAPS{BEST};

end