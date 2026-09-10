function MAPH = maph2m_fit(M1, M2, M3, P, B)
% Computes the second-order MAPH[m] fitting the given ordinary moments 
% (of order up to three), the class probabilities (always fitted exactly)
% and the backward moments,
%
% Input
% - M1,M2,M3: moments of the inter-arrival times
% - P: class probabilities
% - B: first-order backward moments
% Output
% - MAPH: fitted second-order MAPH[m]

if size(B,1) == 1
    B = B';
end

% fit underlying AMAP(2)
[~,APHS] = aph2_fit(M1, M2, M3);

% fit the MAMAP(2,m) using the underlying AMAP(2) form which produces the
% least error
MAPHS = cell(1,length(APHS));
ERRORS = zeros(1,length(APHS));
for j = 1:length(APHS)    
    [MAPHS{j},fB] = maph2m_fit_multiclass(APHS{j}, P, B);
    ERRORS(j) = sum((fB./B - 1).^2);
end
[~,BEST] = min(ERRORS);

MAPH = MAPHS{BEST};

end