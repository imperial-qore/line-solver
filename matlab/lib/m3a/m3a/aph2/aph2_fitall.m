function APHS = aph2_fitall(M1,M2,M3)
% Fits an acyclic phase type distribution with the given moments.
% The result may be unfeasible.

degentol = 1e-8;

SCV = (M2-M1^2)/M1^2;
% lower bound of M3 for SCV <= 1
M3lb = 3*M1^3*(3*SCV-1+sqrt(2)*(1-SCV)^(3/2));

if SCV <= 1 && abs(M3 - M3lb) < degentol
    % when M3 is close to the lower bound, the sqrt argument is zero
    tmp0 = 0;
else
    tmp0 = M3^2/9 + ((8*M1^3)/3 - 2*M2*M1)*M3 - 3*M1^2*M2^2 + 2*M2^3;
    if tmp0 < 0
        % infeasible: square root of negative element
        % APHS = {};
        % added by GC
        APHS = {aph_fit(M1,M2,M3,2)};
        return;
    end
end
tmp1 = 3*sqrt(tmp0);
tmp2 = M3 - 3*M1*M2;
tmp3 = (6*M2 - 12*M1^2);

% maximum number of solutions
if tmp0 == 0
    % the diagonal elements of D0 are identical
	n = 1;
else
    n = 2;
end

% solution parameters
h1v = zeros(n,1);
h2v = zeros(n,1);
r1v = zeros(n,1);

if n == 1
    h2v = tmp2/tmp3;
    h1v = h2v;
else
    h2v(1) = (tmp2 + tmp1)/tmp3;
    h2v(2) = (tmp2 - tmp1)/tmp3;
    h1v(2) = h2v(1);
    h1v(1) = h2v(2);
end

for j = 1:n
    h1 = h1v(j);
    h2 = h2v(j);
    r1v(j) = (M1 - h1)/h2;
end

APHS = cell(1,n);
idx = 1;
for j = 1:n
    h1 = h1v(j);
    h2 = h2v(j);
    r1 = r1v(j);
    if h1 > 0 && h2 > 0 && r1 >= -degentol && r1 <= (1+degentol)
        % feasible (or almost feasible) solution!
        r1 = max(min(r1,1),0);
        APHS{idx} = aph2_assemble(h1, h2, r1);
        idx = idx + 1;
    end
end

% returns 0, 1 or 2 feasible solutions
APHS = APHS(1:(idx-1));
% added by GC
if isempty(APHS)
    APHS = {aph_fit(M1,M2,M3,2)};    
end
end