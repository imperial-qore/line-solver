function AMAPS = amap2_fitall_gamma(M1,M2,M3,GAMMA)
% Finds the equivalent AMAP(2)s (up to 4) fitting the given set of 
% characteristics. If the characteristics are infeasible, an empty cell
% array is returned (no approximate fitting is performed).
% Input:
% - M1,M2,M3: moments of the inter-arrival times
% - GAMMA: auto-correlation decay rate of the inter-arrival times
% Output:
% - AMAPS: a cell array of feasible AMAP(2) processes

degentol = 1e-8;
r12tol = 1e-6;

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
        AMAPS = {};
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

if n == 1
    h2v = tmp2/tmp3;
    h1v = h2v;
else
    h2v(1) = (tmp2 + tmp1)/tmp3;
    h2v(2) = (tmp2 - tmp1)/tmp3;
    h1v(2) = h2v(1);
    h1v(1) = h2v(2);
end

% check for infeasibility
if min(h2v) <= 0
    AMAPS = {};
    return;
end

if GAMMA > 0
    AMAPS = cell(1,n*2);
else
    AMAPS = cell(1,n);
end

idx = 1;
for j = 1:n
    
    h1 = h1v(j);
    h2 = h2v(j);
    
    if (GAMMA >= 0)

        % FIRST canonical form
    
        z = M1^2*GAMMA^2 + (2*M1*h1 + 2*M1*h2 - 4*h1*h2 - 2*M1^2)*GAMMA + M1^2 - 2*M1*h1 - 2*M1*h2 + h1^2 + 2*h1*h2 + h2^2;
        
        if abs(z) < degentol
            % one solution, argument of square root is practically zero
            r2 = (h1 - M1 + h2 + GAMMA*M1)/(2*h1);
            r1 = (M1 - h1 - M1*r2 + h1*r2)/(h2 - M1*r2);
            % construct MAP if feasible
            if feasible(r1,r2)
                r1 = fix(r1);
                r2 = fix(r2);
                AMAPS{idx} = amap2_assemble(h1, h2, r1, r2, 1);
                idx = idx + 1;
            end
        elseif z > 0
            % two possible values of r2
            r2v = zeros(2,1);
            r2v(1) = (h1 - M1 + h2 - sqrt(z) + GAMMA*M1)/(2*h1);
            r2v(2) = (h1 - M1 + h2 + sqrt(z) + GAMMA*M1)/(2*h1);
            % two possible values of r1
            r1v = zeros(2,1);
            for i = 1:2
                r2 = r2v(i);
                r1v(i) = (M1 - h1 - M1*r2 + h1*r2)/(h2 - M1*r2);
            end
            % assemble
            for i = 1:2
                r1 = r1v(i);
                r2 = r2v(i);
                % construct MAP if feasible
                if feasible(r1,r2)
                    r1 = fix(r1);
                    r2 = fix(r2);
                    AMAPS{idx} = amap2_assemble(h1, h2, r1, r2, 1);
                    idx = idx + 1;
                end
            end
        end
        
    else
    
        % SECOND canonical form

        r2 = (h1 - M1 + h2 + GAMMA*M1)/h1;
        r1 = (r2 + (h1 + h2 - h1*r2)/M1 - 2)/(r2 - 1);
        % construct MAP if feasible
        if feasible(r1,r2)
            r2 = fix(r2);
            r1 = fix(r1);
            AMAPS{idx} = amap2_assemble(h1, h2, r1, r2, 2);
            idx = idx + 1;
        end
        
    end
end

% return feasible solutions
AMAPS = AMAPS(1:(idx-1));

    function q = fix(q)
        q = max(min(q,1),0);
    end

    function feas = feasible(r1, r2)
        feas = isreal(r1) && isreal(r2) && ...
               r1 >= -r12tol && r1 <= 1+r12tol && ...
               r2 >= -r12tol && r2 <= 1+r12tol;
    end

end
