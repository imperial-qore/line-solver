function MMAP = mmap_mixture(PHs, P2)
% Fits a MMAP with m classes using a mixture of m^2 PH-distributions.
% PH{i,j} represents the IAT distribution conditioned
% on the fact that the last arrival was of class i and the next arrival is
% of class j. Currently, the PH distributions are of the second order,
% hence the cross moments of order 1, 2 and 3 are required.
%
% INPUT
% - PHs: PH distributions
% - P2: two-step class transition probabilities
% OUTPUT
% - MMAP: the fitted MMAP
% - PHs: the fitted PH-distributions for each transition

% number of classes
m = size(PHs,1);

MMAP = cell(1,2+m);
MMAP{1} = zeros(2*m^2, 2*m^2);
MMAP{2} = zeros(2*m^2, 2*m^2);
for i = 1:m
    MMAP{2+i} = zeros(2*m^2, 2*m^2);
end

k = 1;
for i = 1:m
    for j = 1:m
        first = (k-1)*2+1;
        last = k*2;
        MMAP{1}(first:last,first:last) = PHs{i,j}{1};
        k = k + 1;
    end
end

k1 = 1;
for i1 = 1:m
    for j1 = 1:m
        k2 = 1;
        for i2 = 1:m
            for j2 = 1:m
                firstr = (k1-1)*2+1;
                lastr  = k1*2;
                firstc = (k2-1)*2+1;
                lastc  = k2*2;
                if j1 == i2
                    Q = (-PHs{i1,j1}{1}) * ones(2,1) * map_pie(PHs{i2,j2});
                    p = P2(i1,i2);
                else
                    Q = zeros(2,2);
                    p = 0;
                end                
                MMAP{2}(firstr:lastr,firstc:lastc) = p * Q;
                MMAP{2+j1}(firstr:lastr,firstc:lastc) = p * Q;
                k2 = k2 + 1;
            end
        end
        k1 = k1 +1;
    end
end

% normalize diagonal
for i = 1:(m^2)
    MMAP{1}(i,i) = 0;
end
for i = 1:(m^2)
    MMAP{1}(i,i) = -sum(MMAP{1}(i,1:end)+MMAP{2}(i,:));
end


end % end function