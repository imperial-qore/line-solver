function [MMAP,PHs] = mmap_mixture_fit(P2,M1,M2,M3)
% Fits a MMAP with m classes using a mixture of m^2 PH-distributions.
% Each PH distribution represents the probability distribution conditioned
% on the fact that the last arrival was of class i and the next arrival is
% of class j. Currently, the PH distributions are of the second order,
% hence the cross moments of order 1, 2 and 3 are required.
%
% INPUT
% - P2: two-step class transition probabilities (use mtrace_sigma2)
% - M1: matrix of first-order cross moments
% - M2: matrix of first-order cross moments
% - M3: matrix of first-order cross moments
% OUTPUT
% - MMAP: the fitted MMAP
% - PHs: the fitted PH-distributions for each transition

% number of classes
m = size(M1,1);

PHs = cell(m,m);
for i = 1:m
    for j = 1:m
        % approximate fitting of second and third moments
        m1 = M1(i,j);
        PHs{i,j} = aph2_fit(m1, M2(i,j), M3(i,j));
        % compare theoretical and fitted moments
        % fprintf('PH{%d,%d}: M1 = %.3e, M2 = %.3e (%.3e), M3 = %.3e (%.3e)\n', i, j, m1, m2, M2(i,j), m3, M3(i,j));
    end
end

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
                    p = P2(i1,i2,j2)/sum(P2(i1,i2,:));
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

M = {M1, M2, M3};
FM = cell(3,1);
FM{1} = mmap_cross_moment(MMAP,1);
FM{2} = mmap_cross_moment(MMAP,2);
FM{3} = mmap_cross_moment(MMAP,3);

% for i = 1:m
%     for j = 1:m
%         %fprintf('MMAP CROSS(%d,%d): ', i, j);
%         for h = 1:3
%           %  fprintf('M%d = %.3e (%.3e)', h, FM{h}(i,j), M{h}(i,j));
%             if h < 3
%                 fprintf(', ');
%             end
%         end
%         fprintf('\n');
%     end
% end

FP2 = mmap_sigma2(MMAP);
% for i = 1:m
%     for j = 1:m
%         %fprintf('MMAP p(%d,%d,:): ', i, j);
%         for h = 1:m
%          %   fprintf('p(%d,%d,%d) = %.3e (%.3e)', i, j, h, P2(i,j,h), FP2(i,j,h));
%             if h < m
%                 fprintf(', ');
%             end
%         end
%         fprintf('\n');
%     end
% end

end % end function