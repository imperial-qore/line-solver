function AN=getAvgArvRFromTput(sn, TN, TH)
% Compute average arrival rate at steady-state
% SN: network struct
% TN: average tput
% TH: tput handles

M = sn.nstations;
R = sn.nclasses;
if ~isempty(TH) && ~isempty(TN)
    AN = zeros(M,R);
    for i=1:M
        for j=1:M
            for k=1:R
                for r=1:R
                    AN(i,k) = AN(i,k) + TN(j,r)*sn.rt((j-1)*R+r, (i-1)*R+k);
                end
            end
        end
    end
else
    AN = [];
end
end