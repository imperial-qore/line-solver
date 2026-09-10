function [ percentileRTs ] = returnPer( vector, Matrix, pers )

% mean
meanRT = sum(-vector/Matrix,2);

c = max(-diag(Matrix));
m = size(Matrix,2);
P_res = Matrix/c+eye(m);

M = sum(vector/(eye(m)-P_res),2);
a0 = sum(vector);
sum_a = a0;
k = 0;
vP = sum(P_res,2);
ak = [];
while abs(sum_a-M) >= 10^-10
    k = k+1;
    ak = [ak,vector*vP];
    sum_a = sum_a+vector*vP;
    vP = P_res*vP;
end
K1 = k;

%% percentiles
percentileRTs = zeros(length(pers),2);
percentileRTs(:,1) = pers;
for p = 1 : length(pers)
    
    if pers(p) < 1-sum(vector)
        percentileRTs(p,2) = 0;
    else
        MaxTime = 3*meanRT;
        flag = 0;
        while flag == 0
            
            %% for the MaxTime
            pM = expm(-c*MaxTime);
            F = pM*a0;
            for k = 1 : K1
                pM = c*MaxTime*pM/k;
                F = F+pM*ak(k);
            end
            per_achieved = 1-F;
            if per_achieved < pers(p)
                MaxTime = MaxTime+0.5*meanRT;
            else
                for t = MaxTime : -0.001 : 0
                    
                    pM = expm(-c*t);
                    F = pM*a0;
                    for k = 1 : K1
                        pM = c*t*pM/k;
                        F = F+pM*ak(k);
                    end
                    CDF_analysis = 1-F;
                    
                    if CDF_analysis < pers(p)
                        temp_percentileRT = t+0.001;
                        break
                    end
                end
                
                percentileRTs(p,2) = temp_percentileRT;
                flag = 1;
            end
        end
    end
    
end
