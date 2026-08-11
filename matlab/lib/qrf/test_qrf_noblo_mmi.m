%% QRF BAS test from code release
%function [UUB,ULB,ID]=runquadraticfcrbas(MAP,N,P,F,MM)
addpath C:\Users\csg\Dropbox\experiments\09_mapqnfc\
mapstr='exp';
P=[0,1; 1,0];
%MAPs={map_exponential(1); map_exponential(1)};
MAPs={map_rand(2); map_exponential(1)};
M = length(MAPs);
N=3; %(F1+1)
F1=2; F=[min([N,F1]),N,N,N];
[XN,QNctmc,UNctmc,EN,p,p1,p1c,p2,MAPQN,MM]=mapqnfc_ezsolve_bas(MAPs,N,P,ones(size(MAPs,1),N),F); 
QNctmc, UNctmc
M=length(MAPs);
for i=1:M
    K(i)=size(MAPs{i}{1},1);
end
r = zeros(M);
for i=1:M
    for j=1:M
        r(i,j)=P(i,j);
    end
end
for i=1:M
    for h=1:size(MAPs{i}{1},1)
        for k=1:size(MAPs{i}{1},1)
            mu(i,h,k)=MAPs{i}{2}(h,k);
        end
    end
end
for i=1:M
    for h=1:size(MAPs{i}{1},1)
        for k=1:size(MAPs{i}{1},1)
            if h==k
                v(i,k,h)=0; 
            else
                v(i,k,h)=MAPs{i}{1}(h,k);
            end
        end
    end
end
MR=size(MM,1);
for m=1:size(MM,1)
    for j=1:M
        BB(m,j)=any(find(MM(m,:)==j));
    end
end
ZZ=[];
for m=1:size(MM,1)
    ZZ(m)=nnz(MM(m,:));
end
ZM=max(ZZ);
MM1=[];
for m=1:size(MM,1)
    if nnz(MM(m,:)) == max(ZZ)-1
        for j=1:M
            M1 = MM(m,:);
            M1(end) = j;
            m1=matchrow(MM,M1);
            MM1(m,j)=m1;
        end
    end
end
%%
% [UNqrf,QNqrf] = qrf_bas_mmi(1,M,MR,MM,MM1,ZZ,ZM,BB,K,F,N,mu,v,r);
% UNex=UNctmc(:)'
% UNap=UNqrf
% QNex=QNctmc(:)'
% QNap=QNqrf/sum(QNqrf)*N
% [UNqrf,QNqrf] = qrf_bas_mem(1,M,MR,MM,MM1,ZZ,ZM,BB,K,F,N,mu,v,r);
% UNex=UNctmc(:)'
% UNap=UNqrf
% QNex=QNctmc(:)'
% QNap=QNqrf/sum(QNqrf)*N
[UNqrf,QNqrf] = qrf_bas_mmi_simple(1,M,MR,BB,K,F,N,mu,v,r);
UNex=sum(UNctmc,2)'
UNap=UNqrf
QNex=sum(QNctmc,2)'
QNap=QNqrf/sum(QNqrf)*N