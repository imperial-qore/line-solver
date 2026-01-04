%{
%{
 % @file pfqn_grnmol.m
 % @brief Normalizing constant using Grundmann-Moeller quadrature.
%}
%}

%{
%{
 % @brief Normalizing constant using Grundmann-Moeller quadrature.
 % @fn pfqn_grnmol(L, N)
 % @param L Service demand matrix.
 % @param N Population vector.
 % @return G Normalizing constant.
%}
%}
function G=pfqn_grnmol(L,N)
[M,R]=size(L);
G=0;
S=ceil(sum(N)-1)/2;
H=zeros(1+S,1);
for i=0:S
    c(1+i)=2*(S-i)+M;
    w(1+i)=2^-(2*S)*(-1)^i*c(1+i)^(2*S+1)/factorial(i)/factorial(i+c(1+i));
    [s,bvec,SD,D]=sprod(M,S-i); bvec=bvec';
    while bvec(1)>=0
        H(1+i) = H(1+i) + prod((((2*bvec+1)/c(1+i))*L).^N);
        [s,bvec]=sprod(s,SD,D); bvec=bvec';
    end
    G = G + w(1+i)*H(1+i);
end
G=G*factorial(sum(N)+M-1)/prod(factorial(N));
end