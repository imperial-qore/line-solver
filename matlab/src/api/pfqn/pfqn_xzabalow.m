function [XN]=pfqn_xzabalow(L,N,Z)
    Ltot=sum(L);    
    XN=N/(Z+Ltot*N);
end