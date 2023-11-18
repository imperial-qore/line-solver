function [XN]=pfqn_xzabaup(L,N,Z)    
    XN=min([1/max(L) N/(sum(L)+Z)]);
end