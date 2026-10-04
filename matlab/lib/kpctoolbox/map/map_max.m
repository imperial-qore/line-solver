function M=map_max(A,B)
% M=map_max(A,B) - MAP of the maximum of two independent MAPs
%
%  Input:
%  A: a MAP in the form of {D0,D1}
%  B: a MAP in the form of {D0,D1}
%
%  Output:
%  M: a MAP in the form of {D0,D1} whose inter-arrival times are
%  distributed as max(X,Y), with X~A and Y~B drawn independently at each
%  arrival epoch from the respective embedded equilibrium distributions
%
%  The phase space is ordered as [(i,j) pairs, B-only phases, A-only
%  phases]: in the first block both A and B are still running, in the
%  second block A has already completed and B is awaited, in the third
%  block B has completed and A is awaited. An arrival is recorded when the
%  second of the two completes, i.e. only out of the last two blocks.
%
na=size(A{1},1);
nb=size(B{1},1);
a=-A{1}*ones(na,1);
b=-B{1}*ones(nb,1);
% pair state (i,j) has index (i-1)*nb+j, matching krons(A{1},B{1})
M0=[
    krons(A{1},B{1}),   kron(a,eye(nb)),    kron(eye(na),b)
    zeros(nb,na*nb),    B{1},               zeros(nb,na)
    zeros(na,na*nb),    zeros(na,nb),       A{1}
    ];
pie = [kron(map_pie(A),map_pie(B)),zeros(1,nb),zeros(1,na)];
% completion rate of the second of the two processes out of each phase; it
% is zero in the pair block, where only the first one can still complete
d = [zeros(na*nb,1); b; a];
M1 = d*pie;
M={M0,M1};
end
