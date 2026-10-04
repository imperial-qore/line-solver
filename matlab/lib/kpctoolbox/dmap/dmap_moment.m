function MOMENTS=dmap_moment(DMAP,ORDERS)
% MOMENTS=dmap_moment(DMAP,ORDERS) - Compute raw moments of
% inter-arrival times for a discrete-time MAP
%
%  Input:
%  DMAP: a D-MAP in the form of {D0,D1}
%  ORDERS: set of moment orders (1=>E[X], 2=>E[X^2], 3=>E[X^3])
%
%  Output:
%  MOMENTS: moments returned in the same order of ORDERS
%
%  NOTE: the previous formula factorial(i)*al*(I-D0)^-i*e returned neither the
%  raw nor the falling-factorial moment for order>=2 (e.g. it gave 4.447 where
%  the true E[X^2] is 2.958). The raw moments below match the JAR
%  implementation and a brute-force sum over P(X=k)=al*D0^(k-1)*D1*e.

D0=DMAP{1};
D1=DMAP{2};
N=size(D0,1);
e=ones(N,1);
P=inv(eye(N)-D0)*D1;
al=dtmc_solve(P);
A=inv(eye(N)-D0);
Ae=A*e;
m1=al*Ae;
for t=1:length(ORDERS)
    i=ORDERS(t);
    if any(isnan(D0(:)))
        MOMENTS(t)=NaN;
    else
        switch i
            case 1
                MOMENTS(t)=m1;
            case 2
                MOMENTS(t)=2*al*A*Ae - m1;
            case 3
                MOMENTS(t)=6*al*A*A*Ae - 6*al*A*Ae + m1;
            otherwise
                error('dmap_moment: raw moments of order > 3 not implemented');
        end
    end
end
end
