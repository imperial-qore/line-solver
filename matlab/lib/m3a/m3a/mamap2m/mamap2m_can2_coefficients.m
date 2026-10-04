function [E,V,Z] = mamap2m_can2_coefficients(h1,h2,r1,r2)
% Returns the coefficients used in the direct and inverse formulas for
% fitting a MAMAP(2,m) in first canonical form (gamma < 0).
% Input:
% - h1,h2,r1,r2: the parameters of the underlying AMAP(2) with gamma < 0
% Output:
% - E: coefficients of p1, p11, F11 and B11
% - V: coefficients to fit {p1,p11,F11} or {p1,p11,B11}
% - Z: denominators

if issym(h1) || issym(h2) || issym(r1) || issym(r2)
    E = sym(zeros(14,1));
    V = sym(zeros(12,1));
    Z = sym(zeros(3,1));
else
    E = zeros(15,1);
    V = zeros(12,1);
    Z = zeros(3,1);
end

E(1) = 1 - 1/(r2*(r1 - 1) - r1 + 2);
E(2) = -(r2 - 1)/(r1*(r2 - 1) - r2 + 2);
E(3) = r2/(r1*(r2 - 1) - r2 + 2);
E(4) = (r2 - 2)/(r1*(r2 - 1) - r2 + 2) - r2 + 2;
E(5) = r2 - r2/(r1*(r2 - 1) - r2 + 2);
E(6) = -(r1*(r2 - 1)^2)/(r1 + r2 - r1*r2 - 2);
E(7) = - r2 - (r2*(2*r2 - 3))/(r1*(r2 - 1) - r2 + 2);
E(8) = r2^2/(r2*(r1 - 1) - r1 + 2);
E(9) = h1 - h1/(r2*(r1 - 1) - r1 + 2);
E(10) = h1*(r2 - 1) - ((r2 - 1)*(2*h1 + h2 - h1*r2))/(r1*(r2 - 1) - r2 + 2);
E(11) = (r2*(2*h1 + h2 - h1*r2))/(r1*(r2 - 1) - r2 + 2) - h1*r2;
E(12) = h2 - h2/(r2*(r1 - 1) - r1 + 2);
E(13) = ((h1 + h2*r1)*(r2 - 1))/(r1 + r2 - r1*r2 - 2);
E(14) = (h2*r2)/(r1*(r2 - 1) - r2 + 2);
V(1) = -(r1 + r2 - r1*r2 - 2)^2;
V(2) = -(r1 + r2 - r1*r2 - 2)*(2*h1 + 2*h2 - h1*r2 - h2*r2 + h2*r1*r2);
V(3) = h2*(2*h1 - h1*r2 + h2*r1)*(r1 + r2 - r1*r2 - 2);
V(4) = (r2 - 1)*(h1 - h2 + h2*r1)^2;
V(5) = (h1 - h2 + h2*r1)*(2*r2 - r1*r2 + r1*r2^2 - r2^2);
V(6) = -(h1*r2 + h2*r2 - h1*r2^2)*(h1 - h2 + h2*r1);
V(7) = -(r1 + r2 - r1*r2 - 2)^2;
V(8) = -(r1 + r2 - r1*r2 - 2)*(2*h1 + 2*h2 - h1*r2 - h2*r2 + h1*r1*r2^2 - h1*r1*r2);
V(9) = h1*(r1 + r2 - r1*r2 - 2)*(2*h2 + h1*r1 - h2*r2 + h1*r1*r2^2 - 2*h1*r1*r2);
V(10) = (r2 - 1)*(h1 - h2 - h1*r1 + h1*r1*r2)^2;
V(11) = -r2*(h1 - h2 - h1*r1 + h1*r1*r2)*(r1 + r2 - r1*r2 - 2);
V(12) = -r2*(h1 + h2 - h1*r2)*(h1 - h2 - h1*r1 + h1*r1*r2);
Z(1) = E(10)*E(12)*E(3)-E(10)*E(14)*E(1)-E(11)*E(12)*E(2)+E(11)*E(13)*E(1)-E(13)*E(3)*E(9)+E(14)*E(2)*E(9);
Z(2) = E(12) * E(2) - E(13) * E(1);
Z(3) = E(10)*E(1) - E(2)*E(9);

end