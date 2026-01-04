function [G,U,Y] = mamap2m_can1_coefficients(h1,h2,r1,r2)
% Returns the coefficients used in the direct and inverse formulas for
% fitting a MAMAP(2,m) in first canonical form (gamma > 0).
% Input:
% - h1,h2,r1,r2: the parameters of the underlying AMAP(2) with gamma > 0
% Output:
% - G: coefficients of p1, p11, F11 and B11
% - U: coefficients to fit {p1,p11,F11} or {p1,p11,B11}
% - Y: denominators

if issym(h1) || issym(h2) || issym(r1) || issym(r2)
    G = sym(zeros(15,1));
    U = sym(zeros(12,1));
    Y = sym(zeros(3,1));
else
    G = zeros(15,1);
    U = zeros(12,1);
    Y = zeros(3,1);
end

G(1) = 1 - r1/(r2*(r1 - 1) + 1);
G(2) = -(r1*(r2 - 1))/(r1*r2 - r2 + 1);
G(3) = (r1*r2)/(r1*r2 - r2 + 1);
G(4) = (r1*(r1 - 1))/(r2*(r1 - 1) + 1) - r1 + 1;
G(5) = -(r1*(r1 - 1)*(r2 - 1)*(r2 - 2))/(r1*r2 - r2 + 1);
G(6) = (r1*r2*(r1 - 1)*(r2 - 1))/(r1*r2 - r2 + 1);
G(7) = (r1^2*(r2 - 1)^2)/(r2*(r1 - 1) + 1);
G(8) = -(r1*r2*(r1 + 1)*(r2 - 1))/(r1*r2 - r2 + 1);
G(10) = (r1*r2^2)/(r1*r2 - r2 + 1);
G(10) = h1 - (h1*r1)/(r2*(r1 - 1) + 1);
G(11) = -(r1*(r2 - 1)*(h1 + h2 - h1*r2))/(r1*r2 - r2 + 1);
G(12) = (r1*r2*(h1 + h2 - h1*r2))/(r1*r2 - r2 + 1);
G(13) = ((h1 + h2*r1)*(r1 - 1)*(r2 - 1))/(r1*r2 - r2 + 1);
G(14) = -(r1*(h1 + h2*r1)*(r2 - 1))/(r1*r2 - r2 + 1);
G(15) = (h2*r1*r2)/(r1*r2 - r2 + 1);
U(1) = (r1*r2 - r2 + 1)^2;
U(2) = -(r1*r2 - r2 + 1)*(2*h1 - h1*r1 - 2*h1*r2 + 3*h2*r1 - h2*r1^2 + h2*r1^2*r2 + h1*r1*r2 - h2*r1*r2);
U(3) = r1*(r2 - 1)*(h1 - h2 + h2*r1)^2;
U(4) = (r1*r2 - r2 + 1)*(h2^2*r1 - h1^2*r2 + h1^2 + h1*h2*r1 - h1*h2*r1*r2);
U(5) = -r1*(r2 - 1)*(r1*r2 - r2 + 1)*(h1 - h2 + h2*r1);
U(6) = r1*(r2 - 1)*(h1 - h1*r2 + h2*r1)*(h1 - h2 + h2*r1);
U(7) = (r1*r2 - r2 + 1)^2;
U(8) = -(r1*r2 - r2 + 1)*(2*h1 - 2*h1*r2 + h2*r1 - h1*r1*r2^2 + h1*r1*r2 + h2*r1*r2);
U(9) = r1*(h2 - h1*r2)^2*(r2 - 1);
U(10) = (r1*r2 - r2 + 1)*(h2^2*r1 - h1^2*r2 + h1^2 + h1*h2*r1 - h1*h2*r1*r2);
U(11) = -r1*(h2 - h1*r2)*(r2 - 1)*(r1*r2 - r2 + 1);
U(12) = r1*(h2 - h1*r2)*(r2 - 1)*(h1 - h1*r2 + h2*r1);
Y(1) = G(1)*G(11)*G(15)-G(1)*G(12)*G(14)-G(2)*G(10)*G(15)+G(2)*G(12)*G(13)+G(3)*G(10)*G(14)-G(3)*G(11)*G(13);
Y(2) = G(3)*G(13) - G(1)*G(15);
Y(3) = G(10)*G(3) - G(12)*G(1);

end