clear;
clc;

%% 1D Element Stiffness

syms x xa xb u1 u2 c1 c2 le E A;
xb = le + xa;

C = [1 xa; 1 xb];
B = inv(C);

ph1(x) = B(1,1) + B(2,1)*x;
ph2(x) = B(1,2) + B(2,2)*x;

ph1(x) = simplify(ph1(x));
ph2(x) = simplify(ph2(x));

y(x) = ph1(x)*u1 + ph2(x)*u2;
yd(x) = diff(y(x));

w1(x) = ph1(x);
wd1(x) = diff(w1(x));

k1(x) = int((E*A)*yd(x)*wd1(x), x, xa, xb);
c1 = coeffs(k1(x), u1);
Ke_Truss(1,1) = c1(1,2);
c2 = coeffs(k1(x), u2);
Ke_Truss(1,2) = c2(1,2);

w2(x) = ph2(x);
wd2(x) = diff(w2(x));

k2(x) = int((E*A)*yd(x)*wd2(x), x, xa, xb);
c1 = coeffs(k2(x), u1);
Ke_Truss(2,1) = c1(1,2);
c2 = coeffs(k2(x), u2);
Ke_Truss(2,2) = c2(1,2);

%% Truss Element Stiffness

% 1 Element has 2 nodes and 2 DoF per node (H & V)
% DoF = [u1 v1 u2 v2]

Ke_local(1,1) = Ke_Truss(1,1);
Ke_local(3,1) = Ke_Truss(2,1);
Ke_local(1,3) = Ke_Truss(1,2);
Ke_local(3,3) = Ke_Truss(2,2);
Ke_local(4,4) = 0;

% Truss element inclined at angle theta(th)

syms th

T = [cos(th) sin(th) 0 0; -sin(th) cos(th) 0 0;
     0 0 cos(th) sin(th); 0 0 -sin(th) cos(th)];
 
Ke_Global(th) = T.'*Ke_local*T;