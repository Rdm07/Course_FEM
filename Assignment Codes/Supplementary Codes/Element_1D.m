clc;
clear;
syms x x1 x2 u1 u2 c1 c2 le E A;
x2 = le + x1;

C = [1 x1; 1 x2];
B = inv(C);

ph1(x) = B(1,1) + B(2,1)*x;
ph2(x) = B(1,2) + B(2,2)*x;

ph1(x) = simplify(ph1(x));
ph2(x) = simplify(ph2(x));

y(x) = ph1(x)*u1 + ph2(x)*u2;
yd(x) = diff(y(x));

v1(x) = ph1(x);
vd1(x) = diff(v1(x));

k1(x) = int((E*A)*yd(x)*vd1(x), x, x1, x2);
c1 = coeffs(k1(x), u1);
Ke_1D(1,1) = c1(1,2);
c2 = coeffs(k1(x), u2);
Ke_1D(1,2) = c2(1,2);

v2(x) = ph2(x);
vd2(x) = diff(v2(x));

k2(x) = int((E*A)*yd(x)*vd2(x), x, x1, x2);
c1 = coeffs(k2(x), u1);
Ke_1D(2,1) = c1(1,2);
c2 = coeffs(k2(x), u2);
Ke_1D(2,2) = c2(1,2);