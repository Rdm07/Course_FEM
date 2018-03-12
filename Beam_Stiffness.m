clc;
clear;
syms x xa xb y1 y2 th1 th2 c1 c2 c3 c4 le E I;
xb = le + xa;

A = [1 xa (xa^2) (xa^3); 0 1 (2*xa) (3*(xa^2)); 1 xb (xb^2) (xb^3); 0 1 (2*xb) (3*(xb^2))];
B = inv(A);
ph1(x) = B(1,1) + B(2,1)*x + B(3,1)*(x^2) + B(4,1)*(x^3);
ph2(x) = B(1,2) + B(2,2)*x + B(3,2)*(x^2) + B(4,2)*(x^3);
ph3(x) = B(1,3) + B(2,3)*x + B(3,3)*(x^2) + B(4,3)*(x^3);
ph4(x) = B(1,4) + B(2,4)*x + B(3,4)*(x^2) + B(4,4)*(x^3);

ph1(x) = simplify(ph1(x));
ph2(x) = simplify(ph2(x));
ph3(x) = simplify(ph3(x));
ph4(x) = simplify(ph4(x));

y(x) = ph1(x)*y1 + ph2(x)*th1 + ph3(x)*y2 + ph4(x)*th2;
yd(x) = diff(y(x));
ydd(x) = diff(yd(x));

w1(x) = ph1(x);
wd1(x) = diff(w1(x));
wdd1(x) = diff(wd1(x));

k1(x) = int((E*I)*ydd(x)*wdd1(x), x, xa, xb);
c1 = coeffs(k1(x), y1);
Ke(1,1) = c1(1,2);
c2 = coeffs(k1(x), th1);
Ke(1,2) = c2(1,2);
c3 = coeffs(k1(x), y2);
Ke(1,3) = c3(1,2);
c4 = coeffs(k1(x), th2);
Ke(1,4) = c4(1,2);

w2(x) = ph2(x);
wd2(x) = diff(w2(x));
wdd2(x) = diff(wd2(x));

k2(x) = int((E*I)*ydd(x)*wdd2(x), x, xa, xb);
c1 = coeffs(k2(x), y1);
Ke(2,1) = c1(1,2);
c2 = coeffs(k2(x), th1);
Ke(2,2) = c2(1,2);
c3 = coeffs(k2(x), y2);
Ke(2,3) = c3(1,2);
c4 = coeffs(k2(x), th2);
Ke(2,4) = c4(1,2);

w3(x) = ph3(x);
wd3(x) = diff(w3(x));
wdd3(x) = diff(wd3(x));

k3(x) = int((E*I)*ydd(x)*wdd3(x), x, xa, xb);
c1 = coeffs(k3(x), y1);
Ke(3,1) = c1(1,2);
c2 = coeffs(k3(x), th1);
Ke(3,2) = c2(1,2);
c3 = coeffs(k3(x), y2);
Ke(3,3) = c3(1,2);
c4 = coeffs(k3(x), th2);
Ke(3,4) = c4(1,2);

w4(x) = ph4(x);
wd4(x) = diff(w4(x));
wdd4(x) = diff(wd4(x));

k4(x) = int((E*I)*ydd(x)*wdd4(x), x, xa, xb);
c1 = coeffs(k4(x), y1);
Ke(4,1) = c1(1,2);
c2 = coeffs(k4(x), th1);
Ke(4,2) = c2(1,2);
c3 = coeffs(k4(x), y2);
Ke(4,3) = c3(1,2);
c4 = coeffs(k4(x), th2);
Ke(4,4) = c4(1,2);

Ke
