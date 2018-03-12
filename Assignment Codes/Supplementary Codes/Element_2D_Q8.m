clc;
clear;

%% Quad-8 Element

% Element Nodes - (x1,y1); (x2,y2); (x3,y3); (x4,y4); (x5,y5); (x6,y6); (x7,y7); (x8,y8)
% 1 Node has 2 DoF (H & V)
% Element DoF - [q1 q2 q3 q4 q5 q6 q7 q8 q9 q10 q11 q12 q13 q14 q15 q16]
% m,n represent zhi and neeta

syms x1 x2 x3 x4 x5 x6 x7 x8 y1 y2 y3 y4 y5 y6 y7 y8 q1 q2 q3 q4 q5 q6 q7 q8 q9 q10 q11 q12 q13 q14 q15 q16 x y E poi m n a b

% Shape Functions

N1 = -0.25.*(1-m)*(1-n)*(1+m+n);
N2 = 0.5*(1+n)*(1-n)*(1-m);
N3 = -0.25.*(1-m)*(1+n)*(1+m-n);
N4 = 0.5*(1+m)*(1-m)*(1+n);
N5 = -0.25.*(1+m)*(1+n)*(1-m-n);
N6 = 0.5*(1+n)*(1-n)*(1+m);
N7 = -0.25.*(1+m)*(1-n)*(1-m+n);
N8 = 0.5*(1+m)*(1-m)*(1-n);

x = N1*x1 + N2*x2 + N3*x3 + N4*x4 + N5*x5 + N6*x6 + N7*x7 + N8*x8;
y = N1*y1 + N2*y2 + N3*y3 + N4*y4 + N5*y5 + N6*y6 + N7*y7 + N8*y8;

xd_m = diff(x,m);
xd_n = diff(x,n);
yd_m = diff(y,m);
yd_n = diff(y,n);

J = [xd_m yd_m; xd_n yd_n];                             % Jacobian

A1 = inv(J)*[diff(N1,m); diff(N1,n)];
N1d_x = A1(1,1);
N1d_y = A1(2,1);

A2 = inv(J)*[diff(N2,m); diff(N2,n)];
N2d_x = A2(1,1);
N2d_y = A2(2,1);

A3 = inv(J)*[diff(N3,m); diff(N3,n)];
N3d_x = A3(1,1);
N3d_y = A3(2,1);

A4 = inv(J)*[diff(N4,m); diff(N4,n)];
N4d_x = A4(1,1);
N4d_y = A4(2,1);

A5 = inv(J)*[diff(N5,m); diff(N5,n)];
N5d_x = A5(1,1);
N5d_y = A5(2,1);

A6 = inv(J)*[diff(N6,m); diff(N6,n)];
N6d_x = A6(1,1);
N6d_y = A6(2,1);

A7 = inv(J)*[diff(N7,m); diff(N7,n)];
N7d_x = A7(1,1);
N7d_y = A7(2,1);

A8 = inv(J)*[diff(N8,m); diff(N8,n)];
N8d_x = A8(1,1);
N8d_y = A8(2,1);


% e = B*q           % e is strain

B(1,1) = N1d_x;
B(1,3) = N2d_x;
B(1,5) = N3d_x;
B(1,7) = N4d_x;
B(1,9) = N5d_x;
B(1,11) = N6d_x;
B(1,13) = N7d_x;
B(1,15) = N8d_x;
B(2,2) = N1d_y;
B(2,4) = N2d_y;
B(2,6) = N3d_y;
B(2,8) = N4d_y;
B(2,10) = N5d_y;
B(2,12) = N6d_y;
B(2,14) = N7d_y;
B(2,16) = N8d_y;
B(3,1) = N1d_y;
B(3,2) = N1d_x;
B(3,3) = N2d_y;
B(3,4) = N2d_x;
B(3,5) = N3d_y;
B(3,6) = N3d_x;
B(3,7) = N4d_y;
B(3,8) = N4d_x;
B(3,9) = N5d_y;
B(3,10) = N5d_x;
B(3,11) = N6d_y;
B(3,12) = N6d_x;
B(3,13) = N7d_y;
B(3,14) = N7d_x;
B(3,15) = N8d_y;
B(3,16) = N8d_x;

D = E.*(1/1+poi).*(1/1-(2.*poi)).*[1-poi poi 0; poi 1-poi 0; 0 0 0.5.*poi];

% Problem with values
% Thickness is unity

E1 = 2.1.*(10.^11);
poi1 = 0.3;

B = subs(B, [E, poi, x1, x2, x3, x4, x5, x6, x7, x8, y1, y2, y3, y4, y5, y6, y7, y8], [E1, poi1, 0, 10, 10, 0, 5, 10, 5, 0, 0, 0, 10, 10, 0, 5, 10, 5]);      % B Matrix for specific element
B = simplify(B);
D = subs(D, [E, poi, x1, x2, x3, x4, x5, x6, x7, x8, y1, y2, y3, y4, y5, y6, y7, y8], [E1, poi1, 0, 10, 10, 0, 5, 10, 5, 0, 0, 0, 10, 10, 0, 5, 10, 5]);
J = subs(J, [E, poi, x1, x2, x3, x4, x5, x6, x7, x8, y1, y2, y3, y4, y5, y6, y7, y8], [E1, poi1, 0, 10, 10, 0, 5, 10, 5, 0, 0, 0, 10, 10, 0, 5, 10, 5]);
J = simplify(J);

SE = B.'*D*B*det(J);
SE = simplify(SE);

for i = 1:16
    for j = 1:16
        Ke_2D_Q8(i,j) = subs(SE(i,j), [m,n], [(-sqrt(1/3)), (-sqrt(1/3))]) + subs(SE(i,j), [m,n], [(sqrt(1/3)), (-sqrt(1/3))]) + subs(SE(i,j), [m,n], [(-sqrt(1/3)), (sqrt(1/3))]) + subs(SE(i,j), [m,n], [(sqrt(1/3)), (sqrt(1/3))]);
        Ke_2D_Q8(i,j) = simplify(Ke_2D_Q8(i,j));
    end
end

