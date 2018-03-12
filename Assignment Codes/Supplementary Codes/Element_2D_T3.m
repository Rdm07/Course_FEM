clc;
clear;

%% Tria-3 CST Element

% Element Nodes - (x1,y1); (x2,y2); (x3,y3)
% 1 Node has 2 DoF (H & v)
% Element DoF - [q1 q2 q3 q4 q5 q6 q7 q8 q9]

syms x1 x2 x3 y1 y2 y3 q1 q2 q3 q4 q5 q6 q7 q8 q9 x y m n E poi

% Points in new coordinate system

m1 = 1; m2 = 0; m3 = 0; m4 = 0;
n1 = 0; n2 = 1; n3 = 0; n4 = 0;

% Shape Functions

N1 = m;
N2 = n;
N3 = 1-m-n;

x = N1.*x1 + N2.*x2 + N3.*x3;
y = N1.*y1 + N2.*y2 + N3.*y3;

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

% e = B*q           % e is strain

B(1,1) = N1d_x;
B(1,3) = N2d_x;
B(1,5) = N3d_x;
B(2,2) = N1d_y;
B(2,4) = N2d_y;
B(2,6) = N3d_y;
B(3,1) = N1d_y;
B(3,2) = N1d_x;
B(3,3) = N2d_y;
B(3,4) = N2d_x;
B(3,5) = N3d_y;
B(3,6) = N3d_x;

D = E*(1/1+poi)*(1/1-(2*poi))*[1-poi poi 0; poi 1-poi 0; 0 0 0.5*poi];

% Problem with values
% Thickness is unity

E1 = 2.1.*(10.^11);
poi1 = 0.3;

B = subs(B, [E, poi, x1, x2, x3, y1, y2, y3], [E1, poi1, 0, 0, 5, 0, 5, 4]);      % B Matrix for specific element

D = subs(D, [E, poi, x1, x2, x3, y1, y2, y3], [E1, poi1, 0, 0, 5, 0, 5, 4]);
J = subs(J, [E, poi, x1, x2, x3, y1, y2, y3], [E1, poi1, 0, 0, 5, 0, 5, 4]);

SE = B.'*D*B*det(J);
SE = simplify(SE);

Ke_2D_T3 = int(int(SE,m,0,1-n),n,0,1);