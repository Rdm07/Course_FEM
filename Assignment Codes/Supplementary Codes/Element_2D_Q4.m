clc;
clear;

%% Quad-4 Element

% Element Nodes - (x1,y1); (x2,y2); (x3,y3); (x4,y4)
% 1 Node has 2 DoF (H & V)
% Element DoF - [q1 q2 q3 q4 q5 q6 q7 q8]
% m,n represent zhi and neeta

syms x1 x2 x3 x4 y1 y2 y3 y4 q1 q2 q3 q4 q5 q6 q7 q8 x y E poi m n a b

% Points in new coordinate system

m1 = -1; m2 = 1; m3 = 1; m4 = -1;
n1 = -1; n2 = -1; n3 = 1; n4 = 1;

% Shape Functions

N1 = 0.25.*(1+(m.*m1)).*(1+(n.*n1));
N2 = 0.25.*(1+(m.*m2)).*(1+(n.*n2));
N3 = 0.25.*(1+(m.*m3)).*(1+(n.*n3));
N4 = 0.25.*(1+(m.*m4)).*(1+(n.*n4));

x = N1.*x1 + N2.*x2 + N3.*x3 + N4.*x4;
y = N1.*y1 + N2.*y2 + N3.*y3 + N4.*y4;

xd_m = diff(x,m);
xd_n = diff(x,n);
yd_m = diff(y,m);
yd_n = diff(y,n);

J = [xd_m yd_m; xd_n yd_n];                             % Jacobian

% N = [N1 0 N2 0 N3 0 N4 0; 0 N1 0 N2 0 N3 0 N4];       % N Matrix

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

% e = B*q           % e is strain

B(1,1) = N1d_x;
B(1,3) = N2d_x;
B(1,5) = N3d_x;
B(1,7) = N4d_x;
B(2,2) = N1d_y;
B(2,4) = N2d_y;
B(2,6) = N3d_y;
B(2,8) = N4d_y;
B(3,1) = N1d_y;
B(3,2) = N1d_x;
B(3,3) = N2d_y;
B(3,4) = N2d_x;
B(3,5) = N3d_y;
B(3,6) = N3d_x;
B(3,7) = N4d_y;
B(3,8) = N4d_x;      

% Problem with values
% Thickness is unity

E1 = 2.1.*(10.^11);
poi1 = 0.3;

B = subs(B, [E, poi, x1, x2, x3, x4, y1, y2, y3, y4], [E1, poi1, 0, 0, 10, 10, 0, 10, 0, 10]);      % B Matrix for specific element

D = E.*(1/1+poi).*(1/1-(2.*poi)).*[1-poi poi 0; poi 1-poi 0; 0 0 0.5.*poi];

D = subs(D, [E, poi, x1, x2, x3, x4, y1, y2, y3, y4], [E1, poi1, 0, 0, 10, 10, 0, 10, 0, 10]);
J = subs(J, [E, poi, x1, x2, x3, x4, y1, y2, y3, y4], [E1, poi1, 0, 0, 10, 10, 0, 10, 0, 10]);

SE = B.'*D*B*det(J);
SE = simplify(SE);

for i = 1:8
    for j = 1:8
        Ke_2D_Q4(i,j) = subs(SE(i,j), [m,n], [(-sqrt(1/3)), (-sqrt(1/3))]) + subs(SE(i,j), [m,n], [(sqrt(1/3)), (-sqrt(1/3))]) + subs(SE(i,j), [m,n], [(-sqrt(1/3)), (sqrt(1/3))]) + subs(SE(i,j), [m,n], [(sqrt(1/3)), (sqrt(1/3))]);
        Ke_2D_Q4(i,j) = simplify(Ke_2D_Q4(i,j));
    end
end
