% Name: 	Rohan Dharmadhikari
% Roll No.: 	154104002

clc;
clear;

%% Tria-3 CST Element

% Element Nodes - (x1,y1); (x2,y2); (x3,y3); (x4,y4)
% 1 Node has 3 DoF (x, y, z)
% Element DoF - [q1 q2 q3 q4 q5 q6 q7 q8 q9 q10 q11 q12]

syms x1 x2 x3 x4 y1 y2 y3 y4 z1 z2 z3 z4 q1 q2 q3 q4 q5 q6 q7 q8 q9 q10 q11 q12 x y z m n p E poi

% Points in new coordinate system

m1 = 1; m2 = 0; m3 = 0; m4 = 0;
n1 = 0; n2 = 1; n3 = 0; n4 = 0;
p1 = 0; p2 = 0; p3 = 1; p4 = 0;

% Shape Functions

N1 = m;
N2 = n;
N3 = p;
N4 = 1 - m - n - p;

x = N1.*x1 + N2.*x2 + N3.*x3 + N4.*x4;
y = N1.*y1 + N2.*y2 + N3.*y3 + N4.*y4;
z = N1.*z1 + N2.*z2 + N3.*z3 + N4.*z4;

xd_m = diff(x,m);
xd_n = diff(x,n);
xd_p = diff(x,p);
yd_m = diff(y,m);
yd_n = diff(y,n);
yd_p = diff(y,p);
zd_m = diff(z,m);
zd_n = diff(z,n);
zd_p = diff(z,p);

J = [xd_m yd_m zd_m; xd_n yd_n zd_n; xd_p yd_p zd_p];                             % Jacobian

A1 = inv(J)*[diff(N1,m); diff(N1,n); diff(N1,p)];
N1d_x = A1(1,1);
N1d_y = A1(2,1);
N1d_z = A1(3,1);

A2 = inv(J)*[diff(N2,m); diff(N2,n); diff(N2,p)];
N2d_x = A2(1,1);
N2d_y = A2(2,1);
N2d_z = A2(3,1);

A3 = inv(J)*[diff(N3,m); diff(N3,n); diff(N3,p)];
N3d_x = A3(1,1);
N3d_y = A3(2,1);
N3d_z = A3(3,1);

A4 = inv(J)*[diff(N4,m); diff(N4,n); diff(N4,p)];
N4d_x = A4(1,1);
N4d_y = A4(2,1);
N4d_z = A4(3,1);

% e = B*q           % e is strain

B(1,1) = N1d_x;
B(1,4) = N2d_x;
B(1,7) = N3d_x;
B(1,10) = N4d_x;
B(2,2) = N1d_y;
B(2,5) = N2d_y;
B(2,8) = N3d_y;
B(2,11) = N4d_y;
B(3,3) = N1d_z;
B(3,6) = N2d_z;
B(3,9) = N3d_z;
B(3,12) = N4d_z;
B(4,1) = N1d_y;
B(4,2) = N1d_x;
B(4,4) = N2d_y;
B(4,5) = N2d_x;
B(4,7) = N3d_y;
B(4,8) = N3d_x;
B(4,10) = N4d_y;
B(4,11) = N4d_x;
B(5,1) = N1d_z;
B(5,3) = N1d_x;
B(5,4) = N2d_z;
B(5,6) = N2d_x;
B(5,7) = N3d_z;
B(5,9) = N3d_x;
B(5,10) = N4d_z;
B(5,12) = N4d_x;
B(6,2) = N1d_z;
B(6,3) = N1d_y;
B(6,5) = N2d_z;
B(6,6) = N2d_y;
B(6,8) = N3d_z;
B(6,9) = N3d_y;
B(6,11) = N4d_z;
B(6,12) = N4d_y;

D = E.*(1/1+poi).*(1/1-(2.*poi)).*[1-poi poi poi 0 0 0; poi 1-poi poi 0 0 0; poi poi 1-poi 0 0 0;
                                   0 0 0 0.5.*poi 0 0; 0 0 0 0 0.5.*poi 0; 0 0 0 0 0 0.5.*poi];

D = subs(D, [E,poi], [(30.*(10.^6)),0.3]);

% SE = B.'*D*B*det(J);
% SE = simplify(SE);

%% Solution

% Splitting the cuboid into 6 tetrahedrons
% Cuboid - [1 2 3 4, 5 6 7 8] - points

ND = [1 2 3 4 5 6 7 8];                 % Nodes
EL = [1 2 3 4 5 6];                     % Elements

% Element 1

X1 = 0; X2 = 2.5; X3 = 0; X4 = 0;
Y1 = 0; Y2 = 0; Y3 = 0; Y4 = 20;
Z1 = 0; Z2 = 0; Z3 = 0.8; Z4 = 0;

B1 = subs(B, [x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4], [X1, X2, X3, X4, Y1, Y2, Y3, Y4, Z1, Z2, Z3, Z4]);
J1 = subs(J, [x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4], [X1, X2, X3, X4, Y1, Y2, Y3, Y4, Z1, Z2, Z3, Z4]);

SE1 = B1.'*D*B1*det(J1);

Ke1 = int(int(int(SE1,m,0,1-n-p),n,0,1-p),p,0,1);

Dof1 = [1 2 3 4 5 6 13 14 15 10 11 12];

% Element 2

X1 = 2.5; X2 = 0; X3 = 2.5; X4 = 0;
Y1 = 0; Y2 = 0; Y3 = 0; Y4 = 20;
Z1 = 0; Z2 = 0.8; Z3 = 0.8; Z4 = 0.8;

B1 = subs(B, [x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4], [X1, X2, X3, X4, Y1, Y2, Y3, Y4, Z1, Z2, Z3, Z4]);
J1 = subs(J, [x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4], [X1, X2, X3, X4, Y1, Y2, Y3, Y4, Z1, Z2, Z3, Z4]);

SE1 = B1.'*D*B1*det(J1);

Ke2 = int(int(int(SE1,m,0,1-n-p),n,0,1-p),p,0,1);

Dof2 = [4 5 6 13 14 15 16 17 18 22 23 24];

% Element 3

X1 = 2.5; X2 = 0; X3 = 0; X4 = 0;
Y1 = 0; Y2 = 20; Y3 = 0; Y4 = 20;
Z1 = 0; Z2 = 0; Z3 = 0.8; Z4 = 0.8;

B1 = subs(B, [x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4], [X1, X2, X3, X4, Y1, Y2, Y3, Y4, Z1, Z2, Z3, Z4]);
J1 = subs(J, [x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4], [X1, X2, X3, X4, Y1, Y2, Y3, Y4, Z1, Z2, Z3, Z4]);

SE1 = B1.'*D*B1*det(J1);

Ke3 = int(int(int(SE1,m,0,1-n-p),n,0,1-p),p,0,1);

Dof3 = [4 5 6 10 11 12 13 14 15 22 23 24];

% Element 4

X1 = 02.5; X2 = 2.5; X3 = 0; X4 = 2.5;
Y1 = 0; Y2 = 20; Y3 = 20; Y4 = 20;
Z1 = 0; Z2 = 0; Z3 = 0; Z4 = 0.8;

B1 = subs(B, [x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4], [X1, X2, X3, X4, Y1, Y2, Y3, Y4, Z1, Z2, Z3, Z4]);
J1 = subs(J, [x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4], [X1, X2, X3, X4, Y1, Y2, Y3, Y4, Z1, Z2, Z3, Z4]);

SE1 = B1.'*D*B1*det(J1);

Ke4 = int(int(int(SE1,m,0,1-n-p),n,0,1-p),p,0,1);

Dof4 = [4 5 6 7 8 9 10 11 12 19 20 21];

% Element 5

X1 = 0; X2 = 2.5; X3 = 2.5; X4 = 0;
Y1 = 20; Y2 = 0; Y3 = 20; Y4 = 20;
Z1 = 0; Z2 = 0.8; Z3 = 0.8; Z4 = 0.8;

B1 = subs(B, [x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4], [X1, X2, X3, X4, Y1, Y2, Y3, Y4, Z1, Z2, Z3, Z4]);
J1 = subs(J, [x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4], [X1, X2, X3, X4, Y1, Y2, Y3, Y4, Z1, Z2, Z3, Z4]);

SE1 = B1.'*D*B1*det(J1);

Ke5 = int(int(int(SE1,m,0,1-n-p),n,0,1-p),p,0,1);

Dof5 = [10 11 12 16 17 18 19 20 21 22 23 24];

% Element 6

X1 = 2.5; X2 = 0; X3 = 2.5; X4 = 2.5;
Y1 = 0; Y2 = 20; Y3 = 0; Y4 = 20;
Z1 = 0; Z2 = 0; Z3 = 0.8; Z4 = 0.8;

B1 = subs(B, [x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4], [X1, X2, X3, X4, Y1, Y2, Y3, Y4, Z1, Z2, Z3, Z4]);
J1 = subs(J, [x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4], [X1, X2, X3, X4, Y1, Y2, Y3, Y4, Z1, Z2, Z3, Z4]);

SE1 = B1.'*D*B1*det(J1);

Ke6 = int(int(int(SE1,m,0,1-n-p),n,0,1-p),p,0,1);

Dof6 = [4 5 6 10 11 12 16 17 18 19 20 21];

% Assembling Global Stiffness

KeG = zeros(24,24);

for i = 1:12
    for j = 1:12
        KeG(Dof1(1,i),Dof1(1,j)) = KeG(i,j) + Ke1(i,j);
    end
end

for i = 1:12
    for j = 1:12
        KeG(Dof2(1,i),Dof2(1,j)) = KeG(i,j) + Ke2(i,j);
    end
end

for i = 1:12
    for j = 1:12
        KeG(Dof3(1,i),Dof3(1,j)) = KeG(i,j) + Ke3(i,j);
    end
end

for i = 1:12
    for j = 1:12
        KeG(Dof4(1,i),Dof4(1,j)) = KeG(i,j) + Ke4(i,j);
    end
end

for i = 1:12
    for j = 1:12
        KeG(Dof1(1,i),Dof1(1,j)) = KeG(i,j) + Ke1(i,j);
    end
end

for i = 1:12
    for j = 1:12
        KeG(Dof5(1,i),Dof5(1,j)) = KeG(i,j) + Ke5(i,j);
    end
end

for i = 1:12
    for j = 1:12
        KeG(Dof6(1,i),Dof6(1,j)) = KeG(i,j) + Ke6(i,j);
    end
end

% Loads And Restraints

RN = [1 2 3 4 5 6 13 14 15 16 17 18];                 % Restrained DoF
n_RN = length(RN);

F = zeros(24,1);
F(24,1) = -600;

KeG1 = KeG;

for i = 1:n_RN
    KeG1(RN(1,i),:) = 0;
    KeG1(:,RN(1,i)) = 0;
end

U = (pinv(KeG1))*F;

F_Final = KeG*U;

R = F_Final - F;

disp('Deflection of Point 3 in z-direction is = ')
U(9)

disp('Deflection of Point 4 in z-direction is = ')
U(12)

disp('Deflection of Point 7 in z-direction is = ')
U(21)

disp('Deflection of Point 8 in z-direction is = ')
U(24)
