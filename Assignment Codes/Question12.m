% Name: 	Rohan Dharmadhikari
% Roll No.: 	154104002


clc;
clear;

tic;

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

N = [N1 0 N2 0 N3 0 N4 0; 0 N1 0 N2 0 N3 0 N4];

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

D = E.*(1/1+poi).*(1/1-(2.*poi)).*[1-poi poi 0; poi 1-poi 0; 0 0 0.5.*poi];

% SE = B.'*D*B*det(J);
% SE = simplify(SE);

%% Solution

% We can consider the given block as 2D block with constant thickness 10mm

t = 10;

% Elements

EL = [1 2 3 4; 3 4 6 5; 4 6 8 7; 7 8 9 10; 10 7 11 12; 11 12 14 13; 12 14 15 16; 15 16 18 17; 16 18 20 19];

for i = 1:9
    for k = 1:4
        Dof(i,2*k) = 2*(EL(i,k));
        Dof(i,2*k - 1) = 2*(EL(i,k)) - 1;
    end
end

Coor = [0,3; 0,0; 7,0; 7,3; 10,0; 10,3; 7,10.5; 10,10.5; 10,13.5; 7,13.5;
        3,10.5; 3,13.5; 0,10.5; 0,13.5; 0,21; 3,21; 0,24; 3,24; 10,21; 10,24];

% Properties

E1 = 70000;
poi1 = 0.3;

D = subs(D, [E, poi], [E1, poi1]);

% Solving

KeG = zeros(40,40);

for i = 1:9
    Bi = subs(B, [x1, x2, x3, x4, y1, y2, y3, y4], [Coor(EL(i,1),1), Coor(EL(i,2),1), Coor(EL(i,3),1), Coor(EL(i,4),1), Coor(EL(i,1),2), Coor(EL(i,2),2), Coor(EL(i,3),2), Coor(EL(i,4),2)]);
    Bi = simplify(Bi);
    Ji = subs(J, [x1, x2, x3, x4, y1, y2, y3, y4], [Coor(EL(i,1),1), Coor(EL(i,2),1), Coor(EL(i,3),1), Coor(EL(i,4),1), Coor(EL(i,1),2), Coor(EL(i,2),2), Coor(EL(i,3),2), Coor(EL(i,4),2)]);
    Ji = simplify(Ji);
    SE = Bi.'*D*Bi*t*det(Ji);
    for j = 1:8
        for k = 1:8
            Ke(j,k) = subs(SE(j,k), [m,n], [(-sqrt(1/3)), (-sqrt(1/3))]) + subs(SE(j,k), [m,n], [(sqrt(1/3)), (-sqrt(1/3))]) + subs(SE(j,k), [m,n], [(-sqrt(1/3)), (sqrt(1/3))]) + subs(SE(j,k), [m,n], [(sqrt(1/3)), (sqrt(1/3))]);
            KeG(Dof(i,j),Dof(i,k)) = KeG(Dof(i,j),Dof(i,k)) + Ke(j,k);
        end
    end
end

% Load

N_I = N.'*N;

% Element 8

Ti1 = [0; 0; 0; 0; 0; -200; 0; -200];
NI1 = subs(N_I, n, 1);
I1 = int(NI1, m, -1, 1);
T1 = (I1*Ti1)*3;

% Element 9

Ti2 = [0; 0; 0; -200; 0; -200; 0; 0];
NI2 = subs(N_I, m, 1);
I2 = int(NI2, n, -1, 1);
T2 = (I2*Ti2)*7;

% Assembling Force Vector

F(34,1) = T1(8);
F(36,1) = T1(6) + T2(4);
F(40,1) = T2(6);

% Restraints

RN = [3 4 5 6 9 10];                 % Restrained DoF
n_RN = length(RN);

KeG1 = KeG;

for i = 1:n_RN
    KeG1(RN(1,i),:) = 0;
    KeG1(:,RN(1,i)) = 0;
end


U = (pinv(KeG1))*F;

F_Final = KeG*U;

R = F_Final - F;

R = double(R);                  % Reactions
U = double(U);                  % Deflections

Xi = Coor(:,1);
Yi = Coor(:,2);

for i = 1:20
    CoorN(i,1) = Coor(i,1) + U(2*1 - 1);
    CoorN(i,2) = Coor(i,2) + U(2*1);
end

XiN = CoorN(:,1);
YiN = CoorN(:,2);

scatter(Xi,Yi)
hold on
scatter(XiN,YiN)

TimeSpent = toc;