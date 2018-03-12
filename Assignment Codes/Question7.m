% Name: 	Rohan Dharmadhikari
% Roll No.: 	154104002


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

%% Solution

EL = [ 1 2 3;2 3 4;3 4 5;4 5 6; 5 6 7; 6 7 8; 7 8 9; 8 9 10; 9 10 11; 10 11 12];

for i = 1:10
    for k = 1:3
        Dof(i,2*k) = 2*(EL(i,k));
        Dof(i,2*k - 1) = 2*(EL(i,k)) - 1;
    end
end

x_c=[0 1 0.84 0.86 0.71 0.74 0.59 0.61 0.48 0.5 0.36 0.4];
y_c=[0 0.1 0.015 0.12 0.045 0.14 0.06 0.165 0.07 0.185  0.095 0.2];

Coor(:,1) = x_c.';
Coor(:,2) = y_c.';

% Properties

E1 = 70000;
poi1 = 0.3;

D = subs(D, [E, poi], [E1, poi1]);

% Solving

KeG = zeros(24,24);

for i = 1:10
    Bi = subs(B, [x1, x2, x3, y1, y2, y3], [Coor(EL(i,1),1), Coor(EL(i,2),1), Coor(EL(i,3),1), Coor(EL(i,1),2), Coor(EL(i,2),2), Coor(EL(i,3),2)]);
    Ji = subs(J, [x1, x2, x3, y1, y2, y3], [Coor(EL(i,1),1), Coor(EL(i,2),1), Coor(EL(i,3),1), Coor(EL(i,1),2), Coor(EL(i,2),2), Coor(EL(i,3),2)]);
    SE = Bi.'*D*Bi*det(Ji);
    SE = simplify(SE);
    Ke = int(int(SE,m,0,1-n),n,0,1);
    for j = 1:6
        for k =1:6
            KeG(Dof(i,j),Dof(i,k)) = KeG(Dof(i,j),Dof(i,k)) + Ke(j,k);
        end
    end
end

%SOLVING

RN = [2];
n_RN=length(RN);

KeGTemp1 = KeG;
for i = 1:n_RN
    KeGTemp1(RN(1,i),:) = 0;
    KeGTemp1(:,RN(1,i)) = 0;
end

P = 20:30;

for i = 1:11
    F(24,i) = P(i);

    U(:,i) = pinv(KeGTemp1)*F(:,i);
    Df(i) = U(24,i);   
end

for i = 1:11
    if Df(i) > 0.1
        Load = P(i);
        break
    end
end
