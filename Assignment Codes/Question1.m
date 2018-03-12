% Name: 	Rohan Dharmadhikari
% Roll No.: 	154104002

%% 1D Element

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

%% Beam

syms x xa xb y1 y2 th1 th2 c1 c2 c3 c4 le E I;
xb = le + xa;

AB = [1 xa (xa^2) (xa^3); 0 1 (2*xa) (3*(xa^2)); 1 xb (xb^2) (xb^3); 0 1 (2*xb) (3*(xb^2))];
B = inv(AB);
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
Ke_Beam(1,1) = c1(1,2);
c2 = coeffs(k1(x), th1);
Ke_Beam(1,2) = c2(1,2);
c3 = coeffs(k1(x), y2);
Ke_Beam(1,3) = c3(1,2);
c4 = coeffs(k1(x), th2);
Ke_Beam(1,4) = c4(1,2);

w2(x) = ph2(x);
wd2(x) = diff(w2(x));
wdd2(x) = diff(wd2(x));

k2(x) = int((E*I)*ydd(x)*wdd2(x), x, xa, xb);
c1 = coeffs(k2(x), y1);
Ke_Beam(2,1) = c1(1,2);
c2 = coeffs(k2(x), th1);
Ke_Beam(2,2) = c2(1,2);
c3 = coeffs(k2(x), y2);
Ke_Beam(2,3) = c3(1,2);
c4 = coeffs(k2(x), th2);
Ke_Beam(2,4) = c4(1,2);

w3(x) = ph3(x);
wd3(x) = diff(w3(x));
wdd3(x) = diff(wd3(x));

k3(x) = int((E*I)*ydd(x)*wdd3(x), x, xa, xb);
c1 = coeffs(k3(x), y1);
Ke_Beam(3,1) = c1(1,2);
c2 = coeffs(k3(x), th1);
Ke_Beam(3,2) = c2(1,2);
c3 = coeffs(k3(x), y2);
Ke_Beam(3,3) = c3(1,2);
c4 = coeffs(k3(x), th2);
Ke_Beam(3,4) = c4(1,2);

w4(x) = ph4(x);
wd4(x) = diff(w4(x));
wdd4(x) = diff(wd4(x));

k4(x) = int((E*I)*ydd(x)*wdd4(x), x, xa, xb);
c1 = coeffs(k4(x), y1);
Ke_Beam(4,1) = c1(1,2);
c2 = coeffs(k4(x), th1);
Ke_Beam(4,2) = c2(1,2);
c3 = coeffs(k4(x), y2);
Ke_Beam(4,3) = c3(1,2);
c4 = coeffs(k4(x), th2);
Ke_Beam(4,4) = c4(1,2);

%% Solution

% Elements - AB, BC, CD
% Stiffness Matrix in local coordinates - [V1 Theta1 V2 Theta2 U1 U2]
Ke = Ke_Beam;
Ke(5,5) = Ke_1D(1,1);
Ke(5,6) = Ke_1D(1,2);
Ke(6,5) = Ke_1D(2,1);
Ke(6,6) = Ke_1D(2,2);

% Transformation Matrix

syms th

T = [cos(th) 0 0 0 -sin(th) 0;0 1 0 0 0 0; 0 0 cos(th) 0 0 -sin(th);
    0 0 0 1 0 0; sin(th) 0 0 0 cos(th) 0; 0 0 -sin(th) 0 0 cos(th)];

Ke_G = T.'*Ke*T;

% Element 1

Dof1 = [1 2 3 4 5 6];
E1 = 0.04385949*(10^11);    % pound-force/ft2
A1 = 15*(1/(12^2));         % ft2
I1 = 305*(1/(12^4));        % ft4
le1 = sqrt((10^2)+(20^2));        % ft
th1 = atan(2);              % Angle of Inclination

Ke1 = subs(Ke_G, [E, I, A, le, th], [E1, I1, A1, le1, th1]);
Ke1 = double(Ke1);

% Element 2 - 1

Dof2 = [3 4 7 8 6 9];
E1 = 0.04385949*(10^11);    % pound-force/ft2
A2 = 7.5*(1/(12^2));        % ft2
I2 = 125*(1/(12^4));        % ft4
le2 = 10;                   % ft
th2 = 0;                    % Angle of Inclination

Ke2 = subs(Ke_G, [E, I, A, le, th], [E1, I2, A2, le2, th2]);
Ke2 = double(Ke2);

% Element 2 - 1

Dof3 = [7 8 10 11 9 12];
E1 = 0.04385949*(10^11);    % pound-force/ft2
A3 = 7.5*(1/(12^2));        % ft2
I3 = 125*(1/(12^4));        % ft4
le3 = 10;                   % ft
th3 = 0;                    % Angle of Inclination

Ke3 = subs(Ke_G, [E, I, A, le, th], [E1, I3, A3, le3, th3]);
Ke3 = double(Ke3);

% Element 4

Dof4 = [10 11 13 14 12 15];
E1 = 0.04385949*(10^11);    % pound-force/ft2
A4 = 15*(1/(12^2));         % ft2
I4 = 305*(1/(12^4));        % ft4
le4 = sqrt((10^2)+(20^2));        % ft
th4 = (pi/2) - atan(2);              % Angle of Inclination

Ke4 = subs(Ke_G, [E, I, A, le, th], [E1, I4, A4, le4, th4]);
Ke4 = double(Ke4);

% Assembling Global Stiffness Matrix

KeG = zeros(15,15);

for i = 1:6
    for j = 1:6
        KeG(Dof1(1,i),Dof1(1,j)) = KeG(i,j) + Ke1(i,j);
    end
end

for i = 1:6
    for j = 1:6
        KeG(Dof2(1,i),Dof2(1,j)) =  KeG(i,j) + Ke2(i,j);
    end
end

for i = 1:6
    for j = 1:6
        KeG(Dof3(1,i),Dof3(1,j)) =  KeG(i,j) + Ke3(i,j);
    end
end

for i = 1:6
    for j = 1:6
        KeG(Dof4(1,i),Dof4(1,j)) =  KeG(i,j) + Ke4(i,j);
    end
end

% Loads And Restraints

RN = [1 2 5 13 14 15];                 % Restrained DoF
URN = [4 5 6 7 8 9];            % Unrestrained DoF
n_RN = length(RN);
n_URN = length(URN);

q1 = 1200;

F2_1(1,1) = int(w1(x)*q1, x, xa, xb);
F2_1(1,1) = simplify(F2_1(1,1));
F2_1(2,1) = int(w2(x)*q1, x, xa, xb);
F2_1(2,1) = simplify(F2_1(2,1));
F2_1(3,1) = int(w3(x)*q1, x, xa, xb);
F2_1(3,1) = simplify(F2_1(3,1));
F2_1(4,1) = int(w4(x)*q1, x, xa, xb);
F2_1(4,1) = simplify(F2_1(4,1));

F2 = subs(F2_1, le, 10);

% Load Vector
F = zeros(15,1);

F(3,1) = F(3,1) + F2(1,1);
F(4,1) = F(4,1) + F2(2,1);
F(7,1) = F(7,1) + F2(3,1) + F2(1,1);
F(8,1) = F(8,1) + F2(4,1) + F2(2,1);
F(10,1) = F(10,1) + F2(3,1);
F(11,1) = F(11,1) + F2(4,1);

KeG1 = KeG;

for i = 1:n_RN
    KeG1(RN(1,i),:) = 0;
    KeG1(:,RN(1,i)) = 0;
end

U = (pinv(KeG1))*F;

F_Final = KeG*U;

R = F_Final - F;

R1 = sqrt(R(1)^2 + R(2)^2);
R2 = sqrt(R(13)^2 + R(15)^2);

disp('Deflection at the centre of BC = ')
U(7)

disp('Reaction at A = ')
R1

disp('Reaction at D = ')
R2