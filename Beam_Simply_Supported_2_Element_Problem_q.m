%% Beam Element Stiffness Matrix

% Element lenght is Le, X Co-ordinate from xa to xb

clc;
clear;
syms x xa xb y1 y2 th1 th2 c1 c2 c3 c4 Le E I
xb = Le + xa;

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

%% Beam Problem - Using 2 Elements

% Simply Supported Beam under 'P' Point load at centre

ND = [1 2 3];               % Nodes
EL = [1 2];                 % Elements
DoF = [1 2 3 4 5 6];        % DoFs
% Moment of Inertia and Young's Modulus already defined as Symbols
RN = [1 5];                 % Restrained DoF
URN = [2 3 4 6];            % Unrestrained DoF
n_RN = length(RN);

syms R1 R2 R3 R4 R5 R6 u1 u2 u3 u4 u5 u6 P q L

F1 = [R1; 0; 0; 0; R5; 0];  % Add any Point Loads in this vector

n_DoF = length(DoF);        % Number of Nodes
n_EL = length(EL);          % Number of Elements
L = 2*Le;

% Element 1

DoF_1 = [1 2 3 4];
for i = 1:4
    for j = 1:4
        KeG1(DoF_1(1,i),DoF_1(1,j)) =  Ke(i,j);
    end
end
KeG1(n_DoF,n_DoF) = 0;      % Take out for last Element

F2_1(1,1) = int(w1(x)*q, x, xa, xb);
F2_1(1,1) = simplify(F2_1(1,1));
F2_1(2,1) = int(w2(x)*q, x, xa, xb);
F2_1(2,1) = simplify(F2_1(2,1));
F2_1(3,1) = int(w3(x)*q, x, xa, xb);
F2_1(3,1) = simplify(F2_1(3,1));
F2_1(4,1) = int(w4(x)*q, x, xa, xb);
F2_1(4,1) = simplify(F2_1(4,1));
F2_1(5,1) = 0;
F2_1(6,1) = 0;

% Element 2

DoF_2 = [3 4 5 6];
for i = 1:4
    for j = 1:4
        KeG2(DoF_2(1,i),DoF_2(1,j)) = Ke(i,j);
    end
end

F2_2(3,1) = int(w1(x)*q, x, xa, xb);
F2_2(3,1) = simplify(F2_2(3,1));
F2_2(4,1) = int(w2(x)*q, x, xa, xb);
F2_2(4,1) = simplify(F2_2(4,1));
F2_2(5,1) = int(w3(x)*q, x, xa, xb);
F2_2(5,1) = simplify(F2_2(5,1));
F2_2(6,1) = int(w4(x)*q, x, xa, xb);
F2_2(6,1) = simplify(F2_2(6,1));

% Global Matrices

KeG = KeG1 + KeG2;

F2 = F2_1 + F2_2;
F = F1 + F2;

% Solving

KeGTemp = KeG;
for i = 1:n_RN
    KeGTemp(RN(1,i),:) = 0;
    KeGTemp(:,RN(1,i)) = 0;
end

U = (pinv(KeGTemp2))*F;
U = simplify(U);                       % Displacement and Rotations

F_Final = KeG*U;

R = F_Final - F2;                      % Reactions - Not neat

for i = 1:n_URN
    R(URN(1,i)) = 0;                   % Reactions
end

U                                      % Output
R                                      % Output