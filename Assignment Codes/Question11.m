% Name: 	Rohan Dharmadhikari
% Roll No.: 	154104002


clc;
clear;

%% Tria-3 CST Element

% Element Nodes - (x1,y1); (x2,y2); (x3,y3)
% 1 Node has 2 DoF (H & v)
% Element DoF - [q1 q2 q3 q4 q5 q6 q7 q8 q9]

syms x1 x2 x3 y1 y2 y3 q1 q2 q3 q4 q5 q6 q7 q8 q9 x y m n E poi

% Shape Functions

N1 = m;
N2 = n;
N3 = 1-m-n;

N = [N1 0 N2 0 N3 0; 0 N1 0 N2 0 N3];

% Element 1

X1 = 60; X2 = 60; X3 = 30;
Y1 = 0; Y2 = 20; Y3 = 0;

% Element 2

x1 = 60; x2 = 30; x3 = 30;
y1 = 20; y2 = 0; y3 = 20;

% Solving

% Over length 1-2

Ti1 = [-.5; 0; -0.5; 0; 0; 0];

N_I1 = N.'*N;
N_I1 = subs(N_I1, m+n, 1);

I1 = int(int(N_I1,n,0,1),m,0,1)*(sqrt((X1-X2)^2 + (Y2-Y1)^2));

T1 = I1*Ti1;

% Over length 1-3

Ti2 = [0; -0.5; 0; 0; 0; -0.5];
N = [N1 0 N2 0 N3 0; 0 N1 0 N2 0 N3];

N_I2 = N.'*N;
N_I2 = subs(N_I2, n, 0);

I2 = int(N_I2,m,0,1)*(sqrt((x1-x3)^2 + (y1-y3)^2));

T2 = I2*Ti2;

% Over length 2-3

N = [N1 0 N2 0 N3 0; 0 N1 0 N2 0 N3];

N_I3 = N.'*N;
N_I3 = subs(N_I3, m, 0);

I3 = int(N_I3,n,0,1);


F1 = T1(1)
F2 = T1(2)
F3 = T1(3) + T2(1)
F4 = T1(4) + T2(2)
F7 = T2(5)
F8 = T2(6)

