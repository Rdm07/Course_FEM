% Name: 	Rohan Dharmadhikari
% Roll No.: 	154104002

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           PROBLEM NO. 8
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PREPARED BY ANANTA SINHA (ROLL NO - 154104003)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   1. DATA PREPARATION::
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   A. USER INPUT (input values required only in this part)
%   a. Geometry of the Beam
nnd = 96;               % number of nodes
nel = 34;               % number of elements
nne = 8;                % number of nodes per elements
nodof = 3;              % number of degree of freedom per node
eldof = nne*nodof;      % number of degree of freedom per element
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Coordinates of each node
filename = 'geom_Q8.xlsx';
geom = xlsread(filename);        
%   Connectivity of each element        
filename = 'connec_Q8.xlsx';
connec = xlsread(filename);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dof = (1:nnd*nodof);
dof = reshape(dof,3,96)';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   b. Material Properties of Frame
E = 1.65e5;                                % Young's modulus(N/mm2)
nu = 0.25;                                 % Poisson's ratio
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   B. DATA MODIFICATION
%   Constitutive D Matrix 
D1 = 1-nu;
D2 = nu;
D3 = 0.5-nu;
D = [D1 D2 D2 0  0   0
     D2 D1 D2 0  0   0
     D2 D2 D1 0  0   0 
     0  0  0  D3 0   0
     0  0  0  0  D3  0
     0  0  0  0  0  D3];
D = E/((1+nu)*(1-2*nu))*D;             
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BC = [37:48 85:96];
BCN = [];                      % Boundary Condition
for i = 1:eldof
  H = [dof(BC(i),:)];
  BCN = [BCN;H];
end
BCN=reshape(BCN,1,72);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Py = 200;                      % Uniformly Distributed Line Load(N/mm)
F = zeros(nnd*nodof,1);
Y = [24 23 23 72 72 71];
Le = [25 25 200 200 25 25];
for i = 1:6
    y1 = (3*Y(i)-1);
    fe = -Py*Le(i)/2;
    F(y1) = F(y1)+fe;
end
 F(BCN,:)=[];
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %  C. ANALYSIS USING TETRAHEDRON ELEMENT
 for i = 1:2
    KG=zeros(nnd*nodof);
    if i == 1 
         connec = connec;
    else
         connec = [connec; 29 30 9 10 33 34 13 14;77 78 57 58 81 82 61 62];
    end
    elem = [];         % Node Connectivity of All Elements
    for j = 1:length(connec)
        j = connec(j,:); 
        elem1 = j;
        elem1 =    [j(1,1) j(1,2) j(1,4) j(1,8);
                    j(1,1) j(1,2) j(1,8) j(1,5);
                    j(1,2) j(1,8) j(1,5) j(1,6);
                    j(1,1) j(1,3) j(1,4) j(1,7);
                    j(1,1) j(1,7) j(1,8) j(1,5);
                    j(1,1) j(1,8) j(1,4) j(1,7)];
        elem = [elem; elem1];
    end             
    for j = 1:6*length(connec)
        n1 = elem(j,1);
        n2 = elem(j,2);
        n3 = elem(j,3);
        n4 = elem(j,4);
        
        geoco = [geom(n1,:); geom(n2,:); geom(n3,:); geom(n4,:)];
        
        x14 = geoco(1,1)-geoco(4,1);
        y14 = geoco(1,2)-geoco(4,2);
        z14 = geoco(1,3)-geoco(4,3);
        
        x24 = geoco(2,1)-geoco(4,1);
        y24 = geoco(2,2)-geoco(4,2);
        z24 = geoco(2,3)-geoco(4,3);
        
        x34 = geoco(3,1)-geoco(4,1);
        y34 = geoco(3,2)-geoco(4,2);
        z34 = geoco(3,3)-geoco(4,3);
        
        J = [x14  y14  z14;             % Jacobian Matrix
             x24  y24  z24;
             x34  y34  z34];    
        dj=det(J);
        Ve=abs(dj)/6;
        A=inv(J);
        A1 = sum(A(1,:));
        A2 = sum(A(2,:));
        A3 = sum(A(3,:));
B =[A(1,1)   0      0    A(1,2)   0      0    A(1,3)   0      0    -A1  0   0 
      0    A(2,1)   0      0    A(2,2)   0      0    A(2,3)   0     0  -A2  0
      0      0    A(3,1)   0      0    A(3,2)   0      0    A(3,3)  0   0  -A3
      0    A(3,1) A(2,1)   0    A(3,2) A(2,2)   0    A(3,3) A(2,3)  0  -A3 -A2
    A(3,1)   0    A(1,1) A(3,2)   0    A(1,2) A(3,3)   0    A(1,3) -A3  0  -A1
    A(2,1) A(1,1)   0    A(2,2) A(1,2)   0    A(2,3) A(1,3)   0    -A2 -A1  0];

    Ke=Ve*B'*D*B;

    loc = [dof(n1,:) dof(n2,:) dof(n3,:) dof(n4,:)];
        for m = 1:length(loc)
            m1 = loc(m);
            for n = 1:length(loc)
                n1 = loc(n);
                KG(m1,n1) = KG(m1,n1) + Ke(m,n);
            end
        end
    end
KG(BCN,:) = [];
KG(:,BCN) = [];
disp = zeros(nnd*nodof,1);
disp([1:108 145:252],1) = (KG)\F;

    if i == 1
        display('SECTION WITH OPENINGS')
    else
        display('SECTION WITHOUT OPENINGS')
    end
    
MAX_DEFLECTION_mm = max(abs(disp))
DISPLACEMENT_UNDER_LOAD_mm = disp([71 68 215 212])

% Determination of axial strains and stresses in global coordinates  
q = disp([loc]);
sigma = D*B*q;
% Calculation of Principal Stress
I1 = sigma(1) + sigma(2) + sigma(3);
I2 = sigma(1)*sigma(2)+sigma(2)*sigma(3)+sigma(3)*sigma(1)-sigma(4)^2-sigma(5)^2-sigma(6)^2;
I3 = sigma(1)*sigma(2)*sigma(3)+2*sigma(4)*sigma(5)*sigma(6)-sigma(1)*sigma(5)^2-sigma(2)*sigma(6)^2-sigma(3)*sigma(4)^2;

a = (I1^2)/3 - I2;
b = -2*(I1/3)^3 + I1*I2/3 - I3;
c = 2*sqrt(a/3);
theta = (1/3)*acos(-3*b/(a*c));

ps1(i) = I1/3 + c*cosd(theta);
ps2(i) = I1/3 + c*cosd(theta + 2*pi/3);
ps3(i) = I1/3 + c*cosd(theta + 4*pi/3);

MAXIMUM_PRINCIPLE_STRESS_MPa_X_DIRECTION = max(ps1)
MAXIMUM_PRINCIPLE_STRESS_MPa_Y_DIRECTION = max(ps2)
MAXIMUM_PRINCIPLE_STRESS_MPa_Z_DIRECTION = max(ps3)
 end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%