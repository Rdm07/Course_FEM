% Name: 	Rohan Dharmadhikari
% Roll No.: 	154104002

%-------------------------
%Q5 
%--------------------------
clear;
clc;

node= [1:1:5];%Number of nodes
ne = [1:1:5];% No of elements
Le =[80.97 48.28 48.28 80.97 35]*10;
%Node connections
Nc = [ 1 2; 2 3; 4 3; 5 4; 2 4 ];
as=1;
ang = [351.119 338.75 21.25 8.88 pi/2];
ro = 10*(6 + 0.5*1);            % Outer Radius
ri = 10*(6 - 0.5*1);            % Inner Radius
A = pi*(ro^2-ri^2);         
I = 0.25*pi*(ro^4-ri^4);

nd = length(node);
ne= length(ne);
dof = nd*2;
dofn = [(2.*node-1)' (2.*node)'];
KG = zeros(dof,dof); 
syms E le I;
E =21*10^5;


for i=1:ne
    
n1=Nc(i,1);
n2=Nc(i,2);

  loc= [dofn(n1,:) dofn(n2,:)];
  nj = length(loc);
  c = cos(ang(i));
  s = sin(ang(i));
  T=[c s 0 0;-s c 0 0;0 0 c s; 0 0 -s c];
  ke= (E*A/Le(i))*[1 0 -1 0; 0 0 0 0;-1 0 1 0;0 0 0 0];
  Ke = T'*ke*T;
for j =1:nj 
    L=loc(j);
     for k =1:nj
        L2= loc(k);
        KG(L,L2)=KG(L,L2)+ Ke(j,k);
        end
end
end

% Force

F = [0 0 0 0 -10000 6000 0 0 0 0];

RN = [1 2 9 10];                % Restrained DoF
URN = [3 4 5 6 7 8];            % Unrestrained DoF
n_RN = length(RN);

KeGTemp1 = KG;
for i = 1:n_RN
    KeGTemp1(RN(1,i),:) = 0;
    KeGTemp1(:,RN(1,i)) = 0;
end

U = pinv(KeGTemp1)*F'% Deflection





