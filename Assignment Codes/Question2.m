% Name: 	Rohan Dharmadhikari
% Roll No.: 	154104002

%%Q2: Quad 4 meshing
%6 elements have been used


clear;
clc;
node= [1:1:12];%Number of nodes
ne = [1:1:6];% No of elements
xo=[3 4.5 2.75 3.75;4.5 3.75 4.7 6;2.75 3.75 1.65 2.4;3.75 2.4 3.2 4.7;1.65 2.4 0 0;2.4 0 0 3.2 ];
yo=[0 0 1.35 1.9;0 1.9 2.5 0;1.35 1.9 2.6 3.75;1.9 3.75 5.196 2.5;2.6 3.75 3 4.2;3.75 4.2 5.196 5.196];
for i=1:6
 
        x = xo(i,:);
        y= yo(i,:);
scatter(x , y);

hold on
    
end


%Node connections
Nc = [ 1 2 3 4; 2 4 11 12; 3 4 5 6; 4 6 10 11; 5 6 7 8; 6 8 9 10]; 
nd = length(node);
ne= length(ne);
dof = nd*2;
dofn = [(2.*node-1)' (2.*node)'];
KG = zeros(dof,dof); 

syms EA le I;

for i=1:ne
    
n1=Nc(i,1);
n2=Nc(i,2);
n3=Nc(i,3);
n4=Nc(i,4);

  loc= [dofn(n1,:) dofn(n2,:) dofn(n3,:) dofn(n4,:)];
  nj = length(loc);
%% Quad-4 Element

syms x1 x2 x3 x4 y1 y2 y3 y4 q1 q2 q3 q4 q5 q6 q7 q8 x y E poi m n ab

m1 = -1; m2 = 1; m3 = 1; m4 = -1;
n1 = -1; n2 = -1; n3 = 1; n4 = 1;

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

E1 = 2.*(10.^11);
poi1 = 0.3;

B = subs(B, [E, poi, x1, x2, x3, x4, y1, y2, y3, y4], [E1, poi1,  xo(i,1), xo(i,2), xo(i,3), xo(i,4),yo(i,1), yo(i,2), yo(i,3), yo(i,4)]);      % B Matrix for specific element

D = E.*(1/1+poi).*(1/1-(2.*poi)).*[1-poi poi 0; poi 1-poi 0; 0 0 0.5.*poi];

D = subs(D, [E, poi, x1, x2, x3, x4, y1, y2, y3, y4], [E1, poi1, xo(i,1), xo(i,2), xo(i,3), xo(i,4),yo(i,1), yo(i,2), yo(i,3), yo(i,4)]);
J = subs(J, [E, poi, x1, x2, x3, x4, y1, y2, y3, y4], [E1, poi1,  xo(i,1), xo(i,2), xo(i,3), xo(i,4),yo(i,1), yo(i,2), yo(i,3), yo(i,4)]);

SE = B.'*D*B*det(J);
SE = simplify(SE);

for i = 1:8
    for j = 1:8
        Ke_2D_Q4(i,j) = subs(SE(i,j), [m,n], [(-sqrt(1/3)), (-sqrt(1/3))]) + subs(SE(i,j), [m,n], [(sqrt(1/3)), (-sqrt(1/3))]) + subs(SE(i,j), [m,n], [(-sqrt(1/3)), (sqrt(1/3))]) + subs(SE(i,j), [m,n], [(sqrt(1/3)), (sqrt(1/3))]);
        Ke_2D_Q4(i,j) = simplify(Ke_2D_Q4(i,j));
    end
end


for j =1:nj 
    L=loc(j);
     for k =1:nj
        L2= loc(k);
        KG(L,L2)=KG(L,L2)+ Ke_2D_Q4(j,k);
        end
end
end


T(24) = 0;
T(18)= -8000;
T(20) = -8000;

RN = [1 2 3 4 13 15 17 23 24];                 % Restrained DoF
URN = [2 3 4 6];
n_RN=length(RN);

KeGTemp = KG;

for i = 1:n_RN
    KeGTemp(RN(1,i),:) = 0;
    KeGTemp(:,RN(1,i)) = 0;
end

U = (pinv(KeGTemp))*T';
f= KG*U
% Displacement and Rotations
k=1;
for i = 1:6
n1=Nc(i,1);
n2=Nc(i,2);
n3=Nc(i,3);
n4=Nc(i,4);
loc= [dofn(n1,:) dofn(n2,:) dofn(n3,:) dofn(n4,:)];
    B = subs(B, [E, poi, x1, x2, x3, x4, y1, y2, y3, y4, m ,n], [E1, poi1,  xo(i,1), xo(i,2), xo(i,3), xo(i,4),yo(i,1), yo(i,2), yo(i,3), yo(i,4),0 , 0]);
    u = U(loc);
        S = todecimal(D*B*u);
        Stress(:,i)=simplify(S);
end
   
Stress_In_first_Element=Stress(:,1)
Stress_In_second_Element=Stress(:,2)
Stress_In_third_Element = Stress(:,3)
Stress_In_fourth_Element = Stress(:,4)
Stress_In_fifth_Element =Stress(:,5)
Stress_In_sixth_Element = Stress(:,6)

