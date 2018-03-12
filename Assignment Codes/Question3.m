% Name: 	Rohan Dharmadhikari
% Roll No.: 	154104002

%-----------------------------------------
%Q3
%-----------------------------------------

clear
clc
node= [1:1:22];%nodes
ne = [1:1:20];% Element
%X and Y coordinates
xo=10^-3*[50 34 50; 34 50 34; 50 34 50; 34 50 34; 50 34 50;34 50 34; 50 34 50;34 50 34;50 34 50;34 50 34; 50 34 50; 34 50 34; 50 34 50; 34 50 34; 50 34 50; 34 50 34; 50 34 50; 34 50 34; 50 34 50; 34 50 34];
yo=10^-3*[0 0 20; 0 20 20; 20 20 40;20 40 40;40 40 60;40 60 60;60 60 80; 60 80 80;80 80 100;80 100 100;100 100 120;100 120 120;120 120 140; 120 140 140;140 140 160;140 160 160;160 160 180; 160 180 180;180 180 200;180 200 200 ];

%Node connections
Nc = [ 1 2 3;2 3 4;3 4 5;4 5 6; 5 6 7; 6 7 8; 7 8 9; 8 9 10; 9 10 11; 10 11 12; 11 12 13; 12 13 14; 13 14 15; 14 15 16; 15 16 17; 16 17 18; 17 18 19;18 19 20; 19 20 21; 20 21 22];
ne= length(ne);%No of element
dof = length(node)*2;
dofn = [(2.*node-1)' (2.*node)'];
KG = zeros(dof,dof); 

syms EA le I;


for i=1:ne
    
n1=Nc(i,1);
n2=Nc(i,2);
n3=Nc(i,3);

  loc= [dofn(n1,:) dofn(n2,:) dofn(n3,:)];
  nj = length(loc);
  
  
%% Tria-3 CST Element

% Element Nodes - (x1,y1); (x2,y2); (x3,y3)
% 1 Node has 2 DoF (H & V)
% Element DoF - [q1 q2 q3 q4 q5 q6]

syms x1 x2 x3 y1 y2 y3 q1 q2 q3 q4 q5 q6 x y E poi m n

J = [(x1 - x3) (y1-y3); (x2 - x3) (y2 - y3)];       % Jacobian

ud = inv(J)*[(q1 - q5);(q3 - q5)];                  

ud_x = ud(1,1);                                     % P.Diff of u wrt x
ud_y = ud(2,1);                                     % P.Diff of u wrt y

vd = inv(J)*[(q2 - q6);(q4 - q6)];

vd_x = vd(1,1);                                     % P.Diff of v wrt x
vd_y = vd(2,1);                                     % P.Diff of v wrt y

x = m*x1+ n *x2+ (1-m-n)*x3; 
r(i) =(xo(i,1)+xo(i,2)+xo(i,3))/3;


% e = B*q           % e is strain

c1 = coeffs(ud_x, q1);
B(1,1) = c1(1,2);
% c2 = coeffs(ud_x, q2);
% B(1,2) = c2(1,2);
c3 = coeffs(ud_x, q3);
B(1,3) = c3(1,2);
% c4 = coeffs(ud_x, q4);
% B(1,4) = c4(1,2);
c5 = coeffs(ud_x, q5);
B(1,5) = c5(1,2);
% c6 = coeffs(ud_x, q6);
% B(1,6) = c6(1,2);

% c1 = coeffs(vd_x, q1);
% B(2,1) = c1(1,2);
c2 = coeffs(vd_x, q2);
B(2,2) = c2(1,2);
% c3 = coeffs(vd_x, q3);
% B(2,3) = c3(1,2);
c4 = coeffs(vd_x, q4);
B(2,4) = c4(1,2);
% c5 = coeffs(vd_x, q5);
% B(2,5) = c5(1,2);
c6 = coeffs(vd_x, q6);
B(2,6) = c6(1,2);

c1 = coeffs(ud_y, q1);
d1 = coeffs(vd_x, q1);
B(3,1) = c1(1,2);
c2 = coeffs(ud_y, q2);
d2 = coeffs(vd_x, q2);
B(3,2) = d2(1,2);
c3 = coeffs(ud_y, q3);
d3 = coeffs(vd_x, q3);
B(3,3) = c3(1,2);
c4 = coeffs(ud_y, q4);
d4 = coeffs(vd_x, q4);
B(3,4) = d4(1,2);
c5 = coeffs(ud_y, q5);
d5 = coeffs(vd_x, q5);
B(3,5) = c5(1,2);
c6 = coeffs(ud_y, q6);
d6 = coeffs(vd_x, q6);
B(3,6) = d6(1,2);
B(4,1) =r(i)/3;
B(4,2)=0;
B(4,3)=r(i)/3;
B(4,4)=0;
B(4,5)=r(i)/3;
D = E*(1/1+poi)*(1/1-(2*poi))*[1-poi poi 0 poi; poi 1-poi 0 poi; 0 0 0.5*(1-2*poi) 0; poi poi 0 1-poi];

%-----% For intergration over domain
 
SE = B.'*D*B;

% Problem with values
% Thickness is unity

E1 = 2.1*(10^11);
poi1 = 0.3;

SE = subs(SE, [E, poi, x1, x2, x3, y1, y2, y3], [E1, poi1, xo(i,1), xo(i,2), xo(i,3),yo(i,1), yo(i,2), yo(i,3) ]);
J = subs(J, [E, poi, x1, x2, x3, y1, y2, y3], [E1, poi1,  xo(i,1), xo(i,2), xo(i,3),yo(i,1), yo(i,2), yo(i,3)]);

x = subs (x, [x1,x2,x3],[xo(i,1), xo(i,2), xo(i,3)]);

Ke_2D_CST = 2*pi*abs(det(J))*r(i)*int(int(SE ,n, 0, 1-m), m, 0, 1);



for j =1:nj 
    L=loc(j);
     for k =1:nj
        L2= loc(k);
        KG(L,L2)=KG(L,L2)+ Ke_2D_CST(j,k);
        end
end
end

%Force Vector
%Force on Seach node = (2*pi*le*1N/mm2)/2
% Forces are on x direction on nodes 1,3,5,7,9,11,13,15,17,19,21 
f=2*0.034*10^6;
F=[f*pi 0 0 0 f*pi 0 0 0 f*pi 0 0 0 f*pi 0 0 0 f*pi 0 0 0 f*pi 0 0 0 f*pi 0 0 0 f*pi 0 0 0 f*pi 0 0 0 f*pi 0 0 0 f*pi 0 0 0 ]; 

%SOLVING
w=1;
RN = [1 2 3 4 41 42 43 44 ];                % Restrained DoF
URN = [3 4 5 6 7  8];                    % Unrestrained DoF
n_RN = length(RN);

KeGTemp1 = KG;
for i = 1:n_RN
    KeGTemp1(RN(1,i),:) = 0;
    KeGTemp1(:,RN(1,i)) = 0;
end

U= pinv(KeGTemp1)*F';

X_dash = [50+U(3) 50+U(7) 50+U(11) 50+U(15) 50+U(19) 50+U(23)  50+U(27)  50+U(31)  50+U(35)  50+U(39) 50+U(43)];
Y_dash = [ 0+U(4) 20+U(8) 40+U(12) 60+U(16) 80+U(20) 100+U(24) 120+U(28) 140+U(32) 160+U(36) 180+U(40) 200+U(44)];
u=[U(3) U(7) U(11) U(15) U(19) U(23)  U(27)  U(31)  U(35)  U(39) U(43)];
n=[0 20 40 60 80 100 120 140 160 180 200];
figure(1)

scatter(u,n)
figure(2)
for i=1:20
 
        x = xo(i,:);
        y= yo(i,:);
scatter(x,y)

hold on
end

