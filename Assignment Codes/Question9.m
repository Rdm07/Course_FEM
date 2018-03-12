% Name: 	Rohan Dharmadhikari
% Roll No.: 	154104002

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           PROBLEM NO. 9
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           PREPARED BY ANANTA SINHA (ROLL NO - 154104003)
%COMMON STAPLER PROBLEM(CASE a b c d)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   1. DATA PREPARATION::
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   A. USER INPUT (input values required only in this part)
%   a. Geometry of Frame
 ND=[1 2 3 4 ];                     %Node numbers
 NC=[1 2;2 3;3 4];                  %Node connectivity
 L1=6.31; L2=12.63; L3=6.31;        %Length of the element
 alf1=pi/2; alf2=0; alf3=(3*pi/2);  %Angle
 EI=623.469;EA=39.585e3;
 dof=3*length(ND);                  %degree of freedom at each node
 KG=zeros(dof,dof);
 f11=[0 -60 -126.28 0 -60 126.28];   % distributed load
 f12=[0 0 -60 -126.28 0 -60 126.28 0];
 
 f21=[0 -60 0 0 -60 0];              %point load
 f22=[0 0 -60 0 0 -60 0 0];
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   2. ANALYSIS OF STRUCTURE

 c1= cos(alf1); c2= cos(alf2); c3= cos(alf3);
 s1= sin(alf1); s2= sin(alf2); s3= sin(alf3);
 
 T1=[c1 s1 0 0 0 0;-s1 c1 0 0 0 0;0 0 1 0 0 0;0 0 0 c1 s1 0;0 0 0 -s1 c1 0;0 0 0 0 0 1];
 ke1=[(EA/L1) 0 0 (-EA/L1) 0 0;0 ((12*EI)/L1.^3) ((6*EI)/L1.^2) 0 (-(12*EI)/L1.^3) ((6*EI)/L1.^2);0 ((6*EI)/L1.^2) ((4*EI)/L1) 0 (-(6*EI)/L1.^2) ((2*EI)/L1);(-EA/L1) 0 0 (EA/L1) 0 0;0 (-(12*EI)/L1.^3) (-(6*EI)/L1.^2) 0 ((12*EI)/L1.^3) (-(6*EI)/L1.^2);0 ((6*EI)/L1.^2) ((2*EI)/L1) 0 (-(6*EI)/L1.^2) ((4*EI)/L1)]; 
 Ke1=T1'*ke1*T1;
 
 T2=[c2 s2 0 0 0 0;-s2 c2 0 0 0 0;0 0 1 0 0 0;0 0 0 c2 s2 0;0 0 0 -s2 c2 0;0 0 0 0 0 1];
 ke2=[(EA/L2) 0 0 (-EA/L2) 0 0;0 ((12*EI)/L2.^3) ((6*EI)/L2.^2) 0 (-(12*EI)/L2.^3) ((6*EI)/L2.^2);0 ((6*EI)/L2.^2) ((4*EI)/L2) 0 (-(6*EI)/L2.^2) ((2*EI)/L2);(-EA/L2) 0 0 (EA/L2) 0 0;0 (-(12*EI)/L2.^3) (-(6*EI)/L2.^2) 0 ((12*EI)/L2.^3) (-(6*EI)/L2.^2);0 ((6*EI)/L2.^2) ((2*EI)/L2) 0 (-(6*EI)/L2.^2) ((4*EI)/L2)];
 Ke2=T2'*ke2*T2;
 
 T3=[c3 s3 0 0 0 0;-s3 c3 0 0 0 0;0 0 1 0 0 0;0 0 0 c3 s3 0;0 0 0 -s3 c3 0;0 0 0 0 0 1];
 ke3=[(EA/L3) 0 0 (-EA/L3) 0 0;0 ((12*EI)/L3.^3) ((6*EI)/L3.^2) 0 (-(12*EI)/L3.^3) ((6*EI)/L3.^2);0 ((6*EI)/L3.^2) ((4*EI)/L3) 0 (-(6*EI)/L3.^2) ((2*EI)/L3);(-EA/L3) 0 0 (EA/L3) 0 0;0 (-(12*EI)/L3.^3) (-(6*EI)/L3.^2) 0 ((12*EI)/L3.^3) (-(6*EI)/L3.^2);0 ((6*EI)/L3.^2) ((2*EI)/L3) 0 (-(6*EI)/L3.^2) ((4*EI)/L3)]; 
 Ke3=T3'*ke3*T3;
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    CALCULATION OF DISPLACEMENT IN GLOBAL COORDINATE
  a1=[1 2 3 4 5 6];
 for i= 1:length(a1);
     b1=a1(i);
     for j=1:length(a1)
         c2=a1(j);
         KG(b1,c2)=KG(b1,c2)+Ke1(i,j);
     end 
 end
 KG;
%  
 a2=[4 5 6 7 8 9];
 for i= 1:length(a2);
     b1=a2(i);
     for j=1:length(a2)
         c2=a2(j);
         KG(b1,c2)=KG(b1,c2)+Ke2(i,j);
     end 
 end
 KG;
%  
 a3=[7 8 9 10 11 12];
 for i= 1:length(a3);
     b1=a3(i);
     for j=1:length(a3)
         c2=a3(j);
         KG(b1,c2)=KG(b1,c2)+Ke3(i,j);
     end 
 end
 KG;
 

%%%%%%% Extracting elements from global stiffness matrix  %%%%%%
%%%%%%%% [case 1 - fixed end] %%%%%%%%%%
  Kg1= zeros(6,6); 
  a6=[4 5 6 7 8 9];
  for i=1:length(a6)
      b1=a6(i);
      for j=1:length(a6)
         c2=a6(j);
         Kg1(i,j)=KG(b1,c2)+Kg1(i,j);
      end 
  end
  Kg1;
%%%%%%%% [case 2 - pinned end] %%%%%%%%%%  
 Kg2= zeros(8,8); 
  a7=[3 4 5 6 7 8 9 12];
  for i=1:length(a7)
      b1=a7(i);
      for j=1:length(a7)
         c2=a7(j);
         Kg2(i,j)=KG(b1,c2)+Kg2(i,j);
      end 
  end
  Kg2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  D. OUTPUT OF RESULTS
%  %%%%%%% calculation of displacements %%%%%%%%%%%%%%%%
fprintf('displacements for fixed end with distributed load(a)')
UR1=inv(Kg1)*f11'        %%%% fixed end with distributed load
fprintf('displacements for fixed end with point load(b)')
UR2=inv(Kg1)*f21'        %%%% fixed end with point load 
fprintf('displacements for pinned end with distributed load(c)')
UR3=inv(Kg2)*f12'        %%%% Pinned end with distributed load
fprintf('displacements for pinned end with point load(d)')
UR4=inv(Kg2)*f22'        %%%% Pinned end with point load
