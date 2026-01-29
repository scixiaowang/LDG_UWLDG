% This program is used to calculate the order of convergence.

% clean the enviroment
clear all;clc

% degree of polynomial
k=3;

% SDC parameters
K=3;
P=3;
theta=1;

% spatial solution region [a,b]
a=0;b=2*pi;

% set the Gauss integration accuracy(number of integration points)
NumGLP=5;  

% number of spatial subdivision 
Np=[10,20,40,80,160];

% spatial step size
h=(b-a)./Np; 

% temporal discretization parameters
t=0;
T=1; 
CFL=1;
dt=CFL*h;
TH=ceil(T./dt);
dt=T./TH;  


iN=length(Np);
if iN==1
    Npth=1;
end

% [L1,L2,L_inf]=UWLDG(a,b,NumGLP,k,K,P,theta,Np,h,t,T,TH,dt,Npth,iN)
% Calculate the error
for i=1:iN
        [UWLDG_L1(i),UWLDG_L2(i),UWLDG_L_inf(i)]= UWLDG(a,b,NumGLP,k,K,P,theta,Np(i),h(i),t,T,TH(i),dt(i),i,iN);
        [LDG_L1(i),LDG_L2(i),LDG_L_inf(i)]= LDG(a,b,NumGLP,k,K,P,theta,Np(i),h(i),t,T,TH(i),dt(i),i,iN);
end

% Calculate the order
for i=2:iN
    UWLDG_order_L1(i)=log(UWLDG_L1(i-1)/UWLDG_L1(i))/log(2);
    UWLDG_order_L2(i)=log(UWLDG_L2(i-1)/UWLDG_L2(i))/log(2);
    UWLDG_order_L_inf(i)=log(UWLDG_L_inf(i-1)/UWLDG_L_inf(i))/log(2);
    LDG_order_L1(i)=log(LDG_L1(i-1)/LDG_L1(i))/log(2);
    LDG_order_L2(i)=log(LDG_L2(i-1)/LDG_L2(i))/log(2);
    LDG_order_L_inf(i)=log(LDG_L_inf(i-1)/LDG_L_inf(i))/log(2);
end



fprintf('LDG_L2 =  ');
fprintf(' %.2e  ', LDG_L2);
fprintf('\n');
fprintf('UWLDG_L2 =  ');
fprintf(' %.2e  ', UWLDG_L2);
fprintf('\n\n');

fprintf('LDG_order_L2 =  ');
fprintf(' %.2f  ', LDG_order_L2);
fprintf('\n');
fprintf('UWLDG_order_L2 =  ');
fprintf(' %.2f  ', UWLDG_order_L2);
fprintf('\n');

loglog(h,LDG_L2,'bo-', 'LineWidth',1, 'MarkerSize',6);
hold on
loglog(h,UWLDG_L2,'r*-', 'LineWidth',1, 'MarkerSize',6);

hold on
wa1=zeros(1,iN);
wa2=zeros(1,iN);
wa1(1)=UWLDG_L2(1);
wa2(1)=LDG_L2(1);
for i=2:iN
    wa1(i)=2^(k+1)\wa1(i-1);
end
for i=2:iN
    wa2(i)=2^(k+1)\wa2(i-1);
end
loglog(h,wa1,'k--');
hold on
loglog(h,wa2,'k--');
xlabel('h');
ylabel('L2 error');
title('LDG vs UWLDG log-log plot');
legend('LDG method','UWLDG method','benchmark','Location','northwest');
