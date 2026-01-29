% clean the enviroment
clear all;clc

% degree of polynomial
k=1;

% SDC
K=1;
P=1;
theta=1;

% spatial solution region [a,b]
coh=1;
a=0;b=coh*2*pi;

% set the Gauss integration accuracy(number of integration points)
NumGLP=5;  

% number of spatial subdivision 
% Np=[2,4,8,16,32,64,128,256,512,1024,2048]
% Np=[5,10,20,40,80,100,360,640,1280,2560]
Np=[5,10,20,40,80,160,320,640,1280,2560]
% spatial step size
h=(b-a)./Np; 

% the initial moment
t=0;

% the end moment
T=1; 
% intial time step size
dt=0.1;

% number of time subdivision 
TH=ceil(T./dt);

% corrected time step size
dt=T./TH;  

% set boundary condition
BC=2;

% penalty term parameters
k1=1;
k2=1;

iN=length(Np);
order_L2=zeros(1,iN);
for i=1:iN
        [L2(i)]= UWLDG(a,b,NumGLP,k,K,P,theta,Np(i),h(i),t,T,TH,dt,BC,k1,k2,i,iN,coh);
end
% loglog(Np,L2,'b:*') % exact solution
% hold on
semilogy(Np,L2,'r:s','LineWidth', 2.2) % exact solution
xlabel('number of N','FontSize', 14);
ylabel('L2 error','FontSize', 12);
xlim([0 2650])

% hold on
% semilogy(Np([2:20:iN]), L2([2:20:iN]),'bd','LineWidth',1,'MarkerSize',5);

% ylim([0,L2(end)]);  % y轴从0开始




function [L1,L2,L_inf]=UWLDG(a,b,NumGLP,k,K,P,theta,Np,h,t,T,TH,dt,BC,k1,k2,Npth,iN,coh)



% calculate the basic function
statement; 

% all midpoints of the spatial subdivision units
X = a+h/2 : h : b-h/2;  

% calculate the total stiffness matrix
Stiffness 

% calculate the initial value
u=initial_value(a,b,h,Np,k,NumGLP,weight,lambda,Ma,phiG,coh); 

% time marching
% SDC_3
NEWSDC

% calculate the L2 error
[L1,L2,L_inf] = L_error(u,h,X,T,Np,k,NumGLP,weight,lambda,phiG,coh);

end