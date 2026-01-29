function [L1,L2,L_inf,u]=UWLDG(a,b,NumGLP,k,K,P,theta,Np,h,t,T,TH,dt,Npth,iN)

% calculate the basic function
statement; 

% all midpoints of the spatial subdivision units
X = a+h/2 : h : b-h/2;  

%% calculate the total stiffness matrix
% construct “assembly matrix”
%j=2,3,...,N-1
Qj_int=diag([0,ones(1,Np-2),0]);
Qj_intR=diag([0,ones(1,Np-2)],1);
Qj_intL=diag([ones(1,Np-2),0],-1);
%j=1
Q1_bry=zeros(Np,Np);
Q1_bry(1,1)=1;
Q1_bryR=zeros(Np,Np);
Q1_bryR(1,2)=1;
Q1_bryL=zeros(Np,Np);
Q1_bryL(1,Np)=1;
%j=N
QN_bry=zeros(Np,Np);
QN_bry(end,end)=1;
QN_bryR=zeros(Np,Np);
QN_bryR(Np,1)=1;
QN_bryL=zeros(Np,Np);
QN_bryL(end,end-1)=1;


% Local mass matrix
M=zeros(k+1,k+1);
for i=1:k+1
    M(i,i)=h/(1+2*(i-1));
end
% Global mass matrix
Ma=kron(Qj_int+Q1_bry+QN_bry,M);


S_int00=zeros(k+1,k+1);
S_int01=zeros(k+1,k+1);
S_int02=zeros(k+1,k+1);

S_LlRl=zeros(k+1,k+1);%v^{-}_{j-1/2}*v^{+}_{j-1/2}
S_LrLr=zeros(k+1,k+1);%v^{-}_{j+1/2}*v^{-}_{j+1/2}
S_LlxRl=zeros(k+1,k+1);%v^{-}_{j-1/2}*vx^{+}_{j-1/2}
S_LrxLr=zeros(k+1,k+1);%v^{-}_{j+1/2}*vx^{-}_{j+1/2}
S_xLlRl=zeros(k+1,k+1);%vx^{-}_{j-1/2}*v^{+}_{j-1/2}
S_xLrLr=zeros(k+1,k+1);%vx^{-}_{j+1/2}*v^{-}_{j+1/2}

S_RlRl=zeros(k+1,k+1);%v^{+}_{j-1/2}*v^{+}_{j-1/2}
S_RrLr=zeros(k+1,k+1);%v^{+}_{j+1/2}*v^{-}_{j+1/2}
S_RlxRl=zeros(k+1,k+1);%v^{+}_{j-1/2}*vx^{+}_{j-1/2}
S_RrxLr=zeros(k+1,k+1);%v^{+}_{j+1/2}*vx^{-}_{j+1/2}
S_xRlRl=zeros(k+1,k+1);%vx^{+}_{j-1/2}*v^{+}_{j-1/2}
S_xRrLr=zeros(k+1,k+1);%vx^{+}_{j+1/2}*v^{-}_{j+1/2}
S_xLrxLr=zeros(k+1,k+1);%vx^{+}_{j-1/2}*vx^{+}_{j-1/2}
S_xRlxRl=zeros(k+1,k+1);%vx^{-}_{j+1/2}*vx^{-}_{j+1/2}

% The matrix for integral terms and flux terms
for i=1:k+1
    for j=1:k+1
        S_int00(i,j) = int_Gauss(phiG(j,:,1),phiG(i,:,1),h,NumGLP,weight);
        S_int01(i,j) = int_Gauss(phiG(j,:,1),phiG(i,:,2),h,NumGLP,weight);
        S_int02(i,j) = int_Gauss(phiG(j,:,1),phiG(i,:,3),h,NumGLP,weight);
        
        S_LlRl(i,j)  = phiGL(j,1)*phiGR(i,1);
        S_LlxRl(i,j) = phiGL(j,1)*phiGR(i,2);
        S_xLlRl(i,j) = phiGL(j,2)*phiGR(i,1);
        S_RlRl(i,j)  = phiGR(j,1)*phiGR(i,1);
        S_RlxRl(i,j) = phiGR(j,1)*phiGR(i,2);
        S_xRlRl(i,j) = phiGR(j,2)*phiGR(i,1);
        S_xRlxRl(i,j)= phiGR(j,2)*phiGR(i,2);

        S_LrLr(i,j)  = phiGL(j,1)*phiGL(i,1);
        S_LrxLr(i,j) = phiGL(j,1)*phiGL(i,2);
        S_xLrLr(i,j) = phiGL(j,2)*phiGL(i,1);
        S_RrLr(i,j)  = phiGR(j,1)*phiGL(i,1);
        S_RrxLr(i,j) = phiGR(j,1)*phiGL(i,2);
        S_xRrLr(i,j) = phiGR(j,2)*phiGL(i,1);
        S_xLrxLr(i,j)= phiGL(j,2)*phiGL(i,2);
    end
end

% Select the flux
aa=1;
bb=1;
aa1=1-aa;
bb1=1-bb;

% Global stiffness matrix 1
%j=2,3,...,N-1
SP=kron(Qj_int,S_int02)...
    +bb*kron(Qj_int,S_xLrLr)+bb1*kron(Qj_intR,S_xRrLr)...
    -bb*kron(Qj_intL,S_xLlRl)-bb1*kron(Qj_int,S_xRlRl)...
    -aa*kron(Qj_int,S_LrxLr)-aa1*kron(Qj_intR,S_RrxLr)...
    +aa*kron(Qj_intL,S_LlxRl)+aa1*kron(Qj_int,S_RlxRl);

%j=1
SP=SP+kron(Q1_bry,S_int02)...
    +bb*kron(Q1_bry,S_xLrLr)+bb1*kron(Q1_bryR,S_xRrLr)...
    -bb*kron(Q1_bryL,S_xLlRl)-bb1*kron(Q1_bry,S_xRlRl)...
    -aa*kron(Q1_bry,S_LrxLr)-aa1*kron(Q1_bryR,S_RrxLr)...
    +aa*kron(Q1_bryL,S_LlxRl)+aa1*kron(Q1_bry,S_RlxRl);
%j=N
SP=SP+kron(QN_bry,S_int02)...
    +bb*kron(QN_bry,S_xLrLr)+bb1*kron(QN_bryR,S_xRrLr)...
    -bb*kron(QN_bryL,S_xLlRl)-bb1*kron(QN_bry,S_xRlRl)...
    -aa*kron(QN_bry,S_LrxLr)-aa1*kron(QN_bryR,S_RrxLr)...
    +aa*kron(QN_bryL,S_LlxRl)+aa1*kron(QN_bry,S_RlxRl);
% Iteration matrix of auxiliary variables1
CP=Ma\SP; 

% Global stiffness matrix 2
%j=2,3,...,N-1
SUt=-kron(Qj_int,S_int02)...
    -bb1*kron(Qj_int,S_xLrLr)-bb*kron(Qj_intR,S_xRrLr)...
    +bb1*kron(Qj_intL,S_xLlRl)+bb*kron(Qj_int,S_xRlRl)...
    +aa1*kron(Qj_int,S_LrxLr)+aa*kron(Qj_intR,S_RrxLr)...
    -aa1*kron(Qj_intL,S_LlxRl)-aa*kron(Qj_int,S_RlxRl);
%j=1
SUt=SUt-kron(Q1_bry,S_int02)...
    -bb1*kron(Q1_bry,S_xLrLr)-bb*kron(Q1_bryR,S_xRrLr)...
    +bb1*kron(Q1_bryL,S_xLlRl)+bb*kron(Q1_bry,S_xRlRl)...
    +aa1*kron(Q1_bry,S_LrxLr)+aa*kron(Q1_bryR,S_RrxLr)...
    -aa1*kron(Q1_bryL,S_LlxRl)-aa*kron(Q1_bry,S_RlxRl);
%j=N
SUt=SUt-kron(QN_bry,S_int02)...
    -bb1*kron(QN_bry,S_xLrLr)-bb*kron(QN_bryR,S_xRrLr)...
    +bb1*kron(QN_bryL,S_xLlRl)+bb*kron(QN_bry,S_xRlRl)...
    +aa1*kron(QN_bry,S_LrxLr)+aa*kron(QN_bryR,S_RrxLr)...
    -aa1*kron(QN_bryL,S_LlxRl)-aa*kron(QN_bry,S_RlxRl);
% Iteration matrix of auxiliary variables2
CUt=Ma\SUt;

% Total iteration matrix
C=CUt*CP;

%% calculate the initial value
u=initial_value(a,b,h,Np,k,NumGLP,weight,lambda,Ma,phiG); 

%% time marching
SDC

%% calculate the L2 error
[L1,L2,L_inf] = L_error(u,h,X,T,Np,k,NumGLP,weight,lambda,phiG);

end