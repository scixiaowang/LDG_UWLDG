[dm,th,Lagrange_coeff]=SDCcoeff(K,P);
th = th(end:-1:1);
LUpre=zeros((k+1)*Np,(k+1)*Np,P);
for i=1:P
As=eye((k+1)*Np)-dm(i)*dt*C;
[Ls, Us, Ps] = lu(sparse(As));
LUpre(:,:,i)=Us\(Ls\Ps);
end

LUcor=zeros((k+1)*Np,(k+1)*Np,P);
for i=1:P
As=eye((k+1)*Np)-dm(i)*dt*theta*C;
[Ls, Us, Ps] = lu(sparse(As));
LUcor(:,:,i)=Us\(Ls\Ps);
end

ut=zeros(Np*(k+1),K+1,P+1);
bu=zeros(Np*(k+1),K+1,P+1);

for i=1:TH
    t=t+dt;
    for l=1:K+1
        ut(:,l,1)=u;
    end
%pre-step
    for m=1:P
        bu(:,1,m+1)=ut(:,1,m);
        ut(:,1,m+1)=LUpre(:,:,m)*bu(:,1,m+1);
    end
%cor-step
    for l=1:K
        for m=1:P
            int=0;
            for li=1:P+1
                int = int + Lagrange_coeff(li,m)*dt*C*ut(:,l,li);
            end
            bu(:,l+1,m+1) = ut(:,l+1,m) - theta*dm(m)*dt*(C*ut(:,l,m+1))+ int;
            ut(:,l+1,m+1) = LUcor(:,:,m)*bu(:,l+1,m+1);
        end
    end
    u=ut(:,end,end);

    clc
    fprintf('Current progress: %s  %%\n',num2str(i/TH*100))
    fprintf('Current number of grids: %i \n',Np)
    fprintf('Total process: %i/%i \n',Npth,iN)
end
clc