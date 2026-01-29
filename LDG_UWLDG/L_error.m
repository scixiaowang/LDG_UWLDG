function [L1,L2,L_inf] = L_error(u,h,X,T,Np,k,NumGLP,weight,lambda,phiG)

uhG = zeros(NumGLP,Np);
ureal = zeros(NumGLP,Np);
uu=zeros(Np,1);

for i = 1:Np
    for j = 1:NumGLP
        ureal(j,i) = exact_sol(X(i) + h/2*lambda(j),T);
    end
end

for i = 1:Np
    for d = 1:k+1
        for j=1:NumGLP
            uhG(j,i) = uhG(j,i) + u((i-1)*(k+1)+d,1)*phiG(d,j,1);
        end
    end
end

uE=abs(uhG-ureal);


L1=0;
for i = 1:Np
    for i1 = 1:NumGLP
        L1 = L1 + h/2*weight(i1)*uE(i1,i);
    end
end

L2=0;
for i = 1:Np
    for i1 = 1:NumGLP
        L2 = L2 + h/2*weight(i1)*uE(i1,i)^2;
    end
end
L2=sqrt(L2);


for i=1:Np
    for j=1:k+1
            uu(i)=uu(i)+u((i-1)*(k+1)+j,1)*phiG(j,(NumGLP+1)/2,1);
    end
end
L_inf=max(uu-exact_sol(X,T)');

end