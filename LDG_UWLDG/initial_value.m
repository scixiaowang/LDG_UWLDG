function uh = initial_value(a,b,h,Np,k,NumGLP,weight,lambda,Ma,phiG)
uh = zeros((k+1)*Np,1);
ureal = zeros(Np,NumGLP);
X = a+h/2 : h : b-h/2;

for i = 1:Np
    for j = 1:NumGLP
        ureal(i,j) = exact_sol(X(i) + 0.5*h*lambda(j),0);
    end
end

for i = 1:Np
    for d = 1:k+1
        uh((k+1)*(i-1)+d,1) = int_Gauss(ureal(i,:),phiG(d,:,1),h,NumGLP,weight);
    end
end


uh = Ma\uh;

end