function [dm,th,Lagrange_coeff]=SDCcoeff(K,P)
if P==1
    lambda=1;
else
    lambda=JacobiGL(0,0,P);
end
dm=zeros(P,1);
if P==1
    dm=1;
else
    for i=1:P
        dm(i)=1/2*(lambda(i+1)-lambda(i));
    end
end
syms x y
th=zeros(P+1,1);
th(1)=0;
for i=1:P
    th(i+1)=th(i)+dm(i);
end
Lagrange_coeff=zeros(P+1,P);
for m=1:P
    for i=1:P+1
        Li=1;
        for j = 1:P+1
            if j ~= i
                Li = Li*(y - (x+th(j))) / (th(i) - th(j));
            end
        end
        Lagrange_coeff(i,m)=int(Li, y, x+th(m),x+th(m+1));
    end
end
end