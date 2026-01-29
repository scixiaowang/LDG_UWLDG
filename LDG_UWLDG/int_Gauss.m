function y=int_Gauss(u,v,h,NumGLP,weight)
y=0;
for i = 1:NumGLP
       y =y + 0.5*h*weight(i)*u(i)*v(i);
end

end