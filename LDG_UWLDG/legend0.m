function y=legend0(k,z)
if k==0
    y=1;
elseif k==1
    y=z;
elseif k==2
    y=(3*z^2)/2 - 1/2;
elseif k==3
    y=(5*z^3)/2 - (3*z)/2;
elseif k==4
    y=(35*z^4)/8 - (15*z^2)/4 + 3/8;
elseif k==5
    y=(63*z^5)/8 - (35*z^3)/4 + (15*z)/8;
end
return