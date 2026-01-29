function y=legend1(k,z)
if k==0
    y=0;
elseif k==1
    y=1;
elseif k==2
    y=3*z;
elseif k==3
    y=(15*z^2)/2 - 3/2;
elseif k==4
    y=(35*z^3)/2 - (15*z)/2;
elseif k==5
    y=(315*z^4)/8 - (105*z^2)/4 + 15/8;
end
return