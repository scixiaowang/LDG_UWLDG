function y=legend2(k,z)
if k==0
    y=0;
elseif k==1
    y=0;
elseif k==2
    y=3;
elseif k==3
    y=15*z;
elseif k==4
    y=(105*z^2)/2 - 15/2;
elseif k==5
    y=(315*z^3)/2 - (105*z)/2;
end
return