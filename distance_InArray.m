function [dlat dlon] = distance_InArray(cent_stla,cent_stlo,stla, stlo)
% This function is used to calculate the distance between stations and reference location, 
% (usually the center of the array) after rotating whole array to equator while refenence 
% location at lat = 0 and long = 0. 
% by Pei-Ying Patty Lin 07/31/2014 

%% =================================================
% determine rotation tensor based on the reference point
%===================================================

phi = deg2rad(cent_stlo);
theta = deg2rad(cent_stla);


x(1)=cos(theta)*cos(phi);
x(2)=cos(theta)*sin(phi);
x(3)=sin(theta);

cphi=pi/2.0-phi;
aphi=pi/2.0+phi;

t(1,1)=cos(phi);
t(1,2)=cos(cphi);
t(1,3)=0.;
t(2,1)=cos(aphi);
t(2,2)=t(1,1);
t(2,3)=0.;
t(3,1)=0.;
t(3,2)=0.;
t(3,3)=1.0;



%% =================================================
%  use stla stlo to do rotation using the tensor
%===================================================


for ista = 1: length(stla)
    phi=deg2rad(stlo(ista));
    theta=deg2rad(stla(ista));


    x(1)= cos(theta)*cos(phi);
    x(2)= cos(theta)*sin(phi);
    x(3)= sin(theta);


    xp = t*x';
    y1 = sqrt(xp(1)*xp(1)+xp(2)*xp(2));

    if (abs(y1) > 1 ) 
       display('program sets abs(arg)=1, and continues ....');
       if ( y1 < -1 ) y1 = -1.0;
       end
       if ( y1 > 1 ) y1 = 1.0;
       end
    end
    lat1 = acos(y1);
    if (xp(3) < 0)
        lat1=-lat1;
    end

    y=xp(1)/y1;
    if (abs(y) > 1 ) 
       display('program sets abs(arg)=1, and continues ....');
       if ( y < -1 ) y = -1.0;
       end
       if ( y > 1 ) y = 1.0;
       end
    end
    lon1=acos(y);
    if (xp(2) < 0) 
        lon1 = -lon1;
    end

    st_la(ista) = rad2deg(lat1);
    st_lo(ista) = rad2deg(lon1);

end

%% =================================================
%  rotate the reference point to Equator about x2 and others accordingly
%===================================================

    phi = deg2rad(cent_stlo);
    theta = deg2rad(cent_stla);


    x(1)=cos(theta)*cos(phi);
    x(2)=cos(theta)*sin(phi);
    x(3)=sin(theta);


    cthe=pi/2.0-theta;
    athe=pi/2.0+theta;

    t1(1,1)=cos(theta);
    t1(1,2)=0.;
    t1(1,3)=cos(cthe);
    t1(2,1)=0.;
    t1(2,2)=1.;
    t1(2,3)=0.;
    t1(3,1)=cos(athe);
    t1(3,2)=0.;
    t1(3,3)=t1(1,1);


    %% =================================================
    %  use st_la st_lo  to do rotation using the tensor
    %===================================================
for ista = 1 : length(stla)
    
    phi=deg2rad(st_lo(ista));
    theta=deg2rad(st_la(ista));

    x(1)=cos(theta)*cos(phi);
    x(2)=cos(theta)*sin(phi);
    x(3)=sin(theta);


    xp = t1*x';

    % ------ recover new lat,lon ----
    y1=sqrt(xp(1)*xp(1)+xp(2)*xp(2));
    if (abs(y1) > 1 ) 
       display('program sets abs(arg)=1, and continues ....');
       if ( y1 < -1 ) y1 = -1.0;
       end
       if ( y1 > 1 ) y1 = 1.0;
       end
    end

    lat1=acos(y1);
    if (xp(3) < 0)
        lat1=-lat1;
    end


    y=xp(1)/y1;
    if (abs(y) > 1 ) 
       display('program sets abs(arg)=1, and continues ....');
       if ( y < -1 ) y = -1.0;
       end
       if ( y > 1 ) y = 1.0;
       end
    end

    lon1=acos(y);
    if (xp(2) < 0) 
        lon1 = -lon1;
    end

    dlat(ista) = rad2deg(lat1);
    dlon(ista) = rad2deg(lon1);
end
return
end







