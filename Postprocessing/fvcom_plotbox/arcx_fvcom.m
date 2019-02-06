function arcx1=arcx_fvcom(xx1,yy1,xx2,yy2,arcx1)
%
%function arcx1= arcx_fvcom(xx1,yy1,xx2,yy2)
%returns deltax (arc length) of a great arc of points (lon=xx1,lat=yy1) to (lon=xx2,lat=yy2) 
%projected on equatorial plane on earth with earth radius assumed to be 6371.0E+3 meter
%
% input:
%   xx1,yy1 --- lon, lat of point 1 on sphere (deg)
%   xx2,yy2 --- lon, lat of point 2 on sphere (deg)
%
% output:
%  arcx1 --- arc length of great arc projected on equatorial plane
%  
   pi=3.14159265;
   deg2rad=pi/180.0;  
   rearth= 6371.0E+3;   %radius of earth (meter)

   if(xx1 == xx2)  %if same longitude 
     arcx1=0.0;       %then no x coordinate difference
   else

     x1=xx1*deg2rad;  % longitude of point 1
     y1=yy1*deg2rad;  % latitude of point 1

     x2=xx2*deg2rad;  % longitude of point 2
     y2=yy2*deg2rad;  % latitude of point 2

     xtmp  = x2-x1;    %angle difference of longitude

     if(xtmp >  pi)
       xtmp = -2*pi+xtmp;  %makes sure angle difference is within [-pi,pi]
     elseif(xtmp < -pi)
       xtmp =  2*pi+xtmp;
     end 

     ty=0.5*(y2+y1);       %mean latitude
     arcx1=rearth*cos(ty)*xtmp;  
   end
   
end 


