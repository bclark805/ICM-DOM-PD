 function arc1=arc_fvcom(xx1,yy1,xx2,yy2)
%
% function arc1=arc_fvcom(xx1,yy1,xx2,yy2)
% calculate the great arc lenth for given two point on the spherical plane
% Assuming sphere radius being 6371.0E+3 meter 
%  input:
%           xx1,yy1,xx2,yy2 :are longitude and latitude of two points (deg)
%  output:
%           arcl :  arc lenth of two points in spherical plane (meter)
%

   pi=3.14159265;
   deg2rad=pi/180.0; 
   rearth=6371.0E+3;   

   x1=xx1*deg2rad;
   y1=yy1*deg2rad;

   x2=xx2*deg2rad;
   y2=yy2*deg2rad;

   xa=cos(y1)*cos(x1);
   ya=cos(y1)*sin(x1);
   za=sin(y1);

   xb=cos(y2)*cos(x2);
   yb=cos(y2)*sin(x2);
   zb=sin(y2);

   ab=sqrt((xb-xa)**2+(yb-ya)**2+(zb-za)**2);
   aob=(2.-ab*ab)/2.;
   aob=acos(aob);
   arcl=rearth*aob;

end

