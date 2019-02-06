function [area1]=fvcom_area(side1,side2,side3)
% function [area1]=fvcom_area(side1,side2,side3)
%  calculate the area of a triangle on a spherical plane
%           on earth assuming earth radius to be  rearth= 6371.0E+3 (meter)
% input:
%     side1,side2 and side3: are 3 arc lenth for one triangle
% output:
%     areal: is area of a triangle on a spherical plane
%
% Code and Documentation:
%
%   Wen Long@PNNL  03/28/2013, Seattle
%
%   Note: adopted from mod_spherical.F in FVCOM model
%
   if(nargin<3)
     error(['oops incorrect number of arguments']);
   end
   rearth= 6371.0E+3;   
   side1=side1/rearth;
   side2=side2/rearth;
   side3=side3/rearth;
   if(side1 == 0. || side2 == 0. || side3 == 0.)
     area1=0.
   else
     psum=0.5*(side1+side2+side3);
     pm=sin(psum)*sin(psum-side1)*sin(psum-side2)*sin(psum-side3);
     pm=sqrt(pm)/(2.0*cos(side1*0.5)*cos(side2*0.5)*cos(side3*0.5));
     qmjc = 2.0*asin(pm);

     area1=rearth*rearth*qmjc;

   end

end
