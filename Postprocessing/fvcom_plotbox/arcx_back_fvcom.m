   function [arcx1]= arcx_back_fvcom(xx1,yy1,xx2,yy2)
%
% function [arcx1]= arcx_back_fvcom(xx1,yy1,xx2,yy2)
%
% calculate the x-arm length of an arc defined 
% by two points 
   pi=3.14159265; 
   deg2rad=pi/180.0; 

   i=0; 
   nx=500;
   xx1=0;
   yy1=0;
   xx2=0;
   yy2=0;
   arcx1=0;
   x1=0;
   y1=0;
   x2=0;
   y2=0;
   ty=0;
   a1=0;
   a2=0;
   b1=0;
   b2=0;
   c1=0;
   c2=0;
   a =0;
   b =0;
   c =0;
   x=zeros([(nx+1),1]);
   y=zeros([(nx+1),1]); 
   xtmp=0;
   if(xx1 == xx2)
     arcx1=0.;
   else
     x1=xx1*deg2rad;
     y1=yy1*deg2rad;
     x2=xx2*deg2rad;
     y2=yy2*deg2rad;
     x(1)=x1;
     y(1)=y1;
     x(nx+1)=x2;
     y(nx+1)=y2;
     xtmp=x(nx+1)-x(1);
     if(xtmp >  pi)
       xtmp = -2*pi+xtmp;
     else if(xtmp < -pi)
       xtmp =  2*pi+xtmp;
     end
     for i=2:nx
       x(i)=x(i-1)+xtmp/float(nx);
     end
     a1=cos(y(1))*cos(x(1))
     a2=cos(y(nx+1))*cos(x(nx+1));
     b1=cos(y(1))*sin(x(1));
     b2=cos(y(nx+1))*sin(x(nx+1));
     c1=sin(y(1));
     c2=sin(y(nx+1));
     a=a1*b2-a2*b1;
     b=b1*c2-b2*c1;
     c=a2*c1-a1*c2;
     for i=2:nx
       y(i)=-b*cos(x(i))-c*sin(x(i));
       y(i)=y(i)/a;
       y(i)=atan(y(i));
     end
     arcx1=0.;
     for i=1:nx
       ty=0.5*(y(i)+y(i+1));
       xtmp=x(i+1)-x(i);
       if(xtmp >  pi)
         xtmp = -2*pi+xtmp;
       else if(xtmp < -pi)
         xtmp =  2*pi+xtmp;
       end
       arcx1=arcx1+rearth*cos(ty)*xtmp;
%       arcx1=arcx1+rearth*cos(ty)*(x(i+1)-x(i))
     end
   end

