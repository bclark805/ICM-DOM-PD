function [aw0,awx,awy]=check_shape_continuity(vx,vy,nv,f)
% function [aw0,awx,awy]=check_shape_continuity(vx,vy,nv,f)
%
%   Calculate the coefficients for a linear function in x-y plane with 3 node values known for a triangular grid.
%   Then do the same for a neighboring triangular grid, and make check if the functions are continuous on the shared edge.
%
%   Formulation:
%
%     	    For a triangle with 3 points located at (x,y)=( vx(1:3), vy(1:3) ), the corresponding values are given as 
%        f(1:3).  A linear function f(x,y) can be constructed based on 
%		f(x',y') = fi + df/dx * x' + df/dy * y'
%        where x'=x-x0, y'=y-y0, with x0,y0 being the centroid of the triangle, x0=mean(vx), y0=mean(vy).
%        fi,df/dx,df/dy are coefficients to be calculated here. 
%
%         fi= aw0(i,1)*f1+aw0(i,2)*f2+aw0(i,3)*f3
%	      aw0(1:n,1:3)  --- coefficients of estimating ceontroid value based on node values of an element
%                               where fi is value at centroid, f1,f2,f3 are values on nodes of a triangular element i
%
%         df/dx = awx(i,1)*f1+awx(i,2)*f2+awx(i,3)*f3
%             awx(1:n,1:3)  --- coefficients of estimating x gradient of function f within an element based
%                               on nodal values of the element
%                               where f1,f2,f3 are values on the nodes of an element i
%
%         df/dy = awx(i,1)*f1+awx(i,2)*f2+awx(i,3)*f3
%             awy(1:n,1:3)  --- coefficients of estimating y gradient of function f within an element based
%                               on nodal values of the element
%                               where f1,f2,f3 are values on the nodes of an element i
%
%   Inputs:
%        vx(1:m)    --- x coordinate of all nodes
%        vy(1:m)    --- y coordinate of all nodes
%        nv(1:n,1:3)--- element to node connectivity (must be clockwise)
%        f(1:n)     --- function values of the nodes
%
%   Outputs: 
%
%       aw0(1:n,1:3)  --- coefficients of estimating ceontroid value based on node values of an element
%                         fi= aw0(i,1)*f1+aw0(i,2)*f2+aw0(i,3)*f3
%                         where fi is value at centroid, f1,f2,f3 are values on nodes of a triangular element i
%       awx(1:n,1:3)  --- coefficients of estimating x gradient of function f within an element based 
%                         on nodal values of the element
%                         df/dx = awx(i,1)*f1+awx(i,2)*f2+awx(i,3)*f3
%                         where f1,f2,f3 are values on the nodes of an element i
%       awy(1:n,1:3)  --- coefficients of estimating y gradient of function f within an element based
%                         on nodal values of the element
%                         df/dy = awx(i,1)*f1+awx(i,2)*f2+awx(i,3)*f3
%                         where f1,f2,f3 are values on the nodes of an element i
%    Example:
%	vx = [0, 1, 1.5, 0.5];
%       vy = [0, 0, 1  , 1  ];
%       nv = [[1, 4, 2], [4, 3, 2]];  %nodes 2,4 are on a shared edge
%        f = [0, 5, 10, 7];
%       check_shape_coninuity(vx,vy,nv,f); 
%
%    Code and documentation
%         Wen Long, PNNL/BSRC @Seattle, 08/20/2013
%

%        xc(1:n)    --- x coord of centroid of elements
%        yc(1:n)    --- y coord of centroid of elements
%        art(1:n)   --- area (m^2) of all elmenents

if nargin <4   %if number of arguments less than 4
               %simply quit and give help
     [ST,I] = dbstack;
     PRETTY_FUNCTION = ST.name;  %get name of current function
     display(PRETTY_FUNCTION)
     help(PRETTY_FUNCTION);
     s=-1;
     return
end

    x1=0;
    x2=0;
    x3=0; 
    y1=0;
    y2=0;
    y3=0;
    ai1=0;
    ai2=0;
    ai3=0;
    bi1=0;
    bi2=0;
    bi3=0;
    ci1=0;
    ci2=0;
    ci3=0;
    i=0;
    j=0;

   n=size(nv,1);      %number of elements
   m=length(vx(:));   %number of nodes
   mf=length(f(:)); 
   mvy=length(vy(:));

   if(n~=2 || m~=4)
     error('oops, this testing program only takes two triangular elments, and totally 4 nodes, one edge should be shared')
   end

   if(mf<m)
    error('oops, array f is too small'); 
   end
   if(mvy~=m)
    error('oops, wrong vy size');
   end

   %find the shared edge, a shared edge has two nodes
    ishare=intersect(nv(1,:),nv(2,:));
   if(length(ishare(:))~=2)
      error('oops, the two elements do not have a shared edge');
   end

   xc=zeros([2,1]);
   yc=zeros([2,1]); 
   art=zeros([2,1]);

   for i=1:n

     %calculate centroid value
     xc(i)=mean(vx(nv(i,:)));
     yc(i)=mean(vy(nv(i,:)));

     %calculate area of element
     art(i) = (vx(nv(i,2)) - vx(nv(i,1))) * (vy(nv(i,3)) - vy(nv(i,1))) -  ...
              (vx(nv(i,3)) - vx(nv(i,1))) * (vy(nv(i,2)) - vy(nv(i,1)));
     art(i) = abs(.5*art(i));

     if(art(i)<=0)
        error(['oops, element i=' int2str(i) ' has zero area']);
     end

   end

   for i=1:n 

     x1=vx(nv(i,1))-xc(i);
     x2=vx(nv(i,2))-xc(i);
     x3=vx(nv(i,3))-xc(i);
     y1=vy(nv(i,1))-yc(i);
     y2=vy(nv(i,2))-yc(i);
     y3=vy(nv(i,3))-yc(i);

     ai1=y2-y3;
     ai2=y3-y1;
     ai3=y1-y2;

     bi1=x3-x2;
     bi2=x1-x3;
     bi3=x2-x1;

     ci1=x2*y3-x3*y2;
     ci2=x3*y1-x1*y3;
     ci3=x1*y2-x2*y1;

     aw0(i,1)=-ci1/2./art(i);
     aw0(i,2)=-ci2/2./art(i);
     aw0(i,3)=-ci3/2./art(i);

     awx(i,1)=-ai1/2./art(i);
     awx(i,2)=-ai2/2./art(i);
     awx(i,3)=-ai3/2./art(i);

     awy(i,1)=-bi1/2./art(i);
     awy(i,2)=-bi2/2./art(i);
     awy(i,3)=-bi3/2./art(i);

  end


  %check whether the reconstructed fucntion value on the shared edge using forumations of each neighboring element is the same (continuous)
 
   %find length of the shared edge

        L=sqrt( (vx(ishare(1))-vx(ishare(2)))^2 ...
	       +(vy(ishare(1))-vy(ishare(2)))^2 ...
	      );
   %find angle of the shared edge
       angle=atan2(vy(ishare(2))-vy(ishare(1)),vx(ishare(2))-vx(ishare(1)));
    
   %sample the shared edge into 100 points and calculat f(x,y) on those 100 points and compare the calculation use coeffcients from each 
   %neighboring element.
       ds=L/100.0;

       s=0:ds:L;

       %find the sampling points on the shared edge

       y=vy(ishare(1))+sin(angle)*s;
       x=vx(ishare(1))+cos(angle)*s;

       %find value of point (x,y) using coefficients of first element

       ie=1;
       xprime=x-xc(ie);
       yprime=y-yc(ie);
       f1=  aw0(ie,1)*f(nv(ie,1)) + aw0(ie,2)*f(nv(ie,2))+aw0(ie,3)*f(nv(ie,3)) ...
         + (awx(ie,1)*f(nv(ie,1)) + awx(ie,2)*f(nv(ie,2))+awx(ie,3)*f(nv(ie,3))).*xprime ...
         + (awy(ie,1)*f(nv(ie,1)) + awy(ie,2)*f(nv(ie,2))+awy(ie,3)*f(nv(ie,3))).*yprime;


       ie=2;
       %find value of point (x,y) using coefficients of second element
       xprime=x-xc(ie);
       yprime=y-yc(ie);
       f2=  aw0(ie,1)*f(nv(ie,1)) + aw0(ie,2)*f(nv(ie,2))+aw0(ie,3)*f(nv(ie,3)) ...
         + (awx(ie,1)*f(nv(ie,1)) + awx(ie,2)*f(nv(ie,2))+awx(ie,3)*f(nv(ie,3))).*xprime ...
         + (awy(ie,1)*f(nv(ie,1)) + awy(ie,2)*f(nv(ie,2))+awy(ie,3)*f(nv(ie,3))).*yprime;
         

       %plot f1(s), f2(s) to show if they are the same

       figure(1); hold on;
       h1=plot(s,f1+1,'r-.');
       h2=plot(s,f2,'g-.');
       h3=plot(s,f1-f2,'k-');
       h4=plot([0,L],[f(ishare(1)),f(ishare(2))],'mo','markersize',15);

       legend(gca,[h1,h2,h3,h4],'interpolation based on element 1 with offset 1','interpolation based on element 2','difference','original');
       title(['Interpolated values on the shared edge with nodes [' int2str(ishare(1)) ',' int2str(ishare(2)) ']']);

end %end of function

