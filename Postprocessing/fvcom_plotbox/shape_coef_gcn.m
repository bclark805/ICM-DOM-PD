function [alpha,a1u,a2u,aw0,awx,awy]=shape_coef_gcn(vx,vy,nv,xc,yc,art,ntve,nbve,nbe,isbce,isonb,spherical)
%   function [alpha,a1u,a2u,aw0,awx,awy]=shape_coef_gcn(vx,vy,nv,xc,yc,art,ntve,nbve,nbe,isbce,isonb,spherical)
%----------------------------------------------------------------------%
%  calculate the coefficient for a linear                              %
%  function on the x-y plane, i.e.:                                    %
%                     r(x,y;phai)=phai_c+cofa1*x+cofa2*y               %
%     isbce(i)=0    cells on the boundary                              %
%     isbce(i)=1    cells in the interior                              %
%----------------------------------------------------------------------%
%
%inputs:
%        vx(1:m)    --- x coordinate of all nodes
%        vy(1:m)    --- y coordinate of all nodes
%        nv(1:n,1:3)--- element to node connectivity
%        xc(1:n)    --- x coord of centroid of elements
%        yc(1:n)    --- y coord of centroid of elements
%        art(1:n)   --- area (m^2) of all elmenents
%        ntve(1:m,1:mx_nbr_elem+1)  --- total number of surrounding elements of a node
%        nbve(1:m,1:mx_nbr_elem+1)  --- element number of surrounding elements of a given node
%                                       counter clockwisely ordered
%        nbe(1:n,1:3)--- elemement numbers of a given element's 3 neighbors
%                          nbe(i,1) stores the interfacing elemnt that bounds edge n2-n3
%                          nbe(i,2)                                                n1-n3
%                          nbe(i,3)                                                n1-n2
%                          where n1,n2,n3 are the node numbers of element i
%                          value is zero if one element edge does not have corresponding neighbor
%        isbce --- type of boundary condition for a given element
%                  0 -- this element is in interior
%                  1 -- this element is on solid boundary
%                  2 -- this elemenent is on open boundary
%                  3 -- this element has two solid boundaries
%        isonb --- type of boundary on a given node
%                  0 -- this node is in interior
%                  1 -- this node is on solid boundary
%                  2 -- this node is on open boundary
%        spherical ---flag for spherical coordinate, 1-spherical, 0-plane
%outputs: 
%         alpha(1:n)    --- angle (orientation) of boundary element's boundary edge
%                           range [-pi, pi] counterclockwise from south 
%
%         a1u(1:n,1:4)  --- coefficients of x gradient based on element and its 3 neighbors
%                           df/dx ~= a1u(i,1)*fi+a1u(i,2)*f1+a1u(i,3)*f2+a1u(i,4)*f3
%                                  where fi is the value at center of element i
%                                        f1~f3 are the values of f at centroids of elemenent i's three neighboring
%                                        elements
%         a2u(1:n,1:4)  --- coefficients of y gradient based on element and its 3 neighbors
%                           df/dy ~= a2u(i,1)*fi+a2u(i,2)*f1+a2u(i,3)*f2+a2u(i,4)*f3
%                                  where fi is the value at center of element i
%                                        f1~f3 are the values of f at centroids of elemenent i's three neighboring
%                                        elements
%         aw0(1:n,1:3)  --- coefficients of estimating ceontroid value  based on node values of an element
%                           fi= aw0(i,1)*f1+aw0(i,2)*f2+aw0(i,3)*f3
%                               where fi is value at centroid, f1,f2,f3 are values on nodes of a triangular element i
%         awx(1:n,1:3)  --- coefficients of estimating x gradient of function f within an element based 
%                           on nodal values of the element
%                           df/dx = awx(i,1)*f1+awx(i,2)*f2+awx(i,3)*f3
%                           where f1,f2,f3 are values on the nodes of an element i
%         awy(1:n,1:3)  --- coefficients of estimating y gradient of function f within an element based
%                           on nodal values of the element
%                           df/dy = awx(i,1)*f1+awx(i,2)*f2+awx(i,3)*f3
%                           where f1,f2,f3 are values on the nodes of an element i
%
%
%    Code and documentation
%         Wen Long, PNNL/BSRC @Seattle, 04/03/2013


%   use all_vars
%#  if  (spherical)
%   use modherical
%#  endif
%   implicit none


    x1=0;
    x2=0;
    x3=0; 

    y1=0;
    y2=0;
    y3=0;

    delt=0;

    ai1=0;
    ai2=0;
    ai3=0;

    bi1=0;
    bi2=0;
    bi3=0;

    ci1=0;
    ci2=0;
    ci3=0;

    deltx=0;
    delty=0;

    temp1=0;

    ang1=0;
    ang2=0;

    b1=0;
    b2=0;

    angle=0;

    i=0;
    ii=0;
    j=0;
    jj=0;
    j1=0;
    j2=0;

if (spherical)
    deg2rad=3.14159265/180.0;
    rearth=6371000.0;
    tpi=deg2rad*rearth;   %length(m) per degree on earth-sphere's great circle

    xxc=0;
    yyc=0;
    xxc1=0;
    yyc1=0;
    xxc2=0;
    yyc2=0;
    xxc3=0;
    yyc3=0;
    side=0;
     ty1=0;
     ty2=0;
   x1_dp=0;
   y1_dp=0;
   x2_dp=0;
   y2_dp=0;

   xtmp1=0;
   xtmp2=0;
   xtmp3=0;

   xtmp11=0;
   xtmp21=0;
   xtmp31=0;

end

%
%---------------interior cells-----------------------------------------%
%

   n=size(nv,1);   %number of elements
   m=size(vx,1);   %number of nodes

   for i=1:n             %loop through all elements
     if(isbce(i) == 0)
       y1 = yc(nbe(i,1))-yc(i);  %y coord difference of first neighbor element's centroid
       y2 = yc(nbe(i,2))-yc(i);  %                      2nd 
       y3 = yc(nbe(i,3))-yc(i);  %                      3rd

if  (spherical)
       x1_dp = xc(i);
       y1_dp = yc(i);
       x2_dp = xc(nbe(i,1));
       y2_dp = yc(nbe(i,1));
       side=arcx_fvcom(x1_dp,y1_dp,x2_dp,y2_dp);        %x coordinate difference of first neighbor element's centroid
       x1=side;

       x2_dp=xc(nbe(i,2));
       y2_dp=yc(nbe(i,2));
       side=arcx_fvcom(x1_dp,y1_dp,x2_dp,y2_dp);        %                                 2nd
       x2=side;            

       x2_dp=xc(nbe(i,3));
       y2_dp=yc(nbe(i,3));
       side=arcx_fvcom(x1_dp,y1_dp,x2_dp,y2_dp);        %                                 3rd
       x3=side;

       y1=tpi*y1;     %y coord difference of 1st neighbor element's centroid
       y2=tpi*y2;     %                      2nd
       y3=tpi*y3;     %                      3rd
else
       x1=xc(nbe(i,1))-xc(i);   %xcoord difference of first neighbor eliement relative to element i
       x2=xc(nbe(i,2))-xc(i);   %                     2nd  
       x3=xc(nbe(i,3))-xc(i);   %                     3rd

end

       %convert to km units before doing the squares to prevent overflow of floating numbers
       x1=x1/1000.0;
       x2=x2/1000.0;
       x3=x3/1000.0;
       y1=y1/1000.0;
       y2=y2/1000.0;
       y3=y3/1000.0;


       delt=(x1*y2-x2*y1)^2+(x1*y3-x3*y1)^2+(x2*y3-x3*y2)^2; 

       delt=delt*1000. ;  

       a1u(i,1)=(y1+y2+y3)*(x1*y1+x2*y2+x3*y3)-  ...
                (x1+x2+x3)*(y1^2+y2^2+y3^2);
       a1u(i,1)=a1u(i,1)/delt;

       a1u(i,2)=(y1^2+y2^2+y3^2)*x1-(x1*y1+x2*y2+x3*y3)*y1;
       a1u(i,2)=a1u(i,2)/delt;

       a1u(i,3)=(y1^2+y2^2+y3^2)*x2-(x1*y1+x2*y2+x3*y3)*y2;
       a1u(i,3)=a1u(i,3)/delt;

       a1u(i,4)=(y1^2+y2^2+y3^2)*x3-(x1*y1+x2*y2+x3*y3)*y3;
       a1u(i,4)=a1u(i,4)/delt;

       a2u(i,1)=(x1+x2+x3)*(x1*y1+x2*y2+x3*y3)- ...
                (y1+y2+y3)*(x1^2+x2^2+x3^2);
       a2u(i,1)=a2u(i,1)/delt;

       a2u(i,2)=(x1^2+x2^2+x3^2)*y1-(x1*y1+x2*y2+x3*y3)*x1;
       a2u(i,2)=a2u(i,2)/delt;

       a2u(i,3)=(x1^2+x2^2+x3^2)*y2-(x1*y1+x2*y2+x3*y3)*x2;
       a2u(i,3)=a2u(i,3)/delt;

       a2u(i,4)=(x1^2+x2^2+x3^2)*y3-(x1*y1+x2*y2+x3*y3)*x3;
       a2u(i,4)=a2u(i,4)/delt;
     end

if  (spherical)

     x1=vx(nv(i,1));
     x2=vx(nv(i,2));
     x3=vx(nv(i,3));

     y1=vy(nv(i,1));
     y2=vy(nv(i,2));
     y3=vy(nv(i,3));

     ai1=tpi*(y2-y3);   %spherical y coordinate difference (meter)
     ai2=tpi*(y3-y1);   %
     ai3=tpi*(y1-y2);

     side=arcx_fvcom(x2,y2,x3,y3);  %spherical x coordinate difference (meter)
     bi1=side                    ;
     side=arcx_fvcom(x3,y3,x1,y1);
     bi2=side;
     side=arcx_fvcom(x1,y1,x2,y2);
     bi3=side;

     x2_dp = xc(i);  %element i'th center corodinate (lon)
     y2_dp = yc(i);  %                               (lat)

     [xxc1,yyc1]=arcc_fvcom(x1,y1,x2_dp,y2_dp);  %mid point of node 1 and centroid
     [xxc2,yyc2]=arcc_fvcom(x2,y2,x2_dp,y2_dp);  %mid point of node 2 and centroid
     [xxc3,yyc3]=arcc_fvcom(x3,y3,x2_dp,y2_dp);  %mid point of node 3 and centroid

     xtmp1  = x1*tpi-xc(i)*tpi;     %x coordinate diffrence in meter of node 1 relative to centroid
     xtmp2  = x2*tpi-xc(i)*tpi;     %                                        2
     xtmp3  = x3*tpi-xc(i)*tpi;     %                                        3

     xtmp11 = x1-xc(i);             %longitude difference of node 1 relative to centroid
     xtmp21 = x2-xc(i);             %                             2
     xtmp31 = x3-xc(i);             %                             3
     
     %limit xccord difference within [-180*tpi 180*tpi], i.e. [-pi*R, pi*R] range, where pi=3.14159265, R is earth radius
     if(xtmp11 >  180.0)
       xtmp1 = -360.0*tpi+xtmp1 ;
     else(xtmp11 < -180.0)
       xtmp1 =  360.0*tpi+xtmp1 ;
     end  
     if(xtmp21 >  180.0)
       xtmp2 = -360.0*tpi+xtmp2 ;
     else(xtmp21 < -180.0)
       xtmp2 =  360.0*tpi+xtmp2 ;
     end  
     if(xtmp31 >  180.0)
       xtmp3 = -360.0*tpi+xtmp3 ;
     else(xtmp31 < -180.0)
       xtmp3 =  360.0*tpi+xtmp3 ;
     end
       
     ci1=xtmp2*(y3-yc(i))*tpi*cos(deg2rad*yyc2)- ...  
         xtmp3*(y2-yc(i))*tpi*cos(deg2rad*yyc3) ;    %dx2*dy3-dx3*dy2 on spherical surface

     ci2=xtmp3*(y1-yc(i))*tpi*cos(deg2rad*yyc3)- ...
         xtmp1*(y3-yc(i))*tpi*cos(deg2rad*yyc1) ;

     ci3=xtmp1*(y2-yc(i))*tpi*cos(deg2rad*yyc1)- ...
         xtmp2*(y1-yc(i))*tpi*cos(deg2rad*yyc2) ;
	     
else

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
end

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

%
%--------boundary cells------------------------------------------------%
%
   for i=1:n
     if(isbce(i) > 1)
       for j=1:4
         a1u(i,j)=0.0;
         a2u(i,j)=0.0;
       end
     elseif(isbce(i) == 1) 
       for j=1:3
         if(nbe(i,j) == 0)
             jj=j;
         end
       end

       j1=jj+1-floor((jj+1)/4)*3;  %the local node number on the boundary
       j2=jj+2-floor((jj+2)/4)*3;  %the next local node number on the boundary for cell i

       x1=vx(nv(i,j1))-xc(i); 
       x2=vx(nv(i,j2))-xc(i);

       y1=vy(nv(i,j1))-yc(i);
       y2=vy(nv(i,j2))-yc(i);

if  (spherical)
%      [xxc,yyc]=arcc_fvcom(vx(nv(i,j1)),vy(nv(i,j1)),xc(i),yc(i))
       ty1=0.5*(vy(nv(i,j1))+yc(i));
%      ty1=yyc

%      [xxc,yyc]=arcc_fvcom(vx(nv(i,j2)),vy(nv(i,j2)),xc(i),yc(i))
       ty2=0.5*(vy(nv(i,j2))+yc(i));
%      ty2=yyc

       xtmp1  = vx(nv(i,j1))*tpi-xc(i)*tpi;
       xtmp2  = vx(nv(i,j2))*tpi-xc(i)*tpi;
       xtmp11 = vx(nv(i,j1))-xc(i);
       xtmp21 = vx(nv(i,j2))-xc(i);
       if(xtmp11 >  180.0)
         xtmp1 = -360.0*tpi+xtmp1;
       elseif(xtmp11 < -180.0)
         xtmp1 =  360.0*tpi+xtmp1;
       end  
       if(xtmp21 >  180.0)
         xtmp2 = -360.0*tpi+xtmp2;
       elseif(xtmp21 < -180.0)
         xtmp2 =  360.0*tpi+xtmp2;
       end

       x1=xtmp1*cos(deg2rad*ty1);
       x2=xtmp2*cos(deg2rad*ty2);
       y1=tpi*y1;
       y2=tpi*y2;
end
       delt=x1*y2-x2*y1;
       b1=(y2-y1)/delt;
       b2=(x1-x2)/delt;

       deltx=vx(nv(i,j1))-vx(nv(i,j2));
       delty=vy(nv(i,j1))-vy(nv(i,j2));

if  (spherical)
       x1_dp=vx(nv(i,j1));
       y1_dp=vy(nv(i,j1));
       x2_dp=vx(nv(i,j2));
       y2_dp=vy(nv(i,j2));
       side=arcx_fvcom(x2_dp,y2_dp,x1_dp,y1_dp);

       deltx=side;
       delty=tpi*delty;
end

       alpha(i)=atan2(delty,deltx);  %alpha is basically the orirentation (angle) of the boundary line
                                     %(i.e. element i's boundary edge, when element i has isbce(i)==1)
                                     %(note that isbce=i means that cell i has only one boundary)
                                     %(hence cell i is not a corner cell)

       alpha(i)=alpha(i)-3.1415926/2.0;

       x1=xc(nbe(i,j1))-xc(i);
       x2=xc(nbe(i,j2))-xc(i);
       y1=yc(nbe(i,j1))-yc(i);
       y2=yc(nbe(i,j2))-yc(i);
if  (spherical)
%      call arcc(xc(nbe(i,j1)),yc(nbe(i,j1)),xc(i),yc(i),xxc,yyc)
       ty1=0.5*(yc(nbe(i,j1))+yc(i));
%      ty1=yyc

%      call arcc(xc(nbe(i,j2)),yc(nbe(i,j2)),xc(i),yc(i),xxc,yyc)
       ty2=0.5*(yc(nbe(i,j2))+yc(i));
%      ty2=yyc

       xtmp1  = xc(nbe(i,j1))*tpi-xc(i)*tpi;
       xtmp2  = xc(nbe(i,j2))*tpi-xc(i)*tpi;
       xtmp11 = xc(nbe(i,j1))-xc(i);
       xtmp21 = xc(nbe(i,j2))-xc(i);
       if(xtmp11 >  180.0)
         xtmp1 = -360.0*tpi+xtmp1;
       elseif(xtmp11 < -180.0)
         xtmp1 =  360.0*tpi+xtmp1;
       end	 
       if(xtmp21 >  180.0)
         xtmp2 = -360.0*tpi+xtmp2;
       else(xtmp21 < -180.0)
         xtmp2 =  360.0*tpi+xtmp2;
       end 

       x1=xtmp1*cos(deg2rad*ty1);
       x2=xtmp2*cos(deg2rad*ty2);
       y1=tpi*y1;
       y2=tpi*y2;
end
       temp1=x1*y2-x2*y1;

       if(abs(temp1)<1.e-6)
         display('shape_f of solid b. c. temp1=0');
         display(['i=' num2str(i) ',jj=' num2str(jj) ',j1=' num2str(j1) ',j2=' num2str(j2) ...
                  'x1=' num2str(x1) ',x2=' num2str(x2) ',y1=' num2str(y1) ',y2=' num2str(y2) ]);
         display(['x1*y2==',num2str(x1*y2)]);
         display(['x2*y1==',num2str(x2*y1)]);
         exit
       end

       a1u(i,1)=0.0;
       a1u(i,jj+1)=0.0;
       a1u(i,j1+1)=0.0;
       a1u(i,j2+1)=0.0;

       a2u(i,1)=0.0;
       a2u(i,jj+1)=0.0;
       a2u(i,j1+1)=0.0;
       a2u(i,j2+1)=0.0;
     end
   end


   ang1=359.9/180.0*3.1415926;
   ang2=-0.1/180.0*3.1415926;

   for i=1:m  %loop through all nodes

     if((isonb(i) ==1)&&(ntve(i)>2))   %for boundary nodes and if there are more than one elements surrounding this node

       angle=alpha(nbve(i,ntve(i)))-alpha(nbve(i,1));  %get the angle difference of last surrounding element and the first
                                                       %surrounding element

       %limit the angle to be within [-pi pi] range
       if(angle>ang1) 
         angle=100000.0;
       elseif(angle>3.1415926) 
         angle=angle-2.0*3.1415926;
       elseif(angle<-3.1415926) 
         angle=angle+2.0*3.1415926;
       elseif(angle<ang2)
         angle=100000.0;
       end

       for j=2:ntve(i)-1  %loop through 2nd to second to the last surrounding elements
         ii=nbve(i,j);    %fetch the surrounding element number 
         if(isbce(ii)~=1) %if the element is not on bc, calculate alpha for that element
                          %by dividing angle into ntve(i)-1.0 shares

            alpha(ii)=alpha(nbve(i,1))+ ...            
                      angle/(ntve(i)-1.0)*(j-1.0) ;
         end
       end

     end
   end

end %end of function

