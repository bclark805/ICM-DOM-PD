function [art_ele,art_tce,art_can]=fvcom_element_area(spherical,ntve,nv,vx,vy,isonb,nbve,nbvt,xc,yc)
%
%function [art_ele,art_tce,art_can]=fvcom_element_area(spherical,ntve,nv,vx,vy,isonb,nbve,nbvt,xc,yc)
%
%    this is used to calculate the area of individual                          
%    triangle based on the three vertex coordinates and also calculate        
%    the sigma-surface area of individual control volume consisted of        
%    triangles with a common node point                                     
%  
% Outputs
%                                                                         
%  art_ele(1:nt)   : area of elements (triangle)                    
%  art_tce(1:mt)   : area of interior cv (for node value integration)
%  art_can(1:mt)   : sum area of all cells around node              
%
% Inputs
%
%  spherical : 1 for spherical coord (lon,lat), 0 for planar
%  vx(1:m)   : x coordinate of nodes
%  vy(1:m)   : y coordinate of nodes
%  nv(1:n,3) : element to node connectivity
%  isonb     : flag of boundary type of a given node
%                  0--this node is in interior
%                  1--this node is on solid boundary
%                  2--this node is on open boundary
%  ntve(1:m) :               : total number of surrounding elements of a node
%  nbve(1:m,1:mx_nbr_elem+1) : element number of surrounding elements of a given node
%  nbvt(1:m,1:mx_nbr_elem+1) : 1, 2, or 3, i.e. the local node id of an surrounding element
%                              that happens to be the same as a given node
%
%  xc(1:n)   : x coordinate of element center
%  yc(1:n)   : y coordinate of element center
%
% Code and documentation:
%
%  Wen Long, PNNL/BSRC, Seattle, 04/02/2013

   m=size(vx(:,1));  %number of nodes
   n=size(nv,1);     %number of elements

   nhn=0;      %number of halo nodes
   nhe=0;      %number of halo elements
   nt=n+nhe;   
   mt=m+nhn; 

   artmax=0;
   arttot=0;
   artmin=0;

   i=0;
   j=0;
   ii=0;
   j1=0;
   j2=0;
   max_nbre=0; 

if (spherical)   
   vx1=0;
   vx2=0;
   vx3=0;

   vy1=0;
   vy2=0;
   vy3=0;

   side1=0;
   side2=0;
   side3=0;

   area1=0;

   xx1=0;
   yy1=0;
   xx2=0;
   yy2=0;
   xxc=0;
   yyc=0;

   x1_dp=0;
   y1_dp=0;
   x2_dp=0;
   y2_dp=0;
   x3_dp=0;
   y3_dp=0;

end

   art  = zeros([nt,1]);
   art1 = zeros([mt,1]);
   art2 = zeros([mt,1]);

   max_nbre = max(ntve(:))+1 ; %maximum number of elements around a node 

   xx=zeros([2*max_nbre+1,1]);
   yy=zeros([2*max_nbre+1,1]); 
 
if (spherical)
   for i=1:nt 

     vx1=vx(nv(i,1));
     vx2=vx(nv(i,2));
     vx3=vx(nv(i,3));

     vy1=vy(nv(i,1));
     vy2=vy(nv(i,2));
     vy3=vy(nv(i,3));

     side1=arc_fvcom(vx2,vy2,vx3,vy3);
     side2=arc_fvcom(vx3,vy3,vx1,vy1);
     side3=arc_fvcom(vx1,vy1,vx2,vy2);

     area1=area_fvcom(side1,side2,side3);

     art(i)=area1;

   end
else
    for i=1:nt 
   	 art(i) = (vx(nv(i,2)) - vx(nv(i,1))) * (vy(nv(i,3)) - vy(nv(i,1))) -  ...
        	  (vx(nv(i,3)) - vx(nv(i,1))) * (vy(nv(i,2)) - vy(nv(i,1)));
    end 
    art    = abs(.5*art);
end    

   artmin = min(art(1:n));
   artmax = max(art(1:n));
   arttot = sum(art(1:n));

   for i=1:m
     if(isonb(i) == 0) 
       for j=1:ntve(i)
         ii=nbve(i,j);
         j1=nbvt(i,j);
         j2=j1+1-floor((j1+1)/4)*3;
if (spherical)
         xx1=vx(nv(ii,j1));
         yy1=vy(nv(ii,j1));
         xx2=vx(nv(ii,j2));
         yy2=vy(nv(ii,j2));

         [xxc,yyc]=arcc_fvcom(xx1,yy1,xx2,yy2); 
         xx(2*j-1)=xxc;
         yy(2*j-1)=yyc ;          
         xx(2*j)=xc(ii);
         yy(2*j)=yc(ii) ;           
else
         xx(2*j-1)=(vx(nv(ii,j1))+vx(nv(ii,j2)))*0.5-vx(i);
         yy(2*j-1)=(vy(nv(ii,j1))+vy(nv(ii,j2)))*0.5-vy(i);
         xx(2*j)=xc(ii)-vx(i);
         yy(2*j)=yc(ii)-vy(i);
end
       end
       xx(2*ntve(i)+1)=xx(1);
       yy(2*ntve(i)+1)=yy(1);

       for j=1:2*ntve(i);
if (spherical)
          x1_dp=xx(j);
          y1_dp=yy(j);
          x2_dp=xx(j+1);
          y2_dp=yy(j+1);
          x3_dp=vx(i);
          y3_dp=vy(i);
          side1=arc_fvcom(x1_dp,y1_dp,x2_dp,y2_dp);
          side2=arc_fvcom(x2_dp,y2_dp,x3_dp,y3_dp);
          side3=arc_fvcom(x3_dp,y3_dp,x1_dp,y1_dp);
          area1=area_fvcom(side1,side2,side3)     ;
          art1(i)=art1(i)+area1            ;
else
          art1(i)=art1(i)+0.5*(xx(j+1)*yy(j)-xx(j)*yy(j+1));
end
       end
       art1(i)=abs(art1(i));
     else
       for j=1:ntve(i)
         ii=nbve(i,j);
         j1=nbvt(i,j);
         j2=j1+1-floor((j1+1)/4)*3;
if(spherical)
         xx1=vx(nv(ii,j1));
         yy1=vy(nv(ii,j1));
         xx2=vx(nv(ii,j2));
         yy2=vy(nv(ii,j2));
         [xxc,yyc]=arcc_fvcom(xx1,yy1,xx2,yy2);
                  
         xx(2*j-1)=xxc;
         yy(2*j-1)=yyc;            
             
         xx(2*j)=xc(ii);
         yy(2*j)=yc(ii);
else
         xx(2*j-1)=(vx(nv(ii,j1))+vx(nv(ii,j2)))*0.5-vx(i);
         yy(2*j-1)=(vy(nv(ii,j1))+vy(nv(ii,j2)))*0.5-vy(i);
         xx(2*j)=xc(ii)-vx(i);
         yy(2*j)=yc(ii)-vy(i);
end	                
       end
        j=ntve(i)+1;
       ii=nbve(i,j-1);
       j1=nbvt(i,ntve(i));
       j2=j1+2-floor((j1+2)/4)*3;
if (spherical)
       xx1 =  vx(nv(ii,j1));
       yy1 =  vy(nv(ii,j1));
       xx2 =  vx(nv(ii,j2));
       yy2 =  vy(nv(ii,j2));
       arcc_fvcom(xx1,yy1,xx2,yy2,xxc,yyc) ;
     
       xx(2*j-1)=xxc;
       yy(2*j-1)=yyc;
          
       xx(2*j)=vx(i);
       yy(2*j)=vy(i);          

       xx(2*j+1)=xx(1);
       yy(2*j+1)=yy(1);

else
       xx(2*j-1)=(vx(nv(ii,j1))+vx(nv(ii,j2)))*0.5-vx(i);
       yy(2*j-1)=(vy(nv(ii,j1))+vy(nv(ii,j2)))*0.5-vy(i);

       xx(2*j)=vx(i)-vx(i);
       yy(2*j)=vy(i)-vy(i);

       xx(2*j+1)=xx(1);
       yy(2*j+1)=yy(1);
end
       for j=1:2*ntve(i)+2
if (spherical)
        x1_dp= xx(j);
        y1_dp= yy(j);
        x2_dp= xx(j+1);
        y2_dp= yy(j+1);
        x3_dp= vx(i);
        y3_dp= vy(i);
        side1=arc_fvcom(x1_dp, y1_dp, x2_dp, y2_dp);
        side2=arc_fvcom(x2_dp, y2_dp, x3_dp, y3_dp);
        side3=arc_fvcom(x3_dp, y3_dp, x1_dp, y1_dp);

        area1=area_fvcom(side1,side2,side3);
        
        art1(i)=art1(i)+area1             ;
else
        art1(i)=art1(i)+0.5*(xx(j+1)*yy(j)-xx(j)*yy(j+1));
end
       end
       art1(i)=abs(art1(i));
     end
   end

   for i=1:m
     art2(i) = sum(art(nbve(i,1:ntve(i))));
   end

%   art(0) = art(1) ;
%   art1(0) = art1(1); 

   if(nt > n)
       art(n+1:nt) = art(n)
   end
   if(mt > m)
      art2(m+1:mt) = art2(m) ;
   end
   if(mt > m)
      art1(m+1:mt) = art1(m) ;
   end

   clear 'xx' 'yy'; 

   art_ele=art;
   art_tce=art1;
   art_can=art2;

   clear 'art' 'art1' 'art2'  ; 

   end

