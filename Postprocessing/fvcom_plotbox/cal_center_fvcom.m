function [xcb,ycb,xcc,ycc,xcd,ycd,xce,yce,xcf,ycf,xcg2,ycg2]=cal_center_fvcom(vx,vy,nv,xc,yc,nbve,nbvt,isonb,m,n)
%
% function [xca,yca,xcb,ycb,xcc,ycc,xcd,ycd,xce,yce,xcf,ycf,xcg2,ycg2]=cal_center_fvcom(vx,vy,nv,xc,yc,nbve,nbvt,isonb,m,n)
%
% Inputs:
%        vx
%        vy
%        nv
%        xc
%        yc
%        nbve
%        nbvt
%        isonb
%        m
%        n
%
% Outputs:
%        xca
%        yca
%        xcb
%        ycb
%        xcc
%        ycc
%        xcd
%        ycd
%        xce
%        yce
%        xcf
%        ycf
%        xcg2
%        ycg2
%
% Code and documentation
%
%    Wen Long@PNNL, Seattle, 03/29/2013
%
      x11=0;
      y11=0;
      x22=0;
      y22=0;
      x33=0;
      y33=0;
     tmp1=0;
     tmp2=0;
       xi=0;
       yi=0;

     txpi=0;
     typi=0;

    xtmp1=0;
     xtmp=0;
    x1_dp=0;
    y1_dp=0;
    x2_dp=0;
    y2_dp=0;

  x11_tmp=0;
  y11_tmp=0;
  x33_tmp=0;
  y33_tmp=0;
      xii=0;
      yii=0;

        i=0;
       i1=0;
       ia=0;
       ib=0;
        j=0;
       j1=0;
       j2=0;
        k=0;
     jtmp=0;

%    nct = nt*3
%    allocate(dltxne(nct,2))       ;dltxne    = zero
%    allocate(dltyne(nct,2))       ;dltyne    = zero
%    allocate(deltux(nt,3))        ;deltux    = zero
%    allocate(deltuy(nt,3))        ;deltuy    = zero
%    allocate(sitau(nt,3))         ;sitau     = zero

     pi=3.14159265;
     deg2rad=pi/180.0;
     rearth= 6371.0E+3; 
     tpi=deg2rad*rearth;%arc length on earth's great circle per degree


     for i=1:m  %loop through all nodes


       val_cos_vy(i)=cos(deg2rad*vy(i)) ; %cosine of latitude, which is projection of R on equator plane

       j=1;
       i1=nbve(i,j);                      %node to element connectivity (j=1:ntve)
       jtmp=nbvt(i,j);                    %local node id of the sharing node of the j'th surrounding element 

       j1=jtmp+1-(jtmp+1)/4*3;            %the nonsharing node (surrounding node) of node i 
       j2=jtmp+2-(jtmp+2)/4*3;            %the nonsharing node (next surrounding node) of node i
 
                                          %spoke mid point 
       x11=0.5*(vx(i)+vx(nv(i1,j1)));     %lon of mid point of edge connecting nodes (i,j1)
       y11=0.5*(vy(i)+vy(nv(i1,j1)));     %lat of mid point of edge connecting nodes (i,j1)

       x22=xc(i1);                        %centroid of the element
       y22=yc(i1);

       x33=0.5*(vx(i)+vx(nv(i1,j2)));     %mid point lon and lat of edge connecting nodes (i,j2)
       y33=0.5*(vy(i)+vy(nv(i1,j2)));     %i.e. next spoke mid point (clockwise relative to (x11,y11))

       x1_dp=vx(i);
       y1_dp=vy(i);

       x2_dp=vx(nv(i1,j1));
       y2_dp=vy(nv(i1,j1));

       %get the center of the arc
       [x11_tmp,y11_tmp]=arcc_fvcom(x2_dp,y2_dp,x1_dp,y1_dp);  %find the mid point using great arc
       xca(i)=x11_tmp;  
       yca(i)=y11_tmp;
       x11=x11_tmp;            %x11, y33 updated using great arc midpoint instead (first spoke mid point)
       y11=y11_tmp;


       x2_dp=vx(nv(i1,j2));
       y2_dp=vy(nv(i1,j2));
       %get the center of the arc
       [x33_tmp,y33_tmp]=arcc_fvcom(x2_dp,y2_dp,x1_dp,y1_dp);
       xcb(i)=x33_tmp;
       ycb(i)=y33_tmp;
       x33=x33_tmp;           %x33,y33 updated using great arc midpoint instead (2nd spoke mid point)
       y33=y33_tmp;

       xtmp  = x33*tpi-x11*tpi; %equatorial distance (arc length on equator ) between the two spoke mid points

       xtmp1 = x33-x11;         %longitude difference

       %bracket xtmp into [-180*tpi +180*tpi] = [-pi*R pi*R] range, where R is earth radius

       if(xtmp1 >  180.0)
         xtmp = -360.0*tpi+xtmp;
       else if(xtmp1 < -180.0)
         xtmp =  360.0*tpi+xtmp;
       end if	 

       txpi=xtmp*cos(deg2rad*vy(i));
       typi=(y11-y33)*tpi;

       %on the boundary
       if(isonb(i) ~= 0) 
         xtmp  = x11*tpi-vx(i)*tpi;
         xtmp1 = x11-vx(i);
         if(xtmp1 >  180.0)
	   xtmp = -360.0*tpi+xtmp;
         else if(xtmp1 < -180.0)
	   xtmp =  360.0*tpi+xtmp;
         end if  
         txpi=xtmp*cos(deg2rad*vy(i));
         typi=(vy(i)-y11)*tpi;
       end if

       for j=2:ntve(i)-1

         i1=nbve(i,j);
         jtmp=nbvt(i,j);
         j1=jtmp+1-(jtmp+1)/4*3;
         j2=jtmp+2-(jtmp+2)/4*3;

         x11=0.5*(vx(i)+vx(nv(i1,j1)));
         y11=0.5*(vy(i)+vy(nv(i1,j1)));
         x22=xc(i1);
         y22=yc(i1);
         x33=0.5*(vx(i)+vx(nv(i1,j2)));
         y33=0.5*(vy(i)+vy(nv(i1,j2)));

         x1_dp=vx(i);
         y1_dp=vy(i);
         x2_dp=vx(nv(i1,j1));
         y2_dp=vy(nv(i1,j1));

         [x11_tmp,y11_tmp]= arcc_fvcom(x2_dp,y2_dp,x1_dp,y1_dp);
         xcc(i,j)=x11_tmp;
	 ycc(i,j)=y11_tmp;

	 x11=x11_tmp;
	 y11=y11_tmp;
         x2_dp=vx(nv(i1,j2));
         y2_dp=vy(nv(i1,j2));

         [x33_tmp,y33_tmp]=arcc_fvcom(x2_dp,y2_dp,x1_dp,y1_dp); 
         xcd(i,j)=x33_tmp;
         ycd(i,j)=y33_tmp;

         x33=x33_tmp;
	 y33=y33_tmp;

       end

       j=ntve(i);
       i1=nbve(i,j);
       jtmp=nbvt(i,j);
       j1=jtmp+1-(jtmp+1)/4*3;
       j2=jtmp+2-(jtmp+2)/4*3;

       x11=0.5*(vx(i)+vx(nv(i1,j1)));
       y11=0.5*(vy(i)+vy(nv(i1,j1)));
       x22=xc(i1);
       y22=yc(i1);
       x33=0.5*(vx(i)+vx(nv(i1,j2)));
       y33=0.5*(vy(i)+vy(nv(i1,j2)));

       x1_dp=vx(i);
       y1_dp=vy(i);
       x2_dp=vx(nv(i1,j1));
       y2_dp=vy(nv(i1,j1));

       [x11_tmp,y11_tmp]=arcc_fvcom(x2_dp,y2_dp,x1_dp,y1_dp);
       xce(i)=x11_tmp;
       yce(i)=y11_tmp;

       x11=x11_tmp;
       y11=y11_tmp;
       x2_dp=vx(nv(i1,j2));
       y2_dp=vy(nv(i1,j2));

       [x33_tmp,y33_tmp]=arcc_fvcom(x2_dp,y2_dp,x1_dp,y1_dp); 
       xcf(i)=x33_tmp;
       ycf(i)=y33_tmp;

     end

     xcg2=zeros([ncv_i,1]);
     ycg2=zeros([ncv_i,1]);

     for i=1:ncv_i
          ia=niec(i,1);
          ib=niec(i,2);
          xi=0.5*(xije(i,1)+xije(i,2));
          yi=0.5*(yije(i,1)+yije(i,2));
          x1_dp=xije(i,1);
          y1_dp=yije(i,1);
          x2_dp=xije(i,2);
          y2_dp=yije(i,2);
          [xii,yii]= arcc_fvcom(x2_dp,y2_dp,x1_dp,y1_dp);
          xcg2(i)=xii;
          ycg2(i)=yii;
     end
end
   
   
