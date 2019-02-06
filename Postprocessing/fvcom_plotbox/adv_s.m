function [f_fw_tce_edge, xflux, xflux_adv, xflux_tce_edge, xflow_tce_edge]= adv_s(spherical,northpole,...
                                                                                 mpdata,semi_implicit,...
                                                              wet_dry,multiprocessor,horzmix,ncv,ntrg ...
                           ,dt1,dz1,u,v,dltxe,dltye,ntsn,nbsn,s1,smean1,vx,vy,art2,viscofh,ncv_i,niec ...
              ,xije,yije,horcon,iobcn,i_obc_n,one_d_model,wts,dz,isonb,art1,point_st_type,inflow_type ...
                                           ,numqbc,sdis,ibfw,bfwdis3,dti,dt,dtfa,node_northarea,kb,nv ...
                                                                     ,TCE_edge_use_uniq,TCE_use_uniq)
%
% function [xflux, xflup_adv, xflux_tce_edge, xflow_tce_edge]=adv_s(spherical,northpole,mpdata,semi_implicit ...
%                                                                   ,wet_dry,multiprocessor,horzmix,ncv,ntrg ...
%                                 ,dt1,dz1,u,v,dltxe,dltye,ntsn,nbsn,s1,smean1,vx,vy,art2,viscofh,ncv_i,niec ...
%                    ,xije,yije,horcon,iobcn,i_obc_n,one_d_model,wts,dz,isonb,art1,point_st_type,inflow_type ...
%                                                 ,numqbc,sdis,ibfw,bfwdis3,dti,dt,dtfa,node_northarea,kb,nv ...
%                                                                            ,TCE_edge_use_uniq,TCE_use_uniq)
%
%   Calculates advection and horizontal diffusion terms for salinity
%
%   inputs:
%     spherical
%     northpole
%     mpdata
%     semi_implicit
%     wet_dry
%     node_northarea
%     multiprocessor
%     one_d_model
%     horzmix
%     ncv
%     ncv_i
%     ntrig
%     dt1(1:nt)
%     u(1:nt,kb)
%     v(1:nt,kb)
%     dltxe(1:nct)
%     dltye(1:nct)
%     ntsn(1:mt)
%     nbsn(1:m,mx_nbr_elem+3)
%     s1(1:mt,kb)
%     smean1(1:mt,kb)
%     vy(1:mt)
%     vx(1:mt)
%     art2(1:mt)
%     viscofh(1:mt,kb)
%     niec(1:nct,2)
%     xije(1:nct,2)
%     yije(1:nct,2)
%     horcon               --- horizontal diffusion coefficient
%     iobcn                --- number of obc nodes
%     i_obc_n(1:iobcn)     --- list of obc nodesa
%     wts(1:mt,kb)
%     dz(1:mt,kb)             delta-sigma value
%%     dzz(1:mt,kb)           delta of intra level sigma
%     dz1(1:nt,1:kb)          delta-sigma value last level is given as (:,kb)=0
%%     dzz1(1:nt,kb)          delta of intra level sigma
%     isonb(1:mt)
%     art1(1:mt)
%     point_st_type   %caculated
%     inflow_type     %node or edge
%     numqbc
%     sdis(1:numqbc)
%     ibfw               ---number of bottom flux
%     bfwdis3(1:ibfw)
%     dti                ---internal timestep, here should be time step between two history outputs
%     dt(1:mt)           ---total depth
%     dtfa(1:mt)         ---adjusted depth for mass conservation
%     kb
%     nv
%     TCE_edge_use_uniq  --- the TCE edge numbers that are going to be used to evaluate xflux_tce***
%     TCE_use_uniq       --- the TCE numbers that are going to be used to evalucate xflux
%     
%   outputs:
%       xflux(1:mt,1:kb)   (salt*m^3/sec)
%       xflux_adv(1:mt,1:kb)
%       xflux_obc                       %debug (not incuded yet)
%       xflux_tce_edge(1:ncv_i,kbm1)
%       xflow_tce_edge(1:ncv_i,kbm1)
%       f_fw_tce_edge(1:ncv_i,kbm1)     %fraction of freshwater flow
%
%   temporary variables:
%
%  	dtij
% 	s1_sf(1:mt,kb)
%  	s1_s=zeros([mt,kb]);    %% temporary salinity in modified upwind
%   	s1_sf=zeros([mt,kb]);   %% temporary salinity in modified upwind
%   	wwws=zeros([mt,kb]);
%   	wwwsf=zeros([mt,kb]);
%   	dtwwws=zeros([mt,1]);
%   	zzzflux=zeros([mt,kb]); %% temporary total flux in corrected part
%   	beta=zeros([mt,kb]);    %% temporary beta coefficient in corrected part
%   	betain=zeros([mt,kb]);  %% temporary beta coefficient in corrected part
%   	betaout=zeros([mt,kb]); %% temporary beta coefficient in corrected part
%   	s1_fresh=zeros([mt,kb]); 
%

% code and documentation 
%
% wen long, pnnl/bsrc, seattle, 04/03/2013
%

  kbm1=kb-1;
  nhe=0;  %number of halo elements
  nhn=0;  %number of halo nodes
  m=size(vx,1);  %number of nodes 
  n=size(nv,1);;  %number of elements
  mt=m+nhn;  %number of total nodes 
  nt=n+nhe;  %number of total elements

  xflux=zeros([mt,kb]); 
  xflux_adv=zeros([mt,kb]);
  pupx=zeros([mt,1]);
  pupy=zeros([mt,1]);
  pvpx=zeros([mt,1]);
  pvpy=zeros([mt,1]);

  pspx=zeros([mt,1]);
  pspy=zeros([mt,1]);
  pspxd=zeros([mt,1]);
  pspyd=zeros([mt,1]);
  viscoff=zeros([mt,1]); 

  dtij=zeros([3*nt,kbm1]);
  uvn =zeros([3*nt,kbm1]); 

   utmp    = 0;
   vtmp    = 0;

   sitai   = 0;
   ffd     = 0;
   ff1     = 0;

   x11     = 0;
   y11     = 0;
   x22     = 0;
   y22     = 0;
   x33     = 0;
   y33     = 0;

   tmp1    = 0;
   tmp2    = 0;
   xi      = 0;
   yi      = 0;
   dxa     = 0;
   dya     = 0;
   dxb     = 0;
   dyb     = 0;
   fij1    = 0;
   fij2    = 0;
   un      = 0;
   txx     = 0;
   tyy     = 0;
   fxx     = 0;
   fyy     = 0;
   viscof  = 0;
   exflux  = 0;
   temp    = 0;
   stpoint = 0;
   fact    = 0;
   fm1     = 0;
   i       = 0;
   i1      = 0;
   i2      = 0;
   ia      = 0;
   ib      = 0;
   j       = 0;
   j1      = 0;
   j2      = 0;
   k       = 0;
   jtmp    = 0;
   jj      = 0;
   ii      = 0;
   s1min   = 0;
   s1max   = 0;
   s2min   = 0;
   s2max   = 0;
if (spherical)
   ty         = 0;
   txpi       = 0;
   typi       = 0;
   xtmp1      = 0;
   xtmp       = 0;
   x1_dp      = 0;
   y1_dp      = 0;
   x2_dp      = 0;
   y2_dp      = 0;
   xii        = 0;
   yii        = 0;
   x11_tmp    = 0;
   y11_tmp    = 0;
   x33_tmp    = 0;
   y33_tmp    = 0;
if (northpole)
    vx1_tmp   = 0;
    vy1_tmp   = 0;
    vx2_tmp   = 0;
    vy2_tmp   = 0;
    txpi_tmp  = 0;
    typi_tmp  = 0;
end
end
if (semi_implicit)
   un1 = 0;
   uvn1=zeros([3*nt,kbm1]);
   dtij1=zeros([3*nt,kbm1]);
end

if (mpdata)
   smin=0;
   smax=0;
   xxxx=0;
   s1_s=zeros([mt,kb]);    %% temporary salinity in modified upwind
   s1_sf=zeros([mt,kb]);   %% temporary salinity in modified upwind
   wwws=zeros([mt,kb]);     
   wwwsf=zeros([mt,kb]);   
   dtwwws=zeros([mt,1]);  
   zzzflux=zeros([mt,kb]); %% temporary total flux in corrected part
   beta=zeros([mt,kb]);    %% temporary beta coefficient in corrected part
   betain=zeros([mt,kb]);  %% temporary beta coefficient in corrected part
   betaout=zeros([mt,kb]); %% temporary beta coefficient in corrected part
   s1_fresh=zeros([mt,kb]);    %% for source term which also bring mass volume
   itera=0;
   ntera=0;
end
%------------------------------------------------------------------------------%

   fact = 0.0;
   fm1  = 1.0;
   if(strcmp(horzmix,'closure'))
     fact = 1.0;
     fm1  = 0.0;
   end
%
%--initialize fluxes-----------------------------------------------------------%
%
   xflux          = zeros([mt,kb]);
   xflux_adv      = zeros([mt,kb]); 

   xflux_tce_edge = zeros([ncv_i,kbm1]); %
   xflow_tce_edge = zeros([ncv_i,kbm1]); %transport
   f_fw_tce_edge  = zeros([ncv_i,kbm1]);

%
%--loop over control volume sub-edges and calculate normal velocity------------%
%
%  tmpp=input('Ctrl-C to stop and debug, enter to continue');
%
%   for i=1:ncv      %loop over all  tce edges
    for itce=1:length(TCE_edge_use_uniq)   %only do a subset of edges to speed up calculation
       i=TCE_edge_use_uniq(itce); 

       i1=ntrg(i) ;  %i1 is the element number of the tce edge
       for k=1:kbm1  %loop all layers
           dtij(i,k) = dt1(i1)*dz1(i1,k);      %dt1 is total depth in previous time step
                                               %dz1 is sigma level step size in previous time step
                                               %for the k'th layer
                       %dtij is then thickness of the k'th layer
           uvn(i,k)  = v(i1,k)*dltxe(i) - u(i1,k)*dltye(i);
                       %u, v are velocity of current time step
                        
                       %uvn is basically (u,v)\dot nds 
                       %where (u,v) is velocity vector at centroid of the element i1
                       %n is the normal vector of the tce edge at the centroid 
                       %ds is the length of the tce edge
                       %i.e. volume flux through the edge (normalized by thick ness of the layer)
                       % (m/s)*m^2/m ==> m^2/sec 
                       %  
                       %
             %
             %this backs out the following orientation
             %
             %               ib o---------------*--------------->o ia              *: element edge (ia,ib)'s mid point
             %                        spoke   2 ^                                 ia: node point
             %                                  |e                                ib: node point
             %                                  |                                  x: element centroid
             %                                  |                                  e: tce edge unit vecotor
             %                          n       |                                  n: normal vector of the tce edge
             %  y ^                    <------  x                %ib is always on left hand side of e
             %    |                           1                  %
             %    |                                                    1: starting point (centroid) of the tce edge
             %    |                                                    2: ending point (mid of edge) of the tce edge
             %    ----->x
             %
             %       if u>0, v=0   ==>  uvn = -u*dltye = -u (y(2)-y(1)) < 0
             %
             %       if u=0, v>0   ==>  uvn =  v*dltxe = v (x(2)-x(1)) =0  %because x(2)==x(1) in this figure
             %
             %       ===> this confirms that uvn is flux across tce edge (1~2) in the n'th direction
             %       i.e. uvn positive is going from ia to ib, with ib always on left hand side of e vector
             %       i.e. positive uvn is flow leaving node ia.
             %                           flow going to node ib.

		if (semi_implicit)
	           dtij1(i,k) = d1(i1)*dz1(i1,k);     %d1 is current time step depth
	                                              %dtij1 is current thickness of the layer k
	           uvn1(i,k) = vf(i1,k)*dltxe(i) - uf(i1,k)*dltye(i);
        	                               %
			                       %uf and vf are velocity of previous step
                                               %
		end
       end
   end
%
%--calculate the advection and horizontal diffusion terms----------------------%
%
   for k=1:kbm1  %loop through all layers
     pspx(:)  = 0.0 ;
     pspy(:)  = 0.0 ;
     pspxd(:) = 0.0 ;
     pspyd(:) = 0.0 ;
%     for i=1:m    %loop through all nodes
      for itce=1:length(TCE_use_uniq)  %only evaluate at chosen TCE's to speed up
       i=TCE_use_uniq(itce); 
       for j=1:ntsn(i)-1  %loop through all surrounding nodes 
                          %note the collection of the nodes around node i  (if i is not on a boundary)
                          %goes back to the starting point, so first and last one collapse
                          %hence only need to loop to ntsn(i)-1

         i1=nbsn(i,j)  ;   %i1 is node number 
         i2=nbsn(i,j+1);   %i2 is the next node number
                          %i, i1, i2 makes an element that is part of the tce around node i

if (wet_dry)
         if(iswetn(i1) == 0 && iswetn(i2) == 1)
          ffd=0.5*(s1(i,k)+s1(i2,k)-smean1(i,k)-smean1(i2,k));
          ff1=0.5*(s1(i,k)+s1(i2,k));
	 elseif(iswetn(i1) == 1 && iswetn(i2) == 0)
          ffd=0.5*(s1(i1,k)+s1(i,k)-smean1(i1,k)-smean1(i,k));
          ff1=0.5*(s1(i1,k)+s1(i,k));
	 elseif(iswetn(i1) == 0 && iswetn(i2) == 0)
          ffd=0.5*(s1(i,k)+s1(i,k)-smean1(i,k)-smean1(i,k));
          ff1=0.5*(s1(i,k)+s1(i,k));
	 else
          ffd=0.5*(s1(i1,k)+s1(i2,k)-smean1(i1,k)-smean1(i2,k));
          ff1=0.5*(s1(i1,k)+s1(i2,k));
	 end
else	 
         ffd=0.5*(s1(i1,k)+s1(i2,k)-smean1(i1,k)-smean1(i2,k));  %calculate 
                                                                 %salinity at the edge center
                                                                 %with edge being the connection of nodes i1, i2
         ff1=0.5*(s1(i1,k)+s1(i2,k));
end	 
	 
if (spherical)
         xtmp  = vx(i2)*tpi-vx(i1)*tpi;
         xtmp1 = vx(i2)-vx(i1);
         if(xtmp1 >  180.0)
           xtmp = -360.0*tpi+xtmp;
         elseif(xtmp1 < -180.0)
           xtmp =  360.0*tpi+xtmp;
         end  

         txpi=xtmp*val_cos_vy(i) ;
         typi=(vy(i1)-vy(i2))*tpi;
if (northpole)
         if(node_northarea(i) == 1)
           vx1_tmp = rearth * cos(vy(i1)*deg2rad) * cos(vx(i1)*deg2rad) ...
                     * 2. /(1.+sin(vy(i1)*deg2rad));
           vy1_tmp = rearth * cos(vy(i1)*deg2rad) * sin(vx(i1)*deg2rad) ...
                     * 2. /(1.+sin(vy(i1)*deg2rad));

           vx2_tmp = rearth * cos(vy(i2)*deg2rad) * cos(vx(i2)*deg2rad) ...
                     * 2. /(1.+sin(vy(i2)*deg2rad));
           vy2_tmp = rearth * cos(vy(i2)*deg2rad) * sin(vx(i2)*deg2rad) ...
                     * 2. /(1.+sin(vy(i2)*deg2rad));

           txpi = (vx2_tmp-vx1_tmp)/(2. /(1.+sin(vy(i)*deg2rad)));
           typi = (vy1_tmp-vy2_tmp)/(2. /(1.+sin(vy(i)*deg2rad)));
  	   if(i ~= node_northpole)
	     txpi_tmp = typi*cos(vx(i)*deg2rad)-txpi*sin(vx(i)*deg2rad);
	     typi_tmp = txpi*cos(vx(i)*deg2rad)+typi*sin(vx(i)*deg2rad);
	     typi_tmp = -typi_tmp;
	    
	     txpi = txpi_tmp;
	     typi = typi_tmp;
	   end  
	 end 
end	 
         pspx(i)=pspx(i)+ff1*typi;
         pspy(i)=pspy(i)+ff1*txpi;
         pspxd(i)=pspxd(i)+ffd*typi;
         pspyd(i)=pspyd(i)+ffd*txpi;
else
         pspx(i)=pspx(i)+ff1*(vy(i1)-vy(i2));    %line integration to get ds/dx 
                                                %wth ff1 being salinity and vy(i1),vy(i2) being 
                                                %the y coordiante of i1 and i2

         pspy(i)=pspy(i)+ff1*(vx(i2)-vx(i1)) ;   %line integration to get ds/dy
                                                

         pspxd(i)=pspxd(i)+ffd*(vy(i1)-vy(i2));  %line integration to get dsmean/dx
         pspyd(i)=pspyd(i)+ffd*(vx(i2)-vx(i1));  %line integration to get dsmean/dy
end
       end %do
       pspx(i)=pspx(i)/art2(i)                ;  %finally average by the area of all elements that surround node i
       pspy(i)=pspy(i)/art2(i);
       pspxd(i)=pspxd(i)/art2(i);
       pspyd(i)=pspyd(i)/art2(i);
     end % do
 
     if(k == kbm1)              %pick out bottom value of the salinity gradient
       for i=1:m
         pfpxb(i) = pspx(i);
         pfpyb(i) = pspy(i);
       end %do
     end

     for i=1:m
       viscoff(i)=viscofh(i,k);
     end % do

     if(k == kbm1)   %pick out the bottom value of the horizontal diffusivity
       ah_bottom(1:m) = horcon*(fact*viscoff(1:m) + fm1);
     end


%     for i=1:ncv_i        %loop through all tce edges
      for itce=1:length(TCE_edge_use_uniq)   %only do a subsect of tce edges to speed up
       i=TCE_edge_use_uniq(itce); 
       ia=niec(i,1)  ;                   %one node of the triangular element edge (spoke) 
                                         %that the tce edge is intersecting with (intersecting cord)
                                         %ia is on the right hand side of the tce edge
       ib=niec(i,2) ;                    %the other node fo the tce spoke
                                         %ib is on the left hand side of the tce edge
                                         %i.e. ia, ib are nodes in counter-clockwise order for the element that the tce edge
                                         %belongs to
       xi=0.5*(xije(i,1)+xije(i,2)) ;    %mid point of the tce edge (a tce edge starts from
                                         %the surrounding element center and ends
                                         %at the spoke center)
       yi=0.5*(yije(i,1)+yije(i,2)) ;

if (spherical)
       x1_dp=xije(i,1);
       y1_dp=yije(i,1);
       x2_dp=xije(i,2);
       y2_dp=yije(i,2);

       xi=xcg2(i);
       xtmp  = xi*tpi-vx(ia)*tpi;
       xtmp1 = xi-vx(ia);
       if(xtmp1 >  180.0)
         xtmp = -360.0*tpi+xtmp;
       elseif(xtmp1 < -180.0)
         xtmp =  360.0*tpi+xtmp;
       end

       dxa=xtmp*val_cos_vy(ia);
       dya=(yi-vy(ia))*tpi;
       xtmp  = xi*tpi-vx(ib)*tpi;
       xtmp1 = xi-vx(ib);
       if(xtmp1 >  180.0)
         xtmp = -360.0*tpi+xtmp;
       elseif(xtmp1 < -180.0)
         xtmp =  360.0*tpi+xtmp;
       end

       dxb=xtmp*val_cos_vy(ib)      ;
       dyb=(yi-vy(ib))*tpi;
else
       dxa=xi-vx(ia);  %local x coordinate of mid point of tce edge  (relative to node ia) 
       dya=yi-vy(ia);  %local y coordinate of mid point of tce edge  (relative to node ia)

       dxb=xi-vx(ib);  %local x and y coordinate of mid point of tce edge (relative to node ib)
       dyb=yi-vy(ib);
end

       fij1=s1(ia,k)+dxa*pspx(ia)+dya*pspy(ia); %salinity at mid point of tce edge
                                                %estimated using one spoke end

       fij2=s1(ib,k)+dxb*pspx(ib)+dyb*pspy(ib); %saliniyt at mid point of tce estimated 
                                                %using the other spoke end

       s1min=min(s1(nbsn(ia,1:ntsn(ia)-1),k)) ;
       s1min=min(s1min, s1(ia,k));    %minimum value all surrounding nodes of ia and at ia

       s1max=max(s1(nbsn(ia,1:ntsn(ia)-1),k));
       s1max=max(s1max, s1(ia,k)) ;   %max of s of ia and around ia

       s2min=min(s1(nbsn(ib,1:ntsn(ib)-1),k));  %min and max around ib
       s2min=min(s2min, s1(ib,k));
       s2max=max(s1(nbsn(ib,1:ntsn(ib)-1),k));
       s2max=max(s2max, s1(ib,k));

       if(fij1 < s1min)
            fij1=s1min;  %this makes sure the linear interpolation does not
                         % go beyond the min and max, around ia 
       end
       if(fij1 > s1max)
          fij1=s1max;
       end
       if(fij2 < s2min)
          fij2=s2min ; %limit to min and max around ib
       end
       if(fij2 > s2max)
          fij2=s2max;
       end
       un=uvn(i,k)  ;  %normal volume flux through the tce edge (at the ceontroid of the element that 
                       %that the tce edge resides in)

if (semi_implicit)
       un1=uvn1(i,k) ; % volume flux of previous time step through the tce edge i over layer k
end

       viscof=horcon*(fact*(viscoff(ia)+viscoff(ib))*0.5 + fm1);

       txx=0.5*(pspxd(ia)+pspxd(ib))*viscof  ;
       tyy=0.5*(pspyd(ia)+pspyd(ib))*viscof;

       fxx=-dtij(i,k)*txx*dltye(i);  %layer thickness * ah*dsmean/dx (see (a.24) of chen 2003, replacing u,v with s)
       fyy= dtij(i,k)*tyy*dltxe(i);  %                               (see (a.25) of chen 2003, replacing u,v with s)

if (~semi_implicit)
       exflux=-un*dtij(i,k)* ...      %calculate transport of salinity across tce 
                             ...        %of layer k (un) is nornal volume flux per unit thickness at the edge i
                             ...       %dtij is thick nesss of the layer
          ((1.0+sign(un))*fij2+(1.0-sign(un))*fij1)*0.5+fxx+fyy;
                                    %if sign(un) is 1, i.e. flux is going out of the tce, then use fij2 
                                    %i.e. use the salinity estimated at the mid point of the tce edge
                                    %using values around node ib (which means ib is inside the tce)

				    %un has unit m^2/sec 
				    % dtij as unit m (thickness), fij has unit psu 
				    % ==> exflux has unit salt*m^3/sec 

else
       exflux=-un*dtij(i,k)* ...
          ((1.0+sign(un))*fij2+(1.0-sign(un))*fij1)*0.5;
       exflux=(1.0-ifceta)*exflux ...
              +ifceta*(-un1*dtij1(i,k)*((1.0+sign(un1))*fij2+(1.0-sign(un1))*fij1)*0.5) ...
              +fxx+fyy;
end


       exflow=-un*dtij(i,k);            %pure flow at this tce edge (m^3/sec)
       xflow_tce_edge(i,k)=exflow;      
       xflux_tce_edge(i,k)=exflux;      %record the flux at tce edge i, layer k


       xflux(ia,k)=xflux(ia,k)+exflux;  %the flux through the tce edge i is 
                                        %added to the two nodes of the spoke
                                        %flux is positive leaving tce of ia and into the tce of ib
                                        %hence exflux is positive from ia to ib
       xflux(ib,k)=xflux(ib,k)-exflux ; %exflux is neagitvely added to xflux(ib), means
                                        %exflux is positivive going into the tce of ib

       xflux_adv(ia,k)=xflux_adv(ia,k)+(exflux-fxx-fyy); %advective flux for tce of node ia
       xflux_adv(ib,k)=xflux_adv(ib,k)-(exflux-fxx-fyy); %advective flux for tce of node ib

             %because exflux (positive) is going from ib to ia
             %and exflux=-un*dtij,with dtij being thickness of the layer (>0), that means un is going from ia to ib
             %yet un is (u,v)\dot \vector(n) with \vector(n) the normal vector of tce edge
             %and \vector(n) is always 90 degrees (counter-clockwise) away from tce edge's along edge vecotr \vector(e)
             %with \vector(e) going from element center to the midpoint of spoke
             %
             %this backs out the following orientation
             %
             %               ib o---------------*---------------o ia               *: element edge (ia,ib)'s mid point
             %                        spoke   2 ^                                 ia: node point
             %                                  |e                                ib: node point
             %                                  |                                  x: element centroid
             %                                  |                                  e: tce edge unit vecotor
             %                          n       |                                  n: normal vector of the tce edge
             %                         <------  x                %ib is always on left hand side of e
             %                                1                  %
             %                                                          1: starting point (centroid) of the tce edge
             %                                                          2: ending point (mid of edge) of the tce edge
             %                                                   %positive uvn (un) gives flux into the tce of ia
             %                                                   %
             % i.e. positive xflux is convergence of flow for node ia
             %                   and  divergencee of flow for node ib


           %calculate fraction of fresh water flux assuming tracer is salinity and maximum salinity is 35
                                                                              %unit: m/s * m^2 = m^3/sec
                 %---Wen Long: How to calculate fresh water flux ?----
                 %   Let maximum salinity be Smax
                 %   volume of bulk water = V (Smax water ) + V( fresh water )                        - eqn 1
                 %   salinity of bulk water = mass of salt / (mass of salt + mass of fresh water)     - eqn 2
                 %
                 %   use eqn(1) and eqn (2) to find V(fresh water) with given salinity of bulk water S, and Smax
                 %
                 %   V_total = V_smax + V_fw                                              eqn 1
                 %         S = m_S / (m_S+m_fw)                                           eqn 2
                 %    m_S  = rho_bulk(S) * V_total - rho_fw * V_fw                        eqn 3
                 %    m_fw = rho_fw*V_fw                                                  eqn 4
                 %    V_Smax =  the volume if mass of S (m_S) is anti-diluted to Smax     eqn 5
                 %           =  m_S/rho_bulk(Smax)
                 %
                 %   ==> V_fw = V_tot -V_smax                                            %use eqn 1
                 %            = V_tot - m_S/rho_bulk(Smax)                               %use eqn 5
                 %            = V_tot - (rho_bulk(S)*V_tot - rho_fw*V_fw)/rho_bulk(Smax) %use eqn 3
                 %
                 %   ==> V_fw = V_tot - rho_bulk(S)*V_tot/rho_bulk(Smax) + rho_fw*V_fw/rho_bulk(Smax)
                 %   ==> V_fw (1 - rho_fw/rho_bulk(Smax)  = [ 1 - rho_bulk(S)/rho_bulk(Smax) ] * V_tot
                 %   ==> V_fw = V_tot* [ 1 - rho_bulk(S) /rho_bulk(Smax) ]/(1 - rho_fw/rho_bulk(Smax))
                 %   ==> f_fw = V_fw/V_tot
                 %            = [ 1 - rho_bulk(S) /rho_bulk(Smax) ]/(1 - rho_fw/rho_bulk(Smax))
                 %
                 %                 [ rho_bulk(Smax) - rho_bulk(S) ]
                 %            =   ----------------------------------                       eqn 6
                 %                 [ rho_bulk(Smax) - rho_fw ]
                 %
                 %
                 %   where V_total == V_tot  is total volume of the bulk water at salinity S, with S between 0 and Smax
                 %         V_Smax            is volume of bulk water with salinity S is anti-diluted to have salinity Smax
                 %         m_S               is mass of salt in bulk water volume V_tot with salinity S
                 %         m_fw              is mass of fresh water in bulk water volume v_tot with salinty S
                 %         V_fw              is volume of fresh water needed to mix with V_Smax water of
                 %                              salinity Smax to produce
                 %                              V_total volume of bulk water with salinity S
                 %         rho_bulk(S_prime) is density of bulk water within given salinity value S_prime
                 %         rho_bulk(Smax)    is density of seawater with salinity equal to Smax
                 %         rho_bulk(S)       is density of seawater with salinity equal to S
                 %         rho_fw = rho_bulk(S_prime=0) is density of fresh water (salinity eqtual to zero)
                 %         f_fw              is fraction of fresh water needed to generate seawater of salinity S by mixing
                 %                           fresh water with seawater of salinity Smax
                 %
                 %  For every V_total volume of water, there will be f_fw fraction of it being freshwater
                 %
                 % ==> fresh water transport = f_fw * bulk water volume transport
                 %                           = f_fw *sum( v_p * dA)         %  unit: (m^3/s)    (eqn6-1)
                 %
                 %      where v_p is velocity perpendicular to transect (m/s)
                 %            dA is cross section area elment size
                 %            sum() means to integrate over the transect area
                 %

                      %calculate fresh water flux (m^3/sec)

                          %get rho_bulk(S)    using equation of state , where S=Tp4
                          %get rho_bulk(Smax) using equation of state,  where Smax is set to 34

                 % Equation of State: Reference: Cushman-Roisin B.,Introduction to Geophysical Fluid Dynamics,
                 %                               page 36, Chapter 3,
                 %                               ISBN 0-13-353301-8, Prentice-Hall Inc., 1994
                 % rho_bulk = rho_0 * [ 1- alpha * (T-T0) + beta *( S-S0) ]  (eqn 7)
                 % where alpha = 1.7E-4  (1/K)
                 %       beta  = 7.6E-4
                 %       T0    = 10 degC
                 %       S0    = 35 (ppt)
                 %       rho0  = 1028 (kg/m^3)

                  salt_tce_edge=((1.0+sign(un))*fij2+(1.0-sign(un))*fij1)*0.5 ;  %salt on tce_edge
                                        %note when un>0, use fij2, when un<-0, use fij1
		  rho_bulk_S    = 1028 * (1 - 1.7E-4 * (10 -10 ) + 7.6E-4 *(salt_tce_edge-35));
                                        %eqn 7 assuming S=Tp4  (mixed)
                  rho_bulk_Smax = 1028 * (1 - 1.7E-4 * (10 -10 ) + 7.6E-4 *(34 -35));
                                        %eqn 7 assuming S=Smax=34 ppt (salt)
                  rho_fw        = 1028 * (1 - 1.7E-4 * (10 -10 ) + 7.6E-4 *(0  -35));
                                        %eqn 7 assuming S= 0 ppt  (fresh)
                                        %get f_fw using  eqn(6) above

                  f_fw_tce_edge(i,k) =( rho_bulk_Smax - rho_bulk_S )./(rho_bulk_Smax - rho_fw) ;


     end% do

if (spherical && northpole)
if (~semi_implicit)
       adv_s_xy_fvcom(xflux,xflux_adv,pspx,pspy,pspxd,pspyd,viscoff,k,0.0);
else
       adv_s_xy_fvcom(xflux,xflux_adv,pspx,pspy,pspxd,pspyd,viscoff,k,ifceta);
end
end

   end %do %%sigma loop

%
%-accumulate fluxes at boundary nodes
%
if (multiprocessor)
%   if(par)call node_match(0,nbn,bn_mlt,bn_loc,bnc,mt,kb,myid,nprocs,xflux,xflux_adv)
end

for k=1:kbm1
    if(iobcn > 0) 
       for i=1:iobcn
         i1=i_obc_n(i);
         xflux_obc(i,k)=xflux_adv(i1,k);  %pick out the advective flux for open boundary node
       end %do
     end
end %do


%--set boundary conditions-for fresh water flux--------------------------------%
%
if (mpdata)
	%   s. hu
	%   using smolarkiewicz, p. k; a fully multidimensional positive definite advection
	%   transport algorithm with small implicit diffusion, journal of computational
	%   physics, 54, 325-362, 1984
	%-----------------------------------------------------------------

	%-----combine all the horizontal flux first-----------------------------------

	%-------fresh water part--------------

	   s1_fresh=s1;
	   if(point_st_type == 'calculated') 
	     if(inflow_type == 'node') 
	       if(numqbc > 0) 
		 for j=1:numqbc
		   jj=inodeq(j);
		   stpoint=sdis(j);
		   for k=1:kbm1
			s1_fresh(jj,k)=sdis(j);
			xflux(jj,k)=xflux(jj,k) - qdis(j)*vqdist(j,k)*stpoint ;   %/dz(jj,k)
					%xflux is discounted by the input flux giving at the node

		   end %do
		 end %do
	       end
	     elseif(inflow_type == 'edge')
	       if(numqbc > 0) 
		 for j=1:numqbc
		   j1=n_icellq(j,1);
		   j2=n_icellq(j,2);
		   stpoint=sdis(j) ;     %%ask liu should this be stpoint1(j1)/stpoint2(j2)
	  		  		 %wen long: distributing the load on to two nodes of the edge
					 %          equally
		   for k=1:kbm1
		     s1_fresh(j1,k)=sdis(j);
		     s1_fresh(j1,k)=sdis(j);
		     xflux(j1,k)=xflux(j1,k)-  ...
				 qdis(j)*rdisq(j,1)*vqdist(j,k)*stpoint ;   %/dz(j1,k)
		     xflux(j2,k)=xflux(j2,k)-  ...
				 qdis(j)*rdisq(j,2)*vqdist(j,k)*stpoint ;   %/dz(j2,k)
		   end %do
		 end %do
	       end
	     end
	   end
	%


	% the horizontal term of advection is neglected here
	   for k=1:kbm1
	     for i=1:m
	       if(isonb(i) == 2)   %node i is on open boundary
		 xflux(i,k)=0.;
	       end
	     end %do
	   end %do
	   
	% initialize variables of mpdata

	   s1_s=zeros([mt,kb]);    %% temporary salinity in modified upwind
	   s1_sf=zeros([mt,kb]);   %% temporary salinity in modified upwind
	   wwws=zeros([mt,kb]);
	   wwwsf=zeros([mt,kb]);
	   dtwwws=zeros([mt,1]);
	   zzzflux=zeros([mt,kb]); %% temporary total flux in corrected part
	   beta=zeros([mt,kb]);    %% temporary beta coefficient in corrected part
	   betain=zeros([mt,kb]);  %% temporary beta coefficient in corrected part
	   betaout=zeros([mt,kb]); %% temporary beta coefficient in corrected part
	   s1_fresh=zeros([mt,kb]); 

	%%   first loop for vertical upwind
	%%   flux including horizontal and vertical upwind
	   for k=1:kbm1
	     for i=1:m
		 dothis=1;
	if (wet_dry)
	       if(iswetn(i)*iswetnt(i) == 1) 
		 dothis=1;
	       else
		 dothis=0;
	       end
	end
	      if(dothis)
		 if(k == 1) 
		   temp = -(wts(i,k+1)-abs(wts(i,k+1)))*s1(i,k)   ...
			  -(wts(i,k+1)+abs(wts(i,k+1)))*s1(i,k+1) ...
			  +(wts(i,k)+abs(wts(i,k)))*s1(i,k)    
		 elseif(k == kbm1) 
		   temp = +(wts(i,k)-abs(wts(i,k)))*s1(i,k-1)     ...
			  +(wts(i,k)+abs(wts(i,k)))*s1(i,k)
		 else
		   temp = -(wts(i,k+1)-abs(wts(i,k+1)))*s1(i,k)   ...
			  -(wts(i,k+1)+abs(wts(i,k+1)))*s1(i,k+1) ...
			  +(wts(i,k)-abs(wts(i,k)))*s1(i,k-1)     ...
			  +(wts(i,k)+abs(wts(i,k)))*s1(i,k)
		 end
		 temp = 0.5*temp 

		 if(k == 1)
		   smax = max(s1(nbsn(i,1:ntsn(i)),k));
		   smin = min(s1(nbsn(i,1:ntsn(i)),k));
		   smax = max(smax,s1(i,k+1),s1(i,k),s1_fresh(i,k));
		   smin = min(smin,s1(i,k+1),s1(i,k),s1_fresh(i,k));
		 elseif(k == kbm1) 
		   smax = max(s1(nbsn(i,1:ntsn(i)),k));
		   smin = min(s1(nbsn(i,1:ntsn(i)),k));
		   smax = max(smax,s1(i,k-1),s1(i,k),s1_fresh(i,k));
		   smin = min(smin,s1(i,k-1),s1(i,k),s1_fresh(i,k));
		 else
		   smax = max(s1(nbsn(i,1:ntsn(i)),k));
		   smin = min(s1(nbsn(i,1:ntsn(i)),k));
		   smax = max(smax,s1(i,k+1),s1(i,k-1),s1(i,k),s1_fresh(i,k));
		   smin = min(smin,s1(i,k+1),s1(i,k-1),s1(i,k),s1_fresh(i,k));
		 end


		 zzzflux(i,k) = temp*(dti/dt(i))/dz(i,k) + xflux(i,k)/art1(i)*(dti/dt(i))/dz(i,k) ;
		 xxxx = zzzflux(i,k)*dt(i)/dtfa(i)+s1(i,k)-s1(i,k)*dt(i)/dtfa(i) ;

		 beta(i,k)=0.5*(1.-sign(xxxx)) * (smax-s1(i,k))/(abs(xxxx)+1.e-10) ...
			  +0.5*(1.-sign(-xxxx)) * (s1(i,k)-smin)/(abs(xxxx)+1.e-10);

		 s1_sf(i,k)=s1(i,k)-min(1.,beta(i,k))*xxxx;
	%if (wet_dry)
	       end %dothis
	%end
	     end %do
	   end %do  %% sigma loop

	%----------------------------------------------------------------------------------------
	   ntera = 4
	   for itera=1:ntera   %% smolaricizw loop 
	     if(itera == 1)
	       wwwsf  = wts
	       s1_s   = s1_sf
	       dtwwws = dt
	     else
	       wwwsf  = wwws
	       s1_s   = s1_sf
	       dtwwws = dtfa
	     end
	     for k=2:kbm1
	       for i=1:m
		 temp=abs(wwwsf(i,k))-dti*(wwwsf(i,k))*(wwwsf(i,k))/dz(i,k)/dtwwws(i)
		 wwws(i,k)=temp*(s1_s(i,k-1)-s1_s(i,k))/(abs(s1_s(i,k-1))+abs(s1_s(i,k))+1.e-14)
	 
		 if(temp < 0.0 || s1_s(i,k) == 0.0) 
		   wwws(i,k)=0. ;
		 end
	       end %do 
	     end %do
	     for i=1:m
	       wwws(i,1)=0.
	     end %do

	     for i=1:m
	       smax = max(s1(nbsn(i,1:ntsn(i)),1));
	       smin = min(s1(nbsn(i,1:ntsn(i)),1));
	       smax = max(smax,s1(i,2),s1(i,1),s1_fresh(i,1));
	       smin = min(smin,s1(i,2),s1(i,1),s1_fresh(i,1));
	 
	       temp=0.5*((wwws(i,2)+abs(wwws(i,2)))*s1_s(i,2))*(dti/dtfa(i))/dz(i,1);
	       betain(i,1)=(smax-s1_s(i,1))/(temp+1.e-10);

	       temp=0.5*((wwws(i,1)+abs(wwws(i,1)))*s1_s(i,1)-        ...
			   (wwws(i,2)-abs(wwws(i,2)))*s1_s(i,1))*(dti/dtfa(i))/dz(i,1);

	       betaout(i,1)=(s1_s(i,1)-smin)/(temp+1.e-10);

	       wwwsf(i,1)=0.5*min(1.,betaout(i,1))*(wwws(i,1)+abs(wwws(i,1))) + ...
			    0.5*min(1.,betain(i,1))*(wwws(i,1)-abs(wwws(i,1)));
	     end %do

	     for k=2:kbm1-1
	       for i=1:m
		 smax = maxval(s1(nbsn(i,1:ntsn(i)),k));
		 smin = minval(s1(nbsn(i,1:ntsn(i)),k));
		 smax = max(smax,s1(i,k+1),s1(i,k-1),s1(i,k),s1_fresh(i,k));
		 smin = min(smin,s1(i,k+1),s1(i,k-1),s1(i,k),s1_fresh(i,k));
	 
		 temp=0.5*((wwws(i,k+1)+abs(wwws(i,k+1)))*s1_s(i,k+1)-  ...
			   (wwws(i,k)-abs(wwws(i,k)))*s1_s(i,k-1))*(dti/dtfa(i))/dz(i,k);
		 betain(i,k)=(smax-s1_s(i,k))/(temp+1.e-10);

		 temp=0.5*((wwws(i,k)+abs(wwws(i,k)))*s1_s(i,k)-        ...
			   (wwws(i,k+1)-abs(wwws(i,k+1)))*s1_s(i,k))*(dti/dtfa(i))/dz(i,k);
		 betaout(i,k)=(s1_s(i,k)-smin)/(temp+1.e-10);

		 wwwsf(i,k)=0.5*min(1.,betain(i,k-1),betaout(i,k))*(wwws(i,k)+abs(wwws(i,k))) + ...
			    0.5*min(1.,betain(i,k),betaout(i,k-1))*(wwws(i,k)-abs(wwws(i,k)));
	       end %do
	     end %do


	     k=kbm1;
	       for i=1:m
		 smax = maxval(s1(nbsn(i,1:ntsn(i)),k));
		 smin = minval(s1(nbsn(i,1:ntsn(i)),k));
		 smax = max(smax,s1(i,k-1),s1(i,k),s1_fresh(i,k));
		 smin = min(smin,s1(i,k-1),s1(i,k),s1_fresh(i,k));
	 
		 temp=0.5*((wwws(i,k+1)+abs(wwws(i,k+1)))*s1_s(i,k+1)-  ...
			   (wwws(i,k)-abs(wwws(i,k)))*s1_s(i,k-1))*(dti/dtfa(i))/dz(i,k);
		 betain(i,k)=(smax-s1_s(i,k))/(temp+1.e-10);

		 temp=0.5*((wwws(i,k)+abs(wwws(i,k)))*s1_s(i,k)-        ...
			   (wwws(i,k+1)-abs(wwws(i,k+1)))*s1_s(i,k))*(dti/dtfa(i))/dz(i,k);
		 betaout(i,k)=(s1_s(i,k)-smin)/(temp+1.e-10);

		 wwwsf(i,k)=0.5*min(1.,betain(i,k-1),betaout(i,k))*(wwws(i,k)+abs(wwws(i,k))) + ...
			    0.5*min(1.,betain(i,k),betaout(i,k-1))*(wwws(i,k)-abs(wwws(i,k)));
	       end %do
	     

	     wwws=wwwsf ;

	     for k=1:kbm1
	       for i=1:m
		   dothis=1;
	if (wet_dry)
		 if(iswetn(i)*iswetnt(i) == 1) 
		   dothis=1;
		 else
		   dothis=0;
		 end
	end
		 if(dothis)
		   if(k == 1) 
		     temp = -(wwws(i,k+1)-abs(wwws(i,k+1)))*s1_s(i,k)   ...
			    -(wwws(i,k+1)+abs(wwws(i,k+1)))*s1_s(i,k+1) ...
			    +(wwws(i,k)+abs(wwws(i,k)))*s1_s(i,k);
		   elseif(k == kbm1) 
		     temp = +(wwws(i,k)-abs(wwws(i,k)))*s1_s(i,k-1)     ...
			    +(wwws(i,k)+abs(wwws(i,k)))*s1_s(i,k);
		   else
		     temp = -(wwws(i,k+1)-abs(wwws(i,k+1)))*s1_s(i,k)   ...
			    -(wwws(i,k+1)+abs(wwws(i,k+1)))*s1_s(i,k+1) ...
			    +(wwws(i,k)-abs(wwws(i,k)))*s1_s(i,k-1)     ...
			    +(wwws(i,k)+abs(wwws(i,k)))*s1_s(i,k);
		   end
		   temp = 0.5*temp;
		   s1_sf(i,k)=(s1_s(i,k)-temp*(dti/dtfa(i))/dz(i,k)) ;
	%if (wet_dry)
		 end %dothis
	%end
	       end %do
	     end %do  %% sigma loop
	   end %do  %% smolarvizw loop
	%--------------------------------------------------------------------------
	% end of smolarkiewicz upwind loop
	%--------------------------------------------------------------------------

end  %end mpdata

if  (~mpdata)
	if  (one_d_model)
	    xflux = 0.0;
	end    

%--------------------------------------------------------------------
%   the central difference scheme in vertical advection
%--------------------------------------------------------------------
   for k=1:kbm1
%     for i=1:m    %loop through all nodes
      for itce=1:length(TCE_use_uniq)  %only evaluate at chosen TCE's to speed up
       i=TCE_use_uniq(itce);

       dothis=1;
if (wet_dry)
       if(iswetn(i)*iswetnt(i) == 1) 
         dothis=0;
       end
end
    if(dothis)
       if(k == 1) 
         temp=-wts(i,k+1)*(s1(i,k)*dz(i,k+1)+s1(i,k+1)*dz(i,k))/   ...
              (dz(i,k)+dz(i,k+1));
       elseif(k == kbm1) 
         temp= wts(i,k)*(s1(i,k)*dz(i,k-1)+s1(i,k-1)*dz(i,k))/(dz(i,k)+dz(i,k-1));
       else
         temp= wts(i,k)*(s1(i,k)*dz(i,k-1)+s1(i,k-1)*dz(i,k))/(dz(i,k)+dz(i,k-1))- ...
               wts(i,k+1)*(s1(i,k)*dz(i,k+1)+s1(i,k+1)*dz(i,k))/(dz(i,k)+dz(i,k+1));
       end

       if(isonb(i) == 2) 
         xflux(i,k)=temp*art1(i);
       else
         xflux(i,k)=xflux(i,k)+temp*art1(i);
       end
%if (wet_dry)
     end  %dothis
%end
     end %do
   end %do  %% sigma loop

%
%--set boundary conditions-for fresh water flux--------------------------------%
%
   if(strcmp(point_st_type, 'calculated')) 
     if(strcmp(inflow_type, 'node')) 
       if(numqbc > 0) 
         for j=1:numqbc
           jj=inodeq(j);
           stpoint=sdis(j);
           for k=1:kbm1
%             xflux(jj,k)=xflux(jj,k) - qdis(j)*vqdist(j,k)*stpoint/dz(jj,k)
             xflux(jj,k)=xflux(jj,k) - qdis(j)*vqdist(j,k)*stpoint;
           end %do
         end %do
       end
     elseif(strcmp(inflow_type, 'edge')) 
       if(numqbc > 0) 
         for j=1:numqbc
           j1=n_icellq(j,1);
           j2=n_icellq(j,2);
           stpoint=sdis(j) %%ask liu should this be stpoint1(j1)/stpoint2(j2)
           for k=1:kbm1
%             xflux(j1,k)=xflux(j1,k)-   &
%                         qdis(j)*rdisq(j,1)*vqdist(j,k)*stpoint/dz(j1,k)
%             xflux(j2,k)=xflux(j2,k)-   &
%                         qdis(j)*rdisq(j,2)*vqdist(j,k)*stpoint/dz(j2,k)
             xflux(j1,k)=xflux(j1,k)-qdis(j)*rdisq(j,1)*vqdist(j,k)*stpoint;
             xflux(j2,k)=xflux(j2,k)-qdis(j)*rdisq(j,2)*vqdist(j,k)*stpoint;
           end %do
         end %do
       end
     end
   end

end %end ~mpdata


%------------- the salinity flux from ground water ----------------------
   if(ibfw > 0)
     for i=1:m
       for j=1:ibfw
         if(i == node_bfw(j))
           xflux(i,kbm1)=xflux(i,kbm1)-bfwdis3(j)*bfwsdis3(j);    %/dz(i,kbm1)
         end
       end %do
     end %do
   end

%--update salinity-------------------------------------------------------------%
%

%  for i=1:m    %loop through all nodes
    for itce=1:length(TCE_use_uniq)  %only evaluate at chosen TCE's to speed up
       i=TCE_use_uniq(itce);

     dothis=1;
if (wet_dry)
     if(iswetn(i)*iswetnt(i) == 1 )
       dothis=1;
     else
       dothis=0;
     end
end
   if(dothis)
     for k=1:kbm1

if (~mpdata)     

       %d(Volume*salt)/dti = - xflux   
       %Volume*(salt2-salt1) = -xflux* dti
       %==> art1*height *(salt2-salt1) = -xflux*dti
       %==>  salt2 = salt1 - xflux*dti/art1/height
       %           = salt1 - xflux*dti/art1/(dt*dz)   %because height = dt*dz, where dt is total depth,
       %                                                       %and dz is delta sigma of the layer
       %xflux has unit salt*velocity*sidearea = salt*m^3/sec 


       sf1(i,k)=(s1(i,k)                              ...   %salinity at previous step
                          -xflux(i,k)/art1(i)         ...   %flux divided by area*dz and multiplied by dt/dtfa
                          *(                          ...   %flux scaled by dti/dt %here dt total depth
                              dti/(dt(i)*dz(i,k))     ...   
                           )                          ...
                )                                     ...   %sign convension: when xflux positive, sf1 < s1
               *(dt(i)/dtfa(i)) ;                            %i.e. xflux is positive leaving tce
                                                            %here dt(i) is previous step depth
else                                                        %dtfa is adjusted depth for mass conservation
                                                            %i.e. dilution due to increase of depth 
                                                            % or concentration due to decrease of depth
                                                            
       sf1(i,k)=s1_sf(i,k);
end              
     end %do

%if(wet_dry)
    else %dothis=0
     for k=1:kbm1
       sf1(i,k)=s1(i,k);
     end %do
    end  %dothis
%end
   end %do

end %end program
