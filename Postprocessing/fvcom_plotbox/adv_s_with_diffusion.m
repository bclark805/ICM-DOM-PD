function [xflow_fw_tce_edge_advection,xflow_fw_tce_edge_diffusion,xflux, xflux_adv, xflux_tce_edge, xflow_tce_edge] = adv_s_with_diffusion(spherical,northpole, ...
                                                                                    mpdata,semi_implicit, ...
		                                                 wet_dry,multiprocessor,horzmix,ncv,ntrg, ...
			       dt1,dz1,u,v,dltxe,dltye,ntsn,nbsn,s1,smean1,vx,vy,art2,viscofh,ncv_i,niec, ...
	          xije,yije,horcon,iobcn,i_obc_n,one_d_model,wts,dz,isonb,art1,point_st_type,inflow_type, ...
	                                       numqbc,sdis,ibfw,bfwdis3,dti,dt,dtfa,node_northarea,kb,nv, ...
	                                                                  TCE_edge_use_uniq,TCE_use_uniq, ...
										    Smax, clockwise_grid)

%
%function [                                       xflow_fw_tce_edge_advection,xflow_fw_tce_edge_diffusion, ...
%             xflux, xflux_adv, xflux_tce_edge, xflow_tce_edge]= adv_s_with_diffusion(spherical,northpole, ...
%                                                                                    mpdata,semi_implicit, ...
%                                                                 wet_dry,multiprocessor,horzmix,ncv,ntrg, ...
%                               dt1,dz1,u,v,dltxe,dltye,ntsn,nbsn,s1,smean1,vx,vy,art2,viscofh,ncv_i,niec, ...
%                  xije,yije,horcon,iobcn,i_obc_n,one_d_model,wts,dz,isonb,art1,point_st_type,inflow_type, ...
%                                              ,numqbc,sdis,ibfw,bfwdis3,dti,dt,dtfa,node_northarea,kb,nv, ...
%                                                                         ,TCE_edge_use_uniq,TCE_use_uniq, ...
%										      Smax,clockwise_grid)
%
%
%Calculates advection and horizontal diffusion terms for salinity
%
%   inputs:
%     spherical			--- flag for spherical coordinate
%     northpole			--- flag for northpole 
%     mpdata			--- flag for mpdata advection scheme (not used)
%     semi_implicit		--- flag for semi_implicit time integration
%     wet_dry			--- flag for wet_dry
%     node_northarea		--- flag for node_northarea
%     multiprocessor		--- flag for parallel comptation (not used)
%     one_d_model		--- flag for one-D model (not used)
%     horzmix			--- flag for horizontal mixing scheme, string of value 'closure' or 'constant'
%     ncv			--- number of tce edges
%     ncv_i			--- number of tce edges in a local process (for parallel domain decomposition)
%     dt1(1:nt)			--- total depth at previous time step (t) at element (1)
%     u(1:nt,kb)		--- u velocity
%     v(1:nt,kb)		--- v velocity
%     dltxe(1:nct)		--- x increment of the tce edges
%     dltye(1:nct)		--- y increment of the tce edges
%     ntsn(1:mt)		--- number of surrounding nodes of a given node
%     nbsn(1:m,mx_nbr_elem+3)	--- list of surrounding nodes of a given node (clockwise or counter-clockwise depending on clockwise_grid)
%     s1(1:mt,kb)		--- salinity (or tracer) at previous time step
%     smean1(1:mt,kb)		--- mean sallinity at previous time step
%     vy(1:mt)			--- x coord of node (vertex)
%     vx(1:mt)			--- y coord of node (vertex)
%     art2(1:mt)		--- area (m^2) of polygon made by all surrounding nodes of a given node 
%     viscofh(1:mt,kbm1)	--- eddy diffusivity (horizontal), [m^2/sec], calculated in viscofh.F
%     niec(1:nct,2)		--- two node numbers of a tce edge's left or right side (left is niec(:,1), right is niec(:,2) for 
%                                       clockwised grid system, left is niec(:,2), right is niec(:,1) for counter-clockwised system
%     xije(1:nct,2)		--- x coordinates of starting and ending points of tce edges
%     yije(1:nct,2)		--- y ""  
%     horcon			--- horizontal diffusion coefficient read from input by data_run.F 
%                                   it is unitless if horzmix=='closure' or diffusivity (m^2/sec) if horzmix is not equal to 'constant'
%     iobcn			--- number of obc nodes
%     i_obc_n(1:iobcn)		--- list of obc nodesa
%     wts(1:mt,kb)		--- pseudo vertical velocity (m/sec) on a node's k level (k=1:KB)
%     dz(1:mt,kb)		--- delta-sigma value
%%    dzz(1:mt,kb)		--- delta of intra level sigma
%     dz1(1:nt,1:kb)            --- delta-sigma value last level is given as (:,kb)=0 
%%    dzz1(1:nt,kb)		--- delta of intra level sigma
%     isonb(1:mt)		--- where a node is on bounary
%     art1(1:mt)		--- area of the element
%     point_st_type		--- flag of point source ="caculated" or "specified"
%     inflow_type		--- flag of river input location "node" or "edge"
%     numqbc			--- number of river flows
%     sdis(1:numqbc)		--- concentration at river input point (psu for salt, [c] for tracers)
%     ibfw			--- number of bottom flux nodes
%     bfwdis3(1:ibfw)		--- bottom flow distribution (not used or set to zero)
%     dti			--- internal timestep, here should be time step between two history outputs
%     dt(1:mt)			--- total depth
%     dtfa(1:mt)		--- adjusted depth for mass conservation
%     kb			--- number of vertical levels
%     nv			--- element to node connectivity
%     TCE_edge_use_uniq		--- the TCE edge numbers that are going to be used to evaluate xflux_tce***
%     TCE_use_uniq		--- the TCE numbers that are going to be used to evalucate xflux
%     Smax			--- maximum value of trace concentration, eg. 34 psu for salinity
%     clockwise_grid		--- flag for clockwised grid (1) or counter-clockwised grid (0)
%     
%   outputs:
%
%       xflux(1:mt,1:kbm1)			    %([c]*m^3/sec, e.g. psu m^3/sec for salt) tracer flux (diverngenc) 
%                                                   %of a node's TCE at layer k (including vertical advection,
%                                                   % horizontal advection and diffusion, but no vertical diffusion)
%
%       xflux_adv(1:mt,1:kbm1)			    %tracer flux ([c]m^3/sec) of layer k of a node's tce (only horizontal advection)
%       xflux_obc				    %debug (not included yet)
%       xflux_tce_edge(1:ncv_i,kbm1)		    %tracer flux on the tce edge (horizontal advection + horizontal diffusion, positive
%                                                   %as leaving TCE ia of the edge
%       xflow_tce_edge(1:ncv_i,kbm1)		    %volume flux on the edge (horizontal advection only)
%       xflow_fw_tce_edge_advection(1:ncv_i,kbm1)   %advective freshwater flux (horizontal advection)
%       xflow_fw_tce_edge_diffusion(1:ncv_i,kbm1)   %diffusive freshwater flux (horizontal diffusion)
%
%   temporary variables:
%
%  	dtij
% 	s1_sf(1:mt,kb)
%  	s1_s=zeros([mt,kb]);			    %% temporary salinity in modified upwind
%   	s1_sf=zeros([mt,kb]);			    %% temporary salinity in modified upwind
%   	wwws=zeros([mt,kb]);
%   	wwwsf=zeros([mt,kb]);
%   	dtwwws=zeros([mt,1]);
%   	zzzflux=zeros([mt,kb]);			    %% temporary total flux in corrected part
%   	beta=zeros([mt,kb]);			    %% temporary beta coefficient in corrected part
%   	betain=zeros([mt,kb]);			    %% temporary beta coefficient in corrected part
%   	betaout=zeros([mt,kb]);			    %% temporary beta coefficient in corrected part
%   	s1_fresh=zeros([mt,kb]); 
%
%
% code and documentation 
%
% wen long, pnnl/bsrc, seattle, 08/28/2013
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

   %if horzmix is chosen to be 'closure' then set fm1=0, and fact=1, i.e. use model calculated eddy viscosity
   %other wise use HORCON (m^2/sec) as horizontal diffusion coefficient
   %when horzmix set to clocure, HORCON should be unitless, and 1??

   fact = 0.0;
   fm1  = 1.0;
   if(strcmp(horzmix,'closure'))
     fact = 1.0;
     fm1  = 0.0;
   end
%
%--initialize fluxes-----------------------------------------------------------%
%
   xflux          = zeros([mt,kb]);			 %tracer divergence (advection + diffusion, horizontal)
   xflux_adv      = zeros([mt,kb]);                      %tracer divergence due to (advection, horizontal)

   xflux_tce_edge = zeros([ncv_i,kbm1]);		 %trace flux (horizontal)
   xflow_tce_edge = zeros([ncv_i,kbm1]);		 %volume flux (horizontal)

%   f_fw_tce_edge  = zeros([ncv_i,kbm1]);		 %freshwater advective fraction of volume flux
   xflow_fw_tce_edge_advection =zeros([ncv_i,kbm1]);     %freshwater advective flux (horizontal)
   xflow_fw_tce_edge_diffusion =zeros([ncv_i,kbm1]);     %freshwater diffusive flux (horizontal)

%loop over control volume sub-edges and calculate normal velocity------------%
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
                     
                     %
		     % uvn = <u,v> \dot n       %where n is the normal vector 
                     %     = u * n_x + v * n_y  %n_x and n_y are the components of the normal vector
                     %   ==> n_x = -dltye
                     %       n_y =  dltxe
                     % ==> n is 90 degress counter-clockwise of vector x-->* 
                     %      

                       %u, v are velocity of current time step
                        
                       %uvn is basically (u,v)\dot n ds 
                       %where (u,v) is velocity vector at centroid of the element i1
                       %n is the normal vector of the tce edge at the centroid 
                       %ds is the length of the tce edge |x*| in the figure below, i.e. length of point 1 to point 2
                       %i.e. volume flux through the edge (normalized by thick ness of the layer)
                       % (m/s)*m^2/m ==> m^2/sec 
                       %  
                       %==> uvn is positive leaving TCE ia, i.e. n is leaving ia (for counter-clockwisely indexd grid)
             %
             %this backs out the following orientation
             %
             %               ib o---------------*--------------->o ia              *: element edge (ia,ib)'s mid point
             %                        spoke   2 ^                                 ia: node point
             %                                  |e                                ib: node point
             %                                  |                                  x: element centroid
             %                                  |                                  e: tce edge unit vecotor
             %                          n       |                                  n: normal vector of the tce edge (pointigng out of tce ia)
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
             %                            flow going to node ib.

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

     for i=1:m    %loop through all nodes
%      for itce=1:length(TCE_use_uniq)  %only evaluate at chosen TCE's to speed up
%       i=TCE_use_uniq(itce); 

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

         pspxd(i)=pspxd(i)+ffd*(vy(i1)-vy(i2));  %line integration to get d(s-smean)/dx
         pspyd(i)=pspyd(i)+ffd*(vx(i2)-vx(i1));  %line integration to get d(s-smean)/dy

         %
         %need to check sign convention of pspx, pspy, pspxd, pspyd here
         %
	 %   for counter-clockwise grid system (nodes of an element are indexed coutner-clockwise in e2n)
         %      pspx = -dS/dx
         %      pspy = -dS/dy
         %      pspxd = -d(S-Smean)/dx
         %      pspyd = -d(S-Smean)/dy
         %

end
	if(~clockwise_grid)  %Wen Long, change pspx to dS/dx, pspy to dS/dy, pspxd to d(S-Smean)/dx, pxpyd to d(S-Smean)/dy
                             %for counter-clockwise grid system
           pspx(i)=-pspx(i)
           pspy(i)=-pspy(i)
           pspxd(i)=-pspxd(i)
           pspyd(i)=-pspyd(i)
        end

       end %do

       pspx(i)=pspx(i)/art2(i);			 %finally average by the area of all elements that surround node i
       pspy(i)=pspy(i)/art2(i);                  %pspx has unit [psu/m]
       pspxd(i)=pspxd(i)/art2(i);                %pspxd also has unit [psu/m]
       pspyd(i)=pspyd(i)/art2(i);

     end % do i
 
     if(k == kbm1)              %pick out bottom value of the salinity gradient
       for i=1:m
         pfpxb(i) = pspx(i);
         pfpyb(i) = pspy(i);
       end %do
     end

     for i=1:m
       viscoff(i)=viscofh(i,k); %horizontal mixing coefficient (m^2/sec)
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
                                         %the surrounding element center of the tce node and ends
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

       txx=0.5*(pspxd(ia)+pspxd(ib))*viscof ; %calculation of d(s-smean)/dx * Am ==> negative diffusion flux of salt (positive towards increasing salt)
       tyy=0.5*(pspyd(ia)+pspyd(ib))*viscof ; %calculation of d(s-smean)/dy * Am ==> negative diffusion flux of salt (positive towards increasing salt)
		
       %sign convention of txx, tyy:
       %
       %        txx = Am*d(S-Smean)/dx
       %        tyy = Am*d(S-Smean)/dy
       %

       fxx=-dtij(i,k)*txx*dltye(i);  %layer thickness *TCE_edge_length*Ah*d(s-smean)/dx (see (a.24) of chen 2003, replacing u,v with s)
       fyy= dtij(i,k)*tyy*dltxe(i);  %                                                  (see (a.25) of chen 2003, replacing u,v with s)
          
                                     %x and y diffusion flux integrated over the side surface of a tce edge (x-*) (positive leaving ia for clockwised
                                     %grid system)
                                     %dtij is the thickness (m) of the tce edge evaluated at the center of the element (one end of the TCE edge)
                                     %dltye is delta y (y coordinate difference) between the midpoint of the element edge (spoke center)
                                     % and the center of the element. 

       %sign convention of fxx and fyy

	%counter clockwised grid:


             %
             %               ib o---------------*---------------o ia               *: element edge (ia,ib)'s mid point
             %                        spoke   2 ^                                 ia: node point
             %                                  |e                                ib: node point
             %                                  |                                  x: element centroid
             %                                  |                                  e: tce edge unit vecotor
             %                          n       |                                  n: normal vector of the tce edge (x-*)
             %                         <------  x                %ib is always on left hand side of e
             %                                1                  %
             %                                                          1: starting point (centroid) of the tce edge
             %                                                          2: ending point (mid of edge) of the tce edge
             %                                                   %positive uvn (un) gives flux into the tce of ia
             %                                                   %
             %
             %
             %
             %
             %                                                              o ib
             %                                                              |
             %                                             n ^              |s
             %                                               |              |p
             %                                               |              |o
             %                                               x------------->*k
             %                                              1            2  |e
             %                                                              |
             %                                                              |
             %                                                              o
             %                                                                ia

             %
             %fxx = -dtij*txx*dltye is positive when dltye>0, txx<0 
             %fyy =  dtij*tyy*dltxe is positive when dltxe>0, tyy>0
             %    for counter-clockwised system, txx=dS/dx ==> dS/dx<0 ==> salt increase from right to left (salt flux going right)
             %        ==> fxx is positive as salt flux toward TCE ia
             %                                   tyy=dS/dy>0 ==> dS/dy>0 ==> salt increases from bottom to top (salt flux going downward)
             %        ==> fyy is positive as salt flux toward TCE ia
             %

	%clockwised grid system:

	     
             %
             %               ia o---------------*---------------o ib               *: element edge (ia,ib)'s mid point
             %                        spoke   2 ^                                 ia: node point
             %                                  |e                                ib: node point
             %                                  |                                  x: element centroid
             %                                  |                                  e: tce edge unit vecotor
             %                          n       |                                  n: normal vector of the tce edge (x-*)
             %                         <------  x                %ib is always on right hand side of e
             %                                1                  %
             %                                                          1: starting point (centroid) of the tce edge
             %                                                          2: ending point (mid of edge) of the tce edge
             %                                                   %positive uvn (un) gives flux into the tce of ia
             %                                                   %
             %
             %
             %
             %
             %                                                              o ia
             %                                                              |
             %                                             n ^              |s
             %                                               |              |p
             %                                               |              |o
             %                                               x------------->*k
             %                                              1            2  |e
             %                                                              |
             %                                                              |
             %                                                              o
             %                                                                ib
             %
 
             %
             %fxx = -dtij*txx*dltye is positive when dltye>0, txx<0
             %fyy =  dtij*tyy*dltxe is positive when dltxe>0, tyy>0
             %    for clockwised system, txx=dS/dx<0   ==> dS/dx<0 ==> salt increases from right to left (salt flux going right)
             %        ==> fxx is positive as salt flux leaving TCE ia
             %                           tyy=dS/dy>0   ==> dS/dy>0 ==> salt increases from bottom to top (salt flux going down)
             %        ==> fyy is positive as sakt flux leaving TCE ia
             %

       %Wen Long, change to make sure fxx, fyy are positive as salt flux leaving TCE ia

       if(~clockwise_grid)  %change sign for counter-clockwised grid system
            fxx=-fxx;
            fyy=-fyy;
       end           

if (~semi_implicit)
     
       if(clockwise_grid)
           %calculate exflux positive as salt flux leaving ia
           exflux = - un*dtij(i,k)*			 ...  %calculate transport of salinity across tce 
				( (1.0+sign(un))*fij2	 ...  %of layer k (un) is nornal volume flux per unit thickness at the edge i
				 +(1.0-sign(un))*fij1	 ...  %dtij is thick nesss of the layer
				)*0.5+fxx+fyy;

             %for clockwised system, n is towards ia, un is positive towards ia
             %hence -un*dtij(i,k) is positive as volume flux leaving ia
             %if un>0, velocity vector is pointing towards ia, hence we should use upwind value of salt at mid point of x* ( ib side) (fij2)
             %if un<0, velocity vector is pointing towards ib, hence we should use upwind value of salt at mid point of x* ( ia side) (fij1)
	      
             %
             %               ia o---------------*---------------o ib               *: element edge (ia,ib)'s mid point
             %                        spoke   2 ^                                 ia: node point
             %                                  |e                                ib: node point
             %                              fij1|fij2                             x: element centroid
             %                                  |                                  e: tce edge unit vecotor
             %                          n       |                                  n: normal vector of the tce edge (x-*)
             %                         <------  x                %ib is always on right hand side of e
             %                                1                  %
             %                                                          1: starting point (centroid) of the tce edge
             %                                                          2: ending point (mid of edge) of the tce edge
             %                                                   %positive uvn (un) gives flux into the tce of ia
             %                                                   %
             %                                                              o ia
             %                                                              |
             %                                             n ^              |s
             %                                               |              |p
             %                                               |    fij1      |o
             %                                               x------------->*k
             %                                              1     fij2   2  |e
             %                                                              |
             %                                                              |
             %                                                              o
             %                                                                ib

       else  %counter-clockwise

           %calculate exflux positive as salt flux leaving ia

            exflux= un*dtij(i,k)* ...		        %calculate transport of salinity across tce
	 	   (  (1.0+sign(un))*fij1	...    %of layer k (un) is nornal volume flux per unit thickness at the edge i
		     +(1.0-sign(un))*fij2	...    %dtij is thick nesss of the layer
		   )*0.5+fxx+fyy;

	     %For counter clockwised grid system, n is towards ib, 
             %hence un is positive towards ib, hence un*dtij is positive as volume flux leaving ia 
             %if un>0, velocity is pointing towards ib, hence we should use salt on the upwind side which is mid point of x* on ia side, i.e. fij1
             %if un<0, velocity is pointing towards ia, hence we should use salt on the upwind side which is mid point of x* on ib side, i.e. fij2
             %
             %
             %               ib o---------------*---------------o ia               *: element edge (ia,ib)'s mid point
             %                        spoke   2 ^                                 ia: node point
             %                                  |e                                ib: node point
             %                             fij2 |fij1                              x: element centroid
             %                                  |                                  e: tce edge unit vecotor
             %                          n       |                                  n: normal vector of the tce edge (x-*)
             %                         <------  x                %ib is always on left hand side of e
             %                                1                  %
             %                                                          1: starting point (centroid) of the tce edge
             %                                                          2: ending point (mid of edge) of the tce edge
             %                                                   %positive uvn (un) gives flux into the tce of ia
             %                                                   %
             %                                                              o ib
             %                                                              |
             %                                             n ^              |s
             %                                               |              |p
             %                                               |    fij2      |o
             %                                               x------------->*k
             %                                              1     fij1   2  |e
             %                                                              |
             %                                                              |
             %                                                              o
             %                                                                ia

       end

else
    
       if(clockwise_grid)  %exflux positive as salt flux leaving ia

           exflux = -un*dtij(i,k)*					...
	             ( (1.0+sign(un))*fij2				...
		      +(1.0-sign(un))*fij1				...
		     )*0.5;

           exflux =  (1.0-ifceta)*exflux				...
	            +ifceta*(-un1*dtij1(i,k)*( (1.0+sign(un1))*fij2	...
					      +(1.0-sign(un1))*fij1	...
					      )*0.5			...
                            )						...
	            +fxx+fyy;
	%   tmptmptmp=input('debugdebugdebug:');


       else           %exflux positive as salt flux leaving ia (swith direction, switch fij1 and fij2)
                      %for positive un, the upwind direction of n is ia (fij1), and positive un is in n 
                      %direction and leaving ia

           exflux =  un*dtij(i,k)*					...
		     ( (1.0+sign(un))*fij1				...
		      +(1.0-sign(un))*fij2				... 
		     )*0.5;  

           exflux =  (1.0-ifceta)*exflux				...
                   + ifceta*( un1*dtij1(i,k)*( (1.0+sign(un1))*fij1	...
			  		      +(1.0-sign(un1))*fij2	...
					    )*0.5			...
			   )						...
                   + fxx+fyy; 

       end
end

       %
       %sign conventions of exflux:
       %
       %   for counter-clockwised system
       %       un*dtij is pointing to same direction of n=(n_x,n_y), and n is pointing towards ib, i.e. n is away from ia
       %	    ==> un*dtij is flow pointing towards ib, i.e. positive value of un*dtij is volume flux towards  ib (leaving ia)
       %                un*dtij(i,k)*((1.0+sign(un))*fij2+(1.0-sign(un))*fij1)*0.5 is positive as salt flux leaving ia
       %
       %  for clockwised grid system
       %       -un*dtij is pointing to opposite direction of n=(n_x,n_y) and n is pointing TOWARDS ia, i.e. n is away from ib
       %            ==> -un*dtij is flow pointing towards ib, i.e. positive value of  -un*dtij is volume flux leaving ia (diverging flow for ia)
       %                 -un*dtij(i,k)*((1.0+sign(un))*fij2+(1.0-sign(un))*fij1)*0.5 is positive as salt flux leaving ia
       %

       if(clockwise_grid)
           exflow=-un*dtij(i,k);    %diverging flow for ia for clockwise grid  
       else
           exflow= un*dtij(i,k);    %also diverging flow for ia for counter-clockwise grid 
       end

       %
       %sign conventions of exflow
       %
       %   for counter-clockwised system
       %       un*dtij is pointing to same direction of n=(n_x,n_y), and n is pointing towards ib, i.e. n is away from ia
       %            ==> un*dtij is flow pointing towards ib, i.e. positive value of un*dtij is volume flux towards  ib (leaving ia)
       %
       %  for clockwised grid system
       %       -un*dtij is pointing to opposite direction of n=(n_x,n_y) and n is pointing TOWARDS ia, i.e. n is away from ib
       %            ==> -un*dtij is flow pointing towards ib, i.e. positive value of  -un*dtij is volume flux leaving ia (diverging flow for ia)
       %

       %total water volume fluxes
	  
          xflow_tce_edge(i,k)=exflow;	    %advective flow flux  (m^3/sec) through edge i of TCE ia (positive as leaving ia)

       %tracer fluxes

	   %total tracer flux
	   xflux_tce_edge(i,k)=exflux;      %total tracer flux [c][m^3]/[sec] through edge i of TCE ia (positive as leaving ia)


      %----debug---
       %output txx,tyy
       %xflux_tce_edge(i,k)=txx;   	%see if txx is right: YES after plotting
       %xflux_tce_edge(i,k)=tyy;   	%see if tyy is correct: YES after plotting
       %xflux_tce_edge(i,k)=dtij(i,k);  %see if dtij is correct: YES after plotting
      %------------



	   xflux(ia,k)=xflux(ia,k)+exflux;  %divergence of tracer for TCE ia (positive as leaving ia) of layer k, unit ([c][m^3]/[sec])
	   xflux(ib,k)=xflux(ib,k)-exflux ; %divergence of tracer for TCE ib (positive as leaving ib) of layer k  unit ([c][m^3]/[sec])
                                            %(because exflux is positive towards ib, i.e. -exflux is positive as leaving ib)

	   %advective tracer flux 

	   xflux_adv(ia,k)=xflux_adv(ia,k)+(exflux-fxx-fyy); %advective flux for tce of node ia (exflux-fxx-fyy is diverging for ia)
           xflux_adv(ib,k)=xflux_adv(ib,k)-(exflux-fxx-fyy); %advective flux for tce of node ib (exflux-fxx-fyy is converging for ib)
        
	% freshwater flux (advection and diffusion)

	      %calculate fraction of fresh water flux assuming tracer is salinity and maximum salinity is 35
	                                                                       %unit: m/s * m^2 = m^3/sec

                 %---Wen Long: How to calculate advective fresh water flux ?----
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
                 %                 rho_bulk(Smax) - rho_fw       rho_bulk(S) - rho_fw
                 %            =   -------------------------  -  -------------------------
                 %                 rho_bulk(Smax) - rho_fw       rho_bulk(Smax) - rho_fw
                 %
                 %		   (rho_bulk(Smax) - rho_fw)/rho_fw        (rho_bulk(S) - rho_fw)/rho_fw
                 %            =   ----------------------------------   -  --------------------------------
                 %                 (rho_bulk(Smax) - rho_fw)/rho_fw        (rho_bulk(Smax) - rho_fw)/rho_fw
                 %
                 %                 Smax       S
                 %            =   -----  - --------
                 %                 Smax      Smax
                 %
                 %            =    (Smax-S)/Smax                                           eqn 7
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
                 %  For every V_total volume of water, there will be f_fw fraction of it being freshwater (advected)
                 %
                 % ==>advective  fresh water transport = f_fw * (bulk water volume transport)
                 %                                     = f_fw * sum( v_p * dA)           %  unit: (m^3/s)    (eqn6-1)
                 %
                 %      where v_p is velocity perpendicular to transect (m/s)
                 %            dA is cross section area elment size
                 %            sum() means to integrate over the transect area
                 %

                 %calculate advective fresh water flux (m^3/sec)

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
                  rho_bulk_Smax = 1028 * (1 - 1.7E-4 * (10 -10 ) + 7.6E-4 *(Smax -35));
										 %eqn 7 assuming S=Smax=Smax ppt (salt)
                  rho_fw        = 1028 * (1 - 1.7E-4 * (10 -10 ) + 7.6E-4 *(0  -35));
										 %eqn 7 assuming S= 0 ppt  (fresh)
										 %get f_fw using  eqn(6) above

                  f_fw_tce_edge =( rho_bulk_Smax - rho_bulk_S )./(rho_bulk_Smax - rho_fw);
		
		%advective freshwater flux

		  xflow_fw_tce_edge_advection(i,k) = f_fw_tce_edge * xflow_tce_edge(i,k);
                  
		%diffusive freshwater flux

                  %for fxx, fyy as horizontal diffusive salt flux (positive leaving ia) through TCE edge i
                  %there will be the following amount of fresh water flux
                  %
                  %  fxx_fw = -fxx/Smax;
                  %  fyy_fw = -fyy/Smax;
                  %
                  % The reason of division by Smax is
                  % to convert from salt flux [psu]*[m^3]/sec to water flux [m^3]/sec
                  % Sign switch is to make fxx_fw to be positive as diffusive freshwater flux leaving ia
                  % See Wen Long notes Eqn (12-236) or follow Robet Hetland 2005 JPO paper, page 1669, eqn(4)
                  %

                      fxx_fw = -fxx/Smax;
                      fyy_fw = -fyy/Smax;

                   xflow_fw_tce_edge_diffusion(i,k)   = fxx_fw+fyy_fw;   %need to return this to calling function
                                                                         %contributon of freshwater diffusive flux
                                                                         %of k'th layer on i'th TCE edge for node ia
                                                                         %and ib, positive as freshwater leaving ia
                                                                         %and towards ib direction through diffusion.
                   

		%total freshwater flux through TCE edge i, layer k (m^3/sec, positive leaving ia)

                   % WLong: this is not calculated to save memory, can be found using the following easily
		   % xflow_fw_tce_edge(i,k)= xflow_fw_tce_edge_advection(i,k)+xflow_fw_tce_edge_diffuse(i,k)  

     end% do i 

if (spherical && northpole)
if (~semi_implicit)
       adv_s_xy_fvcom(xflux,xflux_adv,pspx,pspy,pspxd,pspyd,viscoff,k,0.0);  %Wen Long, still need careful check with sign for 
                                                                             %clockwise vs counter-clockwise grid systems
else
       adv_s_xy_fvcom(xflux,xflux_adv,pspx,pspy,pspxd,pspyd,viscoff,k,ifceta);
end
end

   end %do %sigma loop k

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
         xflux_obc(i,k)=xflux_adv(i1,k);  %pick out the advective flux for open boundary node (positive as leaving node)
       end %do
     end
end %do


%--set boundary conditions-for fresh water flux--------------------------------%
%
if (mpdata) % Wen Long: still need careful check on sign of clockwise vs counter-clockwise grid systems

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
     for i=1:m    %loop through all nodes
%      for itce=1:length(TCE_use_uniq)  %only evaluate at chosen TCE's to speed up
%       i=TCE_use_uniq(itce);

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
%Wen Long: commented out the following for we do not need to update salinity
%          in this flux calculation, and we don't return sf1 back to calling 
%          function
%WLONG 
%WLONG %  for i=1:m    %loop through all nodes
%WLONG     for itce=1:length(TCE_use_uniq)  %only evaluate at chosen TCE's to speed up
%WLONG        i=TCE_use_uniq(itce);
%WLONG 
%WLONG      dothis=1;
%WLONG if (wet_dry)
%WLONG      if(iswetn(i)*iswetnt(i) == 1 )
%WLONG        dothis=1;
%WLONG      else
%WLONG        dothis=0;
%WLONG      end
%WLONG end
%WLONG    if(dothis)
%WLONG      for k=1:kbm1
%WLONG 
%WLONG if (~mpdata)     
%WLONG 
%WLONG        %d(Volume*salt)/dti = - xflux   
%WLONG        %Volume*(salt2-salt1) = -xflux* dti
%WLONG        %==> art1*height *(salt2-salt1) = -xflux*dti
%WLONG        %==> salt2 = salt1 - xflux*dti/art1/height
%WLONG        %          = salt1 - xflux*dti/art1/(dt*dz)   %because height = dt*dz, where dt is total depth,
%WLONG        %                                             %and dz is delta sigma of the layer
%WLONG        %xflux has unit salt*velocity*sidearea = salt*m^3/sec 
%WLONG        %with positive xflux, flow is divergent for tce i, and salinity will be decreasing
%WLONG        %
%WLONG 
%WLONG        sf1(i,k)=(s1(i,k)                              ...   %salinity at previous step
%WLONG                           -xflux(i,k)/art1(i)         ...   %flux divided by area*dz and multiplied by dt/dtfa
%WLONG                           *(                          ...   %flux scaled by dti/dt %here dt total depth
%WLONG                               dti/(dt(i)*dz(i,k))     ...   
%WLONG                            )                          ...
%WLONG                 )                                     ...   %sign convension: when xflux positive, sf1 < s1
%WLONG                *(dt(i)/dtfa(i)) ;                           %i.e. xflux is positive leaving tce i
%WLONG 
%WLONG else                                                       
%WLONG                                                             
%WLONG        sf1(i,k)=s1_sf(i,k);
%WLONG end              
%WLONG      end %do
%WLONG 
%WLONG %if(wet_dry)
%WLONG     else %dothis=0
%WLONG      for k=1:kbm1
%WLONG        sf1(i,k)=s1(i,k);
%WLONG      end %do
%WLONG     end  %dothis
%WLONG %end
%WLONG    end %do
%WLONG 

end %end program
