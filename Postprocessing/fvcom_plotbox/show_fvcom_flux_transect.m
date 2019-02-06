gridfile='chn_island.mat';
load(gridfile);  %grid file
nv=e2n(:,2:4);
vx=xyd_n(:,2);
vy=xyd_n(:,3);
xc=xy_e(:,1);
yc=xy_e(:,2);


iobcn=5;              %temporary boundary nodes
i_obc_n=[1 2 3 4 5];

tri=nv; 
spherical=0;

 [iec,ienode,isbc,isbce,isonb,lisbce_1,lisbce_2,mx_nbr_elem,ne,nbe,nbsn,nbve,nbvt,niec,...
   nisbce_1,nisbce_2,nisbce_3,ntrg,ntsn,ntve,ncv,ncv_i,sitac,sitae,sitau,...
  deltux,deltuy,dltxc,dltxe,dltxyc,dltxye,dltyc,dltye,dltxne,dltyne,xijc,xije,yijc,yije,epor,nflag ... 
 ] =  triangle_grid_edge(nv,vx,vy,xc,yc,iobcn,i_obc_n,spherical); 

%show the areas of elements and around nodes
[art_ele art_tce art_can]=fvcom_element_area(spherical,ntve,nv,vx,vy,isonb,nbve,nbvt,xc,yc);

%
%show the boundary elements' boundary edge angle
%

[alpha,a1u,a2u,aw0,awx,awy]=shape_coef_gcn(vx,vy,nv,xc,yc,art_ele,ntve,nbve,nbe,isbce,isonb,spherical);

%
%show salt fluxes through all tce edges
%

 kb=11;
 kbm1=kb-1;
 
 nhe=0;  %number of halo elements
 nhn=0;  %number of halo nodes

 m=size(vx,1);  %number of nodes
 n=size(nv,1);;  %number of elements

 mt=m+nhn;  %number of total nodes
 nt=n+nhe;  %number of total elements

 northpole=0;
 mpdata=0;
 semi_implicit=0;
 wet_dry=0;
 multiprocessor=0;
 horzmix='closure';
 one_d_model=0;
 point_st_type='calculated';
 inflow_type='node';
 dti=1.0;

 dt1=25*ones(size(xc));
 dt=dt1;
 dtfa=dt;

 dz1=zeros([nt,kb]);
 dz1(1:nt,1:kbm1)=1.0/kbm1;
 dz=dz1;
 
 dz_n=zeros([mt,kb]);
 dz_n(1:mt,1:kbm1)=1.0/kbm1; 
 zeta_n=zeros([mt,1]); 

 zeta_n_1=zeros([mt,1]);

 dt_time=3600.0;  %time interval (in seconds)

 ibfw=0;
 bfwdis3=[];
 
 numqbc=0;
 sdis=[];
 node_northarea=[];

 s1=ones([mt,kb]);

 s1_1=ones([mt,kb]);
 smean1=ones([mt,kb]);
 viscofh=ones([mt,kb]);
 wts=zeros([mt,kb]);
 
 u=ones([nt,kb]);
 v=zeros([nt,kb]);
 horcon=0;

 art=art_ele;
 art1=art_tce;
 art2=art_can;

 [xflux,xflud_adv,xflux_tce_edge]=adv_s(spherical,northpole,mpdata,semi_implicit,wet_dry,...
                                             multiprocessor,horzmix,ncv,ntrg,dt1,dz1,u,v,...
                                           dltxe,dltye,ntsn,nbsn,s1,smean1,vx,vy,art_can,...
                                               viscofh,ncv_i,niec,xije,yije,horcon,iobcn,...
                                                i_obc_n,one_d_model,wts,dz,isonb,art_tce,...
                                              point_st_type,inflow_type,numqbc,sdis,ibfw,...
                                                 bfwdis3,dti,dt,dtfa,node_northarea,kb,nv);

 %get a transect
   xline=[-45000,43000];
   yline=[-500,2800];

   figure;
      trimesh(tri,vx,vy,'color',[0 0 0]); hold on
      %show the TCEs
      line(xije',yije','color',[0 1 1],'linestyle','--');
      %show the TCE vertecies
      plot(xije(:,2),yije(:,2),'k.');  %edge mid points
      plot(xije(:,1),yije(:,1),'ko');  %element centers
      for i=1:100
          text(vx(i),vy(i),num2str(i),'color',[1 0 0]);
      end
      axis([-46000,1000,-500,2800]);


 %loop through all TCE's
 for itce=1:m

           %calculate xp and yp which are the polygons of the TCE's
  	   [xp,yp,i_tce,i_bce]=fvcom_tce_polygon(spherical,ne,ienode,ntve,nv,vx,vy, ...
                                              isonb,nbve,nbvt,xc,yc,ncv,niec,ntrg,itce);
           if(0)
		figure;
	        trimesh(tri,vx,vy,'color',[0 0 0]); hold on
	        %show the TCEs
	        line(xije',yije','color',[0 1 1],'linestyle','--');
	        %show the TCE vertecies
	        plot(xije(:,2),yije(:,2),'k.');  %edge mid points
	        plot(xije(:,1),yije(:,1),'ko');  %element centers
	        for i=1:100
	            text(vx(i),vy(i),num2str(i),'color',[1 0 0]);
	        end
	        axis([-3000,1000,-500,2800]);
	        line(xp,yp,'color','m')
	    	np=length(xp(:));
                %label the tce polygon
	        for i=2:np-1
                     textangle=atan2(yp(i)-vy(itce),xp(i)-vx(itce));
	             text(xp(i)+(max(xp(:))-min(xp(:)))/8.0*cos(textangle) ...
                         ,yp(i)+(max(xp(:))-min(xp(:)))/8.0*sin(textangle) ...
                         ,num2str(i),'color',[1 0 0]);
	        end
                i=np ; %last point goes back to the starting point (avoid overlap of text)
                textangle=atan2(yp(i)-vy(itce),xp(i)-vx(itce));
                text(xp(i)+(max(xp(:))-min(xp(:)))*1/8.0*cos(textangle) ...
                    ,yp(i)+(max(xp(:))-min(xp(:)))*1/8.0*sin(textangle) ...
                    ,['1(' num2str(i) ')'],'color',[1 0 0]);

  		for i=1:np-1
                    itce_edge=i_tce(i);
                    if(itce_edge>0)
	                    ia=niec(itce_edge,1); %get the spoke of the TCE
        	            ib=niec(itce_edge,2);
                            %plot the spoke
                            line([vx(ia),vx(ib)],[vy(ia),vy(ib)],'color','g'); 

                            text_x=(xp(i)+xp(i+1))/2.0;
                            text_y=(yp(i)+yp(i+1))/2.0;
                            textangle=atan2(text_y-vy(itce),text_x-vx(itce));

                            %label the tce_edge
                            text( text_x+(max(xp(:))-min(xp(:)))*1/8.0*cos(textangle),...
                                  text_y+(max(xp(:))-min(xp(:)))*1/8.0*sin(textangle),...
                                  num2str(itce_edge),'color',[0 1 1]);
                    end

                    ibce_edge=i_bce(i);
                    if(ibce_edge>0)
                            text_x=(xp(i)+xp(i+1))/2.0;
                            text_y=(yp(i)+yp(i+1))/2.0;
                            textangle=atan2(text_y-mean(yp(:)),text_x-mean(xp(:)));

                            %label the boundary edge 
                            text( text_x+(max(xp(:))-min(xp(:)))*1/8.0*cos(textangle),...
                                  text_y+(max(xp(:))-min(xp(:)))*1/8.0*sin(textangle),...
                                  num2str(ibce_edge),'color',[1 0 1]);
                   end
		end
          end

          %calculate the intersections of the line with each polygon
          [xi,yi,segments]=polyxline(xp,yp,xline,yline); %get the intersections and polygons

          %calculate the cutout polygons by the segments
          [polycuts]=fvcom_segment_cutout_polygon(tri,vx,vy,xije,yije,itce,xp,yp,i_bce,i_tce,xi,yi,segments);

          %calculate the flux on each of the segments
          NSegs=length(polycuts);

          %calculate flux on each segment
          for iseg=1:NSegs
	       cflux{iseg}=0.0;
               cflux_cutcord{iseg}=0.0;
               if(~isempty(polycuts{iseg}.xpc))  %there is a valid cut

                   line(xp,yp,'color','m')
                   np=length(xp(:));
                   %label the tce polygon
                 %  for i=2:np-1
                 %     textangle=atan2(yp(i)-vy(itce),xp(i)-vx(itce));
                 %     text(xp(i)+(max(xp(:))-min(xp(:)))/5.0*cos(textangle) ...
                 %         ,yp(i)+(max(xp(:))-min(xp(:)))/5.0*sin(textangle) ...
                 %         ,num2str(i),'color',[1 0 0]);
                 %  end
                 %  i=np ; %last point goes back to the starting point (avoid overlap of text)
                 %  textangle=atan2(yp(i)-vy(itce),xp(i)-vx(itce));
                 %  text(xp(i)+(max(xp(:))-min(xp(:)))*1/5.0*cos(textangle) ...
                 %      ,yp(i)+(max(xp(:))-min(xp(:)))*1/5.0*sin(textangle) ...
                 %      ,['1(' num2str(i) ')'],'color',[1 0 0]);

                   for i=1:np-1
                    itce_edge=i_tce(i);
                    if(itce_edge>0)
                            ia=niec(itce_edge,1); %get the spoke of the TCE
                            ib=niec(itce_edge,2);
                            %plot the spoke
                     %       line([vx(ia),vx(ib)],[vy(ia),vy(ib)],'color','g');
                            text_x=(xp(i)+xp(i+1))/2.0;
                            text_y=(yp(i)+yp(i+1))/2.0;
                            textangle=atan2(text_y-vy(itce),text_x-vx(itce));
                            %label the tce_edge
                     %       text( text_x+(max(xp(:))-min(xp(:)))*1/8.0*cos(textangle),...
                     %             text_y+(max(xp(:))-min(xp(:)))*1/8.0*sin(textangle),...
                     %             num2str(itce_edge),'color',[0 1 1]);
                    end
                    ibce_edge=i_bce(i);
                    if(ibce_edge>0)
                            text_x=(xp(i)+xp(i+1))/2.0;
                            text_y=(yp(i)+yp(i+1))/2.0;
                            textangle=atan2(text_y-mean(yp(:)),text_x-mean(xp(:)));
                            %label the boundary edge
                      %      text( text_x+(max(xp(:))-min(xp(:)))*1/8.0*cos(textangle),...
                      %            text_y+(max(xp(:))-min(xp(:)))*1/8.0*sin(textangle),...
                      %            num2str(ibce_edge),'color',[1 0 1]);
                    end
                   end

                   %loop through the cutout polygon edges, collect the fluxes on all the edges
                   xpc=polycuts{iseg}.xpc;
                   ypc=polycuts{iseg}.ypc;

                   line(xpc,ypc,'color','g');  %show the cutout polygon

                   for iedge=1:length(polycuts{iseg}.itce_c(:))  %loop through all the tce edges
                                                                 %that happen to be part of the cutout polygon
                       itce_edge=polycuts{iseg}.itce_c(iedge);   %get the tce edge index
		       if(itce_edge>0) %if 0, the not a tce edge

                               %get the flux on the tce edge
	                       if(iedge==1)  %first edge on the cutout polygon, partial tce edge
                                        %get the proportion due to the cut
					dist_partial= (xpc(iedge)-xpc(iedge+1))^2 + ...
                                                      (ypc(iedge)-ypc(iedge+1))^2 ;
                                        dist_tce    = (xije(itce_edge,1)-xije(itce_edge,2) )^2 + ...
                                                      (yije(itce_edge,1)-yije(itce_edge,2) )^2;
					dist_frac=sqrt(dist_partial/dist_tce);

	                       elseif(iedge==length(polycuts{iseg}.itce_c(:))) %last edge on the cutout polygon 
                                                                               %partial tce edge
					%get the proportion due to the cut
                                        dist_partial= (xpc(iedge)-xpc(iedge+1))^2 + ...
                                                      (ypc(iedge)-ypc(iedge+1))^2 ;
                                        dist_tce    = (xije(itce_edge,1)-xije(itce_edge,2) )^2 + ...
                                                      (yije(itce_edge,1)-yije(itce_edge,2) )^2;
                                        dist_frac=sqrt(dist_partial/dist_tce);
	                       else  %normal full tce edges, collect the flux directly
					dist_frac=1.0;
        	               end

		               %collect the flux according the proportion 
                               %set flux as being positive leaving the tce polygon xp,yp
                               if(niec(itce_edge,1)==itce)  %xflux_tce_edge is into polygon xp,yp, flip the sign
                                     cflux{iseg}= cflux{iseg} - sum(xflux_tce_edge(itce_edge,:)).*dist_frac; 
                               else                         %xflux_tce_edge is entering polygon xp,yp
				     cflux{iseg}= cflux{iseg} + sum(xflux_tce_edge(itce_edge,:)).*dist_frac;
                               end

			       %show the flux on the edges of the cutout polygon

			       %show the flux on the tce edge
                               text_x=(xpc(iedge)+xpc(iedge+1))/2.0;
                               text_y=(ypc(iedge)+ypc(iedge+1))/2.0;
                               textangle=atan2(text_y-vy(itce),text_x-vx(itce));
                               %label flux on the cut out polgyon (in the polygon side) 
                               if(niec(itce_edge,1)==itce)
                                  text( text_x-(max(xpc(:))-min(xpc(:)))*1/8.0*cos(textangle),...
                                        text_y-(max(xpc(:))-min(xpc(:)))*1/8.0*sin(textangle),...
                                        num2str(-sum(xflux_tce_edge(itce_edge,:)).*dist_frac),...
                                        'color',[0 1 1]);
                               else
                                  text( text_x-(max(xpc(:))-min(xpc(:)))*1/8.0*cos(textangle),...
                                        text_y-(max(xpc(:))-min(xpc(:)))*1/8.0*sin(textangle),...
                                        num2str(+sum(xflux_tce_edge(itce_edge,:)).*dist_frac),...
                                        'color',[0 1 1]);
                               end

                       else %must be a boundary edge

			    ibce_edge=polycuts{iseg}.ibce_c(iedge); 
                            if(ibce_edge==0)  %something wrong, must be either a tce edge or a boundary edge
				error(['oops edge type unidentified with the cutout polygon']);
                            else %normal boundary edges

                                 %collect boundary flux if an open boundary, otherwise, flux through this edge is zero

	   		         %check if the boundary is an open boundary

                                 %find the element that ibce_edge is part of

                            end
                       end
                   end

                   %get the rate of change in the cutout polygon
			%current time step
				%get the volume of the full polygon xp,yp
				   art_tce(itce);

                                %get the volume of the cutout polygon xpc,ypc
                                   area_cutout=polygon_area([xpc xpc(1)],[ypc ypc(1)]);  %close polygon by going to starting
                                                                                         %point

				%get the conentration in the box
                                   c_tce=s1(itce,1:kbm1);

				%calculate the total mass
				   mass_tce=sum(c_tce.*area_cutout     ...
                                              .*(h_n(itce)+zeta_n(itce)).*dz_n(itce,1:kbm1));


			 %next time step
                                %get the volume of the full polygon xp,yp
                                   art_tce(itce);

                                %get the volume of the cutout polygon xpc,ypc
                                   area_cutout=polygon_area([xpc xpc(1)],[ypc ypc(1)]);  %close polygon by going to starting
                                                                                         %point

                                %get the conentration in the box
                                   c_tce=s1_1(itce,1:kbm1);

                                %calculate the total mass
                                   mass_tce_1=sum(c_tce.*area_cutout     ...
                                              .*(h_n(itce)+zeta_n_1(itce)).*dz_n(itce,1:kbm1));
                        %
			%calculate the rate of change = [ mass(nextstep) - mass(thisstep)] /dt
                        %

                   %get the flux on the cutting cord using mass balance

                     %
		     %d(mass)/dt =  flux_in -flux_out
                     %===> flux_out =  flux_in     - d(mass)/dt
		     %              = -cflux{iseg} - d(mass)/dt
                     %

                     cflux_cutcord{iseg}=  -cflux{iseg} -(mass_tce_1-mass_tce)/dt_time;

	                     %show the cutting segment
        	              line(segments{iseg}.x,segments{iseg}.y,'color','r');
	                     %show the flux on the cutting segment
	                      text_x=mean(segments{iseg}.x(:));
	                      text_y=mean(segments{iseg}.y(:));
	                      textangle=atan2(text_y-vy(itce),text_x-vx(itce));
                              text( text_x-(max(xpc(:))-min(xpc(:)))*1/8.0*cos(textangle),...
                                    text_y-(max(xpc(:))-min(xpc(:)))*1/8.0*sin(textangle),...
                                    num2str(cflux_cutcord{iseg}),                         ...
                                    'color',[1 0 0]);

               else %not cutting, flux for this segment is zero

                     cflux{iseg}=0.0;
                     cflux_cutcord{iseg}=0.0;

               end %end of valid cut

	  end %end of segments loop

  end  %end itce loop

