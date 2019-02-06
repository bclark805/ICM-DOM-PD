
gridfile='chn_island.mat';

load(gridfile);  %grid file

nv=e2n(:,2:4);
vx=xyd_n(:,2);
vy=xyd_n(:,3);
xc=xy_e(:,1);
yc=xy_e(:,2);

%--debug
%vx(1)=500;
%vy(1)=-200;
%--------

iobcn=5;              %temporary boundary nodes
i_obc_n=[1 2 3 4 5];

tri=nv; 
spherical=0;


i=1;                  %tce node to show
iedge=100;            %edge to show


 [iec,ienode,isbc,isbce,isonb,lisbce_1,lisbce_2,mx_nbr_elem,ne,nbe,nbsn,nbve,nbvt,niec,...
                 nisbce_1,nisbce_2,nisbce_3,ntrg,ntsn,ntve,ncv,ncv_i,sitac,sitae,sitau,...
                     deltux,deltuy,dltxc,dltxe,dltxyc,dltxye,dltyc,dltye,dltxne,dltyne,...
                                                        xijc,xije,yijc,yije,epor,nflag ... 
 ] =  triangle_grid_edge(nv,vx,vy,xc,yc,iobcn,i_obc_n,spherical); 

if(0)
figure(1)

	%show the mesh
	trimesh(tri,vx,vy,'color',[0 0 0]); hold on

	%show the TCE
	line(xije',yije','color',[0 1 1],'linestyle','--');

	%show the TCE vertecies

	plot(xije(:,2),yije(:,2),'k.');  %edge mid points
	plot(xije(:,1),yije(:,1),'k+');  %element centers

	%show a node i and the surrounding elements

	plot(vx(i),vy(i),'rs','markersize',10,'markerfacecolor',[1 1 0])
	for ia=1:ntve(i)+1
	        iele=nbve(i,ia);
	 
	        if(iele>0)  %ia=ntve+1 is only available for interior nodes
	 		plot(xc(nbve(i,ia)),yc(nbve(i,ia)),'kp','markersize',3,'markerfacecolor',[1 0 0]);
			text(xc(nbve(i,ia)),yc(nbve(i,ia)),num2str(ia),'color',[0 1 0],'fontsize',40);
	        end
	end

	%show the local id
	if(0)
		ia=2;
		if ntve(i) >=2
		   iele=nbve(i,ia);  %element number 
		   %the sharing element's local node id
		   text(xc(nbve(i,ia)),yc(nbve(i,ia))+50,num2str(nbvt(i,ia)),'color',[0 0 1],'fontsize',40);
		   for j=1:3
		     text(vx(tri(iele,j)), ...
		          vy(tri(iele,j)), ...
		          num2str(j),      ...
		          'color',[0 0 0] );
		   end
		end
	
	end

	%show the nodes around the node i
	for ia=1:ntsn(i)
	        plot(vx(nbsn(i,ia)),vy(nbsn(i,ia)),'k^','markersize',3,'markerfacecolor',[0 1 0]);
	        text(vx(nbsn(i,ia)),vy(nbsn(i,ia)), num2str(ia),'color',[0 1 0],'fontsize',40);
	end

	%show an edge
	iele1=iec(iedge,1);  %element on one side
	iele2=iec(iedge,2);  %element on the other side

	plot(vx(ienode(iedge,1)),vy(ienode(iedge,1)),'o','markersize',10,'markerfacecolor',[1 1 0]);  %starting point of edge
	plot(vx(ienode(iedge,2)),vy(ienode(iedge,2)),'^','markersize',10,'markerfacecolor',[1 1 0]);  %ending point of edge

	%show the element center on one side of the edge
	plot(xc(iele1),yc(iele1),'s','markersize',10,'markerfacecolor',[0 0 0]);                      %
	
	%show the element on the other side of the edge
	if iele2>0
		plot(xc(iele2),yc(iele2),'p','markersize',10,'markerfacecolor',[0 0 0]);  
	end
end

%show the areas of elements and around nodes
[art_ele art_tce art_can]=fvcom_element_area(spherical,ntve,nv,vx,vy,isonb,nbve,nbvt,xc,yc);


if(0)
figure(2)

	%show the mesh
        trimesh(tri,vx,vy,'color',[0 0 0]); hold on

	%show the TCE
        line(xije',yije','color',[0 1 1],'linestyle','--');

	%show the TCE vertecies

	plot(xije(:,2),yije(:,2),'k.');  %edge mid points
	plot(xije(:,1),yije(:,1),'k+');  %element centers
 
        %show TCE area

        for i=1:10
            text(vx(i),vy(i),num2str(art_tce(i)));
        end
      
end

if(0)
figure(3)

   %show the mesh
	trimesh(tri,vx,vy,'color',[0 0 0]); hold on
 
        %show the TCE
        line(xije',yije','color',[0 1 1],'linestyle','--');

        %show the TCE vertecies

        plot(xije(:,2),yije(:,2),'k.');  %edge mid points
        plot(xije(:,1),yije(:,1),'k+');  %element centers

        %show TCE area

        for i=1:10
            text(vx(i),vy(i),num2str(art_can(i)));
        end

end

%
%show the boundary elements' boundary edge angle
%

[alpha,a1u,a2u,aw0,awx,awy]=shape_coef_gcn(vx,vy,nv,xc,yc,art_ele,ntve,nbve,nbe,isbce,isonb,spherical);

if(0)
%show the mesh
figure(4)
        trimesh(tri,vx,vy,'color',[0 0 0]); hold on

        %show the TCE
        line(xije',yije','color',[0 1 1],'linestyle','--');

        %show the TCE vertecies

        plot(xije(:,2),yije(:,2),'k.');  %edge mid points
        plot(xije(:,1),yije(:,1),'ko');  %element centers

        for i=1:size(nv,1)
              if(isbce(i)>0)  %only show boundary element values
                 text(xc(i),yc(i),num2str(alpha(i)),'color',[1 0 0]);
              end
        end

end

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
 
 ibfw=0;
 bfwdis3=[];
 
 numqbc=0;
 sdis=[];
 node_northarea=[];

 s1=ones([mt,kb]);
 smean1=ones([mt,kb]);
 viscofh=ones([mt,kb]);
 wts=zeros([mt,kb]);
 
 u=ones([nt,kb]);
 v=zeros([nt,kb]);
 horcon=0;

 art=art_ele;
 art1=art_tce;
 art2=art_can;

 [xflux,xflud_adv,xflux_tce_edge]=adv_s(spherical,northpole,mpdata,semi_implicit,wet_dry, ...
                                             multiprocessor,horzmix,ncv,ntrg,dt1,dz1,u,v, ...
                                           dltxe,dltye,ntsn,nbsn,s1,smean1,vx,vy,art_can, ...
                                               viscofh,ncv_i,niec,xije,yije,horcon,iobcn, ...
                                                i_obc_n,one_d_model,wts,dz,isonb,art_tce, ...
                                              point_st_type,inflow_type,numqbc,sdis,ibfw, ...
                                                bfwdis3,dti,dt,dtfa,node_northarea,kb,nv);
if(0)
	figure(5)
        trimesh(tri,vx,vy,'color',[0 0 0]); hold on

        %show the TCE
        line(xije',yije','color',[0 1 1],'linestyle','--');

        %show the TCE vertecies
        plot(xije(:,2),yije(:,2),'k.');  %edge mid points
        plot(xije(:,1),yije(:,1),'ko');  %element centers

        %for i=1:20
        %    text(vx(i),vy(i),num2str(xflux(i,1)),'color',[1 0 0]);
        %end

        %show uvn on all tce edges
        
        for i=1:100
           i1=ntrg(i) ;  %the element number that the i'th tce edge is on

            text(xc(i1)+dltxe(i)/2.0,yc(i1)+dltye(i)/2.0,num2str(i),'color',[0 0 0]);

%            text(xc(i1)+dltxe(i)/2.0,yc(i1)+dltye(i)/2.0,num2str(uvn(i,1)),'color',[1 0 0 ]);
                          %simply show the numbers at the middle of the tce edges
%            text(xc(i1)+dltxe(i)/2.0,yc(i1)+dltye(i)/2.0,num2str(xflux_tce_edge(i,1)),'color',[1 0 0 ]);
                          %simply show the numbers at the middle of the tce edges
        end

 
        %show divergence 

        for i=1:100
            text(vx(i),vy(i),num2str(xflux(i,1)),'color',[0 0 0]);
        end
end

%
% test cross section cutting
%

if(0)
        figure(6)
        trimesh(tri,vx,vy,'color',[0 0 0]); hold on

        %show the TCE
        line(xije',yije','color',[0 1 1],'linestyle','--');

        %show the TCE vertecies
        plot(xije(:,2),yije(:,2),'k.');  %edge mid points
        plot(xije(:,1),yije(:,1),'ko');  %element centers

        for i=1:100
            text(vx(i),vy(i),num2str(i),'color',[1 0 0]);
        end

	axis([-3000,1000,-500,2800]);

end

        %get a line
	xline=[-500,-750];
	yline=[-500,2800];

	%loop through all TCE's
        for itce=6:6

	        %calculate xp and yp which are the polygons of the TCE's
		[xp,yp,i_tce,i_bce]=fvcom_tce_polygon(spherical,ne,ienode,ntve,nv,vx,vy, ...
                                                      isonb,nbve,nbvt,xc,yc,ncv,niec,ntrg,itce);

            if(0) %---debug
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
	             text(xp(i)+(max(xp(:))-min(xp(:)))/5.0*cos(textangle) ...
                         ,yp(i)+(max(xp(:))-min(xp(:)))/5.0*sin(textangle) ...
                         ,num2str(i),'color',[1 0 0]);
	        end
                i=np ; %last point goes back to the starting point (avoid overlap of text)
                textangle=atan2(yp(i)-vy(itce),xp(i)-vx(itce));
                text(xp(i)+(max(xp(:))-min(xp(:)))*1/5.0*cos(textangle) ...
                    ,yp(i)+(max(xp(:))-min(xp(:)))*1/5.0*sin(textangle) ...
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
                            text( text_x+(max(xp(:))-min(xp(:)))*1/5.0*cos(textangle),...
                                  text_y+(max(xp(:))-min(xp(:)))*1/5.0*sin(textangle),...
                                  num2str(itce_edge),'color',[0 1 1]);
                    end

                    ibce_edge=i_bce(i);
                    if(ibce_edge>0)
                            text_x=(xp(i)+xp(i+1))/2.0;
                            text_y=(yp(i)+yp(i+1))/2.0;
                            textangle=atan2(text_y-mean(yp(:)),text_x-mean(xp(:)));

                            %label the boundary edge 
                            text( text_x+(max(xp(:))-min(xp(:)))*1/5.0*cos(textangle),...
                                  text_y+(max(xp(:))-min(xp(:)))*1/5.0*sin(textangle),...
                                  num2str(ibce_edge),'color',[1 0 1]);
                   end
		end
          end  %-----debug

                %calculate the intersections of the line with each polygon
		[xi,yi,segments]=polyxline(xp,yp,xline,yline); %get the intersections and polygons

                %calculate the cutout polygons by the segments
                [polycuts]=fvcom_segment_cutout_polygon(tri,vx,vy,xije,yije,itce,xp,yp,i_bce,i_tce,xi,yi,segments);


        end
