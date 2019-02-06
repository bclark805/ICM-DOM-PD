function [TCE_edge_line,TCE_edge_line_ia,TCE_edge_line_ib, TCE_edge_element, TCE_edge_ele_edge, TCE_edge_flux_sign]=fvcom_select_tce_edge_line(his2flux)
%
%function [TCE_edge_line,TCE_edge_line_ia,TCE_edge_line_ib, TCE_edge_element, TCE_edge_ele_edge, TCE_edge_flux_sign]=fvcom_select_tce_edge_line(his2flux)
%
%    This program calculates the geometries of all nodes, elements and TCE edge lines 
%    then select a TCE edge list that is the shortest path based on two points selected.
%    the outcome of the function is used for constructing casename_flux_map.dat which
%    is an input file for modified FVCOM to print fluxes on those TCE edges for flux analysis
%    Typical case of Example 1) below gives chn_flux_map.dat with the first row stating how many TCE edges
%    the rest of the file list element number, ia point node number, ib point node number and sign to apply
%    for aggregating total flux through the tce edge line in order to have final flux being positive towards
%    right hand side of the edge line looking from the starting edge towards the ending edge of the transect
%
%		    16
%		    441   276   281    -1
%		    441   281   277     1
%		    442   277   281    -1
%		    442   282   277     1
%		    443   277   282    -1
%		    443   282   278     1
%		    444   278   282    -1
%		    444   283   278     1
%		    445   278   283    -1
%		    445   283   279     1
%		    446   279   283    -1
%		    446   284   279     1
%		    447   279   284    -1
%		    447   284   280     1
%		    448   280   284    -1
%		    448   285   280     1
%
%    Input --- his2flux, a structure that stores various information about the fvcom model,
%              the transect and also the output file names
%    Outputs:
%
%    	TCE_edge_line(1:N_count_TCE_edge);        % TCE edge numbers for the TCE edge line transect in FVCOM
%    	TCE_edge_element(1:N_count_TCE_edge);     % Triangle element number for the TCE edge in FVCOM
%    	TCE_edge_ele_edge(1:N_count_TCE_edge);    % Triangle element edge number for the TCE edge in FVCOM
%    	TCE_edge_flux_sign(1:N_count_TCE_edge);   % sign to multiply to the flux in FVCOM to output so that output flux is
%                  	                          % positive to the right hand side of the TCE edge line
%       TCE_edge_line_ia;    % ia point of the TCE edge, it is the left side node of the triangel element edge that the TCE edge is
%                            %    intersecting with.
%       TCE_edge_line_ib;    % ia point of the TCE edge, it is the right side node of the triangel element edge that the TCE edge is
%                            %    intersecting with.
%
%    Dependency: 
%                triangle_grid_edge.m           code to calculate geometry of the grid
%                fvcom_element_area.m           code to calculate areas of the grids
%                shape_coef_gcn.m               code to calculate shape functions of the grids
%                fvcom_grid_iscw.m              code to check if the grid system is clockwisely arranged
%
%    Example 1) use netcdf, no island channel case, with salt, with tide etc
%
%     his2flux.fvcom_mat_grid_file = './CHN_HYD.mat';                         %matlab format of the model grid
%     his2flux.transname         ='cross';                                    %name of the transect
%     his2flux.xt                =[-44363,-44363];  %x coord of the starting and ending points of the transect
%     his2flux.yt                =[-539,3379];      %y coord of the starting and ending points of the transect
%     his2flux.ds_transect       = 1.0;             %step size of discretizing the transect (m)
%     his2flux.t_int             = 1;               %time record interval in processing the history outputs
%                                                   %     1 for every record,2 for every other record,...
%     his2flux.plotoutdir        =['./' his2flux.transname '/'];  %output plt dirname
%
%     [TCE_edge_line,TCE_edge_line_ia,TCE_edge_line_ib, TCE_edge_element, TCE_edge_ele_edge, TCE_edge_flux_sign]=fvcom_select_tce_edge_line(his2flux);
%     edgeline=[TCE_edge_element' TCE_edge_line_ia' TCE_edge_line_ib' TCE_edge_flux_sign']  %use this for constructing model_input/chn_flux_map.dat
%
%   Example 2)
%
%      his2flux.fvcom_mat_grid_file ='./PS2_SPB.mat';
%      his2flux.transname         ='A-A1';           %name of the transect
%      load  trans_loc_A-A1.mat;
%      his2flux.xt                =[xt(1), xt(end)] ; 
%                                        %x coord of the starting and ending points of the transect
%      his2flux.yt                =[yt(1), yt(end)] ;
%                                        %y coord of the starting and ending points of the transect
%      his2flux.ds_transect       = 20.0;           %step size of discretizing the transect (m)
%      his2flux.t_int             = 2;              %time record interval in processing the history outputs
%                                                   %     1 for every record,2 for every other record,...
%      his2flux.plotoutdir        =['./' his2flux.transname '/'];  %output plt dirname
%      [TCE_edge_line,TCE_edge_line_ia,TCE_edge_line_ib, TCE_edge_element, TCE_edge_ele_edge, TCE_edge_flux_sign]=fvcom_select_tce_edge_line(his2flux);
%      edgeline=[TCE_edge_element' TCE_edge_line_ia' TCE_edge_line_ib' TCE_edge_flux_sign']
%
%   Example 3)
%
%      his2flux.fvcom_mat_grid_file ='./PS2_SPB.mat';
%      his2flux.transname         ='B-B1';           %name of the transect
%      load  trans_loc_B-B1.mat;
%      his2flux.xt                =[xt(1), xt(end)] ;
%                                                    %x coord of the starting and ending points of the transect
%      his2flux.yt                =[yt(1), yt(end)] ;
%                                                    %y coord of the starting and ending points of the transect
%      his2flux.ds_transect       = 1.0;            %step size of discretizing the transect (m)
%      his2flux.t_int             = 2;               %time record interval in processing the history outputs
%                                                    %     1 for every record,2 for every other record,...
%      his2flux.plotoutdir        =['./' his2flux.transname '/'];  %output plt dirname
%      [TCE_edge_line,TCE_edge_line_ia,TCE_edge_line_ib, TCE_edge_element, TCE_edge_ele_edge, TCE_edge_flux_sign]=fvcom_select_tce_edge_line(his2flux);
%      edgeline=[TCE_edge_element' TCE_edge_line_ia' TCE_edge_line_ib' TCE_edge_flux_sign']
%
%
% Wen Long @ PNNL, April 27, 2014
%

%obtain input parameters from his2flux structure
%

%     if(nargin<1)
%      error(['oops,missing his2flux argument, quitting ...']);
%     end

     %debug setting
     testcase_number =  0; %1 to 16 for 16 idealized sanity test cases, 0 to turn off testing 
     %testcase_number = 22; %temporary test case for skagit river 
     if(isempty(his2flux.fvcom_mat_grid_file))
       display(['oops, input argument is empty, quitting ...']);
     end

     fieldnames_arg=fieldnames(his2flux) ;

%----check if arguments are there
     have_grid        = ~isempty(find(strcmp(fieldnames_arg,'fvcom_mat_grid_file')==1));
     have_transname   = ~isempty(find(strcmp(fieldnames_arg,'transname')==1));
     have_xt          = ~isempty(find(strcmp(fieldnames_arg,'xt')==1));
     have_yt          = ~isempty(find(strcmp(fieldnames_arg,'yt')==1));
     have_dst	      = ~isempty(find(strcmp(fieldnames_arg,'ds_transect')==1));
     have_tint        = ~isempty(find(strcmp(fieldnames_arg,'t_int')==1));
     have_pltdir      = ~isempty(find(strcmp(fieldnames_arg,'plotoutdir')==1));

     if(~have_grid)
	error(['oops, fvcom_mat_grid_file not given']);
     end
     if(~have_transname)
        error(['oops, transname not given']);
     end
     if(~have_xt || ~have_yt)
        error(['oops, xt or yt not given']);
     end
     if(~have_dst)
        error(['oops, ds_transect not given']);
     end
     if(~have_tint) 
        error(['oops, t_tint not given']);
     end
     if(~have_pltdir)
        error(['oops, plotoutdir not given']);
     end
%----

     fvcom_mat_grid_file  =his2flux.fvcom_mat_grid_file;
     transname          =his2flux.transname;
     xt                 =his2flux.xt;              
     yt                 =his2flux.yt;             
     ds_transect        =his2flux.ds_transect; 
     t_int              =int32(his2flux.t_int); 
     pltdirname         =his2flux.plotoutdir;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify below                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     yesno_tracer_flux  ='y'; %his2flux.yesno_tracer_flux;
     his_fn_int  =1;         %intervals of file processing
     i_par=1;                %choice of variable  (1 means velocity)

     debug_plot=1;
 

     if(~exist(pltdirname,'dir'))
        eval(['!mkdir -p ' pltdirname]);
     else
        eval(['!rm -r -f ' pltdirname]); 
        eval(['!mkdir -p ' pltdirname]);
     end

     if(exist(fvcom_mat_grid_file,'file'))  %if there is grid in matlab format, then load it dirrectly
               load(fvcom_mat_grid_file);   

	    
	%check if grid is clockwise or counter-clockwise
        %if counter-clockwise, then switch to clockwise 

        x_n=xyd_n(:,2);                  %x coordinate of nodes
        y_n=xyd_n(:,3);                  %y coordinate of nodes
        lon_n=lld_n(:,2);                %longitude of nodes
        lat_n=lld_n(:,3);                %latitude of nodes
        h_n=lld_n(:,4);                  %depth of nodes


       if(~fvcom_grid_iscw(e2n(:,2:4),x_n,y_n))  %if not clockwise, change to clockwise
          %switch e2n from counter-clockwise to clockwise
          e2n_tmp=e2n(:,3);e2n(:,3)=e2n(:,4);e2n(:,4)=e2n_tmp;
          clear 'e2n_tmp';
          clockwise_grid=1;
       end

     else                                 %otherwise, load separately from text files 
        display(['not able to find file named ' fvcom_mat_grid_file]);
	error(['oops, please provide matlab format of the FVCOM grid file, it can ' ...
	       'be created using sms2dm2mat.m if you don''t have it']);
     end

     NN=size(xyd_n,1) ;    %number of nodes
     NE=size(xy_e,1) ;     %number of elements

     %obtain history output file (netcdf or medm)

         %plot a grid and let user choose transect

         yesno_matlabfig=1;  %no matlab plots to be made


         %get information of the transect 

         if (isempty(xt)||isempty(yt) )
            error(['xt and yt should not be empty']);
         end

         if(length(xt(:)) ~=2 || length(yt(:)) ~=2)
            error('xt and yt should have length==2'); 
         end

	 s = [0 cumsum(sqrt(diff(xt).^2+diff(yt).^2))' ]';  %distance along the transect
	 if(ds_transect >=0.1*max(s(:)))
	   warning(['the segment length you input is too large']);
	   warning(['setting to default value of ' num2str(max(s(:))) '/100']);
	   ds_transect= max(s(:))/100.0; 
	 end
	 dx=ds_transect;
	 st=[s(1):dx:max(s) max(s)];
	 xp=spline(s,xt,st);                 %sample the transect (proliferation) using ds_transect step size
	 yp=spline(s,yt,st);
	 save(['trans_loc_' transname '.mat'],'xp','yp','xt','yt');
	 Ntrans=length(st);    %number of points to describe the transect
         %plot the transect if a plot has not been made
	 make_transectplot=~exist([ './' transname '_transect.png'],'file');  
	 if(make_transectplot)
	    hmap=figure('visible','off','position',[1,1,1000,1000]);clf
	    set(gcf,'renderer','zbuffer');
	    tricolor(e2n(:,2:4),x_n,y_n,h_n);
	    colorbar;
	    shading interp; hold on;
	    xlabel('x')
	    ylabel('y')
	    %plot the chosen transect using a red line
	    plot(xp,yp,'k.',xt,yt,'r*')
	    grid
	 end
         %chop off transect parts that have nan coordinates
         %the way to do it is to mark xp and yp values as nans if they do not fall in any
         %of the model elements (find the closest element, and check if xp, yp is in the closest element)
         %for those outside of the domain plot them as being gray 
         display(['Plotting transect ...']); 
         x_e=mean(x_n(e2n(:,2:4)),2);  %element coordinate
         y_e=mean(y_n(e2n(:,2:4)),2);
         NE_use=0;  %Number of elements to be used
         NN_use=0;  %Number of nodes to be used
         IE_use=[]; %element ID list to be used
         IN_use=[]; %node ID list to be used
         for i=1:length(xp)
	     isgood(i)=0; %set it to bad first 
   	     for j=1:NE   %check elements one by one
	       %this inpolygon is kind of slow
	       %     if(inpolygon(xp(i),yp(i),x_n(e2n(j,2:4)),y_n(e2n(j,2:4)))==1)
	       %     check if xp, yp is in the any of element
	       %         isgood(i)=1;
	       %     end
	       %use a faster InPolygon() c function instead:
	       [in_on,on,strict] = InPolygon(xp(i),yp(i),x_n(e2n(j,2:4)),y_n(e2n(j,2:4)));
	       if(in_on==1)   %if in or on an element j, then mark as good and 
			      %record the element j and nodes e2n(j,2:4)
		   isgood(i)=1;  %
		   NE_use=NE_use+1;
		   IE_use(NE_use)=j;  %record element number
		   NN_use=NN_use+3;   
		   IN_use(NN_use-2:NN_use)=e2n(j,2:4);  %record nodes of element(j)
		   break;             %quit the j loop once found it good
	       end
	     end
        end
        Ipath_out=find(isgood==0);   %out of domain
        Ipath_in =find(isgood==1);   %in domain
        clear 'distance' 'dsort' 'Isort' 'dmin' 'Imin' 'isgood';

       make_transectplot=~exist([ './' transname '_transect.png'],'file');
       if(make_transectplot)
          plot(xp(Ipath_out),yp(Ipath_out),'o','markerfacecolor',[0.6,0.6,0.6],'markersize',15)   %places out of domain
       end

    %
    %find out unique list of elements and nodes to be used for data processing by sprawling out from existing list
    %

    %--------------
	Nsprawl=1; %rounds of sprawling is set to 1 (can increase if you want to include more data)

	for isprawl=1:Nsprawl
	   %sort out the IE_use and IN_use so that they are unique (to avoid duplicated data points)
	   IE_use_uniq=unique(IE_use);    %these will be the element numbers to be used
	   IN_use_uniq=unique(IN_use);    %these will be the node numbers to be used
	   NE_use_uniq=length(IE_use_uniq);
	   NN_use_uniq=length(IN_use_uniq);

	   %include elements that are connected to the nodes in the list and add them to the list
	   for i=1:NN_use_uniq
	       IE_more= find(e2n(:,2) == IN_use_uniq(i) | e2n(:,3) == IN_use_uniq(i) | e2n(:,4) == IN_use_uniq(i)) ;
	       NE_use = NE_use+length(IE_more);     %add more elements to the list
	       IE_use(NE_use-length(IE_more)+1:NE_use) = IE_more;
	       NN_use = NN_use+length(IE_more)*3;   %add more nodes to the list (nodes of new elements added)
	       IN_use(NN_use-length(IE_more)*3+1:NN_use)=reshape(e2n(IE_more,2:4),[1,length(IE_more)*3]);
	   end
       end

       %final list of unique elements and nodes to be used for data processing
	IE_use_uniq=unique(IE_use);    %these will be the element numbers to be used
	IN_use_uniq=unique(IN_use);    %these will be the node numbers to be used
	NE_use_uniq=length(IE_use_uniq);
	NN_use_uniq=length(IN_use_uniq); 

       clear 'IE_use' 'IN_use' 'NE_use' 'NN_use'; 


       %find the orientation of the path (tangential vector (tx_p, ty_p) and normal vector (nx_p, ny_p))
       %find angle of path (angle_p) relative to x coorinate

       angle_trans=zeros([1, length(yp)]);
       %find angle of transect to East
       angle_trans(1:end-1)=atan2(diff(yp),diff(xp)); %angle of transect line
       angle_trans(end)=angle_trans(end-1); %repeat last one

       tx_p=cos(angle_trans);   %tangential vector along the transect in (x,y) coordinate
       ty_p=sin(angle_trans);
       tx_p(end)=tx_p(end-1);   %repeat last one because last point on transect has same coord as next to last one
       ty_p(end)=ty_p(end-1);

       nx_p=cos(angle_trans+pi/2.0); %normal vector of the transect in (x,y) coordinate
				  %(90 deg counter-clociwise from tangential direction )
       ny_p=sin(angle_trans+pi/2.0);
       nx_p(end)=nx_p(end-1);   %repeat last one because last point on transect has same coord as next to last one
       ny_p(end)=ny_p(end-1);

       if(make_transectplot)
          %show the tangential and normal vectors on the path
          quiver(xp(1:5:end),yp(1:5:end),tx_p(1:5:end),ty_p(1:5:end),'m-','MaxHeadSize',1);  %along the transect  
          quiver(xp(1:5:end),yp(1:5:end),nx_p(1:5:end),ny_p(1:5:end),'k-','MaxHeadSize',1);  %normal to the transect

          trans_width=max(max(xp)-min(xp),max(yp)-min(yp));  %rectangle size of xp, yp domain
     
          %this makes sure axis is shown equal (square box) so that it is not distorted and still
          %has the whole transect in figure 

          axis([mean(xp(:))-1.2*trans_width/2,mean(xp(:))+1.2*trans_width/2,...
   	        mean(yp(:))-1.2*trans_width/2,mean(yp(:))+1.2*trans_width/2]);
         
          title('transect plot'); 

          %save the map and the transect figure
          saveas(gcf,[pltdirname './' transname '_transect.eps'], 'epsc2');
          saveas(gcf,[pltdirname './' transname '_transect.png'], 'png');
          saveas(gcf,[pltdirname './' transname '_transect.fig'], 'fig');

          delete(hmap);
       end

       if(isempty(Ipath_in))
   	   error(['Oops all transect points are out of model domain, quitting ...']);
       end

   %
   %---find out the TCE's and geometries of the model
   %

     nv=e2n(:,2:4);     
     vx=xyd_n(:,2);  %node x 
     vy=xyd_n(:,3);  %node y
     xc=xy_e(:,1);   %element centroid x
     yc=xy_e(:,2);   %element centroid y

     hc=nan*zeros([size(e2n,1),1]);  %depth at element centroid
     for ie=1:size(xy_e,1)
          hc(ie,1)=mean(h_n(e2n(ie,2:4)));  
     end

     %---debug--assuming no open boundary nodes
     iobcn=0;
     i_obc_n=[];
     %---debug--;

     tri=nv;
     spherical=0;

     [iec,ienode,isbc,isbce,isonb,lisbce_1,lisbce_2,mx_nbr_elem,ne,nbe,nbsn,nbve,nbvt,niec,...
                    nisbce_1,nisbce_2,nisbce_3,ntrg,ntsn,ntve,ncv,ncv_i,sitac,sitae,sitau,...
         deltux,deltuy,dltxc,dltxe,dltxyc,dltxye,dltyc,dltye,dltxne,dltyne,xijc,xije,yijc,...
                                                                          yije,epor,nflag ...
      ] =  triangle_grid_edge(nv,vx,vy,xc,yc,iobcn,i_obc_n,spherical);

     %calculate the areas of elements and around nodes
     [art_ele art_tce art_can]=fvcom_element_area(spherical,ntve,nv,vx,vy,isonb,nbve,nbvt,xc,yc);

     %
     %show the boundary elements' boundary edge angle 
     %

     if(clockwise_grid)

	[alpha,a1u,a2u,aw0,awx,awy]=shape_coef_gcn(vx,vy,nv,xc,yc,art_ele,ntve,nbve,nbe,isbce,isonb,spherical);
 	%only good for clockwise grid system

     else   %WLONG: the upper shape_coef_gcn() assumes clockwise grid system
	error(['oops, shape function only works with clockwised grid system, need to develop formula here for counter-clockwised grids']);

     end


     %
     %connstruct a large matrix for Graph showing distance between any two elements a and b, connected by TCE edges
     %

%
%    GG=nan*zeros([NE,NE]);  %NE by NE array
%
%    for iele=1:NE   %loop through all elements
%
%             GG(iele,iele)=0.0 ;  %no distance from self to self
%
%             for j=1:3
%                 ise=nbe(iele,j) ; %neighbor element
%                 if (ise~=0)      %make sure neighbor elemnt is valid
%			for iedge=1:ne  %loop through triange edges
%			    if(  iec(iedge,1)==iele && iec(iedge,2)==ise ||  ...   %find the edge that has elements iele and ise as neighbors
%                                 iec(iedge,1)==ise  && iec(iedge,2)==iele )         %
%				 %get the distance from element iele to middle point of the edge
%                                  r1=sqrt( (xc(iele)-xijc(iedge))^2  ...
%                                          +(yc(iele)-yijc(iedge))^2 ) ; %where xijc, yijc is edge center location
%                                                                        %xc,yc is element center location
%                                 %get the distance from element ise to middle point of the edge
%
%                                  r2=sqrt( (xc(ise) - xijc(iedge))^2 ...
%				          +(yc(ise) - yijc(iedge))^2)
%                                  GG(iele,ise)=r1+r2;
%                                  GG(ise,iele)=r1+r2;
%                            else
%                            end
%			end
%                 else
%		     %do nothing 
%                 end
%             end
%     end
%

     
     %
     %Do the above of connectivity using a different method
     %

     D_elem_list_1=nan*zeros([1,ne]);  %list of elements that are connected (neighboring) to the following list 
     D_elem_list_2=nan*zeros([1,ne]);  %list of destination elements 
     D_edge_list  =nan*zeros([1,ne]);
     W_list=zeros([1,ne]);

     N_edge_connect=0;

     for iedge=1:ne  %loop through all triangle edges
          if(isbc(iedge)==0)      %find the two connecting elements for interior triangle element edges only

                N_edge_connect=N_edge_connect+1;

                % get the left side element and right side element

                ie_left=iec(iedge,2);
                ie_right=iec(iedge,1);

                %length of TCE edge on the left side of the edge 
                r1 = sqrt( (xc(ie_left)-xijc(iedge))^2 ...
                          +(yc(ie_left)-yijc(iedge))^2);
               
                %length of TCE edge on the right side of the edge
                r2 = sqrt( (xc(ie_right) -xijc(iedge))^2 ...
                          +(yc(ie_right) -yijc(iedge))^2);

                %record the connecting elements and weight (length)

                D_elem_list_1(N_edge_connect)= ie_left;
                D_elem_list_2(N_edge_connect)= ie_right;
                D_edge_list(N_edge_connect)  = iedge;
                W_list(N_edge_connect) = r1+r2;
          end
     end
    

    %get the starting and ending element to construct the shortest path connected by an edge line

    if(isempty(xt))
	 hmap=figure('visible','on');clf
	 tricolor(e2n(:,2:4),x_n,y_n,h_n);
	 set(gcf,'renderer','zbuffer');
	 colorbar;
	 shading interp; hold on;
	 xlabel('x')
	 ylabel('y')
	 axis equal;
	 xt=[];
	 yt=[];
	 more_transect_points=true;
	 while(more_transect_points)
	    disp('Zoom in as you like, press enter after zoom and then click on Way Points of transect :')
	    title({'Map to choose transect from.'; ...
		   'Zoom in as you like,  press enter and then left-click on Way Points of transec'; ...
		   'right click to come back to zoom mode again';'middle click to finish choosing points'}, ...
		   'color',[0 0 0],'fontsize',25);
	     %---allow zoom before clicking
	     zoom on; % use mouse button to zoom in or out
		   % Press Enter to get out of the zoom mode.
		   % CurrentCharacter contains the most recent key which was pressed after opening
		   % the figure, wait for the most recent key to become the return/enter key
	     waitfor(gcf,'CurrentCharacter',13);
	     zoom reset;
	     zoom off;
	     while(true)
		 [xt_tmp,yt_tmp,buttonPress]=ginput(1);
		 switch(buttonPress)
		       case 3      %right click will quit this one
			    break;
		       case 2      %quit completely mid-button click
			    more_transect_points=false;
			    break;
		       case 1
			    plot(xt_tmp,yt_tmp,'r*');
			    xt=[xt ;xt_tmp];   %add this point xt_tmp, yt_tmp
			    yt=[yt ;yt_tmp];

                            if(length(xt)==2)
                              more_transect_points=false; 
                            end
		 end
	     end
	     tmp_input=input('please input 1 to end this,0 to continue:')  %WLong: had to repeat twice !! to make it work
	     tmp_input=input('please input 1 to end this,0 to continue:')  %the first time is skipped by matlab by right-clicking

	     if (tmp_input==1)
		more_transect_points=false;
                delete(hmap);
	     else
		delete(hmap); %replot and contiune choosing points
		hmap=figure('visible','on');clf
		set(gcf,'renderer','zbuffer');
		tricolor(e2n(:,2:4),x_n,y_n,h_n);
		colorbar;
		shading interp; hold on;
		plot(xt,yt,'r*');
		axis equal;
	     end
	 end

    end

    %     
    %interpolate through the xt, yt curve and get starting point and ending point that are in the domain
    %and TCE edge line would be based on these two points

         s = [0 cumsum(sqrt(diff(xt).^2+diff(yt).^2))' ]';  %distance along the transect
         if(ds_transect >=0.1*max(s(:)))
           warning(['the segment length you input is too large']);
           warning(['setting to default value of ' num2str(max(s(:))) '/100']);
           ds_transect= max(s(:))/100.0;
         end
         dx=ds_transect;
         st=[s(1):dx:max(s) max(s)];
         xp=spline(s,xt,st);                 %sample the transect (proliferation) using ds_transect step size
         yp=spline(s,yt,st);
         save(['trans_loc_' transname '.mat'],'xp','yp','xt','yt');
         Ntrans=length(st);    %number of points to describe the transect
         %plot the transect if a plot has not been made
         make_transectplot=~exist([ './' transname '_transect.png'],'file');
         make_transectplot=1;
         if(make_transectplot)
            hmap=figure('visible','off','position',[1,1,1000,1000]);clf
            set(gcf,'renderer','zbuffer');
            tricolor(e2n(:,2:4),x_n,y_n,h_n);
            colorbar;
            shading interp; hold on;
            xlabel('x')
            ylabel('y')
            %plot the chosen transect using a red line
            plot(xp,yp,'k.',xt,yt,'r*')
            grid
 
            saveas(gcf,[pltdirname '/' 'selected_transect.fig'],'fig');

            delete(hmap);
         end

    %chop off transect parts that have nan coordinates
    %the way to do it is to mark xp and yp values as nans if they do not fall in any
    %of the model elements (find the closest element, and check if xp, yp is in the closest element)
    %for those outside of the domain plot them as being gray

    display(['Plotting transect ...']);

    %plot the start and end points of the edge line and 
    %highlight the elemnts that are chosen
       hmap=figure('visible','off');clf

         set(gcf,'renderer','zbuffer');
         tricolor(e2n(:,2:4),x_n,y_n,h_n);
         colorbar;
         shading interp; hold on;
         plot(xt,yt,'r*');
         axis equal;
     
 	 %plot the chosen transect using a red line
	 plot(xp,yp,'rs');
   	 grid;

         x_e=mean(x_n(e2n(:,2:4)),2);  %element coordinate
         y_e=mean(y_n(e2n(:,2:4)),2);

         NE_use=0;  %Number of elements to be used
         NN_use=0;  %Number of nodes to be used
         IE_use=[]; %element ID list to be used
         IN_use=[]; %node ID list to be used

         for i=1:length(xp)
             isgood(i)=0; %set it to bad first
             for j=1:NE   %check elements one by one

               %this inpolygon is kind of slow
               %     if(inpolygon(xp(i),yp(i),x_n(e2n(j,2:4)),y_n(e2n(j,2:4)))==1)
               %     check if xp, yp is in the any of element
               %         isgood(i)=1;
               %     end

               %use a faster InPolygon() c function instead:
               [in_on,on,strict] = InPolygon(xp(i),yp(i),x_n(e2n(j,2:4)),y_n(e2n(j,2:4)));
               if(in_on==1)   %if in or on an element j, then mark as good and
                              %record the element j and nodes e2n(j,2:4)
                   isgood(i)=1;  %
                   NE_use=NE_use+1;
                   IE_use(NE_use)=j;  %record element number
                   NN_use=NN_use+3;
                   IN_use(NN_use-2:NN_use)=e2n(j,2:4);  %record nodes of element(j)
                   break;             %quit the j loop once found it good
               end
             end
         end

         Ipath_out=find(isgood==0);   %out of domain
         Ipath_in =find(isgood==1);   %in domain

         plot(x_e(IE_use'),y_e(IE_use'),'ro');
         %plot(x_n(IN_use'),y_n(IN_use'),'m.'); 

         %show the elememnts in blue color (start and ending element of the transect)

         patch(x_n(e2n(IE_use(1)  ,2:4)),y_n(e2n(IE_use(1)  ,2:4)),'b');  
         patch(x_n(e2n(IE_use(end),2:4)),y_n(e2n(IE_use(end),2:4)),'b');
         
         saveas(gcf,[pltdirname '/' 'start_end_element.fig']);


       delete(hmap);


     %calculate the graph of shortest path
     %make sparce matrix
%     GG_sparse=sparse(D_elem_list_1(1:N_edge_connect),D_elem_list_2(1:N_edge_connect),W_list(1:N_edge_connect), NE,NE);
%    [ddist1,ppath,ppred] = graphshortestpath(GG_sparse,IE_use(1),IE_use(end));

%    %switch the source and destination elements, see if the search results are different
%   GG_sparse=sparse(D_elem_list_2(1:N_edge_connect),D_elem_list_1(1:N_edge_connect),W_list(1:N_edge_connect), NE,NE);
%     [ddist2,ppath,ppred] = graphshortestpath(GG_sparse,IE_use(end),IE_use(1));

%    %reconstruct the sparce matrix with both direction
%   %switch the source and destination elements, see if the search results are different

    GG_sparse=sparse([D_elem_list_1(1:N_edge_connect) D_elem_list_2(1:N_edge_connect)], ...
                     [D_elem_list_2(1:N_edge_connect) D_elem_list_1(1:N_edge_connect)], ...
                     [W_list(1:N_edge_connect)        W_list(1:N_edge_connect)], NE,NE);
     [ddist3,ppath,ppred] = graphshortestpath(GG_sparse,IE_use(1),IE_use(end));

   %plot the found path
   hmap=figure('visible','off');clf
       set(gcf,'renderer','zbuffer');
       tricolor(e2n(:,2:4),x_n,y_n,h_n);
       colorbar;
       shading interp; hold on;
       plot(xt,yt,'r*');
       plot(xp,yp,'k*');
       axis equal;

       %show the elememnts in magenta color (start and ending element of the transect)

        patch(x_n(e2n(IE_use(1)  ,2:4)),y_n(e2n(IE_use(1)  ,2:4)),'m');
        patch(x_n(e2n(IE_use(end),2:4)),y_n(e2n(IE_use(end),2:4)),'m');

       for i=2:length(ppath(:))-1;
           iele=ppath(i); %find the element on the path
           if(i==2)
                  patch(x_n(e2n(iele,2:4)),y_n(e2n(iele,2:4)),'r');  %starting
           elseif(i==length(ppath(:))-1)
		  patch(x_n(e2n(iele,2:4)),y_n(e2n(iele,2:4)),'k');
           else
                  patch(x_n(e2n(iele,2:4)),y_n(e2n(iele,2:4)),'y');
           end

       end
      
       title('found path (start from red element, end at yellow element)')
       saveas(gcf,[pltdirname '/' 'found_shortest_element_path.fig'],'fig'); 

   delete(hmap)

   %complete the search of edgeline by listing all the elements and TCE edges from the first to the last element
   %also prepend and append the first and last element if the first or last element is a boundary element

   NE_path=length(ppath(:));  %number of elements on the path
   TCE_edge_line=[];
   TCE_edge_line_ia=[];    %LHS node of the TCE edge 
   TCE_edge_line_ib=[];    %RHS node of the TCE edge     
   TCE_edge_element=[];
   TCE_edge_ele_edge=[];
   TCE_edge_flux_sign=[];
 
   N_count_TCE_edge=0;
   for ie_path=1:NE_path

       %this element and the next element on the path
       iele=ppath(ie_path);      %get the element number

       if(ie_path<NE_path)
          iele_next=ppath(ie_path+1);
       else
          iele_next=0;
       end

       if(ie_path==1)
                %if first element is a boundary element, complete it with a TCE edge that connects to the boundary
                %if the first element is not a boundary element, then search for the TCE edges that connects it to the next element

                if(any(nbe(iele,1:3)==0))                                  %check if any neighbor of element iele is non-existant, 
                                                                           %which means the element iele is on the boundary
                                                                           %find a TCE edge that is in this boundary element, 
                                                                           %and it is NOT the TCE edge that connects this element 
                                                                           %iele with the next element iele_next

                      itce_edge_find=find(ntrg==iele);                     %all 3 tce edges in element iele
                      iele_edge_find=find(iec(:,1)==iele & iec(:,2)==0);   %all boundary triangle edges of this element iele

                      if(length(iele_edge_find>=2))                        %if the element has more than one boundary edges (max 2) 
                                                                           %(3 is isolated element, not allowed in FVCOM)
                         iele_edge_find=iele_edge_find(1);  %only take one of them, for it would be working the same with the other one
                      end 

                      for j=1:length(itce_edge_find)  %find the tce edge that connects to the boundary
 
                          idone=0;

                          if(niec(itce_edge_find(j),1) == ienode(iele_edge_find,1)      && ...  %ia is starting node of the element boundary edge
                             niec(itce_edge_find(j),2) == ienode(iele_edge_find,2))             %ib is ending node of the element boundary edge

                              N_count_TCE_edge=N_count_TCE_edge+1;

                              TCE_edge_line(N_count_TCE_edge)      = itce_edge_find(j);         %get the TCE edge index
                              TCE_edge_flux_sign(N_count_TCE_edge) = -1;                        %positive as to the right hand 
                                                                                                % side of the TCE edge line
                              TCE_edge_element(N_count_TCE_edge)   = iele;                      %index of triangle element that the TCE edge is in 
                              TCE_edge_ele_edge(N_count_TCE_edge)  = iele_edge_find;            %index of triangle element edge
 
                              TCE_edge_line_ia(N_count_TCE_edge)=niec(itce_edge_find(j),1);
                              TCE_edge_line_ib(N_count_TCE_edge)=niec(itce_edge_find(j),2);
                              idone=1;
                              break;  
                          end

                          if(idone==0)
                          if(niec(itce_edge_find(j),1) == ienode(iele_edge_find,2)      && ...   %ia is ending node of the element boundary edge
                             niec(itce_edge_find(j),2) == ienode(iele_edge_find,1))              %ib is starting node of the element boundary edge

                              N_count_TCE_edge=N_count_TCE_edge+1;

                              TCE_edge_line(N_count_TCE_edge)=itce_edge_find(j);     %get the TCE edge index
                              TCE_edge_flux_sign(N_count_TCE_edge) =  -1;  %         %positive as to the right hand side of the TCE edge line
                              TCE_edge_element(N_count_TCE_edge) = iele;             %index of triangle element that the TCE edge is in
                              TCE_edge_ele_edge(N_count_TCE_edge) = iele_edge_find;  %index of triangle element edge

                              TCE_edge_line_ia(N_count_TCE_edge)=niec(itce_edge_find(j),1);
                              TCE_edge_line_ib(N_count_TCE_edge)=niec(itce_edge_find(j),2);
                              break;
                          end
                          end
                      end
                      end

       end

       %regular intenal TCE edges (TCE edges that do not intercept with bounary edges of the domain)
 

       if(ie_path<NE_path)
             if(iele ~=0 && iele_next ~=0)

                 %regular TCE edge
                 %find the element edge r with iec(r,1)=iele and iec(r,2) =iele_next
                 
                 iele_edge_find=find(iec(:,1)==iele & iec(:,2)==iele_next | iec(:,2)==iele & iec(:,1)==iele_next);

                 if(length(iele_edge_find)~=1)
                   error('oops, not able to find a uniqu element edge that connects iele and iele_next'); 
                 end

                 itce_edge_find=find(ntrg==iele);                                % find all tce edges within this element iele
                 

                 %find the TCE edge within  this element that connects with the 
                 %shared element edge between iele and iele_next
                 ifound=0;
                 for j=1:length(itce_edge_find)                                  % check ia ib of the tce_edge to find the one that 
                                                                                 % had the element edge being r
                     if(ifound==0)
                       if(    niec(itce_edge_find(j),1)==ienode(iele_edge_find,1) && niec(itce_edge_find(j),2)==ienode(iele_edge_find,2)  ...
                          ||  niec(itce_edge_find(j),1)==ienode(iele_edge_find,2) && niec(itce_edge_find(j),2)==ienode(iele_edge_find,1))

                         N_count_TCE_edge=N_count_TCE_edge+1;
                         TCE_edge_line(N_count_TCE_edge)=itce_edge_find(j);      % get the TCE edge index
                         TCE_edge_flux_sign(N_count_TCE_edge) =  +1;             % positive as to the right hand side of the TCE edge line
                         TCE_edge_element(N_count_TCE_edge) = iele;
                         TCE_edge_ele_edge(N_count_TCE_edge) = iele_edge_find;
                        
                         TCE_edge_line_ia(N_count_TCE_edge)=niec(itce_edge_find(j),1);
                         TCE_edge_line_ib(N_count_TCE_edge)=niec(itce_edge_find(j),2);
                         ifound=1;
                       end
                     end
                 end

                 %find the TCE edge in the next element
                 itce_edge_find=find(ntrg==iele_next);  %find all tce edges within element iele_next
                 ifound=0;
                 for j=1:length(itce_edge_find)
                    if(ifound==0)
                       if(    niec(itce_edge_find(j),1)==ienode(iele_edge_find,1) && niec(itce_edge_find(j),2)==ienode(iele_edge_find,2) ...
                          ||  niec(itce_edge_find(j),1)==ienode(iele_edge_find,2) && niec(itce_edge_find(j),2)==ienode(iele_edge_find,1))
                              N_count_TCE_edge=N_count_TCE_edge+1;
                              TCE_edge_line(N_count_TCE_edge)=itce_edge_find(j);  % get the TCE edge index
                              TCE_edge_flux_sign(N_count_TCE_edge) =  -1;         % positive as to the right hand side of the TCE edge line
                              TCE_edge_element(N_count_TCE_edge)   = iele_next;
                              TCE_edge_ele_edge(N_count_TCE_edge)  = iele_edge_find;
                               
                              TCE_edge_line_ia(N_count_TCE_edge)=niec(itce_edge_find(j),1);
                              TCE_edge_line_ib(N_count_TCE_edge)=niec(itce_edge_find(j),2);
                              ifound=1;
                      end
                    end
                 end
             end
       end

       if(ie_path==NE_path)

              %if the last element is a boundary element, complete it with a TCE edge that connects to the boundary
              if(any(nbe(iele,1:3)==0))

                  itce_edge_find=find(ntrg==iele);                              % all tce edges in the ending elment iele

                  iele_edge_find=find(iec(:,1)==iele & iec(:,2)==0 | iec(:,1)==0 & iec(:,2)==iele);     
                                                                                % all  boundary edges of the element

                  if(length(iele_edge_find>=2))                                 % if the element has more than one boundary edges (max 2)
                                                                                % (3 is isolated element, not allowed in FVCOM)
                     iele_edge_find=iele_edge_find(1);  %only take one
                  end

                  %find the TCE edge that intersects with the boundary 

                  for j=1:length(itce_edge_find)                                % find the tce edges that connects to the boundary

                     idone=0;
                     if(niec(itce_edge_find(j),1) == ienode(iele_edge_find,1) && ...
                        niec(itce_edge_find(j),2) == ienode(iele_edge_find,2))

                        N_count_TCE_edge                     = N_count_TCE_edge+1;
                        TCE_edge_line(N_count_TCE_edge)      = itce_edge_find(j);     % get the TCE edge index
                        TCE_edge_flux_sign(N_count_TCE_edge) = +1;                    % positive as to the right hand side of the TCE edge line
                        TCE_edge_element(N_count_TCE_edge)   = iele;
                        TCE_edge_ele_edge(N_count_TCE_edge)  = iele_edge_find;

                        TCE_edge_line_ia(N_count_TCE_edge)=niec(itce_edge_find(j),1); 
                        TCE_edge_line_ib(N_count_TCE_edge)=niec(itce_edge_find(j),2); 
                        idone=1;
                        break;
                     end

                     if(idone==0)
                        if(niec(itce_edge_find(j),1) == ienode(iele_edge_find,2) && ...
                           niec(itce_edge_find(j),2) == ienode(iele_edge_find,1))

                           N_count_TCE_edge                     = N_count_TCE_edge+1;
                           TCE_edge_line(N_count_TCE_edge)      = itce_edge_find(j);     % get the TCE edge index
                           TCE_edge_flux_sign(N_count_TCE_edge) = +1;                    % positive as to the right hand side of the TCE edge line
                           TCE_edge_element(N_count_TCE_edge)   = iele;
                           TCE_edge_ele_edge(N_count_TCE_edge)  = iele_edge_find;

                           TCE_edge_line_ia(N_count_TCE_edge)=niec(itce_edge_find(j),1);
                           TCE_edge_line_ib(N_count_TCE_edge)=niec(itce_edge_find(j),2);
                           idone=1;
                           break;
                        end
                     end

                  end
              end
       end
   end

   if(debug_plot) %debug plot all the active TCE edges that will have flux values calculated
          hfig_tce_use_6=figure('position',[1,1, 800, 800],'visible','off'); hold on;
                  %plot the mesh
                  trimesh(tri,vx,vy,'color',[0 0 0]);
                  %show the TCEs
                  line(xije',yije','color',[0 1 1],'linestyle','--');
                  %show the TCE vertices
                  plot(xije(:,2),yije(:,2),'k.');  %edge mid points
                  plot(xije(:,1),yije(:,1),'ko');  %element centers
                  %plot the chosen transect using a red line
                  plot(xt,yt,'r*');
                  grid;
                  %plot the found path
                  for i=1:length(ppath(:));
                        iele=ppath(i); %find the element on the path
                        if(i==1)
	                    patch(x_n(e2n(iele,2:4)),y_n(e2n(iele,2:4)),'r');
                        elseif(i==length(ppath(:)))
                            patch(x_n(e2n(iele,2:4)),y_n(e2n(iele,2:4)),'k');

                        else
   			    patch(x_n(e2n(iele,2:4)),y_n(e2n(iele,2:4)),'y');
                        end
                  end
                  for j=1:N_count_TCE_edge
                      plot(xije(TCE_edge_line(j),1:2),yije(TCE_edge_line(j),1:2),'r-');  %use red
                      text_x=mean(xije(TCE_edge_line(j),1:2));
                      text_y=mean(yije(TCE_edge_line(j),1:2));
                      text(text_x,text_y,num2str(TCE_edge_flux_sign(j)),'color',[0 1 1],'fontsize',15);
                  end

                  title('TCE edges used for the edge line and xflux sign (route from red to black)');
 
		  trans_width=max(max(xp)-min(xp),max(yp)-min(yp));  %rectangle size of xp, yp domain

	          %this makes sure axis is shown equal (square box) so that it is not distorted and still
	          %has the whole transect in figure

	          axis([mean(xp(:))-1.2*trans_width/2,mean(xp(:))+1.2*trans_width/2,...
                  	mean(yp(:))-1.2*trans_width/2,mean(yp(:))+1.2*trans_width/2]);
         
          saveas(gcf,[pltdirname '/' 'hfig_tce_use_6.fig'],'fig');
          delete(hfig_tce_use_6);
    end


 if(debug_plot) %debug plot all the active TCE edges that will have flux values calculated
          hfig_tce_use_7=figure('position',[1,1, 800, 800],'visible','off'); hold on;
                  %plot the mesh
                  trimesh(tri,vx,vy,'color',[0 0 0]);
                  %show the TCEs
                  line(xije',yije','color',[0 1 1],'linestyle','--');
                  %show the TCE vertices
                  plot(xije(:,2),yije(:,2),'k.');  %edge mid points
                  plot(xije(:,1),yije(:,1),'ko');  %element centers
                  %plot the chosen transect using a red line
                  plot(xt,yt,'r*');
                  grid;

                  %plot the found path
                  for i=1:length(ppath(:));
                        iele=ppath(i); %find the element on the path
                        if(i==1)
                            plot(mean(x_n(e2n(iele,2:4))),mean(y_n(e2n(iele,2:4))),'rs');
                        elseif(i==length(ppath(:)))
                            plot(mean(x_n(e2n(iele,2:4))),mean(y_n(e2n(iele,2:4))),'ks');
                        else
                            plot(mean(x_n(e2n(iele,2:4))),mean(y_n(e2n(iele,2:4))),'ys');

                        end
                  end

              
                  %plot the TCE edges that are found for the transect

                  for j=1:N_count_TCE_edge
                      plot(xije(TCE_edge_line(j),1:2),yije(TCE_edge_line(j),1:2),'r-');  %use red

                      text_x=mean(xije(TCE_edge_line(j),1:2));
                      text_y=mean(yije(TCE_edge_line(j),1:2));
                      text(text_x,text_y,num2str(TCE_edge_flux_sign(j)),'color',[0 1 1],'fontsize',15);
%----
                      itce_edge=TCE_edge_line(j);
                      i=itce_edge;
                      ia=niec(i,1) ; %node of ia (right hand side node of tce edge)
                      ib=niec(i,2) ; %node of ib (left  hand side node of tce_edge)

                      %label ia, ib in the triangle (ia,ib, (xije(i,1),yije(i,1)))

                      ia_x=vx(ia);
                      ia_y=vy(ia);

                      ib_x=vx(ib);
                      ib_y=vy(ib);

                      ic_x=xije(i,1);
                      ic_y=yije(i,1);

                      mid_x=(ib_x+ic_x)/2;
                      mid_y=(ib_y+ic_y)/2;
                      mid_angle=atan2(mid_y-ia_y,mid_x-ia_x);
                      mid_length=sqrt((mid_y-ia_y)^2+(mid_x-ia_x)^2);
                      text_x=ia_x+0.1*mid_length*cos(mid_angle);
                      text_y=ia_y+0.1*mid_length*sin(mid_angle);
                      text_angle=mod((mid_angle+360),360);
                      text(text_x,text_y,'ia');

                      mid_x=(ia_x+ic_x)/2;
                      mid_y=(ia_y+ic_y)/2;
                      mid_angle=atan2(mid_y-ib_y,mid_x-ib_x);
                      mid_length=sqrt((mid_y-ib_y)^2+(mid_x-ib_x)^2);
                      text_x=ib_x+0.1*mid_length*cos(mid_angle);
                      text_y=ib_y+0.1*mid_length*sin(mid_angle);
                      text_angle=mod((mid_angle+360),360);
                      text(text_x,text_y,'ib');

                      %also show which direction (normal vector of the tce is (toward ib))
                      n_x = dltye(i);
                      n_y = -dltxe(i);

                      mid_x = (xije(i,1)+xije(i,2))/2.0;
                      mid_y = (yije(i,1)+yije(i,2))/2.0;

                      quiver(mid_x,mid_y,n_x,n_y,'g');
%-----
                  end


                  title({'TCE edges used for the edge line and xflux sign'; ...
                         '(route from red to black, +1 for RHS of transect is same as arrow(ib), -1 otherwise )'});

                  trans_width=max(max(xp)-min(xp),max(yp)-min(yp));  %rectangle size of xp, yp domain

                  %this makes sure axis is shown equal (square box) so that it is not distorted and still
                  %has the whole transect in figure

                  axis([mean(xp(:))-1.3*trans_width/2,mean(xp(:))+1.2*trans_width/3,...
                        mean(yp(:))-1.3*trans_width/2,mean(yp(:))+1.2*trans_width/3]);

                  saveas(gcf,[pltdirname '/' 'hfig_tce_use_7.fig'],'fig');
                  delete(hfig_tce_use_7);
    end


    %finally output TCE_edge_line, TCE_edge_flux_sign

    TCE_edge_line;        % TCE edge numbers for the TCE edge line transect in FVCOM
    TCE_edge_element;     % Triangle element number for the TCE edge in FVCOM
    TCE_edge_ele_edge;    % Triangle element edge number for the TCE edge in FVCOM
    TCE_edge_flux_sign;   % sign to multiply to the flux in FVCOM to output so that output flux is
                          % positive to the right hand side of the TCE edge line
    TCE_edge_line_ia;     % ia point of the TCE edge, it is the left side node of the triangel element edge that the TCE edge is 
                          %    intersecting with.
    TCE_edge_line_ib;     % ia point of the TCE edge, it is the right side node of the triangel element edge that the TCE edge is
                          %    intersecting with.

end  %end of function

