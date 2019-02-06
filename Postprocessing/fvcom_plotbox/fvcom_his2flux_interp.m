function fvcom_his2flux_interp(his2flux)
%
%function fvcom_his2flux_interp(his2flux)
%
%  This function reads fvcom model outputs and then calculates volume flux and tracer flux across a transect 
% using interpolations
%
%    Input --- his2flux, a structure that stores various information about the fvcom model,
%              the transect and also the output file names
%
%    Dependency: fvcom_medm2txt.f90             compiled by fortran if history files are in medm format
%                InPolygon.c                    compiled by mex in matlab
%                mat2tecplot.m                  code to convert matlab data to tecplot
%                mexcdf                          matlab toolbox to load netcdf file if history files are in netcdf
%
%    Platform: Unix/Linux, Mac OS
%
%    Example 1) use netcdf files
%             his2flux.fvcom_mat_grid_file = './chn_island.mat';                      %matlab format of the model grid
%             his2flux.his_dir           = '/home/long075/TK_to_wen_CHN_HYD_FLXTST/CHN_HYD_FLXTST/model_output/netcdf';
%                                                                                     %model history output dir
%             his2flux.hfile_prefix      = 'chn';                                     %model CASENAME
%             his2flux.hfile_ext         = '.nc';                                     %model history output extension,
%                                                                                     %      '.nc' for netcdf
%                                                                                     %      '.dat' for medm
%             his2flux.his_Nz            = 10;                                        %number of vertical layers
%             his2flux.his_fn_start      = 1;                                         %starting history file number
%             his2flux.his_fn_end        = 10;                                        %ending history file number
%             his2flux.fvcom_medm2txt_exec = 'fvcom_medm2txt';                        %executable program for medm2txt
%             his2flux.tracername        = 'salinity';                                %tracer name :'S' or 'T' for medm
%                                                                                     %   'salinity' or 'temp' for netcdf
%             his2flux.transname         ='cross';                                    %name of the transect
%             his2flux.xt                =[-44363,-43079];  %x coord of the starting and ending points of the transect
%             his2flux.yt                =[-539,3379];      %y coord of the starting and ending points of the transect
%             his2flux.ds_transect       = 100.0;           %step size of discretizing the transect (m)
%             his2flux.t_int             = 1;               %time record interval in processing the history outputs
%                                                           %     1 for every record,2 for every other record,...
%             his2flux.plotoutdir        =['./' his2flux.transname '/'];  %output plt dirname
%             fvcom_his2flux_interp(his2flux);
%
%    Example 2) use medm outputs
%
%             his2flux.fvcom_mat_grid_file = './chn_island.mat';                        %matlab format of the model grid
%             his2flux.his_dir           = '/home/long075/TK_to_wen_CHN_HYD_FLXTST/CHN_HYD_FLXTST/model_output/medm';
%                                                                                     %model history output dir
%             his2flux.hfile_prefix      = 'chn_sim';                                 %model CASENAME (add _sim for medm)
%             his2flux.hfile_ext         = '.dat';                                    %model history output extension,
%                                                                                     %      '.nc' for netcdf
%                                                                                     %      '.dat' for medm
%             his2flux.his_Nz            = 10;                                        %number of vertical layers
%             his2flux.his_fn_start      = 1;                                         %starting history file number
%             his2flux.his_fn_end        = 10;                                        %ending history file number
%             his2flux.fvcom_medm2txt_exec = 'fvcom_medm2txt';                        %executable program for medm2txt
%             his2flux.tracername        = 'S';                                       %tracer name :'S' or 'T' for medm
%                                                                                     %   'salinity' or 'temp' for netcdf
%             his2flux.transname         ='cross';                                    %name of the transect
%             his2flux.xt                =[-44363,-43079];  %x coord of the starting and ending points of the transect
%             his2flux.yt                =[-539,3379];      %y coord of the starting and ending points of the transect
%             his2flux.ds_transect       = 100.0;           %step size of discretizing the transect (m)
%             his2flux.t_int             = 1;               %time record interval in processing the history outputs
%                                                           %     1 for every record,2 for every other record,...
%             his2flux.plotoutdir        =['./' his2flux.transname '/'];  %output plt dirname
%             fvcom_his2flux_interp(his2flux);
%
% Wen Long @ PNNL, May  2nd, 2012
%                  Nov 12th, 2012
%                  May  3rd, 2013
%                  May 24th, 2013

%
%obtain input parameters from his2flux structure
%

     if(nargin<1)
      error(['oops,missing his2flux argument, quitting ...']);
     end

     if(isempty(his2flux.fvcom_mat_grid_file))
       display(['oops, input argument is empty, quitting ...']);
     end

     fieldnames_arg=fieldnames(his2flux) ;

%----check if arguments are there
     have_grid        = ~isempty(find(strcmp(fieldnames_arg,'fvcom_mat_grid_file')==1));
     have_hisdir      = ~isempty(find(strcmp(fieldnames_arg,'his_dir')==1));
     have_hpfix       = ~isempty(find(strcmp(fieldnames_arg,'hfile_prefix')==1));
     have_hext        = ~isempty(find(strcmp(fieldnames_arg,'hfile_ext')==1));
     have_Nz          = ~isempty(find(strcmp(fieldnames_arg,'his_Nz')==1));
     have_fn_start    = ~isempty(find(strcmp(fieldnames_arg,'his_fn_start')==1));
     have_fn_end      = ~isempty(find(strcmp(fieldnames_arg,'his_fn_end')==1));
     have_medmexec    = ~isempty(find(strcmp(fieldnames_arg,'fvcom_medm2txt_exec')==1));
%     have_yesnotf     = ~isempty(find(strcmp(fieldnames_arg,'yesno_tracer_flux')==1));
     have_tracername  = ~isempty(find(strcmp(fieldnames_arg,'tracername')==1));
     have_transname   = ~isempty(find(strcmp(fieldnames_arg,'transname')==1));
     have_xt          = ~isempty(find(strcmp(fieldnames_arg,'xt')==1));
     have_yt          = ~isempty(find(strcmp(fieldnames_arg,'yt')==1));
     have_dst         = ~isempty(find(strcmp(fieldnames_arg,'ds_transect')==1));
     have_tint        = ~isempty(find(strcmp(fieldnames_arg,'t_int')==1));
     have_pltdir      = ~isempty(find(strcmp(fieldnames_arg,'plotoutdir')==1));

     if(~have_grid)
        error(['oops, fvcom_mat_grid_file not given']);
     end
     if(~have_hisdir)
        error(['oops, his_dir not given']);
     end
     if(~have_hpfix)
        error(['oops, hfile_prefix not given']);
     end
     if(~have_hext)
        error(['oops, hfile_ext not given']);
     end
     if(~have_Nz)
        error(['oops, his_Nz not given']);
     end
     if(~have_fn_start)
        error(['oops, his_fn_start not given']);
     end
     if(~have_fn_end)
        error(['oops, his_fn_end not given']);
     end
     if(~have_medmexec)
        error(['oops, fvcom_medm2txt_exec not given']);
     end
%     if(~have_yesnotf)
%        error(['oops, yesno_tracer_flux not given']);
%     end
     if(~have_tracername)
        error(['oops, tracername not given']);
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

     fvcom_mat_grid_file  =his2flux.fvcom_mat_grid_file;
     his_dir            =his2flux.his_dir;
     hfile_prefix       =his2flux.hfile_prefix;
     hfile_ext          =his2flux.hfile_ext;
     his_Nz             =his2flux.his_Nz;
     his_fn_start       =his2flux.his_fn_start;
     his_fn_end         =his2flux.his_fn_end;
     fvcom_medm2txt_exec=his2flux.fvcom_medm2txt_exec;
     tracername         =his2flux.tracername;
     transname          =his2flux.transname;
     xt                 =his2flux.xt;
     yt                 =his2flux.yt;
     ds_transect        =his2flux.ds_transect;
     t_int              =int32(his2flux.t_int);
     pltdirname         =his2flux.plotoutdir;



     i_par=1;             %choice of variable
     yesno_tracer_flux='y'; 
                          %'salinity' or 'temp' for netcdf files
     climits=[-1 1];      %min and max value for velocity in the color plots
     zMax = 2.0;          %maximum z coordinate (above MSL) (meter)
     his_fn_int  =1;      %intervals of file processing

     %temporary directory name 
     if(~exist(pltdirname,'dir'))
        eval(['!mkdir -p ' pltdirname]);
     else
        eval(['!rm -r -f ' pltdirname]);
        eval(['!mkdir -p ' pltdirname]);
     end

     if(exist(fvcom_mat_grid_file,'file'))  %if there is grid in matlab format, then load it dirrectly
               load(fvcom_mat_grid_file);   
               x_n=xyd_n(:,2);                  %x coordinate of nodes
               y_n=xyd_n(:,3);                  %y coordinate of nodes
               lon_n=lld_n(:,2);                %longitude of nodes
               lat_n=lld_n(:,3);                %latitude of nodes
               h_n=lld_n(:,4);                  %depth of nodes
     else                                 %otherwise, load separately from text files 
        display(['not able to find file named ' fvcom_mat_grid_file]);
	error(['oops, please provide matlab format of the FVCOM grid file, it can ' ...
	       'be created using sms2dm2mat.m if you don''t have it']);
     end

     NN=size(xyd_n,1) ;    %number of nodes
     NE=size(xy_e,1) ;     %number of elements

     %obtain history output file (netcdf or medm)


             if(~exist(his_dir,'dir'))
               error(['Oops, history output file dir:' his_dir ' does not exist!']);
             end

             if(his_fn_start>his_fn_end)
               error(['Oops ending file number cannot be smaller than starting file number']);
             end

             if(his_fn_int<1 || his_fn_int>his_fn_end-his_fn_start+1)
               error(['Oops wrong file interval, must be of range : 1<= interval <=' int2str(his_fn_end-his_fn_start+1)]); 
             end

             %A typical netcdf output of fvcom model has the following information:
               %netcdf psm_1290 {
               %dimensions:
               %        scalar = 1 ;
               %        node = 9013 ;     %Number of nodes
               %        nele = 13941 ;    %Number of elements
               %        siglay = 10 ;     %Number of sigma layers
               %        siglev = 11 ;     %Number of sigma levels
               %        three = 3 ;
               %        four = 4 ;
               %        obc = 11 ;        %Number of open boundary nodes
               %        obc2 = 11 ;
               %        time = UNLIMITED ; // (72 currently)  %Number of time steps
               %variables:
               %        int nprocs(scalar) ;
               %                nprocs:long_name = "number of processors" ;
               %        int partition(nele) ;
               %                partition:long_name = "partition" ;
               %        float Initial_Density(siglay, node) ;
               %                Initial_Density:long_name = "Initial Density" ;
               %        float x(node) ;
               %                x:long_name = "nodal x-coordinate" ;
               %                x:units = "meters" ;
               %        float y(node) ;
               %                y:long_name = "nodal y-coordinate" ;
               %                y:units = "meters" ;
               %        float lon(node) ;
               %                lon:long_name = "Longitude" ;
               %                lon:standard_name = "longitude" ;
               %                lon:units = "degrees_east" ;
               %        float lat(node) ;
               %                lat:long_name = "Latitude" ;
               %                lat:standard_name = "latitude" ;
               %                lat:units = "degrees_north" ;
               %        float siglay(siglay) ;
               %                siglay:long_name = "Sigma Layers" ;
               %                siglay:standard_name = "ocean_sigma_coordinate" ;
               %                siglay:positive = "up" ;
               %                siglay:valid_min = "-1" ;
               %                siglay:valid_max = "0" ;
               %                siglay:formula_terms = "siglay:siglay eta:zeta depth:depth" ;
               %        float siglay_shift(siglay) ;
               %                siglay_shift:long_name = "Shifted Sigma Layers" ;
               %        float siglev(siglev) ;
               %                siglev:long_name = "Sigma Levels" ;
               %                siglev:standard_name = "ocean_sigma_coordinate" ;
               %                siglev:positive = "up" ;
               %                siglev:valid_min = "-1" ;
               %                siglev:valid_max = "0" ;
               %                siglev:formula_terms = "siglev:siglev eta:zeta depth:depth" ;
               %        float h(node) ;
               %                h:long_name = "Bathymetry" ;
               %                h:units = "meters" ;
               %                h:positive = "down" ;
               %                h:standard_name = "depth" ;
               %                h:grid = "fvcom_grid" ;
               %        float nv(three, nele) ;                                    %element to node mapping
               %                nv:long_name = "nodes surrounding element" ;
               %        float a1u(four, nele) ;
               %                a1u:long_name = "a1u" ;
               %        float a2u(four, nele) ;
               %                a2u:long_name = "a2u" ;
               %        float aw0(three, nele) ;
               %                aw0:long_name = "aw0" ;
               %        float awx(three, nele) ;
               %                awx:long_name = "awx" ;
               %        float awy(three, nele) ;
               %                awy:long_name = "awy" ;
               %        float time(time) ;
               %                time:long_name = "Time" ;
               %                time:units = "seconds after 00:00:00" ;
               %                time:calendar = "none" ;
               %        int iint(time) ;
               %                iint:long_name = "internal mode iteration number" ;
               %        float u(time, siglay, nele) ;
               %                u:long_name = "Eastward Water Velocity" ;
               %                u:units = "meters s-1" ;
               %                u:grid = "fvcom_grid" ;
               %                u:type = "data" ;
               %        float v(time, siglay, nele) ;              
               %                v:long_name = "Northward Water Velocity" ;
               %                v:units = "meters s-1" ;
               %                v:grid = "fvcom_grid" ;
               %                v:type = "data" ;
               %        float ww(time, siglay, nele) ;
               %                ww:long_name = "Upward Water Velocity" ;
               %                ww:units = "meters s-1" ;
               %                ww:grid = "fvcom_grid" ;
               %                ww:type = "data" ;
               %        float wts(time, siglev, node) ;
               %                wts:long_name = "Upward Water at node" ;
               %                wts:units = "meters s-1" ;
               %                wts:grid = "fvcom_grid" ;
               %                wts:type = "data" ;
               %        float uard_obcn(time, obc) ;
               %                uard_obcn:long_name = "UARD at OBC" ;
               %                uard_obcn:units = "meters s-1" ;
               %                uard_obcn:grid = "fvcom_grid" ;
               %                uard_obcn:type = "data" ;
               %        float xflux_obc(time, siglay, obc2) ;
               %                xflux_obc:long_name = "Xflux at OBC" ;
               %                xflux_obc:units = "meters s-1" ;
               %                xflux_obc:grid = "fvcom_grid" ;
               %                xflux_obc:type = "data" ;
               %        float dtfa(time, node) ;
               %                dtfa:long_name = "FSH + Water Height" ;
               %                dtfa:units = "meters" ;
               %                dtfa:positive = "down" ;
               %                dtfa:grid = "fvcom_grid" ;
               %                dtfa:type = "data" ;
               %        float kh(time, siglev, node) ;
               %                kh:long_name = "Turbulent Eddy Diffusivity" ;
               %                kh:units = "meters2 s-1" ;
               %                kh:grid = "fvcom_grid" ;
               %                kh:type = "data" ;
               %        float zeta(time, node) ;                             %Surface elevation
               %                zeta:long_name = "Water Surface Elevation" ;
               %                zeta:units = "meters" ;
               %                zeta:positive = "up" ;
               %                zeta:standard_name = "sea_surface_elevation" ;
               %                zeta:type = "data" ;
               %        float salinity(time, siglay, node) ;                 %Salinity
               %                salinity:long_name = "salinity" ;
               %                salinity:standard_name = "sea_water_salinity" ;
               %                salinity:units = "1e-3" ;
               %                salinity:grid = "fvcom_grid" ;
               %                salinity:type = "data" ;
               %        float temp(time, siglay, node) ;                     %temperature
               %                temp:long_name = "temperature" ;
               %                temp:standard_name = "sea_water_temperature" ;
               %                temp:units = "degrees_C" ;
               %                temp:grid = "fvcom_grid" ;
               %                temp:type = "data" ;
               %

          %A typical medm file contains the following information :

               %     !! ELEMENT BASED VALUES
               %     WRITE(1) IINTT,NGL,MGL,THOUR    !Internal time step, number of elements, 
                                                     %number of nodes and time in hours
               %
               %#    if !defined (TWO_D_MODEL)                                     !3D model 
               %        DO I=1,NGL
               %           WRITE(1) (U(I,K),V(I,K),WW(I,K),KM1(I,K),K=1,KBM1)      !  U, V, WW
               %         END DO
               %     !! NODE BASED VALUES
               %#       if !defined (DYE_RELEASE)
               %           DO I=1,M
               %              WRITE(1) EL(I),(T1(I,K),S1(I,K),RHO1(I,K),K=1,KBM1)  !  Elevation, temperature, 
               %                                                                   %salinity, density
               %           END DO
               %#       else
               %           DO I=1,M
               %              WRITE(1) EL(I),(T1(I,K),S1(I,K),RHO1(I,K),DYE(I,K),K=1,KBM1)
               %           END DO
               %#       endif
               %#    else                           !2D model
               %        DO I=1,NGL
               %           WRITE(1) UA(I),VA(I)     !    U, V
               %         END DO
               %        !! NODE BASED VALUES
               %        DO I=1,M
               %           WRITE(1) EL(I)           !    Elevation
               %        END DO
               %#    endif


 %         varnames={  'H'; 'zeta';  'D'; 'DL'; 'DO';         'LDOC';     'ALG1';   ... 
 %                      'ALG2';     'NH4';     'NO3';     'PO4';       'T';    'S'}; 
 %         varunits={'(m)';  '(m)';'(m)';'(m)';'(mgDO/L)';'(gC/m^3)'; '(gC/m^3)';'(gC/m^3)'; ...
 %                     '(gN/m^3)';'(gN/m^3)';'(gP/m^3)';'(\circC)';'(psu)'};


          %check if files exist
          read_his_nc=false;
          if(strcmp(hfile_ext,'.nc'))
              for nf=his_fn_start:his_fn_int:his_fn_end
                 hfile=[his_dir '/' hfile_prefix '_' num2str(nf,'%04d') hfile_ext];
                 if(~exist(hfile,'file'))
                   error(['Oops, cannot find file ' hfile ]);
                 end
              end
              read_his_nc=true;
          end

          read_his_medm=false; 
          if(strcmp(hfile_ext,'.dat'))
              for nf=his_fn_start:his_fn_int:his_fn_end
                 hfile=[his_dir '/' hfile_prefix  num2str(nf,'%04d') hfile_ext];
                 if(~exist(hfile,'file'))
                   error(['Oops, cannot find file ' hfile ]);
                 end
              end
              read_his_medm=true;
          end

          read_his_matlab=false;
          if(strcmp(hfile_ext,'.mat'))
              for nf=his_fn_start:his_fn_int:his_fn_end
                 hfile=[his_dir '/' hfile_prefix  num2str(nf,'%04d') hfile_ext];
                 if(~exist(hfile,'file'))
                   error(['Oops, cannot find file ' hfile ]);
                 end
              end
              read_his_matlab=true;
          end

         if(~(read_his_matlab || read_his_medm || read_his_nc))
            error(['Oops, wrong hisfile file extension']);
         end


         if(i_par==1)
            calculate_tracerflux=0;
            calculate_tracerflux=(lower(yesno_tracer_flux(1:1))=='y'); 
            if(calculate_tracerflux)
               if(read_his_medm)
			if (tracername ~= 'S' || tracername ~= 'T')
                            error(['tracer name must be ''S'' or ''T'' for medm history outputs']);
                        end
               end
               if(read_his_nc)
                        if (tracername ~= 'salinity' || tracername ~= 'temp')
                            error(['tracer name must be ''salinity'' or ''temp'' for medm history outputs']);
                        end
               end
            end
         end

           if(climits(2)<=climits(1))
              error('input of color value is incorrect, maximum can not be less than minimum');
           end

           zlimit(1)=inf;
           zlimit(2)=-inf;
           zlimit(2)=zMax; 

           %plot a grid and let user choose transect

           yesno_matlabfig=0;  %no matlab plots to be made

          %get information of the transect 

         if (isempty(xt)||isempty(yt) )
            error(['xt and yt should not be empty']);
         end
	 s = [0 cumsum(sqrt(diff(xt).^2+diff(yt).^2))' ]';

	 if(ds_transect >=0.5*max(s(:)))
	   warning(['the segment length you input is too large']);
	   warning(['setting to default value of ' num2str(max(s(:))) '/100']);
	   ds_transect= max(s(:))/100.0; 
	 end

	 dx=ds_transect;
	 st=[s(1):dx:max(s) max(s)];
	 xp=spline(s,xt,st);
	 yp=spline(s,yt,st);
	 save(['trans_loc_' transname '.mat'],'xp','yp','xt','yt');

	 Ntrans=length(st);  
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

       %save the map and the transect figure
       saveas(gcf,[pltdirname './' transname '_transect.eps'], 'epsc2');
       saveas(gcf,[pltdirname './' transname '_transect.png'], 'png');
       saveas(gcf,[pltdirname './' transname '_transect.fig'], 'fig');

       delete(hmap);
    end

    if(isempty(Ipath_in))
	error(['Oops all transect points are out of model domain, quiting ...']);
	exit(1);
    end

    %
    %load history data
    %
    display(['Please wait while loading history output files']);


    %Loop through all files one by one and make the transect plots
    %save data into matlab file for statistics after the time loop is finished

      tdata2d=[];   %2Dimensional tecplot in (r,z) coordinates using ordered surface
      tdata3d=[];   %3Dimensional tecplot in (x,y,z) coordinates using ordered surface
      tdata2d_timeavg=[];  %time average
      tdata3d_timeavg=[];  %time average

      %run mean (to finish coding)
      %tdata2d_runmean=[];
      %tdata3d_runmean=[];

      %get variable information

      switch(i_par)
            case  1  %three components of u,v,w
                tdata2d.Nvar=8;
                tdata2d.varnames={'r','z','zdummy','u','v','w','u_p','v_p'};

                tdata3d.Nvar=9;
                tdata3d.varnames={'x','y','z','u','v','w','u_p','v_p','depth'};

                tdata2d_timeavg.Nvar=8;
                tdata2d_timeavg.varnames={'r','z','zdummy','u','v','w','u_p','v_p'};

                tdata3d_timeavg.Nvar=9;
                tdata3d_timeavg.varnames={'x','y','z','u','v','w','u_p','v_p','depth'};
     end

    %
    %load history data
    %

    %
    %Loop through all files one by one
    %save data into matlab file for statistics after the time loop is finished
    %

    display(['Please wait while loading history output files']);

    it_all=0;  %time record counter
    htime_all=[];

    for nf=his_fn_start:his_fn_int:his_fn_end

               %read history output file
            if(read_his_nc)            %open netcdf file
               hfile=[his_dir '/' hfile_prefix '_' num2str(nf,'%04d') hfile_ext];
               hhis=netcdf(hfile);   %open the history file using mexcdf/mexnc/netcdf_toolbox

              %read all time history of surface elevation (this is memory intensive)
              %read all time history of chosen variable

              htime=hhis{'time'}(:)./86400;    %read time (sec) and convert to days
              zeta=hhis{'zeta'}(:);     %read surface elevation ((m)
              h=hhis{'h'}(:); 
              par2=[];
              par3=[];
              par4=[]; 
              
              switch(i_par)
                    case  1  %three components of u,v,w
                         par1=hhis{'u'}(:);
                         par2=hhis{'v'}(:);
                         par3=hhis{'ww'}(:);

                         parname='velocity';
                         parunit='(m/s)';
                         par4=hhis{tracername}(:); 
                         if(prod(size(par4))==0)
                           error(['oops, not able to read variable ' ...
                                  tracername ' in file ' hfile, ' or it is empty']); 
                         end
                         %get tracer units here
                         %tracerunit=

                    case  2  %temperature
                         par1=hhis{'temp'}(:);
                         parname='temperature';
                         parunit='(\circC)';

                    case  3  %salt
                         par1=hhis{'salinity'}(:);
                         parname='salinity';
                         parunit='(psu)';
                    case  4  %turbulent diffusivity
                         par1=hhis{'kh'}(:);
                         parname='kh';
                         parunit='(m^2/s)'
		    case 5
		 	par1=hhis{myvarname}(:);
			parname=myvarname;
                        if(~isempty(hhis{myvarname}.units(:)))          
   			   parunit=['(' hhis{myvarname}.units(:) ')'];  %get the units of the variable
                        else 
			   parunit='(unknown unit)';
			end
              end
              close(hhis);   %close netcdf file

              %save data into matlab format for later use

               tmp_matfile=[pltdirname  'fvcom_hisdata_' num2str(nf,'%04d') '.mat'];
               save(tmp_matfile,   'par1', ...  %1st variable
                                   'par2', ...  %2nd variable
                                   'par3', ...  %3rd variable
                                   'par4', ...  %4th variable 
                                'parname', ...  %variable name
                                'parunit', ...  %unit of variable
                                   'zeta', ...  %surface elevation
                                    'x_n', ...  %node x coordinate
                                    'y_n', ...  %node y coordinate
                                  'lon_n', ...  %node longitude
                                  'lat_n', ...  %node latitude
                                    'h_n', ...  %depth (still water) (read from grid file)
                                      'h', ...  %depth read from netcdf file
                                   'xy_e', ...  %x,y coordinates of elements
                                    'e2n', ...  %connectivity matrix (element to node)
                                  'htime', ...  %time in days
                    '-mat');

                NT=length(htime);   %number of time records in this file
                display(['successfully saved ' num2str(NT) ' records of data into file: ' tmp_matfile]);

            end

            if(read_his_medm)

               [status, find_fvcom_medm2txt_exec]=system(['!which ' fvcom_medm2txt_exec ' 2>&1 ' '|awk ''{print $2}''']);

               if(strcmp(lower(find_fvcom_medm2txt_exec(1:end-1)),'no'))
                 error(['oops, can''t find executable ' fvcom_medm2txt_exec ' in your system']);
               end


               hfile=[his_dir '/' hfile_prefix  num2str(nf,'%04d') hfile_ext];

               if(i_par==6)  %if sediment concentration from medm outputs
                   hfile_sed=[his_dir '/' hfile_prefix_sed  num2str(nf,'%04d') hfile_ext];
               end

               %read time in days into htime 
               tmp_file_name='thour_tmp_tmp.txt';
               command_line=[fvcom_medm2txt_exec ' ' hfile ' ' tmp_file_name '  thour ' num2str(his_Nz+1)];
               eval(['!' command_line]);
               if(exist(tmp_file_name,'file'))
                    htime=load(tmp_file_name);
                    htime=htime/24.0 ;              %convert from hour to day
                    eval(['!rm -f ' tmp_file_name]);;
               else
                    warning(['Oops, not able to read ' tmp_file_name  ',please check the following command:']);
                    error(['This command failed:' command_line ]);
               end

               %read surface elevation
               command_line=[fvcom_medm2txt_exec ' ' hfile ' zeta_tmp_tmp.txt el ' num2str(his_Nz+1)];
               eval(['!' command_line]);  %call fortran program fvcom_medm2txt to dump the variable el in a txt
               if(exist('zeta_tmp_tmp.txt','file'))  %pick up the text file
                  zeta=load('zeta_tmp_tmp.txt');
                  zeta=zeta';             %transpose to have first dimension being depth
                  eval(['!rm -r -f zeta_tmp_tmp.txt']);  %remove temp file
               else
                  warning(['Oops, not able to read zeta_tmp_tmp.txt,', 'please check the following command:']);
                  error(['This command failed:' command_line ]);
               end

               h=h_n; 
               NT=length(htime);
               par2=[];
               par3=[];
               par4=[]; 

               switch(i_par)
                    case  1  %three components of u,v,w

                         %read u into par1
                         tmp_file_name='u_tmp_tmp.txt';
                         command_line=[fvcom_medm2txt_exec ' ' hfile ' ' tmp_file_name '  u ' num2str(his_Nz+1)];
                         eval(['!' command_line]);  %call fortran program fvcom_medm2txt to dump the variable u in a txt
                         if(exist(tmp_file_name,'file'))
                            par1=load(tmp_file_name);
                            par1=par1';
                            par1=reshape(par1,[1,size(par1)]);  %[1,his_Nz,NGL]
                            eval(['!rm -f ' tmp_file_name]);;
                         else
                            warning(['Oops, not able to read ' tmp_file_name ', please check the following command:'])
                            error(['This command failed:' command_line ]);
                         end
  
                         %read v into par2
                         tmp_file_name='v_tmp_tmp.txt';
                         command_line=[fvcom_medm2txt_exec ' ' hfile ' ' tmp_file_name '  v ' num2str(his_Nz+1)];
                         eval(['!' command_line]);  %call fortran program fvcom_medm2txt to dump the variable v in a txt
                         if(exist(tmp_file_name,'file'))
                            par2=load(tmp_file_name);
                            par2=par2';
                            par2=reshape(par2,[1,size(par2)]);
                            eval(['!rm -f ' tmp_file_name]);;
                         else
                            warning(['Oops, not able to read ' tmp_file_name ', please check the following command:']);
                            error(['This command failed:' command_line ]);
                         end

                         %read w into par3
                         tmp_file_name='w_tmp_tmp.txt';
                         command_line=[fvcom_medm2txt_exec ' ' hfile ' ' tmp_file_name '  w ' num2str(his_Nz+1)];
                         eval(['!' command_line]);  %call fortran program fvcom_medm2txt to dump the variable w in a txt
                         if(exist(tmp_file_name,'file'))
                            par3=load(tmp_file_name);
                            par3=par3';
                            par3=reshape(par3,[1,size(par3)]);
                            eval(['!rm -f ' tmp_file_name]);;
                         else
                            warning(['Oops, not able to read ' tmp_file_name ', please check the following command:']);
                            error(['This command failed:' command_line ]);
                         end

                         parname='velocity';
                         parunit='(m/s)';
 
                         %read tracer into par4 if need to calculate flux through transect
                         tmp_file_name='par_tmp_tmp.txt';
                         command_line=[fvcom_medm2txt_exec ' ' hfile ' ' tmp_file_name ' ' ...
                                       lower(tracername) ' ' num2str(his_Nz+1)];
                         eval(['!' command_line]);  %call fortran program fvcom_medm2txt to dump 
                                                    % the variable tracername in a txt
	                         if(exist(tmp_file_name,'file'))
        	                    par4=load(tmp_file_name);
	                            par4=par4';
	                            par4=reshape(par4,[1,size(par4)]);
	                            eval(['!rm -f ' tmp_file_name]);;
	                         else
	                            warning(['Oops, not able to read ' tmp_file_name  ...
                                             ', please check the following command:']);
	                            error(['This command failed:' command_line ]);
	                         end

                    case  2  %temperature

                         %read t into par1
                         tmp_file_name='t_tmp_tmp.txt';
                         command_line=[fvcom_medm2txt_exec ' ' hfile ' ' tmp_file_name '  t ' num2str(his_Nz+1)];
                         eval(['!' command_line]);  %call fortran program fvcom_medm2txt 
                                                    %to dump the variable el in a txt
                         if(exist(tmp_file_name,'file'))
                            par1=load(tmp_file_name);
                            par1=par1';
                            par1=reshape(par1,[1,size(par1)]);
                            eval(['!rm -f ' tmp_file_name]);;
                         else
                            warning(['Oops, not able to read ' tmp_file_name ', please check the following command:']);
                            error(['This command failed:' command_line ]);
                         end
                         parname='temperature';
                         parunit='(\circC)';

                    case  3  %salt

                         %read s into par1
                         tmp_file_name='s_tmp_tmp.txt';
                         command_line=[fvcom_medm2txt_exec ' ' hfile ' ' tmp_file_name '  s ' num2str(his_Nz+1)];
                         eval(['!' command_line]);  %call fortran program fvcom_medm2txt to 
                                                    %dump the variable el in a txt
                         if(exist(tmp_file_name,'file'))
                            par1=load(tmp_file_name);
                            par1=par1';
                            par1=reshape(par1,[1,size(par1)]);
                            eval(['!rm -f ' tmp_file_name]);;
                         else
                            warning(['Oops, not able to read ' tmp_file_name ', please check the following command:']);
                            error(['This command failed:' command_line ]);
                         end
                         parname='salinity';
                         parunit='(psu)';

                    case 6  %sediment concentration
                            %read sediment concentration into par1
 
                         tmp_file_name='sed_tmp_tmp.txt';
                         command_line=['fvcom_medmsed2txt ' hfile_sed ' ' tmp_file_name ' sed ' ...
                                        num2str(his_Nz+1) ' ' ...
                                        num2str(nsed) ' '     ...
                                        num2str(nbed) ' '     ...
                                        num2str(ised) ];

                         eval(['!' command_line]);  %call fortran program fvcom_medm2txt to 
                                                    %dump the variable el in a txt
                         if(exist(tmp_file_name,'file'))
                            par1=load(tmp_file_name);
                            par1=par1';
                            par1=reshape(par1,[1,size(par1)]);
                            eval(['!rm -f ' tmp_file_name]);;
                         else
                            warning(['Oops, not able to read ' tmp_file_name  ...
                                     ', please check the following command:']);
                            error(['This command failed:' command_line ]);
                         end
                         parname='sediment';
                         parunit='(kg/m^3)';
              end

              %save data into matlab for later use
               tmp_matfile=[pltdirname 'fvcom_hisdata_' num2str(nf,'%04d') '.mat'];
               save(tmp_matfile,   'par1', ...  %1st variable
                                   'par2', ...  %2nd variable
                                   'par3', ...  %3rd variable
                                   'par4', ...  %4th variable
                                'parname', ...  %variable name
                                'parunit', ...  %unit of variable
                                   'zeta', ...  %surface elevation
                                    'x_n', ...  %node x coordinate
                                    'y_n', ...  %node y coordinate
                                  'lon_n', ...  %node longitude
                                  'lat_n', ...  %node latitude
                                    'h_n', ...  %depth (still water) read from grid
                                    'h'  , ...  %depth 
                                   'xy_e', ...  %x,y coordinates of elements
                                    'e2n', ...  %connectivity matrix (element to node)
                                  'htime', ...  %time in days
                    '-mat');
                display(['successfully saved ' num2str(NT) ' records of data into file: ' tmp_matfile]);
            end  %end of reading medm file

   end %end of nf loop

   %
   %Now only read matlab format data 
   %

   %find total number of time records by going through all the files
   %also record the file numbers for for each time record

   NT_all=0;
   htime_global=[];
   filenum_global=[];
   it_local_global=[];

   for nf=his_fn_start:his_fn_int:his_fn_end
       if(read_his_matlab==true) 
         hfile=[his_dir '/' hfile_prefix  num2str(nf,'%04d') hfile_ext];
       else
         hfile=[pltdirname ',/' 'fvcom_hisdata_' num2str(nf,'%04d') '.mat'];
       end
       if(~exist(hfile,'file'))
         error(['Oops, cannot find file ' hfile ]);
       end
       load(hfile);
       his_NT=length(htime);   %number of time records in this file
                               %give references to find file numbers and record numbers
                               %in file for any given record in htime_all

       htime_global=[htime_global htime(:)'];
       filenum_global=[filenum_global nf*ones([1,his_NT])];    %file number for each time record in htime_all
       it_local_global=[it_local_global (1:his_NT)];           %local record number in each file
   end

   NT_global=length(htime_global);
   NT_all=length(htime_global(1:t_int:NT_global);

   %pre-allocate htime_all 
   htime_all=zeros([NT_all,1]); 

   %Loop through all files and then do interpolation onto transect
   it_all=0; 

   for  it_global=1:t_int:NT_global

        it_all=it_all+1;

        %find the history file to load
        nf=filenum_global(it_global);

        if(it_all==1)
              hfile_old=[];
        else
              hfile_old=hfile;
        end

        %Do interpolation onto transect
        %read matlab file of the data saved above
        if(read_his_matlab)
            hfile=[his_dir '/' hfile_prefix  num2str(nf,'%04d') hfile_ext];
        else
            hfile=[pltdirname 'fvcom_hisdata_' num2str(nf,'%04d') '.mat'];
        end
        if(~exist(hfile,'file'))
            error(['not able to find file :' hisfile]); 
        end


        %find the local time record (record number in individual history file) to load

        it_local=it_local_global(it_global);
        it=it_local;

 
            if(~strcmp(hfile,hfile_old)) %if file is new, read new file
                load(hfile);
                NT=length(htime);   %number of time records in this file
                display(['successfully read ' num2str(NT) ' records from file: ' hfile]);

                if(it_all==1)
                        his_NT=NT;  %number of time records in this file
                        his_Nz=size(par1,2);   % number of vertical layers (siglay for u,v,ww,T,S)  (siglev for kh)
                        if(strcmp(parname,'kh'))   %siglev
                           KB=his_Nz;
                        else
                           KB=his_Nz+1;           %siglev=siglay+1
                        end
                        KBM1=KB-1;
                        z=nan.*ones(size(KB,1));
                        for k=1:KB                          %Note that this is hardwired!
                            z(k)= -((k-1.0)/(KB-1.0))^1.5 ; %vertical distortion (sigma levels)
                        end
                end
                display(['using file ' hfile ' for global time step ' num2str(it_global)]);
            else
                display(['using file ' hfile ' for global time step ' num2str(it_global)]);
                %do nothing
            end

            htime_all(it_all)=htime(it);   %store the time


            %interpolate the data onto the transect horizontally layer by layer

            %Calculate the along transect distance as axis along the transect curve in r
            %latlon2meter=110.0;   %should calculate this based on lat,lon at center of grid
            dxp=diff(xp);
            dyp=diff(yp);
            %dr=latlon2meter.*sqrt(dxp.*dxp+dyp.*dyp);  %dr in km 

            dr=sqrt(dxp.*dxp+dyp.*dyp);   %dr in meter
            dr=dr/1000.0;                 %conver to km unit
            r=nan.*ones(his_Nz+1,length(xp));
            for k=1:his_Nz+1
                r(k,1)=0;
                r(k,2:end)=cumsum(dr);    
            end

            if(nf==his_fn_start)
               
               %find min and max vertical cordinate along the transect so that vertical coord scale is fixed
               zlimit(1)=inf;   %min z coord
               hpath=griddata(x_n,y_n,h_n,xp,yp);                  
               for it=1:max(t_int,floor(his_NT/10.0)):his_NT  %in order to save time, only take 10 time records or less
                 zlimit(1)=min(zlimit(1),  min(-hpath(Ipath_in)));
               end 

            end

            Tp1=nan.*ones(KB,length(xp));   %value of parameter
            if(strcmp(parname,'velocity'))
                Tp2=nan.*ones(KB,length(xp));   %2nd component of velocity (v) if parameter is velocity (Northward)
                Tp3=nan.*ones(KB,length(xp));   %3rd component of velocity (w) if parameter is velocity (Upward)
                if(calculate_tracerflux)
                  Tp4=nan.*ones(KB,length(xp)); %tracer to calculate tracer flux through the transect
                end
                u_p=nan.*ones(KB,length(xp));   %velocity component along the transect direction
                v_p=nan.*ones(KB,length(xp));   %velocity component normal (through) the transect
            end

            Tz=nan.*ones(KB,length(xp));   %vertical coordinate
            Tz_tec=zeros(size(Tz));        %nans are good for matlab but not for tecplot
            z=nan.*ones(size(KB,1));
            
            for k=1:KB                     %Note that this is hardwired!  
                z(k)= -((k-1.0)/(KB-1.0))^1.5 ;  %vertical distortion (sigma levels)
            end

            %plot

            if(yesno_matlabfig)
               hfig=figure('visible','off','position',[1,1,1000,800]); hold on;
               set(gcf,'renderer','zbuffer');
            end

               %interpolate the data horizontally layer by layer onto the transect path (xp,yp)
            for k=1:his_Nz  %for kh, his_NZ is KB, for T, S, u,v,w, his_NZ is KB-1

                   %find vertical coordinate (layer or level)

                   if(his_Nz==KB-1)  %sigma layer coordinate (layer is demarcated by two levels)
                       zz(k)=(z(k)+z(k+1))/2.0;   %sigma layer (at center of each layer)
                       zcoord_layer=zeta(it,:)+(zeta(it,:)+h').*zz(k); 
                       Tz(k,:)=griddata(x_n(IN_use_uniq),y_n(IN_use_uniq), ...
                                       zcoord_layer(IN_use_uniq),xp,yp); %z coordinate of k'th layer (center of layer)

                       Tz_tec(k,:)=Tz(k,:);
                       isTznan=isnan(Tz_tec(k,:));  %find locations in Tz_tec that are nans
                                                    %for this locations reinterpolate with nearest value instead 
                                                    %to avoid nans in Tz_tec
                       if(length(isTznan)>=1)
                         Tz_tec(k,isTznan)=griddata(x_n(IN_use_uniq),y_n(IN_use_uniq),...
                                                    zcoord_layer(IN_use_uniq),xp(isTznan),yp(isTznan),'nearest'); 
                       end

                   else                             %sigma level coordinate
                       zcoord_lev=zeta(it,:)+(zeta(it,:)+h').*z(k); 
                       Tz(k,:)=griddata(x_n(IN_use_uniq),y_n(IN_use_uniq), ...
                                        zcoord_lev(IN_use_uniq),xp,yp); %z coordinate of k'th level

                       Tz_tec(k,:)=Tz(k,:);
                       isTznan=isnan(Tz_tec(k,:));  %find locations in Tz_tec that are nans
                                                    %for this locations reinterpolate with nearest value instead
                                                    %to avoid nans in Tz_tec
                       if(length(isTznan)>=1)
                         Tz_tec(k,isTznan)=griddata(x_n(IN_use_uniq),y_n(IN_use_uniq), ...
                                                   zcoord_lev(IN_use_uniq),xp(isTznan),yp(isTznan),'nearest');
                       end

                   end

                   %interpolate parameter onto path (linear interpolation,and nearest value used for extrapolaton) 
                   %         ideally we should use least square method to
                   %         regression into linear function of v(x,y)=ax+by+c with
                   %         given values at element center xy_e surround a node location (x,y)
                   %         If there are 3 elements only surrounding a node, then this solves exactly a,b,c
                   %         If there are more than 3 elements, then it is over specified, and we use least
                   %         square to find a,b,c that minimizes the error on all element centers xy_e
                   %

                   if(strcmp(parname,'velocity'))


                       Tp1(k,:)=griddata(xy_e(IE_use_uniq,1),xy_e(IE_use_uniq,2),par1(it,k,IE_use_uniq),xp,yp);
                       %for velocity, also need to interpolate par2 and par3
                       Tp2(k,:)=griddata(xy_e(IE_use_uniq,1),xy_e(IE_use_uniq,2),par2(it,k,IE_use_uniq),xp,yp);   %v
                       Tp3(k,:)=griddata(xy_e(IE_use_uniq,1),xy_e(IE_use_uniq,2),par3(it,k,IE_use_uniq),xp,yp);   %w


                       isTpnan=isnan(Tp1(k,:));  %find locations in Tp1 that are nans 
                                                 %(which means locations outside of the IE_use_uniq
                                                 %elements, and linear interpolation gives nans)
                                                 %for these nan places, use nearest value as 
                                                 %extrapolation method. This happens when 
                                                 %interpolated locations are within
                                                 %boundary but out of the convex hull of 
                                                 %set of points made of all chosen element centers,
                                                 %because element centers are
                                                 %away (interior) from the boundaries by nodes. 
                                                 %E.g., to obtain velocity on boundary nodes 
                                                 %(positive depth), one needs to 
                       %do extrapolation. 

                       if(length(isTpnan)>=1)
                         Tp1(k,isTpnan)=griddata(xy_e(IE_use_uniq,1),xy_e(IE_use_uniq,2), ...
                                       par1(it,k,IE_use_uniq),xp(isTpnan),yp(isTpnan),'nearest');
                         Tp2(k,isTpnan)=griddata(xy_e(IE_use_uniq,1),xy_e(IE_use_uniq,2), ...
                                       par2(it,k,IE_use_uniq),xp(isTpnan),yp(isTpnan),'nearest');
                         Tp3(k,isTpnan)=griddata(xy_e(IE_use_uniq,1),xy_e(IE_use_uniq,2), ...
                                       par3(it,k,IE_use_uniq),xp(isTpnan),yp(isTpnan),'nearest');
                       end

                       if(length(Ipath_out)>0)  %for places out of domain, still set to nans 
                                                %so that matlab will not plot them
                          Tp1(k,Ipath_out)=nan; %and tecplot can set it to blank value of -999999999.0
                          Tp2(k,Ipath_out)=nan;
                          Tp3(k,Ipath_out)=nan;
                       end

                      %also rotate velocity (u,v) = (TP1, TP2) to along and perpendicular to the path

                      u_p(k,:)=Tp1(k,:).*tx_p+Tp2(k,:).*ty_p;    %project onto tangential direction (inner product)
                      v_p(k,:)=Tp1(k,:).*nx_p+Tp2(k,:).*ny_p;    %project onto normal direction

                      %if calculate tracer flux through the transect, then get ther tracer concentration, flux will
                      %then be caluclated based on v_p*concentration integrated over the transect
                      Tp4(k,:)=griddata(x_n(IN_use_uniq),y_n(IN_use_uniq),par4(it,k,IN_use_uniq),xp,yp); 
                                                %tracer concentration
                      isTpnan=isnan(Tp4(k,:));
                      if(length(isTpnan)>=1)    %fix interpolation if interpolated to be nan (need extrapolcation)
                          Tp4(k,isTpnan)=griddata(x_n(IN_use_uniq),y_n(IN_use_uniq),par4(it,k,IN_use_uniq), ...
                                                  xp(isTpnan),yp(isTpnan),'nearest');
                      end
                      if(length(Ipath_out)>0)   %set to nan for out of domain section of the transect
                          Tp4(k,Ipath_out)=nan;
                      end
                   else
                       Tp1(k,:)=griddata(x_n(IN_use_uniq),y_n(IN_use_uniq),par1(it,k,IN_use_uniq),xp,yp); 
                        %u, T, S, kh
                       if(length(Ipath_out)>0)
                          Tp1(k,Ipath_out)=nan;
                       end
                   end
            end

            %bed level values are created so that movie pictures are not jumping up and down 
            %because of thickness change over time

            if(KB==his_Nz+1)  %if data is on sigma layer, append bottom value as bottom level (z=-h)
                   k=KB; 
                   zcoord_bottom= zeta(it,:)+(zeta(it,:)+h').*z(k);
                   Tz(k,:)=griddata(x_n(IN_use_uniq),y_n(IN_use_uniq), zcoord_bottom(IN_use_uniq),xp,yp); 
                                                %z coordinate of bottom
                   Tz_tec(k,:)=Tz(k,:);
                   isTznan=isnan(Tz_tec(k,:));  %find locations in Tz_tec that are nans
                                                %for this locations reinterpolate with nearest value instead
                                                %to avoid nans in Tz_tec
                   if(length(isTznan)>=1)
                      Tz_tec(k,isTznan)=griddata(x_n(IN_use_uniq),y_n(IN_use_uniq), ... 
                        zcoord_bottom(IN_use_uniq),xp(isTznan),yp(isTznan),'nearest');
                   end
                   Tp1(his_Nz+1,:)=Tp1(his_Nz,:);        %repeat bottom level's data
                   if(strcmp(parname,'velocity'))
                       Tp2(his_Nz+1,:)=Tp2(his_Nz,:);
                       Tp3(his_Nz+1,:)=Tp3(his_Nz,:);    %this is vertical velocity on bed level,
                                                         %maybe we should set it to zero for flat bottom
                                                         %however, when there is a slope, then 
                                                         %vertical velocity can be non-zero, as 
                                                         %the only physical argument is that normal 
                                                         %velocity should be zero for slip boundary condition
                       u_p(his_Nz+1,:)=u_p(his_Nz,:);    %velocity along the transect (horizontally)
                       v_p(his_Nz+1,:)=v_p(his_Nz,:);    %velocity perpendicular to the transect (horizontally)
                       Tp4(his_Nz+1,:)=Tp4(his_Nz,:);    %tracer concentration (default is salinity)
                   end
            end

            %calculate flux of chosen velocity or tracer
            if(strcmp(parname,'velocity'))
                   %find cross sectional area by taking difference in z and then in r dierction
                   Tz_new=nan.*ones(KB+1,length(xp));
                   %Tz_new is vertical coordinate of sigma level
                   zcoord_lev=zeta(it,:)+(zeta(it,:)+h').*z(1);       %add surface level (surface elevation)
                   Tz_new(1,:)=griddata(x_n(IN_use_uniq),y_n(IN_use_uniq), ...
                               zcoord_lev(IN_use_uniq) ,xp,yp); %z coordinate of k'th level, k=1
                   Tz_new(2:end-1,:)=Tz(1:end-1,:)+diff(Tz,1,1)/2.0;  %off set dowward by half grid, k=2 to KB
                   Tz_new(end,:)=Tz(end,:);       %make the last two levels the same (bottom level), k=KB+1
                   Tz_new(end-1,:)=Tz(end,:);     %
                   r_new=nan*ones([size(r,1),size(r,2)+1]);          %r=ones(his_Nz+1,length(xp)), 
                                                                     %r_new has dimension his_Nz+1 by length(xp)+1
                   r_new(:,1)=r(:,1)+(r(:,1)-r(:,2))/2.0;            %off set to left by half grid
                   r_new(:,2:end-1)=r(:,1:end-1)+diff(r,1,2)/2.0;    %off set right by half grid
                   r_new(:,end)=r(:,end)+(r(:,end)-r(:,end-1))/2.0;  %off set right by half grid

                   dz_trans=-diff(Tz_new,1,1);     %first derivative in z (row) direction. %last row of dz_trans is zero
                                                   %size of dz_trans is same as v_p  
                   dr_trans=diff(r_new,1,2)*1000;  %first derivatinve in along transect direction (column) 
                                                   % (*1000 to convert from km to m)
                                                   %last two columns of dr_trans is incorrect
                                                   %size of dr_trans is same as v_p

                   dr_trans(:,end-2)=dr_trans(:,end-3);  %repeat the last two columns
                   dr_trans(:,end)=dr_trans(:,end-3);
                   dA_trans=dz_trans(1:end-1,:).*dr_trans(1:end-1,:);  %area calculation (m^2) size KB-1 by length(xp) 
                   if(calculate_tracerflux)  %flux of tracer
                      flux_trans=dA_trans.*v_p(1:end-1,:).*Tp4(1:end-1,:);  %unit: (m^3/sec) * 
                                             %(tracer volume concentration (tracer/m^3))
                   else 
                      flux_trans=dA_trans.*v_p(1:end-1,:);                 %velocity volume  flux through each window
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

				rho_bulk_S    = 1028 * (1 - 1.7E-4 * (10 -10 ) + 7.6E-4 *(Tp4-35));
                                          %eqn 7 assuming S=Tp4  (mixed)
                                rho_bulk_Smax = 1028 * (1 - 1.7E-4 * (10 -10 ) + 7.6E-4 *(34 -35)); 
                                          %eqn 7 assuming S=Smax=34 ppt (salt)
                                rho_fw        = 1028 * (1 - 1.7E-4 * (10 -10 ) + 7.6E-4 *(0  -35));
                                          %eqn 7 assuming S= 0 ppt  (fresh)
                          %get f_fw using  eqn(6) above
                                f_fw= ( rho_bulk_Smax - rho_bulk_S )./(rho_bulk_Smax - rho_fw) ;
                          %get fresh water flux (m^3/sec)
				flux_trans_fw=flux_trans.*f_fw(1:end-1,:);    %eqn (6-1) 
                   end

                   Igood=~isnan(flux_trans(:));     %find places that are not nans, when part of transect
                                                    % is outside of model domain
                   flux_total=sum(flux_trans(Igood));   %total flux is summation of all 
                                                        %good grids on transect (m^3/sec)
                   %initialize flux_all if it_all =1
                   if(it_all==1)
                     flux_all=zeros([NT_all,1]); 
                   end
                   flux_all(it_all)=flux_total;                         %overall flux
                   if(~calculate_tracerflux)
                           if(it_all==1)  %initialize flux_all_fw if it_all==1
                              flux_all_fw=zeros([NT_all,1]); 
                           end
                           flux_total_fw=sum(flux_trans_fw(Igood));     %total fresh water flux (m^3/s)
                           flux_all_fw(it_all)=flux_total_fw;           %overall fresh water flux (m^3/s)
                   end

                   %
                   %calculate vertical distribution of the flux 
                   %Because cross section has variable depth, we can not use sigma coordinates to calculate
                   %vertical distribution of flux (horzontally integrated normal velocity) 
                   %we have to switch to z-coorinate griding vertically
                   %

                   Tz_min=  zlimit(1)  ;   %min(-hpath(Ipath_in));         
                                           %minimum vertical coordinate on the transect
                   Tz_max=  zlimit(2)  ;   %max(zcoord_lev(Ipath_in));     
                                           %maximum vertical coordinate on the transect (surface elevation)
                   
                   dz_vdis=(Tz_max-Tz_min)/(his_Nz);     %divide the vertical distance by same
                                                         %number of layers as his_Nz
                   Tz_vert=(Tz_min:dz_vdis:Tz_max);      %evenly distributed vertical coordinate 
                   Nz_vert=length(Tz_vert);
                   Tz_vert2D=repmat(Tz_vert',[1,size(r,2)]);
                   r_vert2D =repmat(r(1,:),[Nz_vert,1]);

                   Tz_tmp=Tz(1:end-1,:);                 %same size as dA_trans flux_trans
                   r_tmp =r(1:end-1,:); 
                   flux_vert=zeros([1,Nz_vert]); 

                   if(~calculate_tracerflux)
                       flux_vert_fw=zeros([1,Nz_vert]);
                   end

                   %collect cross sectional flux on each interval of Tz_vert
                   for k=1:Nz_vert-1
                       %find locations in transect that are between
                       %Tz_vert(k) and Tz_vert(K+1)
                       Ifind=find(Tz_tmp>Tz_vert(k) & Tz_tmp <=Tz_vert(k+1) );   
                       if(length(Ifind>0))
                           %assign found values of flux_trans to vertical distribution 
                           flux_vert_thislayer=flux_trans(Ifind);
                           %exclude nan values in transect (where parts of transect is out of model domain)
                           Igood=~isnan(flux_vert_thislayer);
                           flux_vert(k)=sum(flux_vert_thislayer(Igood)) ; %total flux of this layer 
                                                                          %(m^3/sec)*(concentration)
                           flux_vert(k)=flux_vert(k)/dz_vdis            ; %distribution vertically 
                                                                          %(m^2/sec)*(concentration)
                          if(~calculate_tracerflux)
                             flux_vert_thislayer_fw=flux_trans_fw(Ifind);
                             Igood=~isnan(flux_vert_thislayer_fw);
	                     flux_vert_fw(k)=sum(flux_vert_thislayer_fw(Igood)) ; %total flux of this layer (m^3/sec)
        	             flux_vert_fw(k)=flux_vert_fw(k)/dz_vdis            ; %distribution vertically (m^2/sec)
                          end
                       else
                           flux_vert(k)=0.0; 
                           if(~calculate_tracerflux)
                              flux_vert_fw(k)=0.0;
                           end
                       end
                   end

                   %collect all vertical flux distributions in time
                   if(it_all==1) 
                        flux_vert_all = zeros([NT_all,Nz_vert]); 
                   end
                   flux_vert_all(it_all,:)=flux_vert; 
                   if(~calculate_tracerflux)
                         if(it_all==1)
                            flux_vert_all_fw = zeros([NT_all,Nz_vert]);  %initalize flux_vert_all_fw
                         end
                         flux_vert_all_fw(it_all,:)=flux_vert_fw;
                   end
            else
                  %do nothing
            end

            display(['plotting ', parname, ' at time step: ', num2str(it)]);
 
            %plot the transect based on parameter choice
             r_plot=r;     r_plot(:,Ipath_out)=nan;
             Tz_plot=Tz;  Tz_plot(:,Ipath_out)=nan;

            %A) make the plot in matlab format and save the plot data into .mat files 
            if(strcmp(parname,'velocity'))

                 %plot velocity for each time step 
                 %1) plot in matlab format

                 if(yesno_matlabfig)
                   %plot velocity with vectors
                   subplot(2,1,1)
                      hold on;
                      %Wen Long customized quiver plotting program with u in km/s and v in m/s
                       quiverc_ukms_vms(r,Tz,u_p/1000.0,Tp3,1.2); colorbar; %plot vectors 
                                                                            %u velocity is scaled to km/sec
                       if(climits(2)>0)
                         caxis([0,climits(2)]);  %for quiverc_ukms_vms plot, color is based on 
                                                 %size of vector, and there is 
                       else                      %no negative value. so color axis should start from 0 
                         caxis([0,max(sqrt(u_p(:).^2+Tp3(:).^2))]);
                       end

                      hbottom=plot(r_plot(end,:)',Tz_plot(end,:)','k-');    %plot bottom with thick line
                      set(hbottom,'linewidth',3)
 
                      %make backoground black and boarders blue
                      set(gca, 'color', [0 0 0],'Xcolor','b','Ycolor','b');
                      set(gcf, 'color', [0 0 0]);

                      ylabel('Z(m)');
                      set(gca,'xlim',[0,max(r(:))]);
                      set(gca,'ylim',[zlimit(1),zlimit(2)+2]);  %offset zlimit(2) by 2 meter
                      title(['(u,w) on ' transname ' at time = ' num2str(htime(it),'%f8') ...
                            '(days)'],'fontsize',15,'color',[0 0 1]);

                   subplot(2,1,2)
                      pcolor(r,Tz,v_p);shading interp; colorbar
                      if(diff(climits(1:2))>0)
                         caxis(climits(1:2));
                      end
                      hold on;
                      plot(r_plot(:,:)',Tz_plot(:,:)','k-');
                      hbottom=plot(r_plot(end,:)',Tz_plot(end,:)','k-');  %plot bottom with thick line
                      set(hbottom,'linewidth',3)
                      xlabel('Transect distance (km)');ylabel('depth (m)');
                      ylabel('Z(m)');
                      set(gca,'xlim',[0,max(r(:))]);
                      set(gca,'ylim',[zlimit(1),zlimit(2)+2]);  %offset zlimit(2)  by 2 meter

                      set(gca, 'color', [0 0 0],'Xcolor','b','Ycolor','b');
                      set(gcf, 'color', [0 0 0]);

                      if(calculate_tracerflux)
                          title({['normal velocity on ' transname ' at time = ' num2str(htime(it),'%f8') '(days)']; ...
      	                     ['Total flux Q\_transect=' num2str(flux_total,'%f20') '(m3/s)*(conc.)'] } ...
                              ,'fontsize',15,'color',[0 0 1]);
		      else
                          title({['normal velocity on ' transname ' at time = ' num2str(htime(it),'%f8') '(days)']; ...
                             ['Total flux Q\_transect=' num2str(flux_total,'%f20') '(m3/s)'];  ...
                             ['Total fresh waterflux Q\_transect\_fw=' num2str(flux_total_fw,'%f20') '(m3/s)']; } ...
                              ,'fontsize',15,'color',[0 0 1]);
                      end
 
                    %debug
                    % display(['Debugging pause, enter to continue']); 
                    % pause

                    %plot vertical distribution of flow 
                    hflux_vert=figure('visible','off','position',[1,1,400,800]);
                    set(gcf,'renderer','zbuffer');
 
                    if(calculate_tracerflux)   %tracer flux
	                    hbar_vert=barh((Tz_vert(1:end-1)+Tz_vert(2:end))/2.0,flux_vert(1:end-1));
                                                       %horizontal bar chart
	                    title({[tracername ' flux per unit transect depth ((m^3/sec)*(conc.)/m)'];...
                                   [' at time = ' num2str(htime(it),'%f8') '(days)'] } ...
                            ,'fontsize',15);
	                    xlabel('flux per unit depth ((m^3/sec)*(conc)./m)','fontsize',15);
	                    ylabel('z(m)','fontsize',15);
	                    hflux_figfilename=[pltdirname parname,'_trans_', transname, '_', tracername, ...
                                               '_Qprofile_', num2str(it_all,'%010d'),'.png'];
	                    saveas(gcf,hflux_figfilename,'png');
	                    hflux_figfilename=[pltdirname parname,'_trans_', transname, '_', tracername, ...
                                               '_Qprofile_', num2str(it_all,'%010d'),'.fig'];
	                    saveas(gcf,hflux_figfilename,'fig');
		    else  %volume flux of bulk sea water  and fresh water volume flux
                          subplot(2,1,1)
         	            hbar_vert=barh((Tz_vert(1:end-1)+Tz_vert(2:end))/2.0,flux_vert(1:end-1));
                                                        %horizontal bar chart
                            title({['Flux per unit transect depth ((m^3/sec)/m)'];...
                              [' at time = ' num2str(htime(it),'%f8') '(days)'] } ...
                            ,'fontsize',15);
                            xlabel('flux per unit depth ((m^3/sec)./m)','fontsize',15);
                            ylabel('z(m)','fontsize',15);

			  subplot(2,1,2)
                            hbar_vert_fw=barh((Tz_vert(1:end-1)+Tz_vert(2:end))/2.0,flux_vert_fw(1:end-1));
                                                         %horizontal bar chart
                            title({['Freshwater flux per unit transect depth ((m^3/sec)/m)'];...
                                         [' at time = ' num2str(htime(it),'%f8') '(days)'] } ...
                            ,'fontsize',15);
                            xlabel('flux per unit depth ((m^3/sec)/m)','fontsize',15);
                            ylabel('z(m)','fontsize',15);
                          hflux_figfilename=[pltdirname parname,'_trans_', transname, ...
                                             '_Qprofile_fw_', num2str(it_all,'%010d'),'.png'];
                          saveas(gcf,hflux_figfilename,'png');
                          hflux_figfilename=[pltdirname parname,'_trans_', transname, ...
                                             '_Qprofile_fw_', num2str(it_all,'%010d'),'.fig'];
                          saveas(gcf,hflux_figfilename,'fig');
                    end
                    delete(hflux_vert);
                 end

                 %2)save plot data in matlab format 

                 if(calculate_tracerflux)   %tracer flux
                     
                    %save the data into matlab as well so that it can be replotted easily
                     tmpfilename=[pltdirname tracername '_' parname,'_trans_', transname, ...
                                 '_',num2str(nf,'%06d'),'_' num2str(it,'%010d'),'.mat'];
                     save(tmpfilename,'r',          ...   % distance along transect (km)
                                      'Tz',         ...   % vertical coordinate (m)
                                      'u_p',        ...   % velocity along the transect in transect surface
                                      'Tp3',        ...   % upward velocity (w) on the transect in transect surface
                                      'Tp4',        ...   % tracer
                                      'climits',    ...   % color limits
                                      'zlimit',     ...   % vertical coordinate limits
                                      'r_plot',     ...   % distance
                                      'Tz_plot',    ...   % vertical oordinate
                                      'transname',  ...   % name of transect
                                      'parname',    ...   % name of parameter
                                      'tracername', ...   % name of tracer
                                      'htime',      ...   % time (days)
                                      'it',         ...   % index in time array
                                      'flux_total', ...   % flow flux through the transect (m^3/sec)*(conc.)
                                      '-mat');

                 else %volume flux and fresh water flux
			 %save the data into matlab as well so that it can be replotted easily
                         tmpfilename=[pltdirname parname,'_trans_', transname, ...
                                  '_',num2str(nf,'%06d'),'_' num2str(it,'%010d'),'_fw.mat'];
                         save(tmpfilename,'r',          ...   % distance along transect (km)
                                      'Tz',         ...   % vertical coordinate (m)
                                      'u_p',        ...   % velocity along the transect in transect surface
                                      'Tp3',        ...   % upward velocity (w) on the transect in transect surface
                                      'Tp4',        ...   % tracer
                                      'climits',    ...   % color limits
                                      'zlimit',     ...   % vertical coordinate limits
                                      'r_plot',     ...   % distance
                                      'Tz_plot',    ...   % vertical oordinate
                                      'transname',  ...   % name of transect
                                      'parname',    ...   % name of parameter
                                      'tracername', ...   % name of tracer is S
                                      'htime',      ...   % time (days)
                                      'it',         ...   % index in time array
                                      'flux_total', ...   % flow flux through the transect (m^3/sec)
                                      'flux_total_fw',... % freshwater flux *m^3/sec)
                                      '-mat');
                 end

            else  %scalar only

                 %plot scalar only
                 %1)plot in matlab format 
                 if(yesno_matlabfig)
                   %plot pcolor if scalar
                   pcolor(r,Tz,Tp1);shading interp; colorbar; 
                   hold on;
                   plot(r_plot(:,:)',Tz_plot(:,:)','k-');
                   hbottom=plot(r_plot(end,:)',Tz_plot(end,:)','k-');  %plot bottom with thick line
                   set(hbottom,'linewidth',3)
                   xlabel('Transect distance (km)','color',[1 1 1]);ylabel('depth (m)','color',[1,1,1]);
                   ylabel('Z(m)');
                   set(gca,'xlim',[0,max(r(:))]);
                   set(gca,'ylim',[zlimit(1),zlimit(2)+2]);  %offset zlimit(2)  by 2 meter
                   title([parname parunit ' on ' transname  ' at time = ' num2str(htime(it),'%f8') ...
                         '(days)'],'fontsize',15,'color',[0 0 1]);

                   if(diff(climits(1:2))>0)
                      caxis(climits(1:2));
                   end
                   set(gca, 'color', [0 0 0],'Xcolor','b','Ycolor','b');
                   set(gcf, 'color', [0 0 0]);
                 end

                 %2) save the plot data into matlab data 
                   tmpfilename=[pltdirname parname,'_trans_', transname, '_', num2str(nf,'%06d'),...
                               '_',num2str(it,'%010d'),'.mat'];
                   save(tmpfilename,'r',         ...   % distance along transect (km)
                                    'Tz',        ...   % vertical coordinate (m)
                                    'Tp1',       ...   % parameter value 
                                    'climits',   ...   % color limits
                                    'zlimit',    ...   % vertical coordinate limits
                                    'r_plot',    ...   % distance
                                    'Tz_plot',   ...   % vertical oordinate
                                    'transname', ...   % name of transect
                                    'parname',   ...   % name of the parameter
                                    'htime',     ...   % time (days)
                                    'it',        ...   % index in time array
                                    '-mat');
            end

	    %B) Save the matlab plot into .png and .fig files 

            if(yesno_matlabfig)
	               %save plots
	               figfilename=[pltdirname parname,'_trans_', transname, ...
                                      '_',num2str(nf,'%06d'),'_', num2str(it,'%010d'),'.png'];
	               saveas(gcf,figfilename,'png');
	               figfilename=[pltdirname parname,'_trans_', transname, ...
                                      '_',num2str(nf,'%06d'),'_', num2str(it,'%010d'),'.fig'];
	               saveas(gcf,figfilename,'fig');
	               %clear figure
	               delete(hfig);
            end

           %C) plot in tecplot format

            %save the transect into tecplot data 
            if(strcmp(parname,'velocity'))  %save velocity transects
                   tdata2d.surfaces(1).zonename=[parname ' transect ' transname, ...
                                      ' 2D zone file number ' num2str(nf,'%06d') ...
                                      ' time reocrd ' num2str(it_all,'%8d')];
                   tdata2d.surfaces(1).x=r'*1000.0;        %'x' is in meter
                   tdata2d.surfaces(1).y=Tz_tec';          %'y' is actually vertical coordinate here and it is in meter
                   tdata2d.surfaces(1).z=zeros(size(r'));  %dummy z value
                   tdata2d.surfaces(1).v(1,:,:)=Tp1'; %'u' (m/s)
                   tdata2d.surfaces(1).v(2,:,:)=Tp2'; %'v' (m/s)
                   tdata2d.surfaces(1).v(3,:,:)=Tp3'; %'w' (m/s) 
                   tdata2d.surfaces(1).v(4,:,:)=u_p';      %'u_p' (m/s)  %horizontal veloity along transect 
                   tdata2d.surfaces(1).v(5,:,:)=v_p';      %'v_p' (m/s)  %horizontal velocity permendicular to transect
                   tdata2d.surfaces(1).varloc=0;           %variables are provided nodal
                   tdata2d.surfaces(1).order=3;            %surface is definded on (x,y) intecplot, i.e. ('r','z') here
                   tdata2d.surfaces(1).solutiontime=htime(it); 
                   
                   %need to worry about nan values in x, y, z, v here when
                   %the transects run out of domain (replace nans with
                   %-999999999.0 so that tecplot can blank them out
                   
                   inan_x=find(isfinite(tdata2d.surfaces(1).x)==0);
                   tdata2d.surfaces(1).x(inan_x)=-999999999.0;  %give big negative value so that we can blank it
                   inan_y=find(isfinite(tdata2d.surfaces(1).y)==0);
                   tdata2d.surfaces(1).y(inan_y)=-999999999.0;  %give big negative value so that we can blank it

                   inan_v=find(isfinite(tdata2d.surfaces(1).v)==0);    %
                   tdata2d.surfaces(1).v(inan_v)=-999999999.0;
                   tdata2d.surfaces(1).v(1,inan_x)=-999999999.0;
                   tdata2d.surfaces(1).v(1,inan_y)=-999999999.0;
                   tdata2d.surfaces(1).v(2,inan_x)=-999999999.0;
                   tdata2d.surfaces(1).v(2,inan_y)=-999999999.0;
                   tdata2d.surfaces(1).v(3,inan_x)=-999999999.0;
                   tdata2d.surfaces(1).v(3,inan_y)=-999999999.0;
                   tdata2d.surfaces(1).v(4,inan_x)=-999999999.0;
                   tdata2d.surfaces(1).v(4,inan_y)=-999999999.0;                   
                   tdata2d.surfaces(1).v(5,inan_x)=-999999999.0;
                   tdata2d.surfaces(1).v(5,inan_y)=-999999999.0;

                   %3D tecplot in x, y, z coordinates of the transect (x, y, z are 2D ordered arrays)
                   tdata3d.surfaces(1).zonename=[parname ' transect ' transname, ...
                                          ' 3D zone file number ' num2str(nf,'%06d') ...
                                          ' time record ' num2str(it_all,'%8d')];
                   tdata3d.surfaces(1).x=repmat(reshape(xp,[length(xp),1]),[1,size(r,1)]); 
                   %first dimension is number of points along xp
                   tdata3d.surfaces(1).y=repmat(reshape(yp,[length(yp),1]),[1,size(r,1)]);
                   tdata3d.surfaces(1).z=Tz_tec';   %Tz has size [Nz+1, length(r)]; have to transpose it to [length(r),Nz+1]
                   tdata3d.surfaces(1).v(1,:,:)=Tp1'; % 'u' (m/s)
                   tdata3d.surfaces(1).v(2,:,:)=Tp2'; % 'v' (m/s)
                   tdata3d.surfaces(1).v(3,:,:)=Tp3'; % 'w' (m/s)
                   tdata3d.surfaces(1).v(4,:,:)=u_p'; % 'u_p' (m/s) horizontal velocity component that is along transect
                   tdata3d.surfaces(1).v(5,:,:)=v_p'; % 'v_p' (m/s) horizontal velocity component that is 
                                                      %perpendicular to transect
                   tdata3d.surfaces(1).v(6,:,:)=Tz';  %last variable in v is -depth
                   tdata3d.surfaces(1).order=2;       %surface is IK-odered (I is first dimension of x, y, z, v, K 
                                                      %is 2nd dimension of of x, y, z, v)
                                                      %surfaces are defined on (x,z) plane
                   tdata3d.surfaces(1).solutiontime=htime(it);
                   %need to worry abut nan values in x, y, z, v
                   inan_x=find(isfinite(tdata3d.surfaces(1).x)==0);
                   tdata3d.surfaces(1).x(inan_x)=-999999999.0;  %give big negative value so that we can blank it
                   inan_y=find(isfinite(tdata3d.surfaces(1).y)==0);
                   tdata3d.surfaces(1).y(inan_y)=-999999999.0;  %give big negative value so that we can blank it
                   inan_z=find(isfinite(tdata3d.surfaces(1).z)==0);
                   tdata3d.surfaces(1).z(inan_z)=-999999999.0;   
                   inan_v=find(isfinite(tdata3d.surfaces(1).v)==0);
                   tdata3d.surfaces(1).v(inan_v)=-999999999.0;
                   tdata3d.surfaces(1).v(1,inan_x)=-999999999.0;
                   tdata3d.surfaces(1).v(1,inan_y)=-999999999.0;
                   tdata3d.surfaces(1).v(1,inan_z)=-999999999.0;                   
                   tdata3d.surfaces(1).v(2,inan_x)=-999999999.0;
                   tdata3d.surfaces(1).v(2,inan_y)=-999999999.0;
                   tdata3d.surfaces(1).v(2,inan_z)=-999999999.0;                   
                   tdata3d.surfaces(1).v(3,inan_x)=-999999999.0;
                   tdata3d.surfaces(1).v(3,inan_y)=-999999999.0;
                   tdata3d.surfaces(1).v(3,inan_z)=-999999999.0;                   
                   tdata3d.surfaces(1).v(4,inan_x)=-999999999.0;
                   tdata3d.surfaces(1).v(4,inan_y)=-999999999.0;
                   tdata3d.surfaces(1).v(4,inan_z)=-999999999.0;                   
                   tdata3d.surfaces(1).v(5,inan_x)=-999999999.0;
                   tdata3d.surfaces(1).v(5,inan_y)=-999999999.0;
                   tdata3d.surfaces(1).v(5,inan_z)=-999999999.0;                   
                   tdata3d.surfaces(1).v(6,inan_x)=-999999999.0;
                   tdata3d.surfaces(1).v(6,inan_y)=-999999999.0;
                   tdata3d.surfaces(1).v(6,inan_z)=-999999999.0;                   
                   
                   %for 3D transect tecplot, also add bathymtry as a zone
                   %ideally we do not have to have this for every time step
                   %but that then makes it hard to do movie with time steps
                   
                   tdata3d.FEsurfaces(1).zonename=['PSWQM 2.0 grid file number '  ...
                                   num2str(nf,'%06d') ' time record ' num2str(it_all,'%8d')];
                   tdata3d.FEsurfaces(1).x=x_n;   %x of all nodes
                   tdata3d.FEsurfaces(1).y=y_n;   %y of all nodes
                   tdata3d.FEsurfaces(1).z=-h_n;  %z of all nodes is given as -h
                   tdata3d.FEsurfaces(1).order=3; %order 3, surface is defined on (x, y) coordinates
                   tdata3d.FEsurfaces(1).e2n=e2n(:,2:4);
                   tdata3d.FEsurfaces(1).v(1,:)=zeros(size(h_n')); %dummy values as place holder for u
                   tdata3d.FEsurfaces(1).v(2,:)=zeros(size(h_n')); %dummy values as place holder for v
                   tdata3d.FEsurfaces(1).v(3,:)=zeros(size(h_n')); %dummy values as place holder for w
                   tdata3d.FEsurfaces(1).v(4,:)=zeros(size(h_n')); %dummy values as place holder for u_p
                   tdata3d.FEsurfaces(1).v(5,:)=zeros(size(h_n')); %dummy values as place holder for v_p
                   tdata3d.FEsurfaces(1).v(6,:)=-h_n';         %-depth give sthe depth on this surfacs
                   tdata3d.FEsurfaces(1).solutiontime=htime(it);
                   
            else %save scalar transects
                   
                   tdata2d.surfaces(1).zonename=[parname ' transect ' transname, ...
                                      ' 2D zone file number ' num2str(nf,'%06d') ... 
                                      ' time reocrd ' num2str(it_all,'%8d')];
                   tdata2d.surfaces(1).x=r'*1000.0;        %'x' is in meter
                   tdata2d.surfaces(1).y=Tz_tec';          %'y' is actually vertical coordinate here and it is in meter
                   tdata2d.surfaces(1).z=zeros(size(r'));  %dummy z value
                   tdata2d.surfaces(1).v(1,:,:)=Tp1';      %parname (T, S, kh)
                   tdata2d.surfaces(1).varloc=0;           %variables are provided nodal
                   tdata2d.surfaces(1).order=3;            %surface is definded on (x,y) intecplot, i.e. ('r','z') here
                   tdata2d.surfaces(1).solutiontime=htime(it); 
                   
                   %need to worry about nan values in x, y, z, v here when
                   %the transects run out of domain (replace nans with
                   %-999999999.0 so that tecplot can blank them out
                   
                   inan_x=find(isfinite(tdata2d.surfaces(1).x)==0);
                   tdata2d.surfaces(1).x(inan_x)=-999999999.0;  %give big negative value so that we can blank it
                   inan_y=find(isfinite(tdata2d.surfaces(1).y)==0);
                   tdata2d.surfaces(1).y(inan_y)=-999999999.0;  %give big negative value so that we can blank it

                   inan_v=find(isfinite(tdata2d.surfaces(1).v)==0);     %
                   tdata2d.surfaces(1).v(inan_v)=-999999999.0;
                   tdata2d.surfaces(1).v(1,inan_x)=-999999999.0;
                   tdata2d.surfaces(1).v(1,inan_y)=-999999999.0;

                   %3D tecplot in x, y, z coordinates of the transect (x, y, z are 2D ordered arrays)
                   tdata3d.surfaces(1).zonename=[parname ' transect ' transname, ...
                                      ' 3D zone file number ' num2str(nf,'%06d') ...
                                      ' time record ' num2str(it,'%8d')];
                   tdata3d.surfaces(1).x=repmat(reshape(xp,[length(xp),1]),[1,size(r,1)]); 
                                                      %first dimension is number of points along xp
                   tdata3d.surfaces(1).y=repmat(reshape(yp,[length(yp),1]),[1,size(r,1)]);
                   tdata3d.surfaces(1).z=Tz_tec';     %Tz has size [Nz+1, length(r)]; have to transpose 
                                                      %it to [length(r),Nz+1]
                   tdata3d.surfaces(1).v(1,:,:)=Tp1'; %parname (T, S, kh)
                   tdata3d.surfaces(1).v(2,:,:)=Tz_tec';  %last variable in v is -depth
                   tdata3d.surfaces(1).order=2;       %surface is IK-odered (I is first dimension of x, y, z,
                                                      %v, K is 2nd dimension of of x, y, z, v)
                                                      %surfaces are defined on (x,z) plane
                   tdata3d.surfaces(1).solutiontime=htime(it);

                   %need to worry abut nan values in x, y, z, v
                   inan_x=find(isfinite(tdata3d.surfaces(1).x)==0);
                   tdata3d.surfaces(1).x(inan_x)=-999999999.0;  %give big negative value so that we can blank it
                   inan_y=find(isfinite(tdata3d.surfaces(1).y)==0);
                   tdata3d.surfaces(1).y(inan_y)=-999999999.0;  %give big negative value so that we can blank it
                   inan_z=find(isfinite(tdata3d.surfaces(1).z)==0);
                   tdata3d.surfaces(1).z(inan_z)=-999999999.0;  %give big negative value so that we can blank it

                   inan_v=find(isfinite(tdata3d.surfaces(1).v)==0);
                   tdata3d.surfaces(1).v(inan_v)=-999999999.0;
                   tdata3d.surfaces(1).v(1,inan_x)=-999999999.0;
                   tdata3d.surfaces(1).v(1,inan_y)=-999999999.0;
                   tdata3d.surfaces(1).v(1,inan_z)=-999999999.0;                   
                   tdata3d.surfaces(1).v(2,inan_x)=-999999999.0;
                   tdata3d.surfaces(1).v(2,inan_y)=-999999999.0;
                   tdata3d.surfaces(1).v(2,inan_z)=-999999999.0;                   
                   
                   %for 3D transect tecplot, also add bathymtry as a zone
                   %ideally we do not have to have this for every time step
                   %but that then makes it hard to do movie with time steps
                   
                   tdata3d.FEsurfaces(1).zonename=['PSWQM 2.0 grid time record file number ' ... 
                                                   num2str(nf,'%06d') ' time record' num2str(it,'%8d')];
                   tdata3d.FEsurfaces(1).x=x_n;   %x of all nodes
                   tdata3d.FEsurfaces(1).y=y_n;   %y of all nodes
                   tdata3d.FEsurfaces(1).z=-h_n;  %z of all nodes is given as -h
                   tdata3d.FEsurfaces(1).order=3; %order 3, surface is defined on (x, y) coordinates
                   tdata3d.FEsurfaces(1).e2n=e2n(:,2:4);
                   tdata3d.FEsurfaces(1).v(1,:)=zeros(size(h_n')); %dummy values as place holder for parname (T,S, kh)
                   tdata3d.FEsurfaces(1).v(2,:)=-h_n';         %-depth give sthe depth on this surfacs
                   tdata3d.FEsurfaces(1).solutiontime=htime(it);
                   
            end

            %save the tecplot data into matlab file first, will load them after the nf loop
             tmpfilename=[pltdirname parname,'_trans_', transname, '_tecplot_', ... 
                          num2str(nf,'%06d') '_' num2str(it,'%010d'),'.mat'];
             save(tmpfilename,'tdata2d', ...
                              'tdata3d', ...
                              '-mat');

     end  %end it_global loop


 %Aggregrate all tecplot plots of individual time records into a big tecplot history file
 %so that it can be played as an animation in tecplot

     tdata2d_all=[];
     tdata3d_all=[];
     it_all=0;
     
     for it_global=1:t_int:NT_global

         it_all=it_all+1;

         nf=filenum_global(it_global);
         it_local=it_local_global(it_global);
         it=it_local;
         %find the tecplot data stored in matlab to load

         display(['making tecplot history file, nf=', num2str(nf,'%06d')]);

         tmpfilename=[pltdirname parname,'_trans_', transname, '_tecplot_', ...
                            num2str(nf,'%06d') '_' num2str(it,'%010d'),'.mat'];
         load(tmpfilename); %load tdata2d and tdata3d
         if(it_all==1)
            tdata2d_all.Nvar=tdata2d.Nvar;
            tdata2d_all.varnames=tdata2d.varnames;
         end

         tdata2d_all.surfaces(it_all)=tdata2d.surfaces;
         if(it_all==1)
            tdata3d_all.Nvar=tdata3d.Nvar;
            tdata3d_all.varnames=tdata3d.varnames;
         end
         tdata3d_all.surfaces(it_all)=tdata3d.surfaces;
         tdata3d_all.FEsurfaces(it_all)=tdata3d.FEsurfaces;
     end

     %output the transect time history into tecplot format
     tecfile2d=[pltdirname parname,'_trans_', transname,'_', num2str(his_fn_start,'%06d'), ...
                                   '_', num2str(his_fn_end,'%06d'), '_2D.plt'];
     mat2tecplot(tdata2d_all,tecfile2d);

     tecfile3d=[pltdirname parname,'_trans_', transname,'_', num2str(his_fn_start,'%06d'), ...
                                   '_', num2str(his_fn_end,'%06d'), '_3D.plt'];
     mat2tecplot(tdata3d_all,tecfile3d);

 %Create a time average for the 2D and 3D zones 
     tdata2d_timeavg.surfaces(1).zonename=[parname ' transect ' transname, ' 2D zone time average from ' ...
                                           num2str(htime_all(1),'%8.2f') ' to ' num2str(htime_all(end),'%8.2f')];
     tdata2d_timeavg.surfaces(1).x=zeros(size(tdata2d.surfaces(1).x));
     tdata2d_timeavg.surfaces(1).y=zeros(size(tdata2d.surfaces(1).y));
     tdata2d_timeavg.surfaces(1).z=zeros(size(tdata2d.surfaces(1).z));
     tdata2d_timeavg.surfaces(1).v=zeros(size(tdata2d.surfaces(1).v));
     tdata2d_timeavg.surfaces(1).varloc=0;           %variables are provided nodal
     tdata2d_timeavg.surfaces(1).order=3;            %surface is definded on (x,y) intecplot, i.e. ('r','z') here

     tdata3d_timeavg.surfaces(1).zonename=[parname ' transect ' transname, ' 3D zone time average from ' ...
                                           num2str(htime_all(1),'%8.2f') ' to ' num2str(htime_all(end),'%8.2f')];
     tdata3d_timeavg.surfaces(1).x=zeros(size(tdata3d.surfaces(1).x));
     tdata3d_timeavg.surfaces(1).y=zeros(size(tdata3d.surfaces(1).y));
     tdata3d_timeavg.surfaces(1).z=zeros(size(tdata3d.surfaces(1).z));
     tdata3d_timeavg.surfaces(1).v=zeros(size(tdata3d.surfaces(1).v));
     tdata3d_timeavg.surfaces(1).order=2;    %surface is IK-odered (I is first dimension of x, y, z, 
                                             %v, K is 2nd dimension of of x, y, z, v)
                                             %surfaces are defined on (x,z) plane
     %add bathymetry zone
     tdata3d_timeavg.FEsurfaces(1).zonename=['PSWQM 2.0 grid'];
     tdata3d_timeavg.FEsurfaces(1).x=x_n;   %x of all nodes
     tdata3d_timeavg.FEsurfaces(1).y=y_n;   %y of all nodes
     tdata3d_timeavg.FEsurfaces(1).z=-h_n;  %z of all nodes is given as -h
     tdata3d_timeavg.FEsurfaces(1).order=3; %order 3, surface is defined on (x, y) coordinates
     tdata3d_timeavg.FEsurfaces(1).e2n=e2n(:,2:4);
     for ivar=1:tdata3d_timeavg.Nvar-3-1   %first 3 is x, y,z, then rest of them has be filled up with dummy valuses
        tdata3d_timeavg.FEsurfaces(1).v(ivar+3,:)=zeros(size(h_n'));%dummy values as place holder for parname (T,S, kh)
     end
     tdata3d_timeavg.FEsurfaces(1).v(tdata3d_timeavg.Nvar,:)=-h_n'; %-depth give sthe depth on this surfacs (last variable)

     %do time average
     for it=1:length(htime_all)
         tdata2d_timeavg.surfaces(1).x=tdata2d_timeavg.surfaces(1).x+tdata2d_all.surfaces(it).x;
         tdata2d_timeavg.surfaces(1).y=tdata2d_timeavg.surfaces(1).y+tdata2d_all.surfaces(it).y;
         tdata2d_timeavg.surfaces(1).z=tdata2d_timeavg.surfaces(1).z+tdata2d_all.surfaces(it).z;
         tdata2d_timeavg.surfaces(1).v=tdata2d_timeavg.surfaces(1).v+tdata2d_all.surfaces(it).v;

         tdata3d_timeavg.surfaces(1).x=tdata3d_timeavg.surfaces(1).x+tdata3d_all.surfaces(it).x;
         tdata3d_timeavg.surfaces(1).y=tdata3d_timeavg.surfaces(1).y+tdata3d_all.surfaces(it).y;
         tdata3d_timeavg.surfaces(1).z=tdata3d_timeavg.surfaces(1).z+tdata3d_all.surfaces(it).z;
         tdata3d_timeavg.surfaces(1).v=tdata3d_timeavg.surfaces(1).v+tdata3d_all.surfaces(it).v;
     end

     tdata2d_timeavg.surfaces(1).x=tdata2d_timeavg.surfaces(1).x/length(htime_all);
     tdata2d_timeavg.surfaces(1).y=tdata2d_timeavg.surfaces(1).y/length(htime_all);
     tdata2d_timeavg.surfaces(1).z=tdata2d_timeavg.surfaces(1).z/length(htime_all);
     tdata2d_timeavg.surfaces(1).v=tdata2d_timeavg.surfaces(1).v/length(htime_all);

     tdata3d_timeavg.surfaces(1).x=tdata3d_timeavg.surfaces(1).x/length(htime_all);
     tdata3d_timeavg.surfaces(1).y=tdata3d_timeavg.surfaces(1).y/length(htime_all);
     tdata3d_timeavg.surfaces(1).z=tdata3d_timeavg.surfaces(1).z/length(htime_all);
     tdata3d_timeavg.surfaces(1).v=tdata3d_timeavg.surfaces(1).v/length(htime_all);

     %write time average to tecplot binary
     tecfile2d_timeavg=[pltdirname parname,'_trans_', transname,'_',...
             num2str(his_fn_start,'%06d'), '_', num2str(his_fn_end,'%06d'), '_2D_timeavg.plt'];
     mat2tecplot(tdata2d_timeavg,tecfile2d_timeavg);
     tecfile3d_timeavg=[pltdirname parname,'_trans_', transname,'_',...
             num2str(his_fn_start,'%06d'), '_', num2str(his_fn_end,'%06d'), '_3D_timeavg.plt'];
     mat2tecplot(tdata3d_timeavg,tecfile3d_timeavg);

     %Create run mean of the history and save in tecplot format 
        %to add calculation of runmean here ...
        %to save run mean in tecplot here

     %Finally save everything about the tecplot data  in a matlab file for future inspection
     save([pltdirname parname '_' transname '_calculation_' num2str(his_fn_start,'%06d') ... 
                              '_' num2str(his_fn_end,'%06d') '.mat'],                    ...
                       'tdata2d_all','tdata3d_all','tdata2d_timeavg','tdata3d_timeavg');
     clear 'tdata2d' 'tdata3d' 'tdata2d_timeavg' 'tdata3d_timeavg' ;
     clear 'tdata2d_all', tdata3d_all';

 %
 %plot flux calculation as time series and vertical distribution in both matlab figure, tecplot figure formats
 %
     if(strcmp(parname,'velocity'))
       %A) plot figures in matlab format if chosen so
         if(yesno_matlabfig)
	 %1) plot and save time series of total flux records
	  hflux=figure('visible','off','position',[1,1,1000,600]);
	  set(gcf,'renderer','zbuffer');

	  plot(htime_all,flux_all,'b-'); 
	  title({['Total flux Q\_transect (m^3/s)*(conc.)']},'fontsize',15);
	  ylabel(['Q (m^3/s)*(conc.)'],'fontsize',15);
	  xlabel('Time(days)','fontsize',15); 
	  figfilename=[pltdirname parname,'_trans_', transname, '_', ...
            num2str(his_fn_start,'%06d'), '_', num2str(his_fn_end,'%06d'),'_flux','.png'];
	  saveas(gcf,figfilename,'png');

	  figfilename=[pltdirname parname,'_trans_', transname, '_', ...
            num2str(his_fn_start,'%06d'), '_', num2str(his_fn_end,'%06d'),'_flux','.fig'];
	  saveas(gcf,figfilename,'fig');

	  delete(hflux);
	 %2) plot vertical distribution of time averaged flux
	  hflux_vert=figure('visible','off','position',[1,1,400,800]);
	  set(gcf,'renderer','zbuffer');

	  hbar_vert=barh((Tz_vert(1:end-1)+Tz_vert(2:end))/2.0, ...
                   mean(flux_vert_all(:,1:end-1),1));  %horizontal bar chart of time averaged flux
	  title({['Time averaged flux per unit transect depth ((m^3/sec)*(conc.)/m)'];...
		 ['time=[' num2str(htime_all(1),'%f8')  ','  num2str(htime_all(end),'%f8') '] (days)'] } ...
		 ,'fontsize',15);
	  xlabel('<flux per unit depth> ((m^3/sec)*(conc.)/m)','fontsize',15);
	  ylabel('z(m)','fontsize',15);
	  hflux_figfilename=[pltdirname parname,'_trans_', transname, '_', ...
                  num2str(his_fn_start,'%06d'), '_', num2str(his_fn_end,'%06d'),'_Qprof_avg','.png'];
	  saveas(gcf,hflux_figfilename,'png');
	  hflux_figfilename=[pltdirname parname,'_trans_', transname, '_', ...
                  num2str(his_fn_start,'%06d'), '_', num2str(his_fn_end,'%06d'),'_Qprof_avg','.fig'];
	  saveas(gcf,hflux_figfilename,'fig');

	  delete(hflux_vert);

	  if(~calculate_tracerflux)  %plot fresh water flux time series and vertical distribution

		 %plot and save time series of total flux
		  hflux=figure('visible','off','position',[1,1,1000,600]);
		  set(gcf,'renderer','zbuffer');

		  plot(htime_all,flux_all_fw,'b-');
		  title({['Fresh water flux Q\_transect (m^3/s)']},'fontsize',15);
		  ylabel(['Q (m^3/s)'],'fontsize',15);
		  xlabel('Time(days)','fontsize',15);
		  figfilename=[pltdirname parname,'_trans_', transname, '_', ...
                    num2str(his_fn_start,'%06d'), '_', num2str(his_fn_end,'%06d'),'_fw_flux','.png'];
		  saveas(gcf,figfilename,'png');

		  figfilename=[pltdirname parname,'_trans_', transname, '_', ...
                    num2str(his_fn_start,'%06d'), '_', num2str(his_fn_end,'%06d'),'_fw_flux','.fig'];
		  saveas(gcf,figfilename,'fig');

		  delete(hflux);

		 %plot vertical distribution of time averaged flux
		  hflux_vert=figure('visible','off','position',[1,1,400,800]);
		  set(gcf,'renderer','zbuffer');

		  hbar_vert=barh((Tz_vert(1:end-1)+Tz_vert(2:end))/2.0, ...
                      mean(flux_vert_all_fw(:,1:end-1),1));  %horizontal bar chart of time averaged flux
		  title({['Time averaged fresh water flux per unit transect depth ((m^3/sec)/m)'];...
			 ['time=[' num2str(htime_all(1),'%f8')  ','  num2str(htime_all(end),'%f8') '] (days)'] } ...
			 ,'fontsize',15);
		  xlabel('<flux per unit depth> ((m^3/sec)/m)','fontsize',15);
		  ylabel('z(m)','fontsize',15);

		  hflux_figfilename=[pltdirname parname,'_trans_', transname, '_',...
                          num2str(his_fn_start,'%06d'), '_', num2str(his_fn_end,'%06d'),'_fw_Qprof_avg','.png'];
		  saveas(gcf,hflux_figfilename,'png');
		  hflux_figfilename=[pltdirname parname,'_trans_', transname, '_',...
                          num2str(his_fn_start,'%06d'), '_', num2str(his_fn_end,'%06d'),'_fw_Qprof_avg','.fig'];
		  saveas(gcf,hflux_figfilename,'fig');
		  delete(hflux_vert);
	  end
	 end  %end of matlab plot of flux 

        %B) save the tracer flux or freshwater flux into a matlab file (both time series and vertical profile)
	  if(calculate_tracerflux) %tracer flux
		  %save the flux data into a matlab file 
		  datafilename=[pltdirname parname,'_trans_', transname, '_',tracername,'_Qprof.mat'];

		  readmetxt={'htime_all: time in days'; ...
			     'Tz_vert: vertical coordinate (m) of vertical distribution of tracer flux';...
			     'flux_vert_all: time series of vertical distribution of tracer flux through the transect';
			     'flux_all : time series of tracer flux integrated over the whole transect';
			     'transname: name of the transect';
			     'xp, yp: x and y cordiantes of the trasnect';
			     'Ipath_in: indices of points in xp, yp that are in the domain';
			     'Ipath_out: indices of points in xp, yp that are out of the domain';};

		 save(datafilename,'htime_all','transname','flux_all','Tz_vert','flux_vert_all', ...
                                   'Ipath_in','Ipath_out','xp','yp','readmetxt','-mat');

	 else  %volume and freshwater volume flux

		  %save the flux data into a matlab file
		  datafilename=[pltdirname parname,'_trans_', transname, '_fw_Qprof.mat'];

		  readmetxt={'htime_all: time in days'; ...
			     'Tz_vert: vertical coordinate (m) of vertical distribution of flux';...
			     'flux_vert_all: time series of vertical distribution of flow through the transect';
			     'flux_all : time series of flow integrated over the whole transect';
 	    'flux_vert_all_fw: time series of vertical distribution of fresh water flow through the transect';
			     'flux_all_fw : time series of fresh water flux integrated over the whole transect';
			     'transname: name of the transect';
			     'xp, yp: x and y cordiantes of the trasnect';
			     'Ipath_in: indices of points in xp, yp that are in the domain';
			     'Ipath_out: indices of points in xp, yp that are out of the domain';};

		 save(datafilename,'htime_all','transname','flux_all','flux_all_fw','Tz_vert','flux_vert_all',...
                                   'flux_vert_all_fw', ...
				   'Ipath_in','Ipath_out','xp','yp','readmetxt','-mat');
	 end

       %C) save the results into tecplot format

         %time series
	 tdata_TS=[] ;      % time series 
	 tdata_Prf=[];      % time series of vertical profile
	 tdata_Prf_mean=[]; % mean profile (time averagedrecords)
	 
	 tdata_TS.Nvar=2;
	 tdata_TS.varnames={'time','flux'};
	 tdata_TS.lines(1).x=htime_all;  %time 
	 tdata_TS.lines(1).y=flux_all;   %flux (m^3*(conc.)/sec)

	 if(~calculate_tracerflux)
	    tdata_TS.lines(2).x=htime_all;     %time
	    tdata_TS.lines(2).y=flux_all_fw;   %fresh water flux (m^3/sec)
	 end
	 
	 tecfile_TS=[pltdirname parname,'_trans_', transname, '_', ...
                num2str(his_fn_start,'%06d'), '_', num2str(his_fn_end,'%06d'),'_flux_timeseries.plt'];
	 mat2tecplot(tdata_TS,tecfile_TS);
	 
	 tdata_Prf.Nvar=2;
	 tdata_Prf.varnames={'z','flux'};
	 for it=1:length(htime_all)
	    if(calculate_tracerflux)
		    tdata_Prf.lines(it).zonename=[tracername ' flux zone of transect ' ...
                         transname ' at time ' num2str(htime_all(it),'%8.4f')];
	    else
		    tdata_Prf.lines(it).zonename=['volume flux zone of transect ' ...
                         transname ' at time ' num2str(htime_all(it),'%8.4f')];
	    end
	    tdata_Prf.lines(it).x=Tz_vert';             %transpose from row to column % coordinate (vertical)
	    tdata_Prf.lines(it).y=flux_vert_all(it,:)'; %transpose from it's row to column  %flux at this time
	    tdata_Prf.lines(it).solutiontime=htime_all(it);  
	 end

	 if(calculate_tracerflux)  %tracer flux
	     tecfile_Prf=[pltdirname parname,'_trans_', transname, '_', tracername, '_', ...
                               num2str(his_fn_start,'%06d'), ...
			  '_', num2str(his_fn_end,'%06d') '_fluxprofile_timeseries.plt'];
	 else  %volume flux
	     tecfile_Prf=[pltdirname parname,'_trans_', transname, '_volume_', ...
                               num2str(his_fn_start,'%06d'), ...
			  '_', num2str(his_fn_end,'%06d') '_fluxprofile_timeseries.plt'];

	 end 
	 mat2tecplot(tdata_Prf,tecfile_Prf);

	 %fresh water volume flux
	 if(~calculate_tracerflux)   %also do fresh water flux
	     tdata_Prf_fw=[];        % time series of vertical profile
	     tdata_Prf_fw.Nvar=2;
	     tdata_Prf_fw.varnames={'z','flux'};
	     for it=1:length(htime_all)
	       tdata_Prf_fw.lines(it).zonename=['freshwater flux zone of transect ' ...
                     transname ' at time ' num2str(htime_all(it),'%8.4f')];
	       tdata_Prf_fw.lines(it).x=Tz_vert';                %transpose from row to column % coordinate (vertical)
	       tdata_Prf_fw.lines(it).y=flux_vert_all_fw(it,:)'; %transpose from it's row to column  %flux at this time
	       tdata_Prf_fw.lines(it).solutiontime=htime_all(it);
	     end
	     tecfile_Prf_fw=[pltdirname parname,'_trans_', transname, '_volume_fw_', ...
                             num2str(his_fn_start,'%06d'), ...
			     '_', num2str(his_fn_end,'%06d') '_fluxprofile_timeseries.plt'];
	     mat2tecplot(tdata_Prf_fw,tecfile_Prf_fw);
	 end

         %vertical profile distribution
	 %time average to get mean profile

	 tdata_Prf_mean.Nvar=2;
	 tdata_Prf_mean.varnames={'z','flux'};

	 if(calculate_tracerflux) %tracer flux
		 tdata_Prf_mean.lines(1).zonename=[tracername ' mean flux zone transect ' ...
                 transname ' during time from ' ... 
		 num2str(htime_all(1),'%8.4f') ' to ' num2str(htime_all(end),'%8.4f')];
	 else  %volume flux
		 tdata_Prf_mean.lines(1).zonename=['volume time mean flux zone transect ' ...
                 transname ' during time from ' ...
	         num2str(htime_all(1),'%8.4f') ' to ' num2str(htime_all(end),'%8.4f')];

	 end
	 tdata_Prf_mean.lines(1).x=((Tz_vert(1:end-1)+Tz_vert(2:end))/2.0)';   %coordinate (vertical) (m)
									       %transpose from row to column 
	 tdata_Prf_mean.lines(1).y=(mean(flux_vert_all(:,1:end-1),1))';        %flux after time average (m^3*(conc.)/sec)
									       %transpose from it's row to column
	 if(calculate_tracerflux)
		tecfile_Prf_mean=[pltdirname parname,'_trans_', transname,'_'   ...
                                  tracername '_', num2str(his_fn_start,'%06d'), ...
				   '_' num2str(his_fn_end,'%06d'), '_fluxprofile_mean.plt'];
	 else
		tecfile_Prf_mean=[pltdirname parname,'_trans_', transname,'_volume_', ...
                                  num2str(his_fn_start,'%06d'), ...
				   '_' num2str(his_fn_end,'%06d'), '_fluxprofile_mean.plt'];
	 end

	 mat2tecplot(tdata_Prf_mean,tecfile_Prf_mean);

	 %fresh water volume flux if not calculating traceer flux
	 if(~calculate_tracerflux)  %fresh water time  mean flux vertical distribution

	   tdata_Prf_mean_fw=[];        % time series of vertical profile
	   tdata_Prf_mean_fw.Nvar=2;
	   tdata_Prf_mean_fw.varnames={'z','flux'};
	   tdata_Prf_mean_fw.lines(1).zonename=['fresh water mean flux zone transect ' ...
                                            transname ' during time from ' ...
					    num2str(htime_all(1),'%8.4f') ' to ' ...
                                            num2str(htime_all(end),'%8.4f')];
	   tdata_Prf_mean_fw.lines(1).x=((Tz_vert(1:end-1)+Tz_vert(2:end))/2.0)';   %coordinate (vertical) (m)
						       %transpose from row to column
	   tdata_Prf_mean_fw.lines(1).y=(mean(flux_vert_all(:,1:end-1),1))';   %flux after time average (m^3*(conc.)/sec)
									       %transpose from it's row to column
	   tecfile_Prf_mean_fw=[pltdirname parname,'_trans_', 
                                transname,'_volume_fw_', num2str(his_fn_start,'%06d'), ...
			        '_' num2str(his_fn_end,'%06d'), '_fluxprofile_mean.plt'];

	   mat2tecplot(tdata_Prf_mean_fw,tecfile_Prf_mean_fw);

	 end
		   
	 %clear 'tdata_TS' 'tdata_Prf' 'tdata_Prf_mean'

     end %end of "if(strcmp(parname,'velocity'))" block

     display(['Success, plots are made in folder ' pltdirname]);


