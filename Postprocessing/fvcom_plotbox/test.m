%function fvcom_his2flux_fv(his2flux)
%
% fvcom_his2flux_fv(his2flux)
%
%    This program reads fvcom model outputs and then calculates volume flux and tracer flux across a transect.
%
%    Input --- his2flux, a structure that stores various information about the fvcom model,
%              the transect and also the output file names
%
%    Dependency: fvcom_medm2txt.f90             compiled by fortran if history files are in medm format
%                InPolygon.c                    compiled by mex in matlab
%                polyxline.m                    code to calculate intersection of a polygon with a line
%                triangle_grid_edge.m           code to calculate geometry of the grid
%                fvcom_element_area.m           code to calculate areas of the grids
%                shape_coef_gcn.m               code to calculate shape functions of the grids
%                fvcom_tce_polygon.m            code to clculate TCE polygons of the FVCOM model grid
%                fvcom_segment_cutout_polygon.m code to calculate the cutout polygons by a straight line going
%                                                    through all TCE polygons
%                adv_s.m                        code to calculate flux on all TCE edges
%                mat2tecplot.m                  code to convert matlab data to tecplot
%                mexcdf                          matlab toolbox to load netcdf file if history files are in netcdf
%
%    Platform: Unix/Linux, Mac OS
%
%    Example 1) use netcdf files
%             his2flux.fvcom_mat_grid_file = './chn_island.mat';                        %matlab format of the model grid
%             his2flux.his_dir           = '/home/long075/hitach/TK_to_Wen_CHN-HYD_6_2013_HCC_12_D_25_SD_25_M2K1_poseidon';
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
%             fvcom_his2flux_fv(his2flux); 
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
%             fvcom_his2flux_fv(his2flux);
%
%    Example 3) use netcdf, no island channel case, with salt, with tide etc
%
%             his2flux.fvcom_mat_grid_file = './CHN_HYD.mat';                        %matlab format of the model grid
%             his2flux.his_dir           = '/home/long075/mypsfvm/trunk/example/channel-with-salt-tide/model_output/netcdf';';
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
%             his2flux.xt                =[-44363,-44363];  %x coord of the starting and ending points of the transect
%             his2flux.yt                =[-539,3379];      %y coord of the starting and ending points of the transect
%             his2flux.ds_transect       = 100.0;           %step size of discretizing the transect (m)
%             his2flux.t_int             = 1;               %time record interval in processing the history outputs
%                                                           %     1 for every record,2 for every other record,...
%             his2flux.plotoutdir        =['./' his2flux.transname '/'];  %output plt dirname
%             fvcom_his2flux_fv(his2flux);
%
%
% Wen Long @ PNNL, May  3rd, 2013
%

%
%obtain input parameters from his2flux structure
%

%     if(nargin<1)
%      error(['oops,missing his2flux argument, quitting ...']);
%     end

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
     have_dst	      = ~isempty(find(strcmp(fieldnames_arg,'ds_transect')==1));
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
%----

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
			if (~(strcmp(tracername,'S') || strcmp(tracername, 'T')))
                            error(['tracer name must be ''S'' or ''T'' for medm history outputs']);
                        end
               end
               if(read_his_nc)
                        if (~(strcmp(tracername, 'salinity') || strcmp( tracername, 'temp')))
                            error(['tracer name must be ''salinity'' or ''temp'' for netcdf history outputs']);
                        end
               end
            end
         end


         %plot a grid and let user choose transect

         yesno_matlabfig=0;  %no matlab plots to be made


         %get information of the transect 

         if (isempty(xt)||isempty(yt) )
            error(['xt and yt should not be empty']);
         end

         if(length(xt(:)) ~=2 || length(yt(:)) ~=2)
            error('xt and yt should have length==2'); 
         end

	 s = [0 cumsum(sqrt(diff(xt).^2+diff(yt).^2))' ]';  %distance along the transect

	 if(ds_transect >=0.5*max(s(:)))
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

     [alpha,a1u,a2u,aw0,awx,awy]=shape_coef_gcn(vx,vy,nv,xc,yc,art_ele,ntve,nbve,nbe,isbce,isonb,spherical);

     %---find out the cutout polygons of the transect 
      xline=[xp(1),xp(end)];
      yline=[yp(1),yp(end)];

     for itce=1:NN
          %calculate xpp and ypp which are the polygons of the TCE's
          [xpp{itce},ypp{itce},i_tce{itce},i_bce{itce}]=fvcom_tce_polygon(spherical,ne,ienode,ntve,nv,vx,vy, ...
                                                                    isonb,nbve,nbvt,xc,yc,ncv,niec,ntrg,itce);

          %calculate the intersections of the line with each polygon
          [xi{itce},yi{itce},segments{itce}]=polyxline(xpp{itce},ypp{itce},xline,yline); 

          %calculate the cutout polygons of each segment
          [polycuts{itce}]=fvcom_segment_cutout_polygon(tri,vx,vy,xije,yije,itce,...
                            xpp{itce},ypp{itce},i_bce{itce},i_tce{itce},xi{itce},yi{itce},segments{itce});
     end

     %
     %Calculate tce edges that are going to be used in evaluating fluxes at the cutout polygons
     %This will label tce edges on which fluxes will be calculated so that we 
     %do not have to calculate fluxes on tce edges that are irrelevant to the transect
     %This will speed up the calculation in  adv_s() function
     %
    
    TCE_edge_use=[]; 

    for itce=1:NN   %loop through all tce edges
        NSegs=length(polycuts{itce});  %number of segmetns of the tranect cut by itce'th polygon
	for iseg=1:NSegs   %loop all segments
	     if(~isempty(polycuts{itce}{iseg}.xpc)) %if there is a valid cutout polygon
                 for iedge=1:length(polycuts{itce}{iseg}.itce_c(:))  %loop through all the tce edges
                                                                               %that happen to be part of the cutout polygon
                     itce_edge=polycuts{itce}{iseg}.itce_c(iedge);   %get the tce edge index
  	             if(itce_edge~=0)  %valid tce_edge, if zero then it is a boundary element edge
		        TCE_edge_use=[TCE_edge_use itce_edge];
                     end
                 end
  	     end
        end
    end
 
    %sort in ascending order and remove duplication
    TCE_edge_use_uniq=unique(sort(TCE_edge_use,'ascend'));

    if(debug_plot) %debug plot all the active TCE eges that will have flux values calculated
	  hfig_tce_use=figure('position',[1,1, 800, 800],'visible','on'); hold on;
         	  %plot the mesh
 	          trimesh(tri,vx,vy,'color',[0 0 0]);
        	  %show the TCEs
	          line(xije',yije','color',[0 1 1],'linestyle','--');
	          %show the TCE vertices
	          plot(xije(:,2),yije(:,2),'k.');  %edge mid points
	          plot(xije(:,1),yije(:,1),'ko');  %element centers
	  	  %plot the transect
		  line(xt,yt,'color',[1 0 0]);
		  %plot the active TCE edges for inspection
	          for itce=1:length(TCE_edge_use_uniq(:))
	              itce_edge=TCE_edge_use_uniq(itce);
	              line(xije(itce_edge,:)',yije(itce_edge,:)','color',[0 1 0]);
	          end
        	  tttt_tmp=input('ctrl-c to quit, enter to continue'); 
          delete(hfig_tce_use); 
    end

%    TCE_edge_use_uniq=(1:1:ncv_i);  %test with all TCE edges

    %
    %Also collect the TCE's that are going to be used
    %

     TCE_use=[];

    for itce=1:NN
        NSegs=length(polycuts{itce});
        for iseg=1:NSegs
            if(~isempty(polycuts{itce}{iseg}.xpc))  %there is a valid cut
              TCE_use=[TCE_use itce];
            end
        end
    end
    
    %sort in ascending order and remove duplication
    TCE_use_uniq=unique(sort(TCE_use,'ascend'));

    if(debug_plot) %debug plot all the active TCE eges that will have flux values calculated
          hfig_tce_use=figure('position',[1,1, 800, 800],'visible','on'); hold on;
                  %plot the mesh
                  trimesh(tri,vx,vy,'color',[0 0 0]);
                  %show the TCEs
                  line(xije',yije','color',[0 1 1],'linestyle','--');
                  %show the TCE vertices
                  plot(xije(:,2),yije(:,2),'k.');  %edge mid points
                  plot(xije(:,1),yije(:,1),'ko');  %element centers
                  %plot the transect
                  line(xt,yt,'color',[1 0 0]);
                  %plot the active TCE edges for inspection
                  for itce_use=1:length(TCE_use_uniq(:))
                      itce=TCE_use_uniq(itce_use);
                      plot(vx(itce),vy(itce),'ro','markersize',10,'markerfacecolor',[1 0 0]);
                  end
                  tttt_tmp=input('ctrl-c to quit, enter to continue');
          delete(hfig_tce_use);
    end


     %run mean (to finish coding)
     %tdata2d_runmean=[];
     %tdata3d_runmean=[];

     %get variable information

     switch(i_par)
            case  1  %three components of u,v,w
                tdata2d_fvf.Nvar=9;
                tdata2d_fvf.varnames={'r','z','zdummy',...
                                       'flow_up','flow_up_fw','transport', ...      %piston velocity and transport
                                       'flow_prfd','flow_prfd_fw','tflux_prfd'}; %flow or flux desnity (per unit depth)
                tdata3d_fvf.Nvar=9;
                tdata3d_fvf.varnames={'x','y','z',                             ...  %piston velocity and transport
                                      'flow_up','flow_up_fw','transport',      ...  %flow or flux density (per unit depth)
                                      'flow_prfd','flow_prfd_fw','tflux_prfd', ...  %flow or flux density
                                      'depth'};
     end

    %
    %load history data
    %

    %Loop through all files one by one
    %save data into matlab file for statistics after the time loop is finished

    display(['Please wait while loading history output files']);

    it_all=0;    %time record counter
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
              par5=[];
              
              switch(i_par)
                    case  1  %three components of u,v,w
                         par1=hhis{'u'}(:);
                         par2=hhis{'v'}(:);
                         par3=hhis{'ww'}(:);

                         parname='velocity';
                         parunit='(m/s)';
                         par4=hhis{tracername}(:); 

			 par5=hhis{'w'}(:); 

                         if(prod(size(par4))==0)
                           error(['oops, not able to read variable ' ...
                                  tracername ' in file ' hfile, ' or it is empty']); 
                         end
                          
                         

                         %
                         %get tracer units here
                         %tracerunit = 
                         %

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
		    case  5
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

               tmp_matfile=[pltdirname './'  'fvcom_hisdata_' num2str(nf,'%04d') '.mat'];
               save(tmp_matfile,   'par1', ...  %1st variable
                                   'par2', ...  %2nd variable
                                   'par3', ...  %3rd variable
                                   'par4', ...  %4th variable 
                                   'par5', ...  %5th variable
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
               tmp_matfile=[pltdirname './'  'fvcom_hisdata_' num2str(nf,'%04d') '.mat'];
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

   end

   %
   %Now only read matlab format data 
   %

   %find total number of time records
   %also record the file numbers for for each time record
   htime_global=[]; 
   filenum_global=[];
   it_local_global=[];

   for nf=his_fn_start:1:his_fn_end
          if(read_his_matlab==true) 
            hfile=[his_dir '/' hfile_prefix  num2str(nf,'%04d') hfile_ext];
          else
            hfile=[pltdirname './' 'fvcom_hisdata_' num2str(nf,'%04d') '.mat'];
          end
          if(~exist(hfile,'file'))
            error(['Oops, cannot find file ' hfile ]);
          end
          load(hfile);
          his_NT=length(htime);   %number of time records in this file

         %give references to find file numbers and record numbers in file for any given record in htime_all
          htime_global=[htime_global htime(:)'];
          filenum_global=[filenum_global nf*ones([1,his_NT])];    %file number for each time record in htime_all
          it_local_global=[it_local_global (1:his_NT)];           %local record number in each file
     
   end
 
   NT_global=length(htime_global); 

   it_all=0;
   for it_global=1:t_int:NT_global

       it_all=it_all+1; 

       %find the history file to load
       nf=filenum_global(it_global);

       %read matlab file of the data saved above
            if(it_all==1)
              hfile_old=[];
            else
              hfile_old=hfile;
            end

            if(read_his_matlab)
                hfile=[his_dir '/' hfile_prefix  num2str(nf,'%04d') hfile_ext];
            else
                hfile=[pltdirname './' 'fvcom_hisdata_' num2str(nf,'%04d') '.mat'];
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

            %plot

            if(yesno_matlabfig)
              hfig=figure('visible','off','position',[1,1,1000,800]); hold on;
              set(gcf,'renderer','zbuffer');
            end

               %calculate flux 

	       %----debug----
		    nhe=0;  %number of halo elements
                    nhn=0;  %number of halo nodes

                    m=size(vx,1);  %number of nodes
                    n=size(nv,1);;  %number of elements
                    mt=m+nhn;  %number of total nodes
                    nt=n+nhe;  %number of total elements

                    northpole=0;               %whether to use northpole
		    mpdata=0;                  %whether to use mpdata
                    semi_implicit=0;           %whether to use semi_implicit 
                    wet_dry=0;                 %whether to use wetting and drying
                    multiprocessor=0;          %where to use multiple processors
                    horzmix='closure';         %horzontal mixing method
                    one_d_model=0;             %whether it is one-d model
                    point_st_type='calculated';%point source of t and s 
                    inflow_type='node';        %inflow type
                    dti=1.0;                   %internal model time step 

                    kb=KB;
                    KBM1=KB-1; 
                    kbm1=KBM1;
                    dz1=zeros([nt,kb]);        %vertical sigma coordinate step size 
                                               %for elements for each layer
                    for k=1:kbm1
	                    dz1(1:nt,k)=  z(k)-z(k+1) ;
		    end
                    dz=dz1;

                    dz_n=zeros([mt,kb]);
                    for k=1:kbm1
                        dz_n(1:mt,k)=  z(k)-z(k+1) ;
                    end

                    if(it_all>=2)
                        zeta_n_1=zeta_n;    %zeta at node for previous time step
                    else
                        zeta_n_1=zeta(it,:)';
                    end
                    zeta_n=zeta(it,:)';     %zeta at node for current step


                    zeta_e=nan*zeros(size(xc));    %element center surface elevation (current step)
                    zeta_e_1=nan*zeros(size(xc));  %" "                              (previous step)

                     

                    for ie=1:length(xc(:))
                        zeta_e(ie)=mean(zeta(it,nv(ie,:)));         %take mean of  nodal values (current step)
                        zeta_e_1(ie)=mean(zeta_n_1(nv(ie,:)));   %                           (previous step)
                    end

		    d=h_n+zeta_n;            %total depth at node (current step)
                    d1=hc+zeta_e;            %total depth at element (current step)

                    dt=h_n+zeta_n_1;         %total depth at node at previous step
                    dt1=hc+zeta_e_1;         %total depth at element at previous step
                    
		    dtfa=dt;                 %total depth at node after adjustment

                   % dt_time=3600.0;         %time interval (in seconds)
                    if(it_all>=2)
                        dt_time=(htime_all(it_all)-htime_all(it_all-1))*86400.0;  %days to seconds
                    else
                        dt_time=inf; 
                    end

                    ibfw=0;                    %bottom flux
                    bfwdis3=[];                %distribution of bottom flux
                    numqbc=0;                  %number of point sources
                    sdis=[];                   %point source distribution
                    node_northarea=[];	       %???

                    s1=zeros([mt,kb]);         %tracer concentration at nodes
                    s1_1=zeros([mt,kb]);       %tracer concentration at nodes previous step
                    if(it_all>=2)
                       s1_1=s1;           
                    else
                       s1_1(:,1:KBM1)=squeeze(par4(it,1:KBM1,:))';  %get the tracer concentration
                    end

                    s1(:,1:KBM1)=squeeze(par4(it,1:KBM1,:))'; 

	            %--debug---set salinity to 1.0
%                      s1(:,:)=32.0*ones(size(s1));
%                      s1_1(:,:)=s1(:,:); 
                    %----------

                    smean1=zeros([mt,kb]);     %vertical mean tracer concentration
                    smean1=s1;                 %these is only used for calculating diffusion of 
                                               %baroclinic salinity, i.e.  diffusion of  
                                               %(tracer - depth average of tracer)

                    viscofh=zeros([mt,kb]);    %horizontal viscosity  %turbulent diffusivity
                    wts=zeros([mt,kb]);        %vertical velocity

                    u=zeros([nt,kb]);          %horizontal velocity u and v
                    v=zeros([nt,kb]);         
                    u(:,1:KBM1)=squeeze(par1(it,1:KBM1,:))';  %u velocity (m/sec) at time step it in file

                    %---debug----
                       %reverse u direction in the bottom layers
                       %u(:,1:4)   = 1.0*ones(size(u(:,1:4)   ));
                       %u(:,5:KBM1)=-1.0*ones(size(u(:,5:KBM1)));
                    %------------

                    v(:,1:KBM1)=squeeze(par2(it,1:KBM1,:))';  %v velocity (m/sec) at time step it in file

		    %----debug---%set v to zero
                    %v(:,1:KBM1)=zeros(size(v(:,1:KBM1)));
                    %------------

                    horcon=0;                  %parameter in casename_run.dat that controls
                                               %horizontal diffusivity. HORCON  = 2.000E-1
                                               %final diffusivity in code = horcon*viscofh
                    art =art_ele;
                    art1=art_tce;
                    art2=art_can;
	       %-------------

               %-------------hardwired temp variables-------
               nmqbc=0;                     %number of rivers
               qdis=[];                     %river flux at current time
               vqdis=[];                    %!!DISCHARGE VERTICAL DISTRIBUTION
               mean_flow=0;                 %!!flag to have mean flow calculated
               nmfcell=0;                   %mean flow cell
               node_mfcell=[];              %nodes of mean flow
               ibfw=0;                      %number of nodes with ground water flow
               rofvros=1000.0/1023.0;       %!RATIO OF THE DENSITY OF FRESH AND SEA WATER 1000./1023.
               qevap3=zeros([m,1]);         %surface evaporation for internal mode
               qprec3=zeros([m,1]);         %surface precipitation for internal mode
               bfwdis=[];                   %!GROUNDWATER FLUX AT CURRENT TIME
               node_bfw=[];                 %node numbers of bottom fresh water
               mfqdis=[];                   %
               rdismf=[];                   %ratio:  RDISMF(I,1)=ART1(I1)/(ART1(I1)+ART1(I2))
               mfdist=[];                   %
               iswetnt=ones([m,1]);         %NODE POROSITY AT NODES FOR TIME N-1 INTERNAL
               iswetn =ones([m,1]);         %NODE POROSITY AT NODES FOR TIME N
               iswetc =ones([m,1]);         %CELL POROSITY AT CELLS FOR TIME N
               uf=u;			    %
               vf=v;                        %
               ifceta=0     ;               %IMPLICIT FACTOR SELECTED
               %-------------------------------------------

               %----
               vqdist=[];
               rdisq=[];
               elf=zeta_n;
               %----

               %calculate pseudo vertical velocity
               [w,wts]=fvcom_wpseudo(ncv,niec,ntrg,dt1,dz1,u,v,dltye,dltxe,d1,uf,vf,...
                               spherical,northpole,wet_dry, semi_implicit,one_d_model,isonb,numqbc,...
                               qdis,vqdist,...
                               nmfcell,node_mfcell,node_bfw,rdisq,mfqdis,rdismf,mfdist,qevap3,qprec3,ibfw,...
                               iswetnt,iswetn,iswetc,nv,ifceta,kb,kbm1,inflow_type,mean_flow,n,m,rofvros,art1,dz,d,dt,dti,elf,h);

               %save w to w_all for plotting

               if(it_all==1)
                 w_all=[];        %w calculated by matlab here      %on levels 
                 w_fvcom_all=[];  %w calculated by fvcom internally %on levels
                 ww_all=[];       %ww calculated by fvcom   %on layers
		 zeta_n_all=[];   %zeta calculated by fvcom
                 u_all=[];
                 v_all=[];
               end

               w_all=[w_all;reshape(w,[1,n,KB])]; 
               w_fvcom_all =[w_fvcom_all; reshape(squeeze(par5(it,:,:))',[1,n,KB])];
               ww_all=[ww_all; reshape(squeeze(par3(it,:,:))',[1,n,KBM1])];
               u_all =[ u_all; reshape(squeeze(par1(it,:,:))',[1,n,KBM1])];
               v_all =[ v_all; reshape(squeeze(par2(it,:,:))',[1,n,KBM1])];

               zeta_n_all=[zeta_n_all; zeta(it,:)]; 

                %----plot w------against model simulation results

% iline=find(y_e< 1600 & y_e>1400);
% xiline=x_e(iline);
% yiline=y_e(iline);
% xiline2D=repmat(xiline,[1,KB]);
% zetailine_1=zeta_e_1(iline);
% ziline=zeros(100,11);
% for k=1:KB
%       ziline(:,k)=zetailine+z(k)*(zetailine+hiline);
% end
% uiline=u(iline,:);
% viline=v(iline,:);
% wiline=w(iline,:);
% ww=squeeze(par3(1,:,:))';
% wwiline=ww(iline,:);
% zetailine_1=zeta_e_1(iline);
% zeta_e_all=zeros([72,800]);
% for itt=1:72
%  for ie=1:800
% zeta_e_all(itt,ie)=mean(zeta(itt,e2n(ie,2:4)));
% end
%end

                %------------------------------------------------

end

tttt=input('ctrl-c to debug'); 

for it_global=1:t_int:NT_global

               %calculate fluxes using the cutout polygons, for this particular vertical layer
               [f_fw_tce_edge, xflux,xflux_adv,xflux_tce_edge,xflow_tce_edge]=adv_s(spherical,northpole,...
                                                                           mpdata,semi_implicit,wet_dry,...
                                                            multiprocessor,horzmix,ncv,ntrg,dt1,dz1,u,v,...
                                                          dltxe,dltye,ntsn,nbsn,s1,smean1,vx,vy,art_can,...
                                                              viscofh,ncv_i,niec,xije,yije,horcon,iobcn,...
                                                               i_obc_n,one_d_model,wts,dz,isonb,art_tce,...
                                                             point_st_type,inflow_type,numqbc,sdis,ibfw,...
                                                               bfwdis3,dti,dt,dtfa,node_northarea,kb,nv,...
                                                               TCE_edge_use_uniq,TCE_use_uniq);
	       %---debug----
               if(debug_plot)
	          figure;
                  trimesh(tri,vx,vy,'color',[0 0 0]); hold on
                  %show the TCEs
                  line(xije',yije','color',[0 1 1],'linestyle','--');
                  %show the TCE vertices
                  plot(xije(:,2),yije(:,2),'k.');  %edge mid points
                  plot(xije(:,1),yije(:,1),'ko');  %element centers
                  %for i=1:100
                  %    text(vx(i),vy(i),num2str(i),'color',[1 0 0]);
                  %end
                  axis([-46000,1000,-500,2800]);
               end
               %-----------

               %collect flux through the transect
               for itce=1:m
                   NSegs=length(polycuts{itce});
                   %calculate flux on each segment
                   for iseg=1:NSegs
                       cflow{itce}{iseg}=0.0;                  %total volume flux  (m^3/sec)
                       cflow_prf{itce}{iseg}=0.0*zeros([1,KBM1]);
                       cflow_cutcord{itce}{iseg}=0.0;
                       cflow_cutcord_prf{itce}{iseg}=0.0*zeros([1,KBM1]);  %vertical profile

                       cflow_fw{itce}{iseg}=0.0;               %fresh water volume flux (m^3/sec)
                       cflow_fw_prf{itce}{iseg}=0.0*zeros([1,KBM1]);
                       cflow_cutcord_fw{itce}{iseg}=0.0;
                       cflow_cutcord_fw_prf{itce}{iseg}=0.0*zeros([1,KBM1]);
                                                               %salt water volume flux = 
                                                               %   total volume flux - fresh water volume  flux 

                       tflux{itce}{iseg}=0.0;                  %tracer flux  [c]*m^3/sec 
                       tflux_prf{itce}{iseg}=0.0*zeros([1,KBM1]);
                       tflux_cutcord{itce}{iseg}=0.0;
                       tflux_cutcord_prf{itce}{iseg}=0.0*zeros([1,KBM1]);

  		       if(~isempty(polycuts{itce}{iseg}.xpc))  %there is a valid cut

		  if(debug_plot)
			   line(xpp{itce},ypp{itce},'color','m')
			   np=length(xpp{itce}(:)); 

			 %label the tce polygon
			 %  for i=2:np-1
			 %     textangle=atan2(ypp{itce}(i)-vy(itce),xpp{itce}(i)-vx(itce));
			 %     text(xpp{itce}(i)+(max(xpp{itce}(:))-min(xpp{itce}(:)))/5.0*cos(textangle) ...
			 %         ,ypp{itce}(i)+(max(xpp{itce}(:))-min(xpp{itce}(:)))/5.0*sin(textangle) ...
			 %         ,num2str(i),'color',[1 0 0]);
			 %  end
			 %  i=np ; %last point goes back to the starting point (avoid overlap of text)
			 %  textangle=atan2(ypp{itce}(i)-vy(itce),xpp{itce}(i)-vx(itce));
			 %  text(xpp{itce}(i)+(max(xpp{itce}(:))-min(xpp{itce}(:)))*1/5.0*cos(textangle) ...
			 %      ,ypp{itce}(i)+(max(xpp{itce}(:))-min(xpp{itce}(:)))*1/5.0*sin(textangle) ...
			 %      ,['1(' num2str(i) ')'],'color',[1 0 0]);

			   for i=1:np-1
			    itce_edge=i_tce{itce}(i);
			    if(itce_edge>0)
				    ia=niec(itce_edge,1); %get the spoke of the TCE
				    ib=niec(itce_edge,2);
				    %plot the spoke
			     %       line([vx(ia),vx(ib)],[vy(ia),vy(ib)],'color','g');
				    text_x=(xpp{itce}(i)+xpp{itce}(i+1))/2.0;
				    text_y=(ypp{itce}(i)+ypp{itce}(i+1))/2.0;
				    textangle=atan2(text_y-vy(itce),text_x-vx(itce));
				    %label the tce_edge
			     %       text( text_x+(max(xpp{itce}(:))-min(xpp{itce}(:)))*1/8.0*cos(textangle),...
			     %             text_y+(max(xpp{itce}(:))-min(xpp{itce}(:)))*1/8.0*sin(textangle),...
			     %             num2str(itce_edge),'color',[0 1 1]);
			    end
			    ibce_edge=i_bce{itce}(i);
			    if(ibce_edge>0)
				    text_x=(xpp{itce}(i)+xpp{itce}(i+1))/2.0;
				    text_y=(ypp{itce}(i)+ypp{itce}(i+1))/2.0;
				    textangle=atan2(text_y-mean(ypp{itce}(:)),text_x-mean(xpp{itce}(:)));
				    %label the boundary edge
			      %      text( text_x+(max(xpp{itce}(:))-min(xpp{itce}(:)))*1/8.0*cos(textangle),...
			      %            text_y+(max(xpp{itce}(:))-min(xpp{itce}(:)))*1/8.0*sin(textangle),...
			      %            num2str(ibce_edge),'color',[1 0 1]);
			    end
			   end
                  end
			   %loop through the cutout polygon edges, collect the fluxes on all the edges
			   xpc=polycuts{itce}{iseg}.xpc;
			   ypc=polycuts{itce}{iseg}.ypc;

                           if(debug_plot)
				   line(xpc,ypc,'color','g');  %show the cutout polygon
                           end

			   for iedge=1:length(polycuts{itce}{iseg}.itce_c(:))  %loop through all the tce edges
				  					       %that happen to be part of the cutout polygon
			       itce_edge=polycuts{itce}{iseg}.itce_c(iedge);   %get the tce edge index
			       if(itce_edge>0) %if 0, then not a tce edge

				       %get the flux on the tce edge
				       if(iedge==1)  %first edge on the cutout polygon, partial tce edge
						%get the proportion due to the cut
						dist_partial= (xpc(iedge)-xpc(iedge+1))^2 + ...
							      (ypc(iedge)-ypc(iedge+1))^2 ;
						dist_tce    = (xije(itce_edge,1)-xije(itce_edge,2) )^2 + ...
							      (yije(itce_edge,1)-yije(itce_edge,2) )^2;
						dist_frac=sqrt(dist_partial/dist_tce);

				       elseif(iedge==length(polycuts{itce}{iseg}.itce_c(:))) %last edge on the cutout polygon
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
				       %set flux as being positive leaving the tce polygon xpp,ypp
				       if(niec(itce_edge,1)==itce)  %xflux_tce_edge is into polygon xpp,ypp, flip the sign
                                             cflow{itce}{iseg}= cflow{itce}{iseg}              ...
                                                               - sum(xflow_tce_edge(itce_edge,:)).*dist_frac;
                                             cflow_prf{itce}{iseg}= cflow_prf{itce}{iseg}      ...
                                                               - (xflow_tce_edge(itce_edge,:)).*dist_frac;

                                             cflow_fw{itce}{iseg}=cflow_fw{itce}{iseg}           ...
                                                               - sum(xflow_tce_edge(itce_edge,:) ...
                                                                    .*f_fw_tce_edge(itce_edge,:)).*dist_frac;
                                             cflow_fw_prf{itce}{iseg}=cflow_fw_prf{itce}{iseg}   ...
                                                               - (xflow_tce_edge(itce_edge,:)    ...
                                                                    .*f_fw_tce_edge(itce_edge,:)).*dist_frac;

					     tflux{itce}{iseg}= tflux{itce}{iseg}              ...
                                                               - sum(xflux_tce_edge(itce_edge,:)).*dist_frac;
                                             tflux_prf{itce}{iseg}= tflux_prf{itce}{iseg}      ...
                                                               -    (xflux_tce_edge(itce_edge,:)).*dist_frac;
				       else                         %xflux_tce_edge is entering polygon xpp,ypp
					     cflow{itce}{iseg}= cflow{itce}{iseg}              ...
                                                               + sum(xflow_tce_edge(itce_edge,:)).*dist_frac;
                                             cflow_prf{itce}{iseg}= cflow_prf{itce}{iseg}      ...
                                                               + (xflow_tce_edge(itce_edge,:)).*dist_frac;

                                             cflow_fw{itce}{iseg}=cflow_fw{itce}{iseg}           ...
                                                               + sum(xflow_tce_edge(itce_edge,:) ...
                                                                    .*f_fw_tce_edge(itce_edge,:)).*dist_frac;
                                             cflow_fw_prf{itce}{iseg}=cflow_fw_prf{itce}{iseg}   ...
                                                               + (xflow_tce_edge(itce_edge,:)    ...
                                                                    .*f_fw_tce_edge(itce_edge,:)).*dist_frac;

                                             tflux{itce}{iseg}= tflux{itce}{iseg}      ...
                                                               + sum(xflux_tce_edge(itce_edge,:)).*dist_frac;
                                             tflux_prf{itce}{iseg}= tflux_prf{itce}{iseg}      ...
                                                                  + (xflux_tce_edge(itce_edge,:)).*dist_frac;
				       end

				       %show the flux on the edges of the cutout polygon
                                    if(debug_plot)
				       %show the flux on the tce edge
				       text_x=(xpc(iedge)+xpc(iedge+1))/2.0;
				       text_y=(ypc(iedge)+ypc(iedge+1))/2.0;
				       textangle=atan2(text_y-vy(itce),text_x-vx(itce));
				       %label flux on the cut out polgyon (in the polygon side)
				       if(niec(itce_edge,1)==itce)
					 % text( text_x-(max(xpc(:))-min(xpc(:)))*1/8.0*cos(textangle),...
					 %	text_y-(max(xpc(:))-min(xpc(:)))*1/8.0*sin(textangle),...
					 %	num2str(-sum(xflux_tce_edge(itce_edge,:)).*dist_frac),...
					 %	'color',[0 1 1]);

                                          text( text_x-(max(xpc(:))-min(xpc(:)))*1/8.0*cos(textangle),...
                                               text_y-(max(xpc(:))-min(xpc(:)))*1/8.0*sin(textangle),...
                                               num2str(-(xflux_tce_edge(itce_edge,8)).*dist_frac),...
                                               'color',[0 1 1]);
				       else
					 % text( text_x-(max(xpc(:))-min(xpc(:)))*1/8.0*cos(textangle),...
					 %	text_y-(max(xpc(:))-min(xpc(:)))*1/8.0*sin(textangle),...
					 %	num2str(+sum(xflux_tce_edge(itce_edge,:)).*dist_frac),...
					 %	'color',[0 1 1]);
                                          text( text_x-(max(xpc(:))-min(xpc(:)))*1/8.0*cos(textangle),...
                                                text_y-(max(xpc(:))-min(xpc(:)))*1/8.0*sin(textangle),...
                                               num2str(+(xflux_tce_edge(itce_edge,8)).*dist_frac),...
                                               'color',[0 1 1]);
				       end

%                                       tmp_debug=input('ctrl-c to stop 0');

                                    end

			       else %must be a boundary edge

				    ibce_edge=polycuts{itce}{iseg}.ibce_c(iedge);
				    if(ibce_edge==0)  %something wrong, must be either a tce edge or a boundary edge
					error(['oops edge type unidentified with the cutout polygon']);
				    else %normal boundary edges

					 %collect boundary flux if an open boundary, otherwise,
                                         % flux through this edge is zero

					 %check if the boundary is an open boundary

					 %find the element that ibce_edge is part of
				    end
			       end
			   end

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
			 %   where V_total == V_tot  is total volume of the bulk water at salinity S, 
                         %                              with S between 0 and Smax
			 %         V_Smax            is volume of bulk water with salinity S is anti-diluted
                         %                              to have salinity Smax
			 %         m_S               is mass of salt in bulk water volume V_tot with salinity S
			 %         m_fw              is mass of fresh water in bulk water volume v_tot with salinty S
			 %         V_fw              is volume of fresh water needed to mix with V_Smax water of
			 %                              salinity Smax to produce
			 %                              V_total volume of bulk water with salinity S
			 %         rho_bulk(S_prime) is density of bulk water within given salinity value S_prime
			 %         rho_bulk(Smax)    is density of seawater with salinity equal to Smax
			 %         rho_bulk(S)       is density of seawater with salinity equal to S
			 %         rho_fw = rho_bulk(S_prime=0) is density of fresh water (salinity eqtual to zero)
			 %         f_fw              is fraction of fresh water needed to generate seawater 
                         %                              of salinity S by mixing
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

			      %calculate fresh water flux (m^3/sec)
				  %get rho_bulk(S)    using equation of state , where S=Tp4
				  %get rho_bulk(Smax) using equation of state,  where Smax is set to 34
	 		 %
			 % Equation of State: Reference: Cushman-Roisin B.,Introduction to Geophysical Fluid Dynamics,
			 %                               page 36, Chapter 3,
			 %                               ISBN 0-13-353301-8, Prentice-Hall Inc., 1994
			 % rho_bulk = rho_0 * [ 1- alpha * (T-T0) + beta *( S-S0) ]  (eqn 7)
			 % where alpha = 1.7E-4  (1/K)
			 %       beta  = 7.6E-4
			 %       T0    = 10 degC
			 %       S0    = 35 (ppt)
			 %       rho0  = 1028 (kg/m^3)
                         %

			 %get the rate of change in the cutout polygon
			   %current time step
				%get the volume of the full polygon xpp,ypp
				   art_tce(itce);
				%get the volume of the cutout polygon xpc,ypc
				   area_cutout=polygon_area([xpc xpc(1)],[ypc ypc(1)]); 
                                                     %close polygon by going to starting point
				%get the conentration in the box for all layers
				   c_tce=s1(itce,1:kbm1);
				%calculate the total mass in the cutout polygon
				   mass_tce=(c_tce.*area_cutout     ...
					      .*(h_n(itce)+zeta_n(itce)).*dz_n(itce,1:kbm1));
                                   volume_tce=(area_cutout     ...
                                              .*(h_n(itce)+zeta_n(itce)).*dz_n(itce,1:kbm1));
                                   %calculate fraction of fresh water volume in the cutout polygon
                                   rho_bulk_S    = 1028 * (1 - 1.7E-4 * (10 -10 ) + 7.6E-4 *(c_tce-35));
                                          %eqn 7 assuming S=c_tce  (mixed)
                                   rho_bulk_Smax = 1028 * (1 - 1.7E-4 * (10 -10 ) + 7.6E-4 *(34 -35));
                                          %eqn 7 assuming S=Smax=34 ppt (salt)
                                   rho_fw        = 1028 * (1 - 1.7E-4 * (10 -10 ) + 7.6E-4 *(0  -35));
                                          %eqn 7 assuming S= 0 ppt  (fresh)
                                    %get f_fw using  eqn(6) above
                                   f_fw= ( rho_bulk_Smax - rho_bulk_S )./(rho_bulk_Smax - rho_fw) ;

			   %next time step
				%get the volume of the full polygon xpp,ypp
				   art_tce(itce);
				%get the volume of the cutout polygon xpc,ypc
				   area_cutout=polygon_area([xpc xpc(1)],[ypc ypc(1)]); 
                                                      %close polygon by going to starting
						      %point
				%get the conentration in the box
				   c_tce=s1_1(itce,1:kbm1);
				%calculate the total mass
				   mass_tce_1=(c_tce.*area_cutout     ...
					      .*(h_n(itce)+zeta_n_1(itce)).*dz_n(itce,1:kbm1)); 
                                   volume_tce_1=(area_cutout     ...
                                              .*(h_n(itce)+zeta_n_1(itce)).*dz_n(itce,1:kbm1));
                                   %calculate fraction of fresh water volume in the cutout polygon
                                    rho_bulk_S    = 1028 * (1 - 1.7E-4 * (10 -10 ) + 7.6E-4 *(c_tce-35));
                                          %eqn 7 assuming S=c_tce  (mixed)
                                    rho_bulk_Smax = 1028 * (1 - 1.7E-4 * (10 -10 ) + 7.6E-4 *(34 -35));
                                          %eqn 7 assuming S=Smax=34 ppt (salt)
                                    rho_fw        = 1028 * (1 - 1.7E-4 * (10 -10 ) + 7.6E-4 *(0  -35));
                                          %eqn 7 assuming S= 0 ppt  (fresh)
                                    %get f_fw using  eqn 6
                                    f_fw_1= ( rho_bulk_Smax - rho_bulk_S )./(rho_bulk_Smax - rho_fw) ;

			     %get the flux on the cutting cord using mass balance

                             %tracer flux:
			     %d(mass)/dt =  flux_in -flux_out
			     %===> flux_out =  flux_in     - d(mass)/dt
			     %              = -tflux{itce}{iseg} - d(mass)/dt
			      tflux_cutcord{itce}{iseg}= -tflux{itce}{iseg} ...
                                                         -sum(mass_tce_1-mass_tce)/dt_time; %([c]*m^3/sec)
                              tflux_cutcord_prf{itce}{iseg}= -tflux_prf{itce}{iseg} ...
                                                         -(mass_tce_1-mass_tce)./dt_time; %([c]*m^3/sec)

                              %total volume flux
                              %d(Volume)/dt = flux_in - flux_out
                              %==> flux_out = flux_in - d(volume)/dt
                              cflow_cutcord{itce}{iseg}= -cflow{itce}{iseg} ...
                                                         -sum(volume_tce_1-volume_tce)/dt_time; %(m^3/sec)
                              cflow_cutcord_prf{itce}{iseg}= -cflow_prf{itce}{iseg} ...
                                                         -(volume_tce_1-volume_tce)/dt_time; %(m^3/sec)

                              %fresh water volume flux
                              %d(Volume*f_fw)/dt = flux_in_fw -flux_out_fw
                              %==> flux_out_fw = flux_in_fw - d(Volume*f_fw)/dt
                              cflow_cutcord_fw{itce}{iseg} = -cflow_fw{itce}{iseg} ...
                                         -sum(volume_tce_1.*f_fw_1-volume_tce.*f_fw)/dt_time;   %(m^3/sec)
                              cflow_cutcord_fw_prf{itce}{iseg} = -cflow_fw_prf{itce}{iseg} ...
                                         -(volume_tce_1.*f_fw_1-volume_tce.*f_fw)/dt_time;   %(m^3/sec)

                             if(debug_plot)
			     %show the cutting segment
			      line(segments{itce}{iseg}.x,segments{itce}{iseg}.y,'color','r');
			     %show the flux on the cutting segment
			      text_x=mean(segments{itce}{iseg}.x(:));
			      text_y=mean(segments{itce}{iseg}.y(:));
			      textangle=atan2(text_y-vy(itce),text_x-vx(itce));
			      text( text_x-(max(xpc(:))-min(xpc(:)))*1/8.0*cos(textangle),...
				    text_y-(max(xpc(:))-min(xpc(:)))*1/8.0*sin(textangle),...
				    num2str(tflux_cutcord{itce}{iseg}),                   ...
				    'color',[1 0 0]);
                             end

		       else %not cutting, flux for this segment is zero
			     cflow{itce}{iseg}=0.0;
			     cflow_cutcord{itce}{iseg}=0.0;
                             cflow_cutcord_fw{itce}{iseg}=0.0; 

                             tflux{itce}{iseg}=0.0;
                             tflux_cutcord{itce}{iseg}=0.0;
		       end %end of valid cut
                   end %iseg loop
               end %itce loop

%               tmp_debug=input('ctrl-c to stop 1');

               %summarize the flux ( sum up all the segments), re-order them along the line

               cflow_use=[];         %the flow of each segment on the transect
               cflow_prf_use=[];     %flow with vertical profile
               cflow_fw_use=[];      %fresh water flow
               cflow_fw_prf_use=[];  %fresh water flow with vertical profile
               tflux_use=[];         %tracer flux
               tflux_prf_use=[];     %tracer flux with vertical profile

               trasect_seg_list=[];  %keep the list of segments that went in 
                                     %   so that we can order them along the transect line
                                     %   and be able to make cross sectinal profile plots later
               iseg_tot=0;
               for itce=1:m
                   NSegs=length(polycuts{itce});
                   %calculate flux on each segment
                   for iseg=1:NSegs
                       if(~isempty(polycuts{itce}{iseg}.xpc))  %there is a valid cut
                           iseg_tot=iseg_tot+1;
                           %water volume flux  (m^3/sec) (sum up all the segments)
                           cflow_use=[cflow_use; cflow_cutcord{itce}{iseg}];
                           cflow_prf_use=[cflow_prf_use; cflow_cutcord_prf{itce}{iseg}];

                           %fresh water volume flux (m^3/sec)
                           cflow_fw_use=[cflow_fw_use; cflow_cutcord_fw{itce}{iseg}];
                           cflow_fw_prf_use=[cflow_fw_prf_use; cflow_cutcord_fw_prf{itce}{iseg}];

                           %total tracer flux ([c]*m^3/sec)
                           tflux_use=[tflux_use; tflux_cutcord{itce}{iseg}];
                           tflux_prf_use=[tflux_prf_use; tflux_cutcord_prf{itce}{iseg}];

                           %also note down the starting and ending point of each segment
                           transect_seg_list{iseg_tot}.x=segments{itce}{iseg}.x ;
                                        %stating and ending point x coord of the transect
                           transect_seg_list{iseg_tot}.y=segments{itce}{iseg}.y ; 
                                        %stating and ending point y coord of the transect

                           %also note down the surface and depth of the starting and ending points of each segment
                           transect_seg_list{iseg_tot}.zeta=nan*zeros([2,1]); 
                           transect_seg_list{iseg_tot}.h=nan*zeros([2,1]);

                           %loop through the cutout polygon edges, collect the fluxes on all the edges
                           xpc=polycuts{itce}{iseg}.xpc;
                           ypc=polycuts{itce}{iseg}.ypc;

			   if(debug_plot)
	                           line(xpc,ypc,'color','g');  %show the cutout polygon
                           end

                           for iedge=1:length(polycuts{itce}{iseg}.itce_c(:))  %loop through all the tce edges
                                                                               %that happen to be part of the cutout
                                                                               %polygon
                               itce_edge=polycuts{itce}{iseg}.itce_c(iedge);   %get the tce edge index
                               if(itce_edge>0) %if 0, the not a tce edge
                                       if(iedge==1)  %first edge on the cutout polygon, partial tce edge
					%depth at the start point of the tce edge
					 ie_tce=ntrg(itce_edge); %element number of the tce edge

                                         %depth at start point of tce edge, which is element center
					 htce_start= mean(   h_n(nv(ie_tce,:)));
					 ztce_start= mean(zeta_n(nv(ie_tce,:)));  %starting point is element center
  					   			                  %average using nodes of the element

					%depth at the end point of the tce edge (mid point of spoke) 
					 ia=niec(itce_edge,1);  %spoke node number 1
					 ib=niec(itce_edge,2);  %spoke node number 2
                                         htce_end  = (h_n(ia)+h_n(ib))/2.0;      %depth and zeta at mid of spoke
					 ztce_end  = (zeta_n(ia)+zeta_n(ib))/2.0;

                                        %distance of interection point to start and ending points of the tce edge
					 dist_to_start=sqrt((xpc(iedge)-xije(itce_edge,1))^2 + ...
							    (ypc(iedge)-yije(itce_edge,1))^2); 
					 dist_to_end  =sqrt((xpc(iedge)-xije(itce_edge,2))^2 + ...
							       (ypc(iedge)-yije(itce_edge,2))^2);

                                        %interpolate elevation based on distance
					 transect_seg_list{iseg_tot}.zeta(2) = ( ztce_start*dist_to_end  + ...
					    				         ztce_end  *dist_to_start) ...
							/(dist_to_start+dist_to_end) ; 
                                         
                                         transect_seg_list{iseg_tot}.h(2)    = ( htce_start*dist_to_end  + ...
                                                                                 htce_end  *dist_to_start) ...
                                                        /(dist_to_start+dist_to_end) ;

                                       elseif(iedge==length(polycuts{itce}{iseg}.itce_c(:))) 
                                                                 %last edge on the cutout polygon
                                                                 %partial tce edge
					%depth at the start point of the tce edge
					 ie_tce=ntrg(itce_edge); %element number of the tce edge

                                        %depth at start point of tce edge, which is element center
                                         htce_start= mean(   h_n(nv(ie_tce,:)));
                                         ztce_start= mean(zeta_n(nv(ie_tce,:)));  %starting point is element center
                                                                                  %average using nodes of the
					%depth at the end point of the tce edge (mid point of element edge)
					 ia=niec(itce_edge,1);  %spoke node number 1
					 ib=niec(itce_edge,2);  %spoke node number 2

                                         htce_end  = (h_n(ia)   +    h_n(ib))/2.0; 
					 ztce_end  = (zeta_n(ia)+ zeta_n(ib))/2.0;  %depth at the midpoint of the spoke

					 dist_to_start=sqrt((xpc(iedge)-xije(itce_edge,1))^2 + ...
							    (ypc(iedge)-yije(itce_edge,1))^2);
					 dist_to_end  =sqrt((xpc(iedge)-xije(itce_edge,2))^2 + ...
	 					            (ypc(iedge)-yije(itce_edge,2))^2);

					 transect_seg_list{iseg_tot}.zeta(1) = ( ztce_start*dist_to_end  + ...
									         ztce_end  *dist_to_start) ...
							/(dist_to_start+dist_to_end) ; 
                                         transect_seg_list{iseg_tot}.h(1) = ( htce_start*dist_to_end  + ...
                                                                              htce_end  *dist_to_start) ...
                                                        /(dist_to_start+dist_to_end) ;

                                       else     %normal full tce edges,do nothing
                                         %do nothing   
                                       end
                               else %must be a boundary edge
                                    ibce_edge=polycuts{itce}{iseg}.ibce_c(iedge);
                                    if(ibce_edge==0)  %something wrong, must be either a tce edge or a boundary edge
                                        error(['oops edge type unidentified with the cutout polygon']);
                                    else %normal boundary edges
					if(iedge==1)  %first edge of the cutout polygon
                                           inode1=ienode(ibce_edge,1);  %one end of the bce edge
                                           inode2=ienode(ibce_edge,2);  %other end of the bce edge

                                           znode_bce_1=zeta_n(inode1);  %zeta of inode1
                                           znode_bce_2=zeta_n(inode2);  %zeta of inode2
                                           hnode_bce_1=h_n(inode1);     %h of inode1
                                           hnode_bce_2=h_n(inode2);     %h of inode2

                                           %length of edge
                                           dist_tce    = sqrt((vx(inode1)-vx(inode2))^2 + ...
                                                              (vy(inode1)-vy(inode2))^2);
                                           %linear interpolation to the starting point of the cutout polgyon 
                                           %(xpc(iedge),ypc(iedge))
                                           dist_node1= sqrt((xpc(iedge)-vx(inode1))^2 + ...
                                                            (ypc(iedge)-vy(inode1))^2) ;
                                           dist_node2= sqrt((xpc(iedge)-vx(inode2))^2 + ...
                                                            (ypc(iedge)-vy(inode2))^2) ;
                                           transect_seg_list{iseg_tot}.zeta(2)=(znode_bce_1*dist_node2+ ...
                                                                                znode_bce_2*dist_node1)/dist_tce; 
                                           transect_seg_list{iseg_tot}.h(2)=(hnode_bce_1*dist_node2+ ...
                                                                             hnode_bce_2*dist_node1)/dist_tce;

                                        elseif(iedge==length(polycuts{itce}{iseg}.itce_c(:))) 
                                           %last edge of the cutout polgyon
                                           inode1=ienode(ibce_edge,1);   %one end of the bce edge
                                           inode2=ienode(ibce_edge,2);   %other end of the bce edge
                                           znode_bce_1= zeta_n(inode1);  %zeta of inode1
                                           znode_bce_2= zeta_n(inode2);  %zeta depth of inode2
                                           hnode_bce_1=h_n(inode1);  %depth of inode1
                                           hnode_bce_2=h_n(inode2);  %depth of inode2

                                           dist_tce  = sqrt((vx(inode1)-vx(inode2))^2 + ...
                                                            (vy(inode1)-vy(inode2))^2);
 			                   dist_node1= sqrt((xpc(iedge)-vx(inode1))^2 + ...
                                                           (ypc(iedge)-vy(inode1))^2) ;
                                           dist_node2= sqrt((xpc(iedge)-vx(inode2))^2 + ...
                                                           (ypc(iedge)-vy(inode2))^2) ;
                                           transect_seg_list{iseg_tot}.zeta(1)=(znode_bce_1*dist_node2+ ...
                                                                                znode_bce_2*dist_node1)/dist_tce;
                                           transect_seg_list{iseg_tot}.h(1)   =(hnode_bce_1*dist_node2+ ...
                                                                                hnode_bce_2*dist_node1)/dist_tce; 
                                        else %normal full tce edges which also a boundary edge
                                             %do nothing
                                        end %iedge if
                                    end %ibce_edge if
                               end %itce_edge if
                           end %iedge loop
                       end  %valid cut
                   end %iseg loop
	       end %itce loop

%              tmp_debug=input('ctrl-c to stop 2');

               %sort the segments according to locations
               %----##################################################################################
               %           NOTE ON SORTING: SORTING using distance works if the transect
               %                                    is a straight line, otherwise, we have to
               %                                    use a parametric equation for a curved transect 
               %                                    to sort points on it.
               %----##################################################################################  
               NSegs_use=iseg_tot;
               sseg=zeros([NSegs_use,1]);
               xseg=zeros([NSegs_use,1]);
               yseg=zeros([NSegs_use,1]);

               %distance of each segment mid point to starting point of the transect
               for iseg=1:NSegs_use
                   xseg(iseg)= mean(transect_seg_list{iseg}.x(:));  %x coord of the mid point of the segment 
                   yseg(iseg)= mean(transect_seg_list{iseg}.y(:));  %y coord of the mid point of the segment
                   sseg(iseg)= sqrt((xseg(iseg) - xline(1))^2+ ...   
                                    (yseg(iseg) - yline(1))^2);     %distance to start point
               end

               [sseg_sort, II_sort]= sort( sseg,'ascend');          %sort for points on a straight line

               %plot the segments and label the segment number and also the depth 
               if(debug_plot)   %debug
                    figure; hold on;
                    for iseg=1:NSegs_use
                	   line(transect_seg_list{iseg}.x(:),transect_seg_list{iseg}.y(:),'color','r');
	                   plot(mean(transect_seg_list{iseg}.x(:)),mean(transect_seg_list{iseg}.y(:)),'ko');
        	           text(mean(transect_seg_list{iseg}.x(:)),mean(transect_seg_list{iseg}.y(:)), ...
	                        ['i=',num2str(iseg), ',h=',num2str(mean(transect_seg_list{iseg}.h(:)))]);
                    end
               end

               %calculate flux per segment length
               for iseg=1:NSegs_use
                   seg_length=sqrt((transect_seg_list{iseg}.x(2)-transect_seg_list{iseg}.x(1))^2 + ...
                                   (transect_seg_list{iseg}.y(2)-transect_seg_list{iseg}.y(1))^2); 

                   %depth integrated flux per unit segment length
                   cflow_use_mean(iseg)    =    cflow_use(iseg)/seg_length;  %segment mean flow (m^3/sec/m ==> m^2/sec) 
                   cflow_fw_use_mean(iseg) = cflow_fw_use(iseg)/seg_length;  
                   tflux_use_mean(iseg)    =    tflux_use(iseg)/seg_length;  %[c]m^3/sec/m == [c]m^2/sec

                   %vertical profile  (layer integrated flux) for  each layer per unit segment length
                   cflow_prf_use_mean(iseg,:)    = cflow_prf_use(iseg,:)/seg_length;    %m^3/sec/m
                   cflow_fw_prf_use_mean(iseg,:) = cflow_fw_prf_use(iseg,:)/seg_length; %m^3/sec/m
                   tflux_prf_use_mean(iseg,:)    = tflux_prf_use(iseg,:)/seg_length;    %[c]m^3/sec/m
               end

               %show distribution of the flux along the transect
               if(debug_plot)  %debug
	               figure; hold on;
	               plot(sseg, cflow_use_mean,'r*');
	               plot(sseg, cflow_fw_use_mean,'bo');
	               plot(sseg, tflux_use_mean,'m^'); 
               end

               %vertical distribution of the flux as function of layers
               %---divide by vertical thickness of the layer
               %   and then sum up all the segments of the transect

               zseg=zeros([NSegs_use,KBM1]);
               dz_seg_avg=zeros([1,KBM1]);

               for iseg=1:NSegs_use
                   seg_length=sqrt((transect_seg_list{iseg}.x(2)-transect_seg_list{iseg}.x(1))^2 + ...
                                   (transect_seg_list{iseg}.y(2)-transect_seg_list{iseg}.y(1))^2);
                   %get the vertical thickness of each layer (average thickness of the two ends of the segment)
                   dz_seg=zeros([1,KBM1]);
                   for k=1:KBM1
                       %thickness of the segment of layer k
                       dz_seg(k)=(z(k)-z(k+1))*( ...    %z is sigma coordinate
                                (transect_seg_list{iseg}.zeta(1) + transect_seg_list{iseg}.h(1)) + ...
                                (transect_seg_list{iseg}.zeta(2) + transect_seg_list{iseg}.h(2)))*0.5; 

		       dz_seg_avg(k)=dz_seg_avg(k)+dz_seg(k); 
                       %vertical coordinate of the segment of layer k
                       zseg(iseg,k)= (transect_seg_list{iseg}.zeta(1) + transect_seg_list{iseg}.zeta(2))*0.5 ... %surface
                                     +((transect_seg_list{iseg}.zeta(1) + transect_seg_list{iseg}.h(1)) + ...  
                                       (transect_seg_list{iseg}.zeta(2) + transect_seg_list{iseg}.h(2))   ...
                                      )*0.5 ...   %total depth
                                     *(z(k)+z(k+1))*0.5;    %sigma coord at center of layer


                   end

                   %vertical profile of flux density (flow or tracer flux per vertical layer thickness)
                   cflow_prfd_use(iseg,:)    =    cflow_prf_use(iseg,:)./dz_seg;     %m^3/sec/m       %flow
                   cflow_fw_prfd_use(iseg,:) = cflow_fw_prf_use(iseg,:)./dz_seg;     %m^3/sec/m       %fresh water
                   tflux_prfd_use(iseg,:)    =    tflux_prf_use(iseg,:)./dz_seg;     %[c]m^3/sec/m    %tracer flux

                   %piston velocity (flux per area)
                   cflow_prf_piston(iseg,:)    =    cflow_prfd_use(iseg,:)/seg_length; %m/sec (basically normal velocity)
                   cflow_fw_prf_piston(iseg,:) = cflow_fw_prfd_use(iseg,:)/seg_length; %m/sec (basically normal velocity)
                   tflux_prf_piston(iseg,:)    =    tflux_prfd_use(iseg,:)/seg_length; %[c]*m/sec   %transport

               end

               dz_seg_avg(:)=dz_seg_avg(:)./NSegs_use;                                 %averaged thickness of layrs

               %
               %collect the transect integrated flux 
               %i.e. integrate the flux along the transect to get total flow (m^3/sec), total fresh water flow (m^3/sec)
               %and also total tracer flux [c]*m^3/sec
               %

               if(it_all==1)
                  cflow_all    = sum(cflow_use(1:NSegs_use));
                  cflow_all_fw = sum(cflow_fw_use(1:NSegs_use));
                  tflux_all    = sum(tflux_use(1:NSegs_use));
               else
                  cflow_all    = [cflow_all    ; sum(   cflow_use(1:NSegs_use))];      %append one number for new row
                  cflow_all_fw = [cflow_all_fw ; sum(cflow_fw_use(1:NSegs_use))];
                  tflux_all    = [tflux_all    ; sum(   tflux_use(1:NSegs_use))];
               end

               cflow_prfd_tmp    = zeros([1,KBM1]);
               cflow_fw_prfd_tmp = zeros([1,KBM1]);
               tflux_prfd_tmp    = zeros([1,KBM1]);
               for k=1:KBM1
                   cflow_prfd_tmp(k)   =   sum(   cflow_prf_use(1:NSegs_use,k))/dz_seg_avg(k);
                   cflow_fw_prfd_tmp(k)=   sum(cflow_fw_prf_use(1:NSegs_use,k))/dz_seg_avg(k);
                   tflux_prfd_tmp(k)   =   sum(   tflux_prf_use(1:NSegs_use,k))/dz_seg_avg(k);
               end

               %
               %collect the flux density (flux per averaged layer thickness) integrated over segments
               %

               if(it_all==1)
                   cflow_prfd    =[cflow_prfd_tmp];
                   cflow_fw_prfd =[ cflow_fw_prfd_tmp ];
                   tflux_prfd    =[tflux_prfd_tmp];
               else
                   cflow_prfd = [ cflow_prfd; cflow_prfd_tmp];            %append a new row
                   cflow_fw_prfd = [cflow_fw_prfd ; cflow_fw_prfd_tmp];
                   tflux_prfd = [ tflux_prfd; tflux_prfd_tmp ];
               end

               %output the flux to tecplot
               tdata2d_fvf=[] ; %finite volume flux 2d tecplot
               tdata3d_fvf=[] ; %                   3d tecplot

	       switch(i_par)
	            case  1 
	                tdata2d_fvf.Nvar=9;
	                tdata2d_fvf.varnames={'r','z','zdummy',...
                                       'flow_up','flow_up_fw','transport',  ...    %piston velocity and transport
                                       'flow_prfd','flow_prfd_fw','tflux_prfd'}; %flow or flux desnity (per unit depth)
	                tdata3d_fvf.Nvar=10;
	                tdata3d_fvf.varnames={'x','y','z',                     ...  %piston velocity and transport
                                      'flow_up','flow_up_fw','transport',      ...  %flow or flux density (per unit depth)
                                      'flow_prfd','flow_prfd_fw','tflux_prfd', ...  %flow or flux density
                                      'depth'};
	       end

               if(strcmp(parname,'velocity'))  %save velocity transects
                   tdata2d_fvf.surfaces(1).zonename=[parname ' transect ' transname, ...
                                      ' 2D zone file number ' num2str(nf,'%06d') ...
                                      ' time reocrd ' num2str(it,'%8d')];
                   tdata2d_fvf.surfaces(1).x=repmat(sseg_sort,[1,KBM1]);   %'x' is in meter as distance on the segment
                   tdata2d_fvf.surfaces(1).y=zseg(II_sort,:);             %'y' is actually vertical coordinate (m)
                   tdata2d_fvf.surfaces(1).z=zeros([NSegs_use,KBM1]);      %dummy z value

                   %piston velocity and transport
                   tdata2d_fvf.surfaces(1).v(1,:,:) = cflow_prf_piston(II_sort,:);    %piston velocity u_p (m/sec)
                   tdata2d_fvf.surfaces(1).v(2,:,:) = cflow_fw_prf_piston(II_sort,:); %fresh water piston velocity (m/sec)  
                   tdata2d_fvf.surfaces(1).v(3,:,:) = tflux_prf_piston(II_sort,:);    %tracer flux transport [c]*m/sec
                   %flux density
                   tdata2d_fvf.surfaces(1).v(4,:,:) = cflow_prfd_use(II_sort,:);      %m^3/sec/m
                   tdata2d_fvf.surfaces(1).v(5,:,:) = cflow_fw_prfd_use(II_sort,:);   %m^3/sec/m
                   tdata2d_fvf.surfaces(1).v(6,:,:) = tflux_prfd_use(II_sort,:);      %[c]*m^3/sec/m

                   tdata2d_fvf.surfaces(1).varloc=0;       %variables are provided nodal
                   tdata2d_fvf.surfaces(1).order=3;        %surface is definded on (x,y) in tecplot, 
                                                           %i.e. ('sseg','zseg') here

                   tdata2d_fvf.surfaces(1).solutiontime=htime(it);

                   %need to worry about nan values in x, y, z, v here when
                   %the transects run out of domain (replace nans with
                   %-999999999.0 so that tecplot can blank them out)

                   inan_x=find(isfinite(tdata2d_fvf.surfaces(1).x)==0);
                   tdata2d_fvf.surfaces(1).x(inan_x)=-999999999.0;       %give big negative value so that we can blank it
                   inan_y=find(isfinite(tdata2d_fvf.surfaces(1).y)==0);
                   tdata2d_fvf.surfaces(1).y(inan_y)=-999999999.0;       %give big negative value so that we can blank it
                   inan_v=find(isfinite(tdata2d_fvf.surfaces(1).v)==0);  %
                   tdata2d_fvf.surfaces(1).v(inan_v)=-999999999.0;
                   tdata2d_fvf.surfaces(1).v(1,inan_x)=-999999999.0;
                   tdata2d_fvf.surfaces(1).v(1,inan_y)=-999999999.0;
                   tdata2d_fvf.surfaces(1).v(2,inan_x)=-999999999.0;
                   tdata2d_fvf.surfaces(1).v(2,inan_y)=-999999999.0;
                   tdata2d_fvf.surfaces(1).v(3,inan_x)=-999999999.0;
                   tdata2d_fvf.surfaces(1).v(3,inan_y)=-999999999.0;
                   tdata2d_fvf.surfaces(1).v(4,inan_x)=-999999999.0;
                   tdata2d_fvf.surfaces(1).v(4,inan_y)=-999999999.0;
                   tdata2d_fvf.surfaces(1).v(5,inan_x)=-999999999.0;
                   tdata2d_fvf.surfaces(1).v(5,inan_y)=-999999999.0;
                   tdata2d_fvf.surfaces(1).v(6,inan_x)=-999999999.0;
                   tdata2d_fvf.surfaces(1).v(6,inan_y)=-999999999.0;

                   %3D tecplot in x, y, z coordinates of the transect (x, y, z are 2D ordered arrays)
                   tdata3d_fvf.surfaces(1).zonename=[           parname ' transect ' transname, ...
                                                 ' 3D zone file number ' num2str(nf,'%06d') ...
                                                 ' time record ' num2str(it,'%8d')];
                   tdata3d_fvf.surfaces(1).x=repmat(xseg,[1,KBM1]);    
                   %first dimension is number of points along xp
                   tdata3d_fvf.surfaces(1).y=repmat(yseg,[1,KBM1]);
                   tdata3d_fvf.surfaces(1).z=zseg; 

                   %piston velocity 
                   tdata3d_fvf.surfaces(1).v(1,:,:) = cflow_prf_piston(II_sort,:);    %piston velocity u_p (m/sec)
                   tdata3d_fvf.surfaces(1).v(2,:,:) = cflow_fw_prf_piston(II_sort,:); %fresh water piston velocity (m/sec)
                   tdata3d_fvf.surfaces(1).v(3,:,:) = tflux_prf_piston(II_sort,:);    %tracer flux transport [c]*m/sec
                   %flux density
                   tdata3d_fvf.surfaces(1).v(4,:,:) = cflow_prfd_use(II_sort,:);      %m^3/sec/m
                   tdata3d_fvf.surfaces(1).v(5,:,:) = cflow_fw_prfd_use(II_sort,:);   %m^3/sec/m
                   tdata3d_fvf.surfaces(1).v(6,:,:) = tflux_prfd_use(II_sort,:);   
                   tdata3d_fvf.surfaces(1).v(7,:,:) = zseg;  %last variable in v is -depth
                   tdata3d_fvf.surfaces(1).order=2;          %surface is IK-odered (I is first dimension of x, y, z, v, K
                                                             %is 2nd dimension of of x, y, z, v)
                                                             %surfaces are defined on (x,z) plane
                   tdata3d_fvf.surfaces(1).solutiontime=htime(it);
                   %need to worry abut nan values in x, y, z, v
                   inan_x=find(isfinite(tdata3d_fvf.surfaces(1).x)==0);
                   tdata3d_fvf.surfaces(1).x(inan_x)=-999999999.0;  %give big negative value so that we can blank it
                   inan_y=find(isfinite(tdata3d_fvf.surfaces(1).y)==0);
                   tdata3d_fvf.surfaces(1).y(inan_y)=-999999999.0;  %give big negative value so that we can blank it
                   inan_z=find(isfinite(tdata3d_fvf.surfaces(1).z)==0);
                   tdata3d_fvf.surfaces(1).z(inan_z)=-999999999.0;
                   inan_v=find(isfinite(tdata3d_fvf.surfaces(1).v)==0);
                   tdata3d_fvf.surfaces(1).v(inan_v)=-999999999.0;
                   tdata3d_fvf.surfaces(1).v(1,inan_x)=-999999999.0;
                   tdata3d_fvf.surfaces(1).v(1,inan_y)=-999999999.0;
                   tdata3d_fvf.surfaces(1).v(1,inan_z)=-999999999.0;
                   tdata3d_fvf.surfaces(1).v(2,inan_x)=-999999999.0;
                   tdata3d_fvf.surfaces(1).v(2,inan_y)=-999999999.0;
                   tdata3d_fvf.surfaces(1).v(2,inan_z)=-999999999.0;
                   tdata3d_fvf.surfaces(1).v(3,inan_x)=-999999999.0;
                   tdata3d_fvf.surfaces(1).v(3,inan_y)=-999999999.0;
                   tdata3d_fvf.surfaces(1).v(3,inan_z)=-999999999.0;
                   tdata3d_fvf.surfaces(1).v(4,inan_x)=-999999999.0;
                   tdata3d_fvf.surfaces(1).v(4,inan_y)=-999999999.0;
                   tdata3d_fvf.surfaces(1).v(4,inan_z)=-999999999.0;
                   tdata3d_fvf.surfaces(1).v(5,inan_x)=-999999999.0;
                   tdata3d_fvf.surfaces(1).v(5,inan_y)=-999999999.0;
                   tdata3d_fvf.surfaces(1).v(5,inan_z)=-999999999.0;
                   tdata3d_fvf.surfaces(1).v(6,inan_x)=-999999999.0;
                   tdata3d_fvf.surfaces(1).v(6,inan_y)=-999999999.0;
                   tdata3d_fvf.surfaces(1).v(6,inan_z)=-999999999.0;
                   tdata3d_fvf.surfaces(1).v(7,inan_y)=-999999999.0;
                   tdata3d_fvf.surfaces(1).v(7,inan_z)=-999999999.0;

                   %for 3D transect tecplot, also add bathymtry as a zone
                   %ideally we do not have to have this for every time step
                   %but that then makes it hard to do movie with time steps

                   tdata3d_fvf.FEsurfaces(1).zonename=['PSWQM 2.0 grid file number '  ...
                                   num2str(nf,'%06d') ' time record ' num2str(it,'%8d')];
                   tdata3d_fvf.FEsurfaces(1).x=x_n;   %x of all nodes
                   tdata3d_fvf.FEsurfaces(1).y=y_n;   %y of all nodes
                   tdata3d_fvf.FEsurfaces(1).z=-h_n;  %z of all nodes is given as -h
                   tdata3d_fvf.FEsurfaces(1).order=3; %order 3, surface is defined on (x, y) coordinates
                   tdata3d_fvf.FEsurfaces(1).e2n=e2n(:,2:4);
                   tdata3d_fvf.FEsurfaces(1).v(1,:)=zeros(size(h_n')); %dummy values as place holder for u
                   tdata3d_fvf.FEsurfaces(1).v(2,:)=zeros(size(h_n')); %dummy values as place holder for v
                   tdata3d_fvf.FEsurfaces(1).v(3,:)=zeros(size(h_n')); %dummy values as place holder for w
                   tdata3d_fvf.FEsurfaces(1).v(4,:)=zeros(size(h_n')); %dummy values as place holder for u_p
                   tdata3d_fvf.FEsurfaces(1).v(5,:)=zeros(size(h_n')); %dummy values as place holder for v_p
                   tdata3d_fvf.FEsurfaces(1).v(6,:)=zeros(size(h_n')); %dummy values as place holder for v_p
                   tdata3d_fvf.FEsurfaces(1).v(7,:)=-h_n';             %-depth gives the depth on this surface
                   tdata3d_fvf.FEsurfaces(1).solutiontime=htime(it);

                  %save the tecplot data into matlab file first, will load them after the nf loop
                   tmpfilename=[pltdirname './' parname,'_trans_', transname, '_fvf_tecplot_', ...
                                num2str(nf,'%06d') '_' num2str(it,'%010d'),'.mat'];

                   save(tmpfilename,'tdata2d_fvf', ...
                                    'tdata3d_fvf', ...
                                    '-mat');

%                   tmp_debug=input('ctrl-c to stop 3');

               end  %end 'velocity' if 
     end     %end of time loop (it_global)

%
%Aggregate all tecplot plots of individual time records into a big tecplot history file
%so that it can be played as an animation in tecplot
%

     tdata2d_all=[];  %2D and 3D zones that have all the time history
     tdata3d_all=[];

     tdata2d_timeavg=[];  %time average of the 2D and 3D plots
     tdata3d_timeavg=[];

     it_all=0;

     for it_global=1:t_int:NT_global
         it_all=it_all+1;
         %find the tecplot data stored in matlab to load
         nf=filenum_global(it_global);
         it_local=it_local_global(it_global);
         it=it_local;
         tmpfilename=[pltdirname './'  parname,'_trans_', transname, '_fvf_tecplot_', ...
                            num2str(nf,'%06d') '_' num2str(it,'%010d'),'.mat'];
         load(tmpfilename); %load tdata2d_fvf and tdata3d_fvf
         if(it_all==1)
             tdata2d_all.Nvar=tdata2d_fvf.Nvar;
             tdata2d_all.varnames=tdata2d_fvf.varnames;
         end
         tdata2d_all.surfaces(it_all)=tdata2d_fvf.surfaces;
         if(it_all==1)
           tdata3d_all.Nvar=tdata3d_fvf.Nvar;
           tdata3d_all.varnames=tdata3d_fvf.varnames;
         end
         tdata3d_all.surfaces(it_all)=tdata3d_fvf.surfaces;
         tdata3d_all.FEsurfaces(it_all)=tdata3d_fvf.FEsurfaces;

         %fill out variable information for tdata2d_timeavg and tdata3d_timeavg
         if(it_all==1)
	         tdata2d_timeavg.Nvar=tdata2d_fvf.Nvar;
	         tdata2d_timeavg.varnames=tdata2d_fvf.varnames;
	         tdata3d_timeavg.Nvar=tdata3d_fvf.Nvar;
	         tdata3d_timeavg.varnames=tdata3d_fvf.varnames;
         end
     end

     %output the transect time history into tecplot format
     tecfile2d=[pltdirname './' parname,'_trans_', transname,'_', num2str(his_fn_start,'%06d'), ...
                                   '_', num2str(his_fn_end,'%06d'), '_fvflux_2D_contour.plt'];
     mat2tecplot(tdata2d_all,tecfile2d);
     tecfile3d=[pltdirname './' parname,'_trans_', transname,'_', num2str(his_fn_start,'%06d'), ...
                                   '_', num2str(his_fn_end,'%06d'), '_fvflux_3D_contour.plt'];
     mat2tecplot(tdata3d_all,tecfile3d);

%     tmp_debug=input('ctrl-c to stop 4');


  %
  %Create a time average for the 2D and 3D zones
  %

     tdata2d_timeavg.surfaces(1).zonename=[parname ' transect ' transname, ' 2D zone time average from ' ...
                                           num2str(htime_all(1),'%8.2f') ' to ' num2str(htime_all(end),'%8.2f')];
     tdata2d_timeavg.surfaces(1).x=zeros(size(tdata2d_fvf.surfaces(1).x));
     tdata2d_timeavg.surfaces(1).y=zeros(size(tdata2d_fvf.surfaces(1).y));
     tdata2d_timeavg.surfaces(1).z=zeros(size(tdata2d_fvf.surfaces(1).z));
     tdata2d_timeavg.surfaces(1).v=zeros(size(tdata2d_fvf.surfaces(1).v));
     tdata2d_timeavg.surfaces(1).varloc=0;   %variables are provided nodal
     tdata2d_timeavg.surfaces(1).order=3;    %surface is definded on (x,y) in tecplot, i.e. ('r','z') here

     tdata3d_timeavg.surfaces(1).zonename=[parname ' transect ' transname, ' 3D zone time average from ' ...
                                           num2str(htime_all(1),'%8.2f') ' to ' num2str(htime_all(end),'%8.2f')];
     tdata3d_timeavg.surfaces(1).x=zeros(size(tdata3d_fvf.surfaces(1).x));
     tdata3d_timeavg.surfaces(1).y=zeros(size(tdata3d_fvf.surfaces(1).y));
     tdata3d_timeavg.surfaces(1).z=zeros(size(tdata3d_fvf.surfaces(1).z));
     tdata3d_timeavg.surfaces(1).v=zeros(size(tdata3d_fvf.surfaces(1).v));
     tdata3d_timeavg.surfaces(1).order=2;    %surface is IK-odered (I is first dimension of x, y, z,
                                             %v, K is 2nd dimension of of x, y, z, v)
                                             %surfaces are defined on (x,z) plane
     %add bathymetry zone
     tdata3d_timeavg.FEsurfaces(1).zonename=['FVCOM grid and bathymetry'];
     tdata3d_timeavg.FEsurfaces(1).x=x_n;   %x of all nodes
     tdata3d_timeavg.FEsurfaces(1).y=y_n;   %y of all nodes
     tdata3d_timeavg.FEsurfaces(1).z=-h_n;  %z of all nodes is given as -h
     tdata3d_timeavg.FEsurfaces(1).order=3; %order 3, surface is defined on (x, y) coordinates
     tdata3d_timeavg.FEsurfaces(1).e2n=e2n(:,2:4);
     for ivar=1:tdata3d_timeavg.Nvar-3-1    %first 3 is x, y,z, then rest of them has to be filled up with dummy valuses
        tdata3d_timeavg.FEsurfaces(1).v(ivar+3,:)=zeros(size(h_n'));%dummy values as place holder for parname (T,S, kh)
     end
     tdata3d_timeavg.FEsurfaces(1).v(tdata3d_timeavg.Nvar,:)=-h_n'; %-depth gives the depth on this surfacs (last variable)

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
     tecfile2d_timeavg=[pltdirname './' parname,'_trans_', transname,'_',...
             num2str(his_fn_start,'%06d'), '_', num2str(his_fn_end,'%06d'), '_fvflux_2D_contour_timeavg.plt'];
     mat2tecplot(tdata2d_timeavg,tecfile2d_timeavg);
     tecfile3d_timeavg=[pltdirname './' parname,'_trans_', transname,'_',...
             num2str(his_fn_start,'%06d'), '_', num2str(his_fn_end,'%06d'), '_fvflux_3D_contour_timeavg.plt'];
     mat2tecplot(tdata3d_timeavg,tecfile3d_timeavg);


     %Finally save everything about the tecplot data in a matlab file for future inspection
     save([pltdirname './' parname '_trans_' transname '_' num2str(his_fn_start,'%06d') ...
                              '_' num2str(his_fn_end,'%06d') '_fvflux_calculation.mat'],                ...
           'tdata2d_all','tdata3d_all','tdata2d_timeavg','tdata3d_timeavg');

%     tmp_debug=input('ctrl-c to stop 5');

     clear 'tdata2d_fvf' 'tdata3d_fvf' 'tdata2d_timeavg' 'tdata3d_timeavg' ;
     clear 'tdata2d_all' 'tdata3d_all';

  %
  %time series of total flux (transect cross-sectional integration)
  %

         tdata_TS=[] ;   
         tdata_TS.Nvar=2;
         tdata_TS.varnames={'time','flux'};
         tdata_TS.lines(1).zonename=[transname ' total flow time series'];
         tdata_TS.lines(1).x=htime_all;  %time
         tdata_TS.lines(1).y=cflow_all;      %flux (m^3/sec)
         tdata_TS.lines(2).zonename=[transname ' fresh water flow time series'];
         tdata_TS.lines(2).x=htime_all;
         tdata_TS.lines(2).y=cflow_all_fw;   %flux (m^3/sec) of fresh water
         tdata_TS.lines(3).zonename=[transname ' tracer flux time series'];
         tdata_TS.lines(3).x=htime_all;
         tdata_TS.lines(3).y=tflux_all;      %tracer flux [c]*(m^3/sec)

         tecfile_TS=[pltdirname './' parname,'_trans_', transname, '_', ...
                num2str(his_fn_start,'%06d'), '_', num2str(his_fn_end,'%06d'),'_fvflux_timeseries.plt'];
         mat2tecplot(tdata_TS,tecfile_TS);
  
%         tmp_debug=input('ctrl-c to stop 6');

         clear 'tdata_TS';

  %
  %time series of flux profile
  %
         tdata_Prf=[];  
         tdata_Prf.Nvar=2;
         tdata_Prf.varnames={'z','flux'};
         for it=1:length(htime_all)
             tdata_Prf.lines(it).zonename=['volume flux zone of transect ' ...
                                           transname];
             tdata_Prf.lines(it).x=(1:KBM1);
             tdata_Prf.lines(it).y=cflow_prfd(it,:)';         %transpose from it's row to column
             tdata_Prf.lines(it).solutiontime=htime_all(it);  %flux at this time
         end
         tecfile_Prf=[pltdirname './' parname,'_trans_', transname, '_', ...
                               num2str(his_fn_start,'%06d'), ...
                          '_', num2str(his_fn_end,'%06d') '_fvflux_profile_timeseries_volume.plt'];
         mat2tecplot(tdata_Prf,tecfile_Prf);

         

%         tmp_debug=input('ctrl-c to stop 7');

         %fresh water volume flux
         tdata_Prf=[];   
         tdata_Prf.Nvar=2;
         tdata_Prf.varnames={'z','flux'};
         for it=1:length(htime_all)
             tdata_Prf.lines(it).zonename=['fresh water volume flux zone of transect ' ...
                         transname];
             tdata_Prf.lines(it).x=(1:KBM1);
             tdata_Prf.lines(it).y=cflow_fw_prfd(it,:)';      %transpose from it's row to column
             tdata_Prf.lines(it).solutiontime=htime_all(it);  %flux at this time
         end

         tecfile_Prf=[pltdirname './' parname,'_trans_', transname, '_',...
                               num2str(his_fn_start,'%06d'), ...
                          '_', num2str(his_fn_end,'%06d') '_fvflux_profile_timeseries_freshwater.plt'];
         mat2tecplot(tdata_Prf,tecfile_Prf);

%         tmp_debug=input('ctrl-c to stop 8');


         %tracer flux
         tdata_Prf=[]; 
         tdata_Prf.Nvar=2;
         tdata_Prf.varnames={'z','flux'};
         for it=1:length(htime_all)
             tdata_Prf.lines(it).zonename=['tracer flux zone of transect ' ...
                         transname];
             tdata_Prf.lines(it).x=(1:KBM1);
             tdata_Prf.lines(it).y=tflux_prfd(it,:)';         %transpose from it's row to column
             tdata_Prf.lines(it).solutiontime=htime_all(it);  %
         end
         tecfile_Prf=[pltdirname './' parname,'_trans_', transname,'_', ...
                               num2str(his_fn_start,'%06d'), ...
                          '_', num2str(his_fn_end,'%06d'),   ...
                          '_fvflux_profile_timeseries_'      ...
                          tracername '.plt'];

         mat2tecplot(tdata_Prf,tecfile_Prf);

%        tmp_debug=input('ctrl-c to stop 9');

%         clear 'tdata_Prf';
  %
  %time average to get mean profile
  %
         tdata_Prf_mean=[];
         tdata_Prf_mean.Nvar=2;
         tdata_Prf_mean.varnames={'z','flux'};

         %volume flux
         tdata_Prf_mean.lines(1).zonename=[tracername ' mean volume flux zone transect ' ...
                                            transname ' during time from '        ...
                        num2str(htime_all(1),'%8.4f') ' to ' num2str(htime_all(end),'%8.4f')];
         tdata_Prf_mean.lines(1).x=(1:KBM1)';                             %vertical coordinate as layers
         tdata_Prf_mean.lines(1).y=(mean(cflow_prfd(:,1:end),1))';        %average the profile over time

         %fresh water volume flux
         tdata_Prf_mean.lines(2).zonename=[tracername ' mean fresh water volume flux zone transect ' ...
                                            transname ' during time from '        ...
                        num2str(htime_all(1),'%8.4f') ' to ' num2str(htime_all(end),'%8.4f')];
         tdata_Prf_mean.lines(2).x=(1:KBM1)';                             %vertical coordinate as layers
         tdata_Prf_mean.lines(2).y=(mean(cflow_fw_prfd(:,1:end),1))';     %average the profile over time

         %tracer flux 
         tdata_Prf_mean.lines(3).zonename=[tracername ' mean tracer flux zone transect ' ...
                                            transname ' during time from '        ...
                        num2str(htime_all(1),'%8.4f') ' to ' num2str(htime_all(end),'%8.4f')];
         tdata_Prf_mean.lines(3).x=(1:KBM1)';                             %vertical coordinate as layers
         tdata_Prf_mean.lines(3).y=(mean(tflux_prfd(:,1:end),1))';        %average the profile over time

         tecfile_Prf_mean=[pltdirname './' parname,'_trans_', transname,'_',   ...
                                    num2str(his_fn_start,'%06d'), ...
                                   '_' num2str(his_fn_end,'%06d'), '_fvflux_profile_timeavg.plt'];

         %produce the tecplot file
         mat2tecplot(tdata_Prf_mean,tecfile_Prf_mean);

%         tmp_debug=input('ctrl-c to stop 10');

%         clear 'tdata_Prf_mean';

         save([pltdirname './' 'prfile_timeseries.mat'], ...
              'cflow_prfd','cflow_fw_prfd','tflux_prfd','htime_all','tdata_Prf', ...
              'KBM1','his_fn_start','his_fn_end','tracername','transname','pltdirname', ...
              'tdata_Prf_mean' );

 display(['Success, plots are made in folder ' pltdirname]); 


