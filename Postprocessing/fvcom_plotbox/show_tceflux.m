%simple matlab code to show the TCE flux 

%read the tce flux 
tflux=load('tceflux_output.dat');

%TCE flux is stored in the tflux with the following fortran statements

%
%	   WRITE(107,'(F12.5,1X,I8,1X,I8,1X,I8,1X,I8,1X,4000(F12.5,1X))') THOUR,                                 &
%				   		                  I_edgeline_ele_GL(I),                          &
%							          I_edgeline_ele_edge_ia_GL(I),                  &
%							          I_edgeline_ele_edge_ib_GL(I),                  &
%							          I_edgeline_flux_sign_GL(I),                    &
%							          (EDGELINE_FLUX_SALT_GL(I,K),     K=1,KBM1),    &
%							          (EDGELINE_FLUX_SALTWATER_GL(I,K),K=1,KBM1),    &
%							          (EDGELINE_FLUX_FRSHWATER_GL(I,K),K=1,KBM1),    &
%							          (EDGELINE_FLUX_WATER_GL(I,K),    K=1,KBM1 )
%

  %THOUR ele ia ib sign   salt(1:KBM1) saltwater(1:KBM1) freshwater(1:KBM1) water(1:KBM1)

   KBM1= (size(tflux,2)-5) /4; %number of layers

   t_hour=tflux(:,1);  %time in hours

   t_day= t_hour/24.0;  %time in days

   t_day_uniq=unique(t_day);
   t_hour_uniq=unique(t_hour);

   NT=length(t_hour_uniq); 

   NEL_GL = length(t_hour)/length(t_hour_uniq);  %number of edge segments on the chosen edgeline

   ele      =    tflux(1:NEL_GL,2); 
   ia       =    tflux(1:NEL_GL,3);
   ib       =    tflux(1:NEL_GL,4);
   flx_sign =    tflux(1:NEL_GL,5);

 %loop over time and plot the depth profile of salt flux 

%  hfig=figure('visible','on'); 
  saltflux_timeseries=[];
  frshflux_timeseries=[];
  waterflux_timeseries=[];
  swtrflux_timeseries=[];

  for it=1:NT
%      hold on;
      this_hr=t_hour_uniq(it);
      this_d =t_day_uniq(it);
      II=find(t_hour==this_hr);


      %get the salt flux of this time for all the segments
      ipar=1;

      jstart=6+(ipar-1)*KBM1;
      jend  =jstart+KBM1-1;

      saltflux=tflux(II,jstart:jend); %fetch the time and the columns for saltflux
      %multiply the sign to get flux as positive to the right hand side of the edge line
      for itce=1:NEL_GL
          saltflux(itce,:)=saltflux(itce,:)*flx_sign(itce);
      end
      saltflux_sum=sum(saltflux(:,:),1);  %sum up all edge segments to get vertical profile
%     plot(saltflux_sum,-(1:1:KBM1));  %plot vertical profile as from surface(-1) to bottom (-KBM1) layers
%     title(['salt flux (psu*m^3/sec) at time tday=',num2str(this_d)]);
      saltflux_timeseries=[saltflux_timeseries; saltflux_sum]; 

      %salt water flux
      jstart=6+(ipar-1)*KBM1;
      jend  =jstart+KBM1-1;

      swtrflux=tflux(II,jstart:jend); %fetch the time and the columns for swtrflux
      %multiply the sign to get flux as positive to the right hand side of the edge line
      for itce=1:NEL_GL
          swtrflux(itce,:)=swtrflux(itce,:)*flx_sign(itce);
      end
      swtrflux_sum=sum(swtrflux(:,:),1);  %sum up all edge segments to get vertical profile
%     plot(swtrflux_sum,-(1:1:KBM1));  %plot vertical profile as from surface(-1) to bottom (-KBM1) layers
%     title(['salt water flux (m^3/sec) at time tday=',num2str(this_d)]);
      swtrflux_timeseries=[swtrflux_timeseries; swtrflux_sum];

      %freshwater flux
      ipar=3;
      jstart=6+(ipar-1)*KBM1;
      jend  =jstart+KBM1-1;

      frshflux=tflux(II,jstart:jend); %fetch the time and the columns for frshflux
      %multiply the sign to get flux as positive to the right hand side of the edge line
      for itce=1:NEL_GL
          frshflux(itce,:)=frshflux(itce,:)*flx_sign(itce);
      end
      frshflux_sum=sum(frshflux(:,:),1);  %sum up all edge segments to get vertical profile
%     plot(frshflux_sum,-(1:1:KBM1));  %plot vertical profile as from surface(-1) to bottom (-KBM1) layers
%     title(['freshwater (m^3/sec) at time tday=',num2str(this_d)]);
      frshflux_timeseries=[frshflux_timeseries; frshflux_sum];

      %water volume flux
      ipar=4;  
      jstart=6+(ipar-1)*KBM1;
      jend  =jstart+KBM1-1;
      waterflux=tflux(II,jstart:jend); %fetch the time and the columns for waterflux
      %multiply the sign to get flux as positive to the right hand side of the edge line
      for itce=1:NEL_GL
          waterflux(itce,:)=waterflux(itce,:)*flx_sign(itce);
      end
      waterflux_sum=sum(waterflux(:,:),1);  %sum up all edge segments to get vertical profile
%     plot(waterflux_sum,-(1:1:KBM1));  %plot vertical profile as from surface(-1) to bottom (-KBM1) layers
%     title(['total water flux (m^3/sec) at time tday=',num2str(this_d)]);
      waterflux_timeseries=[waterflux_timeseries; waterflux_sum];

%     pause(1)
%     hold off;

  end

  %plot time average of the flux

   hold on;  
   plot(mean(saltflux_timeseries ,1),-(1:1:KBM1),'r-');

   figure('visible','on'); hold on;
   h1=plot(mean(frshflux_timeseries ,1),-(1:1:KBM1),'bo-');
   h2=plot(mean(swtrflux_timeseries ,1),-(1:1:KBM1),'ms-');
   h3=plot(mean(waterflux_timeseries,1),-(1:1:KBM1),'g^-');
   legend([h1, h2,h3],'freshwater','saltwater','totalwater');
   xlabel('flux (m^3/s)');
   ylabel('-layer'); 
   title('water flux (m^3/sec)');

