
%find the line (vertical transect) to plot

 iline=find(y_e< 1600 & y_e>1400);
 NX=length(iline);

 xiline=x_e(iline);
 yiline=y_e(iline);
 
 xiline2D_level=repmat(xiline,[1,KB]);
 xiline2D_layer=repmat(xiline,[1,KBM1]);
 yiline2D_level=repmat(yiline,[1,KB]);
 yiline2D_layer=repmat(yiline,[1,KBM1]);

 hiline=h_e(iline);

%zetailine_1=zeta_e_1(iline);

 ziline2D_level=zeros(NX,KB);
 ziline2D_layer=zeros(NX,KBM1);

%
%loop in time to make plots of vertical transect
%

 tdata2d=[];
 tdata3d=[];

 tdata2d_w.Nvar=9;
 tdata2d_w.varnames={'x',          'y',          'z',       ...
                     'w_x',        'w_y',        'w_z',          ...
                     'w_fvcom_x',  'w_fvcom_y',  'w_fvcom_z',    ...
                    };
 
 tdata3d_w.Nvar=10;
 tdata3d_w.varnames={'x','y','z',
                     'w_x',        'w_y',        'w_z',          ...
                     'w_fvcom_x',  'w_fvcom_y',  'w_fvcom_z',    ...
		     'zlevel'};

 tdata2d_uvww.Nvar=12;
 tdata2d_uvww.varnames={'x',          'y',          'z',    ...
                        'u_x',        'u_y',        'u_z',       ...
                        'v_x',        'v_y',        'v_z',       ...
                        'ww_fvcom_x', 'ww_fvcom_y', 'ww_fvcom_z'
                       };

 tdata3d_uvww.Nvar=13;
 tdata3d_uvww.varnames={'x','y','z',
                        'u_x',        'u_y',        'u_z',       ...
                        'v_x',        'v_y',        'v_z',       ...
                        'ww_fvcom_x', 'ww_fvcom_y', 'ww_fvcom_z',...
                        'zlayer'};

for it_all=1:NT_global

  %get surface elevation at node
   zeta_n_tmp=squeeze(zeta_n_all(it_all,:)); %zeta on node
   
  %get zeta at element center (averaging) of all 3 nodes of an element
   zeta_e_tmp=zeros([1,n]); 

   for ie=1:n
       zeta_e_tmp(ie)=mean(zeta_n_tmp(e2n(ie,2:4)));
   end

  %get zeta on the transect
   zetailine=zeta_e_tmp(iline);
   
  %get vertical coordinate

    for k=1:KB
       ziline2D_level(:,k)=zetailine+z(k)*(zetailine+hiline);
    end

    for k=1:KBM1
       ziline2D_layer(:,k)=zetailine+(z(k)+z(k+1))*0.5*(zetailine+hiline);
    end

  %get vertical velocities
     wiline=w_all(it_all,iline,:);              %pseudo-vertical velocity by matlab

     wiline_fvcom=w_fvcom_all(it_all,iline,:);  %psuedo-vertical velocity by fvcom

     wwiline=ww_all(it_all,iline,:);            %true vertical velocity
 
  %get horizontal velocities (u,v)

     uiline = u_all(it_all,iline,:);
     viline = v_all(it_all,iline,:); 

   %make plots
     
   %add to tecplot

    %2D plot showing w and w_fvcom (for comparison purpose)
    tdata2d_w.surfaces(it_all).zonename=['w-pseudo at time step ' num2str(it_all)];
    tdata2d_w.surfaces(it_all).x=xiline2D_level;
    tdata2d_w.surfaces(it_all).y=yiline2D_level;
    tdata2d_w.surfaces(it_all).z=ziline2D_level;
    tdata2d_w.surfaces(it_all).order=2;         %surface defined on x and z
    tdata2d_w.surfaces(it_all).v(1,:,:)=zeros(size(wiline))    ;%w_x ==0
    tdata2d_w.surfaces(it_all).v(2,:,:)=zeros(size(wiline))    ;%w_y ==0
    tdata2d_w.surfaces(it_all).v(3,:,:)=w                      ;%w_z ==w
    tdata2d_w.surfaces(it_all).v(4,:,:)=zeros(size(wiline))    ;%w_fvcom_x==0
    tdata2d_w.surfaces(it_all).v(5,:,:)=zeros(size(wiline))    ;%w_fvcom_y==0
    tdata2d_w.surfaces(it_all).v(6,:,:)=w_fvcom                ;%w_fvcom_z==w_fvcom
    tdata2d_w.surfaces(it_all).solutiontime=htime_all(it_all); 

    %2D plot showing u,v,ww
    tdata2d_uvww.surfaces(it_all).zonename=['u,v,ww at time step ' num2str(it_all)]; 
    tdata2d_uvww.surfaces(it_all).x=xiline2D_layer;
    tdata2d_uvww.surfaces(it_all).y=yiline2D_layer;
    tdata2d_uvww.surfaces(it_all).z=ziline2D_layer;
    tdata2d_uvww.surfaces(it_all).order=2;      %surface defined on x and z
    tdata2d_uvww.surfaces(it_all).v(1,:,:)=uiline                 ;%u_x ==u
    tdata2d_uvww.surfaces(it_all).v(2,:,:)=zeros(size(uiline))    ;%u_y ==0
    tdata2d_uvww.surfaces(it_all).v(3,:,:)=zeros(size(uiline))    ;%u_z ==0
    tdata2d_uvww.surfaces(it_all).v(4,:,:)=zeros(size(viline))    ;%v_x ==0
    tdata2d_uvww.surfaces(it_all).v(5,:,:)=viline                 ;%v_y ==v
    tdata2d_uvww.surfaces(it_all).v(6,:,:)=zeros(size(viline))    ;%v_z ==0
    tdata2d_uvww.surfaces(it_all).v(7,:,:)=zeros(size(wwiline))   ;%ww_x ==0
    tdata2d_uvww.surfaces(it_all).v(8,:,:)=zeros(size(wwiiine))   ;%ww_y ==0
    tdata2d_uvww.surfaces(it_all).v(9,:,:)=wwiline                ;%ww_z ==ww 
    tdata2d_uvww.surfaces(it_all).solutiontime=htime_all(it_all); 

    %3D plot showing w and w_fvcom
    tdata3d_w.surfaces(it_all).zonename=['w-pseudo transect at time step ' num2str(it_all)];
    tdata3d_w.surfaces(it_all).x=xiline2D_level;
    tdata3d_w.surfaces(it_all).y=yiline2D_level;
    tdata3d_w.surfaces(it_all).z=ziline2D_level;
    tdata3d_w.surfaces(it_all).order=2;      %surface defined on x, z
    tdata3d_w.surfaces(it_all).v(1,:,:)=zeros(size(wiline))    ;%w_x ==0
    tdata3d_w.surfaces(it_all).v(2,:,:)=zeros(size(wiline))    ;%w_y ==0
    tdata3d_w.surfaces(it_all).v(3,:,:)=w                      ;%w_z ==w
    tdata3d_w.surfaces(it_all).v(4,:,:)=zeros(size(wiline))    ;%w_fvcom_x==0
    tdata3d_w.surfaces(it_all).v(5,:,:)=zeros(size(wiline))    ;%w_fvcom_y==0
    tdata3d_w.surfaces(it_all).v(6,:,:)=w_fvcom                ;%w_fvcom_z==w_fvcom
    tdata3d_w.surfaces(it_all).v(7,:,:)=ziline2D_level;          %height to color the vertical transect
    tdata2d_w.surfaces(it_all).solutiontime=htime_all(it_all); 
   
    tdata3d_w.FEsurfaces(it_all).zonename=['FVCOM grid at time step ' num2str(it_all)];
    tdata3d_w.FEsurfaces(it_all).x=x_n;
    tdata3d_w.FEsurfaces(it_all).y=y_n;
    tdata3d_w.FEsurfaces(it_all).z=-h_n;
    tdata3d_w.FEsurfaces(it_all).order=3;    %surface defined on x, y
    tdata3d_w.FEsurfaces(it_all).e2n=e2n(:,2:4); 
    tdata3d_w.FEsurfaces(it_all).v(1,:)=zeros(size(h_n'));
    tdata3d_w.FEsurfaces(it_all).v(2,:)=zeros(size(h_n'));
    tdata3d_w.FEsurfaces(it_all).v(3,:)=zeros(size(h_n'));
    tdata3d_w.FEsurfaces(it_all).v(4,:)=zeros(size(h_n'));
    tdata3d_w.FEsurfaces(it_all).v(5,:)=zeros(size(h_n'));
    tdata3d_w.FEsurfaces(it_all).v(6,:)=zeros(size(h_n'));
    tdata3d_w.FEsurfaces(it_all).v(7,:)=-h_n; 
    tdata3d_w.FEsurfaces(it_all).solutiontime=htime_all(it_all); 

    %3D plot showing u,v,ww
    tdata3d_uvww.surfaces(it_all).zonename=['u,v,ww transect at time step ' num2str(it_all)];
    tdata3d_uvww.surfaces(it_all).x        =xiline2D_layer;
    tdata3d_uvww.surfaces(it_all).y        =yiline2D_layer;
    tdata3d_uvww.surfaces(it_all).z        =ziline2D_layer;
    tdata3d_uvww.surfaces(it_all).order    =2;   %surface defined on x, z 
    tdata3d_uvww.surfaces(it_all).v(1,:,:) =uiline                 ;%u_x ==u
    tdata3d_uvww.surfaces(it_all).v(2,:,:) =zeros(size(uiline))    ;%u_y ==0
    tdata3d_uvww.surfaces(it_all).v(3,:,:) =zeros(size(uiline))    ;%u_z ==0
    tdata3d_uvww.surfaces(it_all).v(4,:,:) =zeros(size(viline))    ;%v_x ==0
    tdata3d_uvww.surfaces(it_all).v(5,:,:) =viline                 ;%v_y ==v
    tdata3d_uvww.surfaces(it_all).v(6,:,:) =zeros(size(viline))    ;%v_z ==0
    tdata3d_uvww.surfaces(it_all).v(7,:,:) =zeros(size(wwiline))   ;%ww_x ==0
    tdata3d_uvww.surfaces(it_all).v(8,:,:) =zeros(size(wwiiine))   ;%ww_y ==0
    tdata3d_uvww.surfaces(it_all).v(9,:,:) =wwiline                ;%ww_z ==ww
    tdata3d_uvww.surfaces(it_all).v(10,:,:)=ziline2D_layer;        ;%height used to color the vertical transect (vertical surface)
    tdata3d_uvww.surfaces(it_all).solutiontime=htime_all(it_all); 

    tdata3d_uvww.FEsurfaces(it_all).zonename=['FVCOM grid at time step ' num2str(it_all)];
    tdata3d_uvww.FEsurfaces(it_all).x      =x_n;
    tdata3d_uvww.FEsurfaces(it_all).y      =y_n;
    tdata3d_uvww.FEsurfaces(it_all).z      =-h_n;
    tdata3d_uvww.FEsurfaces(it_all).order  =3;    %surface defined on x, y
    tdata3d_uvww.FEsurfaces(it_all).e2n    =e2n(:,2:4);
    tdata3d.uvww.FEsurfaces(it_all).v(1,:) =zeros(size(h_n'));  %fake variables that do not exist
    tdata3d.uvww.FEsurfaces(it_all).v(2,:) =zeros(size(h_n'));  %on the bottom of the ocean
    tdata3d.uvww.FEsurfaces(it_all).v(3,:) =zeros(size(h_n'));
    tdata3d.uvww.FEsurfaces(it_all).v(4,:) =zeros(size(h_n'));
    tdata3d.uvww.FEsurfaces(it_all).v(5,:) =zeros(size(h_n'));
    tdata3d.uvww.FEsurfaces(it_all).v(6,:) =zeros(size(h_n'));
    tdata3d.uvww.FEsurfaces(it_all).v(7,:) =zeros(size(h_n'));
    tdata3d.uvww.FEsurfaces(it_all).v(8,:) =zeros(size(h_n'));
    tdata3d.uvww.FEsurfaces(it_all).v(9,:) =zeros(size(h_n'));
    tdata3d.uvww.FEsurfaces(it_all).v(10,:)=-h_n; 
    tdata3d.uvww.FEsurfaces(it_all).solutiontime=htime_all(it_all); 

 end

 %file names for tecplot
 tecfile2d_w    ='transect2d_w.plt';
 tecfile2d_uvww ='transect2d_uvww.plt';
 tecfile3d_w    ='transect3d_w.plt';
 tecfile3d_uvww ='transect3d_uvww.plt';
 
 %write the tecplot data 
 mat2tecplot(tdata2d_w,tecfile2d_w);
 mat2tecplot(tdata2d_uvww,tecfile2d_uvww);
 mat2tecplot(tdata3d_w,tecfile3d_w);
 mat2tecplot(tdata3d_uvww,tecfile3d_uvww);

 %


