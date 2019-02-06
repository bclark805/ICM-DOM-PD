function [h]=show_fvcom_horizontal(filename,nt,parname,layer,cmin,cmax,tecfile)
%function [h]=show_fvcom_horizontal(filename,nt,parname,layer [,cmin,cmax, tecfile])
%Simple function to show scalar variables at a given layer from fvcom
%netcdf history outputs
%
%filename ---- string,   netcdf file name, e.g. 'psm_0012.nc'
%nt       ---- integer,  time record to choose, if nt=0, then we will plot a non-time changing field
%parname  ---- string,   scalar variable name, e.g. 'temp','salinity' etc, if 'uvmag', then it will
%                        try to get sqrt(u^2+v^2); if 'uv', then it will try to get u,v both in tecplot file
%                        in addition to 'uvmag'
%
%layer    ---- integer,  layer to choose to plot, e.g. 10 for bottom, 1 for
%                        surface, if it is 0, then we will plot a 2D array such as depth or surface elevation
%cmin     ---- minimum value to show in color (optional)
%cmax     ---- maximum value to show in color (optional)
%tecfile  ---- a tecplot file name to also output the picture to tecplot format (optional)
%h        ---- returned handle to the figure 
%
% when cmin, cmax not given, we will pick them up from data 
%                            automatically
%
%e.g. :   
%       show_fvcom_horizontal('psm_0012.nc',20,'salinity',1)
%
% Wen Long, PNNL,  wen.long@pnnl.gov      02/03/2012
%

if nargin<5
   cmin=[]; 
end
if nargin<6
   cmax=[];
end

if(nargin==7)
   do_tecplot=1; 
else
   do_tecplot=0;
   tecfile=[]; 
end

if(~exist(filename,'file'))
    error('Oops can not find file to open');
end

hisnc=netcdf(filename);   %open netcdf

%get number of sigma layers
siglay=hisnc{'siglay'}(:);
Nz=length(siglay);
if layer< 0 || layer >Nz
   error(['Oops, layer is out of bounds [0,' int2str(Nz) ']']);
end

%get number of time steps
time=hisnc{'time'}(:);
Nt=length(time);
if nt <0 || nt> Nt
   error(['Oops, time record # nt is out of bounds [0,' int2str(Nt) ']']);
end

lon=hisnc{'lon'}(:);
lat=hisnc{'lat'}(:);
nv=hisnc{'nv'}(:)';            %get element node indices triplets

if(layer>=1)  %has vertical dimension
   if(nt>=1)    %has time dimension
       switch(parname)  
            case {'uvmag','uv'}  %get u, v components and then calculate magnitude
               uu=hisnc{'u'}(nt,layer,:);
	       vv=hisnc{'v'}(nt,layer,:);
	       v=sqrt( uu.^2+vv.^2);
            otherwise
               v=hisnc{parname}(nt,layer,:);  %fetch the nt'th record, layer
       end
   else         %no time dimension

       switch(parname)
           case {'uvmag','uv'}  %get magnitude of current u and v, uvmag=sqrt(u^2+v^2)
               uu=hisnc{'u'}(layer,:);
               vv=hisnc{'v'}(layer,:);
               v=sqrt(uu.^2+vv.^2);
           otherwise
               v=hisnc{parname}(layer,:);  %no time dimension
       end
   end
else          %no vertical dimension
   if(nt>=1)     %has time dimension
      switch (parname)   %get magnitude of current u and v, uvmag=sqrt(u^2+v^2)
          case {'uvmag','uv'}
              uu=hisnc{'u'}(nt,:);
              vv=hisnc{'v'}(nt,:);
              v=sqrt( uu.^2+vv.^2);
     	  otherwise
              v=hisnc{parname}(nt,:);  
      end
   else          %no time dimension
      switch(parname)  %get magnitude of current u and v, uvmag=sqrt(u^2+v^2)
	case {'uvmag','uv'}
          uu=hisnc{'u'}(:);
          vv=hisnc{'v'}(:);
          v=sqrt(uu.^2+vv.^2);
	otherwise
          v=hisnc{parname}(:);  %no time dimension
      end
   end
end

if isempty(cmin) || isempty(cmax)
   hfig=tricolor(nv,lon,lat,squeeze(v));  %plot tricolor
else
   hfig=tricolor(nv,lon,lat,squeeze(v),cmin,cmax);
end

%hold on;
%trimesh(nv,lon,lat,'color','k');
hold off;

h=hfig;  %return handle of the figure



%generate tecplot of the picture
if(do_tecplot || ~isempty(tecfile))
    tdata=[];
    if(layer>=1)
        tdata.Nvar=3+1;  %first 3 are x,y,z, and last one is the parameter
        tdata.varnames=[{'x','y','z'},parname];
    else
        tdata.NVar=3+1;    %first 3 are x,y,z and the last one is the parameter
                             %however since no z dimension, z will be filled out with zeros
        tdata.varnames=[{'x','y','z'},parname];
    end

    %to have vector u, v 
    switch(parname)
        case {'uvmag','uv'}
            tdata.Nvar=3+3;
	    tdata.varnames=[{'x','y','z'},'u','v','uvmag'];
        otherwise

	    %do nothing
    end


    if(nt>=1)
        tdata.FEsurfaces(1).zonename=[filename ' record number ' num2str(nt)]; 
    else
        tdata.FEsurfaces(1).zonename=[filename ' ' parname]; 
    end
    tdata.FEsurfaces(1).x=lon;
    tdata.FEsurfaces(1).y=lat; 

    if(layer>=1)  %has vertical dimension  give true z coordinates
       if(nt>=1)        %get surface elevation
        zeta=hisnc{'zeta'}(nt,:); 

        %give error if not able to have sane zeta values
        %... to do ...

       else             %if no time dimension for the parameter parname, then use zero as surface elevation
        zeta=zeros(size(lon)); 
       end

       %get depth
        depth=hisnc{'h'}(:); 
       %give error if not able to have sane h values
       % ... to do ...

       %get layer elevation (layer is layer number) 
       %siglay is sigma coordinate value of range [-1,0], when siglay= -1, zcoord=-depth -zeta+zeta   = -depth --> bottom
       %                                                  when siglay= 0,  zcoord=0*(depth+zeta)+zeta = zeta   --> surface
        zcoord=siglay(layer)*(depth(:)+zeta(:))+zeta(:); 

        tdata.FEsurfaces(1).z=zcoord(:); 

    else  %no vertical dimension, give dummy z coordinates

        tdata.FEsurfaces(1).z=zeros([length(lon(:)),1]);   %zeros of size [NN,1], where NN is number of nodes
    end

    tdata.FEsurfaces(1).e2n=nv;

    if(length(v(:))==size(nv,1))   %variable is on element
        tdata.FEsurfaces(1).varloc=1;   %data is at cell center (element center)
    else
        tdata.FEsurfaces(1).varloc=0;   %data is nodal
    end

    switch(parname)

	case {'uvmag','uv'}
            tdata.FEsurfaces(1).v(1,:)=uu(:)';   %u component
            tdata.FEsurfaces(1).v(2,:)=vv(:)';   %v component
            tdata.FEsurfaces(1).v(3,:)=v(:)';    %velocity magnitude

	otherwise

            tdata.FEsurfaces(1).v(1,:)=v(:)';  
    end

    if(layer>=1)
       tdata.FEsurfaces(1).order=4;  %order 4, surface is defined on (x,y,z) coordinates, i.e. v=f(x,y,z)
    else
       tdata.FEsurfaces(1).order=3;  %order 3, 2D surface defined on (x,y), i.e. z=z(x,y), v=v(x,y)
    end

    if(nt>=1)  %has time dimension then give solution time to the tecplot zone
        tdata.FEsurfaces(1).solutiontime=time(nt)/86400;   %get time in days
    else
        tdata.FEsurfaces(1).solutiontime=[]; 
    end

    %output to the tecplot file

    mat2tecplot(tdata,tecfile); 

end

close(hisnc);
