function [c,h]=show_fvcom_contour(filename,nt,parname,layer,cval);
%function [c,h]=show_fvcom_contour(filename,nt,parname,layer[,cval] );
%Simple function to show scalar variables at a given layer from fvcoms
%
%netcdf history outputs
%
%filename ---- string,   netcdf file name, e.g. 'psm_0012.nc'
%nt       ---- integer,  time record to choose
%parname  ---- string,   scalar variable name, e.g. 'temp','salinity' etc
%layer    ---- integer,  layer to choose to plot, e.g. 10 for bottom, 1 for
%                        surface
%cval     ---- vector: one-d array of contour values to show,
%         ---- N :number of contours to show (default is 10)
%
%c        ---- output handle of contours
%h        ---- returned handle to the figure 
%
% when cmin, cmax not given, we will pick them up from data 
%                            automatically
%e.g. :   
%       show_fvcom_contour('psm_0012.nc',20,'salinity',1,(1:1:30))
%
% Wen Long, PNNL,  wen.long@pnnl.gov      02/03/2012
%

if nargin<5
   cval=10;   %number of contours to show by default is 10
end

if(~exist(filename,'file'))
    error('Oops can not find file to open');
end

hisnc=netcdf(filename);   %open netcdf

%get number of sigma layers
siglay=hisnc{'siglay'}(:);
Nz=length(siglay);
if layer< 1 || layer >Nz
   error(['Oops, layer is out of bounds [1,' int2str(Nz) ']']);
end

%get number of time steps
time=hisnc{'time'}(:);
Nt=length(time);
if nt <1 || nt> Nt
   error(['Oops, time record # nt is out of bounds [1,' int2str(Nt) ']']);
end

lon=hisnc{'lon'}(:);
lat=hisnc{'lat'}(:);
nv=hisnc{'nv'}(:)';            %get element node indices triplets
v=hisnc{parname}(nt,layer,:);  %fetch the nt'th record, layer
%[c, h]=tricontour(nv,[lon(:) lat(:)],squeeze(v(:)),cval);  %plot tricontour

 %[htmp,h,c]=triplot2([lon(:) lat(:) v(:)],nv, v(:), 'surface',cval);

   whos lon lat v nv v cval 
  [htmp,h,c]=triplot2([lon(:) lat(:) v(:)],nv, v(:), 'contour',cval);

 
%hold on;
%trimesh(nv,lon,lat,'color','k');
hold off;
close(hisnc);
