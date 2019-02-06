function [plot1]=velocity_vectors(cpbmatfile,ncfile,nt,parname,layer)
%simple function to plot velocity vectors
%reads in grid file  and
%over a colored map of a given parameter 
%have grid variables in workspace

%Blake Clark UMCES Oct 2014

% gridfile ---- string, sms   **.2dm grid file
% filename ---- string,   netcdf file name, e.g. 'psm_0012.nc'
%  nt       ---- integer,  time record to choose, if nt=0, then we will plot a non-time changing field
%  parname  ---- string,   scalar variable name, e.g. 'temp','salinity' etc, if 'uvmag', then it will
%                          try to get sqrt(u^2+v^2); if 'uv', then it will try to get u,v both in tecplot file
%                          in addition to 'uvmag'
%  
%  layer    ---- integer,  layer to choose to plot, e.g. 10 for bottom, 1 for
%                          surface, if it is 0, then we will plot a 2D array such as depth or surface elevation

%sms2dm2mat(gridfile,'cpb_grd_v1.mat',0,'18-S')
load(cpbmatfile)

u=ncread(ncfile,'u');
v=ncread(ncfile,'v');

show_fvcom_horizontal(ncfile,nt,parname,layer)
hold on
u=u(:,layer,nt);
v=v(:,layer,nt);
quivers(xy_e(:,1),xy_e(:,2),u,v,2,1,'m/s','r')

end