%Script to plot Transects of DOC and CDOM absorbance
%in the RHode RIver Using data collected across multiple years
%and ebb tides that are decreasing master than the median ebb tide flow
% of -0.0216 m/hr

% to compare to data from the RR and also do visual comparisons with
% previous research done at the Kirkpatrick Marsh and RR

% B Clark, UMCES, March 2016

% Updated B Clark, UMCES, January 2019

%make sure this is run in the directory with the netcdf output file from
%the ICM model

%Load in a file that has the area (m^2) surrounding each node defined
% [node_area]=importdata('*/RhodeFVCOM_2005/Data/grid/node_CVs_area.dat');
node_file=input('What is the full path and file name to the Node Area file? --->  ','s');
[node_area]=importdata(node_file);

%Load in a file that has the DOC concentration and absorbance data from the
% load('*/RhodeFVCOM_2005/Data/RR_doc_transects.mat');% transect data from RR
Transect_file=input('What is the full path and file name to the Rhode River DOC transect data? -->  ','s');
load(Transect_file);

%load in the grid definition file

grid_file=input('What is the full path and file name to the grid definition *.mat file? --> ','s');
load(grid_file);
%
%load in the csv file that has the spectral optical characteristics
% ICM_input=importdata('*/RhodeFVCOM_2005/Data/wqm_kei_photoDEG_v8.csv');
optics_file=input('What is the full path and file name to the spectrally defined optics? --->   ','s');
ICM_input=importdata(optics_file);
lambda=ICM_input(:,1);
astar_cdom1 = ICM_input(:,3);
astar_cdom2 = ICM_input(:,4);
astar_cdom3 = ICM_input(:,5);

%find where 440 nm is
id_440=find(lambda==440);

first_day=input('What is the first day of the model, i.e. April 1 is 91 -->  ');

outdir='transects';
mkdir(outdir);
%%
%cpb_stations = [7355;7224;7196;7129;7071;6984;6875;6744;6434;5856];%7656;

%
nc=netcdf('rhd_0001.nc');

%find the model nodes closest to the transect coordinates
cpb_stations=dsearchn([lld_n(:,3) lld_n(:,2)],ts1_coords);

%     LPOC= nc{'LPOC'}(:);
%     RPOC= nc{'RPOC'}(:);
%     ALG1=nc{'B1'}(:);
%     ALG2=nc{'B2'}(:);
%
%load in all DOC data from the model
WC_CDOC1=nc{'WC_CDOC1'}(:);
disp('Loaded CDOC1');
WC_NCDOC1=nc{'WC_NCDOC1'}(:);
disp('Loaded NCDOC1');
WC_CDOC2=nc{'WC_CDOC2'}(:);
disp('Loaded CDOC2');
WC_NCDOC2=nc{'WC_NCDOC2'}(:);
disp('Loaded NCDOC2');
WC_CDOC3=nc{'WC_CDOC3'}(:);
disp('Loaded CDOC3');
WC_NCDOC3=nc{'WC_NCDOC3'}(:);
disp('Loaded NCDOC3');

% 7540 is the node at the middle of the creek, get the depth
creek_depth_in= nc{'depth'}(:,7540);
%set up the time
mytime=nc{'time'}(:);
ICM_int=(mytime(2)-mytime(1))./86400;
% get the model time, in days
time_vector=(mytime./86400)+first_day;

%find all ebbtides
%get the tidal velocity
dzdx=diff(creek_depth_in);
% now find all ebbtides
ebbtides_in=find(dzdx<0.0);
%%
% set up a transect with distance increasing from the marsh
%get some plotting coordinates for our surface analysis
xplot = [0 (cumsum(sqrt(diff(xyd_n(cpb_stations,2)).^2 +diff(xyd_n(cpb_stations,3)).^2)))'];
%find all ebbtides less than the median ebbtide
ebbtide=find(diff(creek_depth_in)<median(dzdx(ebbtides_in)));
%get total DOC
tDOC=WC_CDOC1+WC_NCDOC1+WC_CDOC2+WC_NCDOC2+WC_CDOC3+WC_NCDOC3;

station_doc = squeeze(tDOC(:,2,cpb_stations));

%the first day of the model run);
%find the Times for July
july_ints = find((time_vector)>182 & (time_vector)< 212);%

july_ints_ebbtide=intersect(july_ints,ebbtide);

gcf=subplot(2,1,1);

% now plot the model DOC and the transect data
plot(xplot,station_doc(july_ints_ebbtide,:),'k');
hold on
plot(xplot,TS1_DOC,'cd','markersize',10)
xlabel('Distance (m)');ylabel('DOC (mg C l^-^1)');

a440_1= astar_cdom1(id_440).*(WC_CDOC1(july_ints_ebbtide,2,cpb_stations));
a440_2= astar_cdom2(id_440).*(WC_CDOC2(july_ints_ebbtide,2,cpb_stations));
a440_3= astar_cdom3(id_440).*(WC_CDOC3(july_ints_ebbtide,2,cpb_stations));
%
%now plot the absorbance at 440 nm (m)
subplot(2,1,2)
plot(xplot,squeeze(a440_1+a440_2+a440_3),'k')
hold on;
plot(xplot,ts1_a440_avg ,'cd','markersize',10)
xlabel('Distance (m)');ylabel('a440 nm (m^-^1)')

% save the files
saveas(gcf,[outdir '/DOC_transect2','.png']);
saveas(gcf,[outdir '/DOC_transect2','.fig']);
saveas(gcf,[outdir '/DOC_transect2','.eps']);



%Finished



