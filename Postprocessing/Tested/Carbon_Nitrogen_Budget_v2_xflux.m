%Carbon Budget calculation for ICM with netcdf outputs

% B Clark UMCES HPL Nov 2018
% Updated B Clark UMCES, Jan 2019
%
%ALl non-concentration variables in the budget (e.g. right had side of
%equation)
% are integrated internally in ICM and loaded from the hisdata*.mat files

%====**********
% First run the wqmhist2matlab_budget_v2 script to get the hisdata*.mat
% files that extract and cut up all of the integrated budget terms
%

%load in the files that define the nodes in the Tributary
%the river discharge file
%and the watershed inputs of the water quality variables
load('/discover/nobackup/projects/yukonriver/RhodeRiver/Water_Quality/CBAY_DATA/RhodeRiver_nodes.mat');
% load(input('What is the full path and file name to Tributary Node Information? --> ','s'));
 load('/discover/nobackup/projects/yukonriver/RhodeRiver/Water_Quality/CBAY_DATA/RR_discharge_2005.mat');
% load(input('What is the full path and file name to the River discharge file? --> ','s'));
load('/discover/nobackup/projects/yukonriver/RhodeRiver/Water_Quality/CBAY_DATA/RR_Watershed_WQ_2005.mat');
% load(input('What is the full path and file name to the Watershed Water Quality Input data? --> ','s'));

%load the sediment definition file where marsh nodes are flagged with 1
marsh_file= '/discover/nobackup/projects/yukonriver/RhodeRiver/Water_Quality/CBAY_DATA/seddom_grid.dat';
% marsh_file=input('What is the full path and file name to the sediment definition file? --> ','s');
marsh_data=importdata(marsh_file);
% load the grid file
grid_file='/discover/nobackup/projects/yukonriver/RhodeRiver/Water_Quality/CBAY_DATA/rhoderiver_grid_v1_2.mat';
% grid_file=input('What is the full path and file name to the grid definition *.mat file? --> ','s');
load(grid_file);
%load the node area file
[node_area]=importdata('/discover/nobackup/projects/yukonriver/RhodeRiver/Water_Quality/CBAY_DATA/node_CVS_area.dat');%get our nodal area for marsh flux calcs;
% node_file=input('What is the full path and file name to the Node Area file? --->  ','s');
% [node_area]=importdata(node_file);
%load the file with the different polygons from the Rhode River
%each polygon has a different set of nodes with a decreasing area of the
%Rhode River
load('/discover/nobackup/projects/yukonriver/RhodeRiver/Water_Quality/CBAY_DATA/polygon_work.mat')
% load(input('What is the full path and file name of the polygon definition file? --> ','s'));

% first_day =input('What is the starting day of the year for this run (i.e. April 1 is 91 ---> ');
%
% Use this flag depending what region to set up the budget for, 1 being the
% whole tributary decreasing size of estuary area to 5
nodes_flag=input('What polygon do you want to build a budget for? 1-5 for trib, 0 for mainstem --->  ')
%get the rest of the information
wqm_hisdir = input('What is full path to the WQM history directory? --->  ','s');
% numfile=input('How many files do you want to create a budget for? ---> ')

% nameseg=input('What is the name segment for the FVCOM files ---> ','s');
% nameseg='rhd'
%SET UP EMPTY ARRAYS FOR ALL THE BUDGET VARIABLES

fname=[wqm_hisdir '/netcdf/rhd_0001.nc']

%pull out the vertical coordinate distribution
%this is the FVCOM hydrodynamics file, used to load in the vertical
%coordinate system

nc=netcdf(fname);
% s_hybrid=nc{'siglev'}(:);
% dz=diff(s_hybrid,1);
% for i = 1:nlayers;
%     layer_thickness(i,:)=dz(i,:).*xyd_n(:,4)'*-1;
% end

myinfo=ncinfo(fname);

ngrids=myinfo.Dimensions(2).Length;
nlayers=myinfo.Dimensions(4).Length
ntimes=myinfo.Dimensions(end).Length;

%%
%how many days was the model run for?

outdir=['budget_plots_poly' num2str(nodes_flag)];

mkdir(outdir);
%
%load marsh definition file and pull out the marsh nodes
%and the area at each marsh node
marsh.nodes= find(marsh_data(:,3)>0);
marsh.area=node_area(marsh.nodes,2);

%this pulls out the nodes we want, depending on the polygon
if (~nodes_flag)
    estuary.nodes=1:size(lld_n);%RR_nodes_out(find(lld_n(RR_nodes_out,3) > 38.875));
    river_id=8;
elseif(nodes_flag==1)
    estuary.nodes=RR_nodes_out(find(lld_n(RR_nodes_out,3) > 38.872));
    river_id=8;
elseif(nodes_flag==2)
    estuary.nodes=poly1_nodes;
    river_id=7;
elseif(nodes_flag==3)  
    estuary.nodes=poly2_nodes;
    river_id=6;
elseif(nodes_flag==4)  
    estuary.nodes=poly3_nodes;
    river_id=4;
elseif(nodes_flag==5)  
    estuary.nodes=poly4_nodes;
    river_id=3;    
end
%

% this is the entire area including the marsh

estuary.nodes=estuary.nodes;
estuary.area=node_area(estuary.nodes,2);

total_estuary_area=sum(estuary.area);
total_marsh_area=sum(marsh.area);

% get the nodes

% this is only the subtidal area, excluding the marsh
%use in the sediment flux calculations
estuary.nodesonly=setdiff(estuary.nodes,marsh.nodes);
estuary.areaonly=node_area(estuary.nodesonly,2);

total_estuary_areaonly=sum(estuary.areaonly);

%plot the polygons
tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3),lld_n(:,4));
hold on

plot(lld_n(marsh.nodes,2),lld_n(marsh.nodes,3),'ro')
plot(lld_n(estuary.nodes,2),lld_n(estuary.nodes,3),'kd')
saveas(gcf,[outdir '/Marsh_map.png'])
%
%
% zLayer=diff(s_hybrid,1)*-1;% the amount of water column for each layer, adds up to 1;
tic
disp('Getting Depth and Time Information');
% nc = netcdf([fname]);
DEPTH=zeros(ntimes,ngrids);
% DEPTH=nc{'depth'}(:);
time_vector_in=nc{'time'}(:);
time_vector=(time_vector_in./86400);
time_vector_hr=time_vector.*24;

disp('Got time and depth data')
toc
ICM_int=diff(time_vector_in(1:2))/86400;%get the output time period in days
ndays=ntimes.*ICM_int;
%%
tic
disp('getting river data')
%
%algal fractionation, from the input file
%change this if it is changed in the wqm_algae file
FCD11=1.0;FCD12=0.1;FCD13=0.05
FCD21=1.0;FCD22=0.1;FCD23=0.05
FNCD11=0.3;FNCD12=0.45;FNCD13=0.0;
FNCD21=0.3;FNCD22=0.45;FNCD23=0.0;
FCND11=0.2;FCND12=0.05;FCND13=0.0;
FCND21=0.2;FCND22=0.05;FCND23=0.0;
FNCND11=0.6;FNCND12=0.15;FNCND13=0.0;
FNCND21=0.6;FNCND22=0.15;FNCND23=0.0;
% FCND11=0.2;FCND12=0.05;FCND13=0.0;
% FCND21=0.2;FCND22=0.05;FCND23=0.0;
% FNCND11=0.6;FNCND12=0.15;FNCND13=0.0;
% FNCND21=0.6;FNCND22=0.15;FNCND23=0.0;
%Riverine_budget

%total algal fractions
FCD1=0.4;FCD2=0.55;FCD3=0.05;
FND1=0.4;FND2=0.55;FND3=0.05;

%hydrolysis fractionation;
FHDRC1=0.35;FHDRC2=0.60;FHDRC3=0.05;
FHDRN1=0.35;FHDRN1=0.6;FHDRN3=0.05;

%fractionation of photobleaching
FPDC31 =0.5;  FPDC32=0.25 ;  FPDC33=0.25;

FPDC21 =0.5 ; FPDC22 =0.25;  FPDC23=0.25;

FPDC11=0.5;   FPDC12 =0.25;  FPDC13=0.25;

% get the discharge data for the time period and interpolate to the
% discharge frequency of 1 hour for easier plotting and manipulation
% first_id=find((91+3)/ICM_int==length(mydischarge));
river_time=(1:length(mydischarge))-24*3;

%match up the model times with the discharge times
matching_times=dsearchn(river_time',time_vector_hr);%
%get the discharge times at the model time frequency for each river
cut_discharge=mydischarge(matching_times,1:river_id);
% discharge_Myint=(cut_discharge(2:round(ICM_int*24):end,:));
discharge_Myint=cut_discharge;
river.times=Final_RR_forcing(1:45:end,2);
%
Final_RR_forcing=Final_RR_forcing(:,1:river_id);

% will take the concentration and multiple my discharge g m^-3 x m^3 s^-1 *
% 3600 secs hr^-1 then integrate to get kg
river.CDOC1.all=interp1(river.times,Final_RR_forcing(10:45:end,:),time_vector_hr,'linear')
river.NCDOC1.all=interp1(river.times,Final_RR_forcing(11:45:end,:),time_vector_hr,'linear');
river.LPOC.all=interp1(river.times,Final_RR_forcing(12:45:end,:),time_vector_hr,'linear');
river.RPOC.all=interp1(river.times,Final_RR_forcing(13:45:end,:),time_vector_hr,'linear');
river.CDOC2.all=interp1(river.times,Final_RR_forcing(34:45:end,:),time_vector_hr,'linear');
river.NCDOC2.all=interp1(river.times,Final_RR_forcing(35:45:end,:),time_vector_hr,'linear');
river.CDOC3.all=interp1(river.times,Final_RR_forcing(36:45:end,:),time_vector_hr,'linear');
river.NCDOC3.all=interp1(river.times,Final_RR_forcing(37:45:end,:),time_vector_hr,'linear');

% now integrate and get one time series for each constituent
river.CDOC1.ts=nansum(river.CDOC1.all.*discharge_Myint.*86400.*ICM_int./1000,2);% kg per ICM_int period
river.CDOC1.int=nansum(river.CDOC1.ts);
river.NCDOC1.ts=nansum(river.NCDOC1.all.*discharge_Myint.*86400.*ICM_int./1000,2);% kg per ICM_int period
river.NCDOC1.int=nansum(river.NCDOC1.ts);
river.NCDOC2.ts=nansum(river.NCDOC2.all.*discharge_Myint.*86400.*ICM_int./1000,2);% kg per ICM_int period
river.NCDOC2.int=nansum(river.NCDOC2.ts);
river.CDOC2.ts=nansum(river.CDOC2.all.*discharge_Myint.*86400.*ICM_int./1000,2);% kg per ICM_int period
river.CDOC2.int=nansum(river.CDOC2.ts);
river.NCDOC3.ts=nansum(river.NCDOC3.all.*discharge_Myint.*86400.*ICM_int./1000,2);% kg per ICM_int period
river.NCDOC3.int=nansum(river.NCDOC3.ts);
river.CDOC3.ts=nansum(river.CDOC3.all.*discharge_Myint.*86400.*ICM_int./1000,2);% kg per ICM_int period
river.CDOC3.int=nansum(river.CDOC3.ts);

river.LPOC.ts=nansum(river.LPOC.all.*discharge_Myint.*86400.*ICM_int./1000,2);% kg per ICM_int period
river.LPOC.int=nansum(river.LPOC.ts);
river.RPOC.ts=nansum(river.RPOC.all.*discharge_Myint.*86400.*ICM_int./1000,2);% kg per ICM_int period
river.RPOC.int=nansum(river.RPOC.ts);

river.TDOC.ts=river.CDOC1.ts+river.CDOC2.ts+river.CDOC3.ts+river.NCDOC1.ts+river.NCDOC2.ts+river.NCDOC3.ts;
river.TDOC.int=river.CDOC1.int+river.CDOC2.int+river.CDOC3.int+river.NCDOC1.int+river.NCDOC2.int+river.NCDOC3.int;

river.POC.ts=river.LPOC.ts+river.RPOC.ts;
river.POC.int=river.LPOC.int+river.RPOC.int;

%Nitrogen now
% will take the concentration and multiple my discharge g m^-3 x m^3 s^-1 *
% 3600 secs hr^-1 then integrate to get kg
river.CDON1.all=interp1(river.times,Final_RR_forcing(17:45:end,:),time_vector_hr,'linear','extrap');
river.NCDON1.all=interp1(river.times,Final_RR_forcing(18:45:end,:),time_vector_hr,'linear','extrap');
river.LPON.all=interp1(river.times,Final_RR_forcing(19:45:end,:),time_vector_hr,'linear','extrap');
river.RPON.all=interp1(river.times,Final_RR_forcing(20:45:end,:),time_vector_hr,'linear','extrap');
river.CDON2.all=interp1(river.times,Final_RR_forcing(38:45:end,:),time_vector_hr,'linear','extrap');
river.NCDON2.all=interp1(river.times,Final_RR_forcing(39:45:end,:),time_vector_hr,'linear','extrap');
river.CDON3.all=interp1(river.times,Final_RR_forcing(40:45:end,:),time_vector_hr,'linear','extrap');
river.NCDON3.all=interp1(river.times,Final_RR_forcing(41:45:end,:),time_vector_hr,'linear','extrap');

river.NH4.all=interp1(river.times,Final_RR_forcing(14:45:end,:),time_vector_hr,'linear','extrap');
river.NO3.all=interp1(river.times,Final_RR_forcing(15:45:end,:),time_vector_hr,'linear','extrap');

% now integrate and get one time series for each constituent

river.CDON1.ts=nansum(river.CDON1.all.*discharge_Myint.*86400.*ICM_int./1000,2);% kg per ICM_int period
river.CDON1.int=nansum(river.CDON1.ts);
river.NCDON1.ts=nansum(river.NCDON1.all.*discharge_Myint*86400*ICM_int./1000,2);% kg per ICM_int period
river.NCDON1.int=nansum(river.NCDON1.ts);
river.NCDON2.ts=nansum(river.NCDON2.all.*discharge_Myint*86400*ICM_int./1000,2);% kg per ICM_int period
river.NCDON2.int=nansum(river.NCDON2.ts);
river.CDON2.ts=nansum(river.CDON2.all.*discharge_Myint*86400*ICM_int./1000,2);% kg per ICM_int period
river.CDON2.int=nansum(river.CDON2.ts);
river.NCDON3.ts=nansum(river.NCDON3.all.*discharge_Myint*86400*ICM_int./1000,2);% kg per ICM_int period
river.NCDON3.int=nansum(river.NCDON3.ts);
river.CDON3.ts=nansum(river.CDON3.all.*discharge_Myint*86400*ICM_int./1000,2);% kg per ICM_int period
river.CDON3.int=nansum(river.CDON3.ts);

river.LPON.ts=nansum(river.LPON.all.*discharge_Myint*86400*ICM_int./1000,2);% kg per ICM_int period
river.LPON.int=nansum(river.LPON.ts);
river.RPON.ts=nansum(river.RPON.all.*discharge_Myint*86400*ICM_int./1000,2);% kg per ICM_int period
river.RPON.int=nansum(river.RPON.ts);

river.NH4.ts=nansum(river.NH4.all.*discharge_Myint*86400*ICM_int./1000,2);% kg per ICM_int period
river.NH4.int=nansum(river.NH4.ts);
river.NO3.ts=nansum(river.NO3.all.*discharge_Myint*86400*ICM_int./1000,2);% kg per ICM_int period
river.NO3.int=nansum(river.NO3.ts);

river.TDON.ts=river.CDON1.ts+river.CDON2.ts+river.CDON3.ts+river.NCDON1.ts+river.NCDON2.ts+river.NCDON3.ts;
river.TDON.int=river.CDON1.int+river.CDON2.int+river.CDON3.int+river.NCDON1.int+river.NCDON2.int+river.NCDON3.int;

river.PON.ts=river.LPON.ts+river.RPON.ts;

river.PON.int=river.LPON.int+river.RPON.int;

toc
disp('got river data')
%% get TCE data

%get the Flux across the transect data from tceoutput_flux.mat
Read_TCE_FLUX

CDOC1fx=sum(myfluxes(:,:,9),2);
CDOC2fx=sum(myfluxes(:,:,33),2);
CDOC3fx=sum(myfluxes(:,:,35),2);
NCDOC1fx=sum(myfluxes(:,:,10),2);
NCDOC2fx=sum(myfluxes(:,:,34),2);
NCDOC3fx=sum(myfluxes(:,:,36),2);
DOC1fx=CDOC1fx+NCDOC1fx;
DOC2fx=CDOC2fx+NCDOC2fx;
DOC3fx=CDOC3fx+NCDOC3fx;
tDOCfx=DOC1fx+DOC2fx+DOC3fx;

CDON1fx=sum(myfluxes(:,:,16),2);
CDON2fx=sum(myfluxes(:,:,37),2);
CDON3fx=sum(myfluxes(:,:,39),2);
NCDON1fx=sum(myfluxes(:,:,17),2);
NCDON2fx=sum(myfluxes(:,:,38),2);
NCDON3fx=sum(myfluxes(:,:,40),2);
DON1fx=CDON1fx+NCDON1fx;
DON2fx=CDON2fx+NCDON2fx;
DON3fx=CDON3fx+NCDON3fx
tDONfx=DON1fx+DON2fx+DON3fx;

LPOCfx=sum(myfluxes(:,:,11),2);
RPOCfx=sum(myfluxes(:,:,12),2)
POCfx=LPOCfx+RPOCfx;
LPONfx=sum(myfluxes(:,:,18),2);
RPONfx=sum(myfluxes(:,:,19),2)
PONfx=LPONfx+RPONfx;

NH4fx=sum(myfluxes(:,:,13),2);
NO3fx=sum(myfluxes(:,:,14),2);

B1fx=sum(myfluxes(:,:,4),2);
B2fx=sum(myfluxes(:,:,5),2);
%% get the sediment flux data; only need to open the first file

load('hisdata_00001.mat')

%POC
JPOC1=hisdata(end,:,strcmp(varnames,'JPOC1'));
JPOC2=hisdata(end,:,strcmp(varnames,'JPOC2'));
JPOC3=hisdata(end,:,strcmp(varnames,'JPOC3'));

%PON
JPON1=hisdata(end,:,strcmp(varnames,'JPON1'));
JPON2=hisdata(end,:,strcmp(varnames,'JPON2'));
JPON3=hisdata(end,:,strcmp(varnames,'JPON3'));

%DOC
JCDOC1=hisdata(end,:,strcmp(varnames,'JCDOC1'));
JCDOC2=hisdata(end,:,strcmp(varnames,'JCDOC2'));
JCDOC3=hisdata(end,:,strcmp(varnames,'JCDOC3'));
JNCDOC1=hisdata(end,:,strcmp(varnames,'JNCDOC1'));
JNCDOC2=hisdata(end,:,strcmp(varnames,'JNCDOC2'));
JNCDOC3=hisdata(end,:,strcmp(varnames,'JNCDOC3'));

%DON
JCDON1=hisdata(end,:,strcmp(varnames,'JCDON1'));
JCDON2=hisdata(end,:,strcmp(varnames,'JCDON2'));
JCDON3=hisdata(end,:,strcmp(varnames,'JCDON3'));
JNCDON1=hisdata(end,:,strcmp(varnames,'JNCDON1'));
JNCDON2=hisdata(end,:,strcmp(varnames,'JNCDON2'));
JNCDON3=hisdata(end,:,strcmp(varnames,'JNCDON3'));

%NH4
JNH4=hisdata(end,:,strcmp(varnames,'JNH4'));

%NO3
JNO3=hisdata(end,:,strcmp(varnames,'JNO3'));

%% now the water column data for every other budget term
for ilay =  1 : nlayers;
        nlayer =ilay
    
    load(['hisdata_' num2str(nlayer,'%05d') '.mat'])
    
    algdoc(ilay,:)=abs(hisdata(end,:,(strcmp(varnames,'AlgDOC'))));
    algpoc(ilay,:)=abs(hisdata(end,:,(strcmp(varnames,'AlgPOC'))));
    denit(ilay,:)=abs(hisdata(end,:,(strcmp(varnames,'DenitDOC'))));
    hdrlpoc(ilay,:)=abs(hisdata(end,:,(strcmp(varnames,'HDRLPOC'))));
    hdrrpoc(ilay,:)=abs(hisdata(end,:,(strcmp(varnames,'HDRRPOC'))));
    coagc(ilay,:)=abs(hisdata(end,:,(strcmp(varnames,'COAGC'))));
    mnldoc1(ilay,:)=abs(hisdata(end,:,(strcmp(varnames,'MNLDOC1'))));
    mnldoc2(ilay,:)=abs(hisdata(end,:,(strcmp(varnames,'MNLDOC2'))));
    mnldoc3(ilay,:)=abs(hisdata(end,:,(strcmp(varnames,'MNLDOC3'))));
  
    %change in water column concentration of DOC
    DTWCDOC1(ilay,:)=(hisdata(end,:,strcmp(varnames,'DTWCDOC1')));
    DTWCDOC2(ilay,:)=(hisdata(end,:,strcmp(varnames,'DTWCDOC2')));
    DTWCDOC3(ilay,:)=(hisdata(end,:,strcmp(varnames,'DTWCDOC3')));
    DTWNCDOC1(ilay,:)=(hisdata(end,:,strcmp(varnames,'DTWNCDOC1')));
    DTWNCDOC2(ilay,:)=(hisdata(end,:,strcmp(varnames,'DTWNCDOC2')));
    DTWNCDOC3(ilay,:)=(hisdata(end,:,strcmp(varnames,'DTWNCDOC3')));
    algdon(ilay,:)=abs(hisdata(end,:,(strcmp(varnames,'AlgDON'))));
    
    algpon(ilay,:)=abs(hisdata(end,:,(strcmp(varnames,'AlgPON'))));
    algnh4(ilay,:)=abs(hisdata(end,:,(strcmp(varnames,'AlgNH4'))));
    algno3(ilay,:)=abs(hisdata(end,:,(strcmp(varnames,'AlgNO3')))); 
    denno3(ilay,:)=abs(hisdata(end,:,(strcmp(varnames,'DenitNO3'))));
    nt(ilay,:)=abs(hisdata(end,:,(strcmp(varnames,'NitrifNH4'))));
    denitn(ilay,:)=abs(hisdata(end,:,(strcmp(varnames,'DenitDON'))));
    hdrlpon(ilay,:)=abs(hisdata(end,:,(strcmp(varnames,'HDRLPON'))));
    hdrrpon(ilay,:)=abs(hisdata(end,:,(strcmp(varnames,'HDRRPON'))));
    coagn(ilay,:)=abs(hisdata(end,:,(strcmp(varnames,'COAGN'))));
    mnldon1(ilay,:)=abs(hisdata(end,:,(strcmp(varnames,'MNLDON1'))));
    mnldon2(ilay,:)=abs(hisdata(end,:,(strcmp(varnames,'MNLDON2'))));
    mnldon3(ilay,:)=abs(hisdata(end,:,(strcmp(varnames,'MNLDON3'))));
    
     %change in water column concentration of DON
    DTWCDON1(ilay,:)=(hisdata(end,:,strcmp(varnames,'DTWCDON1')));
    DTWCDON2(ilay,:)=(hisdata(end,:,strcmp(varnames,'DTWCDON2')));
    DTWCDON3(ilay,:)=(hisdata(end,:,strcmp(varnames,'DTWCDON3')));
    DTWNCDON1(ilay,:)=(hisdata(end,:,strcmp(varnames,'DTWNCDON1')));
    DTWNCDON2(ilay,:)=(hisdata(end,:,strcmp(varnames,'DTWNCDON2')));
    DTWNCDON3(ilay,:)=(hisdata(end,:,strcmp(varnames,'DTWNCDON3')));
    
    %change in water column concentration of PON
    DTRPON(ilay,:)=(hisdata(end,:,strcmp(varnames,'DTRPON')));
    DTLPON(ilay,:)=(hisdata(end,:,strcmp(varnames,'DTLPON')));
    
    %change in water column concentration of POC
    DTRPOC(ilay,:)=(hisdata(end,:,strcmp(varnames,'DTRPOC')));
    DTLPOC(ilay,:)=(hisdata(end,:,strcmp(varnames,'DTLPOC')));
    
    %change in water column concentration of NH4
    DTNH4(ilay,:)=hisdata(end,:,strcmp(varnames,'DTNH4'));
    %change in water column concentration of NO3
    DTNO3(ilay,:)=hisdata(end,:,strcmp(varnames,'DTNO3'));
  
    %photochemical pathways of DOC
%     dpd32(ilay,:)=hisdata(end,:,:,(strcmp(varnames,'DPD32')));
    dpd31(ilay,:)=abs(hisdata(end,:,(strcmp(varnames,'DPD31'))));
    dpd30(ilay,:)=abs(hisdata(end,:,(strcmp(varnames,'DPD30'))));
    dpd3N(ilay,:)=abs(hisdata(end,:,(strcmp(varnames,'DPD3N'))));
    dpd21(ilay,:)=abs(hisdata(end,:,(strcmp(varnames,'DPD21'))));
    dpd20(ilay,:)=abs(hisdata(end,:,(strcmp(varnames,'DPD20'))));
    dpd2N(ilay,:)=abs(hisdata(end,:,(strcmp(varnames,'DPD2N'))));    
    dpd10(ilay,:)=abs(hisdata(end,:,(strcmp(varnames,'DPD10'))));
    dpd1N(ilay,:)=abs(hisdata(end,:,(strcmp(varnames,'DPD1N')))); 
    
    %photochemical pathways of DON
%     dpd32(ilay,:)=hisdata(end,:,:,(strcmp(varnames,'DPD32')));
    dpd31N(ilay,:)=abs(hisdata(end,:,(strcmp(varnames,'DPD31N'))));
    dpd30N(ilay,:)=abs(hisdata(end,:,(strcmp(varnames,'DPD30N'))));
    dpd3NN(ilay,:)=abs(hisdata(end,:,(strcmp(varnames,'DPD3NN'))));
    dpd21N(ilay,:)=abs(hisdata(end,:,(strcmp(varnames,'DPD21N'))));
    dpd20N(ilay,:)=abs(hisdata(end,:,(strcmp(varnames,'DPD20N'))));
    dpd2NN(ilay,:)=abs(hisdata(end,:,(strcmp(varnames,'DPD2NN'))));    
    dpd10N(ilay,:)=abs(hisdata(end,:,(strcmp(varnames,'DPD10N'))));
    dpd1NN(ilay,:)=abs(hisdata(end,:,(strcmp(varnames,'DPD1NN'))));
      
end


%%
disp('Calculating the Carbon budgets');
% DOC1
%marsh sediment
budgets.DOC1(1).name='Marsh JDOC1';
budgets.DOC1(1).units='kg C';
budgets.DOC1(1).value=sum(JCDOC1(marsh.nodes)+JNCDOC1(marsh.nodes));
budgets.DOC1(1).areal=budgets.DOC1(1).value./total_marsh_area./ndays;%get the areal average per annum
budgets.DOC1(1).type='Source';
%estuary sediment
budgets.DOC1(2).name='Estuary JDOC1';
budgets.DOC1(2).units='kg C';
budgets.DOC1(2).value=sum(JCDOC1(estuary.nodesonly)+JNCDOC1(estuary.nodesonly));
budgets.DOC1(2).areal=budgets.DOC1(1).value./total_estuary_areaonly./ndays;%get the areal average per annum
budgets.DOC1(2).type='Source';
%algae production
budgets.DOC1(3).name='Estuary AlgDOC1';
budgets.DOC1(3).units='kg C';
budgets.DOC1(3).value=sum(sum(algdoc(:,estuary.nodes)))*FCD1;
budgets.DOC1(3).areal=budgets.DOC1(3).value./total_estuary_area./ndays;
budgets.DOC1(3).type='Source';
%hydrolysis of POC --> DOC
budgets.DOC1(4).name='Estuary HdrPOC1';
budgets.DOC1(4).units='kg C';
budgets.DOC1(4).value=sum(sum(hdrlpoc(:,estuary.nodes)+hdrrpoc(:,estuary.nodes)))*FHDRC1;
budgets.DOC1(4).areal=budgets.DOC1(4).value./total_estuary_area./ndays;
budgets.DOC1(4).type='Source';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%denitrification of DOC
budgets.DOC1(5).name='Denitrification of DOC';
budgets.DOC1(5).units='kg C';
budgets.DOC1(5).value=sum(sum(denit(:,estuary.nodes)));
budgets.DOC1(5).areal=budgets.DOC1(5).value./total_estuary_area./ndays;
budgets.DOC1(5).type='Sink';
%remineralization of DOC
budgets.DOC1(6).name='Remineralization of DOC';
budgets.DOC1(6).units='kg C';
budgets.DOC1(6).value=sum(sum(mnldoc1(:,estuary.nodes)));
budgets.DOC1(6).areal=budgets.DOC1(6).value./total_estuary_area./ndays;
budgets.DOC1(6).type='Sink';
%photochemical remineralization
budgets.DOC1(7).name='Photochemical Remin.';
budgets.DOC1(7).units= 'kg C'
budgets.DOC1(7).value=sum(sum(dpd10(:,estuary.nodes)));
budgets.DOC1(7).areal=budgets.DOC1(7).value./total_estuary_area./ndays;
budgets.DOC1(7).type='Sink';
%riverine input
budgets.DOC1(8).name='Riverine Inputs';
budgets.DOC1(8).units= 'kg C';
budgets.DOC1(8).ts=river.CDOC1.ts+river.NCDOC1.ts;% riverine inputs
budgets.DOC1(8).value=river.CDOC1.int+river.NCDOC1.int;% riverine inputs
budgets.DOC1(8).type='Source';
%photochem production
budgets.DOC1(9).name='Photochem Prod. DOC1';
budgets.DOC1(9).units= 'kg C';
budgets.DOC1(9).value=sum(sum(dpd21(:,estuary.nodes)...
                                        +  dpd31(:,estuary.nodes)...
                                        +  dpd3N(:,estuary.nodes).*FPDC31 ...
                                        +  dpd2N(:,estuary.nodes).*FPDC21...
                                        +  dpd1N(:,estuary.nodes).*FPDC11));
                                    
budgets.DOC1(9).areal=budgets.DOC1(9).value./total_estuary_area./ndays;                                   
budgets.DOC1(9).type='Source';
%photochemical loss of CDOC1
budgets.DOC1(10).name='Photochem Loss DOC1';
budgets.DOC1(10).units= 'kg C';
budgets.DOC1(10).value=sum(sum(dpd1N(:,estuary.nodes)));
budgets.DOC1(10).areal=budgets.DOC1(10).value./total_estuary_area./ndays;
budgets.DOC1(10).type='Sink';
%Coagulation of CDOC1
budgets.DOC1(11).name='Coagulation DOC1';
budgets.DOC1(11).units= 'kg C';
budgets.DOC1(11).value=0.0;
budgets.DOC1(11).areal=0.0;
budgets.DOC1(11).type='Sink';
%change in concentration of DOC % dt is negative if increases
budgets.DOC1(12).name='Delta DOC1';
budgets.DOC1(12).units= 'kg C';
budgets.DOC1(12).value=sum(sum(DTWCDOC1(:,estuary.nodes)+DTWNCDOC1(:,estuary.nodes)));
budgets.DOC1(12).areal=budgets.DOC1(12).value./total_estuary_area./ndays;
budgets.DOC1(12).type='Source';
 
myvals=[];
mynames={};

 %change the sign and extract the data for plotting
for i = 1 : length(budgets.DOC1);
    if(strcmp(budgets.DOC1(i).type,'Source'));
      myvals(i)=budgets.DOC1(i).value;
    else
      myvals(i)=budgets.DOC1(i).value*-1;
    end
    mynames{i}=budgets.DOC1(i).name;
end

% get the difference for the export to the mainstem (e.g. attributed to
% physics)
%Updated by B Clark March 2020 to use the real flux, output is positive out
%of the river, so need to change the sign
myvals(length(budgets.DOC1)+1)=sum(DOC1fx)./1e3.*-1;

mynames{length(budgets.DOC1)+1}='MainStem'
mypercent=[];

% Now change the sign of the terms based on source vs sink
sinks=myvals(myvals(1:end-2)<0);
sources=myvals(myvals(1:end-2)>=0);
mynames_out.DOC1=mynames;
for i = 1: length(myvals);
    mypercent(i)=myvals(i)./nansum(sources);
    mynames{i}=[mynames{i} ' ' num2str(round(mypercent(i)*100,3)) '%'];
end
myvals_out.DOC1=myvals;

% Below will make a bunch of plots
% make vector of ones to blow up the pie graph
explode1=ones(1,length(myvals(myvals>0)));
explode2=ones(1,length(myvals(myvals<0)));

figure;
pie(myvals(myvals>0),explode1,mynames(myvals>0));

saveas(gcf,[outdir '/DOC1_budget1'],'png');
saveas(gcf,[outdir '/DOC1_budget1'],'fig');
saveas(gcf,[outdir '/DOC1_budget1'],'eps');

figure;
pie(abs(myvals(myvals<0)),explode2,mynames(myvals<0));

saveas(gcf,[outdir '/DOC1_budget2'],'png');
saveas(gcf,[outdir '/DOC1_budget2'],'fig');
saveas(gcf,[outdir '/DOC1_budget2'],'eps');

figure;
b=bar(myvals,'FaceColor','flat')
xticklabels(mynames)
xtickangle(45)
%b.CData(myvals<0,:)=repmat([1 0 1],length(find(myvals<0)),1);
ylabel('DOC1 (kg C)');

saveas(gcf,[outdir '/DOC1_budget3'],'png');
saveas(gcf,[outdir '/DOC1_budget3'],'fig');
saveas(gcf,[outdir '/DOC1_budget3'],'eps');
%% DOC 2
%

% DOC2
%marsh sediment
budgets.DOC2(1).name='Marsh JDOC2';
budgets.DOC2(1).units='kg C';
budgets.DOC2(1).value=sum(JCDOC2(marsh.nodes)+JNCDOC2(marsh.nodes));
budgets.DOC2(1).areal=budgets.DOC2(1).value./total_marsh_area./ndays;%get the areal average per annum
budgets.DOC2(1).type='Source';
%estuary sediment
budgets.DOC2(2).name='Estuary JDOC2';
budgets.DOC2(2).units='kg C';
budgets.DOC2(2).value=sum(JCDOC2(estuary.nodesonly)+JNCDOC2(estuary.nodesonly));
budgets.DOC2(2).areal=budgets.DOC2(1).value./total_estuary_areaonly./ndays;%get the areal average per annum
budgets.DOC2(2).type='Source';
%algae production
budgets.DOC2(3).name='Estuary AlgDOC2';
budgets.DOC2(3).units='kg C';
budgets.DOC2(3).value=sum(sum(algdoc(:,estuary.nodes)))*FCD2;
budgets.DOC2(3).areal=budgets.DOC2(3).value./total_estuary_area./ndays;
budgets.DOC2(3).type='Source';
%hydrolysis of POC --> DOC
budgets.DOC2(4).name='Estuary HdrPOC2';
budgets.DOC2(4).units='kg C';
budgets.DOC2(4).value=sum(sum(hdrlpoc(:,estuary.nodes)+hdrrpoc(:,estuary.nodes)))*FHDRC2;
budgets.DOC2(4).areal=budgets.DOC2(4).value./total_estuary_area./ndays;
budgets.DOC2(4).type='Source';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%denitrification of DOC
budgets.DOC2(5).name='Denitrification of DOC';
budgets.DOC2(5).units='kg C';
budgets.DOC2(5).value=0.0;
budgets.DOC2(5).areal=0.0;
budgets.DOC2(5).type='Sink';
%remineralization of DOC
budgets.DOC2(6).name='Remineralization of DOC';
budgets.DOC2(6).units='kg C';
budgets.DOC2(6).value=sum(sum(mnldoc2(:,estuary.nodes)));
budgets.DOC2(6).areal=budgets.DOC2(6).value./total_estuary_area./ndays;
budgets.DOC2(6).type='Sink';
%photochemical remineralization
budgets.DOC2(7).name='Photochemical Remin.';
budgets.DOC2(7).units= 'kg C'
budgets.DOC2(7).value=sum(sum(dpd20(:,estuary.nodes)));
budgets.DOC2(7).areal=budgets.DOC2(7).value./total_estuary_area./ndays;
budgets.DOC2(7).type='Sink';
%riverine input
budgets.DOC2(8).name='Riverine Inputs';
budgets.DOC2(8).units= 'kg C';
budgets.DOC2(8).ts=river.CDOC2.ts+river.NCDOC2.ts;% riverine inputs
budgets.DOC2(8).value=river.CDOC2.int+river.NCDOC2.int;% riverine inputs
budgets.DOC2(8).type='Source';
%photochem production
budgets.DOC2(9).name='Photochem Prod. DOC2';
budgets.DOC2(9).units= 'kg C';
budgets.DOC2(9).value=sum(sum(dpd3N(:,estuary.nodes).*FPDC32 ...
                                        +  dpd2N(:,estuary.nodes).*FPDC22...
                                        +  dpd1N(:,estuary.nodes).*FPDC12));
                                    
budgets.DOC2(9).areal=budgets.DOC2(9).value./total_estuary_area./ndays;                                   
budgets.DOC2(9).type='Source';
%photochemical loss of CDOC2
budgets.DOC2(10).name='Photochem Loss DOC2';
budgets.DOC2(10).units= 'kg C';
budgets.DOC2(10).value=sum(sum(dpd2N(:,estuary.nodes)+dpd21(:,estuary.nodes)));
budgets.DOC2(10).areal=budgets.DOC2(10).value./total_estuary_area./ndays;
budgets.DOC2(10).type='Sink';
%Coagulation of CDOC2
budgets.DOC2(11).name='Coagulation DOC3';
budgets.DOC2(11).units= 'kg C';
budgets.DOC2(11).value=0.0;
budgets.DOC2(11).areal=0.0;
budgets.DOC2(11).type='Sink';

%change in concentration of DOC % dt is negative if increases
budgets.DOC2(12).name='Delta DOC2';
budgets.DOC2(12).units= 'kg C';
budgets.DOC2(12).value=sum(sum(DTWCDOC2(:,estuary.nodes)+DTWNCDOC2(:,estuary.nodes)));
budgets.DOC2(12).areal=budgets.DOC2(12).value./total_estuary_area./ndays;
budgets.DOC2(12).type='Source';
 
myvals=[];
mynames={};

 %change the sign and extract the data for plotting
for i = 1 : length(budgets.DOC2);
    if(strcmp(budgets.DOC2(i).type,'Source'));
      myvals(i)=budgets.DOC2(i).value;
    else
      myvals(i)=budgets.DOC2(i).value*-1;
    end
    mynames{i}=budgets.DOC2(i).name;
end

% get the difference for the export to the mainstem (e.g. attributed to
% physics)
myvals(length(budgets.DOC2)+1)=sum(DOC2fx)./1e3.*-1;;

mynames{length(budgets.DOC2)+1}='MainStem'
mypercent=[];

% Now change the sign of the terms based on source vs sink
sinks=myvals(myvals(1:end-2)<0);
sources=myvals(myvals(1:end-2)>=0);
mynames_out.DOC2=mynames;
for i = 1: length(myvals);
    mypercent(i)=myvals(i)./nansum(sources);
    mynames{i}=[mynames{i} ' ' num2str(round(mypercent(i)*100,3)) '%'];
end
myvals_out.DOC2=myvals;

% Below will make a bunch of plots
% make vector of ones to blow up the pie graph
explode1=ones(1,length(myvals(myvals>0)));
explode2=ones(1,length(myvals(myvals<0)));

figure;
pie(myvals(myvals>0),explode1,mynames(myvals>0));

saveas(gcf,[outdir '/DOC2_budget1'],'png');
saveas(gcf,[outdir '/DOC2_budget1'],'fig');
saveas(gcf,[outdir '/DOC2_budget1'],'eps');

figure;
pie(abs(myvals(myvals<0)),explode2,mynames(myvals<0));

saveas(gcf,[outdir '/DOC2_budget2'],'png');
saveas(gcf,[outdir '/DOC2_budget2'],'fig');
saveas(gcf,[outdir '/DOC2_budget2'],'eps');


figure;
b=bar(myvals,'FaceColor','flat')
xticklabels(mynames)
xtickangle(45)
%b.CData(myvals<0,:)=repmat([1 0 1],length(find(myvals<0)),1);
ylabel('DOC2 (kg C)');

saveas(gcf,[outdir '/DOC2_budget3'],'png');
saveas(gcf,[outdir '/DOC2_budget3'],'fig');
saveas(gcf,[outdir '/DOC2_budget3'],'eps');
%% DOC3 
%

% DOC3
%marsh sediment
budgets.DOC3(1).name='Marsh JDOC3';
budgets.DOC3(1).units='kg C';
budgets.DOC3(1).value=sum(JCDOC3(marsh.nodes)+JNCDOC3(marsh.nodes));
budgets.DOC3(1).areal=budgets.DOC3(1).value./total_marsh_area./ndays;%get the areal average per annum
budgets.DOC3(1).type='Source';
%estuary sediment
budgets.DOC3(2).name='Estuary JDOC3';
budgets.DOC3(2).units='kg C';
budgets.DOC3(2).value=sum(JCDOC3(estuary.nodesonly)+JNCDOC3(estuary.nodesonly));
budgets.DOC3(2).areal=budgets.DOC3(1).value./total_estuary_areaonly./ndays;%get the areal average per annum
budgets.DOC3(2).type='Source';
%algae production
budgets.DOC3(3).name='Estuary AlgDOC3';
budgets.DOC3(3).units='kg C';
budgets.DOC3(3).value=sum(sum(algdoc(:,estuary.nodes)))*FCD3;
budgets.DOC3(3).areal=budgets.DOC3(3).value./total_estuary_area./ndays;
budgets.DOC3(3).type='Source';
%hydrolysis of POC --> DOC
budgets.DOC3(4).name='Estuary HdrPOC3';
budgets.DOC3(4).units='kg C';
budgets.DOC3(4).value=sum(sum(hdrlpoc(:,estuary.nodes)+hdrrpoc(:,estuary.nodes)))*FHDRC3;
budgets.DOC3(4).areal=budgets.DOC3(4).value./total_estuary_area./ndays;
budgets.DOC3(4).type='Source';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%denitrification of DOC
budgets.DOC3(5).name='Denitrification of DOC';
budgets.DOC3(5).units='kg C';
budgets.DOC3(5).value=0.0;
budgets.DOC3(5).areal=0.0;
budgets.DOC3(5).type='Sink';
%remineralization of DOC
budgets.DOC3(6).name='Remineralization of DOC';
budgets.DOC3(6).units='kg C';
budgets.DOC3(6).value=sum(sum(mnldoc3(:,estuary.nodes)));
budgets.DOC3(6).areal=budgets.DOC3(6).value./total_estuary_area./ndays;
budgets.DOC3(6).type='Sink';
%photochemical remineralization
budgets.DOC3(7).name='Photochemical Remin.';
budgets.DOC3(7).units= 'kg C'
budgets.DOC3(7).value=sum(sum(dpd30(:,estuary.nodes)));
budgets.DOC3(7).areal=budgets.DOC3(7).value./total_estuary_area./ndays;
budgets.DOC3(7).type='Sink';
%riverine input
budgets.DOC3(8).name='Riverine Inputs';
budgets.DOC3(8).units= 'kg C';
budgets.DOC3(8).ts=river.CDOC3.ts+river.NCDOC3.ts;% riverine inputs
budgets.DOC3(8).value=river.CDOC3.int+river.NCDOC3.int;% riverine inputs
budgets.DOC3(8).type='Source';
%photochem production
budgets.DOC3(9).name='Photochem Prod. DOC3';
budgets.DOC3(9).units= 'kg C';
budgets.DOC3(9).value=sum(sum(dpd3N(:,estuary.nodes).*FPDC33 ...
                                        +  dpd2N(:,estuary.nodes).*FPDC23...
                                        +  dpd1N(:,estuary.nodes).*FPDC13));
                                    
budgets.DOC3(9).areal=budgets.DOC3(9).value./total_estuary_area./ndays;                                   
budgets.DOC3(9).type='Source';
%photochemical loss of CDOC3
budgets.DOC3(10).name='Photochem Loss DOC3';
budgets.DOC3(10).units= 'kg C';
budgets.DOC3(10).value=sum(sum(dpd3N(:,estuary.nodes)+dpd31(:,estuary.nodes)));
budgets.DOC3(10).areal=budgets.DOC3(10).value./total_estuary_area./ndays;
budgets.DOC3(10).type='Sink';
%Coagulation of CDOC3
budgets.DOC3(11).name='Coagulation DOC3';
budgets.DOC3(11).units= 'kg C';
budgets.DOC3(11).value=sum(sum(coagc(:,estuary.nodes)));
budgets.DOC3(11).areal=budgets.DOC3(11).value./total_estuary_area./ndays;
budgets.DOC3(11).type='Sink';
%change in concentration of DOC % dt is negative if increases
budgets.DOC3(12).name='Delta DOC3';
budgets.DOC3(12).units= 'kg C';
budgets.DOC3(12).value=sum(sum(DTWCDOC3(:,estuary.nodes)+DTWNCDOC3(:,estuary.nodes)));
budgets.DOC3(12).areal=budgets.DOC3(12).value./total_estuary_area./ndays;
budgets.DOC3(12).type='Source';
 
myvals=[];
mynames={};

 %change the sign and extract the data for plotting
for i = 1 : length(budgets.DOC3);
    if(strcmp(budgets.DOC3(i).type,'Source'));
      myvals(i)=budgets.DOC3(i).value;
    else
      myvals(i)=budgets.DOC3(i).value*-1;
    end
    mynames{i}=budgets.DOC3(i).name;
end

% get the difference for the export to the mainstem (e.g. attributed to
% physics)
myvals(length(budgets.DOC3)+1)=sum(DOC3fx)./1e3.*-1;;

mynames{length(budgets.DOC3)+1}='MainStem'
mypercent=[];

% Now change the sign of the terms based on source vs sink
sinks=myvals(myvals(1:end-2)<0);
sources=myvals(myvals(1:end-2)>=0);
mynames_out.DOC3=mynames;
for i = 1: length(myvals);
    mypercent(i)=myvals(i)./nansum(sources);
    mynames{i}=[mynames{i} ' ' num2str(round(mypercent(i)*100,3)) '%'];
end
myvals_out.DOC3=myvals;

% Below will make a bunch of plots
% make vector of ones to blow up the pie graph
explode1=ones(1,length(myvals(myvals>0)));
explode2=ones(1,length(myvals(myvals<0)));

figure;
pie(myvals(myvals>0),explode1,mynames(myvals>0));

saveas(gcf,[outdir '/DOC3_budget1'],'png');
saveas(gcf,[outdir '/DOC3_budget1'],'fig');
saveas(gcf,[outdir '/DOC3_budget1'],'eps');

figure;
pie(abs(myvals(myvals<0)),explode2,mynames(myvals<0));

saveas(gcf,[outdir '/DOC3_budget2'],'png');
saveas(gcf,[outdir '/DOC3_budget2'],'fig');
saveas(gcf,[outdir '/DOC3_budget2'],'eps');

figure;
b=bar(myvals,'FaceColor','flat')
xticklabels(mynames)
xtickangle(45)
%b.CData(myvals<0,:)=repmat([1 0 1],length(find(myvals<0)),1);
ylabel('DOC3 (kg C)');

saveas(gcf,[outdir '/DOC3_budget3'],'png');
saveas(gcf,[outdir '/DOC3_budget3'],'fig');
saveas(gcf,[outdir '/DOC3_budget3'],'eps');



    %% total DOC
    
    % Now change the sign of the terms based on source vs sink
    myvals=(myvals_out.DOC1(1:8)+myvals_out.DOC2(1:8)+myvals_out.DOC3(1:8));
    myvals=[myvals 0.0 0.0 myvals_out.DOC1(11:13)+myvals_out.DOC2(11:13)+myvals_out.DOC3(11:13)];
    
    sinks=myvals(myvals(1:end-2)<0);
    sources=myvals(myvals(1:end-2)>=0);
    mynames_out.tDOC=mynames;
    mypercent=[];
    mynames={};
    for i = 1: length(myvals);
        mypercent(i)=myvals(i)./nansum(sources);
        mynames{i}=[num2str(round(mypercent(i)*100,3)) '%'];
    end
    myvals_out.tDOC=myvals;
    

        figure;
        bar((myvals_out.DOC1+myvals_out.DOC2+myvals_out.DOC3)./1E3);
        hold on
        bar(([myvals_out.DOC1' myvals_out.DOC2' myvals_out.DOC3'])./1E3);
        legend('tDOC','DOC_1','DOC_2','DOC_3');
        xticklabels(mynames);
        xtickangle(45);
        %b.CData(myvals<0,:)=repmat([1 0 1],length(find(myvals<0)),1);
        ylabel('DOC (tons C)');
        
        saveas(gcf,[outdir '/tDOC_barplot'],'png');
        saveas(gcf,[outdir '/tDOC_barplot'],'fig');
        saveas(gcf,[outdir '/tDOC_barplot'],'eps');

%% POC

% POC
%marsh sediment
budgets.POC(1).name='Marsh JPOC';
budgets.POC(1).units='kg C';
budgets.POC(1).value=sum(JPOC1(marsh.nodes)+JPOC2(marsh.nodes)+JPOC3(marsh.nodes));
budgets.POC(1).areal=budgets.POC(1).value./total_marsh_area./ndays;%get the areal average per annum
budgets.POC(1).type='Sink';
%estuary sediment
budgets.POC(2).name='Estuary JPOC';
budgets.POC(2).units='kg C';
budgets.POC(2).value=sum(JPOC1(estuary.nodesonly)+JPOC2(estuary.nodesonly)+JPOC3(estuary.nodesonly));
budgets.POC(2).areal=budgets.POC(1).value./total_estuary_areaonly./ndays;%get the areal average per annum
budgets.POC(2).type='Sink';
%algae production
budgets.POC(3).name='Estuary AlgPOC';
budgets.POC(3).units='kg C';
budgets.POC(3).value=sum(sum(algpoc(:,estuary.nodes)));
budgets.POC(3).areal=budgets.POC(3).value./total_estuary_area./ndays;
budgets.POC(3).type='Source';
%hydrolysis of POC --> DOC
budgets.POC(4).name='Estuary HdrPOC2';
budgets.POC(4).units='kg C';
budgets.POC(4).value=sum(sum(hdrlpoc(:,estuary.nodes)+hdrrpoc(:,estuary.nodes)))*FHDRC2;
budgets.POC(4).areal=budgets.POC(4).value./total_estuary_area./ndays;
budgets.POC(4).type='Sink';

%riverine input
budgets.POC(5).name='Riverine Inputs';
budgets.POC(5).units= 'kg C';
budgets.POC(5).ts=river.LPOC.ts+river.RPOC.ts;% riverine inputs
budgets.POC(5).value=river.LPOC.int+river.RPOC.int;% riverine inputs
budgets.POC(5).type='Source';
%Coagulation of POC
budgets.POC(6).name='Coagulation POC';
budgets.POC(6).units= 'kg C';
budgets.POC(6).value=sum(sum(coagc(:,estuary.nodes)));
budgets.POC(6).areal=budgets.POC(6).value./total_estuary_area./ndays;
budgets.POC(6).type='Source';

%change in concentration of DON % dt is negative if increases
budgets.POC(7).name='Delta POC';
budgets.POC(7).units= 'kg C';
budgets.POC(7).value=sum(sum(DTLPOC(:,estuary.nodes)+DTRPOC(:,estuary.nodes)));
budgets.POC(7).areal=budgets.POC(7).value./total_estuary_area./ndays;
budgets.POC(7).type='Source';
 
myvals=[];
mynames={};

 %change the sign and extract the data for plotting
for i = 1 : length(budgets.POC);
    if(strcmp(budgets.POC(i).type,'Source'));
      myvals(i)=budgets.POC(i).value;
    else
      myvals(i)=budgets.POC(i).value*-1;
    end
    mynames{i}=budgets.POC(i).name;
end

% get the difference for the export to the mainstem (e.g. attributed to
% physics)
myvals(length(budgets.POC)+1)=sum(POCfx)./1e3.*-1;

mynames{length(budgets.POC)+1}='MainStem'
mypercent=[];

% Now change the sign of the terms based on source vs sink
sinks=myvals(myvals(1:end-2)<0);
sources=myvals(myvals(1:end-2)>=0);
mynames_out.POC=mynames;
for i = 1: length(myvals);
    mypercent(i)=myvals(i)./nansum(sources);
    mynames{i}=[mynames{i} ' ' num2str(round(mypercent(i)*100,3)) '%'];
end
myvals_out.POC=myvals;

% Below will make a bunch of plots
% make vector of ones to blow up the pie graph
explode1=ones(1,length(myvals(myvals>0)));
explode2=ones(1,length(myvals(myvals<0)));

figure;
pie(myvals(myvals>0),explode1,mynames(myvals>0));

saveas(gcf,[outdir '/POC_budget1'],'png');
saveas(gcf,[outdir '/POC_budget1'],'fig');
saveas(gcf,[outdir '/POC_budget1'],'eps');

figure;
pie(abs(myvals(myvals<0)),explode2,mynames(myvals<0));

saveas(gcf,[outdir '/POC_budget2'],'png');
saveas(gcf,[outdir '/POC_budget2'],'fig');
saveas(gcf,[outdir '/POC_budget2'],'eps');


figure;
b=bar(myvals,'FaceColor','flat')
xticklabels(mynames)
xtickangle(45)
%b.CData(myvals<0,:)=repmat([1 0 1],length(find(myvals<0)),1);
ylabel('POC (kg C)');

saveas(gcf,[outdir '/POC_budget3'],'png');
saveas(gcf,[outdir '/POC_budget3'],'fig');
saveas(gcf,[outdir '/POC_budget3'],'eps');

%% DON

disp('Calculating the Carbon budgets');
% DON1
%marsh sediment
budgets.DON1(1).name='Marsh JDON1';
budgets.DON1(1).units='kg N';
budgets.DON1(1).value=sum(JCDON1(marsh.nodes)+JNCDON1(marsh.nodes));
budgets.DON1(1).areal=budgets.DON1(1).value./total_marsh_area./ndays;%get the areal average per annum
budgets.DON1(1).type='Source';
%estuary sediment
budgets.DON1(2).name='Estuary JDON1';
budgets.DON1(2).units='kg N';
budgets.DON1(2).value=sum(JCDON1(estuary.nodesonly)+JNCDON1(estuary.nodesonly));
budgets.DON1(2).areal=budgets.DON1(1).value./total_estuary_areaonly./ndays;%get the areal average per annum
budgets.DON1(2).type='Source';
%algae production
budgets.DON1(3).name='Estuary AlgDON1';
budgets.DON1(3).units='kg N';
budgets.DON1(3).value=sum(sum(algdon(:,estuary.nodes)))*FND1;
budgets.DON1(3).areal=budgets.DON1(3).value./total_estuary_area./ndays;
budgets.DON1(3).type='Source';
%hydrolysis of PON --> DON
budgets.DON1(4).name='Estuary HdrPON1';
budgets.DON1(4).units='kg N';
budgets.DON1(4).value=sum(sum(hdrlpon(:,estuary.nodes)+hdrrpon(:,estuary.nodes)))*FHDRC1;
budgets.DON1(4).areal=budgets.DON1(4).value./total_estuary_area./ndays;
budgets.DON1(4).type='Source';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%denitrification of DON
budgets.DON1(5).name='Denitrification of DON';
budgets.DON1(5).units='kg N';
budgets.DON1(5).value=sum(sum(denitn(:,estuary.nodes)));
budgets.DON1(5).areal=budgets.DON1(5).value./total_estuary_area./ndays;
budgets.DON1(5).type='Sink';
%remineralization of DON
budgets.DON1(6).name='Remineralization of DON';
budgets.DON1(6).units='kg N';
budgets.DON1(6).value=sum(sum(mnldon1(:,estuary.nodes)));
budgets.DON1(6).areal=budgets.DON1(6).value./total_estuary_area./ndays;
budgets.DON1(6).type='Sink';
%photochemical remineralization
budgets.DON1(7).name='Photochemical Remin.';
budgets.DON1(7).units= 'kg N'
budgets.DON1(7).value=sum(sum(dpd10N(:,estuary.nodes)));
budgets.DON1(7).areal=budgets.DON1(7).value./total_estuary_area./ndays;
budgets.DON1(7).type='Sink';
%riverine input
budgets.DON1(8).name='Riverine Inputs';
budgets.DON1(8).units= 'kg N';
budgets.DON1(8).ts=river.CDON1.ts+river.NCDON1.ts;% riverine inputs
budgets.DON1(8).value=river.CDON1.int+river.NCDON1.int;% riverine inputs
budgets.DON1(8).type='Source';
%photochem production
budgets.DON1(9).name='Photochem Prod. DON1';
budgets.DON1(9).units= 'kg N';
budgets.DON1(9).value=sum(sum(dpd21N(:,estuary.nodes)...
                                        +  dpd31N(:,estuary.nodes)...
                                        +  dpd3NN(:,estuary.nodes).*FPDC31 ...
                                        +  dpd2NN(:,estuary.nodes).*FPDC21...
                                        +  dpd1NN(:,estuary.nodes).*FPDC11));
                                    
budgets.DON1(9).areal=budgets.DON1(9).value./total_estuary_area./ndays;                                   
budgets.DON1(9).type='Source';
%photochemical loss of CDON1
budgets.DON1(10).name='Photochem Loss DON1';
budgets.DON1(10).units= 'kg N';
budgets.DON1(10).value=sum(sum(dpd1NN(:,estuary.nodes)));
budgets.DON1(10).areal=budgets.DON1(10).value./total_estuary_area./ndays;
budgets.DON1(10).type='Sink';
%Coagulation of CDON1
budgets.DON1(11).name='Coagulation DON1';
budgets.DON1(11).units= 'kg N';
budgets.DON1(11).value=0.0;
budgets.DON1(11).areal=0.0;
budgets.DON1(11).type='Sink';
%change in concentration of DON % dt is negative if increases
budgets.DON1(12).name='Delta DON1';
budgets.DON1(12).units= 'kg N';
budgets.DON1(12).value=sum(sum(DTWCDON1(:,estuary.nodes)+DTWNCDON1(:,estuary.nodes)));
budgets.DON1(12).areal=budgets.DON1(12).value./total_estuary_area./ndays;
budgets.DON1(12).type='Source';
 
myvals=[];
mynames={};

 %change the sign and extract the data for plotting
for i = 1 : length(budgets.DON1);
    if(strcmp(budgets.DON1(i).type,'Source'));
      myvals(i)=budgets.DON1(i).value;
    else
      myvals(i)=budgets.DON1(i).value*-1;
    end
    mynames{i}=budgets.DON1(i).name;
end

% get the difference for the export to the mainstem (e.g. attributed to
% physics)
myvals(length(budgets.DON1)+1)=sum(DON1fx)./1e3.*-1;

mynames{length(budgets.DON1)+1}='MainStem'
mypercent=[];

% Now change the sign of the terms based on source vs sink
sinks=myvals(myvals(1:end-2)<0);
sources=myvals(myvals(1:end-2)>=0);
mynames_out.DON1=mynames;
for i = 1: length(myvals);
    mypercent(i)=myvals(i)./nansum(sources);
    mynames{i}=[mynames{i} ' ' num2str(round(mypercent(i)*100,3)) '%'];
end
myvals_out.DON1=myvals;

% Below will make a bunch of plots
% make vector of ones to blow up the pie graph
explode1=ones(1,length(myvals(myvals>0)));
explode2=ones(1,length(myvals(myvals<0)));

figure;
pie(myvals(myvals>0),explode1,mynames(myvals>0));

saveas(gcf,[outdir '/DON1_budget1'],'png');
saveas(gcf,[outdir '/DON1_budget1'],'fig');
saveas(gcf,[outdir '/DON1_budget1'],'eps');

figure;
pie(abs(myvals(myvals<0)),explode2,mynames(myvals<0));

saveas(gcf,[outdir '/DON1_budget2'],'png');
saveas(gcf,[outdir '/DON1_budget2'],'fig');
saveas(gcf,[outdir '/DON1_budget2'],'eps');


figure;
b=bar(myvals,'FaceColor','flat')
xticklabels(mynames)
xtickangle(45)
%b.CData(myvals<0,:)=repmat([1 0 1],length(find(myvals<0)),1);
ylabel('DON1 (kg N)');

saveas(gcf,[outdir '/DON1_budget3'],'png');
saveas(gcf,[outdir '/DON1_budget3'],'fig');
saveas(gcf,[outdir '/DON1_budget3'],'eps');
%% DON 2
%

% DON2
%marsh sediment
budgets.DON2(1).name='Marsh JDON2';
budgets.DON2(1).units='kg N';
budgets.DON2(1).value=sum(JCDON2(marsh.nodes)+JNCDON2(marsh.nodes));
budgets.DON2(1).areal=budgets.DON2(1).value./total_marsh_area./ndays;%get the areal average per annum
budgets.DON2(1).type='Source';
%estuary sediment
budgets.DON2(2).name='Estuary JDON2';
budgets.DON2(2).units='kg N';
budgets.DON2(2).value=sum(JCDON2(estuary.nodesonly)+JNCDON2(estuary.nodesonly));
budgets.DON2(2).areal=budgets.DON2(1).value./total_estuary_areaonly./ndays;%get the areal average per annum
budgets.DON2(2).type='Source';
%algae production
budgets.DON2(3).name='Estuary AlgDON2';
budgets.DON2(3).units='kg N';
budgets.DON2(3).value=sum(sum(algdon(:,estuary.nodes)))*FND2;
budgets.DON2(3).areal=budgets.DON2(3).value./total_estuary_area./ndays;
budgets.DON2(3).type='Source';
%hydrolysis of PON --> DON
budgets.DON2(4).name='Estuary HdrPON2';
budgets.DON2(4).units='kg N';
budgets.DON2(4).value=sum(sum(hdrlpon(:,estuary.nodes)+hdrrpon(:,estuary.nodes)))*FHDRC2;
budgets.DON2(4).areal=budgets.DON2(4).value./total_estuary_area./ndays;
budgets.DON2(4).type='Source';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%denitrification of DON
budgets.DON2(5).name='Denitrification of DON';
budgets.DON2(5).units='kg N';
budgets.DON2(5).value=0.0;
budgets.DON2(5).areal=0.0;
budgets.DON2(5).type='Sink';
%remineralization of DON
budgets.DON2(6).name='Remineralization of DON';
budgets.DON2(6).units='kg N';
budgets.DON2(6).value=sum(sum(mnldon2(:,estuary.nodes)));
budgets.DON2(6).areal=budgets.DON2(6).value./total_estuary_area./ndays;
budgets.DON2(6).type='Sink';
%photochemical remineralization
budgets.DON2(7).name='Photochemical Remin.';
budgets.DON2(7).units= 'kg N'
budgets.DON2(7).value=sum(sum(dpd20N(:,estuary.nodes)));
budgets.DON2(7).areal=budgets.DON2(7).value./total_estuary_area./ndays;
budgets.DON2(7).type='Sink';
%riverine input
budgets.DON2(8).name='Riverine Inputs';
budgets.DON2(8).units= 'kg N';
budgets.DON2(8).ts=river.CDON2.ts+river.NCDON2.ts;% riverine inputs
budgets.DON2(8).value=river.CDON2.int+river.NCDON2.int;% riverine inputs
budgets.DON2(8).type='Source';
%photochem production
budgets.DON2(9).name='Photochem Prod. DON2';
budgets.DON2(9).units= 'kg N';
budgets.DON2(9).value=sum(sum(dpd3NN(:,estuary.nodes).*FPDC32 ...
                                        +  dpd2NN(:,estuary.nodes).*FPDC22...
                                        +  dpd1NN(:,estuary.nodes).*FPDC12));
                                    
budgets.DON2(9).areal=budgets.DON2(9).value./total_estuary_area./ndays;                                   
budgets.DON2(9).type='Source';
%photochemical loss of CDON2
budgets.DON2(10).name='Photochem Loss DON2';
budgets.DON2(10).units= 'kg N';
budgets.DON2(10).value=sum(sum(dpd2NN(:,estuary.nodes)+dpd21N(:,estuary.nodes)));
budgets.DON2(10).areal=budgets.DON2(10).value./total_estuary_area./ndays;
budgets.DON2(10).type='Sink';
%Coagulation of CDON2
budgets.DON2(11).name='Coagulation DON3';
budgets.DON2(11).units= 'kg N';
budgets.DON2(11).value=0.0;
budgets.DON2(11).areal=0.0;
budgets.DON2(11).type='Sink';

%change in concentration of DON % dt is negative if increases
budgets.DON2(12).name='Delta DON2';
budgets.DON2(12).units= 'kg N';
budgets.DON2(12).value=sum(sum(DTWCDON2(:,estuary.nodes)+DTWNCDON2(:,estuary.nodes)));
budgets.DON2(12).areal=budgets.DON2(12).value./total_estuary_area./ndays;
budgets.DON2(12).type='Source';
 
myvals=[];
mynames={};

 %change the sign and extract the data for plotting
for i = 1 : length(budgets.DON2);
    if(strcmp(budgets.DON2(i).type,'Source'));
      myvals(i)=budgets.DON2(i).value;
    else
      myvals(i)=budgets.DON2(i).value*-1;
    end
    mynames{i}=budgets.DON2(i).name;
end

% get the difference for the export to the mainstem (e.g. attributed to
% physics)
myvals(length(budgets.DON2)+1)=sum(DON2fx)./1e3.*-1;

mynames{length(budgets.DON2)+1}='MainStem'
mypercent=[];

% Now change the sign of the terms based on source vs sink
sinks=myvals(myvals(1:end-2)<0);
sources=myvals(myvals(1:end-2)>=0);
mynames_out.DON2=mynames;
for i = 1: length(myvals);
    mypercent(i)=myvals(i)./nansum(sources);
    mynames{i}=[mynames{i} ' ' num2str(round(mypercent(i)*100,3)) '%'];
end
myvals_out.DON2=myvals;

% Below will make a bunch of plots
% make vector of ones to blow up the pie graph
explode1=ones(1,length(myvals(myvals>0)));
explode2=ones(1,length(myvals(myvals<0)));

figure;
pie(myvals(myvals>0),explode1,mynames(myvals>0));

saveas(gcf,[outdir '/DON2_budget1'],'png');
saveas(gcf,[outdir '/DON2_budget1'],'fig');
saveas(gcf,[outdir '/DON2_budget1'],'eps');

figure;
pie(abs(myvals(myvals<0)),explode2,mynames(myvals<0));

saveas(gcf,[outdir '/DON2_budget2'],'png');
saveas(gcf,[outdir '/DON2_budget2'],'fig');
saveas(gcf,[outdir '/DON2_budget2'],'eps');


figure;
b=bar(myvals,'FaceColor','flat')
xticklabels(mynames)
xtickangle(45)
%b.CData(myvals<0,:)=repmat([1 0 1],length(find(myvals<0)),1);
ylabel('DON2 (kg N)');

saveas(gcf,[outdir '/DON2_budget3'],'png');
saveas(gcf,[outdir '/DON2_budget3'],'fig');
saveas(gcf,[outdir '/DON2_budget3'],'eps');
%% DON3 
%

% DON3
%marsh sediment
budgets.DON3(1).name='Marsh JDON3';
budgets.DON3(1).units='kg N';
budgets.DON3(1).value=sum(JCDON3(marsh.nodes)+JNCDON3(marsh.nodes));
budgets.DON3(1).areal=budgets.DON3(1).value./total_marsh_area./ndays;%get the areal average per annum
budgets.DON3(1).type='Source';
%estuary sediment
budgets.DON3(2).name='Estuary JDON3';
budgets.DON3(2).units='kg N';
budgets.DON3(2).value=sum(JCDON3(estuary.nodesonly)+JNCDON3(estuary.nodesonly));
budgets.DON3(2).areal=budgets.DON3(1).value./total_estuary_areaonly./ndays;%get the areal average per annum
budgets.DON3(2).type='Source';
%algae production
budgets.DON3(3).name='Estuary AlgDON3';
budgets.DON3(3).units='kg N';
budgets.DON3(3).value=sum(sum(algdon(:,estuary.nodes)))*FND3;
budgets.DON3(3).areal=budgets.DON3(3).value./total_estuary_area./ndays;
budgets.DON3(3).type='Source';
%hydrolysis of POC --> DON
budgets.DON3(4).name='Estuary HdrPON3';
budgets.DON3(4).units='kg N';
budgets.DON3(4).value=sum(sum(hdrlpon(:,estuary.nodes)+hdrrpon(:,estuary.nodes)))*FHDRC3;
budgets.DON3(4).areal=budgets.DON3(4).value./total_estuary_area./ndays;
budgets.DON3(4).type='Source';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%denitrification of DON
budgets.DON3(5).name='Denitrification of DON';
budgets.DON3(5).units='kg N';
budgets.DON3(5).value=0.0;
budgets.DON3(5).areal=0.0;
budgets.DON3(5).type='Sink';
%remineralization of DON
budgets.DON3(6).name='Remineralization of DON';
budgets.DON3(6).units='kg N';
budgets.DON3(6).value=sum(sum(mnldon3(:,estuary.nodes)));
budgets.DON3(6).areal=budgets.DON3(6).value./total_estuary_area./ndays;
budgets.DON3(6).type='Sink';
%photochemical remineralization
budgets.DON3(7).name='Photochemical Remin.';
budgets.DON3(7).units= 'kg N'
budgets.DON3(7).value=sum(sum(dpd30N(:,estuary.nodes)));
budgets.DON3(7).areal=budgets.DON3(7).value./total_estuary_area./ndays;
budgets.DON3(7).type='Sink';
%riverine input
budgets.DON3(8).name='Riverine Inputs';
budgets.DON3(8).units= 'kg N';
budgets.DON3(8).ts=river.CDON3.ts+river.NCDON3.ts;% riverine inputs
budgets.DON3(8).value=river.CDON3.int+river.NCDON3.int;% riverine inputs
budgets.DON3(8).type='Source';
%photochem production
budgets.DON3(9).name='Photochem Prod. DON3';
budgets.DON3(9).units= 'kg N';
budgets.DON3(9).value=sum(sum(dpd3NN(:,estuary.nodes).*FPDC33 ...
                                        +  dpd2NN(:,estuary.nodes).*FPDC23...
                                        +  dpd1NN(:,estuary.nodes).*FPDC13));
                                    
budgets.DON3(9).areal=budgets.DON3(9).value./total_estuary_area./ndays;                                   
budgets.DON3(9).type='Source';
%photochemical loss of CDON3
budgets.DON3(10).name='Photochem Loss DON3';
budgets.DON3(10).units= 'kg N';
budgets.DON3(10).value=sum(sum(dpd3NN(:,estuary.nodes)+dpd31N(:,estuary.nodes)));
budgets.DON3(10).areal=budgets.DON3(10).value./total_estuary_area./ndays;
budgets.DON3(10).type='Sink';
%Coagulation of CDON3
budgets.DON3(11).name='Coagulation DON3';
budgets.DON3(11).units= 'kg N';
budgets.DON3(11).value=sum(sum(coagn(:,estuary.nodes)));
budgets.DON3(11).areal=budgets.DON3(11).value./total_estuary_area./ndays;
budgets.DON3(11).type='Sink';
%change in concentration of DON % dt is negative if increases
budgets.DON3(12).name='Delta DON3';
budgets.DON3(12).units= 'kg N';
budgets.DON3(12).value=sum(sum(DTWCDON3(:,estuary.nodes)+DTWNCDON3(:,estuary.nodes)));
budgets.DON3(12).areal=budgets.DON3(12).value./total_estuary_area./ndays;
budgets.DON3(12).type='Source';
 
myvals=[];
mynames={};

 %change the sign and extract the data for plotting
for i = 1 : length(budgets.DON3);
    if(strcmp(budgets.DON3(i).type,'Source'));
      myvals(i)=budgets.DON3(i).value;
    else
      myvals(i)=budgets.DON3(i).value*-1;
    end
    mynames{i}=budgets.DON3(i).name;
end

% get the difference for the export to the mainstem (e.g. attributed to
% physics)
myvals(length(budgets.DON3)+1)=sum(DON3fx)./1e3.*-1;

mynames{length(budgets.DON3)+1}='MainStem'
mypercent=[];

% Now change the sign of the terms based on source vs sink
sinks=myvals(myvals(1:end-2)<0);
sources=myvals(myvals(1:end-2)>=0);
mynames_out.DON3=mynames;
for i = 1: length(myvals);
    mypercent(i)=myvals(i)./nansum(sources);
    mynames{i}=[mynames{i} ' ' num2str(round(mypercent(i)*100,3)) '%'];
end
myvals_out.DON3=myvals;

% Below will make a bunch of plots
% make vector of ones to blow up the pie graph
explode1=ones(1,length(myvals(myvals>0)));
explode2=ones(1,length(myvals(myvals<0)));

figure;
pie(myvals(myvals>0),explode1,mynames(myvals>0));

saveas(gcf,[outdir '/DON3_budget1'],'png');
saveas(gcf,[outdir '/DON3_budget1'],'fig');
saveas(gcf,[outdir '/DON3_budget1'],'eps');

figure;
pie(abs(myvals(myvals<0)),explode2,mynames(myvals<0));

saveas(gcf,[outdir '/DON3_budget2'],'png');
saveas(gcf,[outdir '/DON3_budget2'],'fig');
saveas(gcf,[outdir '/DON3_budget2'],'eps');


figure;
b=bar(myvals,'FaceColor','flat')
xticklabels(mynames)
xtickangle(45)
%b.CData(myvals<0,:)=repmat([1 0 1],length(find(myvals<0)),1);
ylabel('DON3 (kg N)');

saveas(gcf,[outdir '/DON3_budget3'],'png');
saveas(gcf,[outdir '/DON3_budget3'],'fig');
saveas(gcf,[outdir '/DON3_budget3'],'eps');



%%total DON

% Now change the sign of the terms based on source vs sink
 myvals=(myvals_out.DON1+myvals_out.DON2+myvals_out.DON3);

sinks=myvals(myvals(1:end-2)<0);
sources=myvals(myvals(1:end-2)>=0);
mynames_out.tDON=mynames;

for i = 1: length(myvals);
    mypercent(i)=myvals(i)./nansum(sources);
    mynames{i}=[num2str(round(mypercent(i)*100,3)) '%'];
end
myvals_out.tDON=myvals;

figure;
bar((myvals_out.DON1+myvals_out.DON2+myvals_out.DON3)./1E3);
hold on
bar(([myvals_out.DON1' myvals_out.DON2' myvals_out.DON3'])./1E3);
legend('tDON','DON_1','DON_2','DON_3')
xticklabels(mynames);
xtickangle(45);
%b.CData(myvals<0,:)=repmat([1 0 1],length(find(myvals<0)),1);
ylabel('DON3 (kg N)');

saveas(gcf,[outdir '/tDON_barplot'],'png');
saveas(gcf,[outdir '/tDON_barplot'],'fig');
saveas(gcf,[outdir '/tDON_barplot'],'eps');

%% PON

% PON
%marsh sediment
budgets.PON(1).name='Marsh JPON';
budgets.PON(1).units='kg C';
budgets.PON(1).value=sum(JPON1(marsh.nodes)+JPON2(marsh.nodes)+JPON3(marsh.nodes));
budgets.PON(1).areal=budgets.PON(1).value./total_marsh_area./ndays;%get the areal average per annum
budgets.PON(1).type='Sink';
%estuary sediment
budgets.PON(2).name='Estuary JPON';
budgets.PON(2).units='kg C';
budgets.PON(2).value=sum(JPON1(estuary.nodesonly)+JPON2(estuary.nodesonly)+JPON3(estuary.nodesonly));
budgets.PON(2).areal=budgets.PON(1).value./total_estuary_areaonly./ndays;%get the areal average per annum
budgets.PON(2).type='Sink';
%algae production
budgets.PON(3).name='Estuary AlgPON';
budgets.PON(3).units='kg C';
budgets.PON(3).value=sum(sum(algpon(:,estuary.nodes)));
budgets.PON(3).areal=budgets.PON(3).value./total_estuary_area./ndays;
budgets.PON(3).type='Source';
%hydrolysis of PON --> DON
budgets.PON(4).name='Estuary HdrPON2';
budgets.PON(4).units='kg C';
budgets.PON(4).value=sum(sum(hdrlpon(:,estuary.nodes)+hdrrpon(:,estuary.nodes)))*FHDRC2;
budgets.PON(4).areal=budgets.PON(4).value./total_estuary_area./ndays;
budgets.PON(4).type='Sink';

%riverine input
budgets.PON(5).name='Riverine Inputs';
budgets.PON(5).units= 'kg C';
budgets.PON(5).ts=river.LPON.ts+river.RPON.ts;% riverine inputs
budgets.PON(5).value=river.LPON.int+river.RPON.int;% riverine inputs
budgets.PON(5).type='Source';
%Coagulation of PON
budgets.PON(6).name='Coagulation PON';
budgets.PON(6).units= 'kg C';
budgets.PON(6).value=sum(sum(coagn(:,estuary.nodes)));
budgets.PON(6).areal=budgets.PON(6).value./total_estuary_area./ndays;
budgets.PON(6).type='Source';

%change in concentration of DON % dt is negative if increases
budgets.PON(7).name='Delta PON';
budgets.PON(7).units= 'kg C';
budgets.PON(7).value=sum(sum(DTLPON(:,estuary.nodes)+DTRPON(:,estuary.nodes)));
budgets.PON(7).areal=budgets.PON(7).value./total_estuary_area./ndays;
budgets.PON(7).type='Source';
 
myvals=[];
mynames={};

 %change the sign and extract the data for plotting
for i = 1 : length(budgets.PON);
    if(strcmp(budgets.PON(i).type,'Source'));
      myvals(i)=budgets.PON(i).value;
    else
      myvals(i)=budgets.PON(i).value*-1;
    end
    mynames{i}=budgets.PON(i).name;
end

% get the difference for the export to the mainstem (e.g. attributed to
% physics)
myvals(length(budgets.PON)+1)=sum(PONfx)./1e3.*-1;

mynames{length(budgets.PON)+1}='MainStem';
mypercent=[];

% Now change the sign of the terms based on source vs sink
sinks=myvals(myvals(1:end-2)<0);
sources=myvals(myvals(1:end-2)>=0);
mynames_out.PON=mynames;
for i = 1: length(myvals);
    mypercent(i)=myvals(i)./nansum(sources);
    mynames{i}=[mynames{i} ' ' num2str(round(mypercent(i)*100,3)) '%'];
end
myvals_out.PON=myvals;

% Below will make a bunch of plots
% make vector of ones to blow up the pie graph
explode1=ones(1,length(myvals(myvals>0)));
explode2=ones(1,length(myvals(myvals<0)));

figure;
pie(myvals(myvals>0),explode1,mynames(myvals>0));

saveas(gcf,[outdir '/PON_budget1'],'png');
saveas(gcf,[outdir '/PON_budget1'],'fig');
saveas(gcf,[outdir '/PON_budget1'],'eps');

figure;
pie(abs(myvals(myvals<0)),explode2,mynames(myvals<0));

saveas(gcf,[outdir '/PON_budget2'],'png');
saveas(gcf,[outdir '/PON_budget2'],'fig');
saveas(gcf,[outdir '/PON_budget2'],'eps');

figure;
b=bar(myvals,'FaceColor','flat')
xticklabels(mynames)
xtickangle(45)
%b.CData(myvals<0,:)=repmat([1 0 1],length(find(myvals<0)),1);
ylabel('PON (kg C)');

saveas(gcf,[outdir '/PON_budget3'],'png');
saveas(gcf,[outdir '/PON_budget3'],'fig');
saveas(gcf,[outdir '/PON_budget3'],'eps');



%% NH4

% NH4
%marsh sediment
budgets.NH4(1).name='Marsh JNH4';
budgets.NH4(1).units='kg N';
budgets.NH4(1).value=sum(JNH4(marsh.nodes));
budgets.NH4(1).areal=budgets.NH4(1).value./total_marsh_area./ndays;%get the areal average per annum
budgets.NH4(1).type='Source';
%estuary sediment
budgets.NH4(2).name='Estuary JNH4';
budgets.NH4(2).units='kg N';
budgets.NH4(2).value=sum(JNH4(estuary.nodesonly));
budgets.NH4(2).areal=budgets.NH4(1).value./total_estuary_areaonly./ndays;%get the areal average per annum
budgets.NH4(2).type='Source';
%algae uptake
budgets.NH4(3).name='Estuary AlgNH4';
budgets.NH4(3).units='kg N';
budgets.NH4(3).value=sum(sum(algnh4(:,estuary.nodes)));
budgets.NH4(3).areal=budgets.NH4(3).value./total_estuary_area./ndays;
budgets.NH4(3).type='Sink';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%denitrification of DON
budgets.NH4(4).name='Denitrification of DON';
budgets.NH4(4).units='kg N';
budgets.NH4(4).value=sum(sum(denitn(:,estuary.nodes)));
budgets.NH4(4).areal=budgets.NH4(4).value./total_estuary_area./ndays;
budgets.NH4(4).type='Source';
%remineralization of DON
budgets.NH4(5).name='Remineralization of DON';
budgets.NH4(5).units='kg N';
budgets.NH4(5).value=sum(sum(mnldon1(:,estuary.nodes)+mnldon2(:,estuary.nodes)+mnldon3(:,estuary.nodes)));
budgets.NH4(5).areal=budgets.NH4(5).value./total_estuary_area./ndays;
budgets.NH4(5).type='Source';
%photochemical remineralization
budgets.NH4(6).name='Photochemical Remin.';
budgets.NH4(6).units= 'kg N'
budgets.NH4(6).value=sum(sum(dpd30N(:,estuary.nodes)));
budgets.NH4(6).areal=budgets.NH4(6).value./total_estuary_area./ndays;
budgets.NH4(6).type='Sink';
%riverine input
budgets.NH4(7).name='Riverine Inputs';
budgets.NH4(7).units= 'kg N';
budgets.NH4(7).ts=river.NH4.ts;% riverine inputs
budgets.NH4(7).value=river.NH4.int;% riverine inputs
budgets.NH4(7).type='Source';

%Nitrification loss of NH4
budgets.NH4(8).name='Nitrification of NH4';
budgets.NH4(8).units='kg N';
budgets.NH4(8).value=sum(sum(nt(:,estuary.nodes)));
budgets.NH4(8).areal=budgets.NH4(8).value./total_estuary_area./ndays;
budgets.NH4(8).type='Sink';
%change in concentration of DON % dt is negative if increases
budgets.NH4(9).name='Delta NH4';
budgets.NH4(9).units= 'kg N';
budgets.NH4(9).value=sum(sum(DTNH4(:,estuary.nodes)));
budgets.NH4(9).areal=budgets.NH4(9).value./total_estuary_area./ndays;
budgets.NH4(9).type='Source';
 
myvals=[];
mynames={};

 %change the sign and extract the data for plotting
for i = 1 : length(budgets.NH4);
    if(strcmp(budgets.NH4(i).type,'Source'));
      myvals(i)=budgets.NH4(i).value;
    else
      myvals(i)=budgets.NH4(i).value*-1;
    end
    mynames{i}=budgets.NH4(i).name;
end

% get the difference for the export to the mainstem (e.g. attributed to
% physics)
myvals(length(budgets.NH4)+1)=sum(NH4fx)./1e3.*-1;

mynames{length(budgets.NH4)+1}='MainStem'
mypercent=[];

% Now change the sign of the terms based on source vs sink
sinks=myvals(myvals(1:end-2)<0);
sources=myvals(myvals(1:end-2)>=0);
mynames_out.NH4=mynames;
for i = 1: length(myvals);
    mypercent(i)=myvals(i)./nansum(sources);
    mynames{i}=[mynames{i} ' ' num2str(round(mypercent(i)*100,3)) '%'];
end
myvals_out.NH4=myvals;

% Below will make a bunch of plots
% make vector of ones to blow up the pie graph
explode1=ones(1,length(myvals(myvals>0)));
explode2=ones(1,length(myvals(myvals<0)));

figure;
pie(myvals(myvals>0),explode1,mynames(myvals>0));

saveas(gcf,[outdir '/NH4_budget1'],'png');
saveas(gcf,[outdir '/NH4_budget1'],'fig');
saveas(gcf,[outdir '/NH4_budget1'],'eps');

figure;
pie(abs(myvals(myvals<0)),explode2,mynames(myvals<0));

saveas(gcf,[outdir '/NH4_budget2'],'png');
saveas(gcf,[outdir '/NH4_budget2'],'fig');
saveas(gcf,[outdir '/NH4_budget2'],'eps');

figure;
b=bar(myvals,'FaceColor','flat')
xticklabels(mynames)
xtickangle(45)
%b.CData(myvals<0,:)=repmat([1 0 1],length(find(myvals<0)),1);
ylabel('NH4 (kg N)');

saveas(gcf,[outdir '/NH4_budget3'],'png');
saveas(gcf,[outdir '/NH4_budget3'],'fig');
saveas(gcf,[outdir '/NH4_budget3'],'eps');



%% NO3

% NO3
%marsh sediment
budgets.NO3(1).name='Marsh JNO3';
budgets.NO3(1).units='kg N';
budgets.NO3(1).value=sum(JNO3(marsh.nodes));
budgets.NO3(1).areal=budgets.NO3(1).value./total_marsh_area./ndays;%get the areal average per annum
budgets.NO3(1).type='Source';
%estuary sediment
budgets.NO3(2).name='Estuary JNO3';
budgets.NO3(2).units='kg N';
budgets.NO3(2).value=sum(JNO3(estuary.nodesonly));
budgets.NO3(2).areal=budgets.NO3(1).value./total_estuary_areaonly./ndays;%get the areal average per annum
budgets.NO3(2).type='Source';
%algae uptake
budgets.NO3(3).name='Estuary AlgNO3';
budgets.NO3(3).units='kg N';
budgets.NO3(3).value=sum(sum(algno3(:,estuary.nodes)));
budgets.NO3(3).areal=budgets.NO3(3).value./total_estuary_area./ndays;
budgets.NO3(3).type='Sink';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%denitrification of NO3
budgets.NO3(4).name='Denitrification of NO3';
budgets.NO3(4).units='kg N';
budgets.NO3(4).value=sum(sum(denno3(:,estuary.nodes)));
budgets.NO3(4).areal=budgets.NO3(4).value./total_estuary_area./ndays;
budgets.NO3(4).type='Sink';

%riverine input
budgets.NO3(5).name='Riverine Inputs';
budgets.NO3(5).units= 'kg N';
budgets.NO3(5).ts=river.NO3.ts;% riverine inputs
budgets.NO3(5).value=river.NO3.int;% riverine inputs
budgets.NO3(5).type='Source';

%Nitrification Production of NO3
budgets.NO3(6).name='Nitrification of NO3';
budgets.NO3(6).units='kg N';
budgets.NO3(6).value=sum(sum(nt(:,estuary.nodes)));
budgets.NO3(6).areal=budgets.NO3(6).value./total_estuary_area./ndays;
budgets.NO3(6).type='Source';
%change in concentration of DON % dt is negative if increases
budgets.NO3(7).name='Delta NO3';
budgets.NO3(7).units= 'kg N';
budgets.NO3(7).value=sum(sum(DTNO3(:,estuary.nodes)));
budgets.NO3(7).areal=budgets.NO3(7).value./total_estuary_area./ndays;
budgets.NO3(7).type='Source';
 
myvals=[];
mynames={};

 %change the sign and extract the data for plotting
for i = 1 : length(budgets.NO3);
    if(strcmp(budgets.NO3(i).type,'Source'));
      myvals(i)=budgets.NO3(i).value;
    else
      myvals(i)=budgets.NO3(i).value*-1;
    end
    mynames{i}=budgets.NO3(i).name;
end

% get the difference for the export to the mainstem (e.g. attributed to
% physics)
myvals(length(budgets.NO3)+1)=sum(NO3fx)./1e3.*-1;

mynames{length(budgets.NO3)+1}='MainStem'
mypercent=[];

% Now change the sign of the terms based on source vs sink
sinks=myvals(myvals(1:end-2)<0);
sources=myvals(myvals(1:end-2)>=0);
mynames_out.NO3=mynames;
for i = 1: length(myvals);
    mypercent(i)=myvals(i)./nansum(sources);
    mynames{i}=[mynames{i} ' ' num2str(round(mypercent(i)*100,3)) '%'];
end
myvals_out.NO3=myvals;

% Below will make a bunch of plots
% make vector of ones to blow up the pie graph
explode1=ones(1,length(myvals(myvals>0)));
explode2=ones(1,length(myvals(myvals<0)));

figure;
pie(myvals(myvals>0),explode1,mynames(myvals>0));

saveas(gcf,[outdir '/NO3_budget1'],'png');
saveas(gcf,[outdir '/NO3_budget1'],'fig');
saveas(gcf,[outdir '/NO3_budget1'],'eps');

figure;
pie(abs(myvals(myvals<0)),explode2,mynames(myvals<0));

saveas(gcf,[outdir '/NO3_budget2'],'png');
saveas(gcf,[outdir '/NO3_budget2'],'fig');
saveas(gcf,[outdir '/NO3_budget2'],'eps');

figure;
b=bar(myvals,'FaceColor','flat')
xticklabels(mynames)
xtickangle(45)
%b.CData(myvals<0,:)=repmat([1 0 1],length(find(myvals<0)),1);
ylabel('NO3 (kg N)');

saveas(gcf,[outdir '/NO3_budget3'],'png');
saveas(gcf,[outdir '/NO3_budget3'],'fig');
saveas(gcf,[outdir '/NO3_budget3'],'eps');




% save

save(['budgets_poly' num2str(nodes_flag) '.mat'],'budgets','myvals_out');














