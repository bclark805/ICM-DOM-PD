%Nitrogen Budget calculation for ICM with netcdf outputs

% B Clark UMCES HPL Nov 2018
% Updated B Clark UMCES, Jan 2019
%
%ALl non-concentration variables in the budget (e.g. right had side of
%equation)
% are integrated internally in ICM and loaded from the hisdata*.mat files

%load in the files that define the nodes in the Tributary
%the river discharge file
%and the watershed inputs of the water quality variables
% load('*/RhodeFVCOM_2005/Data/RhodeRiver_nodes.mat');
load(input('What is the full path and file name to Tributary Node Information? --> ','s'));
% load('/data/users/bclark/CBAY_DATA/RR_discharge_2005.mat');
load(input('What is the full path and file name to the River discharge file? --> ','s'));
% load('*/RhodeFVCOM_2005/Data/RR_Watershed_WQ_2005.mat');
load(input('What is the full path and file name to the Watershed Water Quality Input data? --> ','s'));

%load the sediment definition file where marsh nodes are flagged with 1
% marsh_file= '/data/users/bclark/CBAY_DATA/seddom_grid.dat';
marsh_file=input('What is the full path and file name to the sediment definition file? --> ','s');
marsh_data=importdata(marsh_file);
% load the grid file
% grid_file='/data/users/bclark/chesfvm/cpb/Water_Quality/Rhode_river_PD/rhoderiver_grid_v1_2.mat';
grid_file=input('What is the full path and file name to the grid definition *.mat file? --> ','s');
load(grid_file);
%load the node area file
% [node_area]=importdata('/data/users/bclark/CBAY_DATA/node_CVS_area.dat');%get our nodal area for marsh flux calcs;
node_file=input('What is the full path and file name to the Node Area file? --->  ','s');
[node_area]=importdata(node_file);
%load the file with the different polygons from the Rhode River
%each polygon has a different set of nodes with a decreasing area of the
%Rhode River
% load('/data/users/bclark/CBAY_DATA/polygon_work.mat');
load(input('What is the full path and file name of the polygon definition file? --> ','s'));

first_day =input('What is the starting day of the year for this run (i.e. April 1 is 91 ---> ');

%this is the FVCOM hydrodynamics file, used to load in the vertical
%coordinate system
ncfile=input('What is the full path and file name to the first netcdf hydrodynamics file? --->  ','s');
nc=netcdf(ncfile);
s_hybrid=nc{'s_hybrid'}(:);
dz=diff(s_hybrid,1);
for i = 1:10;
    layer_thickness(i,:)=dz(i,:).*xyd_n(:,4)'*-1;
end
%
% Use this flag depending what region to set up the budget for, 1 being the
% whole tributary decreasing size of estuary area to 5
nodes_flag=input('What polygon do you want to build a budget for? 1-5 for trib, 0 for mainstem --->  ')
%get the rest of the information
wqm_hisdir = input('What is full path to the WQM history directory with the netcdf file? --->  ','s');
numfile=input('How many files do you want to create a budget for? ---> ')
% ICM_int = input('What is the interval that the data was output at, in days?  ---> ')
% ncfilein = input('What is the full path and file name of the FVCOM nc file used to run the WQM? ---> ','s');
% nc=netcdf(ncfilein);
nameseg=input('What is the name segment for the FVCOM files ---> ','s');

%SET UP EMPTY ARRAYS FOR ALL THE BUDGET VARIABLES

fname=[wqm_hisdir '/rhd_0001.nc']

myinfo=ncinfo(fname);

ngrids=myinfo.Dimensions(4).Length;
nlayers=myinfo.Dimensions(3).Length;
ntimes=myinfo.Dimensions(end).Length;

%%
%how many days was the model run for?

%
mkdir('budget_plots');
%
%load marsh definition file
marsh.nodes= find(marsh_data(:,3)>0);
marsh.area=node_area(marsh.nodes,2);

if (nodes_flag)
   estuary.nodes=1:size(lld_n);%RR_nodes_out(find(lld_n(RR_nodes_out,3) > 38.875));
else
    estuary.nodes=RR_nodes_out(find(lld_n(RR_nodes_out,3) > 38.875));
end
% 

% this is the entire area including the marsh
% estuary.nodes=poly4_nodes;

estuary.nodes=estuary.nodes;
estuary.area=node_area(estuary.nodes,2);
% get the nodes

% this is only the subtidal area, excluding the marsh
%use in the sediment flux calculations
    estuary.nodesonly=setdiff(estuary.nodes,marsh.nodes);
    estuary.areaonly=node_area(estuary.nodesonly,2);
    
tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3),lld_n(:,4));
hold on
plot(lld_n(estuary.nodes,2),lld_n(estuary.nodes,3),'kd')
plot(lld_n(marsh.nodes,2),lld_n(marsh.nodes,3),'ro')
saveas(gcf,'budget_plots/Marsh_map.png')
%
%    
 nlayers=10;
 
 zLayer=diff(s_hybrid,1)*-1;% the amount of water column for each layer, adds up to 1;
tic
    nc = netcdf([fname]);
    DEPTH=zeros(ntimes,ngrids);
    
    DEPTH=nc{'depth'}(:);
    time_vector_in=nc{'time'}(:);
    time_vector=(time_vector_in(2:end)./86400)+first_day;
    time_vector_hr=time_vector.*24;
    ndays=round(time_vector(end))-first_day;
    disp('Got time and depth data')
    toc

%%
tic
disp('getting river data')
%
%algal fractionation
FCD11=1.0;FCD12=0.1;FCD13=0.05
FCD21=1.0;FCD22=0.1;FCD23=0.05
FNCD11=0.3;FNCD12=0.45;FNCD13=0.0;
FNCD21=0.3;FNCD22=0.45;FNCD23=0.0;
FCND11=0.2;FCND12=0.05;FCND13=0.0;
FCND21=0.2;FCND22=0.05;FCND23=0.0;
FNCND11=0.6;FNCND12=0.15;FNCND13=0.0;
FNCND21=0.6;FNCND22=0.15;FNCND23=0.0;

%total algal fractions
FCD1=0.4;FCD2=0.55;FCD3=0.05;
FND1=0.8;FND2=0.2;FND3=0.0;

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

matching_times=dsearchn(river_time',time_vector_hr);
%
cut_discharge=mydischarge(matching_times,1:8);
% discharge_Myint=(cut_discharge(2:round(ICM_int*24):end,:));
discharge_Myint=cut_discharge;
river.times=Final_RR_forcing(1:45:end,2)+first_day*24;
%
Final_RR_forcing=Final_RR_forcing(:,1:8);

% will take the concentration and multiple my discharge g m^-3 x m^3 s^-1 *
% 3600 secs hr^-1 then integrate to get kg
river.CDON1.all=interp1(river.times,Final_RR_forcing(17:45:end,:),time_vector_hr,'pchip');
river.NCDON1.all=interp1(river.times,Final_RR_forcing(18:45:end,:),time_vector_hr,'pchip');
river.LPON.all=interp1(river.times,Final_RR_forcing(19:45:end,:),time_vector_hr,'pchip');
river.RPON.all=interp1(river.times,Final_RR_forcing(20:45:end,:),time_vector_hr,'pchip');
river.CDON2.all=interp1(river.times,Final_RR_forcing(38:45:end,:),time_vector_hr,'pchip');
river.NCDON2.all=interp1(river.times,Final_RR_forcing(39:45:end,:),time_vector_hr,'pchip');
river.CDON3.all=interp1(river.times,Final_RR_forcing(40:45:end,:),time_vector_hr,'pchip');
river.NCDON3.all=interp1(river.times,Final_RR_forcing(41:45:end,:),time_vector_hr,'pchip');

river.NH4.all=interp1(river.times,Final_RR_forcing(14:45:end,:),time_vector_hr,'pchip');
river.NO3.all=interp1(river.times,Final_RR_forcing(15:45:end,:),time_vector_hr,'pchip');

% now integrate and get one time series for each constituent

river.CDON1.ts=nansum(river.CDON1.all.*discharge_Myint*86400*ICM_int./1000,2);% kg per ICM_int period
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

 tic
JNH4=zeros(ntimes,ngrids);
JNO3=zeros(ntimes,ngrids)    ;
JPOC=zeros(ntimes,ngrids)    ;
JPON=zeros(ntimes,ngrids)    ;

JWCDON1=zeros(ntimes,ngrids)    ;
JWCDON2=zeros(ntimes,ngrids)    ;
JWCDON3=zeros(ntimes,ngrids)    ;
JWNCDON1=zeros(ntimes,ngrids)    ;
JWNCDON2=zeros(ntimes,ngrids)    ;
JWNCDON3=zeros(ntimes,ngrids)    ;

LPON=zeros(ntimes,nlayers,ngrids)    ;
RPON=zeros(ntimes,nlayers,ngrids)    ;

NO3=zeros(ntimes,nlayers,ngrids)    ;
NH4=zeros(ntimes,nlayers,ngrids)    ;

WC_CDON1=zeros(ntimes,nlayers,ngrids)    ;
WC_NCDON1=zeros(ntimes,nlayers,ngrids)    ;
WC_CDON2=zeros(ntimes,nlayers,ngrids)    ;
WC_NCDON2=zeros(ntimes,nlayers,ngrids)    ;
WC_CDON3=zeros(ntimes,nlayers,ngrids)    ;
WC_NCDON3=zeros(ntimes,nlayers,ngrids)    ;

toc
disp('allocated memory for all the vars');

 tic
 disp('assigning all values')
 %Nitrogen
 disp('Getting sediment data')
    JNH4=nc{'JNH4'}(:);    
    marsh.JNH4.all=JNH4(2:end,marsh.nodes);  
    estuary.JNH4.all=JNH4(2:end,estuary.nodesonly);   
    clear JNh4
    
    JNO3=nc{'JNO3'}(:);   
    marsh.JNO3.all=JNO3(2:end,marsh.nodes);  
    estuary.JNO3.all=JNO3(2:end,estuary.nodesonly);  
    clear JNO3
    
    JPON=nc{'JNIN'}(:);
    marsh.JPON.all=JPON(2:end,marsh.nodes);  
    estuary.JPON.all=JPON(2:end,estuary.nodesonly);
    clear JPON
    
    JDON1=nc{'JWCDON1'}(:)+nc{'JWNCDON1'}(:);
    marsh.JDON1.all=JDON1(2:end,marsh.nodes);  
    estuary.JDON1.all=JDON1(2:end,estuary.nodesonly);
     clear JDON1
     
    JDON2=nc{'JWCDON2'}(:)+nc{'JWNCDON2'}(:); 
    marsh.JDON2.all=JDON2(2:end,marsh.nodes);  
    estuary.JDON2.all=JDON2(2:end,estuary.nodesonly);  
    clear JDON2

    JDON3=nc{'JWCDON3'}(:)+nc{'JWNCDON3'}(:);     
    marsh.JDON3.all=JDON3(2:end,marsh.nodes);  
    estuary.JDON3.all=JDON3(2:end,estuary.nodesonly); 
    clear JDON3
     
%set up sediment variables for the budget

    for i = 1 : length(marsh.area);
        
        marsh.JNH4.flux(:,i)=marsh.JNH4.all(:,i).*marsh.area(i)./1E3 % multiple by the area of each node and convert to kg, now in kg /node d^-1
        marsh.JNO3.flux(:,i)=marsh.JNO3.all(:,i).*marsh.area(i)./1E3;
        marsh.JPON.flux(:,i)=marsh.JPON.all(:,i).*marsh.area(i)./1E3;
        marsh.JDON1.flux(:,i)=marsh.JDON1.all(:,i).*marsh.area(i)./1E3;   % DOC and DON are in units of g m^-2 d^-1, convert to kg / node/ day
        marsh.JDON2.flux(:,i)=marsh.JDON2.all(:,i).*marsh.area(i)./1E3;
        marsh.JDON3.flux(:,i)=marsh.JDON3.all(:,i).*marsh.area(i)./1E3;

    end
    
    for i = 1 : length(estuary.areaonly);
        
        estuary.JNH4.flux(:,i)=estuary.JNH4.all(:,i).*estuary.areaonly(i)./1E3; % multiple by the area of each node and convert to kg, now in kg /node d^-1
        estuary.JNO3.flux(:,i)=estuary.JNO3.all(:,i).*estuary.areaonly(i)./1E3;
        estuary.JPON.flux(:,i)=estuary.JPON.all(:,i).*estuary.areaonly(i)./1E3;
         estuary.JDON1.flux(:,i)=estuary.JDON1.all(:,i).*estuary.areaonly(i)./1E3;   % DOC and DON are in units of g m^-2 d^-1, convert to kg / node/ day
        estuary.JDON2.flux(:,i)=estuary.JDON2.all(:,i).*estuary.areaonly(i)./1E3;
        estuary.JDON3.flux(:,i)=estuary.JDON3.all(:,i).*estuary.areaonly(i)./1E3;
 
    end


    %now add up all the sediment fluxes to get the total coming across the
    %marsh and coming across the estuary.  the net estuary will subtract
    %the marsh
    
    marsh.JNH4.total1=nansum(marsh.JNH4.flux,2);% total flux per time step 
    marsh.JNH4.total2=nansum(marsh.JNH4.total1.*ICM_int); % total flux over entire run
    marsh.JNH4.areal=marsh.JNH4.total2./nansum(marsh.area)/ndays*1E6;
    marsh.JNO3.total1=nansum(marsh.JNO3.flux,2);
    marsh.JNO3.total2=nansum(marsh.JNO3.total1.*ICM_int);
    marsh.JNO3.areal=marsh.JNO3.total2./nansum(marsh.area)/ndays*1E6;
    
    marsh.JPON.total1=nansum(marsh.JPON.flux,2);% total flux per time step 
    marsh.JPON.total2=nansum(marsh.JPON.total1.*ICM_int); % total flux over entire run
    marsh.JPON.areal=marsh.JPON.total2./nansum(marsh.area)/ndays*1E6;
        
    marsh.JDON1.total1=nansum(marsh.JDON1.flux,2);% total flux per time step 
    marsh.JDON1.total2=nansum(marsh.JDON1.total1.*ICM_int); % total flux over entire run 
    marsh.JDON1.areal=marsh.JDON1.total2./nansum(marsh.area)/ndays*1E6;% mg N m^-s d^-1
    marsh.JDON2.total1=nansum(marsh.JDON2.flux,2);% total flux per time step 
    marsh.JDON2.total2=nansum(marsh.JDON2.total1.*ICM_int); % total flux over entire run 
    marsh.JDON2.areal=marsh.JDON1.total2./nansum(marsh.area)/ndays*1E6; % mg N m^-s d^-1   
    marsh.JDON3.total1=nansum(marsh.JDON3.flux,2);% total flux per time step 
    marsh.JDON3.total2=nansum(marsh.JDON3.total1.*ICM_int); % total flux over entire run
    marsh.JDON3.areal=marsh.JDON3.total2./nansum(marsh.area)/ndays*1E6; % mg N m^-s d^-1   
    marsh.JTDON.total1=marsh.JDON1.total1+marsh.JDON2.total1+marsh.JDON3.total1;
    marsh.JTDON.total2=marsh.JDON1.total2+marsh.JDON2.total2+marsh.JDON3.total2;  
    marsh.JTDON.areal=(marsh.JDON1.total2+marsh.JDON2.total2+marsh.JDON3.total2)/nansum(marsh.area)/ndays*1E6;  % mg N m^-s d^-1   
    
   % now estuary sediment water column flux, must subtract the marsh
    estuary.JNH4.total1=nansum(estuary.JNH4.flux,2);%-marsh.JNH4.total1;% total flux per time step 
    estuary.JNH4.total2=nansum(estuary.JNH4.total1.*ICM_int); % total flux over entire run
    estuary.JNH4.areal=estuary.JNH4.total2./nansum(estuary.areaonly)/ndays*1E6;
    
    estuary.JNO3.total1=nansum(estuary.JNO3.flux,2);%;-marsh.JNO3.total1;
    estuary.JNO3.total2=nansum(estuary.JNO3.total1.*ICM_int);
    estuary.JNO3.areal=estuary.JNO3.total2./nansum(estuary.areaonly)/ndays*1E6;
        
    estuary.JPON.total1=nansum(estuary.JPON.flux,2);%;-marsh.JPON.total1;% total flux per time step 
    estuary.JPON.total2=nansum(estuary.JPON.total1.*ICM_int); % total flux over entire run
    estuary.JPON.areal=estuary.JPON.total2./nansum(estuary.areaonly)/ndays*1E6;
    
    estuary.JDON1.total1=nansum(estuary.JDON1.flux,2);%-marsh.JDON1.total1;% total flux per time step 
    estuary.JDON1.total2=nansum(estuary.JDON1.total1.*ICM_int); % total flux over entire run 
    estuary.JDON1.areal=estuary.JDON1.total2./nansum(estuary.areaonly)/ndays*1E6;
    
    estuary.JDON2.total1=nansum(estuary.JDON2.flux,2)-marsh.JDON2.total1;% total flux per time step 
    estuary.JDON2.total2=nansum(estuary.JDON2.total1.*ICM_int); % total flux over entire run 
    estuary.JDON2.areal=estuary.JDON1.total2./nansum(estuary.areaonly)/ndays*1E6;   
    
    estuary.JDON3.total1=nansum(estuary.JDON3.flux,2);%-marsh.JDON3.total1;% total flux per time step 
    estuary.JDON3.total2=nansum(estuary.JDON3.total1.*ICM_int); % total flux over entire run
    estuary.JDON3.areal=estuary.JDON3.total2./nansum(estuary.areaonly)/ndays*1E6;
    
    estuary.JTDON.total1=estuary.JDON1.total1+estuary.JDON2.total1+estuary.JDON3.total1;%-marsh.JTDON.total1;
    estuary.JTDON.total2=estuary.JDON1.total2+estuary.JDON2.total2+estuary.JDON3.total2;  
    estuary.JTDON.areal=(estuary.JDON1.total2+estuary.JDON2.total2+estuary.JDON3.total2)/nansum(estuary.areaonly)/ndays/1E6; 
    

    disp('got Sediment data');
    toc

    tic
    disp('getting nitrogen data');

    LPON=nc{'LPON'}(:);
    RPON= nc{'RPON'}(:);
    
    NH4= nc{'NH4'}(:);
    NO3= nc{'NO3'}(:);
    
    WC_CDON1=nc{'WC_CDON1'}(:);
    WC_NCDON1=nc{'WC_NCDON1'}(:);
    WC_CDON2=nc{'WC_CDON2'}(:);
    WC_NCDON2=nc{'WC_NCDON2'}(:);
    WC_CDON3=nc{'WC_CDON3'}(:);
    WC_NCDON3=nc{'WC_NCDON3'}(:);
    
    toc   
  % open up hisdata files with the integrated budget output
disp('getting the budget data from hisdata');
for ilay =  1 : nlayers;
        nlayer =ilay
    
    load(['hisdata_' num2str(nlayer,'%05d') '.mat'])

    %these are already integrated now, units of g N m^-3
    algdon(ilay,:)=abs(hisdata(end,estuary.nodes,(strcmp(varnames,'AlgDON'))));
    
    algpon(ilay,:)=abs(hisdata(end,estuary.nodes,(strcmp(varnames,'AlgPON'))));
    algnh4(ilay,:)=abs(hisdata(end,estuary.nodes,(strcmp(varnames,'AlgNH4'))));
    algno3(ilay,:)=abs(hisdata(end,estuary.nodes,(strcmp(varnames,'AlgNO3')))); 
    denno3(ilay,:)=abs(hisdata(end,estuary.nodes,(strcmp(varnames,'DenitNO3'))));
    nt(ilay,:)=abs(hisdata(end,estuary.nodes,(strcmp(varnames,'NitrifNH4'))));
    denitn(ilay,:)=abs(hisdata(end,estuary.nodes,(strcmp(varnames,'DenitDON'))));
    hdrlpon(ilay,:)=abs(hisdata(end,estuary.nodes,(strcmp(varnames,'HDRLPON'))));
    hdrrpon(ilay,:)=abs(hisdata(end,estuary.nodes,(strcmp(varnames,'HDRRPON'))));
    coagn(ilay,:)=abs(hisdata(end,estuary.nodes,(strcmp(varnames,'COAGN'))));
    mnldon1(ilay,:)=abs(hisdata(end,estuary.nodes,(strcmp(varnames,'MNLDON1'))));
    mnldon2(ilay,:)=abs(hisdata(end,estuary.nodes,(strcmp(varnames,'MNLDON2'))));
    mnldon3(ilay,:)=abs(hisdata(end,estuary.nodes,(strcmp(varnames,'MNLDON3'))));

    dpd30N(ilay,:)=abs(hisdata(end,estuary.nodes,(strcmp(varnames,'DPD30N'))));

    dpd20N(ilay,:)=abs(hisdata(end,estuary.nodes,(strcmp(varnames,'DPD20N'))));

    dpd10N(ilay,:)=abs(hisdata(end,estuary.nodes,(strcmp(varnames,'DPD10N'))));

end   
disp('now filling up the structure arrays with full water column data')
for ilay=1:nlayers;    
    ilay
%     %%% Nitrogen % all units are g N m^-3 d^-1
    
    estuary.lpon(ilay).all=(squeeze(LPON(2:end,ilay,estuary.nodes)));
    estuary.rpon(ilay).all=(squeeze(RPON(2:end,ilay,estuary.nodes)));
    
    estuary.dtrpon(ilay).all=diff((squeeze(RPON(1:end,ilay,estuary.nodes))),1,1);
    estuary.dtlpon(ilay).all=diff((squeeze(LPON(1:end,ilay,estuary.nodes))),1,1);

    estuary.dtno3(ilay).all=diff((squeeze(NO3(1:end,ilay,estuary.nodes))),1,1);
    
    estuary.dtnh4(ilay).all=diff((squeeze(NH4(1:end,ilay,estuary.nodes))),1,1);

    
    %%% % all units are g C m^-3 d^-1
    
     estuary.dtdon(ilay).all=diff((squeeze(WC_CDON1(1:end,ilay,estuary.nodes)))+ ...
                                  (squeeze(WC_NCDON1(1:end,ilay,estuary.nodes))) + ...
                                (squeeze(WC_CDON2(1:end,ilay,estuary.nodes)))+ ...
                                  (squeeze(WC_NCDON2(1:end,ilay,estuary.nodes))) + ...
                               (squeeze(WC_CDON3(1:end,ilay,estuary.nodes)))+ ...
                                  (squeeze(WC_NCDON3(1:end,ilay,estuary.nodes))),1,1);
                            
     estuary.dtdon1(ilay).all=diff((squeeze(WC_CDON1(1:end,ilay,estuary.nodes)))+ ...
                                  (squeeze(WC_NCDON1(1:end,ilay,estuary.nodes))),1,1);
    estuary.dtdon2(ilay).all=diff((squeeze(WC_CDON2(1:end,ilay,estuary.nodes)))+ ...
                                  (squeeze(WC_NCDON2(1:end,ilay,estuary.nodes))),1,1);
    estuary.dtdon3(ilay).all=diff((squeeze(WC_CDON3(1:end,ilay,estuary.nodes)))+ ...
                                  (squeeze(WC_NCDON3(1:end,ilay,estuary.nodes))),1,1);                          

        estuary.depth(ilay).all=zLayer(ilay,estuary.nodes).*(squeeze(DEPTH(2:end,estuary.nodes)));% LAYER THICKNESS
    estuary.volume(ilay).all=estuary.depth(ilay).all.*repmat(estuary.area,1,ntimes-1)';
  
end

% now we can do the budgeting, integrating every node by estuary.volume to
% get the mass d^-1, question is what is the best way to represent the data
% will also integrate over the entire model period to get total mass   
% NH4
%
%marsh sediment
budgets.NH4(1).name='Marsh JNH4';
budgets.NH4(1).units='kg N';
budgets.NH4(1).value=marsh.JNH4.total2;
budgets.NH4(1).ts=marsh.JNH4.total1;
budgets.NH4(1).type='Source';
%estuary sediment
budgets.NH4(2).name='Estuary JNH4';
budgets.NH4(2).units='kg N';
budgets.NH4(2).value=estuary.JNH4.total2;
budgets.NH4(2).ts=estuary.JNH4.total1;
budgets.NH4(2).type='Source';
%algae production
budgets.NH4(3).name='Estuary AlgNH4';
budgets.NH4(3).units='kg N';
budgets.NH4(3).ts=zeros(ntimes-1,1);
budgets.NH4(3).type='Sink';
%nitrification loss of NH4
budgets.NH4(4).name='Nitrification Loss';
budgets.NH4(4).units= 'kg N';
budgets.NH4(4).ts=zeros(ntimes-1,1);
budgets.NH4(4).type='Sink';
%riverine input
budgets.NH4(5).name='Riverine Inputs';
budgets.NH4(5).units= 'kg N';
budgets.NH4(5).ts=river.NH4.ts;
budgets.NH4(5).value=river.NH4.int;% riverine inputs
budgets.NH4(5).type='Source';
%Remineralization of DON
budgets.NH4(6).name='Remineralization of DON';
budgets.NH4(6).units= 'kg N';
budgets.NH4(6).ts=zeros(ntimes-1,1);
budgets.NH4(6).type='Source';
%PhotoRemineralization of DON
budgets.NH4(7).name='PhotoRemineralization of DON';
budgets.NH4(7).units= 'kg N';
budgets.NH4(7).ts=zeros(ntimes-1,1);%budgets.TDOC(7).ts*0.075;
budgets.NH4(7).value=0.0;%budgets.TDOC(7).value*0.075;
budgets.NH4(7).type='Source';
%denitrification of DON
budgets.NH4(8).name='Denitrification of DON';
budgets.NH4(8).units='kg N';
budgets.NH4(8).ts=zeros(ntimes-1,1);
budgets.NH4(8).type='Source';
%Delta NH4
budgets.NH4(9).name='Delta NH4';
budgets.NH4(9).units='kg N';
budgets.NH4(9).ts=zeros(ntimes-1,1);
budgets.NH4(9).type='Source';

% have to add up all the layers
for ilay=1:nlayers;    
    %algae
    budgets.NH4(3).layers(ilay)=nansum(algnh4(ilay,:).*estuary.volume(ilay).all(end,:),2)/1000;
    %Nitrification of NH4--> NO3
    budgets.NH4(4).layers(ilay)=nansum(nt(ilay,:).*estuary.volume(ilay).all(end,:),2)/1000;                              
    % remin of DON ---> NH4                                
    budgets.NH4(6).layers(ilay)=nansum((mnldon1(ilay,:) + mnldon2(ilay,:)...
                                            +mnldon3(ilay,:)).*estuary.volume(ilay).all(end,:),2)/1000;
     %photochemical remineralization
    budgets.NH4(7).layers(ilay)=nansum((dpd30N(ilay,:)+dpd20N(ilay,:)+dpd10N(ilay,:)).*estuary.volume(ilay).all(end,:),2)/1000;                                         
        %coagulation loss of DOC      
                                        
    % remin of DON ---> NH4 via denitrification                              
    budgets.NH4(8).layers(ilay)=nansum(denitn(ilay,:).*estuary.volume(ilay).all(end,:),2)/1000;
    
    budgets.NH4(9).ts=budgets.NH4(9).ts+nansum(estuary.dtnh4(ilay).all.*estuary.volume(ilay).all,2)/1000;
    
end

 budgets.NH4(3).value=nansum(budgets.NH4(3).layers);
 budgets.NH4(4).value=nansum(budgets.NH4(4).layers);
 budgets.NH4(6).value=nansum(budgets.NH4(6).layers);
 budgets.NH4(7).value=nansum(budgets.NH4(7).layers);
 budgets.NH4(8).value=nansum(budgets.NH4(8).layers);
 budgets.NH4(9).value=nansum(budgets.NH4(9).ts.*ICM_int);

myvals=[];
mynames={};
myts=[];
 %change the sign and extract the data for plotting
for i = 1 : length(budgets.NH4);
    if(strcmp(budgets.NH4(i).type,'Source'));
      myvals(i)=budgets.NH4(i).value;
      myts(i,:)=budgets.NH4(i).ts;
    else
      myvals(i)=budgets.NH4(i).value*-1;
      myts(i,:)=budgets.NH4(i).ts*-1;
    end
    hold on
    mynames{i}=budgets.NH4(i).name;
end

myvals(length(budgets.NH4)+1)=myvals(9)-nansum(myvals(1:8));
mynames{length(budgets.NH4)+1}='MainStem'
mypercent=[];
sinks=myvals(myvals(1:end-2)<0);
sources=myvals(myvals(1:end-2)>=0);

for i = 1: length(myvals);
    mypercent(i)=myvals(i)./nansum(abs(sources))
    mynames{i}=[mynames{i} ' ' num2str(round(mypercent(i)*100,3)) '%']
end

myts(length(budgets.NH4)+1,:)=myts(9,:)-nansum(myts(1:8,:),1);
myts_out.NH4=myts;

myvals_out.NH4=myvals;
mynames_out.NH4=mynames;
% make vector of ones to blow up the pie graph
explode1=ones(1,length(myvals(myvals>0)));
explode2=ones(1,length(myvals(myvals<0)));
figure;
pie(myvals(myvals>0),explode1,mynames(myvals>0));
saveas(gcf,'budget_plots/NH4_budget1','png');
saveas(gcf,'budget_plots/NH4_budget1','fig');
saveas(gcf,'budget_plots/NH4_budget1','eps');

figure;
pie(abs(myvals(myvals<0)),explode2,mynames(myvals<0));
saveas(gcf,'budget_plots/NH4_budget2','png');
saveas(gcf,'budget_plots/NH4_budget2','fig');
saveas(gcf,'budget_plots/NH4_budget2','eps');

figure;
b=bar(myvals,'FaceColor','flat')
xticklabels(mynames)
xtickangle(45)
b.CData(myvals<0,:)=repmat([1 0 1],length(find(myvals<0)),1);
ylabel('NH4 (kg N)');
saveas(gcf,'budget_plots/NH4_budget3','png');
saveas(gcf,'budget_plots/NH4_budget3','fig');
saveas(gcf,'budget_plots/NH4_budget3','eps');
colormap hsv;
cmap=colormap;

[B,Ids]=sort(abs(myvals),'descend');
figure;

color_Vector=1:round(64/length(myvals)):64;
for i = 1 : length(myvals)
     mycolor=cmap(color_Vector(i),:);
     area(time_vector,(myts(Ids(i),:))','FaceAlpha',0.4,'EdgeAlpha',0.0,'FaceColor',mycolor);
     hold on
     
end
legend(mynames{Ids});
xlabel('Days in 2005');ylabel('NH4 (g N d^-^1)');

saveas(gcf,'budget_plots/NH4_budget4','png');
saveas(gcf,'budget_plots/NH4_budget4','fig');
saveas(gcf,'budget_plots/NH4_budget4','eps');
%
% NO3
%marsh sediment flux
budgets.NO3(1).name='Marsh JNO3'
budgets.NO3(1).units='kg N';
budgets.NO3(1).value=marsh.JNO3.total2;
budgets.NO3(1).ts=marsh.JNO3.total1;
budgets.NO3(1).type='Source';
%estuary sediment flux
budgets.NO3(2).name='Estuary JNO3'
budgets.NO3(2).units='kg N';
budgets.NO3(2).value=estuary.JNO3.total2;
budgets.NO3(2).ts=estuary.JNO3.total1;
budgets.NO3(2).type='Source';
%algae uptake
budgets.NO3(3).name='Estuary AlgNO3';
budgets.NO3(3).units='kg N';
budgets.NO3(3).ts=zeros(ntimes-1,1);
budgets.NO3(3).type='Sink';
%nitrification Production of NO3
budgets.NO3(4).name='Nitrification Production';
budgets.NO3(4).units= 'kg N';
budgets.NO3(4).ts=zeros(ntimes-1,1);
budgets.NO3(4).type='Source';
%riverine input
budgets.NO3(5).name='Riverine Inputs';
budgets.NO3(5).units= 'kg N';
budgets.NO3(5).ts=river.NO3.ts;
budgets.NO3(5).value=river.NO3.int;% riverine inputs
budgets.NO3(5).type='Source';
%denitrification of NO3
budgets.NO3(6).name='Denitrification of NO3';
budgets.NO3(6).units='kg N';
budgets.NO3(6).ts=zeros(ntimes-1,1);
budgets.NO3(6).type='Sink';
%Delta NO3
budgets.NO3(7).name='Delta NO3';
budgets.NO3(7).units='kg N';
budgets.NO3(7).ts=zeros(ntimes-1,1);
budgets.NO3(7).type='Source';


% have to add up all the layers
for ilay=1:nlayers;    
    %algae
    budgets.NO3(3).layers(ilay)=nansum(algno3(ilay,:).*estuary.volume(ilay).all(end,:),2)/1000;
    %Nitrification of NH4--> NO3
    budgets.NO3(4).layers(ilay)=nansum(nt(ilay,:).*estuary.volume(ilay).all(end,:),2)/1000;
    %Denitrification of NO3--> N2
    budgets.NO3(6).layers(ilay)=nansum(denno3(ilay,:).*estuary.volume(ilay).all(end,:),2)/1000; 
    
    budgets.NO3(7).ts=budgets.NO3(7).ts+nansum(estuary.dtno3(ilay).all.*estuary.volume(ilay).all,2)/1000; 
end

 budgets.NO3(3).value=nansum(budgets.NO3(3).layers);
 budgets.NO3(4).value=nansum(budgets.NO3(4).layers);
 budgets.NO3(6).value=nansum(budgets.NO3(6).layers);
 budgets.NO3(7).value=nansum(budgets.NO3(7).layers);
 
 
myvals=[];
mynames={};
myts=[];

 %change the sign and extract the data for plotting
for i = 1 : length(budgets.NO3);
    if(strcmp(budgets.NO3(i).type,'Source'));
      myvals(i)=budgets.NO3(i).value;
      myts(i,:)=budgets.NO3(i).ts
    else
      myvals(i)=budgets.NO3(i).value*-1;
      myts(i,:)=budgets.NO3(i).ts*-1;
    end
    hold on
    mynames{i}=budgets.NO3(i).name;
      
end

figure;
myvals(length(budgets.NO3)+1)=myvals(7)-nansum(myvals(1:6));
mynames{length(budgets.NO3)+1}='MainStem'
mypercent=[];
sinks=myvals(myvals(1:end-2)<0);
sources=myvals(myvals(1:end-2)>=0);
mynames_out.NO3=mynames;
for i = 1: length(myvals);
    mypercent(i)=myvals(i)./nansum(abs(sources))
    mynames{i}=[mynames{i} ' ' num2str(round(mypercent(i)*100,3)) '%']
end

myts(length(budgets.NO3)+1,:)=myts(7,:)-nansum(myts(1:6,:),1);
myvals_out.NO3=myvals;
myts_out.NO3=myts;

% make vector of ones to blow up the pie graph
explode1=ones(1,length(myvals(myvals>0)));
explode2=ones(1,length(myvals(myvals<0)));
figure;
pie(myvals(myvals>0),explode1,mynames(myvals>0));
saveas(gcf,'budget_plots/NO3_budget1','png');
saveas(gcf,'budget_plots/NO3_budget1','fig');
saveas(gcf,'budget_plots/NO3_budget1','eps');
figure;
pie(abs(myvals(myvals<0)),explode2,mynames(myvals<0));
saveas(gcf,'budget_plots/NO3_budget2','png');
saveas(gcf,'budget_plots/NO3_budget2','fig');
saveas(gcf,'budget_plots/NO3_budget2','eps');

figure;
b=bar(myvals,'FaceColor','flat')
xticklabels(mynames)
xtickangle(45)
b.CData(myvals<0,:)=repmat([1 0 1],length(find(myvals<0)),1);
ylabel('NO3 (kg N)');
saveas(gcf,'budget_plots/NO3_budget3','png');
saveas(gcf,'budget_plots/NO3_budget3','fig');
saveas(gcf,'budget_plots/NO3_budget3','eps');
colormap hsv;
cmap=colormap;

[B,Ids]=sort(abs(myvals),'descend');
figure;
color_Vector=1:round(64/length(myvals)):64;
for i = 1 : length(myvals)
     mycolor=cmap(color_Vector(i),:);
     area(time_vector,(myts(Ids(i),:))','FaceAlpha',0.4,'EdgeAlpha',0.0,'FaceColor',mycolor);
     hold on
     
end
legend(mynames{Ids});
xlabel('Days in 2005');ylabel('NO3 (g N d^-^1)');

saveas(gcf,'budget_plots/NO3_budget4','png');
saveas(gcf,'budget_plots/NO3_budget4','fig');
saveas(gcf,'budget_plots/NO3_budget4','eps');
%
%DON first
%marsh sediment
budgets.TDON(1).name='Marsh JDON';
budgets.TDON(1).units='kg C';
budgets.TDON(1).value=marsh.JTDON.total2;
budgets.TDON(1).ts=marsh.JTDON.total1;
budgets.TDON(1).type='Source';
%estuary sediment
budgets.TDON(2).name='Estuary JDON';
budgets.TDON(2).units='kg C';
budgets.TDON(2).value=estuary.JTDON.total2;
budgets.TDON(2).ts=estuary.JTDON.total1;
budgets.TDON(2).type='Source';
%algae production
budgets.TDON(3).name='Estuary AlgDON';
budgets.TDON(3).units='kg C';
budgets.TDON(3).ts=zeros(ntimes-1,1);
budgets.TDON(3).type='Source';
%hydrolysis of POC --> DON
budgets.TDON(4).name='Estuary HdrPON';
budgets.TDON(4).units='kg N';
budgets.TDON(4).ts=zeros(ntimes-1,1);
budgets.TDON(4).type='Source';
%denitrification of DON
budgets.TDON(5).name='Denitrification of DON';
budgets.TDON(5).units='kg N';
budgets.TDON(5).ts=zeros(ntimes-1,1);
budgets.TDON(5).type='Sink';
%remineralization of DON
budgets.TDON(6).name='Remineralization of DON';
budgets.TDON(6).units='kg N';
budgets.TDON(6).ts=zeros(ntimes-1,1);
budgets.TDON(6).type='Sink';
% commented out photoremin because the value is currently not output
%PhotoRemineralization of DON
budgets.TDON(7).name='PhotoRemineralization of DON';
budgets.TDON(7).units= 'kg N';
budgets.TDON(7).ts=zeros(ntimes-1,1);%budgets.TDOC(7).ts*0.075;
budgets.TDON(7).value=0.0;%budgets.TDOC(7).value*0.075;
budgets.TDON(7).type='Sink';
%riverine input
budgets.TDON(8).name='Riverine Inputs';
budgets.TDON(8).units= 'kg N';
budgets.TDON(8).ts=river.TDON.ts;% riverine inputs
budgets.TDON(8).value=river.TDON.int;% riverine inputs
budgets.TDON(8).type='Source';
 %Coagulation of DOC
budgets.TDON(9).name='Coagulation of DON';
budgets.TDON(9).units='kg N';
budgets.TDON(9).ts=zeros(ntimes-1,1);
budgets.TDON(9).type='Sink';
 %Delta DON
budgets.TDON(10).name='Delta DON';
budgets.TDON(10).units='kg N';
budgets.TDON(10).ts=zeros(ntimes-1,1);
budgets.TDON(10).type='Source';

% have to add up all the layers
for ilay=1:nlayers;   
  
    %algae
    budgets.TDON(3).layers(ilay)=nansum(algdon(ilay,:).*estuary.volume(ilay).all(end,:),2)/1000;
    %hydrolysis of POC --> DOC
    budgets.TDON(4).layers(ilay)=nansum((hdrlpon(ilay,:)+hdrrpon(ilay,:))...
                                    .*estuary.volume(ilay).all(end,:),2)/1000;
    %denitrification loss of DOC                     
    budgets.TDON(5).layers(ilay)=nansum(denitn(ilay,:).*estuary.volume(ilay).all(end,:),2)/1000;   
        %photochemical remineralization
    budgets.DOC3(7).layers(ilay)=nansum((dpd30N(ilay,:)+dpd20N(ilay,:)+dpd10N(ilay,:)).*estuary.volume(ilay).all(end,:),2)/1000;                                         
        %coagulation loss of DOC    
    
    budgets.TDON(6).layers(ilay)=nansum((mnldon1(ilay,:)+mnldon2(ilay,:)...
                                                        +mnldon3(ilay,:)).*estuary.volume(ilay).all(end,:),2)/1000;

     budgets.TDON(9).layers(ilay)=nansum(coagn(ilay).*estuary.volume(ilay).all(end,:),2)/1000;
    
     budgets.TDON(10).ts=budgets.TDON(10).ts+nansum(estuary.dtdon(ilay).all.*estuary.volume(ilay).all,2)/1000;  

% %                                           
end



 budgets.TDON(3).value=nansum(budgets.TDON(3).layers);
 budgets.TDON(4).value=nansum(budgets.TDON(4).layers);
 budgets.TDON(5).value=nansum(budgets.TDON(5).layers);
 budgets.TDON(6).value=nansum(budgets.TDON(6).layers);
 budgets.TDON(9).value=nansum(budgets.TDON(9).layers);
budgets.TDON(10).value=nansum(budgets.TDON(10).ts.*ICM_int);
 

 mynames={};
 myvals=[];
 myts=[];
 
 %change the sign and extract the data for plotting
for i = 1 : length(budgets.TDON);
    if(strcmp(budgets.TDON(i).type,'Source'));
      myvals(i)=budgets.TDON(i).value;
      myts(i,:)=budgets.TDON(i).ts;
    else
      myvals(i)=budgets.TDON(i).value*-1;
      myts(i,:)=budgets.TDON(i).ts*-1;
    end
    hold on
    mynames{i}=budgets.TDON(i).name;
end


figure;
myvals(length(budgets.TDON)+1)=myvals(10)-nansum(myvals(1:9));
mynames{length(budgets.TDON)+1}='MainStem'
mypercent=[];
sinks=myvals(myvals(1:end-2)<0);
sources=myvals(myvals(1:end-2)>=0);
mynames_out.DON=mynames;

for i = 1: length(myvals);
    mypercent(i)=myvals(i)./nansum(abs(sources))
    mynames{i}=[mynames{i} ' ' num2str(round(mypercent(i)*100,3)) '%']
end

myts(length(budgets.TDON)+1,:)=myts(10,:)-nansum(myts(1:9,:),1);
myts_out.DON=myts;
myvals_out.DON=myvals;

% make vector of ones to blow up the pie graph
explode1=ones(1,length(myvals(myvals>0)));
explode2=ones(1,length(myvals(myvals<0)));
figure;
pie(myvals(myvals>0),explode1,mynames(myvals>0));

saveas(gcf,'budget_plots/DON_budget1','png');
saveas(gcf,'budget_plots/DON_budget1','fig');
saveas(gcf,'budget_plots/DON_budget1','eps');

figure;
pie(abs(myvals(myvals<0)),explode2,mynames(myvals<0));

saveas(gcf,'budget_plots/DON_budget2','png');
saveas(gcf,'budget_plots/DON_budget2','fig');
saveas(gcf,'budget_plots/DON_budget2','eps');

figure;
b=bar(myvals,'FaceColor','flat')
xticklabels(mynames)
xtickangle(45)
b.CData(myvals<0,:)=repmat([1 0 1],length(find(myvals<0)),1);
ylabel('DON (kg N)');

saveas(gcf,'budget_plots/DON_budget3','png');
saveas(gcf,'budget_plots/DON_budget3','fig');
saveas(gcf,'budget_plots/DON_budget3','eps');
colormap hsv;
cmap=colormap;
[B,Ids]=sort(abs(myvals),'descend');
figure;
color_Vector=1:round(64/length(myvals)):64;
for i = 1 : length(myvals)
     mycolor=cmap(color_Vector(i),:);
     h=area(time_vector,(myts(Ids(i),:))','FaceAlpha',0.4,'EdgeAlpha',0.0,'FaceColor',mycolor);
     hold on
     
end
legend(mynames{Ids});
xlabel('Days in 2005');ylabel('DON (g N d^-^1)');

saveas(gcf,'budget_plots/DON_budget4','png');
saveas(gcf,'budget_plots/DON_budget4','fig');
saveas(gcf,'budget_plots/DON_budget4','eps');
%
% PON
%
%marsh sediment
budgets.PON(1).name='Marsh JPON';
budgets.PON(1).units='kg N';
budgets.PON(1).value=marsh.JPON.total2;
budgets.PON(1).ts=marsh.JPON.total1;
budgets.PON(1).type='Sink';
%estuary sediment
budgets.PON(2).name='Estuary JPON';
budgets.PON(2).units='kg N';
budgets.PON(2).value=estuary.JPON.total2;
budgets.PON(2).ts=estuary.JPON.total1;
budgets.PON(2).type='Sink';
%algae production
budgets.PON(3).name='Estuary AlgPON';
budgets.PON(3).units='kg N';
budgets.PON(3).ts=zeros(ntimes-1,1);
budgets.PON(3).type='Source';
%hydrolysis of PON --> DON
budgets.PON(4).name='Estuary HdrPON';
budgets.PON(4).units='kg N';
budgets.PON(4).ts=zeros(ntimes-1,1);
budgets.PON(4).type='Sink';
%riverine input
budgets.PON(5).name='Riverine Inputs';
budgets.PON(5).units= 'kg C';
budgets.PON(5).ts=river.LPON.ts+river.RPON.ts;% riverine inputs
budgets.PON(5).value=river.LPON.int+river.RPON.int;% riverine inputs
budgets.PON(5).type='Source';
%coagulation of DON
budgets.PON(6).name='Coagulation DON'
budgets.PON(6).units = 'kg N'
budgets.PON(6).ts=budgets.TDON(9).ts;
budgets.PON(6).value=budgets.TDON(9).value
budgets.PON(6).type='Source'
%Delta PON
budgets.PON(7).name='Delta PON'
budgets.PON(7).units = 'kg N'
budgets.PON(7).ts=zeros(ntimes-1,1)
budgets.PON(7).type='Source'

% have to add up all the layers
for ilay=1:nlayers;    
    %algae
    budgets.PON(3).layers(ilay)=nansum(algpon(ilay,:).*estuary.volume(ilay).all(end,:))/1000;
    %hydrolysis of PON --> DON
    budgets.PON(4).layers(ilay)=nansum((hdrlpon(ilay,:)+hdrrpon(ilay,:))...
                                    .*estuary.volume(ilay).all(end,:))/1000;
                                
    budgets.PON(7).ts=budgets.PON(7).ts+nansum((estuary.dtlpon(ilay).all+estuary.dtrpon(ilay).all).*estuary.volume(ilay).all,2)/1000;
                                
                                                             
end

 budgets.PON(3).value=nansum(budgets.PON(3).layers);
 budgets.PON(4).value=nansum(budgets.PON(4).layers);
 budgets.PON(7).value=nansum(budgets.PON(7).ts.*ICM_int);
 
myvals=[];
mynames={};
myts=[];
 %change the sign and extract the data for plotting
for i = 1 : length(budgets.PON);
    if(strcmp(budgets.PON(i).type,'Source'));
      myvals(i)=budgets.PON(i).value;
      myts(i,:)=budgets.PON(i).ts;
    else
      myvals(i)=budgets.PON(i).value*-1;
      myts(i,:)=budgets.PON(i).ts*-1;
    end
    hold on
    mynames{i}=budgets.PON(i).name;
end

myvals(length(budgets.PON)+1)=myvals(7)-nansum(myvals(1:6));

mynames{length(budgets.PON)+1}='MainStem'
mypercent=[];
sinks=myvals(myvals(1:end-2)<0);
sources=myvals(myvals(1:end-2)>=0);

mynames_out.PON=mynames;

for i = 1: length(myvals);
    mypercent(i)=myvals(i)./nansum(abs(sources))
    mynames{i}=[mynames{i} ' ' num2str(round(mypercent(i)*100,3)) '%']
end

myts(length(budgets.PON)+1,:)=myts(7,:)-nansum(myts(1:6,:),1);

myvals_out.PON=myvals;
myts_out.PON=myts;
% make vector of ones to blow up the pie graph
explode1=ones(1,length(myvals(myvals>0)));
explode2=ones(1,length(myvals(myvals<0)));
figure;
pie(myvals(myvals>0),explode1,mynames(myvals>0));
saveas(gcf,'budget_plots/PON_budget1','png');
saveas(gcf,'budget_plots/PON_budget1','fig');
saveas(gcf,'budget_plots/PON_budget1','eps');
figure;
pie(abs(myvals(myvals<0)),explode2,mynames(myvals<0));
saveas(gcf,'budget_plots/PON_budget2','png');
saveas(gcf,'budget_plots/PON_budget2','fig');
saveas(gcf,'budget_plots/PON_budget2','eps');

figure;
b=bar(myvals,'FaceColor','flat')
xticklabels(mynames)
xtickangle(45)
b.CData(myvals<0,:)=repmat([1 0 1],length(find(myvals<0)),1);
ylabel('PON(kg N)');
saveas(gcf,'budget_plots/PON_budget3','png');
saveas(gcf,'budget_plots/PON_budget3','fig');
saveas(gcf,'budget_plots/PON_budget3','eps');
colormap hsv;
cmap=colormap;

[B,Ids]=sort(abs(myvals),'descend');
figure;
color_Vector=1:round(64/length(myvals)):64;
for i = 1 : length(myvals)
     mycolor=cmap(color_Vector(i),:);
     area(time_vector,(myts(Ids(i),:))','FaceAlpha',0.4,'EdgeAlpha',0.0,'FaceColor',mycolor);
     hold on
     
end
legend(mynames{Ids});
xlabel('Days in 2005');ylabel('PON (g N  d^-^1)');

saveas(gcf,'budget_plots/PON_budget4','png');
saveas(gcf,'budget_plots/PON_budget4','fig');
saveas(gcf,'budget_plots/PON_budget4','eps');


save('Nitrogen_budget_out.mat','myvals_out','mynames_out','myts_out')


%Finished 

