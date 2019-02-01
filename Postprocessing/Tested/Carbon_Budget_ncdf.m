%Carbon Budget calculation for ICM with netcdf outputs

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
% load('*/RhodeFVCOM_2005/Data/RR_discharge_2005.mat');
load(input('What is the full path and file name to the River discharge file? --> ','s'));
% load('*/RhodeFVCOM_2005/Data/RR_Watershed_WQ_2005.mat');
load(input('What is the full path and file name to the Watershed Water Quality Input data? --> ','s'));

%load the sediment definition file where marsh nodes are flagged with 1
% marsh_file= '*/RhodeFVCOM_2005/Data/seddom_grid.dat';
marsh_file=input('What is the full path and file name to the sediment definition file? --> ','s');
marsh_data=importdata(marsh_file);
% load the grid file
% grid_file='*/RhodeFVCOM_2005/Data/rhoderiver_grid_v1_2.mat';
grid_file=input('What is the full path and file name to the grid definition *.mat file? --> ','s');
load(grid_file);
%load the node area file
% [node_area]=importdata('*/RhodeFVCOM_2005/Data*/RhodeFVCOM_2005/Data/node_CVS_area.dat');%get our nodal area for marsh flux calcs;
node_file=input('What is the full path and file name to the Node Area file? --->  ','s');
[node_area]=importdata(node_file);
%load the file with the different polygons from the Rhode River
%each polygon has a different set of nodes with a decreasing area of the
%Rhode River
% load('*/RhodeFVCOM_2005/Data/polygon_work.mat
load(input('What is the full path and file name of the polygon definition file? --> ','s'));

first_day =input('What is the starting day of the year for this run (i.e. April 1 is 91 ---> ');
%
% Use this flag depending what region to set up the budget for, 1 being the
% whole tributary decreasing size of estuary area to 5
nodes_flag=input('What polygon do you want to build a budget for? 1-5 for trib, 0 for mainstem --->  ')
%get the rest of the information
wqm_hisdir = input('What is full path to the WQM history directory with the netcdf file? --->  ','s');
numfile=input('How many files do you want to create a budget for? ---> ')

nameseg=input('What is the name segment for the FVCOM files ---> ','s');
nlayers=input('How many vertical layers are in the model? ---> ');

%SET UP EMPTY ARRAYS FOR ALL THE BUDGET VARIABLES

fname=[wqm_hisdir '/rhd_0001.nc']

%pull out the vertical coordinate distribution
%this is the FVCOM hydrodynamics file, used to load in the vertical
%coordinate system

nc=netcdf(fname);
s_hybrid=nc{'siglev'}(:);
dz=diff(s_hybrid,1);
for i = 1:nlayers;
    layer_thickness(i,:)=dz(i,:).*xyd_n(:,4)'*-1;
end


myinfo=ncinfo(fname);

ngrids=myinfo.Dimensions(2).Length;
nlayers=myinfo.Dimensions(4).Length
ntimes=myinfo.Dimensions(end).Length;

%%
%how many days was the model run for?

mkdir('budget_plots');
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
    estuary.nodes=RR_nodes_out(find(lld_n(RR_nodes_out,3) > 38.875));
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
% get the nodes

% this is only the subtidal area, excluding the marsh
%use in the sediment flux calculations
estuary.nodesonly=setdiff(estuary.nodes,marsh.nodes);
estuary.areaonly=node_area(estuary.nodesonly,2);

%plot the polygons
tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3),lld_n(:,4));
hold on
plot(lld_n(estuary.nodes,2),lld_n(estuary.nodes,3),'kd')
plot(lld_n(marsh.nodes,2),lld_n(marsh.nodes,3),'ro')
saveas(gcf,'budget_plots/Marsh_map.png')
%
%
zLayer=diff(s_hybrid,1)*-1;% the amount of water column for each layer, adds up to 1;
tic
disp('Getting Depth and Time Information');
% nc = netcdf([fname]);
DEPTH=zeros(ntimes,ngrids);
DEPTH=nc{'depth'}(:);
time_vector_in=nc{'time'}(:);
time_vector=(time_vector_in(2:end)./86400)+first_day;
time_vector_hr=time_vector.*24;
ndays=round(time_vector(end))-first_day;
disp('Got time and depth data')
toc
ICM_int=diff(time_vector_in(1:2))/86400;%get the output time period in days

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

%match up the model times with the discharge times
matching_times=dsearchn(river_time',time_vector_hr);%
%get the discharge times at the model time frequency for each river
cut_discharge=mydischarge(matching_times,1:river_id);
% discharge_Myint=(cut_discharge(2:round(ICM_int*24):end,:));
discharge_Myint=cut_discharge;
river.times=Final_RR_forcing(1:45:end,2)+first_day*24;
%
Final_RR_forcing=Final_RR_forcing(:,1:river_id);

% will take the concentration and multiple my discharge g m^-3 x m^3 s^-1 *
% 3600 secs hr^-1 then integrate to get kg
river.CDOC1.all=interp1(river.times,Final_RR_forcing(10:45:end,:),time_vector_hr,'pchip')
river.NCDOC1.all=interp1(river.times,Final_RR_forcing(11:45:end,:),time_vector_hr,'pchip');
river.LPOC.all=interp1(river.times,Final_RR_forcing(12:45:end,:),time_vector_hr,'pchip');
river.RPOC.all=interp1(river.times,Final_RR_forcing(13:45:end,:),time_vector_hr,'pchip');
river.CDOC2.all=interp1(river.times,Final_RR_forcing(34:45:end,:),time_vector_hr,'pchip');
river.NCDOC2.all=interp1(river.times,Final_RR_forcing(35:45:end,:),time_vector_hr,'pchip');
river.CDOC3.all=interp1(river.times,Final_RR_forcing(36:45:end,:),time_vector_hr,'pchip');
river.NCDOC3.all=interp1(river.times,Final_RR_forcing(37:45:end,:),time_vector_hr,'pchip');

% now integrate and get one time series for each constituent
river.CDOC1.ts=nansum(river.CDOC1.all.*discharge_Myint*86400*ICM_int./1000,2);% kg per ICM_int period
river.CDOC1.int=nansum(river.CDOC1.ts);
river.NCDOC1.ts=nansum(river.NCDOC1.all.*discharge_Myint*86400*ICM_int./1000,2);% kg per ICM_int period
river.NCDOC1.int=nansum(river.NCDOC1.ts);
river.NCDOC2.ts=nansum(river.NCDOC2.all.*discharge_Myint*86400*ICM_int./1000,2);% kg per ICM_int period
river.NCDOC2.int=nansum(river.NCDOC2.ts);
river.CDOC2.ts=nansum(river.CDOC2.all.*discharge_Myint*86400*ICM_int./1000,2);% kg per ICM_int period
river.CDOC2.int=nansum(river.CDOC2.ts);
river.NCDOC3.ts=nansum(river.NCDOC3.all.*discharge_Myint*86400*ICM_int./1000,2);% kg per ICM_int period
river.NCDOC3.int=nansum(river.NCDOC3.ts);
river.CDOC3.ts=nansum(river.CDOC3.all.*discharge_Myint*86400*ICM_int./1000,2);% kg per ICM_int period
river.CDOC3.int=nansum(river.CDOC3.ts);

river.LPOC.ts=nansum(river.LPOC.all.*discharge_Myint*86400*ICM_int./1000,2);% kg per ICM_int period
river.LPOC.int=nansum(river.LPOC.ts);
river.RPOC.ts=nansum(river.RPOC.all.*discharge_Myint*86400*ICM_int./1000,2);% kg per ICM_int period
river.RPOC.int=nansum(river.RPOC.ts);

river.TDOC.ts=river.CDOC1.ts+river.CDOC2.ts+river.CDOC3.ts+river.NCDOC1.ts+river.NCDOC2.ts+river.NCDOC3.ts;
river.TDOC.int=river.CDOC1.int+river.CDOC2.int+river.CDOC3.int+river.NCDOC1.int+river.NCDOC2.int+river.NCDOC3.int;

river.POC.ts=river.LPOC.ts+river.RPOC.ts;
river.POC.int=river.LPOC.int+river.RPOC.int;

toc
disp('got river data')


%%

disp('assigning all values')
tic
disp('Getting Sediment Flux data');
%Sediment-water column organic carbon fluxes
JPOC=nc{'JCIN'}(:);
marsh.JPOC.all=JPOC(2:end,marsh.nodes);
estuary.JPOC.all=JPOC(2:end,estuary.nodesonly);
clear JPOC

JDOC1=nc{'JWCDOC1'}(:)+nc{'JWNCDOC1'}(:);
marsh.JDOC1.all=JDOC1(2:end,marsh.nodes);
estuary.JDOC1.all=JDOC1(2:end,estuary.nodesonly);
clear JDOC1

JDOC2=nc{'JWCDOC2'}(:)+nc{'JWNCDOC2'}(:);
marsh.JDOC2.all=JDOC2(2:end,marsh.nodes);
estuary.JDOC2.all=JDOC2(2:end,estuary.nodesonly);
clear JDOC2

JDOC3=nc{'JWCDOC3'}(:)+nc{'JWNCDOC3'}(:);
marsh.JDOC3.all=JDOC3(2:end,marsh.nodes);
estuary.JDOC3.all=JDOC3(2:end,estuary.nodesonly);
clear JDOC3

     
%set up sediment variables for the budget
    for i = 1 : length(marsh.area);
        marsh.JPOC.flux(:,i)=marsh.JPOC.all(:,i).*marsh.area(i)./1E3;
        marsh.JDOC1.flux(:,i)=marsh.JDOC1.all(:,i).*marsh.area(i)./1E3;   % DOC and DON are in units of g m^-2 d^-1, convert to kg / node/ day
        marsh.JDOC2.flux(:,i)=marsh.JDOC2.all(:,i).*marsh.area(i)./1E3;
        marsh.JDOC3.flux(:,i)=marsh.JDOC3.all(:,i).*marsh.area(i)./1E3;
    end
    
    for i = 1 : length(estuary.areaonly);
        estuary.JPOC.flux(:,i)=estuary.JPOC.all(:,i).*estuary.areaonly(i)./1E3;
        estuary.JDOC1.flux(:,i)=estuary.JDOC1.all(:,i).*estuary.areaonly(i)./1E3;   % DOC and DON are in units of g m^-2 d^-1, convert to kg / node/ day
        estuary.JDOC2.flux(:,i)=estuary.JDOC2.all(:,i).*estuary.areaonly(i)./1E3;
        estuary.JDOC3.flux(:,i)=estuary.JDOC3.all(:,i).*estuary.areaonly(i)./1E3; 
    end


    %now add up all the sediment fluxes to get the total coming across the
    %marsh and coming across the estuary.  the net estuary will subtract
    %the marsh
    
    marsh.JPOC.total1=nansum(marsh.JPOC.flux,2);% total flux per time step 
    marsh.JPOC.total2=nansum(marsh.JPOC.total1.*ICM_int); % total flux over entire run
    marsh.JPOC.areal=marsh.JPOC.total2./nansum(marsh.area)/ndays*1E6; % average areal flux per day mg C m^-2 d^-1
        
    marsh.JDOC1.total1=nansum(marsh.JDOC1.flux,2);% total flux per time step 
    marsh.JDOC1.total2=nansum(marsh.JDOC1.total1.*ICM_int); % total flux over entire run 
    marsh.JDOC1.areal=marsh.JDOC1.total2./nansum(marsh.area)/ndays*1E6; % mg C m^-2 d^-1
    
    marsh.JDOC2.total1=nansum(marsh.JDOC2.flux,2);% total flux per time step 
    marsh.JDOC2.total2=nansum(marsh.JDOC2.total1.*ICM_int); % total flux over entire run 
    marsh.JDOC2.areal=marsh.JDOC2.total2./nansum(marsh.area)/ndays*1E6;% mg C m^-s d^-1
    
    marsh.JDOC3.total1=nansum(marsh.JDOC3.flux,2);% total flux per time step 
    marsh.JDOC3.total2=nansum(marsh.JDOC3.total1.*ICM_int); % total flux over entire run
    marsh.JDOC3.areal=marsh.JDOC3.total2./nansum(marsh.area)/ndays*1E6;% mg C m^-s d^-1
    % total DOC
    marsh.JTDOC.total1=marsh.JDOC1.total1+marsh.JDOC2.total1+marsh.JDOC3.total1;
    marsh.JTDOC.total2=marsh.JDOC1.total2+marsh.JDOC2.total2+marsh.JDOC3.total2;
    marsh.JTDOC.areal=(marsh.JDOC1.total2+marsh.JDOC2.total2+marsh.JDOC3.total2)/nansum(marsh.area)/ndays*1E6;% mg C m^-2 d^-1
     
   % now estuary sediment water column flux, must subtract the marsh
    
    estuary.JPOC.total1=nansum(estuary.JPOC.flux,2);%;-marsh.JPOC.total1;% total flux per time step 
    estuary.JPOC.total2=nansum(estuary.JPOC.total1.*ICM_int); % total flux over entire run
    estuary.JPOC.areal=estuary.JPOC.total2./nansum(estuary.areaonly)/ndays*1E6;
    
    estuary.JDOC1.total1=nansum(estuary.JDOC1.flux,2);%;-marsh.JDOC1.total1;% total flux per time step 
    estuary.JDOC1.total2=nansum(estuary.JDOC1.total1.*ICM_int); % total flux over entire run 
    estuary.JDOC1.areal=estuary.JDOC1.total2./nansum(estuary.areaonly)/ndays*1E6;
    
    estuary.JDOC2.total1=nansum(estuary.JDOC2.flux,2);%-marsh.JDOC2.total1;% total flux per time step
    estuary.JDOC2.total2=nansum(estuary.JDOC2.total1.*ICM_int); % total flux over entire run 
    estuary.JDOC2.areal=estuary.JDOC2.total2./nansum(estuary.areaonly)/ndays*1E6;
    
    estuary.JDOC3.total1=nansum(estuary.JDOC3.flux,2);%-marsh.JDOC3.total1;% total flux per time step 
    estuary.JDOC3.total2=nansum(estuary.JDOC3.total1.*ICM_int); % total flux over entire run
    estuary.JDOC3.areal=estuary.JDOC3.total2./nansum(estuary.areaonly)/ndays*1E6;
    
    estuary.JTDOC.total1=estuary.JDOC1.total1+estuary.JDOC2.total1+estuary.JDOC3.total1;%-marsh.JTDOC.total1;
    estuary.JTDOC.total2=estuary.JDOC1.total2+estuary.JDOC2.total2+estuary.JDOC3.total2;
    estuary.JTDOC.areal=(estuary.JDOC1.total2+estuary.JDOC2.total2+estuary.JDOC3.total2)/nansum(estuary.areaonly)/ndays*1E6;

    disp('got Sediment data');
    toc
%% Now do the water column variables
   tic
    disp('getting carbon data');
    
    LPOC= nc{'LPOC'}(:);
    RPOC= nc{'RPOC'}(:);      
    WC_CDOC1=nc{'WC_CDOC1'}(:);
    WC_NCDOC1=nc{'WC_NCDOC1'}(:);
    WC_CDOC2=nc{'WC_CDOC2'}(:);
    WC_NCDOC2=nc{'WC_NCDOC2'}(:);
    WC_CDOC3=nc{'WC_CDOC3'}(:);
    WC_NCDOC3=nc{'WC_NCDOC3'}(:);    
disp('got carbon data')
    toc
    disp('assigned all variables');
clear nc 

% open up hisdata files with the integrated budget output
disp('getting the budget data from hisdata');
for ilay =  1 : nlayers;
        nlayer =ilay
    
    load(['hisdata_' num2str(nlayer,'%05d') '.mat'])
    
    algdoc(ilay,:)=abs(hisdata(end,estuary.nodes,(strcmp(varnames,'AlgDOC'))));
    algpoc(ilay,:)=abs(hisdata(end,estuary.nodes,(strcmp(varnames,'AlgPOC'))));
    denit(ilay,:)=abs(hisdata(end,estuary.nodes,(strcmp(varnames,'DenitDOC'))));
    hdrlpoc(ilay,:)=abs(hisdata(end,estuary.nodes,(strcmp(varnames,'HDRLPOC'))));
    hdrrpoc(ilay,:)=abs(hisdata(end,estuary.nodes,(strcmp(varnames,'HDRRPOC'))));
    coagc(ilay,:)=abs(hisdata(end,estuary.nodes,(strcmp(varnames,'COAGC'))));
    mnldoc1(ilay,:)=abs(hisdata(end,estuary.nodes,(strcmp(varnames,'MNLDOC1'))));
    mnldoc2(ilay,:)=abs(hisdata(end,estuary.nodes,(strcmp(varnames,'MNLDOC2'))));
    mnldoc3(ilay,:)=abs(hisdata(end,estuary.nodes,(strcmp(varnames,'MNLDOC3'))));
    
%     dpd32(ilay,:)=hisdata(end,:,:,(strcmp(varnames,'DPD32')));
    dpd31(ilay,:)=abs(hisdata(end,estuary.nodes,(strcmp(varnames,'DPD31'))));
    dpd30(ilay,:)=abs(hisdata(end,estuary.nodes,(strcmp(varnames,'DPD30'))));
    dpd3N(ilay,:)=abs(hisdata(end,estuary.nodes,(strcmp(varnames,'DPD3N'))));
    dpd21(ilay,:)=abs(hisdata(end,estuary.nodes,(strcmp(varnames,'DPD21'))));
    dpd20(ilay,:)=abs(hisdata(end,estuary.nodes,(strcmp(varnames,'DPD20'))));
    dpd2N(ilay,:)=abs(hisdata(end,estuary.nodes,(strcmp(varnames,'DPD2N'))));    

    dpd10(ilay,:)=abs(hisdata(end,estuary.nodes,(strcmp(varnames,'DPD10'))));
    dpd1N(ilay,:)=abs(hisdata(end,estuary.nodes,(strcmp(varnames,'DPD1N'))));   
    
    
end  
tic
disp('now filling up the structure arrays with full water column data')

for ilay=1:nlayers;    
    ilay
%Carbon

    %calculate the change in concentration over time, dC/dt (g C m^-3
    %timestep^-1)
    estuary.dtdoc(ilay).all=diff((squeeze(WC_CDOC1(1:end,ilay,estuary.nodes)))+ ...
                                  (squeeze(WC_NCDOC1(1:end,ilay,estuary.nodes))) + ...
                                (squeeze(WC_CDOC2(1:end,ilay,estuary.nodes)))+ ...
                                  (squeeze(WC_NCDOC2(1:end,ilay,estuary.nodes))) + ...
                               (squeeze(WC_CDOC3(1:end,ilay,estuary.nodes)))+ ...
                                  (squeeze(WC_NCDOC3(1:end,ilay,estuary.nodes))),1,1);
                            
     estuary.dtdoc1(ilay).all=diff((squeeze(WC_CDOC1(1:end,ilay,estuary.nodes)))+ ...
                                  (squeeze(WC_NCDOC1(1:end,ilay,estuary.nodes))),1,1);
    estuary.dtdoc2(ilay).all=diff((squeeze(WC_CDOC2(1:end,ilay,estuary.nodes)))+ ...
                                  (squeeze(WC_NCDOC2(1:end,ilay,estuary.nodes))),1,1);
    estuary.dtdoc3(ilay).all=diff((squeeze(WC_CDOC3(1:end,ilay,estuary.nodes)))+ ...
                                  (squeeze(WC_NCDOC3(1:end,ilay,estuary.nodes))),1,1);                              
                                                          

    estuary.lpoc(ilay).all=(squeeze(LPOC(2:end,ilay,estuary.nodes)));
    estuary.rpoc(ilay).all=(squeeze(RPOC(2:end,ilay,estuary.nodes)));
    
    estuary.dtrpoc(ilay).all=diff((squeeze(RPOC(1:end,ilay,estuary.nodes))),1,1);
    estuary.dtlpoc(ilay).all=diff((squeeze(LPOC(1:end,ilay,estuary.nodes))),1,1);
    
    % calculates the volume of water at each node for each time step and
    % each depth, it is time varying due to change in tidal elevation
    estuary.depth(ilay).all=zLayer(ilay,estuary.nodes).*(squeeze(DEPTH(2:end,estuary.nodes)));% LAYER THICKNESS
    estuary.volume(ilay).all=estuary.depth(ilay).all.*repmat(estuary.area,1,ntimes-1)';
  
end
toc

% now we can do the budgeting, integrating every node by estuary.volume to
% get the mass d^-1
% will also integrate over the entire model period to get total mass 
%%
disp('Calculating the Carbon budgets');
% DOC1
%marsh sediment
budgets.DOC1(1).name='Marsh JDOC1';
budgets.DOC1(1).units='kg C';
budgets.DOC1(1).value=marsh.JDOC1.total2;
budgets.DOC1(1).ts=marsh.JDOC1.total1;
budgets.DOC1(1).type='Source';
%estuary sediment
budgets.DOC1(2).name='Estuary JDOC1';
budgets.DOC1(2).units='kg C';
budgets.DOC1(2).value=estuary.JDOC1.total2;
budgets.DOC1(2).ts=estuary.JDOC1.total1;
budgets.DOC1(2).type='Source';
%algae production
budgets.DOC1(3).name='Estuary AlgDOC1';
budgets.DOC1(3).units='kg C';
budgets.DOC1(3).ts=zeros(ntimes-1,1);
budgets.DOC1(3).type='Source';
%hydrolysis of POC --> DOC
budgets.DOC1(4).name='Estuary HdrPOC1';
budgets.DOC1(4).units='kg C';
budgets.DOC1(4).ts=zeros(ntimes-1,1);
budgets.DOC1(4).type='Source';
%denitrification of DOC
budgets.DOC1(5).name='Denitrification of DOC';
budgets.DOC1(5).units='kg C';
budgets.DOC1(5).ts=zeros(ntimes-1,1);
budgets.DOC1(5).type='Sink';
%remineralization of DOC
budgets.DOC1(6).name='Remineralization of DOC';
budgets.DOC1(6).units='kg C';
budgets.DOC1(6).ts=zeros(ntimes-1,1);
budgets.DOC1(6).type='Sink';
%photochemical remineralization
budgets.DOC1(7).name='Photochemical Remin.';
budgets.DOC1(7).units= 'kg C'
budgets.DOC1(7).ts=zeros(ntimes-1,1);
budgets.DOC1(7).type='Sink';
%riverine input
budgets.DOC1(8).name='Riverine Inputs';
budgets.DOC1(8).units= 'kg C';
budgets.DOC1(8).ts=river.CDOC1.ts+river.NCDOC1.ts;% riverine inputs
budgets.DOC1(8).value=river.CDOC1.int+river.NCDOC1.int;% riverine inputs
budgets.DOC1(8).type='Source';
%  %Coagulation of DOC
% budgets.DOC1(9).name='Coagulation of DOC1';
% budgets.DOC1(9).units='kg C';
% budgets.DOC1(9).ts=zeros(ntimes-1,1);
% budgets.DOC1(9).type='Sink';
%photochemical production of CDOC1
budgets.DOC1(9).name='Photochem Prod. DOC1';
budgets.DOC1(9).units= 'kg C';
budgets.DOC1(9).ts=zeros(ntimes-1,1);
budgets.DOC1(9).type='Source';
%photochemical loss of CDOC1
budgets.DOC1(10).name='Photochem Loss DOC1';
budgets.DOC1(10).units= 'kg C';
budgets.DOC1(10).ts=zeros(ntimes-1,1);
budgets.DOC1(10).type='Sink';
%change in concentration of DOC % dt is negative if increases
budgets.DOC1(11).name='Delta DOC1';
budgets.DOC1(11).units= 'kg C';
budgets.DOC1(11).ts=zeros(ntimes-1,1);
budgets.DOC1(11).type='Source';


% have to add up all the layers
% and also convert from grams to kgs
for ilay=1:nlayers;     
    %algae
    budgets.DOC1(3).layers(ilay)=nansum(algdoc(ilay,:).*FCD1.*estuary.volume(ilay).all(end,:),2)/1000;
    %hydrolysis of POC --> DOC
    budgets.DOC1(4).layers(ilay)=nansum((hdrlpoc(ilay,:)+hdrrpoc(ilay,:)).*FHDRC1...
                                    .*estuary.volume(ilay).all(end,:),2)/1000;
    %denitrification loss of DOC                     
    budgets.DOC1(5).layers(ilay)=nansum(denit(ilay,:).*estuary.volume(ilay).all(end,:),2)/1000;   
    
    budgets.DOC1(6).layers(ilay)=nansum(mnldoc1(ilay,:).*estuary.volume(ilay).all(end,:),2)/1000;
    %photochemical remineralization
    budgets.DOC1(7).layers(ilay)=nansum(dpd10(ilay,:).*estuary.volume(ilay).all(end,:),2)/1000;                                         

    %Photochemical production
    budgets.DOC1(9).layers(ilay)=nansum((dpd21(ilay,:)+ ...
                                                    dpd31(ilay,:)...
                                                 +  dpd3N(ilay,:).*FPDC31 ...
                                                 +  dpd2N(ilay,:).*FPDC21...
                                                 +  dpd1N(ilay,:).*FPDC11).*estuary.volume(ilay).all(end,:),2)/1000;
%     
    % %        %photochemical bleaching
     budgets.DOC1(10).layers(ilay)=nansum(dpd1N(ilay,:).*estuary.volume(ilay).all(end,:),2)/1000;        
%      budgets.DOC1(9).ts=budgets.DOC1(9).ts+nansum(estuary.coagc(ilay).all.*estuary.volume(ilay).all,2)/1000;
    
     budgets.DOC1(11).ts=budgets.DOC1(11).ts+nansum(estuary.dtdoc1(ilay).all.*estuary.volume(ilay).all,2)/1000;  
   
% %                                           
end

%integrate everyting now
 budgets.DOC1(3).value=nansum(budgets.DOC1(3).layers);
 budgets.DOC1(4).value=nansum(budgets.DOC1(4).layers);
 budgets.DOC1(5).value=nansum(budgets.DOC1(5).layers);
 budgets.DOC1(6).value=nansum(budgets.DOC1(6).layers);
 budgets.DOC1(7).value=nansum(budgets.DOC1(7).layers);
 budgets.DOC1(9).value=nansum(budgets.DOC1(9).layers);
 budgets.DOC1(10).value=nansum(budgets.DOC1(10).layers);
 budgets.DOC1(11).value=nansum(budgets.DOC1(11).ts.*ICM_int);
 
myvals=[];
mynames={};
myts=[];

 %change the sign and extract the data for plotting
for i = 1 : length(budgets.DOC1);
    if(strcmp(budgets.DOC1(i).type,'Source'));
      myvals(i)=budgets.DOC1(i).value;
      myts(i,:)=budgets.DOC1(i).ts';
    else
      myvals(i)=budgets.DOC1(i).value*-1;
      myts(i,:)=budgets.DOC1(i).ts'*-1;
    end
    mynames{i}=budgets.DOC1(i).name;
end

% get the difference for the export to the mainstem (e.g. attributed to
% physics)
myvals(length(budgets.DOC1)+1)=myvals(11)-nansum(myvals(1:10));

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
myts(length(budgets.DOC1)+1,:)=myts(11,:)-nansum(myts(1:10,:),1);
myts_out.DOC1=myts;

% Below will make a bunch of plots
% make vector of ones to blow up the pie graph
explode1=ones(1,length(myvals(myvals>0)));
explode2=ones(1,length(myvals(myvals<0)));

figure;
pie(myvals(myvals>0),explode1,mynames(myvals>0));

saveas(gcf,'budget_plots/DOC1_budget1','png');
saveas(gcf,'budget_plots/DOC1_budget1','fig');
saveas(gcf,'budget_plots/DOC1_budget1','eps');

figure;
pie(abs(myvals(myvals<0)),explode2,mynames(myvals<0));

saveas(gcf,'budget_plots/DOC1_budget2','png');
saveas(gcf,'budget_plots/DOC1_budget2','fig');
saveas(gcf,'budget_plots/DOC1_budget2','eps');

colormap hsv;
cmap=colormap;
figure;
b=bar(myvals,'FaceColor','flat')
xticklabels(mynames)
xtickangle(45)
b.CData(myvals<0,:)=repmat([1 0 1],length(find(myvals<0)),1);
ylabel('DOC1 (kg C)');

saveas(gcf,'budget_plots/DOC1_budget3','png');
saveas(gcf,'budget_plots/DOC1_budget3','fig');
saveas(gcf,'budget_plots/DOC1_budget3','eps');
[B,Ids]=sort(abs(myvals),'descend');
figure;
color_Vector=1:round(64/length(myvals)):64;
for i = 1 : length(myvals)
     mycolor=cmap(color_Vector(i),:);
     area(time_vector,(myts_out.DOC1(Ids(i),:))','FaceAlpha',0.4,'EdgeAlpha',0.0,'FaceColor',mycolor);
     hold on
     
end
legend(mynames{Ids});
xlabel('Days in 2005');ylabel('DOC1 (g C d^-^1)');

saveas(gcf,'budget_plots/DOC1_budget4','png');
saveas(gcf,'budget_plots/DOC1_budget4','fig');
saveas(gcf,'budget_plots/DOC1_budget4','eps');
% DOC2
%marsh sediment
budgets.DOC2(1).name='Marsh JDOC2';
budgets.DOC2(1).units='kg C';
budgets.DOC2(1).value=marsh.JDOC2.total2;
budgets.DOC2(1).ts=marsh.JDOC2.total1;
budgets.DOC2(1).type='Source';
%estuary sediment
budgets.DOC2(2).name='Estuary JDOC2';
budgets.DOC2(2).units='kg C';
budgets.DOC2(2).value=estuary.JDOC2.total2;
budgets.DOC2(2).ts=estuary.JDOC2.total1;
budgets.DOC2(2).type='Source';
%algae production
budgets.DOC2(3).name='Estuary AlgDOC2';
budgets.DOC2(3).units='kg C';
budgets.DOC2(3).ts=zeros(ntimes-1,1);
budgets.DOC2(3).type='Source';
%hydrolysis of POC --> DOC
budgets.DOC2(4).name='Estuary HdrPOC1';
budgets.DOC2(4).units='kg C';
budgets.DOC2(4).ts=zeros(ntimes-1,1);
budgets.DOC2(4).type='Source';
%remineralization of DOC
budgets.DOC2(5).name='Remineralization of DOC';
budgets.DOC2(5).units='kg C';
budgets.DOC2(5).ts=zeros(ntimes-1,1);
budgets.DOC2(5).type='Sink';
%photochemical remineralization
budgets.DOC2(6).name='Photochemical Remin.';
budgets.DOC2(6).units= 'kg C'
budgets.DOC2(6).ts=zeros(ntimes-1,1);
budgets.DOC2(6).type='Sink';
%riverine input
budgets.DOC2(7).name='Riverine Inputs';
budgets.DOC2(7).units= 'kg C';
budgets.DOC2(7).ts=river.CDOC2.ts+river.NCDOC2.ts;% riverine inputs
budgets.DOC2(7).value=river.CDOC2.int+river.NCDOC2.int;% riverine inputs
budgets.DOC2(7).type='Source';
%Photochem production of DOC2
budgets.DOC2(8).name='Photochem Prod. DOC2';
budgets.DOC2(8).units= 'kg C';
budgets.DOC2(8).ts=zeros(ntimes-1,1);
budgets.DOC2(8).type='Source';
%Photochem loss of DOC2
budgets.DOC2(9).name='Photochem Loss DOC2';
budgets.DOC2(9).units= 'kg C';
budgets.DOC2(9).ts=zeros(ntimes-1,1);
budgets.DOC2(9).type='Sink';
%change in concentration of DOC
budgets.DOC2(10).name='Delta DOC2';
budgets.DOC2(10).units= 'kg C';
budgets.DOC2(10).ts=zeros(ntimes-1,1);
budgets.DOC2(10).type='Source';

% have to add up all the layers
for ilay=1:nlayers;   
    
        %algae
    budgets.DOC2(3).layers(ilay)=nansum(algdoc(ilay,:).*FCD2.*estuary.volume(ilay).all(end,:),2)/1000;
    %hydrolysis of POC --> DOC
    budgets.DOC2(4).layers(ilay)=nansum((hdrlpoc(ilay,:)+hdrrpoc(ilay,:)).*FHDRC2...
                                    .*estuary.volume(ilay).all(end,:),2)/1000;

    budgets.DOC2(5).layers(ilay)=nansum(mnldoc2(ilay,:).*estuary.volume(ilay).all(end,:),2)/1000;
    %photochemical remineralization
    budgets.DOC2(6).layers(ilay)=nansum(dpd20(ilay,:).*estuary.volume(ilay).all(end,:),2)/1000;                                         

    %Photochemical production
    budgets.DOC2(8).layers(ilay)=nansum(( dpd3N(ilay,:).*FPDC32 ...
                                                 +  dpd2N(ilay,:).*FPDC22...
                                                 +  dpd1N(ilay,:).*FPDC12).*estuary.volume(ilay).all(end,:),2)/1000;

    % %        %photochemical bleaching
     budgets.DOC2(9).layers(ilay)=nansum((dpd2N(ilay,:)+dpd21(ilay,:)).*estuary.volume(ilay).all(end,:),2)/1000;        
    
     budgets.DOC2(10).ts=budgets.DOC2(10).ts+nansum(estuary.dtdoc2(ilay).all.*estuary.volume(ilay).all,2)/1000;  
                                       
end

 budgets.DOC2(3).value=nansum(budgets.DOC2(3).layers);
 budgets.DOC2(4).value=nansum(budgets.DOC2(4).layers);
 budgets.DOC2(6).value=nansum(budgets.DOC2(6).layers);
 budgets.DOC2(8).value=nansum(budgets.DOC2(8).layers);
 budgets.DOC2(5).value=nansum(budgets.DOC2(5).layers); 
 budgets.DOC2(9).value=nansum(budgets.DOC2(9).layers);
 budgets.DOC2(10).value=nansum(budgets.DOC2(10).ts.*ICM_int);
 
myvals=[];
mynames={};
myts=[];

 %change the sign and extract the data for plotting
for i = 1 : length(budgets.DOC2);
    if(strcmp(budgets.DOC2(i).type,'Source'));
      myvals(i)=budgets.DOC2(i).value;
      myts(i,:)=budgets.DOC2(i).ts;
    else
      myvals(i)=budgets.DOC2(i).value*-1;
      myts(i,:)=budgets.DOC2(i).ts*-1;
    end
    mynames{i}=budgets.DOC2(i).name;
end

myvals(length(budgets.DOC2)+1)=myvals(10)-nansum(myvals(1:9));

mynames{length(budgets.DOC2)+1}='MainStem'
mypercent=[];
sinks=myvals(myvals(1:end-2)<0);
sources=myvals(myvals(1:end-2)>=0);
mynames_out.DOC2=mynames;

for i = 1: length(myvals);
    mypercent(i)=myvals(i)./nansum(abs(sources));
    mynames{i}=[mynames{i} ' ' num2str(round(mypercent(i)*100,3)) '%'];
end

myvals_out.DOC2=myvals;
myts(length(budgets.DOC2)+1,:)=myts(10,:)-nansum(myts(1:9,:),1);

myts_out.DOC2=myts;

% make vector of ones to blow up the pie graph
explode1=ones(1,length(myvals(myvals>0)));
explode2=ones(1,length(myvals(myvals<0)));

figure;
pie(myvals(myvals>0),explode1,mynames(myvals>0));

saveas(gcf,'budget_plots/DOC2_budget1','png');
saveas(gcf,'budget_plots/DOC2_budget1','fig');
saveas(gcf,'budget_plots/DOC2_budget1','eps');

figure;
pie(abs(myvals(myvals<0)),explode2,mynames(myvals<0));

saveas(gcf,'budget_plots/DOC2_budget2','png');
saveas(gcf,'budget_plots/DOC2_budget2','fig');
saveas(gcf,'budget_plots/DOC2_budget2','eps');

colormap hsv;
cmap=colormap;
figure;
b=bar(myvals,'FaceColor','flat')
xticklabels(mynames)
xtickangle(45)
b.CData(myvals<0,:)=repmat([1 0 1],length(find(myvals<0)),1);
ylabel('DOC2 (kg C)');

saveas(gcf,'budget_plots/DOC2_budget3','png');
saveas(gcf,'budget_plots/DOC2_budget3','fig');
saveas(gcf,'budget_plots/DOC2_budget3','eps');

[B,Ids]=sort(abs(myvals),'descend');
figure;

color_Vector=1:round(64/length(myvals)):64;
for i = 1 : length(myvals)
     mycolor=cmap(color_Vector(i),:);
     area(time_vector,(myts_out.DOC2(Ids(i),:))','FaceAlpha',0.4,'EdgeAlpha',0.0,'FaceColor',mycolor);
     hold on
     
end
legend(mynames{Ids});
xlabel('Days in 2005');ylabel('DOC2 (g C d^-^1)');

saveas(gcf,'budget_plots/DOC2_budget4','png');
saveas(gcf,'budget_plots/DOC2_budget4','fig');
saveas(gcf,'budget_plots/DOC2_budget4','eps');
% DOC3
%marsh sediment
budgets.DOC3(1).name='Marsh JDOC3';
budgets.DOC3(1).units='kg C';
budgets.DOC3(1).value=marsh.JDOC3.total2;
budgets.DOC3(1).ts=marsh.JDOC3.total1;
budgets.DOC3(1).type='Source';
%estuary sediment
budgets.DOC3(2).name='Estuary JDOC3';
budgets.DOC3(2).units='kg C';
budgets.DOC3(2).value=estuary.JDOC3.total2;
budgets.DOC3(2).ts=estuary.JDOC3.total1;
budgets.DOC3(2).type='Source';
%algae production
budgets.DOC3(3).name='Estuary AlgDOC3';
budgets.DOC3(3).units='kg C';
budgets.DOC3(3).ts=zeros(ntimes-1,1);
budgets.DOC3(3).type='Source';
%hydrolysis of POC --> DOC
budgets.DOC3(4).name='Estuary HdrPOC1';
budgets.DOC3(4).units='kg C';
budgets.DOC3(4).ts=zeros(ntimes-1,1);
budgets.DOC3(4).type='Source';
%remineralization of DOC
budgets.DOC3(5).name='Remineralization of DOC';
budgets.DOC3(5).units='kg C';
budgets.DOC3(5).ts=zeros(ntimes-1,1);
budgets.DOC3(5).type='Sink';
%photochemical remineralization
budgets.DOC3(6).name='Photochemical Remin.';
budgets.DOC3(6).units= 'kg C'
budgets.DOC3(6).ts=zeros(ntimes-1,1);
budgets.DOC3(6).type='Sink';
%riverine input
budgets.DOC3(7).name='Riverine Inputs';
budgets.DOC3(7).units= 'kg C';
budgets.DOC3(7).ts=river.CDOC3.ts+river.NCDOC3.ts;% riverine inputs
budgets.DOC3(7).value=river.CDOC3.int+river.NCDOC3.int;% riverine inputs
budgets.DOC3(7).type='Source';
%  %Coagulation of DOC
budgets.DOC3(8).name='Coagulation of DOC3';
budgets.DOC3(8).units='kg C';
budgets.DOC3(8).ts=zeros(ntimes-1,1);
budgets.DOC3(8).type='Sink';
%photochemical loss of DOC3
budgets.DOC3(9).name='Photochem Prod DOC3';
budgets.DOC3(9).units= 'kg C';
budgets.DOC3(9).ts=zeros(ntimes-1,1);
budgets.DOC3(9).type='Source';
%photochemical production of DOC3
budgets.DOC3(10).name='Photochem Loss DOC3';
budgets.DOC3(10).units= 'kg C';
budgets.DOC3(10).ts=zeros(ntimes-1,1);
budgets.DOC3(10).type='Sink';
%change in concentration of DOC
budgets.DOC3(11).name='Delta DOC3';
budgets.DOC3(11).units= 'kg C';
budgets.DOC3(11).ts=zeros(ntimes-1,1);
budgets.DOC3(11).type='Source';

% have to add up all the layers
for ilay=1:nlayers;   
      %algae
    budgets.DOC3(3).layers(ilay)=nansum(algdoc(ilay,:).*FCD3.*estuary.volume(ilay).all(end,:),2)/1000;
    %hydrolysis of POC --> DOC
    budgets.DOC3(4).layers(ilay)=nansum((hdrlpoc(ilay,:)+hdrrpoc(ilay,:)).*FHDRC3...
                                    .*estuary.volume(ilay).all(end,:),2)/1000;

    
    budgets.DOC3(5).layers(ilay)=nansum(mnldoc3(ilay,:).*estuary.volume(ilay).all(end,:),2)/1000;
    %photochemical remineralization
    budgets.DOC3(6).layers(ilay)=nansum(dpd30(ilay,:).*estuary.volume(ilay).all(end,:),2)/1000;                                         
        %coagulation loss of DOC                     
    budgets.DOC3(8).layers(ilay)=nansum(coagc(ilay,:).*estuary.volume(ilay).all(end,:),2)/1000;      
    %Photochemical production
    budgets.DOC3(9).layers(ilay)=nansum(( dpd3N(ilay,:).*FPDC33 ...
                                                 +  dpd2N(ilay,:).*FPDC23...
                                                 +  dpd1N(ilay,:).*FPDC13).*estuary.volume(ilay).all(end,:),2)/1000;
%   photochemical loss  
    budgets.DOC3(10).layers(ilay)=nansum( (dpd31(ilay,:) ...
                                                 +  dpd30(ilay,:)).*estuary.volume(ilay).all(end,:),2)/1000;

   %change in concentration over time
     budgets.DOC3(11).ts=budgets.DOC3(11).ts+nansum(estuary.dtdoc3(ilay).all.*estuary.volume(ilay).all,2)/1000;  

% %                                           
end

 budgets.DOC3(3).value=nansum(budgets.DOC3(3).layers);
 budgets.DOC3(4).value=nansum(budgets.DOC3(4).layers);
 budgets.DOC3(6).value=nansum(budgets.DOC3(6).layers);
 budgets.DOC3(8).value=nansum(budgets.DOC3(8).layers);
 budgets.DOC3(5).value=nansum(budgets.DOC3(5).layers); 
 budgets.DOC3(9).value=nansum(budgets.DOC3(9).layers);
 budgets.DOC3(10).value=nansum(budgets.DOC3(10).layers);
 budgets.DOC3(11).value=nansum(budgets.DOC3(11).ts.*ICM_int);
 
myvals=[];
mynames={};
myts=[];

 %change the sign and extract the data for plotting
for i = 1 : length(budgets.DOC3);
    if(strcmp(budgets.DOC3(i).type,'Source'));
      myvals(i)=budgets.DOC3(i).value;
      myts(i,:)=budgets.DOC3(i).ts;
    else
      myvals(i)=budgets.DOC3(i).value*-1;
      myts(i,:)=budgets.DOC3(i).ts*-1;
    end
    mynames{i}=budgets.DOC3(i).name;
end

 myvals(length(budgets.DOC3)+1)=myvals(11)-nansum(myvals(1:10));

mynames{length(budgets.DOC3)+1}='MainStem'
mypercent=[];
sinks=myvals(myvals(1:end-2)<0);
sources=myvals(myvals(1:end-2)>=0);
mynames_out.DOC3=mynames;

for i = 1: length(myvals);
    mypercent(i)=myvals(i)./nansum(abs(sources));
    mynames{i}=[mynames{i} ' ' num2str(round(mypercent(i)*100,3)) '%'];
end

myvals_out.DOC3=myvals;
myts(length(budgets.DOC3)+1,:)=myts(11,:)-nansum(myts(1:10,:),1);

myts_out.DOC3=myts;

% make vector of ones to blow up the pie graph
explode1=ones(1,length(myvals(myvals>0)));
explode2=ones(1,length(myvals(myvals<0)));

figure;
pie(myvals(myvals>0),explode1,mynames(myvals>0));

saveas(gcf,'budget_plots/DOC3_budget1','png');
saveas(gcf,'budget_plots/DOC3_budget1','fig');
saveas(gcf,'budget_plots/DOC3_budget1','eps');

figure;
pie(abs(myvals(myvals<0)),explode2,mynames(myvals<0));

saveas(gcf,'budget_plots/DOC3_budget2','png');
saveas(gcf,'budget_plots/DOC3_budget2','fig');
saveas(gcf,'budget_plots/DOC3_budget2','eps');

colormap hsv;
cmap=colormap;
figure;
b=bar(myvals,'FaceColor','flat')
xticklabels(mynames)
xtickangle(45)
b.CData(myvals<0,:)=repmat([1 0 1],length(find(myvals<0)),1);
ylabel('DOC3 (kg C)');

saveas(gcf,'budget_plots/DOC3_budget3','png');
saveas(gcf,'budget_plots/DOC3_budget3','fig');
saveas(gcf,'budget_plots/DOC3_budget3','eps');

[B,Ids]=sort(abs(myvals),'descend');
figure;

color_Vector=1:round(64/length(myvals)):64;
for i = 1 : length(myvals)
     mycolor=cmap(color_Vector(i),:);
     area(time_vector,(myts_out.DOC3(Ids(i),:))','FaceAlpha',0.4,'EdgeAlpha',0.0,'FaceColor',mycolor);
     hold on
     
end
legend(mynames{Ids});
xlabel('Days in 2005');ylabel('DOC3 (g C d^-^1)');

saveas(gcf,'budget_plots/DOC3_budget4','png');
saveas(gcf,'budget_plots/DOC3_budget4','fig');
saveas(gcf,'budget_plots/DOC3_budget4','eps');
%
% DOC
%marsh sediment
budgets.TDOC(1).name='Marsh JDOC';
budgets.TDOC(1).units='kg C';
budgets.TDOC(1).value=marsh.JDOC1.total2+marsh.JDOC2.total2+marsh.JDOC3.total2;
budgets.TDOC(1).ts=marsh.JDOC1.total1+marsh.JDOC2.total1+marsh.JDOC3.total1;
budgets.TDOC(1).type='Source';
%estuary sediment
budgets.TDOC(2).name='Estuary JDOC';
budgets.TDOC(2).units='kg C';
budgets.TDOC(2).value=estuary.JDOC1.total2+estuary.JDOC2.total2+estuary.JDOC3.total2;
budgets.TDOC(2).ts=estuary.JDOC1.total1+estuary.JDOC2.total1+estuary.JDOC3.total1;
budgets.TDOC(2).type='Source';
%algae production
budgets.TDOC(3).name='Estuary AlgDOC';
budgets.TDOC(3).units='kg C';
budgets.TDOC(3).value=budgets.DOC1(3).value+budgets.DOC2(3).value+budgets.DOC3(3).value;
budgets.TDOC(3).ts=+budgets.DOC1(3).ts+budgets.DOC2(3).ts+budgets.DOC3(3).ts;
budgets.TDOC(3).type='Source';
%hydrolysis of POC --> DOC
budgets.TDOC(4).name='Estuary HdrPOC';
budgets.TDOC(4).units='kg C';
budgets.TDOC(4).value=budgets.DOC1(4).value+budgets.DOC2(4).value+budgets.DOC3(4).value;
budgets.TDOC(4).ts=budgets.DOC1(4).ts+budgets.DOC2(4).ts+budgets.DOC3(4).ts;
budgets.TDOC(4).type='Source';
%denitrification of DOC
budgets.TDOC(5).name='Denitrification of DOC';
budgets.TDOC(5).units='kg C';
budgets.TDOC(5).value=budgets.DOC1(5).value;
budgets.TDOC(5).ts=budgets.DOC1(5).ts;
budgets.TDOC(5).type='Sink';
%remineralization of DOC
budgets.TDOC(6).name='Remineralization of DOC';
budgets.TDOC(6).units='kg C';
budgets.TDOC(6).value=budgets.DOC1(6).value+budgets.DOC2(5).value+budgets.DOC3(5).value;
budgets.TDOC(6).ts=budgets.DOC1(6).ts+budgets.DOC2(5).ts+budgets.DOC3(5).ts;
budgets.TDOC(6).type='Sink';
%photochemical remineralization
budgets.TDOC(7).name='Photochemical Remin.';
budgets.TDOC(7).units= 'kg C'
budgets.TDOC(7).value=budgets.DOC1(7).value+budgets.DOC2(6).value+budgets.DOC3(6).value;
budgets.TDOC(7).ts=budgets.DOC1(7).ts+budgets.DOC2(6).ts+budgets.DOC3(6).ts;
budgets.TDOC(7).type='Sink';
%riverine input
budgets.TDOC(8).name='Riverine Inputs';
budgets.TDOC(8).units= 'kg C';
budgets.TDOC(8).value=budgets.DOC1(8).value+budgets.DOC2(7).value+budgets.DOC3(7).value;
budgets.TDOC(8).ts=budgets.DOC1(8).ts+budgets.DOC2(7).ts+budgets.DOC3(7).ts;
budgets.TDOC(8).type='Source';
 %Coagulation of DOC
budgets.TDOC(9).name='Coagulation of DOC';
budgets.TDOC(9).value=budgets.DOC3(8).value;
budgets.TDOC(9).ts=budgets.DOC3(8).ts;
budgets.TDOC(9).type='Sink';
%change in concentration of DOC
budgets.TDOC(10).name='Delta DOC';
budgets.TDOC(10).units= 'kg C';
budgets.TDOC(10).value=budgets.DOC1(11).value+budgets.DOC2(10).value+budgets.DOC3(11).value;
budgets.TDOC(10).ts=budgets.DOC1(11).ts+budgets.DOC2(10).ts+budgets.DOC3(11).ts;
budgets.TDOC(10).type='Source';

myvals=[];
mynames={};
myts=[];

 %change the sign and extract the data for plotting
for i = 1 : length(budgets.TDOC);
    if(strcmp(budgets.TDOC(i).type,'Source'));
      myvals(i)=budgets.TDOC(i).value;
      myts(i,:)=budgets.TDOC(i).ts;
    else
      myvals(i)=budgets.TDOC(i).value*-1;
      myts(i,:)=budgets.TDOC(i).ts*-1;
    end
    mynames{i}=budgets.TDOC(i).name;
end

myvals(length(budgets.TDOC)+1)=myvals(10)-nansum(myvals(1:9));

mynames{length(budgets.TDOC)+1}='MainStem'
mypercent=[];
sinks=myvals(myvals(1:end-2)<0);
sources=myvals(myvals(1:end-2)>=0);
mynames_out.DOC=mynames;

for i = 1: length(myvals);
    mypercent(i)=myvals(i)./nansum(abs(sources));
    mynames{i}=[mynames{i} ' ' num2str(round(mypercent(i)*100,3)) '%'];
end

myvals_out.DOC=myvals;
myts(length(budgets.TDOC)+1,:)=myts(10,:)-nansum(myts(1:9,:),1);

myts_out.DOC=myts;

% make vector of ones to blow up the pie graph
explode1=ones(1,length(myvals(myvals>0)));
explode2=ones(1,length(myvals(myvals<0)));

figure;
pie(myvals(myvals>0),explode1,mynames(myvals>0));

saveas(gcf,'budget_plots/DOC_budget1','png');
saveas(gcf,'budget_plots/DOC_budget1','fig');
saveas(gcf,'budget_plots/DOC_budget1','eps');

figure;
pie(abs(myvals(myvals<0)),explode2,mynames(myvals<0));

saveas(gcf,'budget_plots/DOC_budget2','png');
saveas(gcf,'budget_plots/DOC_budget2','fig');
saveas(gcf,'budget_plots/DOC_budget2','eps');

colormap hsv;
cmap=colormap;
figure;
b=bar(myvals,'FaceColor','flat')
xticklabels(mynames)
xtickangle(45)
b.CData(myvals<0,:)=repmat([1 0 1],length(find(myvals<0)),1);
ylabel('DOC (kg C)');

saveas(gcf,'budget_plots/DOC_budget3','png');
saveas(gcf,'budget_plots/DOC_budget3','fig');
saveas(gcf,'budget_plots/DOC_budget3','eps');

[B,Ids]=sort(abs(myvals),'descend');
figure;

color_Vector=1:round(64/length(myvals)):64;
for i = 1 : length(myvals)
     mycolor=cmap(color_Vector(i),:);
     area(time_vector,(myts_out.DOC(Ids(i),:))','FaceAlpha',0.4,'EdgeAlpha',0.0,'FaceColor',mycolor);
     hold on
     
end
legend(mynames{Ids});
xlabel('Days in 2005');ylabel('DOC (g C d^-^1)');

saveas(gcf,'budget_plots/DOC_budget4','png');
saveas(gcf,'budget_plots/DOC_budget4','fig');
saveas(gcf,'budget_plots/DOC_budget4','eps');
%% POC
budgets.POC(1).name='marsh JPOC';
budgets.POC(1).units='kg C';
budgets.POC(1).value=marsh.JPOC.total2;
budgets.POC(1).ts=marsh.JPOC.total1;
budgets.POC(1).type='Sink';
%estuary sediment
budgets.POC(2).name='Estuary JPOC';
budgets.POC(2).units='kg C';
budgets.POC(2).value=estuary.JPOC.total2;
budgets.POC(2).ts=estuary.JPOC.total1;
budgets.POC(2).type='Sink';
%algae production
budgets.POC(3).name='Estuary AlgPOC';
budgets.POC(3).units='kg C';
budgets.POC(3).ts=zeros(ntimes-1,1);
budgets.POC(3).type='Source';
%hydrolysis of POC --> DOC
budgets.POC(4).name='Estuary HdrPOC';
budgets.POC(4).units='kg C';
budgets.POC(4).ts=zeros(ntimes-1,1);
budgets.POC(4).type='Sink';
%riverine input
budgets.POC(5).name='Riverine Inputs';
budgets.POC(5).units= 'kg C';
budgets.POC(5).ts=river.LPOC.ts+river.RPOC.ts;% riverine inputs
budgets.POC(5).value=river.LPOC.int+river.RPOC.int;% riverine inputs
budgets.POC(5).type='Source';
%coagulation of DOC
budgets.POC(6).name='Coagulation DOC'
budgets.POC(6).units = 'kg C'
budgets.POC(6).ts=budgets.TDOC(9).ts;
budgets.POC(6).value=budgets.TDOC(9).value
budgets.POC(6).type='Source'
%Delta POC
budgets.POC(7).name='Delta POC'
budgets.POC(7).units = 'kg C'
budgets.POC(7).ts=zeros(ntimes-1,1)
budgets.POC(7).type='Source'

% have to add up all the layers
for ilay=1:nlayers;    
    %algae
    budgets.POC(3).layers(ilay)=nansum(algpoc(ilay,:).*estuary.volume(ilay).all(end,:),2)/1000;
    %hydrolysis of POC --> DOC
    budgets.POC(4).layers(ilay)=nansum((hdrlpoc(ilay,:)+hdrrpoc(ilay,:))...
                                    .*estuary.volume(ilay).all(end,:),2)/1000;
  budgets.POC(7).ts=budgets.POC(7).ts+nansum((estuary.dtlpoc(ilay).all+estuary.dtrpoc(ilay).all).*estuary.volume(ilay).all,2)/1000;                               
end

 budgets.POC(3).value=nansum(budgets.POC(3).layers);
 budgets.POC(4).value=nansum(budgets.POC(4).layers);
 budgets.POC(7).value=nansum(budgets.POC(7).ts.*ICM_int);
 
myvals=[];
mynames={};
myts=[];
 %change the sign and extract the data for plotting
for i = 1 : length(budgets.POC);
    if(strcmp(budgets.POC(i).type,'Source'));
      myvals(i)=budgets.POC(i).value;
      myts(i,:)=budgets.POC(i).ts;
    else
      myvals(i)=budgets.POC(i).value*-1;
      myts(i,:)=budgets.POC(i).ts*-1;
    end
    hold on
    mynames{i}=budgets.POC(i).name;
end

myvals(length(budgets.POC)+1)=myvals(7)-nansum(myvals(1:6));
% mynames{length(budgets.POC)+1}=['Main Stem ' num2str(round(myvals(end)./nansum(abs(sources))*100)) '%'];
mynames{length(budgets.POC)+1}='MainStem'
mypercent=[];
sinks=myvals(myvals(1:end-2)<0);
sources=myvals(myvals(1:end-2)>=0);
mynames_out.POC=mynames;

for i = 1: length(myvals);
    mypercent(i)=myvals(i)./nansum(abs(sources))
    mynames{i}=[mynames{i} ' ' num2str(round(mypercent(i)*100,3)) '%']
end

myvals_out.POC=myvals;
myts(length(budgets.POC)+1,:)=myts(7,:)-nansum(myts(1:6,:),1);
myts_out.POC=myts;
% make vector of ones to blow up the pie graph
explode1=ones(1,length(myvals(myvals>0)));
explode2=ones(1,length(myvals(myvals<0)));
figure;
pie(myvals(myvals>0),explode1,mynames(myvals>0));
saveas(gcf,'budget_plots/POC_budget1','png');
saveas(gcf,'budget_plots/POC_budget1','fig');
saveas(gcf,'budget_plots/POC_budget1','eps');
figure;
pie(abs(myvals(myvals<0)),explode2,mynames(myvals<0));
saveas(gcf,'budget_plots/POC_budget2','png');
saveas(gcf,'budget_plots/POC_budget2','fig');
saveas(gcf,'budget_plots/POC_budget2','eps');
figure;
b=bar(myvals,'FaceColor','flat')
xticklabels(mynames)
xtickangle(45)
b.CData(myvals<0,:)=repmat([1 0 1],length(find(myvals<0)),1);
ylabel('POC (kg C)');
saveas(gcf,'budget_plots/POC_budget3','png');
saveas(gcf,'budget_plots/POC_budget3','fig');
saveas(gcf,'budget_plots/POC_budget3','eps');
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
xlabel('Days in 2005');ylabel('POC (g C m^-^3 d^-^1)');
saveas(gcf,'budget_plots/POC_budget4','png');
saveas(gcf,'budget_plots/POC_budget4','fig');
saveas(gcf,'budget_plots/POC_budget4','eps');

save('Carbon_budget_out.mat','myvals_out','mynames_out','myts_out')















