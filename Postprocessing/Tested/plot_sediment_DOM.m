% script to plot the sediment DOM flux at various statinos
%make a *.mat file with a list of node numbers to plot at

% B Clark February, 2019
%UMCES, HPL



matfilein = input('What is the full path and file name to the *.mat grid file? ---> ','s');
load(matfilein);
wqm_hisdir = input('What is full path to the WQM history directory? --->  ','s');
nlayers = input('How many layers are in the model? --> ');
outdir = 'nutrient_plots';
ICM_int = input('What is the interval that the data was output at, in days?  ---> ')
% ncfilein = input('What is the full path and file name of the FVCOM nc file used to run the WQM? ---> ','s');
% nc=netcdf(ncfilein);
load(input('What is the full path and file name to the sediment stations .mat? ---> ','s'));
% load('/data/users/bclark/CBAY_DATA/Sediment_stations.mat');

% CBP_struct = input('What is the full path and file name to the .mat file with the CBP data structure ---> ','s');
% load(CBP_struct);
first_day = input('What is the first day of the model? --->  ');

domflag = input('Is this the WQM with the new DOM or old DOM formulas? enter 0 for no and 1 for yes ---->')

numfile=input('How many files do you want to compare? ---> ');
nameseg=input('What is the prefix for the ICM files? ---> ','s');

%coordinate system

fname = [wqm_hisdir '/' nameseg,'_0001.nc'];
nc=netcdf(fname);
s_hybrid=nc{'siglev'}(:);
dz=diff(s_hybrid,1)*-1;


JWNCDOC1_mod=[];
JWNCDOC2_mod=[];
JWNCDOC3_mod=[];
JWCDOC1_mod=[];
JWCDOC2_mod=[];
JWCDOC3_mod=[];
depth_mod=[];

mkdir('Sediment_Plots');
for idir=1 : numfile;   
    numfile_count = [1:numfile];   
    kstr = num2str(numfile_count(idir),'%04d');    
    fname = [wqm_hisdir '/' nameseg,'_',kstr,'.nc'];
    nc = netcdf([fname]);
    idir
 
    var_in=nc{'JWCDOC1'}(:);
    JWCDOC1_mod=[JWCDOC1_mod; var_in(:,unique_sedstations)];
    var_in=nc{'JWCDOC2'}(:);
    JWCDOC2_mod=[JWCDOC2_mod; var_in(:,unique_sedstations)];
    var_in=nc{'JWCDOC3'}(:);
    JWCDOC3_mod=[JWCDOC3_mod; var_in(:,unique_sedstations)];

    var_in=nc{'JWNCDOC1'}(:);
    JWNCDOC1_mod=[JWNCDOC1_mod; var_in(:,unique_sedstations)];
    var_in=nc{'JWNCDOC2'}(:);
    JWNCDOC2_mod=[JWNCDOC2_mod; var_in(:,unique_sedstations)];
    var_in=nc{'JWNCDOC3'}(:);
    JWNCDOC3_mod=[JWNCDOC3_mod; var_in(:,unique_sedstations)];

end

for istat=1:length(unique_sedstations);
    mystat=unique_sedstations(istat);
    subplot(4,2,1)
    tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3),-h_n);
    hold on
    plot(lld_n(mystat,2),lld_n(mystat,3),'kd');
    subplot(4,2,2);
    plot(JWCDOC1_mod(:,istat)+JWCDOC2_mod(:,istat)+JWCDOC3_mod(:,istat)+...
        JWNCDOC1_mod(:,istat)+JWNCDOC2_mod(:,istat)+JWNCDOC3_mod(:,istat));
    xlabel('Time');ylabel('tDOC (g C m^-^2 d^-^1');
    subplot(4,2,3);
    plot(JWCDOC1_mod(:,istat))
     xlabel('Time');ylabel('JCDOC1 (g C m^-^2 d^-^1)');
    subplot(4,2,4);
    plot(JWCDOC2_mod(:,istat))
     xlabel('Time');ylabel('JCDOC2 (g C m^-^2 d^-^1)'); 
   subplot(4,2,5);
    plot(JWCDOC3_mod(:,istat))
     xlabel('Time');ylabel('JCDOC3 (g C m^-^2 d^-^1)');
    subplot(4,2,6);
    plot(JWNCDOC1_mod(:,istat))
     xlabel('Time');ylabel('JNCDOC1 (g C m^-^2 d^-^1)');
    subplot(4,2,7);
    plot(JWNCDOC2_mod(:,istat))
     xlabel('Time');ylabel('JNCDOC2 (g C m^-^2 d^-^1)');     
    subplot(4,2,8);
    plot(JWNCDOC3_mod(:,istat))
     xlabel('Time');ylabel('JNCDOC3 (g C m^-^2 d^-^1)');     
     
     
   saveas(gcf,['Sediment_Plots/DOM_station_' num2str(unique_sedstations(istat)) '.png']);
   saveas(gcf,['Sediment_Plots/DOM_station_' num2str(unique_sedstations(istat)) '.fig']);  
   close all
end

disp('Finished');































