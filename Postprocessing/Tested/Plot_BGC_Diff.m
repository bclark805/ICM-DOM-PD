%script to pull out and plot the average difference between

% the +M and -M scenarios for Ammonium (NH4)
%Nitrate (NO3), NPP, and PAR

%this is ideally plotted using tecplot, 

%but here will just be contoured using tricolor

% B Clark, UMCES, HPL, January 2019


nc_RR=netcdf(input('What is the full path and file name to the +M netcdf output file? -->','s'));
nc_NM=netcdf(input('What is the full path and file name to the -M netcdf output file? -->','s'));

% load(input('What is the full path and file name to the *.mat grid file? ---> ','s');
load('/Users/bclark/Desktop/RhodeFVCOM_2005/Data/rhoderiver_grid_v1_2.mat');

tec_flag=input('Do you want to make a tecplot output file as well? 0 for no and 1 for yes --->');
%pull out the data
disp('Getting NPP')
NPP_RR=nc_RR{'NPP'}(:);
NPP_NM=nc_NM{'NPP'}(:);

disp('Getting NH4')
NH4_RR=nc_RR{'NH4'}(:);
NH4_NM=nc_NM{'NH4'}(:);

disp('Getting NO3')
NO3_RR=nc_RR{'NO3'}(:);
NO3_NM=nc_NM{'NO3'}(:);

disp('Getting PAR')
PAR_RR=nc_RR{'PAR'}(:);
PAR_NM=nc_NM{'PAR'}(:);
%%
%get the average difference for each variable
mean_NPP_diff=(mean(NPP_RR-NPP_NM))';
mean_NH4_diff=mean(squeeze(NH4_RR(:,5,:))-squeeze(NH4_NM(:,5,:)),1)';
mean_NO3_diff=mean(squeeze(NO3_RR(:,5,:))-squeeze(NO3_NM(:,5,:)),1)';
mean_PAR_diff=mean(squeeze(PAR_RR(:,5,:))-squeeze(PAR_NM(:,5,:)),1)';

%plot a 4 panel for the RHode River
%convert to mg m^-3
subplot(2,2,1);
tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3),mean_NH4_diff.*1000,0,60);
axis([-76.565 -76.492 38.853 38.905]);

subplot(2,2,2);
tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3),mean_NO3_diff.*1000,-1,12);
axis([-76.565 -76.492 38.853 38.905]);

subplot(2,2,3);
tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3),mean_NPP_diff.*1000,-150,150);
axis([-76.565 -76.492 38.853 38.905]);

subplot(2,2,4);
tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3),mean_PAR_diff,-2.4,0);
axis([-76.565 -76.492 38.853 38.905]);

saveas(gcf,'BGC_Diffed_4pane.eps')
saveas(gcf,'BGC_Diffed_4pane.fig')
%set up a tecplot data structure
%this works wtih Tecplot 360
if tec_flag
    tdata.varnames={'x','y','z','NH4','NO3','NPP','PAR'};
    tdata.FEsurfaces(1).x=xyd_n(:,2);  %2nd column in xyd_n is x
    tdata.FEsurfaces(1).y=xyd_n(:,3);  %3rd                    y
    tdata.FEsurfaces(1).order=3;       %surface defined on (x,y) coord
    tdata.FEsurfaces(1).e2n=e2n(:,2:4);%each row(element) has 3 node numbers
    
    tdata.FEsurfaces(1).v(1,:)=mean_NH4_diff.*1000;
    tdata.FEsurfaces(1).v(2,:)=mean_NO3_diff.*1000;
    tdata.FEsurfaces(1).v(3,:)=mean_NPP_diff.*1000;
    tdata.FEsurfaces(1).v(4,:)=mean_PAR_diff;
    
    mat2tecplot(tdata,'BGC_diffed_4panel.plt');
    
end
    
    







