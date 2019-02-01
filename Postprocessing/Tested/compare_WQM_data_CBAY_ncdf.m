% script to plot the WQM output at model nodes
%closest to the stations where Chesapeake Bay Program Data
%has been sampled
%then make a comparison to the data and run some basic statistics
%
%The data structure to compare to with is created by downloading and
%loading in Chesapeake Bay Program data, and then running the Matlab script
%
% Read_CBP.m to sort the data into a nice, tidy data structure.

% B Clark sep 2015
%UMCES, HPL
%updated B Clark UMCES, January 2019
%requres a few inputs but also for the
%ICM netcdf output to be in one BIG file (can be easily modified to use
%many files, but this is very slow)
%
%
%
%=========================================================================%
%
%
%
matfilein = input('What is the full path and file name to the *.mat grid file? ---> ','s');
load(matfilein);
wqm_hisdir = input('What is full path to the WQM history directory? --->  ','s');
nlayers = input('How many layers are in the model? --> ');

ICM_int = input('What is the interval that the data was output at, in days?  ---> ')

%load(*/RhodeFVCOM_2005/Data/
CBP_struct = input('What is the full path and file name to the .mat file with the CBP data structure ---> ','s');
load(CBP_struct);
first_day = input('What is the first day of the model? --->  ');

numfile=input('How many files do you want to compare? ---> ');
nameseg=input('What is the prefix for the ICM files? ---> ','s');

clear nc
%
outdir = 'nutrient_plots';
if(~exist(outdir));
    mkdir(outdir);
end

% get the station indices by matching the CBP data coordinates
% with the FVCOM grid coordinates
stations_mod = dsearchn([lld_n(:,2) lld_n(:,3)],([stations.lon ; stations.lat])'); %find the model nodes that are close to the wqm station points

tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3), h_n);
hold on
for i = 1 : length(stations_mod);
    text(lld_n(stations_mod(i),2),lld_n(stations_mod(i),3),stations(i).name,'fontsize',24)
    
end

saveas(gcf,'AllStations_map','png');

alg1_mod=[];
alg2_mod=[];
DO_mod=[];
NH4_mod=[];
NO3_mod=[];
depth1=[];
DON_mod=[];
DOC_mod=[];
POC_mod=[];
PON_mod=[];
SALT_mod=[];
TEMP_mod=[];
KD_mod=[];
mod_time=[];
numfile_count=1;
idir=1;
kstr = num2str(numfile_count(idir),'%04d');
fname = [wqm_hisdir '/' nameseg,'_',kstr,'.nc'];

%pull out the vertical coordinate distribution
%this is the FVCOM hydrodynamics file, used to load in the vertical
%coordinate system
% ncfile=input('What is the full path and file name to the first netcdf hydrodynamics file? --->  ','s');
nc=netcdf(fname);
s_hybrid=nc{'siglev'}(:);
% dz=diff(s_hybrid,1);
%get the fraction of the water column that each layer is at
kb=11;
for k = 1 : kb;
s_hybrid(k,1:length(xyd_n))=-((k-1)/(kb-1))^1.5;
end
for i = 2:11;
layer_depth(i-1,:)=mean(s_hybrid(i-1:i,:),1)*-1;
end
% for i = 2:nlayers+1;
%     layer_depth(i,:)=mean(s_hybrid(i-1:i,:),1);
% end
%%
idir
% get out all of the data, can't just do stations because must be at a
% constant interval
disp('Getting all of the Water Quality Variables');
disp('Getting Algae, Oxygen and Inorganic Nitrogen');
alg1_in= nc{'B1'}(:);
alg2_in  = nc{'B2'}(:);
DO_in  = nc{'DOXG'}(:);
NH4_in  = nc{'NH4'}(:);
NO3_in = nc{'NO3'}(:);

depth_in= nc{'depth'}(:);
disp('Getting DON');
DON_in = nc{'WC_CDON1'}(:) + ...
    nc{'WC_NCDON1'}(:) + ...
    nc{'WC_CDON2'}(:) + ...
    nc{'WC_NCDON2'}(:) + ...
    nc{'WC_CDON3'}(:) + ...
    nc{'WC_NCDON3'}(:);
disp('Getting DOC');
DOC_in = nc{'WC_CDOC1'}(:) + ...
    nc{'WC_NCDOC1'}(:) + ...
    nc{'WC_CDOC2'}(:) + ...
    nc{'WC_NCDOC2'}(:) + ...
    nc{'WC_CDOC3'}(:) + ...
    nc{'WC_NCDOC3'}(:);

disp('Getting POM');
POC_in= nc{'LPOC'}(:) + nc{'RPOC'}(:);
PON_in= nc{'LPON'}(:) + nc{'RPON'}(:);


SALT_in = nc{'salinity'}(:);
TEMP_in  = nc{'temp'}(:);
KD_in = nc{'KD'}(:);

mod_time=[mod_time;nc{'time'}(:)./86400];
%now scale to the stations;
alg1_mod= alg1_in(:,:,stations_mod);
alg2_mod= alg2_in(:,:,stations_mod);
DO_mod=DO_in(:,:,stations_mod);
NH4_mod= NH4_in(:,:,stations_mod);
NO3_mod= NO3_in(:,:,stations_mod);

depth1=depth_in(:,stations_mod);
DON_mod= DON_in(:,:,stations_mod);
DOC_mod= DOC_in(:,:,stations_mod);

PON_mod= PON_in(:,:,stations_mod)+(alg1_mod+alg2_mod)./5.68;% Use Redfield to convert algae C to N
POC_mod= POC_in(:,:,stations_mod)+alg1_mod+alg2_mod;

SALT_mod= SALT_in(:,:,stations_mod);
TEMP_mod= TEMP_in(:,:,stations_mod);

KD_mod=[KD_mod ; KD_in(:,:,stations_mod)];
clear nc

starting_date=first_day; data_int=ICM_int;

chla_mod = (alg1_mod./50+alg2_mod./50).*1000  ;% convert mg C/l to ug chl/l based on chl:C ratio from model

mod_time=mod_time+starting_date;

[ntimes,nlayers,nstations]=size(POC_mod);

%This loop builds a grid for the contour plots y dimension
for iz=1:nlayers;
    depth(:,iz,:)=depth1.*repmat(layer_depth(iz,stations_mod),ntimes,1);
end

%
rounded_mod_times=round(mod_time);

% [nlayers,ntimes,nstations]=size(chla_mod);
%%
% Now loop over all the stations and make the comparisons
for istat=1:length(stations_mod);
    istat
    
    statstring = ['station' num2str(istat)];
    mkdir(statstring)
    figure
    tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3),-h_n);
    hold on
    plot(lld_n(stations_mod(istat),2),lld_n(stations_mod(istat),3),'*','color','k','markersize',24);
    
    saveas(gcf,[statstring '/station_map'],'png');

    myvar_id=[];
    ivar=0;
    while isempty(myvar_id)
        ivar=ivar+1;
        if(strcmp(stations(istat).variables(ivar).name,'CHLA'))
            myvar_id=ivar              ; % change this ID based on where the variable is located within the CBP data structure
            % generated by read_cbp.m
        end
    end
    
    % pull out the observational data for this station
    [unique_dates,unique_ids]=unique(stations(istat).variables(myvar_id).date);
    cbp_station_date= stations(istat).variables(myvar_id).date ;
    cbp_station_value= stations(istat).variables(myvar_id).value;
    cbp_station_depths=stations(istat).variables(myvar_id).depth;
    
    % now for each unique date, find corresponding values and take the
    % get the matching model times for statistics and plotting;
    [matching_mod_times,ia,ib]=intersect(rounded_mod_times,unique_dates);
    % now plot date by date
    figure;
    
    X_out=[];
    Y_out=[];
    X_out_mean=[];
    Y_out_mean=[];
    unique_dates_out=[];
    
    subplot(2,2,[1,2]);
    
 %This contours the model output for this station   
    contourf(repmat(mod_time,1,nlayers),depth(:,:,istat)*-1,chla_mod(:,:,istat),'LineStyle','none')
    hold on
    [cmin,cmax]=caxis;
    xlabel('Days in 2005');ylabel('Depth (m)')
    
    for iday=1:length(ib)
        % get the good day from the data back out so we can get
        % each value and depth
        myIds=find(cbp_station_date==unique_dates(ib(iday))); 
        myDates=cbp_station_date(myIds);
        myDepths=cbp_station_depths(myIds);
        myvalues=cbp_station_value(myIds);
        
        % get rid of duplicates;
        [myUniqueDepths,Zids]=unique(myDepths);
        myUniqueValues=myvalues(Zids);
        myUniqueDates=unique(myDates);
        
        % now match up the model
        model_depthIds=dsearchn(depth(ia(iday),:,istat)',myUniqueDepths)  ;
        modelDepths=depth(ia(iday),model_depthIds,istat);
        modelValues=chla_mod(ia(iday),model_depthIds,istat);
        
        % now plot a scatter of the observations over the model
        for iplot=1:length(modelDepths);
            scatter(myUniqueDates,modelDepths(iplot)*-1,100,myUniqueValues(iplot),'filled','markeredgecolor','k');
        end
        %get the values for each depth and each time in one
        %long vector to compare directly to eachother
        X_out=[X_out;myUniqueValues];
        Y_out=[Y_out;modelValues'];
        
        % take the depth averaged quantity
        X_out_mean=[X_out_mean;mean(myUniqueValues)];
        Y_out_mean=[Y_out_mean;mean(modelValues)];
        
        unique_dates_out=[unique_dates_out;myUniqueDates];
        
    end
    h= colorbar;
    h.Label.String='chla (ug chla l^-^1)'   ;
    %Get the statistics
    [MEF,r,WMS,RMSE,MPE,RI,axis_min,axis_max]=model_statistics(X_out,Y_out)   ;
    caxis([axis_min axis_max]);
    title([ 'Chl a ' , '  MEF =' num2str(MEF),', r = ',num2str(r) ', RMSE = ' num2str(RMSE) , ', MPE = ',num2str(MPE)  ', WMS = ' num2str(WMS),', RI = ',num2str(RI)]);
    subplot(2,2,3);
    
    plot(X_out,Y_out,'kd','markersize',8);
    axis([axis_min axis_max axis_min axis_max]);
    
    xlabel('Observed Chl a (ug l^-^1)');ylabel('Modeled Chl a (ug l^-^1)')  ;
    xlabel('Days in 2005');ylabel('Depth (m)')  ;
    
    subplot(2,2,4);
    %     histogram(Y_out,10,'normalization','probability');
    %     hold on
    %     histogram(X_out,10,'normalization','probability');
    g1={};g2={};g3={};
    VAR_cat=reshape(chla_mod(:,:,istat),1,nlayers*ntimes);
    g1(1:length(X_out))={'Observed'};% get cell arrays of our group names
    g2(1:length(Y_out))={'Modeled'};
    g3(1:length(VAR_cat))={'Modeled-ALL'};
    grouping=[g1 g2 g3];
    boxplot([X_out ; Y_out; VAR_cat'],grouping);
    
    ylabel('Chl a (ug l^-^1)');%ylabel('Probability')
    saveas(gcf,[statstring '/chla_depth'],'png');
    saveas(gcf,[statstring '/chla_depth'],'eps');
    saveas(gcf,[statstring '/chla_depth'],'fig');

    figure;     % depth averaged
    subplot(2,2,[1,2]);
    [MEF,r,WMS,RMSE,MPE,RI,axis_min,axis_max]=model_statistics(X_out_mean,Y_out_mean)   ;
    
    plot(mod_time,mean(chla_mod(:,:,istat),2),'k-o');
    hold on
    plot(unique_dates_out,X_out_mean,'r-*');
    title([ 'Depth Averaged Chl a ' ,'  MEF =' num2str(MEF),', r = ',num2str(r) ', RMSE = ' num2str(RMSE) , ', MPE = ',num2str(MPE)  ', WMS = ' num2str(WMS),', RI = ',num2str(RI)]);
    xlabel('Days in 2005');ylabel('Chl a (ug l^-^1)');
    legend('Modeled','Observed');
    subplot(2,2,3);
    plot(X_out_mean,Y_out_mean,'rd','markersize',10);
    axis([axis_min axis_max axis_min axis_max]);
    xlabel('Observed Chl a (ug l^-^1)');ylabel('Modeled Chl a (ug l^-^1)')  ;
    
    subplot(2,2,4);
    %     histogram(Y_out,10,'normalization','probability');
    %     hold on
    %     histogram(X_out,10,'normalization','probability');
    g1={};g2={};g3={};
    VAR_cat=reshape(mean(chla_mod(:,:,istat),2),1,ntimes);
    g1(1:length(X_out_mean))={'Observed'};% get cell arrays of our group names
    g2(1:length(Y_out_mean))={'Modeled'};
    g3(1:length(VAR_cat))={'Modeled-ALL'};
    grouping=[g1 g2 g3];
    boxplot([X_out_mean ; Y_out_mean; VAR_cat'],grouping);
    
    ylabel('Chl a (ug l^-^1)');%ylabel('Probability')
    %     xlabel('Chl a (ug l^-^1)');ylabel('Probability')
    
    saveas(gcf,[statstring '/chla_depth_stats'],'png');
    saveas(gcf,[statstring '/chla_depth_stats'],'eps');
    saveas(gcf,[statstring '/chla_depth_stats'],'fig');
    %
    %========================= Dissolved Oxygen ============================
    
    myvar_id=[];
    ivar=0;
    while isempty(myvar_id);
        ivar=ivar+1;
        if(strcmp(stations(istat).variables(ivar).name,'DO'));
            myvar_id=ivar              ; % change this ID based on where the variable is located within the CBP data structure
            % generated by read_cbp.m
        end
    end
    % pull out the data for this station
    [unique_dates,unique_ids]=unique(stations(istat).variables(myvar_id).date);
    cbp_station_date= stations(istat).variables(myvar_id).date ;
    cbp_station_value= stations(istat).variables(myvar_id).value;
    cbp_station_depths=stations(istat).variables(myvar_id).depth;
    % now for each unique date, find corresponding values and take the
    
    % get the matching model times for statistics and plotting;
    [matching_mod_times,ia,ib]=intersect(rounded_mod_times,unique_dates);
    
    % now plot date by date
    figure;
    
    X_out=[];
    Y_out=[];
    X_out_mean=[];
    Y_out_mean=[];
    unique_dates_out=[];
    subplot(2,2,[1,2]);
    contourf(repmat(mod_time,1,nlayers),depth(:,:,istat)*-1,DO_mod(:,:,istat),'LineStyle','none')
    hold on
    [cmin,cmax]=caxis;
    
    for iday=1:length(ib);
        
        % get the good day from the data back out so we can get
        % each value and depth
        myIds=find(cbp_station_date==unique_dates(ib(iday)));
        
        myDates=cbp_station_date(myIds);
        myDepths=cbp_station_depths(myIds);
        myvalues=cbp_station_value(myIds);
        
        % get rid of duplicates;
        [myUniqueDepths,Zids]=unique(myDepths);
        myUniqueValues=myvalues(Zids);
        myUniqueDates=unique(myDates);
        
        % now match up the model
        model_depthIds=dsearchn(depth(ia(iday),:,istat)',myUniqueDepths)  ;
        modelDepths=depth(ia(iday),model_depthIds,istat);
        modelValues=DO_mod(ia(iday),model_depthIds,istat);
        
        for iplot=1:length(modelDepths);
            scatter(myUniqueDates,modelDepths(iplot)*-1,100,myUniqueValues(iplot),'filled','markeredgecolor','k');
        end
        
        X_out=[X_out;myUniqueValues];
        Y_out=[Y_out;modelValues'];
        
        % take the depth averaged quantity
        X_out_mean=[X_out_mean;mean(myUniqueValues)];
        Y_out_mean=[Y_out_mean;mean(modelValues)];
        unique_dates_out=[unique_dates_out;myUniqueDates];
        
    end
    h= colorbar;
    [MEF,r,WMS,RMSE,MPE,RI,axis_min,axis_max]=model_statistics(X_out,Y_out)   ;
    caxis([axis_min axis_max]);
    h.Label.String='DO (mg O_2 l^-^1)'       ;
    
    
    title([ 'DO ' , '  MEF =' num2str(MEF),', r = ',num2str(r) ', RMSE = ' num2str(RMSE) , ', MPE = ',num2str(MPE)  ', WMS = ' num2str(WMS),', RI = ',num2str(RI)]);
    xlabel('Days in 2005');ylabel('Depth (m)')  ;
    
    subplot(2,2,3)
    plot(X_out,Y_out,'kd','markersize',8);
    axis([axis_min axis_max axis_min axis_max]);
    xlabel('Observed DO (mg O_2 l^-^1)');ylabel('Modeled DO (mg O_2 l^-^1)')  ;
    
    subplot(2,2,4);
    %     histogram(Y_out,10,'normalization','probability');
    %     hold on
    %     histogram(X_out,10,'normalization','probability');
    g1={};g2={};g3={};
    VAR_cat=reshape(DO_mod(:,:,istat),1,nlayers*ntimes);
    g1(1:length(X_out))={'Observed'};% get cell arrays of our group names
    g2(1:length(Y_out))={'Modeled'};
    g3(1:length(VAR_cat))={'Modeled-ALL'};
    grouping=[g1 g2 g3];
    boxplot([X_out ; Y_out; VAR_cat'],grouping);
    ylabel('DO (mg O_2 l^-^1)');%ylabel('Probability')
    %     xlabel('Chl a (ug l^-^1)');ylabel('Probability')
    
    saveas(gcf,[statstring '/DO_depth'],'png');
    saveas(gcf,[statstring '/DO_depth'],'eps');
    saveas(gcf,[statstring '/DO_depth'],'fig');
    
    figure;     % depth averaged
    
    [MEF,r,WMS,RMSE,MPE,RI,axis_min,axis_max]=model_statistics(X_out_mean,Y_out_mean)   ;
    subplot(2,2,[1,2]);
    plot(mod_time,mean(DO_mod(:,:,istat),2),'k-o');
    hold on
    plot(unique_dates_out,X_out_mean,'r-*');
    title([ 'Depth Averaged DO ' ,  '  MEF =' num2str(MEF),', r = ',num2str(r) ', RMSE = ' num2str(RMSE) , ', MPE = ',num2str(MPE)  ', WMS = ' num2str(WMS),', RI = ',num2str(RI)]);
    xlabel('Days in 2005');ylabel('Observed DO (mg O_2 l^-^1)');
    legend('Modeled','Observed');
    subplot(2,2,3);
    plot(X_out_mean,Y_out_mean,'rd','markersize',10);
    axis([axis_min axis_max axis_min axis_max]);
    xlabel('Observed DO (mg O_2 l^-^1)');ylabel('Modeled DO (mg O_2 l^-^1)')  ;
    subplot(2,2,4);
    %     histogram(Y_out,10,'normalization','probability');
    %     hold on
    %     histogram(X_out,10,'normalization','probability');
    g1={};g2={};g3={};
    VAR_cat=reshape(mean(DO_mod(:,:,istat),2),1,ntimes);
    g1(1:length(X_out_mean))={'Observed'};% get cell arrays of our group names
    g2(1:length(Y_out_mean))={'Modeled'};
    g3(1:length(VAR_cat))={'Modeled-ALL'};
    grouping=[g1 g2 g3];
    boxplot([X_out_mean ; Y_out_mean; VAR_cat'],grouping);
    ylabel('DO (mg O_2 l^-^1)');%ylabel('Probability')
    %     xlabel('Chl a (ug l^-^1)');ylabel('Probability')
    saveas(gcf,[statstring '/DO_depth_stats'],'png');
    saveas(gcf,[statstring '/DO_depth_stats'],'eps');
    saveas(gcf,[statstring '/DO_depth_stats'],'fig');
    
    %
    %========================= Nitrate ============================
    
    myvar_id=[];
    ivar=0;
    while isempty(myvar_id);
        ivar=ivar+1;
        if(strcmp(stations(istat).variables(ivar).name,'NO23F'))
            myvar_id=ivar              ; % change this ID based on where the variable is located within the CBP data structure
            % generated by read_cbp.m
        end
    end
    
    % pull out the data for this station
    [unique_dates,unique_ids]=unique(stations(istat).variables(myvar_id).date);
    cbp_station_date= stations(istat).variables(myvar_id).date ;
    cbp_station_value= stations(istat).variables(myvar_id).value;
    cbp_station_depths=stations(istat).variables(myvar_id).depth;
    % now for each unique date, find corresponding values and take the
    
    % get the matching model times for statistics and plotting;
    [matching_mod_times,ia,ib]=intersect(rounded_mod_times,unique_dates);
    
    
    % now plot date by date
    figure;
    
    X_out=[];
    Y_out=[];
    X_out_mean=[];
    Y_out_mean=[];
    unique_dates_out=[];
    
    subplot(2,2,[1,2])
    contourf(repmat(mod_time,1,nlayers),depth(:,:,istat)*-1,NO3_mod(:,:,istat),'LineStyle','none')
    hold on
    [cmin,cmax]=caxis;
    
    for iday=1:length(ib);
        
        % get the good day from the data back out so we can get
        % each value and depth
        myIds=find(cbp_station_date==unique_dates(ib(iday)));
        
        myDates=cbp_station_date(myIds);
        myDepths=cbp_station_depths(myIds);
        myvalues=cbp_station_value(myIds);
        
        % get rid of duplicates;
        [myUniqueDepths,Zids]=unique(myDepths);
        myUniqueValues=myvalues(Zids);
        myUniqueDates=unique(myDates);
        
        % now match up the model
        model_depthIds=dsearchn(depth(ia(iday),:,istat)',myUniqueDepths)  ;
        modelDepths=depth(ia(iday),model_depthIds,istat);
        modelValues=NO3_mod(ia(iday),model_depthIds,istat);
        
        for iplot=1:length(modelDepths);
            scatter(myUniqueDates,modelDepths(iplot)*-1,100,myUniqueValues(iplot),'filled','markeredgecolor','k');
        end
        X_out=[X_out;myUniqueValues];
        Y_out=[Y_out;modelValues'];
        
        
        % take the depth averaged quantity
        X_out_mean=[X_out_mean;mean(myUniqueValues)];
        Y_out_mean=[Y_out_mean;mean(modelValues)];
        
        unique_dates_out=[unique_dates_out;myUniqueDates];
    end
    h= colorbar;
    [MEF,r,WMS,RMSE,MPE,RI,axis_min,axis_max]=model_statistics(X_out,Y_out)   ;
    caxis([axis_min axis_max]);
    h.Label.String='NO_3^- (mg N l^-^1)'       ;
    
    title([ 'NO_3^- ' , '  MEF =' num2str(MEF),', r = ',num2str(r) ', RMSE = ' num2str(RMSE) , ', MPE = ',num2str(MPE)  ', WMS = ' num2str(WMS),', RI = ',num2str(RI)]);
    xlabel('Days in 2005');ylabel('NO_3^- (mg N l^-^1)')
    
    subplot(2,2,3);
    plot(X_out,Y_out,'kd','markersize',8);
    axis([axis_min axis_max axis_min axis_max]);
    xlabel('Observed NO_3^- (mg N l^-^1)');ylabel('Modeled NO_3^- (mg N l^-^1)');
    subplot(2,2,4);
    %     histogram(Y_out,10,'normalization','probability');
    %     hold on
    %     histogram(X_out,10,'normalization','probability');
    g1={};g2={};g3={};
    VAR_cat=reshape(NO3_mod(:,:,istat),1,nlayers*ntimes);
    g1(1:length(X_out))={'Observed'};% get cell arrays of our group names
    g2(1:length(Y_out))={'Modeled'};
    g3(1:length(VAR_cat))={'Modeled-ALL'};
    grouping=[g1 g2 g3];
    boxplot([X_out ; Y_out; VAR_cat'],grouping);
    
    ylabel('NO_3(mg N l^-^1)');%ylabel('Probability')
    %     xlabel('Chl a (ug l^-^1)');ylabel('Probability')
    saveas(gcf,[statstring '/NO3_depth'],'png');
    saveas(gcf,[statstring '/NO3_depth'],'eps');
    saveas(gcf,[statstring '/NO3_depth'],'fig');
    
    
    figure;     % depth averaged
    
    [MEF,r,WMS,RMSE,MPE,RI,axis_min,axis_max]=model_statistics(X_out_mean,Y_out_mean)   ;
    subplot(2,2,[1,2]);
    plot(mod_time,mean(NO3_mod(:,:,istat),2),'k-o');
    hold on
    plot(unique_dates_out,X_out_mean,'r-*');
    title([ 'Depth Averaged NO_3^- ' , '  MEF =' num2str(MEF),', r = ',num2str(r) ', RMSE = ' num2str(RMSE) , ', MPE = ',num2str(MPE)  ', WMS = ' num2str(WMS),', RI = ',num2str(RI)]);
    xlabel('Days in 2005');ylabel('NO_3^- (mg N l^-^1)');
    legend('Modeled','Observed');
    subplot(2,2,3);
    plot(X_out_mean,Y_out_mean,'rd','markersize',10);
    axis([axis_min axis_max axis_min axis_max]);
    xlabel('Observed NO_3^- (mg N l^-^1)');ylabel('Modeled NO_3^- (mg N l^-^1)');
    subplot(2,2,4);
    %     histogram(Y_out,10,'normalization','probability');
    %     hold on
    %     histogram(X_out,10,'normalization','probability');
    g1={};g2={};g3={};
    VAR_cat=reshape(mean(NO3_mod(:,:,istat),2),1,ntimes);
    g1(1:length(X_out_mean))={'Observed'};% get cell arrays of our group names
    g2(1:length(Y_out_mean))={'Modeled'};
    g3(1:length(VAR_cat))={'Modeled-ALL'};
    grouping=[g1 g2 g3];
    boxplot([X_out_mean ; Y_out_mean; VAR_cat'],grouping);
    ylabel('NO_3(mg N l^-^1)');%ylabel('Probability')
    %     xlabel('Chl a (ug l^-^1)');ylabel('Probability')
    
    saveas(gcf,[statstring '/NO3_depth_stats'],'png');
    saveas(gcf,[statstring '/NO3_depth_stats'],'eps');
    saveas(gcf,[statstring '/NO3_depth_stats'],'fig');
    %
    %========================= Ammonium ============================
    myvar_id=[];
    ivar=0;
    while isempty(myvar_id);
        ivar=ivar+1;
        if(strcmp(stations(istat).variables(ivar).name,'NH4F'));
            myvar_id=ivar              ; % change this ID based on where the variable is located within the CBP data structure
            % generated by read_cbp.m
        end
    end
    % pull out the data for this station
    [unique_dates,unique_ids]=unique(stations(istat).variables(myvar_id).date);
    cbp_station_date= stations(istat).variables(myvar_id).date ;
    cbp_station_value= stations(istat).variables(myvar_id).value;
    cbp_station_depths=stations(istat).variables(myvar_id).depth;
    % now for each unique date, find corresponding values and take the
    
    % get the matching model times for statistics and plotting;
    [matching_mod_times,ia,ib]=intersect(rounded_mod_times,unique_dates);
    
    
    % now plot date by date
    figure;
    
    X_out=[];
    Y_out=[];
    X_out_mean=[];
    Y_out_mean=[];
    unique_dates_out=[];
    
    subplot(2,2,[1,2]);
    contourf(repmat(mod_time,1,nlayers),depth(:,:,istat)*-1,NH4_mod(:,:,istat),'LineStyle','none');
    hold on;
    [cmin,cmax]=caxis;
    
    for iday=1:length(ib);
        
        % get the good day from the data back out so we can get
        % each value and depth
        myIds=find(cbp_station_date==unique_dates(ib(iday)));
        
        myDates=cbp_station_date(myIds);
        myDepths=cbp_station_depths(myIds);
        myvalues=cbp_station_value(myIds);
        
        % get rid of duplicates;
        [myUniqueDepths,Zids]=unique(myDepths);
        myUniqueValues=myvalues(Zids);
        myUniqueDates=unique(myDates);
        
        % now match up the model
        model_depthIds=dsearchn(depth(ia(iday),:,istat)',myUniqueDepths)  ;
        modelDepths=depth(ia(iday),model_depthIds,istat);
        modelValues=NH4_mod(ia(iday),model_depthIds,istat);
        
        for iplot=1:length(modelDepths);
            scatter(myUniqueDates,modelDepths(iplot)*-1,100,myUniqueValues(iplot),'filled','markeredgecolor','k');
        end
        X_out=[X_out;myUniqueValues];
        Y_out=[Y_out;modelValues'];
        
        
        % take the depth averaged quantity
        X_out_mean=[X_out_mean;mean(myUniqueValues)];
        Y_out_mean=[Y_out_mean;mean(modelValues)];
        
        unique_dates_out=[unique_dates_out;myUniqueDates];
        
    end
    h= colorbar;
    [MEF,r,WMS,RMSE,MPE,RI,axis_min,axis_max]=model_statistics(X_out,Y_out)   ;
    caxis([axis_min axis_max])
    h.Label.String='NH_4^+ (mg N l^-^1)'       ;
    
    
    
    title([ 'NH_4^+ ' , '  MEF =' num2str(MEF),', r = ',num2str(r) ', RMSE = ' num2str(RMSE) , ', MPE = ',num2str(MPE)  ', WMS = ' num2str(WMS),', RI = ',num2str(RI)]);
    xlabel('Days in 2005');ylabel('Depth (m)')  ;
    
    subplot(2,2,3)
    plot(X_out,Y_out,'kd','markersize',8);
    axis([axis_min axis_max axis_min axis_max]);
    xlabel('Observed NH_4^+ (mg N l^-^1)');ylabel('Modeled NH_4^+ (mg N l^-^1)')  ;
    subplot(2,2,4);
    %     histogram(Y_out,10,'normalization','probability');
    %     hold on
    %     histogram(X_out,10,'normalization','probability');
    g1={};g2={};g3={};
    VAR_cat=reshape(NH4_mod(:,:,istat),1,nlayers*ntimes);
    g1(1:length(X_out))={'Observed'};% get cell arrays of our group names
    g2(1:length(Y_out))={'Modeled'};
    g3(1:length(VAR_cat))={'Modeled-ALL'};
    grouping=[g1 g2 g3];
    boxplot([X_out ; Y_out; VAR_cat'],grouping);
    
    ylabel('NH_4(mg N l^-^1)');%ylabel('Probability')
    %     xlabel('Chl a (ug l^-^1)');ylabel('Probability')
    saveas(gcf,[statstring '/NH4_depth'],'png');
    saveas(gcf,[statstring '/NH4_depth'],'eps');
    saveas(gcf,[statstring '/NH4_depth'],'fig');
    
    figure;     % depth averaged
    
    [MEF,r,WMS,RMSE,MPE,RI,axis_min,axis_max]=model_statistics(X_out_mean,Y_out_mean)   ;
    subplot(2,2,[1,2]);
    plot(mod_time,mean(NH4_mod(:,:,istat),2),'k-o');
    hold on
    plot(unique_dates_out,X_out_mean,'r-*');
    title([ 'Depth Averaged NH_4^+' , '  MEF =' num2str(MEF),', r = ',num2str(r) ', RMSE = ' num2str(RMSE) , ', MPE = ',num2str(MPE)  ', WMS = ' num2str(WMS),', RI = ',num2str(RI)]);
    xlabel('Days in 2005');ylabel('NH_4^+ (mg N l^-^1)')  ;
    legend('Modeled','Observed')
    subplot(2,2,3);
    plot(X_out_mean,Y_out_mean,'rd','markersize',10);
    axis([axis_min axis_max axis_min axis_max]);
    xlabel('Observed NH_4^+ (mg N l^-^1)');ylabel('Modeled NH_4^+ (mg N l^-^1)')  ;
    subplot(2,2,4);
    %     histogram(Y_out,10,'normalization','probability');
    %     hold on
    %     histogram(X_out,10,'normalization','probability');
    g1={};g2={};g3={};
    VAR_cat=reshape(mean(NH4_mod(:,:,istat),2),1,ntimes);
    g1(1:length(X_out_mean))={'Observed'};% get cell arrays of our group names
    g2(1:length(Y_out_mean))={'Modeled'};
    g3(1:length(VAR_cat))={'Modeled-ALL'};
    grouping=[g1 g2 g3];
    boxplot([X_out_mean ; Y_out_mean; VAR_cat'],grouping);
    
    ylabel('NH_4(mg N l^-^1)');%ylabel('Probability')
    %     xlabel('Chl a (ug l^-^1)');ylabel('Probability')
    
    saveas(gcf,[statstring '/NH4_depth_stats'],'png');
    saveas(gcf,[statstring '/NH4_depth_stats'],'eps');
    saveas(gcf,[statstring '/NH4_depth_stats'],'fig');
    
    %========================= DON ============================
    
    myvar_id=[];
    ivar=0;
    while isempty(myvar_id)
        ivar=ivar+1;
        if(strcmp(stations(istat).variables(ivar).name,'DON'))
            myvar_id=ivar              ; % change this ID based on where the variable is located within the CBP data structure
            % generated by read_cbp.m
        end
    end
    % pull out the data for this station
    [unique_dates,unique_ids]=unique(stations(istat).variables(myvar_id).date);
    cbp_station_date= stations(istat).variables(myvar_id).date ;
    cbp_station_value= stations(istat).variables(myvar_id).value;
    cbp_station_depths=stations(istat).variables(myvar_id).depth;
    % now for each unique date, find corresponding values and take the
    
    % get the matching model times for statistics and plotting;
    [matching_mod_times,ia,ib]=intersect(rounded_mod_times,unique_dates);
     
    % now plot date by date
    figure;
    
    X_out=[];
    Y_out=[];
    X_out_mean=[];
    Y_out_mean=[];
    unique_dates_out=[];
    
    subplot(2,2,[1,2]);
    contourf(repmat(mod_time,1,nlayers),depth(:,:,istat)*-1,DON_mod(:,:,istat),'LineStyle','none')
    hold on
    [cmin,cmax]=caxis;
    
    for iday=1:length(ib)
        
        % get the good day from the data back out so we can get
        % each value and depth
        myIds=find(cbp_station_date==unique_dates(ib(iday)));
        
        myDates=cbp_station_date(myIds);
        myDepths=cbp_station_depths(myIds);
        myvalues=cbp_station_value(myIds);
        
        % get rid of duplicates;
        [myUniqueDepths,Zids]=unique(myDepths);
        myUniqueValues=myvalues(Zids);
        myUniqueDates=unique(myDates);
        
        % now match up the model
        model_depthIds=dsearchn(depth(ia(iday),:,istat)',myUniqueDepths)  ;
        modelDepths=depth(ia(iday),model_depthIds,istat);
        modelValues=DON_mod(ia(iday),model_depthIds,istat);
        
        for iplot=1:length(modelDepths);
            scatter(myUniqueDates,modelDepths(iplot)*-1,100,myUniqueValues(iplot),'filled','markeredgecolor','k');
        end
        X_out=[X_out;myUniqueValues];
        Y_out=[Y_out;modelValues'];
        
        % take the depth averaged quantity
        X_out_mean=[X_out_mean;mean(myUniqueValues)];
        Y_out_mean=[Y_out_mean;mean(modelValues)];
        
        unique_dates_out=[unique_dates_out;myUniqueDates];
    end
    h= colorbar;
    [MEF,r,WMS,RMSE,MPE,RI,axis_min,axis_max]=model_statistics(X_out,Y_out)   ;
    caxis([axis_min axis_max]);
    h.Label.String='DON (mg N l^-^1)'       ;
    
    title([ 'DON ' ,  '  MEF =' num2str(MEF),', r = ',num2str(r) ', RMSE = ' num2str(RMSE) , ', MPE = ',num2str(MPE)  ', WMS = ' num2str(WMS),', RI = ',num2str(RI)]);
    xlabel('Days in 2005');ylabel('Depth (m)')
    
    subplot(2,2,3);
    plot(X_out,Y_out,'kd','markersize',8);
    axis([axis_min axis_max axis_min axis_max]);;
    %     title([ 'DON ' , 'MEF =' num2str(MEF) ', RMSE = ' num2str(RMSE)   ' WMS = ' num2str(WMS)]);
    xlabel('Observed DON (mg N l^-^1)');ylabel('Modeled DON (mg N l^-^1)')  ;
    subplot(2,2,4);
    %     histogram(Y_out,10,'normalization','probability');
    %     hold on
    %     histogram(X_out,10,'normalization','probability');
    g1={};g2={};g3={};
    VAR_cat=reshape(DON_mod(:,:,istat),1,nlayers*ntimes);
    g1(1:length(X_out))={'Observed'};% get cell arrays of our group names
    g2(1:length(Y_out))={'Modeled'};
    g3(1:length(VAR_cat))={'Modeled-ALL'};
    grouping=[g1 g2 g3];
    boxplot([X_out ; Y_out; VAR_cat'],grouping);
    
    ylabel('DON (mg N l^-^1)');%ylabel('Probability')
    %     xlabel('Chl a (ug l^-^1)');ylabel('Probability')
    saveas(gcf,[statstring '/DON_depth'],'png');
    saveas(gcf,[statstring '/DON_depth'],'eps');
    saveas(gcf,[statstring '/DON_depth'],'fig');
    
    
    figure;     % depth averaged
    [MEF,r,WMS,RMSE,MPE,RI,axis_min,axis_max]=model_statistics(X_out_mean,Y_out_mean)   ;
    subplot(2,2,[1,2]);
    plot(mod_time,mean(DON_mod(:,:,istat),2),'k-o');
    hold on
    plot(unique_dates_out,X_out_mean,'r-*');
    title([ 'Depth Averaged DON ' , '  MEF =' num2str(MEF),', r = ',num2str(r) ', RMSE = ' num2str(RMSE) , ', MPE = ',num2str(MPE)  ', WMS = ' num2str(WMS),', RI = ',num2str(RI)]);
    xlabel('Days in 2005');ylabel('DON(mg N l^-^1)')
    legend('Modeled','Observed')
    subplot(2,2,3);
    plot(X_out_mean,Y_out_mean,'rd','markersize',10);
    axis([axis_min axis_max axis_min axis_max]);
    xlabel('Observed DON (mg N l^-^1)');ylabel('Modeled DON (mg N l^-^1)')  ;
    subplot(2,2,4);
    %     histogram(Y_out,10,'normalization','probability');
    %     hold on
    %     histogram(X_out,10,'normalization','probability');
    g1={};g2={};g3={};
    VAR_cat=reshape(mean(DON_mod(:,:,istat),2),1,ntimes);
    g1(1:length(X_out_mean))={'Observed'};% get cell arrays of our group names
    g2(1:length(Y_out_mean))={'Modeled'};
    g3(1:length(VAR_cat))={'Modeled-ALL'};
    grouping=[g1 g2 g3];
    boxplot([X_out_mean ; Y_out_mean; VAR_cat'],grouping);
    
    ylabel('DON (mg N l^-^1)');%ylabel('Probability')
    %     xlabel('Chl a (ug l^-^1)');ylabel('Probability')
    
    saveas(gcf,[statstring '/DON_depth_stats'],'png');
    saveas(gcf,[statstring '/DON_depth_stats'],'eps');
    saveas(gcf,[statstring '/DON_depth_stats'],'fig');
    %
    %========================= DOC ============================
    myvar_id=[];
    ivar=0;
    while isempty(myvar_id)
        ivar=ivar+1;
        if(strcmp(stations(istat).variables(ivar).name,'DOC'));
            
            myvar_id=ivar              ; % change this ID based on where the variable is located within the CBP data structure
            % generated by read_cbp.
            break
        end
    end
    CtoN=7.2
    % pull out the data for this station
    
    if(isempty(stations(istat).variables(myvar_id).value));  % pull out the DON data and use a C to N ratio if there is no data
        myvar_id=[];
        ivar=0;
        while isempty(myvar_id);
            ivar=ivar+1;
            if(strcmp(stations(istat).variables(ivar).name,'DON'));
                myvar_id=ivar              ; % change this ID based on where the variable is located within the CBP data structure
                % generated by read_cbp.m
            end
        end
        [unique_dates,unique_ids]=unique(stations(istat).variables(myvar_id).date);
        cbp_station_date= stations(istat).variables(myvar_id).date ;
        cbp_station_value= stations(istat).variables(myvar_id).value*CtoN;
        cbp_station_depths=stations(istat).variables(myvar_id).depth;
        % now for each unique date, find corresponding values and take the
    else
        
        [unique_dates,unique_ids]=unique(stations(istat).variables(myvar_id).date);
        cbp_station_date= stations(istat).variables(myvar_id).date ;
        cbp_station_value= stations(istat).variables(myvar_id).value;
        cbp_station_depths=stations(istat).variables(myvar_id).depth;
        % now for each unique date, find corresponding values and take the
    end
    % get the matching model times for statistics and plotting;
    [matching_mod_times,ia,ib]=intersect(rounded_mod_times,unique_dates);
    
    
    % now plot date by date
    figure;
    X_out=[];
    Y_out=[];
    X_out_mean=[];
    Y_out_mean=[];
    unique_dates_out=[];
    
    subplot(2,2,[1,2]);
    contourf(repmat(mod_time,1,nlayers),depth(:,:,istat)*-1,DOC_mod(:,:,istat),'LineStyle','none');
    hold on
    [cmin,cmax]=caxis;
    
    for iday=1:length(ib);
        
        % get the good day from the data back out so we can get
        % each value and depth
        myIds=find(cbp_station_date==unique_dates(ib(iday)));
        
        myDates=cbp_station_date(myIds);
        myDepths=cbp_station_depths(myIds);
        myvalues=cbp_station_value(myIds);
        
        % get rid of duplicates;
        [myUniqueDepths,Zids]=unique(myDepths);
        myUniqueValues=myvalues(Zids);
        myUniqueDates=unique(myDates);
        
        % now match up the model
        model_depthIds=dsearchn(depth(ia(iday),:,istat)',myUniqueDepths)  ;
        modelDepths=depth(ia(iday),model_depthIds,istat);
        modelValues=DOC_mod(ia(iday),model_depthIds,istat);
        
        for iplot=1:length(modelDepths);
            scatter(myUniqueDates,modelDepths(iplot)*-1,100,myUniqueValues(iplot),'filled','markeredgecolor','k');
        end
        X_out=[X_out;myUniqueValues];
        Y_out=[Y_out;modelValues'];
        
        % take the depth averaged quantity
        X_out_mean=[X_out_mean;mean(myUniqueValues)];
        Y_out_mean=[Y_out_mean;mean(modelValues)];
        
        unique_dates_out=[unique_dates_out;myUniqueDates];
        
    end
    h= colorbar;
    [MEF,r,WMS,RMSE,MPE,RI,axis_min,axis_max]=model_statistics(X_out,Y_out)   ;
    caxis([axis_min axis_max])
    h.Label.String='DOC (mg C l^-^1)'       ;
    
    title([ 'DOC ' ,  '  MEF =' num2str(MEF),', r = ',num2str(r) ', RMSE = ' num2str(RMSE) , ', MPE = ',num2str(MPE)  ', WMS = ' num2str(WMS),', RI = ',num2str(RI)]);
    xlabel('Days in 2005');ylabel('Depth (m)')  ;
    subplot(2,2,3)
    plot(X_out,Y_out,'kd','markersize',8);
    axis([axis_min axis_max axis_min axis_max]);
    %     title([ 'DOC ' , 'r =' num2str(r) ', RMSE = ' num2str(RMSE)   ' WMS = ' num2str(WMS)]);
    xlabel('Observed DOC (mg C l^-^1)');ylabel('Modeled DOC (mg C l^-^1)')  ;
    subplot(2,2,4);
    %     histogram(Y_out,10,'normalization','probability');
    %     hold on
    %     histogram(X_out,10,'normalization','probability');
    g1={};g2={};g3={};
    VAR_cat=reshape(DOC_mod(:,:,istat),1,nlayers*ntimes);
    g1(1:length(X_out))={'Observed'};% get cell arrays of our group names
    g2(1:length(Y_out))={'Modeled'};
    g3(1:length(VAR_cat))={'Modeled-ALL'};
    grouping=[g1 g2 g3];
    boxplot([X_out ; Y_out; VAR_cat'],grouping);
    
    ylabel('DOC (mg C l^-^1)');%ylabel('Probability')
    %     xlabel('Chl a (ug l^-^1)');ylabel('Probability')
    saveas(gcf,[statstring '/DOC_depth'],'png');
    saveas(gcf,[statstring '/DOC_depth'],'eps');
    saveas(gcf,[statstring '/DOC_depth'],'fig');
    
    figure;     % depth averaged
    [MEF,r,WMS,RMSE,MPE,RI,axis_min,axis_max]=model_statistics(X_out_mean,Y_out_mean)   ;
    subplot(2,2,[1,2]);
    plot(mod_time,mean(DOC_mod(:,:,istat),2),'k-o');
    hold on
    plot(unique_dates_out,X_out_mean,'r-*');
    title([ 'Depth Averaged DOC ' , '  MEF =' num2str(MEF),', r = ',num2str(r) ', RMSE = ' num2str(RMSE) , ', MPE = ',num2str(MPE)  ', WMS = ' num2str(WMS),', RI = ',num2str(RI)]);
    xlabel('Days in 2005');ylabel('DOC(mg C l^-^1)');
    legend('Modeled','Observed');
    subplot(2,2,3);
    plot(X_out_mean,Y_out_mean,'rd','markersize',10);
    axis([axis_min axis_max axis_min axis_max]);
    xlabel('Observed DOC (mg C l^-^1)');ylabel('Modeled DOC (mg C l^-^1)')  ;
    subplot(2,2,4);
    %     histogram(Y_out,10,'normalization','probability');
    %     hold on
    %     histogram(X_out,10,'normalization','probability');
    g1={};g2={};g3={};
    VAR_cat=reshape(mean(DOC_mod(:,:,istat),2),1,ntimes);
    g1(1:length(X_out_mean))={'Observed'};% get cell arrays of our group names
    g2(1:length(Y_out_mean))={'Modeled'};
    g3(1:length(VAR_cat))={'Modeled-ALL'};
    grouping=[g1 g2 g3];
    boxplot([X_out_mean ; Y_out_mean; VAR_cat'],grouping);
    
    ylabel('DOC (mg C l^-^1)');%ylabel('Probability')
    
    saveas(gcf,[statstring '/DOC_depth_stats'],'png');
    saveas(gcf,[statstring '/DOC_depth_stats'],'eps');
    saveas(gcf,[statstring '/DOC_depth_stats'],'fig');
    %
    %========================= Salinity ============================
    
    myvar_id=[];
    ivar=0;
    while isempty(myvar_id);
        ivar=ivar+1;
        if(strcmp(stations(istat).variables(ivar).name,'SALINITY'));
            
            myvar_id=ivar              ; % change this ID based on where the variable is located within the CBP data structure
            % generated by read_cbp.m
        end
    end
    % pull out the data for this station
    [unique_dates,unique_ids]=unique(stations(istat).variables(myvar_id).date);
    cbp_station_date= stations(istat).variables(myvar_id).date ;
    cbp_station_value= stations(istat).variables(myvar_id).value ;
    cbp_station_depths=stations(istat).variables(myvar_id).depth;
    % now for each unique date, find corresponding values and take the
    
    % get the matching model times for statistics and plotting;
    [matching_mod_times,ia,ib]=intersect(rounded_mod_times,unique_dates);
    
    
    % now plot date by date
    figure;
    X_out=[];
    Y_out=[];
    X_out_mean=[];
    Y_out_mean=[];
    unique_dates_out=[];
    
    subplot(2,2,[1,2]);
    contourf(repmat(mod_time,1,nlayers),depth(:,:,istat)*-1,SALT_mod(:,:,istat),'LineStyle','none');
    hold on
    [cmin,cmax]=caxis;
    
    for iday=1:length(ib);
        
        % get the good day from the data back out so we can get
        % each value and depth
        myIds=find(cbp_station_date==unique_dates(ib(iday)));
        
        myDates=cbp_station_date(myIds);
        myDepths=cbp_station_depths(myIds);
        myvalues=cbp_station_value(myIds);
        
        % get rid of duplicates;
        [myUniqueDepths,Zids]=unique(myDepths);
        myUniqueValues=myvalues(Zids);
        myUniqueDates=unique(myDates);
        
        % now match up the model
        model_depthIds=dsearchn(depth(ia(iday),:,istat)',myUniqueDepths)  ;
        modelDepths=depth(ia(iday),model_depthIds,istat);;
        modelValues=SALT_mod(ia(iday),model_depthIds,istat);
        
        for iplot=1:length(modelDepths);
            scatter(myUniqueDates,modelDepths(iplot)*-1,100,myUniqueValues(iplot),'filled','markeredgecolor','k');
        end
        X_out=[X_out;myUniqueValues];
        Y_out=[Y_out;modelValues'];
        
        % take the depth averaged quantity
        X_out_mean=[X_out_mean;mean(myUniqueValues)];
        Y_out_mean=[Y_out_mean;mean(modelValues)];
        
        unique_dates_out=[unique_dates_out;myUniqueDates];
        
        
    end
    h= colorbar;
    [MEF,r,WMS,RMSE,MPE,RI,axis_min,axis_max]=model_statistics(X_out,Y_out)   ;
    caxis([axis_min axis_max]);
    h.Label.String='Salinity'       ;
    
    title([ 'Salinity' ,  '  MEF =' num2str(MEF),', r = ',num2str(r) ', RMSE = ' num2str(RMSE) , ', MPE = ',num2str(MPE)  ', WMS = ' num2str(WMS),', RI = ',num2str(RI)]);
    xlabel('Days in 2005');ylabel('Depth (m)')  ;
    
    subplot(2,2,3)
    plot(X_out,Y_out,'kd','markersize',8);
    axis([axis_min axis_max axis_min axis_max]);
    %     title([ 'Salinity ' , 'r =' num2str(r) ', RMSE = ' num2str(RMSE)   ' WMS = ' num2str(WMS)]);
    xlabel('Observed Salinity');ylabel('Modeled Salinity ');
    subplot(2,2,4);
    %     histogram(Y_out,10,'normalization','probability');
    %     hold on
    %     histogram(X_out,10,'normalization','probability');
    g1={};g2={};g3={};
    VAR_cat=reshape(SALT_mod(:,:,istat),1,nlayers*ntimes);
    g1(1:length(X_out))={'Observed'};% get cell arrays of our group names
    g2(1:length(Y_out))={'Modeled'};
    g3(1:length(VAR_cat))={'Modeled-ALL'};
    grouping=[g1 g2 g3];
    boxplot([X_out ; Y_out; VAR_cat'],grouping);
    
    ylabel('Salinity');%ylabel('Probability')
    
    saveas(gcf,[statstring '/Salt_depth'],'png');
    saveas(gcf,[statstring '/Salt_depth'],'eps');
    saveas(gcf,[statstring '/Salt_depth'],'fig');
    
    figure;     % depth averaged
    
    [MEF,r,WMS,RMSE,MPE,RI,axis_min,axis_max]=model_statistics(X_out_mean,Y_out_mean)   ;
    subplot(2,2,[1,2]);
    plot(mod_time,mean(SALT_mod(:,:,istat),2),'k-o');
    hold on
    plot(unique_dates_out,X_out_mean,'r-*');
    title([ 'Depth Averaged Salinity ' ,  '  MEF =' num2str(MEF),', r = ',num2str(r) ', RMSE = ' num2str(RMSE) , ', MPE = ',num2str(MPE)  ', WMS = ' num2str(WMS),', RI = ',num2str(RI)]);
    xlabel('Days in 2005');ylabel('Salinity ');
    legend('Modeled','Observed')
    subplot(2,2,3);
    plot(X_out_mean,Y_out_mean,'rd','markersize',10);
    axis([axis_min axis_max axis_min axis_max]);
    xlabel('Observed Salinity');ylabel('Modeled Salinity')  ;
    subplot(2,2,4);
    %     histogram(Y_out,10,'normalization','probability');
    %     hold on
    %     histogram(X_out,10,'normalization','probability');
    g1={};g2={};g3={};
    VAR_cat=reshape(mean(SALT_mod(:,:,istat),2),1,ntimes);
    g1(1:length(X_out_mean))={'Observed'};% get cell arrays of our group names
    g2(1:length(Y_out_mean))={'Modeled'};
    g3(1:length(VAR_cat))={'Modeled-ALL'};
    grouping=[g1 g2 g3];
    boxplot([X_out_mean ; Y_out_mean; VAR_cat'],grouping);
    
    ylabel('Salinity');%ylabel('Probability')
    
    saveas(gcf,[statstring '/Salt_depth_stats'],'png');
    saveas(gcf,[statstring '/Salt_depth_stats'],'eps');
    saveas(gcf,[statstring '/Salt_depth_stats'],'fig');
    %========================= Temperature ============================
    
    myvar_id=[];
    ivar=0;
    while isempty(myvar_id);
        ivar=ivar+1;
        if(strcmp(stations(istat).variables(ivar).name,'WTEMP'));
            
            myvar_id=ivar              ; % change this ID based on where the variable is located within the CBP data structure
            % generated by read_cbp.m
        end
    end
    % pull out the data for this station
    [unique_dates,unique_ids]=unique(stations(istat).variables(myvar_id).date);
    cbp_station_date= stations(istat).variables(myvar_id).date ;
    cbp_station_value= stations(istat).variables(myvar_id).value ;
    cbp_station_depths=stations(istat).variables(myvar_id).depth;
    % now for each unique date, find corresponding values and take the
    
    % get the matching model times for statistics and plotting;
    [matching_mod_times,ia,ib]=intersect(rounded_mod_times,unique_dates);
    
    
    % now plot date by date
    figure;
    X_out=[];
    Y_out=[];
    X_out_mean=[];
    Y_out_mean=[];
    unique_dates_out=[];
    
    subplot(2,2,[1,2]);
    contourf(repmat(mod_time,1,nlayers),depth(:,:,istat)*-1,TEMP_mod(:,:,istat),'LineStyle','none');
    hold on
    [cmin,cmax]=caxis;
    
    for iday=1:length(ib)
        
        % get the good day from the data back out so we can get
        % each value and depth
        myIds=find(cbp_station_date==unique_dates(ib(iday)));
        
        myDates=cbp_station_date(myIds);
        myDepths=cbp_station_depths(myIds);
        myvalues=cbp_station_value(myIds);
        
        % get rid of duplicates;
        [myUniqueDepths,Zids]=unique(myDepths);
        myUniqueValues=myvalues(Zids);
        myUniqueDates=unique(myDates);
        
        % now match up the model
        model_depthIds=dsearchn(depth(ia(iday),:,istat)',myUniqueDepths)  ;
        modelDepths=depth(ia(iday),model_depthIds,istat);
        modelValues=TEMP_mod(ia(iday),model_depthIds,istat);
        
        for iplot=1:length(modelDepths);
            scatter(myUniqueDates,modelDepths(iplot)*-1,100,myUniqueValues(iplot),'filled','markeredgecolor','k');
        end
        X_out=[X_out;myUniqueValues];
        Y_out=[Y_out;modelValues'];
        
        % take the depth averaged quantity
        X_out_mean=[X_out_mean;mean(myUniqueValues)];
        Y_out_mean=[Y_out_mean;mean(modelValues)];
        
        unique_dates_out=[unique_dates_out;myUniqueDates];
        
        
    end
    h= colorbar;
    [MEF,r,WMS,RMSE,MPE,RI,axis_min,axis_max]=model_statistics(X_out,Y_out)   ;
    caxis([axis_min axis_max]);
    h.Label.String='Temperature'       ;
    
    title([ 'Temperature' ,  '  MEF =' num2str(MEF),', r = ',num2str(r) ', RMSE = ' num2str(RMSE) , ', MPE = ',num2str(MPE)  ', WMS = ' num2str(WMS),', RI = ',num2str(RI)]);
    xlabel('Days in 2005');ylabel('Depth (m)')  ;
    
    subplot(2,2,3)
    plot(X_out,Y_out,'kd','markersize',8);
    axis([axis_min axis_max axis_min axis_max]);
    %     title([ 'Salinity ' , 'r =' num2str(r) ', RMSE = ' num2str(RMSE)   ' WMS = ' num2str(WMS)]);
    xlabel('Observed Temperature');ylabel('Modeled Temperature');
    subplot(2,2,4);
    %     histogram(Y_out,10,'normalization','probability');
    %     hold on
    %     histogram(X_out,10,'normalization','probability');
    g1={};g2={};g3={};
    VAR_cat=reshape(TEMP_mod(:,:,istat),1,nlayers*ntimes);
    g1(1:length(X_out))={'Observed'};% get cell arrays of our group names
    g2(1:length(Y_out))={'Modeled'};
    g3(1:length(VAR_cat))={'Modeled-ALL'};
    grouping=[g1 g2 g3];
    boxplot([X_out ; Y_out; VAR_cat'],grouping);
    
    ylabel('Temperature');%ylabel('Probability')
    
    saveas(gcf,[statstring '/Temp_depth'],'png');
    saveas(gcf,[statstring '/Temp_depth'],'eps');
    saveas(gcf,[statstring '/Temp_depth'],'fig');
    
    figure;     % depth averaged
    
    [MEF,r,WMS,RMSE,MPE,RI,axis_min,axis_max]=model_statistics(X_out_mean,Y_out_mean)   ;
    subplot(2,2,[1,2]);
    plot(mod_time,mean(TEMP_mod(:,:,istat),2),'k-o');
    hold on
    plot(unique_dates_out,X_out_mean,'r-*');
    title([ 'Depth Averaged Temperature ' ,  '  MEF =' num2str(MEF),', r = ',num2str(r) ', RMSE = ' num2str(RMSE) , ', MPE = ',num2str(MPE)  ', WMS = ' num2str(WMS),', RI = ',num2str(RI)]);
    xlabel('Days in 2005');ylabel('Temperature')
    
    subplot(2,2,3);
    plot(X_out_mean,Y_out_mean,'rd','markersize',10);
    axis([axis_min axis_max axis_min axis_max]);
    xlabel('Observed Temperature');ylabel('Modeled Temperature')  ;
    subplot(2,2,4);
    %     histogram(Y_out,10,'normalization','probability');
    %     hold on
    %     histogram(X_out,10,'normalization','probability');
    g1={};g2={};g3={};
    VAR_cat=reshape(mean(TEMP_mod(:,:,istat),2),1,ntimes);
    g1(1:length(X_out_mean))={'Observed'};% get cell arrays of our group names
    g2(1:length(Y_out_mean))={'Modeled'};
    g3(1:length(VAR_cat))={'Modeled-ALL'};
    grouping=[g1 g2 g3];
    boxplot([X_out_mean ; Y_out_mean; VAR_cat'],grouping);
    ylabel('Temperature');%ylabel('Probability')
    
    saveas(gcf,[statstring '/Temp_depth_stats'],'png');
    saveas(gcf,[statstring '/Temp_depth_stats'],'eps');
    saveas(gcf,[statstring '/Temp_depth_stats'],'fig');
    
    %========================= POC ============================
    myvar_id=[];
    myvar_id2=[];
    ivar=0;
    while isempty(myvar_id);
        ivar=ivar+1;
        if(strcmp(stations(istat).variables(ivar).name,'DON'));
            myvar_id=ivar              ; % change this ID based on where the variable is located within the CBP data structure
            % generated by read_cbp.m
        end
    end
    
    ivar=0;
    while isempty(myvar_id2);
        ivar=ivar+1;
        if(strcmp(stations(istat).variables(ivar).name,'TON'));
            myvar_id2=ivar              ; % change this ID based on where the variable is located within the CBP data structure
            % generated by read_cbp.m
        end
    end
    
    % pull out the data for this station
    [unique_dates,unique_ids]=unique(stations(istat).variables(myvar_id).date);
    [unique_dates2,unique_ids2]=unique(stations(istat).variables(myvar_id2).date);
    
    cbp_station_date= stations(istat).variables(myvar_id).date ;
    cbp_station_date2= stations(istat).variables(myvar_id2).date ;
    
    cbp_station_value= stations(istat).variables(myvar_id).value;
    cbp_station_value2= stations(istat).variables(myvar_id2).value;
    
    cbp_station_depths=stations(istat).variables(myvar_id).depth;
    cbp_station_depths2=stations(istat).variables(myvar_id2).depth;
    % now for each unique date, find corresponding values and take the
    
    % get the matching model times for statistics and plotting;
    [matching_mod_times,ia,ib]=intersect(rounded_mod_times,unique_dates);
    
    if(length(unique_dates)~=length(unique_dates2));
        disp('Different days for PON and DON, exiting');
        return
    end
    
    
    % now plot date by date
    figure;
    
    X_out=[];
    Y_out=[];
    X_out_mean=[];
    Y_out_mean=[];
    unique_dates_out=[];
    
    subplot(2,2,[1,2]);
    contourf(repmat(mod_time,1,nlayers),depth(:,:,istat)*-1,POC_mod(:,:,istat),'LineStyle','none');
    hold on
    [cmin,cmax]=caxis;
    
    for iday=1:length(ib)
        
        % get the good day from the data back out so we can get
        % each value and depth
        myIds=find(cbp_station_date==unique_dates(ib(iday)));
        myIds2=find(cbp_station_date2==unique_dates(ib(iday)));
        
        myDates=cbp_station_date(myIds);
        myDepths=cbp_station_depths(myIds);
        
        myvalues=cbp_station_value(myIds);
        myvalues2=cbp_station_value2(myIds2);
        
        
        POC_in=(myvalues2-myvalues)*5.67; % TON-PON multiplied by redfield (g C g N^-1)
        
        % get rid of duplicates;
        [myUniqueDepths,Zids]=unique(myDepths);
        myUniqueValues=POC_in(Zids);
        myUniqueDates=unique(myDates);
        
        % now match up the model
        model_depthIds=dsearchn(depth(ia(iday),:,istat)',myUniqueDepths)  ;
        modelDepths=depth(ia(iday),model_depthIds,istat);
        modelValues=POC_mod(ia(iday),model_depthIds,istat);
        
        for iplot=1:length(modelDepths);
            scatter(myUniqueDates,modelDepths(iplot)*-1,100,myUniqueValues(iplot),'filled','markeredgecolor','k');
        end
        X_out=[X_out;myUniqueValues];
        Y_out=[Y_out;modelValues'];
        
        % take the depth averaged quantity
        X_out_mean=[X_out_mean;mean(myUniqueValues)];
        Y_out_mean=[Y_out_mean;mean(modelValues)];
        
        unique_dates_out=[unique_dates_out;myUniqueDates];
    end
    
    h= colorbar;
    [MEF,r,WMS,RMSE,MPE,RI,axis_min,axis_max]=model_statistics(X_out,Y_out)   ;
    caxis([axis_min axis_max]);
    h.Label.String='POC(mg C l^-^1)'       ;
    
    title([ 'POC ' ,  '  MEF =' num2str(MEF),', r = ',num2str(r) ', RMSE = ' num2str(RMSE) , ', MPE = ',num2str(MPE)  ', WMS = ' num2str(WMS),', RI = ',num2str(RI)]);
    xlabel('Days in 2005');ylabel('Depth (m)');
    
    subplot(2,2,3);
    plot(X_out,Y_out,'kd','markersize',8);
    axis([axis_min axis_max axis_min axis_max]);
    %     title([ 'DON ' , 'MEF =' num2str(MEF) ', RMSE = ' num2str(RMSE)   ' WMS = ' num2str(WMS)]);
    xlabel('Observed POC (mg C l^-^1)');ylabel('Modeled POC (mg c l^-^1)')  ;
    subplot(2,2,4);
    %     histogram(Y_out,10,'normalization','probability');
    %     hold on
    %     histogram(X_out,10,'normalization','probability');
    g1={};g2={};g3={};
    VAR_cat=reshape(POC_mod(:,:,istat),1,nlayers*ntimes);
    g1(1:length(X_out))={'Observed'};% get cell arrays of our group names
    g2(1:length(Y_out))={'Modeled'};
    g3(1:length(VAR_cat))={'Modeled-ALL'};
    grouping=[g1 g2 g3];
    boxplot([X_out ; Y_out; VAR_cat'],grouping);
    
    ylabel('POC (mg C l^-^1)');%ylabel('Probability')
    %     xlabel('Chl a (ug l^-^1)');ylabel('Probability')
    saveas(gcf,[statstring '/POC_depth'],'png');
    saveas(gcf,[statstring '/POC_depth'],'eps');
    saveas(gcf,[statstring '/POC_depth'],'fig');
    
    
    figure;     % depth averaged
    [MEF,r,WMS,RMSE,MPE,RI,axis_min,axis_max]=model_statistics(X_out_mean,Y_out_mean)   ;
    subplot(2,2,[1,2]);
    plot(mod_time,mean(POC_mod(:,:,istat),2),'k-o');
    hold on
    plot(unique_dates_out,X_out_mean,'r-*');
    title([ 'Depth Averaged POC ' , '  MEF =' num2str(MEF),', r = ',num2str(r) ', RMSE = ' num2str(RMSE) , ', MPE = ',num2str(MPE)  ', WMS = ' num2str(WMS),', RI = ',num2str(RI)]);
    xlabel('Days in 2005');ylabel('POC (mg C l^-^1)');
    legend('Modeled','Observed');
    subplot(2,2,3);
    plot(X_out_mean,Y_out_mean,'rd','markersize',10);
    axis([axis_min axis_max axis_min axis_max]);
    xlabel('Observed POC (mg C l^-^1)');ylabel('Modeled POC (mg C l^-^1)')  ;
    subplot(2,2,4);
    %     histogram(Y_out,10,'normalization','probability');
    %     hold on
    %     histogram(X_out,10,'normalization','probability');
    g1={};g2={};g3={};
    VAR_cat=reshape(mean(POC_mod(:,:,istat),2),1,ntimes);
    g1(1:length(X_out_mean))={'Observed'};% get cell arrays of our group names
    g2(1:length(Y_out_mean))={'Modeled'};
    g3(1:length(VAR_cat))={'Modeled-ALL'};
    grouping=[g1 g2 g3];
    boxplot([X_out_mean ; Y_out_mean; VAR_cat'],grouping);
    ylabel('POC (mg C l^-^1)');%ylabel('Probability')
    %      set(gca,'FontWeight','bold','LineWidth',2);
    %     xlabel('Chl a (ug l^-^1)');ylabel('Probability')
    
    saveas(gcf,[statstring '/POC_depth_stats'],'png');
    saveas(gcf,[statstring '/POC_depth_stats'],'eps');
    saveas(gcf,[statstring '/POC_depth_stats'],'fig');
    
    close all

    %========================= PON ============================
    myvar_id=[];
    myvar_id2=[];
    ivar=0;
    while isempty(myvar_id);
        ivar=ivar+1;
        if(strcmp(stations(istat).variables(ivar).name,'DON'));
            myvar_id=ivar              ; % change this ID based on where the variable is located within the CBP data structure
            % generated by read_cbp.m
        end
    end
    
    ivar=0;
    while isempty(myvar_id2);
        ivar=ivar+1;
        if(strcmp(stations(istat).variables(ivar).name,'TON'));
            myvar_id2=ivar              ; % change this ID based on where the variable is located within the CBP data structure
            % generated by read_cbp.m
        end
    end
    
    % pull out the data for this station
    [unique_dates,unique_ids]=unique(stations(istat).variables(myvar_id).date);
    [unique_dates2,unique_ids2]=unique(stations(istat).variables(myvar_id2).date);
    
    cbp_station_date= stations(istat).variables(myvar_id).date ;
    cbp_station_date2= stations(istat).variables(myvar_id2).date ;
    
    cbp_station_value= stations(istat).variables(myvar_id).value;
    cbp_station_value2= stations(istat).variables(myvar_id2).value;
    
    cbp_station_depths=stations(istat).variables(myvar_id).depth;
    cbp_station_depths2=stations(istat).variables(myvar_id2).depth;
    % now for each unique date, find corresponding values and take the
    
    % get the matching model times for statistics and plotting;
    [matching_mod_times,ia,ib]=intersect(rounded_mod_times,unique_dates);
    
    if(length(unique_dates)~=length(unique_dates2));
        disp('Different days for PON and DON, exiting');
        return
    end
    
    
    % now plot date by date
    figure;
    
    X_out=[];
    Y_out=[];
    X_out_mean=[];
    Y_out_mean=[];
    unique_dates_out=[];
    
    subplot(2,2,[1,2]);
    contourf(repmat(mod_time,1,nlayers),depth(:,:,istat)*-1,PON_mod(:,:,istat),'LineStyle','none')
    hold on
    [cmin,cmax]=caxis;
    
    for iday=1:length(ib)
        
        % get the good day from the data back out so we can get
        % each value and depth
        myIds=find(cbp_station_date==unique_dates(ib(iday)));
        myIds2=find(cbp_station_date2==unique_dates(ib(iday)));
        
        myDates=cbp_station_date(myIds);
        myDepths=cbp_station_depths(myIds);
        
        myvalues=cbp_station_value(myIds);
        myvalues2=cbp_station_value2(myIds2);
        
        
        PON_in=(myvalues2-myvalues); % TON-PON
        
        % get rid of duplicates;
        [myUniqueDepths,Zids]=unique(myDepths);
        myUniqueValues=PON_in(Zids);
        myUniqueDates=unique(myDates);
        
        % now match up the model
        model_depthIds=dsearchn(depth(ia(iday),:,istat)',myUniqueDepths)  ;
        modelDepths=depth(ia(iday),model_depthIds,istat);
        modelValues=PON_mod(ia(iday),model_depthIds,istat);
        
        for iplot=1:length(modelDepths);
            scatter(myUniqueDates,modelDepths(iplot)*-1,100,myUniqueValues(iplot),'filled','markeredgecolor','k');
        end
        X_out=[X_out;myUniqueValues];
        Y_out=[Y_out;modelValues'];
        
        % take the depth averaged quantity
        X_out_mean=[X_out_mean;mean(myUniqueValues)];
        Y_out_mean=[Y_out_mean;mean(modelValues)];
        
        unique_dates_out=[unique_dates_out;myUniqueDates];
    end
    
    h= colorbar;
    [MEF,r,WMS,RMSE,MPE,RI,axis_min,axis_max]=model_statistics(X_out,Y_out)   ;
    caxis([axis_min axis_max]);
    h.Label.String='PON(mg C l^-^1)'       ;
    
    title([ 'PON ' ,  '  MEF =' num2str(MEF),', r = ',num2str(r) ', RMSE = ' num2str(RMSE) , ', MPE = ',num2str(MPE)  ', WMS = ' num2str(WMS),', RI = ',num2str(RI)]);
    xlabel('Days in 2005');ylabel('Depth (m)')
    
    subplot(2,2,3);
    plot(X_out,Y_out,'kd','markersize',8);
    axis([axis_min axis_max axis_min axis_max]);;
    %     title([ 'DON ' , 'MEF =' num2str(MEF) ', RMSE = ' num2str(RMSE)   ' WMS = ' num2str(WMS)]);
    xlabel('Observed PON (mg N l^-^1)');ylabel('Modeled PON (mg N l^-^1)')  ;
    subplot(2,2,4);
    %     histogram(Y_out,10,'normalization','probability');
    %     hold on
    %     histogram(X_out,10,'normalization','probability');
    g1={};g2={};g3={};
    VAR_cat=reshape(PON_mod(:,:,istat),1,nlayers*ntimes);
    g1(1:length(X_out))={'Observed'};% get cell arrays of our group names
    g2(1:length(Y_out))={'Modeled'};
    g3(1:length(VAR_cat))={'Modeled-ALL'};
    grouping=[g1 g2 g3];
    boxplot([X_out ; Y_out; VAR_cat'],grouping);
    
    ylabel('POC (mg C l^-^1)');%ylabel('Probability')
    %     xlabel('Chl a (ug l^-^1)');ylabel('Probability')
    saveas(gcf,[statstring '/PON_depth'],'png');
    saveas(gcf,[statstring '/PON_depth'],'eps');
    saveas(gcf,[statstring '/PON_depth'],'fig');
    
    
    figure;     % depth averaged
    [MEF,r,WMS,RMSE,MPE,RI,axis_min,axis_max]=model_statistics(X_out_mean,Y_out_mean)   ;
    subplot(2,2,[1,2]);
    plot(mod_time,mean(PON_mod(:,:,istat),2),'k-o');
    hold on
    plot(unique_dates_out,X_out_mean,'r-*');
    title([ 'Depth Averaged PON ' , '  MEF =' num2str(MEF),', r = ',num2str(r) ', RMSE = ' num2str(RMSE) , ', MPE = ',num2str(MPE)  ', WMS = ' num2str(WMS),', RI = ',num2str(RI)]);
    xlabel('Days in 2005');ylabel('PON (mg N l^-^1)')
    legend('Modeled','Observed')
    subplot(2,2,3);
    plot(X_out_mean,Y_out_mean,'rd','markersize',10);
    axis([axis_min axis_max axis_min axis_max]);
    xlabel('Observed PON (mg N l^-^1)');ylabel('Modeled PON (mg N l^-^1)')  ;
    subplot(2,2,4);
    %     histogram(Y_out,10,'normalization','probability');
    %     hold on
    %     histogram(X_out,10,'normalization','probability');
    g1={};g2={};g3={};
    VAR_cat=reshape(mean(PON_mod(:,:,istat),2),1,ntimes);
    g1(1:length(X_out_mean))={'Observed'};% get cell arrays of our group names
    g2(1:length(Y_out_mean))={'Modeled'};
    g3(1:length(VAR_cat))={'Modeled-ALL'};
    grouping=[g1 g2 g3];
    boxplot([X_out_mean ; Y_out_mean; VAR_cat'],grouping);
    ylabel('POC (mg C l^-^1)');%ylabel('Probability')
    %      set(gca,'FontWeight','bold','LineWidth',2);
    %     xlabel('Chl a (ug l^-^1)');ylabel('Probability')
    
    saveas(gcf,[statstring '/PON_depth_stats'],'png');
    saveas(gcf,[statstring '/PON_depth_stats'],'eps');
    saveas(gcf,[statstring '/PON_depth_stats'],'fig');
    
    close all
 
end

%Finished




