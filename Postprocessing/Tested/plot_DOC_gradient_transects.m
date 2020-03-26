% Plot median concentration gradient int he Rhode River
% THis will generate an 8 panel plot with each panel
%the average dDOC/dx over the Rhode River

% B Clark UMCES HPL, Nov 2018

% load in data
nc_RR=netcdf(input('What is the full path and file name to the +M netcdf output file? -->','s'));
nc_NM=netcdf(input('What is the full path and file name to the -M netcdf output file? -->','s'));
% load('/Users/bclark/Desktop/RhodeFVCOM_2005/Data/rhoderiver_grid_v1_2.mat');
 load(input('What is the full path and file name to the *.mat grid file? ---> ','s'));
% load('/Users/bclark/Desktop/RhodeFVCOM_2005/Data/Transect_IDs.mat');
 load(input('What is the full path and file name to the *.mat transect ID file? ---> ','s'));

%pull out the DOC for each model scenario
WC_CDOC1=nc_RR{'WC_CDOC1'}(:);WC_NCDOC1=nc_RR{'WC_NCDOC1'}(:);WC_CDOC2=nc_RR{'WC_CDOC2'}(:);WC_NCDOC2=nc_RR{'WC_NCDOC2'}(:);WC_CDOC3=nc_RR{'WC_CDOC3'}(:);WC_NCDOC3=nc_RR{'WC_NCDOC3'}(:);
WC_CDOC1_NM=nc_NM{'WC_CDOC1'}(:);WC_NCDOC1_NM=nc_NM{'WC_NCDOC1'}(:);WC_CDOC2_NM=nc_NM{'WC_CDOC2'}(:);WC_NCDOC2_NM=nc_NM{'WC_NCDOC2'}(:);WC_CDOC3_NM=nc_NM{'WC_CDOC3'}(:);WC_NCDOC3_NM=nc_NM{'WC_NCDOC3'}(:);
tDOC_RR=WC_CDOC1+WC_NCDOC1+WC_CDOC2+WC_NCDOC2+WC_CDOC3+WC_NCDOC3;
tDOC_NM=WC_CDOC1_NM+WC_NCDOC1_NM+WC_CDOC2_NM+WC_NCDOC2_NM+WC_CDOC3_NM+WC_NCDOC3_NM;

%%
%plot the average difference in DOC conc and the transect out of the RR
subplot(4,2,1);
tricolor(e2n(:,2:4),lld_n(:,2),lld_n(:,3),mean(mean(tDOC_RR,2),1)-mean(mean(tDOC_NM,2),1));
hold on
plot(lld_n(Idx_unique,2),lld_n(Idx_unique,3),'k*');
axis([-76.565 -76.492 38.853 38.905]);
xplot = [0 (cumsum(sqrt(diff(xyd_n(Idx_unique',2)).^2 +diff(xyd_n(Idx_unique,3)).^2)))'];

dX=diff(xplot);
avg_NCDOC1_RR=squeeze(mean(mean(WC_NCDOC1,2),1));
avg_NCDOC2_RR=squeeze(mean(mean(WC_NCDOC2,2),1));
avg_NCDOC3_RR=squeeze(mean(mean(WC_NCDOC3,2),1));
avg_CDOC1_RR=squeeze(mean(mean(WC_CDOC1,2),1));
avg_CDOC2_RR=squeeze(mean(mean(WC_CDOC2,2),1));
avg_CDOC3_RR=squeeze(mean(mean(WC_CDOC3,2),1));
avg_tDOC_RR=squeeze(mean(mean(tDOC_RR,2),1));
avg_tDOC_NM=squeeze(mean(mean(tDOC_NM,2),1));

avg_NCDOC1_NM=squeeze(mean(mean(WC_NCDOC1_NM,2),1));
avg_NCDOC2_NM=squeeze(mean(mean(WC_NCDOC2_NM,2),1));
avg_NCDOC3_NM=squeeze(mean(mean(WC_NCDOC3_NM,2),1));
avg_CDOC1_NM=squeeze(mean(mean(WC_CDOC1_NM,2),1));
avg_CDOC2_NM=squeeze(mean(mean(WC_CDOC2_NM,2),1));
avg_CDOC3_NM=squeeze(mean(mean(WC_CDOC3_NM,2),1));
avg_tDOC_RR=squeeze(mean(mean(tDOC_RR,2),1));
avg_tDOC_NM=squeeze(mean(mean(tDOC_NM,2),1));

%now calculate gradients, in units of g C m^-3 km^1
dDOC_dx_RR=diff(avg_tDOC_RR(Idx_unique)')./dX*-1000;
dDOC_dx_NM=diff(avg_tDOC_NM(Idx_unique)')./dX*-1000;

dCDOC1_dx_NM=diff(avg_CDOC1_NM(Idx_unique)')./dX*-1000;
dCDOC2_dx_NM=diff(avg_CDOC2_NM(Idx_unique)')./dX*-1000;
dCDOC3_dx_NM=diff(avg_CDOC3_NM(Idx_unique)')./dX*-1000;
dNCDOC3_dx_NM=diff(avg_NCDOC3_NM(Idx_unique)')./dX*-1000;
dNCDOC2_dx_NM=diff(avg_NCDOC2_NM(Idx_unique)')./dX*-1000;
dNCDOC1_dx_NM=diff(avg_NCDOC1_NM(Idx_unique)')./dX*-1000;

dCDOC1_dx_RR=diff(avg_CDOC1_RR(Idx_unique)')./dX*-1000;
dCDOC2_dx_RR=diff(avg_CDOC2_RR(Idx_unique)')./dX*-1000;
dCDOC3_dx_RR=diff(avg_CDOC3_RR(Idx_unique)')./dX*-1000;
dNCDOC3_dx_RR=diff(avg_NCDOC3_RR(Idx_unique)')./dX*-1000;
dNCDOC2_dx_RR=diff(avg_NCDOC2_RR(Idx_unique)')./dX*-1000;
dNCDOC1_dx_RR=diff(avg_NCDOC1_RR(Idx_unique)')./dX*-1000;

%%
subplot(4,2,2)
plot(cumsum(dX)/1000,dDOC_dx_RR)
hold on
plot(cumsum(dX)/1000,dDOC_dx_NM)
area(cumsum(dX)/1000,dDOC_dx_RR-dDOC_dx_NM,'linestyle','none')

subplot(4,2,3)
plot(cumsum(dX)/1000,dCDOC1_dx_RR)
hold on
plot(cumsum(dX)/1000,dCDOC1_dx_NM)
area(cumsum(dX)/1000,dCDOC1_dx_RR-dCDOC1_dx_NM,'linestyle','none')

subplot(4,2,4)
plot(cumsum(dX)/1000,dNCDOC1_dx_RR)
hold on
plot(cumsum(dX)/1000,dNCDOC1_dx_NM)
area(cumsum(dX)/1000,dNCDOC1_dx_RR-dNCDOC1_dx_NM,'linestyle','none')

subplot(4,2,5)
plot(cumsum(dX)/1000,dCDOC2_dx_RR)
hold on
plot(cumsum(dX)/1000,dCDOC2_dx_NM)
area(cumsum(dX)/1000,dCDOC2_dx_RR-dCDOC2_dx_NM,'linestyle','none')

subplot(4,2,6)
plot(cumsum(dX)/1000,dNCDOC2_dx_RR)
hold on
plot(cumsum(dX)/1000,dNCDOC2_dx_NM)
area(cumsum(dX)/1000,dNCDOC2_dx_RR-dNCDOC2_dx_NM,'linestyle','none')


subplot(4,2,7)
plot(cumsum(dX)/1000,dCDOC3_dx_RR)
hold on
plot(cumsum(dX)/1000,dCDOC3_dx_NM)
area(cumsum(dX)/1000,dCDOC3_dx_RR-dCDOC3_dx_NM,'linestyle','none')

subplot(4,2,8)
plot(cumsum(dX)/1000,dNCDOC3_dx_RR)
hold on
plot(cumsum(dX)/1000,dNCDOC3_dx_NM)
area(cumsum(dX)/1000,dNCDOC3_dx_RR-dNCDOC3_dx_NM,'linestyle','none')

saveas(gcf,'DOC_gradient_8panel.fig');
print -painters -depsc DOC_gradient_8panel.eps

disp('Finished!');




