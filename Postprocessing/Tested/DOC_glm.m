% script to make a model of the percent gross DOC production

% and total DOC stock as  function of estuary volume to marsh area

%in the Rhode River with varying sections of the tributary

% B Clark, UMCES, January 2019

%first run the Calculate_DOC_stock.m script to get the total DOC stocks;


%get the percentages from the Carbon_Budget_ncdf.m scripts for the marsh
%contribution at varying sections

%here they are stored in a data structure

%load(*/RhodeFVCOM_2005/Data/EVMA_data.mat);
load(input('What is the full path and file name to thee *.mat file with the EVMA data? ---> ','s'));
%load(*/RhodeFVCOM_2005/Data/RhodeRiver_nodes.mat')
load(input('What is the full path and file name to Tributary Node Information? --> ','s'));
load(input('What is the full path and file name to the RR Polygon file? ---> ','s'));

estuary.nodes=RR_nodes_out(find(lld_n(RR_nodes_out,3) > 38.875));

%get a vector the stock at decreasing EVMA
mystock(1)=(median(sum(DOC_stock_z_RR(:,estuary.nodes),2))-median(sum(DOC_stock_z_NM(:,estuary.nodes),2)))...
    ./median(sum(DOC_stock_z_RR(:,estuary.nodes),2))*100;
mystock(2)=(median(sum(DOC_stock_z_RR(:,poly1_nodes),2))-median(sum(DOC_stock_z_NM(:,poly1_nodes),2)))...
    ./median(sum(DOC_stock_z_RR(:,poly1_nodes),2))*100;
mystock(3)=(median(sum(DOC_stock_z_RR(:,poly2_nodes),2))-median(sum(DOC_stock_z_NM(:,poly2_nodes),2)))...
    ./median(sum(DOC_stock_z_RR(:,poly2_nodes),2))*100;
mystock(4)=(median(sum(DOC_stock_z_RR(:,poly3_nodes),2))-median(sum(DOC_stock_z_NM(:,poly3_nodes),2)))...
    ./median(sum(DOC_stock_z_RR(:,poly3_nodes),2))*100;
mystock(5)=(median(sum(DOC_stock_z_RR(:,poly4_nodes),2))-median(sum(DOC_stock_z_NM(:,poly4_nodes),2)))...
    ./median(sum(DOC_stock_z_RR(:,poly4_nodes),2))*100;

%make a logarithmic model for GDP and fit percentages
mymod_gdp=fitglm(EVMA,mypercent_GDP,'link','log');
gdp_predict=predict(mymod_gdp,(0:100)');
%find the prediction for entire Chesapeake Bay EVMA (59.6 m)
CBAY_gdp=predict(mymod_gdp,59.6);

%now do the same but for the stock

mymod_stock=fitglm(EVMA,mystock,'link','log');
stock_predict=predict(mymod_stock,(0:100)');
CBAY_stock=predict(mymod_stock,59.6);

%now plot them on the same graph
figure;
plot(EVMA,mypercent_GDP,'kd','markerfacecolor','k');
hold on
plot(0:100,gdp_predict,'k');
plot(EVMA,mystock,'d','markerfacecolor',[1.0 0.5 .25],'markeredgecolor','k');
plot(0:100,stock_predict,'color',[1.0 0.5 .25]);
plot(repmat(59.6,1,21),0:20,'k-.')
xlabel('EVMA (m)');ylabel('% Contribution');
saveas(gcf,'GLM_RhodeRiver.fig');
saveas(gcf,'GLM_RhodeRiver.eps');

disp('Finished!');















