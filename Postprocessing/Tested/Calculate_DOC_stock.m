%Script to integrate the DOC in the water column

% of a body of water defined by the nodes

% and then calculate the average integrated DOC stock

% B Clark, Jan 31, 2019

%first pull out all of the DOC 
% load('*/RhodeFVCOM_2005/Data/RhodeRiver_nodes.mat');
load(input('What is the full path and file name to Tributary Node Information? --> ','s'));
nc_RR=netcdf(input('What is the full path and file name to the +M netcdf output file? --->  ','s'));
nc_NM=netcdf(input('What is the full path and file name to the -M netcdf output file? --->  ','s'))
% node_file=('*/RhodeFVCOM_2005/Data/node_CVs_area.dat');
node_file=input('What is the full path and file name to the Node Area file? --->  ','s');
[node_area]=importdata(node_file);
% node_file=('*/RhodeFVCOM_2005/Data/rhoderiver_grid_v1_2.mat');
load(input('what is the full path an file name to the *.mat grid definition file? --> ','s'));

estuary.nodes=RR_nodes_out(find(lld_n(RR_nodes_out,3) > 38.875));
disp('Getting Depth and calculating volume')

%Figure out the volume of each nodal area and each layer
mydepth=nc_RR{'depth'}(:);
s_hybrid=nc_RR{'siglev'}(:);
dz=diff(s_hybrid,1)*-1;
[nlayers,ngrids]=size(dz);
[ntimes,~,~]=size(mydepth);

layer_thickness=permute(permute(repmat(dz,1,1,ntimes),[2 3 1]).*repmat(mydepth',1,1,10),[2 3 1]);
myvolume=layer_thickness.*permute(repmat(repmat(node_area(:,2),1,ntimes),1,1,10),[2 3 1]);

disp('Getting DOC data for +M');
CDOC1_RR=nc_RR{'WC_CDOC1'}(:);
CDOC2_RR=nc_RR{'WC_CDOC2'}(:);
CDOC3_RR=nc_RR{'WC_CDOC3'}(:);
NCDOC1_RR=nc_RR{'WC_NCDOC1'}(:);
NCDOC2_RR=nc_RR{'WC_NCDOC2'}(:);
NCDOC3_RR=nc_RR{'WC_NCDOC3'}(:);

totalDOC_RR=CDOC1_RR+CDOC2_RR+CDOC3_RR+NCDOC1_RR+NCDOC2_RR+NCDOC3_RR;

disp('Getting DOC data for -M');
CDOC1_NM=nc_NM{'WC_CDOC1'}(:);
CDOC2_NM=nc_NM{'WC_CDOC2'}(:);
CDOC3_NM=nc_NM{'WC_CDOC3'}(:);
NCDOC1_NM=nc_NM{'WC_NCDOC1'}(:);
NCDOC2_NM=nc_NM{'WC_NCDOC2'}(:);
NCDOC3_NM=nc_NM{'WC_NCDOC3'}(:);

totalDOC_NM=CDOC1_NM+CDOC2_NM+CDOC3_NM+NCDOC1_NM+NCDOC2_NM+NCDOC3_NM;
%%
disp('Got DOC data, now integrating');

DOC_stock_RR=myvolume.*totalDOC_RR;
DOC_stock_z_RR=squeeze(sum(DOC_stock_RR,2));
DOC_stock_RR_out=sum(DOC_stock_z_RR(:,estuary.nodes),2);
mean_DOC_stock_RR=mean(DOC_stock_RR_out)./1E6
median_DOC_stock_RR=median(DOC_stock_RR_out)./1E6

DOC_stock_NM=myvolume.*totalDOC_NM;
DOC_stock_z_NM=squeeze(sum(DOC_stock_NM,2));
DOC_stock_NM_out=sum(DOC_stock_z_NM(:,estuary.nodes),2);
mean_DOC_stock_NM=mean(DOC_stock_NM_out)./1E6
median_DOC_stock_NM=median(DOC_stock_NM_out)./1E6

mean_volume=mean(sum(sum(myvolume(:,:,estuary.nodes),3),2),1)
save('DOC_stock_work.mat','DOC_stock_NM_out','DOC_stock_RR_out','DOC_stock_z_NM','DOC_stock_z_RR');

disp('Finished!');
























