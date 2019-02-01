% script to calculate oxygen variability in the Rhode River

% B Clark, UMCES, January 2019

%load in the Rhode River nodes and oxygen variables
load(input('What is the full path and file name to Tributary Node Information? --> ','s'));
nc_RR=netcdf(input('What is the full path and file name of the +M netcdf output file ? --> ','s'));
nc_NM=netcdf(input('What is the full path and file name of the -M netcdf output file ? --> ','s'));
marsh_file=input('What is the full path and file name to the sediment definition file? --> ','s');
marsh_data=importdata(marsh_file);

%get the nodes that are in the estuary ONLY
estuary.nodes=RR_nodes_out(find(lld_n(RR_nodes_out,3) > 38.875));
marsh.nodes= find(marsh_data(:,3)>0);
estuary.nodesonly=setdiff(estuary.nodes,marsh.nodes);

%get the oxygen

DOXG_RR=nc{'DOXG'}(:);
DOXG_NM=nc{'DOXG'}(:);

%calculate the average minimum bottom water DO throughout the river
%minimum is in the time dimension
%average is in the space dimension
average_Min_DOXG_RR=mean(squeeze(min(DOXG_RR(:,10,estuary.nodesonly))));
average_Min_DOXG_NM=mean(squeeze(min(DOXG_NM(:,10,estuary.nodesonly))));

%now find the amount of hours each location has spent in hypoxia

for i = 1 : length(estuary.nodesonly);
    hypoxic_hours_RR(i)=length(find(DOXG_RR(:,10,estuary.nodesonly(i))<2.0));
    hypoxic_hours_NM(i)=length(find(DOXG_NM(:,10,estuary.nodesonly(i))<2.0));

end
%average time in hypoxia across the river
average_hypoxic_time_RR=mean(hypoxic_hours_RR);
average_hypoxic_time_NM=mean(hypoxic_hours_NM);

%how many places experience hypoxia
hypoxia_here=find(hypoxic_hours_RR>0);
%percentage of places that experienced hypoxia
percent_hypoxia=length(hypoxia_here)/length(estuary.nodesonly);







