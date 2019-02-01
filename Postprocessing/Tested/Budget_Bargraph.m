%make DOC bargraph

%Quick script to make a bargraph showiong each reactivity class

%and the total DOC using the output from the Carbon_Budget.m script

% B Clark UMCES, Jan 2019


load(input('What is the full path and file name to the *budget*.mat file with the budget structures? --->  ','s'));
DOC1_out=myvals_out.DOC1;
DOC2_out=[];
DOC2_out(1:4)=myvals_out.DOC2(1:4);
DOC2_out(6:12)=myvals_out.DOC2(5:end);
DOC3_out(1:4)=myvals_out.DOC3(1:4);
DOC3_out(5)=0.0;
DOC3_out(6:8)=myvals_out.DOC3(5:7);
DOC3_out(9:12)=myvals_out.DOC3(9:12);;


bar((DOC1_out+DOC2_out+DOC3_out)./1000);% plot bargraph in Tons C
hold on
bar([DOC1_out' DOC2_out' DOC3_out']./1000);
legend('Total DOC','DOC1','DOC2','DOC3');
xticklabels(mynames_out.DOC1);ylabel('DOC (tons C)');
xtickangle(45)