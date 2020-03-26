% SCript to write the initial condition file for
% Restart of ICM
%
%B Clark Oct 2015
%
% wqm_hisdir = input('What is full path to the WQM history directory? --->  ','s');
% nlayers = input('How many layers are in the model? --> ');
% dom_switch = input('Is this the old DOM formula or new DOM formula?  enter 0 for old and 1 for new ---> ')
 load(input('What is the full path and file name to the *.mat grid file? ---> ','s'));
% load(matfilein);
%
% mytime_in = input('What time do you want to write the output file for?  --> ');
% mytime = find(time == mytime_in);
%

nc=netcdf('netcdf/rhd_0001.nc');
Big_array_out(:,:,1)=nc{'temp'}(end,:,:);%T
Big_array_out(:,:,2)=nc{'salinity'}(end,:,:);%T
Big_array_out(:,:,3)=10;%SSI
Big_array_out(:,:,4)=nc{'B1'}(end,:,:); %ALG1
Big_array_out(:,:,5)=nc{'B2'}(end,:,:); %ALG2
Big_array_out(:,:,6)=0.0 ;%AL3
Big_array_out(:,:,7)=0.0 ; % SZ
Big_array_out(:,:,8)=0.0;  % LZ
Big_array_out(:,:,9)=nc{'WC_CDOC1'}(end,:,:); %CDOC1
Big_array_out(:,:,10)=nc{'WC_NCDOC1'}(end,:,:); %NCDOC1

Big_array_out(:,:,11)=nc{'LPOC'}(end,:,:); %LPOC
Big_array_out(:,:,12)=nc{'RPOC'}(end,:,:); %RPOC
Big_array_out(:,:,13)=nc{'NH4'}(end,:,:); %NH4
Big_array_out(:,:,14)=nc{'NO3'}(end,:,:); % %N03
Big_array_out(:,:,15)=0.0; %UREA

Big_array_out(:,:,16)=nc{'WC_CDON1'}(end,:,:); %CDON1
Big_array_out(:,:,17)=nc{'WC_NCDON1'}(end,:,:); %NCDON1

Big_array_out(:,:,18)=nc{'LPON'}(end,:,:); %LPOC; %LPON
Big_array_out(:,:,19)=nc{'RPON'}(end,:,:); %RPOC

Big_array_out(:,:,20)=0.001; %PO4

Big_array_out(:,:,21)=0.001; %CDOP1
Big_array_out(:,:,22)=0.001; %NCDOP1

Big_array_out(:,:,23)=0.01; %LPOP
Big_array_out(:,:,24)=0.01; %RPOP

Big_array_out(:,:,25)=0.0; %PIP
Big_array_out(:,:,26)=0.0; %SOD
Big_array_out(:,:,27)= 0.0;%nc{'SOD'}(end,:,:); %DO2

Big_array_out(:,:,28:32)=0.0; %rest are zeros

Big_array_out(:,:,33)=nc{'WC_CDOC2'}(end,:,:); %CDOC2
Big_array_out(:,:,34)=nc{'WC_NCDOC2'}(end,:,:); %NCDOC2
Big_array_out(:,:,35)=nc{'WC_CDOC3'}(end,:,:);  %CDOC3
Big_array_out(:,:,36)=nc{'WC_NCDOC3'}(end,:,:); %NCDOC3

Big_array_out(:,:,37)=nc{'WC_CDON2'}(end,:,:); %CDOC2
Big_array_out(:,:,38)=nc{'WC_NCDON2'}(end,:,:); %NCDOC2
Big_array_out(:,:,39)=nc{'WC_CDON3'}(end,:,:); %CDOC3
Big_array_out(:,:,40)=nc{'WC_NCDON3'}(end,:,:); %NCDOC3


Big_array_out(:,:,41)=0.001; %CDOC2
Big_array_out(:,:,42)=0.001; %NCDOC2
Big_array_out(:,:,43)=0.001; %CDOC3
Big_array_out(:,:,44)=0.001; %NCDOC3


Big_array_out(11,:,:)=Big_array_out(10,:,:);
[nlayers,ngrids]=size(Big_array_out(:,:,1));
%
%     end
% % end
%
%
% %   Big_array_out(:,nlayers+1,:)=Big_array_out(:,nlayers,:);
%
% end
%%
fid = fopen('RhodeRiver_initial_wq_vert.dat','w');

ngrids=length(lld_n);
nlayers=10;
nvars=44;
for igrid = 1:ngrids
    for jlay=1:nlayers+1;
        for kvars = 1:nvars
            fprintf(fid,'%16.8f',Big_array_out(jlay,igrid,kvars));
        end
        fprintf(fid,'\n');
    end
end


