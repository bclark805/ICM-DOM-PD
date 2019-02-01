         %wqmhis2matlab.m
         %
         %matlab tool to read the history outputs from wqm model and then reseave into matlab format for later use
         %calculate statistics, plot horizontal maps of chosen levels or vertical transects 
         %
         %
         %
         % Wen Long @ PNNL, April 19th, 2012
         %
         % updated by B Clark, UMCES, Jan 2019
         
%          startup

           psm_mat_grid_file = input('What is the full path and file name of the *.mat grid file? -----> ','s');
               load(psm_mat_grid_file);   
               x_n=xyd_n(:,2);                  %x coordinate of nodes
               y_n=xyd_n(:,3);                  %y coordinate of nodes
               lon_n=lld_n(:,2);                %longitude of nodes
               lat_n=lld_n(:,3);                %latitude of nodes
               h_n=lld_n(:,4);                  %depth of nodes

	  
             NN=size(xyd_n,1) ;    %number of nodes
             NE=size(xy_e,1) ;     %number of elements

             
             his_dir=input('What is the full path to the directory with the ICM history.out output files? ---> ','s');
% 
             if(~exist(his_dir,'dir'))
               error(['Oops, history output file dir:' his_dir ' does not exist!']);
             end

           %obtain history output file
             hfile_prefix=input('Please give the history file name prefix (e.g. rhd_history):','s');
             hfile_ext   =input('Please give the history file name extension (e.g. .out):','s');
             his_Nz      =input('please give number of vertical layers in history output (e.g. 10):');
             budget_flag = input('Is this output with or without the budget variables, 0 for no 1 for yes ---> ');
             Nvar_read   =input('please give number of state variables in history output (e.g. 92 for budget, 32 for no budget): '); 
             Nvar        =input('please give number of state variables in history output to use (e.g. 92 for budget, 32 for no budget): ');
                 %number of variables in history file to use in history file

%------------------------------------------------------------------------------------------------------------------%
%
%                             WRITE(UNIT_HIS+K,'(:F8.4,1x,I1,1x,I8/(:8(F12.6,1x),F12.6))')JDAY,4,MGL, &
%				(H_GL(I),I=1,MGL), & !Wen Long output still water depth
%				(EL_GL(I),I=1,MGL), & !Wen Long output surface elevation
%				(D_GL(I),I=1,MGL), & !Wen Long output total depth
%				(-D_GL(I)*ZZ(K),I=1,MGL), & !Wen Long,depth of this layer (from surface)
%				(C2_GL(I,K,27),I=1,MGL), & !DO
%				(C2_GL(I,K,9),I=1,MGL), & !LDOC
%				(C2_GL(I,K,4),I=1,MGL), & !ALG1
%				(C2_GL(I,K,5),I=1,MGL), & !ALG2
%				(C2_GL(I,K,13),I=1,MGL), & !NH4
%				(C2_GL(I,K,14),I=1,MGL), & !NO3
%				(C2_GL(I,K,20),I=1,MGL), & !PO4
%				(T_GL(I,K),I=1,MGL), & !Temperature
%				(S_GL(I,K),I=1,MGL) !Salinity 
%
%-------------------------------------------------------------------------------------------------------------------%

            varnames={  'H'; 'zeta';  'D'; 'DL'; 'DO'; 'NPP';   'ALG1';    'ALG2';  'NH4';     'NO3';     'PO4';  ...
                      'T';    'S';'PAR';'LPOC';'RPOC';...
                      'CDOC1';  'NCDOC1';'CDOC2';  'NCDOC2';'CDOC3';  'NCDOC3';...
                      'CDON1';'NCDON1';'CDON2';'NCDON2';'CDON3';'NCDON3';...
                      'KD';'ISS';'Density';'Tau';};
                      
                  
           varunits={'(m)';  '(m)';'(m)';'(m)';'(mgDO L^-^1)';'(gC m^-^-2 d^-^1)';'(gC m^-^-3)'; '(gC m^-^-3)';'(gN m^-^-3)';'(gN m^-^-3)';'(gP m^-^-3)';  ...
                      '(\circC)';'(psu)';'(E m^-^-2 day^-^1)';'(gC m^-^-3)';'(gC/m^3)';...
                      '(gC m^-^-3)';'(gC m^-^-3)';'(gC m^-^-3)';'(gC m^-^-3)';'(gC m^-^-3)';'(gC m^-^-3)';...
                      '(gN m^-^-3)';'(gN m^-^-3)';'(gN m^-^-3)';'(gN m^-^-3)';'(gN m^-^-3)';'(gN m^-^-3)';...
                      '(m^-^1)'; '(g ISS m^-^3)';'(kg m^-^3)';'(pa)'};
     if(budget_flag)             
% sediment fluxes
                      varnames = [varnames; 'JNH4';'JNO3';'JPOC';'JPON';...
                                            'JDOC1' ;'JDOC2' ;'JDOC3' ;...
                                            'JDON1' ;'JDON2' ;'JDON3' ];
                                        
                      varunits = [varunits ; '(mg N m^-^-2 day^-^1)';'(mg N m^-^-2 day^-^1)';'(mg C m^-^-2 day^-^1)';'(mg N m^-^-2 day^-^1)';....
                                            '(g C m^-^-2 day^-^1)';'(g C m^-^-2 day^-^1)';'(g C m^-^-2 day^-^1)';...
                                            '(g N m^-^-2 day^-^1)';'(g N m^-^-2 day^-^1)';'(g N m^-^-2 day^-^1)'];
          %NITROGEN
            
              varnames=[varnames; 'AlgDON';'AlgPON';'AlgNH4';...
                                    'AlgNO3';'DenitNO3';'NitrifNH4';...
                                  'DenitDON';'HDRLPON';'HDRRPON';...
                                  'COAGN';'MNLDON1';'MNLDON2';...
                                  'MNLDON3';'LPON';'RPON';'DTWDON1'; ...
                                  'DTNH4';'DTNO3';'DTLPON';'DTRPON'];
                              
               varunits=[varunits; '(g N m^-^-3 day^-^1)';'(g N m^-^-3 day^-^1)';'(g N m^-^-3 day^-^1)';...
                                   '(g N m^-^-3 day^-^1)';'(g N m^-^-3 day^-^1)';'(g N m^-^-3 day^-^1)';...
                                   '(g N m^-^-3 day^-^1)';'(g N m^-^-3 day^-^1)';'(g N m^-^-3 day^-^1)';...
                                   '(g N m^-^-3 day^-^1)';'(g N m^-^-3 day^-^1)';'(g N m^-^-3 day^-^1)';...
                                   '(g N m^-^-3 day^-^1)';'(g N m^-^-3 day^-^1)';'(g N m^-^-3 day^-^1)';'(g N m^-^-3 day^-^1)'...
                                   ;'(g N m^-^-3 day^-^1)';'(g N m^-^-3 day^-^1)';'(g N m^-^-3 day^-^1)';'(g N m^-^-3 day^-^1)'] ;                
            % Carbon
              varnames=[varnames; 'AlgDOC';'AlgPOC';'DenitDOC';...
                                  'HDRLPOC';'HDRRPOC';'COAGC';...
                                  'MNLDOC1';'MNLDOC2';'MNLDOC3';'DTWDOC1'...
                                   ;'DTLPOC';'DTRPOC'];
                              
              varunits = [varunits ; '(g C m^-3 d^-1)';'(g C m^-3 d^-1)';'(g C m^-3 d^-1)';...
                                     '(g C m^-3 d^-1)';'(g C m^-3 d^-1)';'(g C m^-3 d^-1)';...
                                     '(g C m^-3 d^-1)';'(g C m^-3 d^-1)';'(g C m^-3 d^-1)';'(g C m^-3 d^-1)'...
                                     ;'(g C m^-3 d^-1)';'(g C m^-3 d^-1)'];
            %photochemistry          
                      varnames = [varnames; 'DPD32';'DPD31';'DPD30';'DPD3N';...
                                            'DPD21';'DPD20';'DPD2N'; 'DPD10';'DPD1N']
                      varunits = [varunits ; '(g C m^-3 d^-1)';'(g C m^-3 d^-1)';'(g C m^-3 d^-1)';...
                                             '(g C m^-3 d^-1)';'(g C m^-3 d^-1)';'(g C m^-3 d^-1)';...
                                             '(g C m^-3 d^-1)';'(g C m^-3 d^-1)';'(g C m^-3 d^-1)'];
                                         
                       varnames= [varnames ; 'DPD32N';'DPD31N';'DPD30N';'DPD3NN';...
                                            'DPD21N';'DPD20N';'DPD2NN'; 'DPD10N';'DPD1NN' ];
                                        
                       varunits= [varunits; '(g N m^-3 d^-1)';'(g N m^-3 d^-1)';'(g N m^-3 d^-1)';...
                                             '(g N m^-3 d^-1)';'(g N m^-3 d^-1)';'(g N m^-3 d^-1)';...
                                             '(g N m^-3 d^-1)';'(g N m^-3 d^-1)';'(g N m^-3 d^-1)'];
                                         
     end
%%
          C2CHLA_ALG1=50.0;    %Carbon to CHLA ratio for ALG1 (mgC/mgCHLA) 
          C2CHLA_ALG2=50.0;    %Carbon to CHLA ratio for ALG2 (mgC/mgCHLA)

           %check if files exist
           for k=1:his_Nz
               hfile=[his_dir '/' hfile_prefix '_' num2str(k,'%05d') hfile_ext];
           %     hfile = [his_dir '/' hfile_prefix hfile_ext];
               if(~exist(hfile,'file'))
                 error(['Oops, cannot find file ' hfile ]);
               end
           end

           %read history output file 
           if(strcmp(hfile_ext,'.mat'))
              read_his=false;
           else
              read_his=true;    %only read history if extension name is not '.mat'
           end
 
           if(read_his)
            for k=1:his_Nz
                  hfile=[his_dir '/' hfile_prefix '_' num2str(k,'%05d') hfile_ext];
                  %hfile = [his_dir '/' hfile_prefix hfile_ext];
                  if(~exist(hfile,'file'))
                    error(['Oops, cannot find file ' hfile ]);
                  end

                  %open file
                  hf_his=fopen(hfile,'r');

                  more_to_read=true;
                  it=0;                  %number of records in time
                  time=[];
                  while(more_to_read)    %keep reading
                        jday_tmp=[];
                        %read time, 4, and number of nodes
                        jday_tmp=fscanf(hf_his,'%f',1);    %read float  (time in days)
                         jday_tmp
                         tmp4_tmp=fscanf(hf_his,'%d',1);    %read integer (4)
                         tmp4_tmp
                         NN_tmp  =fscanf(hf_his,'%d',1);    %read integer (number of nodes)
                         NN_tmp
                        if(size(jday_tmp,1)==0)
                           more_to_read=false;
                           break                       %no more to read, break out of the while loop
                        end

                        it=it+1; 
                        time(it)=jday_tmp;                 %record time
                        display(['reading time step : ',num2str(it), ' of file ', hfile, ' ...']); 
                        %read variables
                        for ivar=1:Nvar_read  
                          
                            for i=1:NN                     %read one variable (all nodes)         
                                tmp_var(i)=fscanf(hf_his,'%e',1);
                          
                            end
			                if(ivar<=Nvar)
                                hisdata(it,:,ivar)=tmp_var;     %put the results in hisdata for this variable
                            end
                        end
                 end
                 NT=it                                     %total number of time records 
                  
                 %close file
                 fclose(hf_his);

                 display(['successfully read ' num2str(NT) ' records from file: ' hfile]);

                 %record this k'th level in a big matlab file

                 tmp_matfile=['hisdata_' num2str(k,'%05d') '.mat'];
                
                 save(tmp_matfile,'hisdata', ...  %history data
                                     'time', ...  %time 
                                       'NT', ...  %number of records
                                 'varnames', ...  %variable names
                                 'varunits', ...  %units of variables
                                      'x_n', ...  %node x coordinate
                                      'y_n', ...  %node y coordinate
                                    'lon_n', ...  %node longitude
                                    'lat_n', ...  %node latitude
                                      'h_n', ...  %depth (still water)
                                     'xy_e', ...  %x,y coordinates of elements
                                     'e2n', ...  %connectivity matrix (element to node)
                                     '-v7.3',... %bigger than 2G file
                       '-mat'); 

                display(['successfully saved ' num2str(NT) ' records of data into file: ' tmp_matfile]);

           end  %end of k loop

         end   %end of read_his if-block
