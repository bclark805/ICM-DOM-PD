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
             Nvar_read   =input('please give number of state variables in history output (e.g. 84 for budget): '); 
             Nvar        =input('please give number of state variables in history output to use (e.g. 84 for budget): ');
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

            varnames={  'H'; 'zeta';  'D'; 'DL'};%1-4
                      
                  
           varunits={'(m)';  '(m)';'(m)';'(m)'}             
% sediment fluxes
                      varnames = [varnames; 'JNH4';'JNO3';...
                                            'JPOC1';'JPOC2';'JPOC3';...
                                            'JPON1';'JPON2';'JPON3';...
                                            'JCDOC1' ;'JNCDOC1' ;...
                                            'JCDOC2' ;'JNCDOC2' ;...
                                            'JCDOC3' ;'JNCDOC3' ;...
                                            'JCDON1' ;'JNCDON1' ;...
                                            'JCDON2' ;'JNCDON2' ;...
                                            'JCDON3' ;'JNCDON3' ] %5-24
                                            
                                        
                      varunits = [varunits ; 'kg N';'kg N';...
                                             'kg C';'kg C';'kg C';...
                                             'kg N';'kg N';'kg N';...
                                             'kg C';'kg C';...
                                             'kg C';'kg C';...
                                             'kg C';'kg C';...
                                             'kg N';'kg N';...
                                             'kg N';'kg N';...
                                             'kg N';'kg N']
                                         
              varnames=[varnames; 'AlgDON';'AlgPON';'AlgNH4';...%27
                                    'AlgNO3';'DenitNO3';'NitrifNH4';...%30
                                  'DenitDON';'HDRLPON';'HDRRPON';...%33
                                  'COAGN';'MNLDON1';'MNLDON2';...%36
                                  'MNLDON3';'LPON';'RPON';...%39
                                  'DTWCDON1';'DTWCDON2';'DTWCDON3'; ...%42
                                  'DTWNCDON1';'DTWNCDON2';'DTWNCDON3'; ...$45
                                  'DTNH4';'DTNO3';...%47
                                  'DTLPON';'DTRPON']; %49
                              
               varunits=[varunits; 'kg N';'kg N';'kg N';...
                                   'kg N';'kg N';'kg N';...
                                   'kg N';'kg N';'kg N';...
                                   'kg N';'kg N';'kg N';...
                                   'kg N';'(g N m^-^-3)';'(g N m^-^-3)';...
                                   'kg N';'kg N';'kg N';...
                                   'kg N';'kg N';'kg N';...
                                   'kg N';'kg N';'kg N';'kg N']
                               
            % Carbon
              varnames=[varnames; 'AlgDOC';'AlgPOC';'DenitDOC';...%52
                                  'HDRLPOC';'HDRRPOC';'COAGC';...%55
                                  'MNLDOC1';'MNLDOC2';'MNLDOC3';...%58
                                  'DTWCDOC1';'DTWCDOC2';'DTWCDOC3';...^61
                                  'DTWNCDOC1';'DTWNCDOC2';'DTWNCDOC3';...%64
                                   'DTLPOC';'DTRPOC'];%66
                              
              varunits = [varunits ; 'kg C';'kg C';'kg C';...
                                     'kg C';'kg C';'kg C';...
                                     'kg C';'kg C';'kg C';...
                                     'kg C';'kg C';'kg C';...
                                     'kg C';'kg C';'kg C';...
                                     'kg C';'kg C'];
            %photochemistry          
                      varnames = [varnames; 'DPD32';'DPD31';'DPD30';'DPD3N';...
                                            'DPD21';'DPD20';'DPD2N'; 'DPD10';'DPD1N'];%75
                                        
                      varunits = [varunits ;'kg C';'kg C';'kg C';...
                                             'kg C';'kg C';'kg C';...
                                             'kg C';'kg C';'kg C'];
                                         
                       varnames= [varnames ; 'DPD32N';'DPD31N';'DPD30N';'DPD3NN';...
                                            'DPD21N';'DPD20N';'DPD2NN'; 'DPD10N';'DPD1NN' ];%84
                                        
                       varunits= [varunits; 'kg N';'kg N';'kg N';...
                                            'kg N';'kg N';'kg N';...
                                            'kg N';'kg N';'kg N'];
                                         

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
