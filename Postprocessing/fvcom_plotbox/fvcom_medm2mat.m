function dat = fvcom_medm2mat(fname,pathname)
%
% dat = fvcom_medm2mat(fname,pathname)
%    
%    Function to read an medm file of FVCOM outputs
%
% Example:
%    dat = fvcom_medm2mat('psm_sim0001.dat','/home/long075/mypsfvm/trunk/example/channel/model_output/medm/')
%
% Code & Documentation
%   
%   Wen Long, Seattle, 05/15/2013
%

dat = [];   % initialize output
opt = 2;    % file loading option

if opt == 1 % slow - initial attempt
    disp(['Option 1: ' fname])
    fid = fopen([pathname fname],'r');

    dat.header          = [];
    datlen              = fread(fid,1,'int32')/4;
    dat.header(1,1:3)   = fread(fid,datlen-1,'int32');
    dat.header(1,4)     = fread(fid,1,'float32');
    fread(fid,1,'int32')/4;

    dat.velocity = [];
    for i = 1:dat.header(1,2)
        datlen                   = fread(fid,1,'int32')/4;
        dat.velocity(i,1:datlen) = fread(fid,datlen,'float32');
        fread(fid,1,'int32')/4;
    end

    dat.hydro    = [];
    for i = 1:dat.header(1,3)
        datlen                = fread(fid,1,'int32')/4;
        dat.hydro(i,1:datlen) = fread(fid,datlen,'float32');
        fread(fid,1,'int32')/4;
    end

    fclose(fid);
elseif opt == 2 % optimized
    disp(['Option 2: ' fname])
    fid = fopen([pathname fname],'r');
        dat_i = fread(fid,'int32');
    fclose(fid);
    fid = fopen([pathname fname],'r');
        dat_f = fread(fid,'float32');
    fclose(fid);
    
    dat.header   = [dat_i(2:4)' dat_f(5)];
    
    numvars      = dat_i(7)/4;
    numvars      = numvars+2;
    dat_i        = dat_i(7:end);
    dat_f        = dat_f(7:end);
    numdat       = dat.header(1,2)*numvars;
    dat_idx      = setdiff(setdiff(1:numdat,1:numvars:numdat),numvars:numvars:numdat);
    dat.velocity = dat_f(dat_idx);
    dat.velocity = reshape(dat.velocity,[numvars-2 dat.header(1,2)])'; 
    
    numvars      = dat_i(numdat+1)/4;
    numvars      = numvars+2;
    dat_i        = dat_i(numdat+1:end);
    dat_f        = dat_f(numdat+1:end);
    numdat       = dat.header(1,3)*numvars;
    dat_idx      = setdiff(setdiff(1:numdat,1:numvars:numdat),numvars:numvars:numdat);
    dat.hydro    = dat_f(dat_idx);
    dat.hydro    = reshape(dat.hydro,[numvars-2 dat.header(1,3)])';
end
