function [xp,yp,xt,yt]=fvcom_transect_create(file_grid,transect_name,ds)
%
%function [xp,yp,xt,yt]=fvcom_transect_create(file_grid,transect_name,ds)
%
% Simple program to create a transect off a fvcom grid (matlab format)
%       by allowing the user to click on points on the intended transect
%
% Outputs:
%	xt,yt  ---x and y coordinates of the points that the user clicked on
%       	    the map
%       xp,yp  ---x and y coordinates of interpolated points on the chosen transct
%
% Inputs: 
%       file_grid      --- string storing file name of the fvcom grid (prepared by sms2dm2mat.m)
%       transect_name  --- string storing name of the transect to be created
%       ds             --- segment length used in the interpolation to get xp, yp from xt,yt
%                          (unit is in x, y units)
%
% Example: 
%            file_grid = '/home/long075/mypsfvm/trunk/example/pugetsound-skagit/forcing/grid/PS2_SPB.mat' ;
%        transect_name = 'A-A1'                                                                           ;
%                   ds = 10.0                                                                             ;
%        [xp,yp,xt,yt] = fvcom_transect_create(file_grid,transect_name,ds)                                ;
%        save(['trans_loc_' transect_name '.mat'],'xp','yp','xt','yt');  %save the results into a .mat file
%
% Dependency: tricolor, cumsum, spline, ginput
%
% Code and Documentation:
%
%   Wen Long @ PNNL, Seattle, WA, 98109
%   03/27/2013
%

  if(nargin <3)
    error('Oops, no input arguments given'); 
  end
  %check if file exist

  if(~exist(file_grid,'file'))
    error(['Oops, not able to find fvcom grid file:' file_grid]); 
    exit;
  end

  %check value of ds

  if(isempty(ds))
    error(['Oops, missing parameter ds']);
  else
    if(ds<=0)
     error('oops, ds argument must be a number greater than zero'); 
    end
  end

  %check transect name
  if(isempty(transect_name))
    error('oops, transect_name argument should not be empty'); 
  end
   

  %load file

   load(file_grid);

   x_n=xyd_n(:,2);  %x coord
   y_n=xyd_n(:,3);  %y coord
   h_n=xyd_n(:,4);  %depth 

  %plot the map and choose the transect
 hmap=figure('visible','on');clf
 tricolor(e2n(:,2:4),x_n,y_n,h_n);
 set(gcf,'renderer','zbuffer');
 colorbar;
 shading interp; hold on;
 xlabel('x')
 ylabel('y')
 axis equal;
 xt=[];
 yt=[];
 more_transect_points=true
 while(more_transect_points)
    disp('Zoom in as you like, press enter after zoom and then click on Way Points of transect :')
    title({'Map to choose transect from.'; ...
	   'Zoom in as you like,  press enter and then left-click on Way Points of transec'; ...
	   'right click to come back to zoom mode again';'middle click to finish choosing points'}, ...
	   'color',[0 0 0],'fontsize',25);
     %---allow zoom before clicking
     zoom on; % use mouse button to zoom in or out
	   % Press Enter to get out of the zoom mode.
	   % CurrentCharacter contains the most recent key which was pressed after opening
	   % the figure, wait for the most recent key to become the return/enter key
%       waitfor(gcf,'currentcharacter',13);
     zoom reset;
     zoom off;
     while(true)
	 [xt_tmp,yt_tmp,buttonPress]=ginput(1);
	 switch(buttonPress)
	       case 3      %right click will quit this one
		    break;
	       case 2      %quit completely mid-button click
		    more_transect_points=false;
		    break;
	       case 1
		    plot(xt_tmp,yt_tmp,'r*');
		    xt=[xt ;xt_tmp];   %add this point xt_tmp, yt_tmp
		    yt=[yt ;yt_tmp];
	 end
     end
     tmp_input=input('please input 1 to end this,0 to continue:')  %WLong: had to repeat twice !! to make it work
     tmp_input=input('please input 1 to end this,0 to continue:')  %the first time is skipped by matlab by right-clicking

     if (tmp_input==1)
	more_transect_points=false;
     else
	delete(hmap); %replot and contiune choosing points
	hmap=figure('visible','on');clf
	set(gcf,'renderer','zbuffer');
	tricolor(e2n(:,2:4),x_n,y_n,h_n);
	colorbar;
	shading interp; hold on;
	plot(xt,yt,'r*');
	axis equal;
     end
  end
 ds_transect=ds;
 s = [0 cumsum(sqrt(diff(xt).^2+diff(yt).^2))' ]';
 if(ds_transect >=0.5*max(s(:)))
   warning(['the segment length you input is too large']);
%   warning(['setting to default value of ' num2str(max(s(:))) '/100']);
%   ds_transect= max(s(:))/100.0;
 end
 dx=ds_transect;
 st=[s(1):dx:max(s) max(s)];
 xp=spline(s,xt,st);
 yp=spline(s,yt,st);
 %plot the chosen transect using a red line
 plot(xp,yp,'k.',xt,yt,'r*'); 
 grid;
 Ntrans=length(st);
end


