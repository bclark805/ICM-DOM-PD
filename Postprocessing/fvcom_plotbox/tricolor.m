function htri=tricolor(tri,x,y,v,cmin,cmax)
%FUNCTION [htri]=TRICOLOR(TRI,X,Y,V,CMIN,CMAX) plots color images of V defined
%at the positions specified by X and Y for triangular
%grids locations defined using indices of X, Y stored in TRI
%
%CMIN, CMAX --minimum value and max value in V that will be used for 
%             showing color
%X --- vector of x coordinate
%Y --- vector of Y coordinate
%V --- vector of Z coordinate (V does not have the same size as X, Y)
%      if V has same length as X, Y, then V will be used as vertex color value
%      to provide color on each face with biliner interpolation. If V has same size 
%      as TRI's first dimension then V will be used to color each face
%
%TRI --- (NE x 3) matrix, each row is a triplet of integers pointing
%        to the indices in X, Y to form a triangular element
%
%htri --- output handle of the patch object of the color image
%
%A colourbar is added on the right side of the figure.
%
% The colorbar strectches from the minimum value of v to its
% maximum in 9 steps (10 values).
%
% The plot is actually consits of patch plots for each triangle
% plane 2D plot. 
%
% Example:
%    [x,y]=meshgrid(1:15,1:15);
%    tri = delaunay(x,y);
%    z = peaks(15);
%    tricolor(tri,x,y,z);
%
%Wen Long, PNNL wen.long@pnnl.gov, Feb 2 2012


%check number of arguments
  if nargin <4
     error('Oops, wrong number of arguments');
  end

% simple check on tri
  if (size(tri,2) ~=3)
     error('Oops, wrong size of tri'); 
  end
  
%get the color
%delete(gca)
map=colormap;
if nargin < 6
  miv=min(v);
  mav=max(v);
  delm=mav-miv;
  if(delm==0)  %make sure color range valid
     delm=abs(mav)/1000;
  end
  if(delm==0)  %again make color range valid when miv=mav=0
     delm=0.001;
  end
  miv=miv-0.01*delm;
  mav=mav+0.01*delm;
else
  miv=cmin;
  mav=cmax;
end  
clrstep = (mav-miv)/(size(map,1)-1) ;

carray=(miv:clrstep:mav);  %array of values for color axis
cindex_array=(1:1:size(map,1));  %index array in map

%---------obsolete and slow----
% % Plot the triangles 
%hold on
%ccvalue=[];
%for ie=1:size(tri,1)  %loop through element
%   %average of 3 nodes to get value for the element
%    vertexdata=[v(tri(ie,1)); v(tri(ie,2)); v(tri(ie,3))];
%    
%    elevalue=(v(tri(ie,1))+v(tri(ie,2))+v(tri(ie,3)))/3.0; 
%   %get color index in the map 
%    cindex=max(1,floor(interp1(carray,cindex_array,elevalue,'nearest')));
%    elecvalue=map(cindex,:); %get color triplet of the patch for the element
%    %patch(x(tri(ie,:)),y(tri(ie,:)),elecvalue); 
%    %patch(x(tri(ie,:)),y(tri(ie,:)),vertexdata); 
%    ccvalue=[ccvalue cindex]; 
%end
%

%---easy and fast----
vcindex=int32(fix(interp1(carray,cindex_array,v,'nearest')));
ifill_down=find(v<=miv);  %find places below miv
ifill_up  =find(v>=mav);  %find places above mav
vcindex(ifill_down)=1;     %use lowerest color
vcindex(ifill_up)=size(map,1);  %use highest color
htri = patch('faces',tri,'vertices',[x(:) y(:)],'facevertexcdata',map(vcindex,:));
%-------------
shading flat;
colorbar;
caxis([miv,mav]);
hold off


