[xxc yyc]=function arcc_fvcom(xx1,yy1,xx2,yy2,xxc,yyc)
%
%[xxc yyc]=function arcc_fvcom(xx1,yy1,xx2,yy2)
%
%calculate lon and lat of middle point C of 
%arc (A,B) with A at (lon,lat)=(xx1,yy1)
%B at  (lon,lat) =(xx2,yy2) on a sphere
%are (A,B) is takein as the great circle of AB on sphere
%midle point C is on middle of great circle AB
% 
% xx1--- lon of point A (deg), [0,360]
% yy1--- lat of point A (deg), [-90, 90]
% 
% xx2--- lon of point B (deg), [0,360]
% yy2--- lat of point B (deg), [-90,90]
%
% xxc--- lon of point C (mid point on great circle AB) (deg), [0,360]
% yyc--- lat of point C (mid point on great circle AB) (deg), [-90, 90]
%

   deg2rad=3.14159265/180.0;

   %find lon and lat in rad for the points
   x1=xx1*deg2rad;
   y1=yy1*deg2rad;

   x2=xx2*deg2rad;
   y2=yy2*deg2rad;

  
   if(abs(yy1-90.0) < 1.0e-6)  %if first point is polar
     xxc = xx2;                       %then take xxc the same as xx2 (degrees) (longitude)
   elseif(abs(yy2-90.0) < 1.0e-6) %if second point is polar
                                        %then take xxc as xx1 (longiude)
     xxc = xx1;
   else     %if neither point is on the poles

	   xxc=cos(y1)*sin(x1)+cos(y2)*sin(x2);   %cos(lat1)*sin(lon1)+cos(lat2)*sin(lon2)

	   xxc=atan2(xxc,cos(y1)*cos(x1)+cos(y2)*cos(x2));  %angle of [cos(lat1)*cos(lon1)+cos(lat2)*cos(lon2), 
                                                           %          cos(lat1)*sin(lon1)+cos(lat2)*sin(lon2)]

	   xxc=xxc/deg2rad;              %convert from rad to deg

	   if(xxc < 0.0) xxc=360.0+xxc ; %get into [0,360] degree range   %xxc would be longitude of the mid point of arc 
                                        %                         between (xx1,yy1) and (xx2,yy2)
                                        %

   end

   yyc=    cos(y1)*cos(y1)+cos(y2)*cos(y2)   ...  %    cos(lat1)*cos(lat1)+cos(lat2)*cos(lat2)
       +2.*cos(y1)*cos(y2)*cos(x1-x2);            %  +2*cos(lat1)*cos(lat2)*cos(lon1-lon2)

   yyc=atan2(sqrt(yyc),sin(y1)+sin(y2));          %angle of [sin(lat1)+sin(lat2),
                                                  %           cos(lat1)*cos(lat1)+cos(lat2)*cos(lat2) 
                                                  %           +2*cos(lat1)*cos(lat2)*cos(lon1-lon2)

	
   %yyc is latitude of mid point of great arc of points connected by (lon=xx1,lat=yy1) and (lon=xx2,lat=yy2) on the sphere

   yyc=90.-yyc/deg2rad   ; % makes sure yyc is of range [-90, 90] (deg)

end

