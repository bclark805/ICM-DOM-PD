function [xx,yy,ie_tce,ie_bce]=fvcom_tce_polygon(spherical,ne,ienode,ntve,nv,vx,vy,isonb,nbve,nbvt,xc,yc,ncv,niec,ntrg,itce)
%
%function [xx,yy,ie_tce,ie_bce]=fvcom_tce_polygon(spherical,ne,ienode,ntve,nv,vx,vy,isonb,nbve,nbvt,xc,yc,ncv,niec,ntrg,itce)
%
%    this is used to calculate the polygon of the itce'th TCE, which is a tracer control element
%    in fvcom grid
%  
% Outputs
%                                                                         
%  xx   : x coordinates of tce itce, 1<=itce<=m, points are stored counter clociwise
%  yy   : y coordinates of tce itce, 1<=itce<=m, points are stored counter clockwise
%  ie_tce(1:length(xx)-1): the tce edge indices for the edges of the tce polygon when
%                              the edge is a tce edge, otherwise value is 0
%  ie_bce(1:length(xx)-1): the element edge indices for the edges of the tce polygon (when itce is a 
%                              a node on the boundary), otherwise value is 0
%
% Inputs
%
%  spherical : 1 for spherical coord (lon,lat), 0 for planar
%  ne        : number of edges
%  ienode(1:ne,1:2) : startiing and ending node numbers of a given edge ie=1...ne
%  vx(1:m)   : x coordinate of nodes
%  vy(1:m)   : y coordinate of nodes
%  nv(1:n,3) : element to node connectivity
%  isonb     : flag of boundary type of a given node
%                  0--this node is in interior
%                  1--this node is on solid boundary
%                  2--this node is on open boundary
%  ntve(1:m) :               : total number of surrounding elements of a node
%  nbve(1:m,1:mx_nbr_elem+1) : element number of surrounding elements of a given node
%  nbvt(1:m,1:mx_nbr_elem+1) : 1, 2, or 3, i.e. the local node id of an surrounding element
%                              that happens to be the same as a given node
%
%  xc(1:n)   : x coordinate of element center
%  yc(1:n)   : y coordinate of element center
%   
%  i         : tce number i (i.e. node number)
%
% Code and documentation:
%
%  Wen Long, PNNL/BSRC, Seattle, 04/22/2013

  m=size(vx(:,1));  %number of nodes
  n=size(nv,1);     %number of elements

  nhn=0;      %number of halo nodes
  nhe=0;      %number of halo elements
  nt=n+nhe;   
  mt=m+nhn; 

  artmax=0;
  arttot=0;
  artmin=0;

  i=0;
  j=0;
  ii=0;
  j1=0;
  j2=0;
  max_nbre=0; 

  if (spherical)   
  	   vx1=0;
	   vx2=0;
	   vx3=0;

	   vy1=0;
	   vy2=0;
	   vy3=0;

         side1=0;
	 side2=0;
	 side3=0;

	 area1=0;

	   xx1=0;
	   yy1=0;
	   xx2=0;
	   yy2=0;
	   xxc=0;
	   yyc=0;

	   x1_dp=0;
	   y1_dp=0;
	   x2_dp=0;
	   y2_dp=0;
	   x3_dp=0;
	   y3_dp=0;
   end

   max_nbre = max(ntve(:))+1 ; %maximum number of elements around a node 

   i=itce;
   if(isonb(i) == 0) %node i is interior
       xx     = zeros([2*ntve(i)+1,1]);
       yy     = zeros([2*ntve(i)+1,1]);
       ie_tce = zeros([2*ntve(i),1]);
       ie_bce = zeros([2*ntve(i),1]);

       %2j-1 and 2j
       for j=1:ntve(i)  %loop through all surrounding elements
           ii=nbve(i,j);  %element number of j'th surrounding element
           j1=nbvt(i,j);  %local node number that node i is for elelment ii
           j2=j1+1-floor((j1+1)/4)*3; %neight local node number in the element to node connectivity

	   if (spherical)
		 xx1=vx(nv(ii,j1));
		 yy1=vy(nv(ii,j1));
		 xx2=vx(nv(ii,j2));
		 yy2=vy(nv(ii,j2));
		 [xxc,yyc]=arcc_fvcom(xx1,yy1,xx2,yy2); 
		 xx(2*j-1)=xxc;
		 yy(2*j-1)=yyc ;          
		 xx(2*j)=xc(ii);
		 yy(2*j)=yc(ii) ;           
	   else
		 xx(2*j-1)=(vx(nv(ii,j1))+vx(nv(ii,j2)))*0.5;  %mid point of the spoke
		 yy(2*j-1)=(vy(nv(ii,j1))+vy(nv(ii,j2)))*0.5;  
                                            %j1 and j2 uniquely identifies an edge
                                            %j1 and j2 and i uniquely identify a TCE edge with ib equal to i
	                                            %we should note down that TCE edge under the name
                                                    %of 2j-1

                                           %j1 and j2 and i uniquely identify a TCE edge with ia equal to i
                                                    %we should note down that TCE edge under the name
                                                    %of 2j

		 xx(2*j)=xc(ii);                               %center of the element
		 yy(2*j)=yc(ii);
	   end

           for itce_edge=1:ncv  %loop though all tce edges
               ia=niec(itce_edge,1);  %left side node of the tce edge
               ib=niec(itce_edge,2);  %right side node of the tce edge
               iele=ntrg(itce_edge);  %element number that the tce edge is in

               %identify the tce edge index in ncv that has end point as the spoke mid point xx(2*j-1),yy(2*j-1)
               %and the tce edge is within element ii

               if(ib==i && iele==ii )
                  ie_tce(2*j)=itce_edge;  %register the tce edge

               end
               %identify the tce edge index in ncv that has starting point as element ii'th centroid xc(ii),yc(ii)
               %i.e. xx(2*j), yy(2*j), and the tce edge is within element ii

               if(ia==i && iele==ii)
                  ie_tce(2*j-1)  =itce_edge;  %register the tce edge
               end
           end

       end
       %repeat for the last point
       xx(2*ntve(i)+1)=xx(1);
       yy(2*ntve(i)+1)=yy(1);         
       %we might also repeast the last TCE edge number
       %ie_tce(2*ntve(i)+1)=ie_tce(1);  %may not be necessary %/*debugging */

   else  %node i is on boundary

       xx     = zeros([2*(ntve(i)+1)+1,1]);
       yy     = zeros([2*(ntve(i)+1)+1,1]);
       ie_tce = zeros([2*(ntve(i)+1)  ,1]);  %tce edge index in niec, ntrg
       ie_bce = zeros([2*(ntve(i)+1)  ,1]);  %edge in dex

       iii=nbve(i,1);                 %get the first surrounding element of node i (counter clockwise for counter-clockwise grids)
       jj1=nbvt(i,1);                 %get the index in e2n for node i in element
                                      % iii=nbve(i,1) (first surrounding element)
       jj2=jj1+1-floor((jj1+1)/4)*3;  %get index of the next node in e2n for that element

       for j=1:ntve(i)  %loop through surrounding elements
         ii=nbve(i,j);  %get the element 
         j1=nbvt(i,j);  %get the index in e2n for node i in element ii
         j2=j1+1-floor((j1+1)/4)*3;  %get index of the next node in  e2n 
	 if(spherical)
		 xx1=vx(nv(ii,j1));
		 yy1=vy(nv(ii,j1));
		 xx2=vx(nv(ii,j2));
		 yy2=vy(nv(ii,j2));
		 [xxc,yyc]=arcc_fvcom(xx1,yy1,xx2,yy2);
			  
		 xx(2*j-1)=xxc;
		 yy(2*j-1)=yyc;            
		     
		 xx(2*j)=xc(ii);
		 yy(2*j)=yc(ii);
	 else
		 xx(2*j-1)=(vx(nv(ii,j1))+vx(nv(ii,j2)))*0.5;  %midpoint of the spoke
		 yy(2*j-1)=(vy(nv(ii,j1))+vy(nv(ii,j2)))*0.5;
		 xx(2*j)=xc(ii);                               %center of the element
		 yy(2*j)=yc(ii);
	 end	                

         for itce_edge=1:ncv  %loop though all tce edges
             ia=niec(itce_edge,1);  %left side node of the tce edge
             ib=niec(itce_edge,2);  %right side node of the tce edge

             iele=ntrg(itce_edge);  %element number that the tce edge is in

             %identify the tce edge index in ncv that has end point as the spoke mid point xx(2*j-1),yy(2*j-1)
             %and the tce edge is within element ii

             if(ib==i && iele==ii)
                ie_tce(2*j)=itce_edge;  %register the tce edge
             end

             %identify the tce edge index in ncv that has starting point as element ii'th centroid xc(ii),yc(ii)
             %i.e. xx(2*j), yy(2*j), and the tce edge is within element ii

             if(ia==i &&  iele==ii)
                ie_tce(2*j-1)  =itce_edge;  %register the tce edge
             end

          end

       end

        j=ntve(i)+1;
       ii=nbve(i,j-1);      %the last surrounding element
       j1=nbvt(i,ntve(i));  %local index in e2n of the element ii that node i is in element ii
       j2=j1+2-floor((j1+2)/4)*3;  %next index in e2n for element ii
                                   %j1,j2 finds the edge of the element ii, and j1,j2 is a boundary edge
       if (spherical)
	       xx1 =  vx(nv(ii,j1));
	       yy1 =  vy(nv(ii,j1));
	       xx2 =  vx(nv(ii,j2));
	       yy2 =  vy(nv(ii,j2));
	       arcc_fvcom(xx1,yy1,xx2,yy2,xxc,yyc) ;

	       xx(2*j-1)=xxc;
	       yy(2*j-1)=yyc;

	       xx(2*j)=vx(i);
	       yy(2*j)=vy(i);          

	       xx(2*j+1)=xx(1);
	       yy(2*j+1)=yy(1);

       else
	       xx(2*j-1)=(vx(nv(ii,j1))+vx(nv(ii,j2)))*0.5;  %mid point of the boundary edge j1,j2
	       yy(2*j-1)=(vy(nv(ii,j1))+vy(nv(ii,j2)))*0.5;  %

	       xx(2*j)=vx(i);                                %back to node i
	       yy(2*j)=vy(i);

	       xx(2*j+1)=xx(1);                              %back to starting point of the TCE polygon
	       yy(2*j+1)=yy(1);
       end

       iii=nbve(i,1);                 %get the first surrounding element of node i (counter clockwise for counter-clockwise grids)
       jj1=nbvt(i,1);                 %get the index in e2n for node i in element
                                      % iii=nbve(i,1) (first surrounding element)
       jj2=jj1+1-floor((jj1+1)/4)*3;  %get index of the next node in e2n for that element

       %find the edge number of the boundary edges
       for ie=1:ne  %loop through all edges
           ienode1=ienode(ie,1);  %starting node number of edge ie
           ienode2=ienode(ie,2);  %ending node number of edge ie
           if(   ienode1==nv(ii,j1) && ienode2==nv(ii,j2)  ...
              || ienode1==nv(ii,j2) && ienode2==nv(ii,j1))
             ie_bce(2*j-1)=ie; 
           end
           if(   ienode1==nv(iii,jj1) && ienode2==nv(iii,jj2)    ...
              || ienode1==nv(iii,jj2) && ienode2==nv(iii,jj1))
             ie_bce(2*j)=ie;
           end
       end

   end

end
