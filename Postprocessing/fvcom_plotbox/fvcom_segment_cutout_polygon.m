function [polycuts]=fvcom_segment_cutout_polygon(tri,vx,vy,xije,yije,itce,xp,yp,i_bce,i_tce,xi,yi,segments)
%function [polycuts]=fvcom_segment_cutout_polygon(tri,vx,vy,xije,yije,itce,xp,yp,i_bce,i_tce,xi,yi,segments)
%
%    Calculates the cutout polygon due to segments created by a line going through an existing polygon xp, yp
%    with the intersection points of the line with the polygon stored in xi, yi
%    and the individual segemnts stored in segments. i_bce,i_tce gives the edge indices of all tce edges
%    and all triangle edges that the original polygon xp yp is made of. i.e. the xp, yp polygon is made of edges
%    either from tce edge pool or from triangle edge pool.
%    Also xp,yp is a tce polygon for node point itce. 
%

% save mycut.mat tri vx vy xije xi yi segments i_bce i_tce xp yp -mat;
% load mycut.mat

polycuts=[];
if(length(xi(:))>2)
%-------

      debug_tmp=0;

      if(debug_tmp)
	figure;
	trimesh(tri,vx,vy,'color',[0 0 0]); hold on
	%show the TCE
	line(xije',yije','color',[0 1 1],'linestyle','--');
	%show the TCE vertecies
	plot(xije(:,2),yije(:,2),'k.');  %edge mid points
	plot(xije(:,1),yije(:,1),'ko');  %element centers
	for i=1:100
	    text(vx(i),vy(i),num2str(i),'color',[1 0 0]);
	end
	axis([-3000,1000,-500,2800]);
      end


	for iseg=1:length(xi)-1

	     %plot the cut-out polygon
             iedge1=segments{iseg}.iedge(1);%edge number of the polygon for the starting
	     iedge2=segments{iseg}.iedge(2);%edge number of the polygon for the starting
					    %point of the segment

                xpc=[];  %cut-out polygon to be constructed
                ypc=[];
                itce_c=[];  %record the edge number of TCE edge that a point on the xpc is part of, 0 if not on any TCE edge
                ibce_c=[];  %record the edge number of a boundary element edge that a point on the xpc is part of, 0 if not on any boudary edge (BCE)
                            %a point must be either on an TCE edge or on a boundary edge BCE, so itce_c and ibce_c can't be both 0.

                search_edges=(1:length(xp));
	    

		if(iedge1*iedge2 >0 && segments{iseg}.inout <=0) %both points are on the polygon xp and the segment is in the polygon

                    iedge=iedge2;
		    ipoint=iedge2;
                    xpc=[segments{iseg}.x(2)];  %starts from point 2
                    ypc=[segments{iseg}.y(2)];

                    itce_c=[i_tce(iedge)];      %index of tce edge that xpc point is part of, 0 if not on tce edges
                    ibce_c=[i_bce(iedge)];      %index of bce edge that xpc point is part of, 0 if not on boundary edges
	            
		    finished_search=0;

                    %loop until finding the first point of the segment
                    while(finished_search==0)

                       iedge=iedge+1;

		       ipoint=ipoint+1;

                       %if this new edge contains point 1 of the segment
                        %starting point of edge iedge
                        %ending point of edge iedge
                        istart=ipoint;
                        iend=istart+1;
                        if(istart>length(xp))
                           istart=mod(istart,length(xp))+1;
                        end
                        if(iend>length(xp))
                           iend=mod(iend,length(xp))+1;
                        end

                        edge_x=[xp(istart), xp(iend)];
                        edge_y=[yp(istart), yp(iend)];

	                iedge_use=mod(iedge,length(xp)-1);
			if(iedge_use==0)
			    iedge_use=length(xp)-1;      %prevent iedge_use from goiing over limit = length(xp)-1
                        end

                        if(pointonline(edge_x,edge_y,segments{iseg}.x(1),segments{iseg}.y(1)))
                             %if this edge of the polygon contains point 1 of the segment
			     %add a new point to the cutout polygon only the new point is distinct from last point

			     if(xp(istart)~=xpc(end) || yp(istart) ~= ypc(end))
				xpc=[xpc,xp(istart)];
				ypc=[ypc,yp(istart)];

	 		        itce_c=[itce_c,  i_tce(iedge_use) ];
                                ibce_c=[ibce_c,  i_bce(iedge_use) ];
			     end
			     if(segments{iseg}.x(1) ~= xpc(end) || segments{iseg}.y(1) ~= ypc(end))
				xpc=[xpc,segments{iseg}.x(1)];
                                ypc=[ypc,segments{iseg}.y(1)];
				itce_c=[itce_c, i_tce(iedge1) ];
                                ibce_c=[ibce_c, i_bce(iedge1) ];
			     end
                             finished_search=1;  %terminate the search
                        else
                             %if this new edge does not contain point 1 of the segment
                             if(xp(istart)~=xpc(end) || yp(istart) ~= ypc(end))
                                xpc=[xpc,xp(istart)];
                                ypc=[ypc,yp(istart)];
                                itce_c=[itce_c,  i_tce(iedge_use) ];
                                ibce_c=[ibce_c,  i_bce(iedge_use) ];
                             end
                             %move on to the next edge
                        end
                    end

		    %remove last edge to make sure number of edges is one less than number of points on xpc
                    itce_c(end)=[];
                    ibce_c(end)=[];

		    %need to remove repeating points and corresponding repeating edegs if any

		    if(debug_tmp)  %---debug----

			%plot the intersection points
			plot(xi,yi,'ro');  %intersection points
			%plot the new half polygon xpc, ypc
			line(xpc,ypc, 'color','y');
			%label the found tce edges of the cutout polygon section
			for i=1:length(xpc(:))-1
			itce_edge=itce_c(i);
			    if(itce_edge>0)
				    text_x=(xpc(i)+xpc(i+1))/2.0;
				    text_y=(ypc(i)+ypc(i+1))/2.0;
				    textangle=atan2(text_y-vy(itce),text_x-vx(itce));
				    %label the tce_edge
				    text( text_x+(max(xpc(:))-min(xpc(:)))*1/5.0*cos(textangle),...
					text_y+(max(xpc(:))-min(xpc(:)))*1/5.0*sin(textangle),...
					num2str(itce_edge),'color',[0 1 1]);
			    end

			    ibce_edge=ibce_c(i);  %if the edge of the cut out polygon is a boundary edge
			    if(ibce_edge>0)       %label the boundary edge
			        text_x=(xpc(i)+xpc(i+1))/2.0;
				text_y=(ypc(i)+ypc(i+1))/2.0;
				textangle=atan2(text_y-mean(ypc(:)),text_x-mean(xpc(:)));
				%label the bce_edge
				text( text_x+(max(xpc(:))-min(xpc(:)))*1/5.0*cos(textangle),...
				      text_y+(max(xpc(:))-min(xpc(:)))*1/5.0*sin(textangle),...
				    num2str(ibce_edge),'color',[1 0 1]);
			    end
			end
		    end		 %---debug----

		    polycuts{iseg}.xpc=xpc;        %points of the cutout polygon 
		    polycuts{iseg}.ypc=ypc;
		    polycuts{iseg}.ibce_c=ibce_c;  %and the TCE or BCE edge id's of each edge on the cutout polygn 
		    polycuts{iseg}.itce_c=itce_c;

		else
		    polycuts{iseg}.xpc=[];      %no polygons found
		    polycuts{iseg}.ypc=[];
		    polycuts{iseg}.ibce_c=[];
		    polycuts{iseg}.itce_c=[];
		end

	        %deal with the case where xpc has only one point, (i.e. when the segment's two points are colloapsing into one point)

		if(length(xpc)<=1)
		    polycuts{iseg}.xpc=[];
                    polycuts{iseg}.ypc=[];
		    polycuts{iseg}.ibce_c=[];
		    polycuts{iseg}.itce_c=[];
		end

	end  %iseg loop

end %if(length(xi(:))>2) 


