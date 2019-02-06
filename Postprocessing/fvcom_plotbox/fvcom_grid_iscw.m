function [cw]= fvcom_grid_iscw(e2n,x_n,y_n)
%function [cw] =  fvcom_grid_iscw(e2n,x_n,y_n)
%
%    Program to check if fvcom grid data stored in e2n(1:n,1:3), x_n(1:m),y_n(1:m) is clockwise or counter-clockwise
%  
%    Inputs: x_n(1:m), y_n(1:m) --- x and y coordinates of nodes. m is number of nodes
%            e2n(1:n,1:3)       --- Element to node connectivity of elements. n is number of elements
%                                   Each element consists of 3 nodes.
%
%    Output: cw=1 if all elements are clockwise
%            cw=0 if all elements are counter-clockwise
%            Code will give error if some of elements are clockwise and some are counter-clockwise
%
%

    %check number of nodes

    m_x=length(x_n(:));
    m_y=length(y_n(:));

    if(m_x ~= m_y)
       error('oops, x_n and y_n should have same size and be 1D');
    end
 
    m=m_x;

    %check number of elements

    n=size(e2n,1);

    n3=size(e2n,2);

    if(n<1)
      error('oops, e2n is empty?');
    end
 
    if(n3~=3)
      error('oops, second dimension size of e2n must be 3');
    end

    %check if e2n is in the range of 1:m

    if(max(e2n(:))> m ||min(e2n(:)) <1)
      error(['oops, e2n values invalid, should be of range 1 to ' int2str(m)]);
    end

    cw_e2n=zeros([n,1]);

    for i1=1:n %check whether clockwise for each element
        if(WLispolycw_WL(x_n(e2n(i1,1:3)),y_n(e2n(i1,1:3))))
           cw_e2n(i1)=1;
        else
           cw_e2n(i1)=0;
        end
    end
    cw_tmp_max=max(cw_e2n(:));
    cw_tmp_min=min(cw_e2n(:));

    if(cw_tmp_max==cw_tmp_min)
	cw=cw_tmp_max;
        return;
    else
	error('oops, the grid''s orient is inconsistent, some elements are clockwise, some elements are coutner clockwise');
    end
  
end
