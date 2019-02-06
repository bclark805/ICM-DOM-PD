function [ww]=wreal(w,u,v,zz1,d,elf,awx,awy,el,kbm1,n,nv,iswetc,wet_dry,et1,dti)
%==============================================================================|
%  compute cartesian vertical velocity                                         |
%==============================================================================|
%  Inputs:
% 	u(1:nt,kb)		- u velocity, on element center and layer center
%  	v(1:nt,kb)              - v velocity, on element center and layer center
%  	w(1:nt,kb)              - pseudo vertical velocity, on element center and levels
%  	zz1(1:nt,kb)            - intra level sigma value, on element centr and layer center
%  	elf(0:mt)               - surface elevation at nodes
%  	d(1:mt)                 - current water depth (total) on nodes
%  	awx(1:nt,3)             - coefficients of estimating x gradient of function f within an element based
% 	                          on nodal values of the element
% 	        	                  df/dx = awx(i,1)*f1+awx(i,2)*f2+awx(i,3)*f3
%       	        	          df/dy = awy(i,1)*f1+awy(i,2)*f2+awy(i,3)*f3
%
%  	awy(1:nt,3)             - coefficients of estimating y gradient of function f within an element based
%                                 on nodal values of the element
%  	el(1:mt)  		- elevation on nodes
%  	kbm1			- number of layers
%  	n			- number of elements
%  	nv(1:nt,1:3)		- e2n
%  	iswetc(1:nt)		- porosity at element center
%  	wet_dry                 - flag for wet_dry, calculation 0-no, 1-yes
%       et1(1:nt)               - surface elevation on element center at previous time step
%       dti                     - time step (sec)
%   Outputs:
%       ww(1:nt,kb)             - vertical velocity at layer center and element center
%

   one_third=1.0/3.0; 
   for i=1:n
	    if(iswetc(i) == 1 )
		     j1=nv(i,1);
		     j2=nv(i,2);
		     j3=nv(i,3);
		     dddx=awx(i,1) * d(j1)+awx(i,2) * d(j2)+awx(i,3)*d(j3);
		     dddy=awy(i,1) * d(j1)+awy(i,2) * d(j2)+awy(i,3)*d(j3);
		     dedx=awx(i,1)*elf(j1)+awx(i,2)*elf(j2)+awx(i,3)*elf(j3);
		     dedy=awy(i,1)*elf(j1)+awy(i,2)*elf(j2)+awy(i,3)*elf(j3);
		     etf1aa=one_third*(el(nv(i,1))+el(nv(i,2))+el(nv(i,3)));
	
		     do k=1:kbm1
			       ww1=0.5(w(i,k)+w(i,k+1))+u(i,k)*(zz1(i,k)*dddx+dedx)+ ...
		                                    v(i,k)*(zz1(i,k)*dddy+dedy);
			       ww2=(zz1(i,k)+1.)*(etf1aa-et1(i))/dti;
			       ww(i,k)=ww1+ww2
		     end
	    else
		     do k=1:kbm1
			      ww(i,k)=0.0
		     end
	    end
   end
   return
end 
