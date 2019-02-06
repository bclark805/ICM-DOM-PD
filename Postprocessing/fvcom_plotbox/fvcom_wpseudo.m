function [w,wts]=fvcom_wpseudo(ncv,niec,ntrg,dt1,dz1,u,v,dltye,dltxe,d1,uf,vf,...
                               spherical,northpole,wet_dry, semi_implicit,one_d_model,isonb,numqbc,qdis,vqdist, ...
                               nmfcell,node_mfcell,node_bfw,rdisq,mfqdis,rdismf,mfdist,qevap3,qprec3,ibfw,...
                               iswetnt,iswetn,iswetc,nv,ifceta,kb,kbm1,inflow_type,mean_flow,n,m,rofvros,art1,dz,d,dt,dti,elf,h) 


%==============================================================================|
%   calculate the sigma coordinate vertical velocity for the 3d mode (omega)   |
%                                                                              |
%   determined from equation:                                                  |
%                                                                              |
%   d/dt(d) + d/dx(ud) + d/dy(ud) = d/sigma(omega)                             |
%==============================================================================|

   xflux = zeros([n,kbm1]);

   w=zeros([n,kb]);
   wts=zeros([m,kb]); 
 
   one_third=1.0/3.0; 
  
   for i=1:ncv
     i1=ntrg(i);
     ia=niec(i,1);
     ib=niec(i,2);

     for k=1:kbm1
        if ~semi_implicit
           dij=dt1(i1)*dz1(i1,k);
           uij=u(i1,k);
           vij=v(i1,k);
           exflux=dij*(-uij*dltye(i)+vij*dltxe(i));
        else
           dij=dt1(i1)*dz1(i1,k);
           dij1=d1(i1)*dz1(i1,k);
           uij=u(i1,k);
           vij=v(i1,k);
           uij1=uf(i1,k);
           vij1=vf(i1,k);
           exflux=(1.0-ifceta)*dij*(-uij*dltye(i)+vij*dltxe(i))+ifceta*dij1*(-uij1*dltye(i)+vij1*dltxe(i));
        end
           xflux(ia,k)=xflux(ia,k)-exflux;
           xflux(ib,k)=xflux(ib,k)+exflux;
      end
   end

   if spherical && northpole
	  if ~semi_implicit
		call fvcom_vertvl_edge_xy(xflux,0.0,kbm1,ncedge_lst,node_northpole);
	  else
		call fvcom_vertvl_edge_xy(xflux,ifceta,kbm1,ncedge_lst,node_northpole);
	  end
   end
   
%-----------------------nullify boundary flux----------------------------------%
% for "tide + meanflow"/"meanflow only" case, this part should be commented out;
% for "tide only" case, this part may be kept.
% however, the effect of this term is small from my experience.
  
  if ~ mean_flow
      for i=1:m
        for k=1:kbm1
            if(isonb(i) == 2)
               xflux(i,k)=0.0  ;
            end
        end
      end
% can be changed to (no if statements)
%     for i=1:iobcn
%         for k=1:kbm1
%             xflux(i_obc_n(i),k)=0.0
%         end
%     end
  end

  if  (one_d_model)
    xflux = 0.0;
  end


%-----------------------fresh water inflow-------------------------------------%
  if(numqbc >= 1) 
     if(inflow_type == 'node') 
       for j=1:numqbc
         jj=inodeq(j)
         for k=1:kbm1
           xflux(jj,k)=xflux(jj,k)-qdis(j)*vqdist(j,k) ;   %/dz(jj,k)
         end
       end
     elseif(inflow_type == 'edge') 
       for j=1:numqbc
         j1=n_icellq(j,1);
         j2=n_icellq(j,2);
         for k=1:kbm1
             xflux(j1,k)=xflux(j1,k)-qdis(j)*rdisq(j,1)*vqdist(j,k)  ;  %/dz1(j1,k)
             xflux(j2,k)=xflux(j2,k)-qdis(j)*rdisq(j,2)*vqdist(j,k)  ;  %/dz1(j2,k)
         end
       end
     end
  end

  if mean_flow
     if (nmfcell > 0)
        for i = 1:nmfcell
            j1= node_mfcell(i,1)
            j2= node_mfcell(i,2)
            for k=1:kbm1
                xflux(j1,k) = xflux(j1,k) - mfqdis(i)*rdismf(i,1)*mfdist(i,k) ;    %/dz1(j1,k)
                xflux(j2,k) = xflux(j2,k) - mfqdis(i)*rdismf(i,2)*mfdist(i,k) ;    %/dz1(j2,k)
             end
         end
      end
  end

%---if no fresh water inflow, omega is zero at free surface and bottom---------%

  wbottom   = zeros([m,1]);
  wts(1:m,kb) = 0.0;
  for i=1:m
     wts(i,1) = (qevap3(i)-qprec3(i))*rofvros    ;    %0.0
     if(ibfw > 0)
       for j=1:ibfw
          if(i == node_bfw(j))
             wbottom(i)= bfwdis3(j)/art1(i);
          end
       end
     end
  end

%--------------------------calculate omega-------------------------------------%

  for i=1:m
     if(iswetnt(i)*iswetn(i) == 1)
        for k=1:kbm1
           wts(i,k+1)=wts(i,k)+xflux(i,k)/art1(i)+dz(i,k)*(d(i)-dt(i))/dti;
        end
     else
        for k=1:kbm1
           wts(i,k+1)=0.0;
        end
     end
  end

%--------------------------adjust omega----------------------------------------%
% improves mass conservation
 
  for i=1:m
     if(abs(wts(i,kb)-wbottom(i)) > 1.0e-8)

%#  if ~mean_flow  %take 2 out if not meanflow, do everything if meanflow
       if((isonb(i) ~= 2 && ~mean_flow) || mean_flow)
%#  endif
        if (~semi_implicit)
	         tmp1=elf(i)*double(kbm1)-(wts(i,kb)-wbottom(i))*dti/dz(i,1);
        else
	         tmp1=el(i)*double(kbm1)-(wts(i,kb)-wbottom(i))*dti/dz(i,1);
        end
        tmp1=tmp1/double(kbm1);
        dtfa(i)=tmp1+h(i);
        for k=2:kb
           wts(i,k)=wts(i,k)-double(k-1)/double(kbm1)*(wts(i,kb)-wbottom(i));
        end
%#  if ~mean_flow
       end
%#  endif
     end
  end
%
%----transfer omega to face center---------------------------------------------%
%
  for i=1:n
      for k=1:kb
          w(i,k) = one_third*(wts(nv(i,1),k)+wts(nv(i,2),k)+wts(nv(i,3),k));
      end
  end

  if (wet_dry)
     for i=1:n
     for k=1:kb
       w(i,k) = double(iswetc(i))*w(i,k);
     end
     end
  end
  return
end
