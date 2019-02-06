function h=pcolor2(X,Y,Z)

if nargin==1
  Z=X;
 [ny,nx]=size(X);
 [X,Y]=meshgrid(1:nx,1:ny);
end

W=ones(2,2)/4;

X2=conv2(X,W);
dx1=(X2(:,3)-X2(:,2));
X2(:,1)=X2(:,2)-dx1;
dxn=(X2(:,end-1)-X2(:,end-2));
X2(:,end)=X2(:,end-1)+dxn;
X2(1,:)=X2(1,:)*2;
X2(end,:)=X2(end,:)*2;
X2;

Y2=conv2(Y,W);
dy1=(Y2(3,:)-Y2(2,:));
Y2(1,:)=Y2(2,:)-dy1;
dyn=(Y2(end-1,:)-Y2(end-2,:));
Y2(end,:)=Y2(end-1,:)+dyn;
Y2(:,1)=Y2(:,1)*2;
Y2(:,end)=Y2(:,end)*2;
Y2;

[ny,nx]=size(Z);
Z2=ones(ny+1,nx+1);
Z2(1:ny,1:nx)=Z;

h=pcolor(X2,Y2,Z2);
