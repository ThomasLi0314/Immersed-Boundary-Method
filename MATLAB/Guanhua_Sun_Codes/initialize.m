%initialize.m
%in seconds,gram,cm
Lx=16
Ly=8
Nx=256
Ny=Nx/Lx*Ly

h=Ly/Ny
hb=h/4
ipx=[(2:Nx),1]
imx=[Nx,(1:(Nx-1))]
ipy=[(2:Ny),1]
imy=[Ny,(1:(Ny-1))]

K=25600000
rho=1    
mu=0.010518
tmax=8
dt=0.000064
clockmax=ceil(tmax/dt)
flow = zeros(clockmax,7);
%geometry;
%initial velocity
u=zeros(Nx,Ny,2);



%[vX,vY] = meshgrid (0:h:L-h);
%quiver (vY,vX,bf(:,:,1),bf(:,:,2),5);

    
xgrid=zeros(Nx,Ny);
ygrid=zeros(Nx,Ny);
for j=0:Nx-1
    xgrid(j+1,:)=j*h;
end

for jj=0:Ny-1
    ygrid(:,jj+1)=jj*h;
end

set(gcf,'double','on')
%contour(xgrid,ygrid,vorticity,values)
%hold on


axis([0,Lx,0,Ly])
%caxis(valminmax)
axis equal
axis manual
drawnow
hold off

