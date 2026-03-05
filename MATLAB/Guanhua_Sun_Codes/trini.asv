%initialize.m
%in seconds,gram,cm
%global angle
Lx=6;
Ly=4.5;
Nx=512;
Ny=Nx/Lx*Ly;

h=Ly/Ny;
hb=h*(2/3);
ipx=[(2:Nx),1];
imx=[Nx,(1:(Nx-1))];
ipy=[(2:Ny),1];
imy=[Ny,(1:(Ny-1))];

K=1000000;
rho=1;    
mu=0.01;
tmax=8;
dt=0.0001;
clockmax=ceil(tmax/dt);
flow = zeros(clockmax,7);
trigeo
%set_plane
%initial velocity
u=zeros(Nx,Ny,2);

xgrid=zeros(Nx,Ny);
ygrid=zeros(Nx,Ny);
for j=0:Nx-1
    xgrid(j+1,:)=j*h;
end

for jj=0:Ny-1
    ygrid(:,jj+1)=jj*h;
end

%set(gcf,'double','on')
%contour(xgrid,ygrid,vorticity,values)
%hold on
% grid = zeros((2/L*N)^2,2);
% for j=0:(2/L*N-1)
%   grid(j*N/2+1:(j+1)*N/2,1)=j*h+1;
%   for k=0:(2/L*N-1)
%   grid(j*N/2+1+k,2)=k*h+1;
%   end
% end

%appro = dsearchn(grid,X);
%Y = grid(appro,:);

%Y=X.*(N/L);
%T=zeros(1,2);
%plot(X(:,1),X(:,2),'ko')
%plot(Y(:,1),Y(:,2),'r.')
%axis([0,Lx,0,Ly])
%caxis(valminmax)
%axis equal
%axis([0,Lx,0,Ly])
%drawnow
%hold off

