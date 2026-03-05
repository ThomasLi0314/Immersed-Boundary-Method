%initialize.m
%in seconds,gram,cm
L=16
N=256
h=L/N
hb=h/4
ip=[(2:N),1]
im=[N,(1:(N-1))]
K=25600000
rho=1    
mu=0.010518
tmax=4
dt=0.064
clockmax=ceil(tmax/dt)
flow = zeros(clockmax,7);
tube;
%initial velocity
u=zeros(N,N,2);



%[vX,vY] = meshgrid (0:h:L-h);
%quiver (vY,vX,bf(:,:,1),bf(:,:,2),5);

    
xgrid=zeros(N,N);
ygrid=zeros(N,N);
for j=0:(N-1)
  xgrid(j+1,:)=j*h;
  ygrid(:,j+1)=j*h;
end

set(gcf,'double','on')
%contour(xgrid,ygrid,vorticity,values)
%hold on

plot(X(:,1),X(:,2),'ko')

axis([0,L,0,L])
%caxis(valminmax)
axis equal
axis manual
drawnow
hold off

