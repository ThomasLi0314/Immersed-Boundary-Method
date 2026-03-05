%initialize.m
L=1.0
N=64
h=L/N
ip=[(2:N),1]
im=[N,(1:(N-1))]
Nb=ceil(pi*(L/2)/(h/2))
dtheta=2*pi/Nb
kp=[(2:Nb),1]
km=[Nb,(1:(Nb-1))]
K=1
rho=1
mu=0.01
tmax=4
dt=0.01
clockmax=ceil(tmax/dt)

for k=0:(Nb-1)
  theta=k*dtheta;
  X(k+1,1)=(L/2)+(L/4)*cos(theta);
  X(k+1,2)=(L/2)+(L/4)*sin(theta);
end

u=zeros(N,N,2);
for j1=0:(N-1)
  x=j1*h;
  u(j1+1,:,2)=sin(2*pi*x/L);
end

vorticity=(u(ip,:,2)-u(im,:,2)-u(:,ip,1)+u(:,im,1))/(2*h);
dvorticity=(max(max(vorticity))-min(min(vorticity)))/5;
values= (-10*dvorticity):dvorticity:(10*dvorticity);
valminmax=[min(values),max(values)];
xgrid=zeros(N,N);
ygrid=zeros(N,N);
for j=0:(N-1)
  xgrid(j+1,:)=j*h;
  ygrid(:,j+1)=j*h;
end

set(gcf,'double','on')
contour(xgrid,ygrid,vorticity,values)
hold on
plot(X(:,1),X(:,2),'ko')
axis([0,L,0,L])
caxis(valminmax)
axis equal
axis manual
drawnow
hold off

