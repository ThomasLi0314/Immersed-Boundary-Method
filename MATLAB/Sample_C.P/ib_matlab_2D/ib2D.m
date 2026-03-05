% ib2D.m
% This script is the main program.

global dt Nb N h rho mu ip im a;
global kp km dtheta K;
initialize
init_a

for clock=1:clockmax
  XX=X+(dt/2)*interp(u,X);
  ff=spread(Force(XX),XX);
  [u,uu]=fluid(u,ff);
  X=X+dt*interp(uu,XX);
  
  %animation:
  vorticity=(u(ip,:,2)-u(im,:,2)-u(:,ip,1)+u(:,im,1))/(2*h);
  contour(xgrid,ygrid,vorticity,values)
  hold on
  plot(X(:,1),X(:,2),'ko')
  axis([0,L,0,L])
  caxis(valminmax)
  axis equal
  axis manual
  drawnow
  hold off
end

