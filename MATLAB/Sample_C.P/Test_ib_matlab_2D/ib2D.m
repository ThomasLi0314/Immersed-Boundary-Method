% ib2D.m
% This script is the main program.

global dt Nb N h rho mu ip im a;
global kp km dtheta K;
initialize
init_a

for clock=1:clockmax
  % XX=X+(dt/2)*interp(u,X);
  XX= X + dt / 2 * interpolation(u, X, Nb, N, N, h, h);
  % ff = spread(Force(XX), XX);
  ff=spreadforce(Force(XX),XX,dtheta,h,h,Nb,N,N);
  [u,uu] = navier_stokes_solver_pIB(u, ff, a, ip, im, ip, im, h, h, dt, mu, rho);
  % [u,uu]=fluid(u,ff);
  % X=X+dt*interp(uu,XX);
  X = X + dt * interpolation(uu, XX, Nb, N, N, h, h);
   
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

