% ib2D.m
% This script is the main program.

global dt Nb N h rho mu ip im a hb;
global kp km dtheta K flow tw u;
initialize2
init_a

for clock=1:clockmax
    t = clock*dt;
    %bf = bforce(fpos,20*pi,t);
    XX=X+(dt/2)*interp(u,X);  
    ff=spread(Force(XX,Z),XX);
    u(N/2,N/2+8:N/2+24,1)= 1;
    [u,uu]=fluid(u,ff);
    %if mod(clock,10)==0   
    %    flow(clock,1)=t;
    %  flow(clock,2)=flux(20*[0.5 0.3125],1);
    %   flow(clock,3)=flux(20*[0.8125 0.25],2);
    %   flow(clock,4)=flux(20*[0.5 0.0625],1);
    %end
    X=X+dt*interp(uu,XX); 
  %animation:
  
  vorticity=(u(ipx,:,2)-u(imx,:,2)-u(:,ipy,1)+u(:,imy,1))/(2*h);
  contour(xgrid,ygrid,vorticity,30);
  

  vu = u(:,:,1);
  vv = u(:,:,2);
  fu = ff(:,:,1);
  fv = ff(:,:,2);
  %quiver (xgrid,ygrid,u(:,:,1),u(:,:,2),3);
  
  %streamline(vX,vY,vu,vv,0.5,0.1);
  title (["t=",t]);
  
  hold on
 
  plot(X(:,1),X(:,2),'ko')

  %caxis(valminmax)
  axis equal
  axis manual
  drawnow
  hold off
  %save ('test.mat',flow);
end
