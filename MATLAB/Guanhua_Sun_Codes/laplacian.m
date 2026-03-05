function w=laplacian(u)
global imx ipx imy ipy h;
w=(u(ipx,:,:)+u(imx,:,:)+u(:,ipy,:)+u(:,imy,:)-4*u)/(h*h);

% im = [N,1:(N-1)] = circular version of i-1
% ip = [2:N,1]     = circular version of i+1
% N  = number of points in each space direction

