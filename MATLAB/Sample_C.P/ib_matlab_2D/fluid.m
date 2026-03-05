function [uuu,uu]=fluid(u,ff)
global a dt rho mu;
w=u-(dt/2)*skew(u)+(dt/(2*rho))*ff;
w=fft(w,[],1);
w=fft(w,[],2);
uu(:,:,1)=a(:,:,1,1).*w(:,:,1)+a(:,:,1,2).*w(:,:,2);
uu(:,:,2)=a(:,:,2,1).*w(:,:,1)+a(:,:,2,2).*w(:,:,2);
uu=ifft(uu,[],2);
uu=real(ifft(uu,[],1));
w=u-dt*skew(uu)+(dt/rho)*ff+(dt/2)*(mu/rho)*laplacian(u);
w=fft(w,[],1);
w=fft(w,[],2);
uuu(:,:,1)=a(:,:,1,1).*w(:,:,1)+a(:,:,1,2).*w(:,:,2);
uuu(:,:,2)=a(:,:,2,1).*w(:,:,1)+a(:,:,2,2).*w(:,:,2);
uuu=ifft(uuu,[],2);
uuu=real(ifft(uuu,[],1));
