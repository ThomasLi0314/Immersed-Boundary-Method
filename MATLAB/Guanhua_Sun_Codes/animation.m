hangle=35;
anglenum=0;
trini;
pdt=0.05;
folder=sprintf('smallera%u',hangle);
name = sprintf('a%ub%u',hangle,anglenum);
T=8;
clockmax=ceil(T/pdt);
oangle=hangle/180*pi;
n1=[sin(oangle) -cos(oangle)]';
n2=[-1 0]';
n3=[sin(oangle) cos(oangle)]';

angle=anglenum/180*pi;
rotation=[cos(angle) -sin(angle);sin(angle) cos(angle)];
n1=(n1'*rotation)'
n2=(n2'*rotation)'
n3=(n3'*rotation)';
sfolder=fullfile(folder,name);
for clock=1:clockmax
         pt=clock*pdt;
         bfname=sprintf('a%ub%ut%g.mat',hangle,anglenum,pt);
         bfname=fullfile(sfolder,bfname);
         xname=sprintf('xa%ub%ut%g.mat',hangle,anglenum,pt);
         xname=fullfile(sfolder,xname);
         uname=sprintf('ua%ub%ut%g.mat',hangle,anglenum,pt);
         uname=fullfile(sfolder,uname);
         load(xname);
         %load(uname);
         load(bfname);
         p1=bf(1:k1,:)*n1;
         p2=bf(k1+1:k1+k2+1,:)*n2;
         p3=bf(k1+k2+2:end,:)*n3;
         p=[p1;p2;p3];
         pc=sum(X.*p,1)/sum(p,1);
         center=sum(X)/length(X);
         plot(X(:,1),X(:,2),'k.-');
         hold on
         plot (center(1),center(2),'ro');
         hold on
         plot(pc(1),pc(2),'k.')
        % hold on
         %vorticity=(u(ipx,:,2)-u(imx,:,2)-u(:,ipy,1)+u(:,imy,1))/(2*h);
         %contour(xgrid,ygrid,vorticity,50);
         
         axis equal
         %axis ([0 Lx 0 Ly])
        
         drawnow
         hold off
end