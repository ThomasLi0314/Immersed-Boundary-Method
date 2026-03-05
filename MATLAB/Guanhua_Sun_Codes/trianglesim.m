%simscircular.m
% This script is the main program.
clear;
clc;
close all;
global dt Nb Nx Ny h rho mu ipx ipy imx imy a b hb center;
global kp km dtheta K flow tw u pw ppp;
%hangle=20;
T=[];
for hangle = [30]
    folder=sprintf('K1e6convergence%u',hangle);
    mkdir (folder);
    for anglenum=[25]
    angle = anglenum/180*pi;
    trini
    init_a
    init_b
    rotation=[cos(angle) -sin(angle);sin(angle) cos(angle)];
    Xt=X-center;
    Xt=Xt*rotation;
    X=Xt+center;
    Z=X;
    
    name = sprintf('a%ub%u',hangle,anglenum);
     sfolder=fullfile(folder,name);
    mkdir(sfolder);
    p_total=zeros(Nx,Ny);
    omega_total=zeros(Nx,Ny);
    u_total=zeros(Nx,Ny,2);
    for clock=1:clockmax
        t = clock*dt
        u(1:Nx/32,:,1)=5;
	    u(1:Nx/32,:,2)=0;
        XX=X+(dt/2)*vec_interp(u,X);
        bf=Force(XX,Z);
        ff=vec_spread(bf,XX);
        %pp=bf(:,1)*cos(angle)+bf(:,2)*sin(angle);
        %center_p=sum(X.*pp)/sum(pp);
        T(clock,1)=t;
        T(clock,2)=torque(XX,Force(XX,Z));      
        T(clock,3:4)=sum(bf)*hb;      
        T(clock,5)=norm(max(XX-Z));
        
        [u,uu]=fluid(u,ff);
%         if clock>6000
%             vorticity=(u(ipx,:,2)-u(imx,:,2)-u(:,ipy,1)+u(:,imy,1))/(2*h);
%             p_total=p_total+ppp;
%             u_total=u_total+u;
%             omega_total=omega_total+vorticity;
%             p_average=p_total/(clock-6000);
%             omega_average=omega_total/(clock-6000);
%             u_average=u_total/(clock-6000);
%             if mod(clock,500)==0
%                 p_average_name=sprintf('p_average%g.mat',t);
%                 p_average_name=fullfile(sfolder,p_average_name);
%                 omega_average_name=sprintf('o_average%g.mat',t);
%                 omega_average_name=fullfile(sfolder,omega_average_name);
%                 u_average_name=sprintf('u_average%g.mat',t);
%                 u_average_name=fullfile(sfolder,u_average_name);
%                 save (p_average_name,'p_average');
%                 save (omega_average_name,'omega_average');
%                 save (u_average_name,'u_average');
%             end
%          end  
%         dvorticity=(max(max(vorticity))-min(min(vorticity)))/60;
%         values= (-17*dvorticity):dvorticity/2:(17*dvorticity);
%         %valminmax=[min(values),max(values)];
%         set(gcf,'double','on')
%         contourf(xgrid,ygrid,vorticity,values)
%         %caxis(valminmax)
%         hold on
%         axis equal
        
        X=X+dt*vec_interp(uu,XX); 
        % plot(X(:,1),X(:,2),'ko')
        % axis([0,Lx,0,Ly])
        % drawnow
       %  if anglenum ==0 
       %  if mod(t,0.05)==0
       %      bfname=sprintf('a%ub%ut%g.mat',hangle,anglenum,t);
       %      bfname=fullfile(sfolder,bfname);
       %      xname=sprintf('xa%ub%ut%g.mat',hangle,anglenum,t);
       %      xname=fullfile(sfolder,xname);
       %      uname=sprintf('ua%ub%ut%g.mat',hangle,anglenum,t);
       %      uname=fullfile(sfolder,uname);
       %      save(bfname,'bf')
        %     save(xname,'X')
             %save(uname,'u')
             %save(bfname,'bf')
        % end
        % end

        
    end
    tname=fullfile(sfolder,name);
    save (tname,'T');
    end
end
