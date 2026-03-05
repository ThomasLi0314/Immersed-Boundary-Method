clear all
%hangle=60;
lfolder='linesofaction';
mkdir(lfolder);
global h hb center
% figure
centerlist=[];
kk=0;
pcmean_list=[];
for hangle = [30]
kk=kk+1;
oangle=hangle/180*pi;
n1=[sin(oangle) cos(oangle)]';
n2=[1 0]';
n3=[sin(oangle) -cos(oangle)]';
folder=sprintf('smallera%u',hangle);
pdt=0.05;
pname=sprintf('pmean%u',hangle);
lname=sprintf('linesa%u',hangle);
lname=fullfile(lfolder,lname);
pname=fullfile(folder,pname);
pclist=zeros(1,4);
linesoa=[];
%pcmean_list=zeros(1,5);
for anglenum = 0	
   angle = anglenum/180*pi;
    trini
    pmean=zeros(Nb,1);
    rotation=[cos(angle) -sin(angle);sin(angle) cos(angle)];
    Xt=X-center;
    Xt=Xt*rotation;
    X=Xt+center;
    n1=(n1'*rotation)'
    n2=(n2'*rotation)'
    n3=(n3'*rotation)';
    num_start=20;
    num_end=120;
    for i=num_start:num_end
        p_start=5;
        p_end=length(X)-5;
        t=i*pdt;
        name=sprintf('a%ub%u',hangle,anglenum);
        sfolder=fullfile(folder,name);
        bfname=sprintf('a%ub%ut%g.mat',hangle,anglenum,t);
        bfname=fullfile(sfolder,bfname);
        load(bfname);
	linesoa(i,1:2)=sum(bf,1);
	linesoa(i,3)=torque(X,bf);
        p1=bf(1:k1,:)*n1;
        p2=bf(k1+1:k1+k2+1,:)*n2;
        p3=bf(k1+k2+2:end,:)*n3;
        p=[p1;p2;p3];
        pmean=pmean+p;
        pc=sum(X.*p,1)/sum(p,1);
	pc1=sum(X(1:k1,:).*p1);
	pc3=sum(X(k1+k2+2:end,:).*p3);
	pcs=(pc1+pc3)/sum(p1+p3);
%	pcs=sum(X(1:k1+k2+1,:).*p(1:k1+k2+1))/sum(p(1:k1+k2+1,:));	
     %   plot(X(:,1),X(:,2),'k.-');
     %   hold on
     %   plot (center(1),center(2),'ro');
     %   hold on
     %   plot(pc(1),pc(2),'k.')
     %   drawnow
     %   axis equal
        pclist(i,:)=[pc pcs];
    end
	
    pmean=pmean/(num_end-num_start+1);
    pcmean=sum(pclist)/(num_end-num_start+1);
    pcmean2=sum(X.*pmean,1)/sum(pmean,1);
    pcmean_list(kk,1:7)=[hangle center pcmean];
    %figure
    %plot(X(:,1),X(:,2),'k.-');
    %hold on
    centerlist=[centerlist;center];
    %plot(center(1),center(2),'o');
    %hold on
%     plot(pcmean2(1)-center(1),pcmean2(2)-center(2),'o');
%     drawnow
%     hold on
    %hold on
    %plot(pcmean(1),pcmean(2),'bo');
   % axis equal
	save(lname,'linesoa');
	save(pname,'pmean');
end
end
