plist=zeros(37,11);
hh=1;
num_start=10000;
u_up=5;
l_tri=0.5;
% figure
%[10:5:20,23:25,30,35,39:41,45:5:60];
for i = [30]

    Tlist=zeros(2,1);
    jj=0;
    for j = [0:36]
    
    h=5*j;
    Tlist(jj+1,1)=h;

        %name = sprintf('smallera%u',i);
	sfolder=sprintf('4.5*6/a%u',i);
        file = sprintf('a%ub%u.mat',i,h);
    
    fullname = fullfile(sfolder,file);
    load (fullname);
     num=size(T(num_start:end,2),1);
    
    Tlist(jj+1,2)=sum(T(num_start:end,2))/num;
    Tlist(jj+1,2)=Tlist(jj+1,2)/(0.25*u_up^2*l_tri^2*sin(2*i*pi/180));
    
    % taverage
   % Tlist(jj+1,2)=t_mean(end,2);
   % Tlist(jj+1,3)=tmm(end,2);
    jj=jj+1;
    end
    plist(:,hh)=Tlist(:,2);
   %  plot (Tlist(:,1),Tlist(:,2),'o-','LineWidth',1.5);
   %   hold on
    % plot (Tlist(:,1),Tlist(:,3),'o-','LineWidth',1.5);
    % hold on 
     hh=hh+1;
end
%save ('smallera6570.mat','plist');
% pp=zeros(37,19);
% for j = 1:37
%     j=38-j;
%     pp(j,:)=sum(plist(1:j,:));
% end
%pp=-pp;
%legend ('a=30','a=35','a=39','a=40','a=41')
   %legend ('a=10','a=15','a=20','a=25','a=30','a=35','a=40','a=45','a=50',...
   %    'a=55','a=60','a=65','a=70','a75','a=80');
