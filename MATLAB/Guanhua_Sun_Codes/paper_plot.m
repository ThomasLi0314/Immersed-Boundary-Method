%plot

%vorticity-zoomed out
% figure
% 
% dvorticity=(max(max(vorticity))-min(min(vorticity)))/60;
% values= (-15*dvorticity):dvorticity/2:(15*dvorticity);
% valminmax=[min(values),max(values)];
% %set(gcf,'double','on')
% contourf(xgrid,ygrid,vorticity,values,'Lines','None')
% caxis([-150,150])
% axis([0 6 0 4.5])
% axis equal
% colormap jet
%  hold on
%  plot(X(:,1),X(:,2),'ko','MarkerSize',2);

%pressure-zoomed in
figure 
p_average=p_average/12.5;
dp=(max(max(p_average))-min(min(p_average)))/40;
values= (-25*dp):dp/10:(25*dp);
valminmax=[min(values),max(values)];
%set(gcf,'double','on')
contourf(xgrid,ygrid,p_average,values,'Lines','None')
caxis([-3,2])
axis equal
axis ([1.6 2.8 1.6 2.8])
colormap jet
hold on
plot(X(:,1),X(:,2),'ko','MarkerSize',4);
colorbar('LineWidth',1,'TickLength',0.02);

%velocity_average_zoomedin
% figure
% xp=[1:8:512];
% yp=[1:8:384];
% plot(X(:,1),X(:,2),'ko','MarkerSize',2);
% hold on
% quiver (xgrid(xp,yp),ygrid(xp,yp),u_average(xp,yp,1),u_average(xp,yp,2),2,'LineWidth',2)
% axis equal

