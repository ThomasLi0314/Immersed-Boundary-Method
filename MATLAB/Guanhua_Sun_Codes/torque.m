%calculate the torque on the body by intergating pressure on the surface
function tor=torque(X,force)
global h hb center;

r=X-center;
%f=force(Yind(:,1),Yind(:,2),:);
f=force.*hb;
crossf=fliplr(f);
crossf(:,2)=-crossf(:,2);
tor=r.*crossf;
tor=sum(sum(tor,2));

%pforce = pres*hb;
%r=((x(:,1)-center(1)).^2+(x(:,2)-center(2)).^2).^(1/2);



%for i = 1:Nb
%    if i<=k1
%        normal = [yA-yO,-xA+xO];
%    elseif i>k1 && i<=k2
%        normal = [yO-yB,xB-xO];
%    elseif i>k2 
%        normal = [yB-yA,-xB+xA];
%    end
%    angle(i)=((x(i,:)-center)*normal.')/(norm(x(i,:)-center)*norm(normal));
%end

%angle=(ones(Nb,1)-angle.^2).^(1/2);

