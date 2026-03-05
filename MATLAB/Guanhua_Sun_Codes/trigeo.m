%triangle 
X=zeros(1,2);
xO=2;
yO=Ly/2;
%two angles
halfangle=pi/180*hangle;
l_tri=0.5;
num=1;
x=xO;
y=yO;
X(1,:)=[x,y];
for k1=1:floor(l_tri/hb)
  num=num+1;
  x=x+cos(halfangle)*hb;
  y=y+sin(halfangle)*hb;
  X(num,1)=x;
  X(num,2)=y;

  
end

lc=2*l_tri*sin(halfangle);

% x=xO;
% y=yO;
 xA=x;
 yA=y;
baseangle=(pi/2-halfangle);
for k2=1:ceil(lc/hb)
  num=num+1;
  y=y-hb;
  X(num,1)=x;
  X(num,2)=y;
  
  
end

 xB=x;
 yB=y;
 
for k3=1:floor(l_tri/hb)
  
  num=num+1;
  x=x+cos(pi/2+baseangle)*hb;
  y=y+sin(pi/2+baseangle)*hb;
  X(num,1)=x;
  X(num,2)=y;
  

  
end


Z=X;
Nb=size(X,1);
Zb=Nb;

%find the centroid of the triangle
center=[xA+xB+xO,yA+yB+yO]/3;
