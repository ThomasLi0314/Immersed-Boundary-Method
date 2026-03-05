%set_plate
xO=2;
yO=2.75;
%two angles
num=1;
x=xO;
y=yO;
X(1,:)=[x,y];

for k1=1:floor(1/hb)
    
  num=num+1;
  y=y-hb;
  X(num,1)=x;
  X(num,2)=y;

  
end

Z=X;
Nb=size(X,1);
Zb=Nb;

%find the centroid of the triangle
center=sum(X)/Nb;