%flow meter
%must measure at prue vertical or horizontal postion and always choose at
%the bottom or the left
%which is not applicable around corners
function Q = flux(f,i,width)
global N h tw u L;


nodes = floor(width/h);
pN =ceil(N*f/16);
xN=pN(1,1);
yN=pN(1,2);
if i==1
    Nyt = yN+nodes;
    x = 16*xN/N*ones(1,nodes);
    y = 16/N*linspace(yN,Nyt,nodes);
    scatter(x,y);
    hold on
    Q = h*sum(u(xN,yN:Nyt,1));
elseif i==2
    Nxr = xN+nodes;
    
    y2 = 16*yN/N*ones(1,nodes);
    x2 = 16/N*linspace(xN,Nxr,nodes);
    scatter(x2,y2);
    Q = h*sum(u(xN:Nxr,yN,2));
    hold on
end   

