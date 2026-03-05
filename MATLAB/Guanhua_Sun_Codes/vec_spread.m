function f=vec_spread(F,X)
% spread F to grid
global h Nx Ny Nb hb;

c=hb/(h*h);
f=zeros(Nx,Ny,2);

s=X/h; % Get body position relative to grid
i=floor(s);
r=s-i;
w=vec_phi1(r(:,1)).*vec_phi2(r(:,2));%Evaluate delta function
w = permute(w, [1,3,2]);


for k=1:Nb
  i1=mod((i(k,1)-1):(i(k,1)+2),Nx)+1; %Find affected cells
  i2=mod((i(k,2)-1):(i(k,2)+2),Ny)+1;
  ww = w(:,:,k);
  f(i1,i2,1)=f(i1,i2,1)+(c*F(k,1))*ww; %Spread force to fluid
  f(i1,i2,2)=f(i1,i2,2)+(c*F(k,2))*ww;
end

