function U=interp(u,X)
global Nb h;
global N;
U=zeros(Nb,2);
for k=1:Nb
  s=X(k,:)/h;
  i=floor(s);
  r=s-i;
  i1=mod((i(1)-1):(i(1)+2),N)+1; 
  i2=mod((i(2)-1):(i(2)+2),N)+1;
  w = twoD_delta_discrete_delta(r(1), r(2));
  % w=phi1(r(1)).*phi2(r(2));
  U(k,1)=sum(sum(w.*u(i1,i2,1)));
  U(k,2)=sum(sum(w.*u(i1,i2,2)));
end

