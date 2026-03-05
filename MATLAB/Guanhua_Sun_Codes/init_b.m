%init_b.m
%This script initializes the array  b 
%that is used in the fluid solver to extract pressure output
b=zeros(Nx,Ny,2);

for m1=0:(Nx-1)
  for m2=0:(Ny-1)
    if~(((m1==0)|(m1==Nx/2))&((m2==0)|(m2==Ny/2)))
      t=[2*pi/Nx;2*pi/Ny].*[m1;m2];
      s=sin(t);
      c = - (i*rho*h/dt)/(s(1)^2+s(2)^2);
      b(m1+1,m2+1,1) = c*s(1);
      b(m1+1,m2+1,2) = c*s(2);
    end
  end
end