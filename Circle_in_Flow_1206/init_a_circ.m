%init_a.m
%This script initializes the array  a  
%that is used in the fluid solver

global a;
a=zeros(Nx,Ny,2,2);
for m1=0:(Nx-1)
  for m2=0:(Ny-1)
    a(m1+1,m2+1,1,1)=1;
    a(m1+1,m2+1,2,2)=1;
  end
end

for m1=0:(Nx-1)
  for m2=0:(Ny-1)
    if~(((m1==0)|(m1==Nx/2))&((m2==0)|(m2==Ny/2)))
      t=2 * pi * [m1 / Nx; m2 / Ny];
      s=sin(t);
      ss=(s*s')/(s'*s);
%     a(m1+1,m2+1,:,:)=a(m1+1,m2+1,:,:)-(s*s')/(s'*s);
      a(m1+1,m2+1,1,1)=a(m1+1,m2+1,1,1)-ss(1,1);
      a(m1+1,m2+1,1,2)=a(m1+1,m2+1,1,2)-ss(1,2);
      a(m1+1,m2+1,2,1)=a(m1+1,m2+1,2,1)-ss(2,1);
      a(m1+1,m2+1,2,2)=a(m1+1,m2+1,2,2)-ss(2,2);
    end
  end
end

for m1=0:(Nx-1)
  for m2=0:(Ny-1)
    t=pi * [m1 / Nx, m2 / Ny];
    s=sin(t);
    matrix = s'*s;
    scaler = dt / 2 * (mu / rho) * 4 / (dx * dx);
    a(m1 +1 , m2 + 1, 1, 1) = a(m1 +1 , m2 + 1, 1, 1)/ (1 + scaler * matrix(1,1));
    a(m1 +1 , m2 + 1, 1, 2) = a(m1 +1 , m2 + 1, 1, 2)/ (1 + scaler * matrix(1,2));
    a(m1 +1 , m2 + 1, 2, 1) = a(m1 +1 , m2 + 1, 2, 1)/ (1 + scaler * matrix(2,1));
    a(m1 +1 , m2 + 1, 2, 2) = a(m1 +1 , m2 + 1, 2, 2)/ (1 + scaler * matrix(2,2)); 
    % a(m1+1,m2+1,:,:)=a(m1+1,m2+1,:,:)...
    %                 /(1+(dt/2)*(mu/rho)*(4/(dx*dx))*(s'*s));
  end
end

