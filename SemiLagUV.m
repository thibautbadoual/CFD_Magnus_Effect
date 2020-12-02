function SemiLagUV(u,v,dt)

%=================================================
%
% Advection routine based on Semi-Lagrangian reconstruction
% using bilinear interpolation
%
%=================================================


global dx;
global dy;
global uadv;
global vadv;

N = size(u,1);
M = size(u,2);

% Matrices where 1 is right, 0 is left or center
Mx2 = sign(sign(u(2:N-1,2:M-1)) + ones(N-2,M-2));
Mx1 = ones(N-2,M-2)-Mx2;

% Matrices where 1 is up, 0 is down or center
My2 = sign(sign(v(2:N-1,2:M-1)) + ones(N-2,M-2));
My1 = ones(N-2,M-2)-My2;

% Matrices of absolute values for u and v
au = abs(u(2:N-1,2:M-1));
av = abs(v(2:N-1,2:M-1));

% Matrices of coefficients respectively central, external, same x, same y
Cc = (dx*ones(N-2,M-2)-au*dt).*(dy*ones(N-2,M-2)-av*dt)/dx/dy;
Ce = dt*dt*au.*av/dx/dy;
Cmx = (dx*ones(N-2,M-2)-dt*au).*av*dt/dx/dy;
Cmy = dt*au.*(dy*ones(N-2,M-2)-dt*av)/dx/dy;

% Computes the advected velocities
uadv(2:N-1,2:M-1)=Cc.*u(2:N-1,2:M-1)+...
    Ce.*(Mx1.*My1.*u(3:N,3:M)+Mx1.*My2.*u(1:N-2,3:M)+Mx2.*My1.*u(3:N,1:M-2)+Mx2.*My2.*u(1:N-2,1:M-2))+...
    Cmx.*(My1.*u(3:N,2:M-1)+My2.*u(1:N-2,2:M-1))+...
    Cmy.*(Mx1.*u(2:N-1,3:M)+Mx2.*u(2:N-1,1:M-2));

vadv(2:N-1,2:M-1)=Cc.*v(2:N-1,2:M-1)+...
    Ce.*(Mx1.*My1.*v(3:N,3:M)+Mx1.*My2.*v(1:N-2,3:M)+Mx2.*My1.*v(3:N,1:M-2)+Mx2.*My2.*v(1:N-2,1:M-2))+...
    Cmx.*(My1.*v(3:N,2:M-1)+My2.*v(1:N-2,2:M-1))+...
    Cmy.*(Mx1.*v(2:N-1,3:M)+Mx2.*v(2:N-1,1:M-2));

end
