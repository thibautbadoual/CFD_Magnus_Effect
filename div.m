function divergence=div(u,v)
% Computes the divergence

global dx
global dy
global N
global M

divergence=zeros(N,M);
divergence(2:N-1,2:M-1) = (u(2:N-1,3:M)-u(2:N-1,1:M-2))/dx/2+(v(3:N,2:M-1)-v(1:N-2,2:M-1))/dy/2;

end