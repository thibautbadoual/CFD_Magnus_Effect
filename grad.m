function grad(phi)
% Compute the gradient of phi

global dx;
global dy;
global gradphix;
global gradphiy;
global N;
global M;

gradphix(2:N-1,2:M-1) = (phi(2:N-1,3:M)-phi(2:N-1,1:M-2))/dx/2;
gradphiy(2:N-1,2:M-1)= (phi(3:N,2:M-1)-phi(1:N-2,2:M-1))/dy/2;

end

