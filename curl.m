function curl(u,v)
% z component of vorticity

global dx
global dy
global N
global M
global vort

vort(2:N-1,2:M-1) = (v(2:N-1,3:M)-v(2:N-1,1:M-2))/dx/2-(u(3:N,2:M-1)-u(1:N-2,2:M-1))/dy/2;

end
